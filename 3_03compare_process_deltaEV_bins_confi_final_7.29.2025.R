setwd("/Users/80030577/Desktop/HiC_analysis/EV_normalization")
load("normalized_ev.100k_multi_samples_06_03_2025.rda")

library(GenomicRanges)
library(IRanges)
library(S4Vectors)

normalized_ev_lists <- list(
  RBL = EV.rbl,
  LCL = EV.lcl,
  GCBC = EV.gcbc,
  MBC = EV.mbc,
  NBC = EV.nbc,
  PC = EV.pc
)

sample_names <- names(normalized_ev_lists)
chroms <- paste0("chr", c(1:22, "X"))
bin_size <- 100000  # 100kb

# Step 1: Create GRanges of Eigenvectors per Sample
ev_to_gr <- function(ev_by_chr, sample_name, bin_size = 100000) {
  gr_list <- lapply(names(ev_by_chr), function(chr) {
    ev_vec <- ev_by_chr[[chr]]
    starts <- seq(0, (length(ev_vec) - 1) * bin_size, by = bin_size)
    ends <- starts + bin_size
    GRanges(
      seqnames = chr,
      ranges = IRanges(start = starts + 1, end = ends),
      eigen = ev_vec,
      sample = sample_name
    )
  })
  unlist(GRangesList(gr_list))
}

gr_by_sample <- lapply(sample_names, function(sname) ev_to_gr(normalized_ev_lists[[sname]], sname))
names(gr_by_sample) <- sample_names

# Step 2: Compute ΔEV for Each Directed Sample Pair
pair_deltas <- list()
for (s1 in sample_names) {
  for (s2 in sample_names) {
    if (s1 != s2) {
      gr1 <- gr_by_sample[[s1]]
      gr2 <- gr_by_sample[[s2]]
      stopifnot(identical(seqnames(gr1), seqnames(gr2)))
      stopifnot(identical(start(gr1), start(gr2)))
      gr_delta <- gr1
      gr_delta$deltaEV <- gr2$eigen - gr1$eigen
      gr_delta$pair <- paste0(s1, "_to_", s2)
      pair_deltas[[paste0(s1, "_to_", s2)]] <- gr_delta
    }
  }
}

# Step 3: Calculate Global ΔEV SD Cutoff
sample_pairs <- combn(sample_names, 2, simplify = FALSE)
all_sds <- sapply(sample_pairs, function(pair) {
  sd(unlist(normalized_ev_lists[[pair[1]]]) - unlist(normalized_ev_lists[[pair[2]]]), na.rm = TRUE)
})
cutoff_value <- 2 * mean(all_sds, na.rm = TRUE)
cat("Global ΔEV cutoff (2×mean SD):", cutoff_value, "\n")

# Helper Functions
get_key <- function(gr) paste0(seqnames(gr), ":", start(gr))
filter_high_conf <- function(gr) abs(gr$deltaEV) > cutoff_value
group_consecutive_sign <- function(gr, sign_col = "sign") {
  signs <- mcols(gr)[[sign_col]]
  rle_signs <- Rle(signs)
  group_ids <- inverse.rle(list(values = seq_along(runLength(rle_signs)), lengths = runLength(rle_signs)))
  gr$group <- group_ids
  merged <- reduce(split(gr, gr$group))
  merged <- unlist(merged)
  merged$sign <- signs[match(names(merged), names(gr))]
  return(merged)
}

# Step 4: Define Comparison Function
compare_transitions <- function(pair1, pair2, label) {
  delta1 <- pair_deltas[[pair1]]
  delta2 <- pair_deltas[[pair2]]
  delta1$sign <- sign(delta1$deltaEV)
  delta2$sign <- sign(delta2$deltaEV)
  
  valid1 <- delta1[filter_high_conf(delta1)]
  valid2 <- delta2[filter_high_conf(delta2)]
  
  common_keys <- intersect(get_key(valid1), get_key(valid2))
  common1 <- valid1[get_key(valid1) %in% common_keys]
  common2 <- valid2[get_key(valid2) %in% common_keys]
  
  stopifnot(length(common1) == length(common2))
  same_sign <- common1$sign == common2$sign
  same_bins <- common1[same_sign]
  diff_bins <- common1[!same_sign]
  
  merged_same <- group_consecutive_sign(same_bins)
  merged_diff <- group_consecutive_sign(diff_bins)
  
  cat("\n===", label, "===\n")
  cat("High-confidence bins in", pair1, ":", length(valid1), "\n")
  cat("High-confidence bins in", pair2, ":", length(valid2), "\n")
  cat("Common bins:", length(common_keys), "\n")
  cat("→ Same direction:", length(same_bins), "bins in", length(merged_same), "regions\n")
  cat("→ Different direction:", length(diff_bins), "bins in", length(merged_diff), "regions\n")
  
  return(list(
    same = same_bins,
    diff = diff_bins,
    merged_same = merged_same,
    merged_diff = merged_diff
  ))
}

# Step 5: Run Comparisons
result_rbl_gcbc_vs_nbc_gcbc <- compare_transitions("RBL_to_GCBC", "NBC_to_GCBC", "RBL→GCBC vs. NBC→GCBC")
result_rbl_lcl_vs_nbc_lcl <- compare_transitions("RBL_to_LCL", "NBC_to_LCL", "RBL→LCL vs. NBC→LCL")
result_rbl_lcl_vs_rbl_nbc <- compare_transitions("RBL_to_LCL", "RBL_to_NBC", "RBL→LCL vs. RBL→NBC")

# Step 6: Save Output
save(
  result_rbl_gcbc_vs_nbc_gcbc,
  result_rbl_lcl_vs_nbc_lcl,
  result_rbl_lcl_vs_rbl_nbc,
  file = "merged_comparison_results.rda"
)

write.csv(as.data.frame(result_rbl_gcbc_vs_nbc_gcbc$diff), "diff_rbl_gcbc_vs_nbc_gcbc.csv")
write.csv(as.data.frame(result_rbl_lcl_vs_nbc_lcl$diff), "diff_rbl_lcl_vs_nbc_lcl.csv")
write.csv(as.data.frame(result_rbl_lcl_vs_rbl_nbc$diff), "diff_rbl_lcl_vs_rbl_nbc.csv")


# Step 7: Save a GRanges object directly to disk as .rds.
# .rds preserves all the genomic ranges info, metadata columns (like sign and deltaEV), and can be loaded back exactly as it was
# This is the best format for TxDb analysis because you can use it directly with findOverlaps() without any conversion.export_gr_rds <- function(gr, path_rds) {
  saveRDS(gr, file = path_rds)
}

#Starts a BED export function.

export_gr_bed <- function(gr, path_bed, name_prefix = "feat") {
  if (length(gr) == 0L) {
    write.table(data.frame(), file = path_bed, sep = "\t", quote = FALSE,
                row.names = FALSE, col.names = FALSE)
    return(invisible(NULL))
  }
  # Converts GRanges to a BED-compatible table:
  df <- data.frame(
    seqnames = as.character(GenomicRanges::seqnames(gr)),
    start    = GenomicRanges::start(gr) - 1L,          # BED is 0-based, half-open
    end      = GenomicRanges::end(gr),                 # end is 1-based in GRanges; OK for BED
    name     = paste0(name_prefix, "_", seq_along(gr)),
    score    = ".",
    strand   = ifelse(is.na(as.character(GenomicRanges::strand(gr))),
                      ".", as.character(GenomicRanges::strand(gr)))
  )
  # Optional: carry sign/deltaEV into the "name" column for quick glance
  # Enhancement: If your GRanges has a sign or deltaEV column, append those values to the name field in the BED.
  # This makes it easier to see direction (sign) and magnitude (deltaEV) in genome browsers.
  if ("sign" %in% names(mcols(gr))) {
    s <- as.integer(mcols(gr)$sign)
    df$name <- paste0(df$name, "_sign", s)
  }
  if ("deltaEV" %in% names(mcols(gr))) {
    # Round for compactness
    d <- round(as.numeric(mcols(gr)$deltaEV), 4)
    df$name <- paste0(df$name, "_dEV", d)
  }
  # Writes the BED file without headers (standard for BED format).
  write.table(df, file = path_bed, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)
}

# Convenience: export both formats with a base name
# Export both formats at once
export_gr_both <- function(gr, base) {
  export_gr_rds(gr, paste0(base, ".rds"))
  export_gr_bed(gr, paste0(base, ".bed"), name_prefix = basename(base))
}

# =========================
# Export SAME/DIFF bins + merged regions for each comparison
# =========================
# result_* objects exist from your compare_transitions() calls:
# - result_rbl_gcbc_vs_nbc_gcbc
# - result_rbl_lcl_vs_nbc_lcl
# - result_rbl_lcl_vs_rbl_nbc

# A little helper to export all pieces for a comparison result
# Export SAME/DIFF bins for each comparison
export_all_for_result <- function(res, tag) {
  # per-bin
  export_gr_both(res$same,        paste0("same_bins_",   tag))
  export_gr_both(res$diff,        paste0("diff_bins_",   tag))
  # merged regions
  export_gr_both(res$merged_same, paste0("regions_same_", tag))
  export_gr_both(res$merged_diff, paste0("regions_diff_", tag))
}

export_all_for_result(result_rbl_gcbc_vs_nbc_gcbc, "rbl_gcbc_vs_nbc_gcbc")
export_all_for_result(result_rbl_lcl_vs_nbc_lcl,   "rbl_lcl_vs_nbc_lcl")
export_all_for_result(result_rbl_lcl_vs_rbl_nbc,   "rbl_lcl_vs_rbl_nbc")

cat("✅ Exported SAME/DIFF bins and merged regions as both .rds (GRanges) and .bed (0-based).\n")
