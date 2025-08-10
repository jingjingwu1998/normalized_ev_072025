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
