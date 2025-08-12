## ================================================================
## EV Δ comparisons + export SAME/DIFF bins/regions (RDS + BED + CSV)
## Robust v2: DF-driven merger, NA-safe, empty-safe, list->GRanges coercion
## ================================================================

## --- I/O roots (adjust if needed) ---
setwd("/Users/80030577/Desktop/HiC_analysis/EV_normalization")
load("normalized_ev.100k_multi_samples_06_03_2025.rda")
setwd("/Users/80030577/Desktop/HiC_analysis/EV_normalization/gene_annotation_latest_8.11.2025")

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(IRanges)
  library(S4Vectors)  # Rle
})

## ---------------------------------------------
## Inputs: per-sample eigenvectors (by chromosome)
## ---------------------------------------------
normalized_ev_lists <- list(
  RBL = EV.rbl,
  LCL = EV.lcl,
  GCBC = EV.gcbc,
  MBC = EV.mbc,
  NBC = EV.nbc,
  PC  = EV.pc
)

sample_names <- names(normalized_ev_lists)
bin_size <- 100000  # 100 kb

## ---------------------------------------------
## Utilities
## ---------------------------------------------

# Coerce anything GRanges-like (even nested lists) into a single GRanges
coerce_to_gr <- function(x) {
  if (inherits(x, "GRanges"))     return(x)
  if (inherits(x, "GRangesList")) return(unlist(x, use.names = FALSE))
  if (is.list(x)) {
    parts <- lapply(x, coerce_to_gr)
    parts <- parts[vapply(parts, function(g) inherits(g, "GRanges") && length(g) > 0, logical(1))]
    if (length(parts) == 0L) return(GenomicRanges::GRanges())
    return(do.call(c, parts))  # c(GRanges, GRanges, ...) -> GRanges
  }
  GenomicRanges::GRanges()
}

# Simple class/length checker
chk <- function(tag, x) {
  cat(sprintf("[chk] %s: class=%s; length=%s\n",
              tag, paste(class(x), collapse = "/"),
              if (is.list(x)) length(x) else length(x)))
}

## ---------------------------------------------
## Step 1: Create GRanges of eigenvectors/sample
## ---------------------------------------------
ev_to_gr <- function(ev_by_chr, sample_name, bin_size = 100000) {
  gr_list <- lapply(names(ev_by_chr), function(chr) {
    ev_vec <- ev_by_chr[[chr]]
    starts <- seq(0, (length(ev_vec) - 1) * bin_size, by = bin_size)
    ends   <- starts + bin_size
    GRanges(
      seqnames = chr,
      ranges   = IRanges(start = starts + 1, end = ends),
      eigen    = ev_vec,
      sample   = sample_name
    )
  })
  gr <- do.call(c, gr_list)
  stopifnot(inherits(gr, "GRanges"))
  gr
}

gr_by_sample <- lapply(sample_names, function(sname) ev_to_gr(normalized_ev_lists[[sname]], sname))
names(gr_by_sample) <- sample_names

## ---------------------------------------------------------
## Step 2: Compute ΔEV for each directed pair (s1 -> s2)
## ---------------------------------------------------------
pair_deltas <- list()
for (s1 in sample_names) {
  for (s2 in sample_names) {
    if (s1 != s2) {
      gr1 <- gr_by_sample[[s1]]
      gr2 <- gr_by_sample[[s2]]
      stopifnot(identical(seqnames(gr1), seqnames(gr2)))
      stopifnot(identical(start(gr1),    start(gr2)))
      gr_delta <- gr1
      gr_delta$deltaEV <- gr2$eigen - gr1$eigen
      gr_delta$pair    <- paste0(s1, "_to_", s2)
      pair_deltas[[paste0(s1, "_to_", s2)]] <- gr_delta
    }
  }
}

## ---------------------------------------------------------
## Step 3: Global ΔEV SD cutoff = 2 * mean SD over pairs
## ---------------------------------------------------------
sample_pairs <- combn(sample_names, 2, simplify = FALSE)
all_sds <- sapply(sample_pairs, function(pair) {
  sd(unlist(normalized_ev_lists[[pair[1]]]) - unlist(normalized_ev_lists[[pair[2]]]),
     na.rm = TRUE)
})
cutoff_value <- 2 * mean(all_sds, na.rm = TRUE)
cat("Global ΔEV cutoff (2×mean SD):", cutoff_value, "\n")

## -------------------
## Helper functions
## -------------------
get_key <- function(gr) paste0(seqnames(gr), ":", start(gr))

# High-confidence filter: drop NAs explicitly
filter_high_conf <- function(gr) {
  !is.na(gr$deltaEV) & abs(gr$deltaEV) > cutoff_value
}

# --- NEW: DF-driven merger (no split/do.call chains) ---
group_consecutive_sign <- function(gr, sign_vec, delta1_vec, delta2_vec) {
  gr <- coerce_to_gr(gr)
  if (length(gr) == 0L) return(gr)
  
  df <- data.frame(
    chr  = as.character(seqnames(gr)),
    start = as.integer(start(gr)),
    end   = as.integer(end(gr)),
    sign  = as.integer(sign_vec),
    d1    = as.numeric(delta1_vec),
    d2    = as.numeric(delta2_vec),
    stringsAsFactors = FALSE
  )
  
  # Keep only complete rows
  keep <- complete.cases(df[, c("chr","start","end","sign","d1","d2")])
  df <- df[keep, , drop = FALSE]
  if (nrow(df) == 0L) return(GRanges())
  
  # Order within chromosome
  o <- order(df$chr, df$start, df$end)
  df <- df[o, , drop = FALSE]
  
  # Identify breaks: new group if chr changes OR sign changes OR not adjacent
  chr_change  <- c(TRUE, df$chr[-1]  != df$chr[-nrow(df)])
  sign_change <- c(TRUE, df$sign[-1] != df$sign[-nrow(df)])
  adj_break   <- c(TRUE, df$start[-1] != (df$end[-nrow(df)] + 1L))
  new_group <- chr_change | sign_change | adj_break
  gid <- cumsum(new_group)
  
  # Aggregate per group
  split_idx <- split(seq_len(nrow(df)), gid)
  
  chr_vec   <- vapply(split_idx, function(ii) df$chr[ii[1]],   character(1))
  start_vec <- vapply(split_idx, function(ii) min(df$start[ii]), integer(1))
  end_vec   <- vapply(split_idx, function(ii) max(df$end[ii]),   integer(1))
  sign_out  <- vapply(split_idx, function(ii) df$sign[ii[1]],    integer(1))
  
  d1_mean <- vapply(split_idx, function(ii) mean(df$d1[ii], na.rm = TRUE), numeric(1))
  d1_min  <- vapply(split_idx, function(ii) min (df$d1[ii], na.rm = TRUE), numeric(1))
  d1_max  <- vapply(split_idx, function(ii) max (df$d1[ii], na.rm = TRUE), numeric(1))
  d2_mean <- vapply(split_idx, function(ii) mean(df$d2[ii], na.rm = TRUE), numeric(1))
  d2_min  <- vapply(split_idx, function(ii) min (df$d2[ii], na.rm = TRUE), numeric(1))
  d2_max  <- vapply(split_idx, function(ii) max (df$d2[ii], na.rm = TRUE), numeric(1))
  n_bins  <- vapply(split_idx, length, integer(1))
  
  out <- GRanges(
    seqnames = chr_vec,
    ranges   = IRanges(start = start_vec, end = end_vec),
    sign          = sign_out,
    deltaEV1_mean = d1_mean,
    deltaEV1_min  = d1_min,
    deltaEV1_max  = d1_max,
    deltaEV2_mean = d2_mean,
    deltaEV2_min  = d2_min,
    deltaEV2_max  = d2_max,
    n_bins        = n_bins
  )
  stopifnot(inherits(out, "GRanges"))
  out
}

## ---------------------------------------------
## Step 4: Pairwise comparison helper (robust)
## ---------------------------------------------
compare_transitions <- function(pair1, pair2, label) {
  delta1 <- pair_deltas[[pair1]]
  delta2 <- pair_deltas[[pair2]]
  
  delta1$sign <- sign(delta1$deltaEV)
  delta2$sign <- sign(delta2$deltaEV)
  
  valid1 <- delta1[filter_high_conf(delta1)]
  valid2 <- delta2[filter_high_conf(delta2)]
  
  keys1 <- get_key(valid1); keys2 <- get_key(valid2)
  common_keys <- intersect(keys1, keys2)
  
  cat("\n===", label, "===\n")
  cat("High-confidence bins in", pair1, ":", length(valid1), "\n")
  cat("High-confidence bins in", pair2, ":", length(valid2), "\n")
  cat("Common bins:", length(common_keys), "\n")
  
  if (length(common_keys) == 0L) {
    cat("→ Same direction: 0 bins in 0 regions\n")
    cat("→ Different direction: 0 bins in 0 regions\n")
    empty <- GRanges()
    return(list(same = empty, diff = empty, merged_same = empty, merged_diff = empty))
  }
  
  idx1 <- match(common_keys, keys1)
  idx2 <- match(common_keys, keys2)
  common1 <- valid1[idx1]
  common2 <- valid2[idx2]
  stopifnot(length(common1) == length(common2))
  
  deltaEV_pair1 <- as.numeric(common1$deltaEV)
  deltaEV_pair2 <- as.numeric(common2$deltaEV)
  s1 <- sign(deltaEV_pair1)
  
  # drop any NA that might remain
  keep <- !is.na(deltaEV_pair1) & !is.na(deltaEV_pair2) & !is.na(s1)
  if (!all(keep)) {
    common1 <- common1[keep]
    common2 <- common2[keep]
    deltaEV_pair1 <- deltaEV_pair1[keep]
    deltaEV_pair2 <- deltaEV_pair2[keep]
    s1 <- s1[keep]
  }
  
  if (length(common1) == 0L) {
    empty <- GRanges()
    cat("→ After NA removal: 0 bins overlap; returning empty.\n")
    return(list(same = empty, diff = empty, merged_same = empty, merged_diff = empty))
  }
  
  same_idx <- s1 == sign(deltaEV_pair2)
  diff_idx <- !same_idx
  
  same_bins <- common1[same_idx]
  diff_bins <- common1[diff_idx]
  
  if (length(same_bins) > 0) {
    same_bins$sign          <- as.integer(s1[same_idx])
    same_bins$deltaEV_pair1 <- as.numeric(deltaEV_pair1[same_idx])
    same_bins$deltaEV_pair2 <- as.numeric(deltaEV_pair2[same_idx])
    same_bins$pair1 <- rep(pair1, length(same_bins))
    same_bins$pair2 <- rep(pair2, length(same_bins))
  }
  if (length(diff_bins) > 0) {
    diff_bins$sign          <- as.integer(s1[diff_idx])
    diff_bins$deltaEV_pair1 <- as.numeric(deltaEV_pair1[diff_idx])
    diff_bins$deltaEV_pair2 <- as.numeric(deltaEV_pair2[diff_idx])
    diff_bins$pair1 <- rep(pair1, length(diff_bins))
    diff_bins$pair2 <- rep(pair2, length(diff_bins))
  }
  
  merged_same <- if (length(same_bins) > 0) {
    group_consecutive_sign(
      same_bins,
      sign_vec   = same_bins$sign,
      delta1_vec = same_bins$deltaEV_pair1,
      delta2_vec = same_bins$deltaEV_pair2
    )
  } else GRanges()
  
  merged_diff <- if (length(diff_bins) > 0) {
    group_consecutive_sign(
      diff_bins,
      sign_vec   = diff_bins$sign,
      delta1_vec = diff_bins$deltaEV_pair1,
      delta2_vec = diff_bins$deltaEV_pair2
    )
  } else GRanges()
  
  # force GRanges
  merged_same <- coerce_to_gr(merged_same)
  merged_diff <- coerce_to_gr(merged_diff)
  
  cat("→ Same direction:", length(same_bins), "bins in", length(merged_same), "regions\n")
  cat("→ Different direction:", length(diff_bins), "bins in", length(merged_diff), "regions\n")
  
  list(
    same        = same_bins,
    diff        = diff_bins,
    merged_same = merged_same,
    merged_diff = merged_diff
  )
}

## -----------------------
## Step 5: Run comparisons
## -----------------------
result_rbl_gcbc_vs_nbc_gcbc <- compare_transitions("RBL_to_GCBC", "NBC_to_GCBC", "RBL→GCBC vs. NBC→GCBC")
result_rbl_lcl_vs_nbc_lcl   <- compare_transitions("RBL_to_LCL",  "NBC_to_LCL",  "RBL→LCL vs. NBC→LCL")
result_rbl_lcl_vs_rbl_nbc   <- compare_transitions("RBL_to_LCL",  "RBL_to_NBC",  "RBL→LCL vs. RBL→NBC")

## ---------------------------------------------
## Step 6: Save R objects (RDA)
## ---------------------------------------------
save(
  result_rbl_gcbc_vs_nbc_gcbc,
  result_rbl_lcl_vs_nbc_lcl,
  result_rbl_lcl_vs_rbl_nbc,
  file = "merged_comparison_results.rda"
)

## ================================================================
## Export: SAME/DIFF bins and merged regions (RDS + BED + CSV)
## ================================================================

export_gr_rds <- function(gr, path_rds) saveRDS(coerce_to_gr(gr), file = path_rds)

# BED writer that always coerces to a single GRanges first
export_gr_bed <- function(gr, path_bed, name_prefix = "region") {
  gr <- coerce_to_gr(gr)
  if (length(gr) == 0L) {
    # Write a truly empty BED (0 lines). Many tools prefer no header/cols.
    writeLines(character(0L), con = path_bed)
    return(invisible(NULL))
  }
  df <- as.data.frame(gr)
  bed <- data.frame(
    chrom = as.character(df$seqnames),
    chromStart = as.integer(df$start) - 1L,  # BED is 0-based, half-open
    chromEnd = as.integer(df$end),
    name = sprintf("%s_%05d", name_prefix, seq_len(nrow(df))),
    score = 0,
    strand = if ("strand" %in% names(df)) as.character(df$strand) else ".",
    stringsAsFactors = FALSE
  )
  write.table(bed, file = path_bed, sep = "\t",
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  invisible(NULL)
}

export_gr_csv <- function(gr, path_csv) {
  gr <- coerce_to_gr(gr)
  if (length(gr) == 0L) {
    # Write an empty CSV with minimal headers for easy reading
    write.csv(data.frame(seqnames=character(),
                         start=integer(),
                         end=integer()),
              path_csv, row.names = FALSE)
  } else {
    write.csv(as.data.frame(gr), path_csv, row.names = FALSE)
  }
}

export_gr_both <- function(gr, base) {
  gr <- coerce_to_gr(gr)  # enforce here as well
  export_gr_rds(gr, paste0(base, ".rds"))
  export_gr_bed(gr, paste0(base, ".bed"), name_prefix = basename(base))
  export_gr_csv(gr, paste0(base, ".csv"))
}

export_all_for_result <- function(res, tag) {
  # Coerce everything up-front (handles empty and list cases)
  res$same        <- coerce_to_gr(res$same)
  res$diff        <- coerce_to_gr(res$diff)
  res$merged_same <- coerce_to_gr(res$merged_same)
  res$merged_diff <- coerce_to_gr(res$merged_diff)
  
  # Quick checks
  chk(paste0(tag, " same"),         res$same)
  chk(paste0(tag, " diff"),         res$diff)
  chk(paste0(tag, " regions_same"), res$merged_same)
  chk(paste0(tag, " regions_diff"), res$merged_diff)
  
  export_gr_both(res$same,         paste0("same_bins_",    tag))
  export_gr_both(res$diff,         paste0("diff_bins_",    tag))
  export_gr_both(res$merged_same,  paste0("regions_same_", tag))
  export_gr_both(res$merged_diff,  paste0("regions_diff_", tag))
  
  cat(sprintf("→ [%s] exported: same_bins=%d, diff_bins=%d, regions_same=%d, regions_diff=%d\n",
              tag, length(res$same), length(res$diff),
              length(res$merged_same), length(res$merged_diff)))
}

## -----------------------
## Step 7: Export all (empty-safe)
## -----------------------
export_all_for_result(result_rbl_gcbc_vs_nbc_gcbc, "rbl_gcbc_vs_nbc_gcbc")
export_all_for_result(result_rbl_lcl_vs_nbc_lcl,   "rbl_lcl_vs_nbc_lcl")
export_all_for_result(result_rbl_lcl_vs_rbl_nbc,   "rbl_lcl_vs_rbl_nbc")

cat("✅ Exported SAME/DIFF bins and merged regions as .rds, .bed, and .csv (empty-safe; ΔEV summaries included in regions)\n")
