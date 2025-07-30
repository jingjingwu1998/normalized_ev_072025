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

########################################################################
######## Step 1: Create GRanges of Eigenvectors per Sample ############
########################################################################

make_ev_gr <- function(ev_by_chr, sample_name, bin_size = 100000) {
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
  names(gr_list) <- names(ev_by_chr)
  return(unlist(GRangesList(gr_list)))
}

# Create GRanges for all samples
gr_by_sample <- lapply(sample_names, function(sname) {
  make_ev_gr(normalized_ev_lists[[sname]], sname)
})
names(gr_by_sample) <- sample_names

str(gr_by_sample, max.level = 1)

########################################################################
####### Step 2: Compute ΔEV for Each Directed Sample Pair #############
########################################################################

pair_deltas <- list()

for (s1 in sample_names) {
  for (s2 in sample_names) {
    if (s1 != s2) {
      gr1 <- gr_by_sample[[s1]]
      gr2 <- gr_by_sample[[s2]]
      
      # Confirm bin alignment
      stopifnot(identical(seqnames(gr1), seqnames(gr2)))
      stopifnot(identical(start(gr1), start(gr2)))
      
      delta <- gr2$eigen - gr1$eigen
      gr_delta <- gr1
      gr_delta$deltaEV <- delta
      gr_delta$pair <- paste0(s1, "_to_", s2)
      pair_deltas[[paste0(s1, "_to_", s2)]] <- gr_delta
    }
  }
}
str(pair_deltas, max.level = 1)

# Step 3: Calculate Global ΔEV SD Cutoff
sample_pairs_to_compare <- combn(sample_names, 2, simplify = FALSE)
all_pairwise_sd_of_differences <- c()
for (pair in sample_pairs_to_compare) {
  ev1 <- unlist(normalized_ev_lists[[pair[1]]])
  ev2 <- unlist(normalized_ev_lists[[pair[2]]])
  all_pairwise_sd_of_differences <- c(
    all_pairwise_sd_of_differences,
    sd(ev1 - ev2, na.rm = TRUE)
  )
}
cutoff_value <- 2 * mean(all_pairwise_sd_of_differences, na.rm = TRUE)
cat("Global ΔEV cutoff (2×mean SD):", cutoff_value, "\n")

# Step 4: Filter High-Confidence Bins
filter_high_confidence <- function(gr) {
  abs(gr$deltaEV) > cutoff_value
}

delta_rbl_lcl <- pair_deltas[["RBL_to_LCL"]]
delta_rbl_nbc <- pair_deltas[["RBL_to_NBC"]]
delta_rbl_lcl$sign <- sign(delta_rbl_lcl$deltaEV)
delta_rbl_nbc$sign <- sign(delta_rbl_nbc$deltaEV)

rbl_lcl_filt <- delta_rbl_lcl[filter_high_confidence(delta_rbl_lcl)]
length(rbl_lcl_filt)
rbl_nbc_filt <- delta_rbl_nbc[filter_high_confidence(delta_rbl_nbc)]
length(rbl_nbc_filt)

# Step 5: Find Common High-Confidence Bins
get_key <- function(gr) paste0(seqnames(gr), ":", start(gr))
rbl_lcl_keys <- get_key(rbl_lcl_filt)
rbl_nbc_keys <- get_key(rbl_nbc_filt)
common_keys <- intersect(rbl_lcl_keys, rbl_nbc_keys)
rbl_lcl_common <- rbl_lcl_filt[rbl_lcl_keys %in% common_keys]
rbl_nbc_common <- rbl_nbc_filt[rbl_nbc_keys %in% common_keys]
length(rbl_lcl_common)
length(rbl_nbc_common)
stopifnot(length(rbl_lcl_common) == length(rbl_nbc_common))

# Step 6: Compare Signs
same_sign <- rbl_lcl_common$sign == rbl_nbc_common$sign
same_bins <- rbl_lcl_common[same_sign]
diff_bins <- rbl_lcl_common[!same_sign]

# Step 7: Group Consecutive Bins by Sign
group_consecutive_sign <- function(gr, sign_col = "sign") {
  signs <- mcols(gr)[[sign_col]]
  rle_signs <- Rle(signs)
  group_ids <- inverse.rle(list(
    values = seq_along(runLength(rle_signs)),
    lengths = runLength(rle_signs)
  ))
  gr$group <- group_ids
  merged <- reduce(split(gr, gr$group))
  merged <- unlist(merged)
  merged$sign <- signs[match(names(merged), names(gr))]
  return(merged)
}

merged_same <- group_consecutive_sign(same_bins)
merged_diff <- group_consecutive_sign(diff_bins)

# Step 8: Output Results
cat("Common High-Confidence Bins:", length(common_keys), "\n")
cat("→ Same Direction:", length(same_bins), "bins in", length(merged_same), "regions\n")
cat("→ Different Direction:", length(diff_bins), "bins in", length(merged_diff), "regions\n")
