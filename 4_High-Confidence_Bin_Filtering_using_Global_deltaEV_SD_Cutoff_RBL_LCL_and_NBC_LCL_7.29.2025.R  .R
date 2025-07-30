## use last script until step 3 line 91
## 3_High-Confidence_Bin_Filtering_using_Global_deltaEV_SD_Cutoff_7.29.2025.R

# Filter High-Confidence Bins
filter_high_confidence <- function(gr) {
  abs(gr$deltaEV) > cutoff_value
}

# Step 1: Create unique bin keys
get_key <- function(gr) paste0(seqnames(gr), ":", start(gr))

# Step 2: Prepare GRanges and deltaEV comparisons for RBL→LCL and NBC→LCL
# Assumes normalized_ev_lists and sample_names already exist
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
  return(unlist(GRangesList(gr_list)))
}

# Step 3: Generate GRanges for each sample
gr_by_sample <- lapply(sample_names, function(sname) {
  make_ev_gr(normalized_ev_lists[[sname]], sname)
})
names(gr_by_sample) <- sample_names

# Step 4: Compute ΔEV GRanges for RBL→LCL and NBC→LCL
gr_rbl <- gr_by_sample[["RBL"]]
gr_nbc <- gr_by_sample[["NBC"]]
gr_lcl <- gr_by_sample[["LCL"]]

delta_rbl_lcl <- gr_rbl
delta_rbl_lcl$deltaEV <- gr_lcl$eigen - gr_rbl$eigen

delta_nbc_lcl <- gr_nbc
delta_nbc_lcl$deltaEV <- gr_lcl$eigen - gr_nbc$eigen

# Step 5: Assign sign direction

delta_rbl_lcl$sign <- sign(delta_rbl_lcl$deltaEV)
delta_nbc_lcl$sign <- sign(delta_nbc_lcl$deltaEV)

# Step 6: Filter finite and high-confidence bins
valid_rbl_lcl <- delta_rbl_lcl[
  is.finite(delta_rbl_lcl$eigen) &
    is.finite(gr_by_sample[["LCL"]]$eigen) &
    abs(delta_rbl_lcl$deltaEV) > cutoff_value
]

valid_nbc_lcl <- delta_nbc_lcl[
  is.finite(delta_nbc_lcl$eigen) &
    is.finite(gr_by_sample[["LCL"]]$eigen) &
    abs(delta_nbc_lcl$deltaEV) > cutoff_value
]

# Step 7: Identify common high-confidence bins
rbl_lcl_keys <- get_key(valid_rbl_lcl)
nbc_lcl_keys <- get_key(valid_nbc_lcl)
common_keys_lcl <- intersect(rbl_lcl_keys, nbc_lcl_keys)

# Step 8: Subset GRanges to common bins
common_rbl_lcl <- valid_rbl_lcl[get_key(valid_rbl_lcl) %in% common_keys_lcl]
common_nbc_lcl <- valid_nbc_lcl[get_key(valid_nbc_lcl) %in% common_keys_lcl]

# Step 9: Compare sign directions for common bins
same_direction_lcl <- common_rbl_lcl$sign == common_nbc_lcl$sign
same_bins_lcl <- common_rbl_lcl[same_direction_lcl]
diff_bins_lcl <- common_rbl_lcl[!same_direction_lcl]

# Step 10: Group consecutive bins with same sign
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

# Step 11: Merge grouped bins
merged_same_lcl <- group_consecutive_sign(same_bins_lcl)
merged_diff_lcl <- group_consecutive_sign(diff_bins_lcl)

# Step 12: Summary
cat("RBL→LCL high-confidence finite bins:", length(rbl_lcl_keys), "\n")
cat("NBC→LCL high-confidence finite bins:", length(nbc_lcl_keys), "\n")
cat("Common high-confidence finite bins:", length(common_keys_lcl), "\n")
cat("Same-direction bins:", length(same_bins_lcl), "\n")
cat("Different-direction bins:", length(diff_bins_lcl), "\n")

# Optional output
head(as.data.frame(merged_same_lcl)[, c("seqnames", "start", "end", "sign")])
