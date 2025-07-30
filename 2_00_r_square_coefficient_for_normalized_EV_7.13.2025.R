# R square calculations

library(dplyr)
# --- Initial Setup ---
setwd('/Users/80030577/Desktop/HiC_analysis/EV_normalization')
getwd()
load("normalized_ev.100k_multi_samples_06_03_2025.rda")

normalized_ev_lists <- list(
  RBL = EV.rbl,
  LCL = EV.lcl,
  GCBC = EV.gcbc,
  MBC = EV.mbc,
  NBC = EV.nbc,
  PC = EV.pc
)

sample_names <- names(normalized_ev_lists)

chroms <- paste0('chr',c(1:22,'X'))

concatenated_normalized_evs <- list()
for (sample_name in sample_names) {
  current_sample_data <- normalized_ev_lists[[sample_name]]
  temp_vec <- c()
  for (chr in chroms) {
    # Assuming all chromosomes exist for all samples.
    temp_vec <- c(temp_vec, current_sample_data[[chr]])
  }
  concatenated_normalized_evs[[sample_name]] <- temp_vec
}

# Assign the concatenated data to my_eigen_data for subsequent analysis
my_eigen_data <- concatenated_normalized_evs
sample_lengths <- sapply(my_eigen_data, length)
print(sample_lengths)

# Verify data points of concatenated_normalized_evs

all_M_values_across_all_pairs_and_chrs <- c()
sample_pairs_to_compare <- combn(sample_names, 2, simplify = FALSE)

for (pair in sample_pairs_to_compare) {
  sample1_name <- pair[1]
  sample2_name <- pair[2]
  
  for (chr in chroms) {
    ev_values_sample1_chr <- normalized_ev_lists[[sample1_name]][[chr]]
    ev_values_sample2_chr <- normalized_ev_lists[[sample2_name]][[chr]]
    M_values_chr <- ev_values_sample2_chr - ev_values_sample1_chr # Sample2 - Sample1
    all_M_values_across_all_pairs_and_chrs <- c(all_M_values_across_all_pairs_and_chrs, M_values_chr)
  }
}
str(all_M_values_across_all_pairs_and_chrs) # 455,640
# 30,376 (bin numbers in each samples) X 15 (pairs number) = 455,640


####################################################################
### Core Data Preparation and Filtering
####################################################################

### Pre-calculate Chromosome Bin Mapping for Location Tracking
cat("--- Pre-calculating Chromosome Bin Mapping ---\n")
chr_bin_map <- data.frame(
  chromosome = character(),
  start_idx = numeric(),
  end_idx = numeric(),
  bin_in_chr_start = numeric(), # First bin number within that chromosome
  stringsAsFactors = FALSE
)
current_global_idx = 1
for (chr in chroms) {
  # Assuming all samples have the same binning structure, use the first sample to map
  chr_length = length(normalized_ev_lists[[sample_names[1]]][[chr]])
  if (chr_length > 0) {
    chr_bin_map <- rbind(chr_bin_map, data.frame(
      chromosome = chr,
      start_idx = current_global_idx,
      end_idx = current_global_idx + chr_length - 1,
      bin_in_chr_start = 1, # Bin numbering starts from 1 for each chr
      stringsAsFactors = FALSE
    ))
    current_global_idx = current_global_idx + chr_length
  }
}
cat("Chromosome bin mapping complete.\n")

# --- Compute the Global Cutoff Value for Highlighting ---
all_pairwise_sd_of_differences <- c()
sample_pairs_to_compare <- combn(sample_names, 2, simplify = FALSE)

for (pair in sample_pairs_to_compare) {
  sample1_name <- pair[1]
  sample2_name <- pair[2]
  
  ev_values_sample1 <- my_eigen_data[[sample1_name]]
  ev_values_sample2 <- my_eigen_data[[sample2_name]]
  
  M_values <- ev_values_sample1 - ev_values_sample2
  sd_of_current_M <- sd(M_values, na.rm = TRUE)
  all_pairwise_sd_of_differences <- c(all_pairwise_sd_of_differences, sd_of_current_M)
}

cutoff_value <- mean(all_pairwise_sd_of_differences, na.rm = TRUE)
cutoff_value <- 2 * cutoff_value # 2 times the mean SD as the cutoff
print(cutoff_value) #0.01885682
cat("-------------------------------------------------------------------\n")
cat("Calculated Cutoff Value for identifying heavily changed bins:", cutoff_value, "\n")
cat("-------------------------------------------------------------------\n\n")

num_samples_N <- length(sample_names) # 6 samples


# These matrices will store the R-squared values
r_squared_all_points_matrix <- matrix(NA, nrow = num_samples_N, ncol = num_samples_N,
                                      dimnames = list(sample_names, sample_names))
r_squared_highlighted_points_matrix <- matrix(NA, nrow = num_samples_N, ncol = num_samples_N,
                                              dimnames = list(sample_names, sample_names))

# Minimum number of points required for meaningful R-squared calculation
min_points_for_lm <- 2 # At least 2 points are needed for lm(), 3 or more for robust fitting

cat("\n--- Calculating Pairwise R-squared Coefficients ---\n")

# Loop through all sample pairs to calculate R-squared
for (i in 1:num_samples_N) {
  s_row_name <- sample_names[i]
  data_row <- my_eigen_data[[s_row_name]]
  
  for (j in 1:num_samples_N) {
    s_col_name <- sample_names[j]
    data_col <- my_eigen_data[[s_col_name]]
    
    # Handle diagonal cases (R-squared is 1.0 for self-comparison)
    if (s_row_name == s_col_name) {
      r_squared_all_points_matrix[s_row_name, s_col_name] <- 1.0
      r_squared_highlighted_points_matrix[s_row_name, s_col_name] <- 1.0
    } else {
      # --- Calculate R-squared for ALL points ---
      # Identify valid (non-NA, non-Inf) indices for both samples
      # (Simplified to rep(TRUE, length(data_row)) if you've confirmed no NA/Inf/NaN)
      valid_indices_full <- !is.na(data_row) & !is.na(data_col) & is.finite(data_row) & is.finite(data_col)
      
      if (sum(valid_indices_full) >= min_points_for_lm) {
        lm_all_points <- lm(data_col[valid_indices_full] ~ data_row[valid_indices_full])
        r_squared_all_points_matrix[s_row_name, s_col_name] <- summary(lm_all_points)$r.squared
      } else {
        r_squared_all_points_matrix[s_row_name, s_col_name] <- NA # Not enough valid points
      }
      
      # --- Calculate R-squared for HIGHLIGHTED points ---
      # First, get M-values to identify highlighted points from the valid subset
      M_val_for_highlight_check <- data_col[valid_indices_full] - data_row[valid_indices_full]
      
      # Identify points beyond the cutoff from the *valid* subset
      points_to_highlight_indices_in_valid_subset <- which(M_val_for_highlight_check > cutoff_value | M_val_for_highlight_check < -cutoff_value)
      
      # Extract original data_row and data_col values for the highlighted indices
      highlighted_data_row_values <- data_row[valid_indices_full][points_to_highlight_indices_in_valid_subset]
      highlighted_data_col_values <- data_col[valid_indices_full][points_to_highlight_indices_in_valid_subset]
      
      if (length(highlighted_data_row_values) >= min_points_for_lm) {
        # Ensure there are enough finite points after subsetting
        valid_lm_subset_indices <- is.finite(highlighted_data_row_values) & is.finite(highlighted_data_col_values)
        
        if (sum(valid_lm_subset_indices) >= min_points_for_lm) {
          # Linear regression between the eigenvector values themselves for highlighted points
          lm_highlighted <- lm(highlighted_data_col_values[valid_lm_subset_indices] ~ highlighted_data_row_values[valid_lm_subset_indices])
          r_squared_highlighted_points_matrix[s_row_name, s_col_name] <- summary(lm_highlighted)$r.squared
        } else {
          r_squared_highlighted_points_matrix[s_row_name, s_col_name] <- NA # Not enough valid points
        }
      } else {
        r_squared_highlighted_points_matrix[s_row_name, s_col_name] <- NA # Assign NA if no highlighted points
      }
    }
  }
}
cat("Pairwise R-squared calculations complete.\n\n")

# You can now print these matrices or use them for heatmaps as discussed previously
print("R-squared matrix for all points:")
print(r_squared_all_points_matrix)

print("R-squared matrix for highlighted points:")
print(r_squared_highlighted_points_matrix)


# Step 2: Cluster

library(Heatplus)
r_squared_matrix_for_heatmap <- as.matrix(r_squared_highlighted_points_matrix)# Define a color palette for the heatmap# You can customize this, e.g., using RColorBrewer# A gradient from white/light yellow to dark blue/purple is common for correlation/R^2# Here, we create a simple gradient from white to blue.
heatmap_colors <- colorRampPalette(c("white", "lightblue", "steelblue", "darkblue"))(100)# Generate the clustered heatmap# We'll save it to a PNG file
png("r_squared_highlighted_heatmap_Heatplus_clustered_7.13.2025.png", width = 800, height = 700, res = 100)
heatmap_plus(
  r_squared_matrix_for_heatmap,
  col = heatmap_colors,                           # Color palette
  col.range = c(0, 1),                            # Force color range from 0 to 1 for R^2
  cluster.rows = TRUE,                            # Cluster rows
  cluster.columns = TRUE,                         # Cluster columns
  scale = "none",                                 # Do not scale the values (R^2 is already scaled 0-1)
  legend = 1,                                     # Show color legend
  cexRow = 1.0,                                   # Font size for row labels
  cexCol = 1.0,                                   # Font size for column labels
  lmat = rbind(c(0, 3, 4), c(2, 1, 0)),           # Layout for dendrograms, heatmap, and legend
  lhei = c(1, 4),                                 # Height ratios for dendrograms and heatmap
  lwid = c(1, 4, 1)                               # Width ratios for row dendrogram, heatmap, and legend
)
dev.off() # Close the PNG device
cat("Generated clustered heatmap: r_squared_heatmap_Heatplus_clustered.png\n")
