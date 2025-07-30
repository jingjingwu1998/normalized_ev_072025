# computing cutoff for MA plot and histogram
# SD adjustment for plotting please visit line 71, 86, 193, 208, 219
setwd('/Users/80030577/Desktop/HiC_analysis/EV_normalization')
getwd()
load("normalized_ev.100k_multi_samples_06_03_2025.rda") 
# List of all NORMALIZED EV objects for easy iteration
normalized_ev_lists <- list(
  RBL = EV.rbl,
  LCL = EV.lcl,
  GCBC = EV.gcbc,
  MBC = EV.mbc,
  NBC = EV.nbc,
  PC = EV.pc
)

# Define sample names based on the normalized data list
sample_names <- names(normalized_ev_lists) 

# Define standard human chromosome names in correct order.
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

my_eigen_data <- concatenated_normalized_evs
sample_lengths <- sapply(my_eigen_data, length) # 'sapply' applies 'length' to each list item.
sample_lengths

#######################################
############ start to edit ############
#######################################

# Part 1: Compute the Global Cutoff Value  #
all_pairwise_sd_of_differences <- c()
sample_pairs_to_compare <- combn(sample_names, 2, simplify = FALSE)


# Loop through each unique pair of samples
for (pair in sample_pairs_to_compare) {
  # Get the names of the two samples in the current pair
  sample1_name <- pair[1]
  sample2_name <- pair[2]
  
  # Extract the eigenvector values for these two samples.
  # We can now assume they have the same length due to the pre-processing step.
  ev_values_sample1 <- my_eigen_data[[sample1_name]]
  ev_values_sample2 <- my_eigen_data[[sample2_name]]
  
  # Calculate the 'Difference' (A value) for each data point
  # A = Value_from_Sample1 - Value_from_Sample2 (without log transformation)
  A_values <- ev_values_sample1 - ev_values_sample2
  
  # Calculate the standard deviation (SD) of these 'A' (Difference) values.
  # 'na.rm = TRUE' tells R to ignore any 'NA' (Not Available) values when calculating SD.
  sd_of_current_A <- sd(A_values, na.rm = TRUE)
  
  # Add this calculated standard deviation to our list
  all_pairwise_sd_of_differences <- c(all_pairwise_sd_of_differences, sd_of_current_A)
}

cutoff_value <- mean(all_pairwise_sd_of_differences, na.rm = TRUE)
cutoff_value <- 2*cutoff_value
cutoff_value # SD = 0.009428408, 2 SD = 0.01885682, 3 SD = 0.02828522


cat("-------------------------------------------------------------------\n")
cat("Calculated Cutoff Value for MA plots:", cutoff_value, "\n")
cat("This cutoff is the mean of the standard deviations of 'Difference (A)' values across all 15 pairwise comparisons.\n")
cat("-------------------------------------------------------------------\n\n")
num_samples_N <- length(sample_names)

# Part 2: The Cutoff Is The Same For All Pairwise Comparisons #


# Part 3: Draw MA Plots in an N x N Matrix to PDF & Collect Highlighted Counts # 

pdf_file_name_matrix <- 'ev.ma_cutoff_at_2_SD_pairwise_plot.pdf'
x_lim_matrix <- c(-0.1, 0.1)  # X-axis (Mean EV) range
y_lim_matrix <- c(-0.15, 0.15) # Y-axis (EV Difference) range
num_samples_N <- length(sample_names)
plot_rows <- num_samples_N
plot_cols <- num_samples_N

# Initialize a data frame to store highlighted counts for the summary plot
highlight_summary_df <- data.frame(
  Pair = character(),
  Highlighted_Count = numeric(),
  Total_Points = numeric(),
  Percentage_Highlighted = numeric(),
  stringsAsFactors = FALSE
)

# Open PDF device for plotting the N x N matrix
pdf(file = pdf_file_name_matrix,
    width = plot_cols * 5, # Calculate total width based on number of columns
    height = plot_rows * 5) # Calculate total height based on number of rows

# Set general plot parameters for the entire matrix
par(mfrow = c(plot_rows, plot_cols),
    font.lab = 2, cex.lab = 1.0,
    mar = c(3, 3, 2, 1),
    oma = c(2, 2, 2, 2),
    mgp = c(1.5, 0.5, 0),
    xaxs = 'i', yaxs = 'i'
)


# Loop through all combinations for the N x N matrix plot
for (i in 1:num_samples_N) {
  s_row_name <- sample_names[i]
  data_row <- my_eigen_data[[s_row_name]]
  
  for (j in 1:num_samples_N) {
    s_col_name <- sample_names[j]
    data_col <- my_eigen_data[[s_col_name]]
    
    # Check if samples are the same (diagonal plots)
    if (s_row_name == s_col_name) {
      # For diagonal plots, show only the sample name
      plot(NA, xlim = x_lim_matrix, ylim = y_lim_matrix, type = "n", xlab = "", ylab = "",
           main = "", cex.main = 1.5, font.main = 2)
      text(mean(x_lim_matrix), mean(y_lim_matrix), labels = s_row_name, cex = 3, col = "gray")
      
    } else {
      # For off-diagonal plots, create the MA plot
      
      # Calculate Mean (X-axis) and Difference (Y-axis) values
      A_val <- (data_row + data_col) / 2 # Mean EV
      M_val <- data_col - data_row     # EV Difference (col - row)
      
      # Identify points beyond the cutoff and count them
      points_to_highlight_indices <- which(M_val > cutoff_value | M_val < -cutoff_value)
      num_highlighted_points <- length(points_to_highlight_indices)
      current_total_points_in_plot <- length(A_val) # Total points in this specific MA plot
      percentage_highlighted <- (num_highlighted_points / current_total_points_in_plot) * 100
      
      # Construct the plot title with the count
      plot_title <- paste0(s_col_name, " vs ", s_row_name, "\nHighlighted: ", num_highlighted_points)
      
      # === FIX: Ensure Pair name is consistently ordered for summary data ===
      # This prevents duplicate entries like "RBL vs LCL" and "LCL vs RBL" in the summary.
      # We sort the two sample names alphabetically before pasting them.
      unique_pair_name <- paste(sort(c(s_row_name, s_col_name)), collapse=" vs ")
      
      # Add to summary data frame ONLY if this unique pair hasn't been added yet.
      # The 'if (i < j)' condition ensures we process each unique comparison only once
      # when adding to the summary data frame, regardless of how the plot title is formatted.
      if (i < j) {
        highlight_summary_df <- rbind(highlight_summary_df, data.frame(
          Pair = unique_pair_name, # Use the consistently ordered name
          Highlighted_Count = num_highlighted_points,
          Total_Points = current_total_points_in_plot,
          Percentage_Highlighted = percentage_highlighted
        ))
      }
      
      # Plot the MA plot
      plot(A_val, M_val,
           pch = 16, cex = 0.25, # Point style
           col = rgb(0, 0, 0, alpha = 0.1), # Semi-transparent black for density
           main = plot_title, # Title with highlighted count
           xlab = "Mean EV", ylab = "EV Difference", # Axis labels
           xlim = x_lim_matrix, # Fixed X-axis limits
           ylim = y_lim_matrix, # Fixed Y-axis limits
           cex.main = 1.0, cex.lab = 0.9, cex.axis = 0.8) # Font sizes
      
      # Add cutoff dash lines (red)
      abline(h = cutoff_value, col = "red", lty = 2, lwd = 2)
      abline(h = -cutoff_value, col = "red", lty = 2, lwd = 2)
      
      # Add a horizontal line at EV Difference = 0
      abline(h = 0, lty = "dashed", col = "blue")
      
      # Re-plot highlighted points in red
      points(A_val[points_to_highlight_indices], M_val[points_to_highlight_indices],
             col = "red", pch = 16, cex = 0.25)
    }
  }
}



# Add an overall title to the entire matrix plot
mtext("Pairwise MA Plots for Normalized EV Data with Cutoff at 2 * SD", side = 3, line = 0, outer = TRUE, cex = 2, font = 2)

# Close the PDF device for the N x N matrix
dev.off()
cat(paste0("\nFull pairwise MA plot matrix generated and saved to '", pdf_file_name_matrix, "'.\n"))

###### to work on ######
#
#
#
#
#
#

### Part 4: Plot Summary of Highlighted Points (Single Bar Per Unique Pair)
pdf_file_name_summary <- 'highlight_summary_plot_with_cutoff_at_2_SD.pdf'
pdf(file = pdf_file_name_summary, width = 12, height = 8)
par(mar = c(8, 5, 4, 2) + 0.1) # Adjust margins: bottom, left, top, right

# === FIX: Ensure uniqueness for the bar plot using duplicated() ===
unique_summary_df <- highlight_summary_df[!duplicated(highlight_summary_df$Highlighted_Count), ]

# Create a bar plot of the percentage of highlighted points
bp <- barplot(unique_summary_df$Percentage_Highlighted,
              names.arg = unique_summary_df$Pair, # Use unique pair names as labels
              ylab = "Percentage of Highlighted Points (%)", # Y-axis label
              main = "Percentage of Highlighted Points per Unique Pairwise MA Plot (2x SD Cutoff)", # Main title
              col = "skyblue", # Bar color
              las = 2, # Rotate X-axis labels vertically for readability
              cex.names = 0.7, # Adjust size of X-axis names
              ylim = c(0, max(unique_summary_df$Percentage_Highlighted) * 1.1), # Extend y-axis slightly above max
              cex.axis = 0.8 # Adjust size of axis numbers
)

# Add text labels on top of each bar showing both count and percentage
text(x = bp,
     y = unique_summary_df$Percentage_Highlighted,
     labels = paste0(unique_summary_df$Highlighted_Count, "\n(", round(unique_summary_df$Percentage_Highlighted, 2), "%)"),
     pos = 3, # Position labels above the bar
     cex = 0.7, # Size of the text labels
     col = "black")

# Close the PDF device for the summary plot
dev.off()
cat(paste0("\nSummary plot of highlighted points saved to '", pdf_file_name_summary, "'.\n"))