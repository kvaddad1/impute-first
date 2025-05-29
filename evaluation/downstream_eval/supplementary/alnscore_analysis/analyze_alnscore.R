# analyze_alnscore.R
# Usage: Rscript analyze_alnscore.R path_to_input_file
#
# This script analyzes alignment score differences.
# The input file should contain these score differences in a single column.
# The column with score differences is expected to be named 'score_diff'.

args <- commandArgs(trailingOnly = TRUE)

# Check if the file path is provided
if (length(args) == 0) {
    stop("No input file provided. Usage: Rscript analyze_alnscore.R path_to_input_file", call. = FALSE)
}

input_file <- args[1]

# Read the data
data <- read.csv(input_file, header = TRUE, sep = ",")

# Extract positive and negative values
positive_values <- data$score_diff[data$score_diff > 0]
negative_values <- data$score_diff[data$score_diff < 0]

# Calculate statistics for positive values
positive_stats <- c(
  Count = length(positive_values),
  Density = length(positive_values) / length(data$score_diff),
  Mean = mean(positive_values),
  Median = median(positive_values),
  Variance = var(positive_values),
  Std_Deviation = sd(positive_values),
  Min = min(positive_values),
  Max = max(positive_values),
  Q1 = quantile(positive_values, 0.25),
  Q3 = quantile(positive_values, 0.75)
)

# Calculate statistics for negative values
negative_stats <- c(
  Count = length(negative_values),
  Density = length(negative_values) / length(data$score_diff),
  Mean = mean(negative_values),
  Median = median(negative_values),
  Variance = var(negative_values),
  Std_Deviation = sd(negative_values),
  Min = min(negative_values),
  Max = max(negative_values),
  Q1 = quantile(negative_values, 0.25),
  Q3 = quantile(negative_values, 0.75)
)

# Print the statistics
cat("Statistics for Positive Differences:\n")
print(positive_stats)

cat("\nStatistics for Negative Differences:\n")
print(negative_stats)
