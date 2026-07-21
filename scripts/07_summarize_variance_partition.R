# Summarize and plot variance-partition results

library(variancePartition)
library(ggplot2)


# Command-line arguments -------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop(
    paste(
      "Usage:",
      "Rscript scripts/07_summarize_variance_partition.R",
      "<results_csv> <analysis_label> <output_prefix>"
    )
  )
}

results_file <- args[1]
analysis_label <- args[2]
output_prefix <- args[3]


# Output directories -----------------------------------------------------

summary_dir <- "results/variance_partition/summaries"
figure_dir <- "results/variance_partition/figures"

dir.create(
  summary_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# Load results -----------------------------------------------------------

message("Loading results: ", results_file)

variance_results <- read.csv(
  results_file,
  row.names = 1,
  check.names = FALSE
)

message(
  "Variance results: ",
  nrow(variance_results),
  " probes x ",
  ncol(variance_results),
  " components"
)


# Validate results -------------------------------------------------------

if (anyNA(variance_results)) {
  stop("Missing values found in variance-partition results.")
}

variance_sums <- rowSums(variance_results)

if (!all(abs(variance_sums - 1) < 1e-6)) {
  stop("Variance fractions do not sum to 1 for every probe.")
}

message("All variance fractions sum to 1.")
message(
  "Components: ",
  paste(colnames(variance_results), collapse = ", ")
)

# Native variancePartition plot ------------------------------------------

variance_plot <- plotVarPart(
  sortCols(variance_results)
) +
  ggtitle(
    paste("Variance Partitioning:", analysis_label)
  ) +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

figure_file <- file.path(
  figure_dir,
  paste0(
    output_prefix,
    "_variance_partition.png"
  )
)

ggsave(
  filename = figure_file,
  plot = variance_plot,
  width = 9,
  height = 6,
  dpi = 300
)

message("Figure saved to: ", figure_file)

# Summarize variance components ------------------------------------------

component_summary <- data.frame(
  component = colnames(variance_results),
  median_variance = apply(
    variance_results,
    2,
    median
  ),
  mean_variance = colMeans(
    variance_results
  ),
  probes_above_50_percent = colSums(
    variance_results > 0.50
  ),
  stringsAsFactors = FALSE
)

component_summary$percent_probes_above_50 <- (
  component_summary$probes_above_50_percent /
    nrow(variance_results)
) * 100

component_summary <- component_summary[
  order(
    component_summary$median_variance,
    decreasing = TRUE
  ),
]

summary_file <- file.path(
  summary_dir,
  paste0(
    output_prefix,
    "_component_summary.csv"
  )
)

write.csv(
  component_summary,
  summary_file,
  row.names = FALSE
)

message("Component summary saved to: ", summary_file)

print(component_summary)
