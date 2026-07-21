# Run variance partition analysis on the combined ImmVar dataset

library(variancePartition)
library(BiocParallel)

# File paths -------------------------------------------------------------

expression_file <- "data/processed/immvar_expression_log2.csv"
metadata_file <- "data/processed/immvar_metadata.csv"

results_dir <- "results/variance_partition"

dir.create(
  results_dir,
  recursive = TRUE,
  showWarnings = FALSE
)



# Load data --------------------------------------------------------------

message("Loading expression data: ", expression_file)

expression <- read.csv(
  expression_file,
  row.names = 1,
  check.names = FALSE
)

message("Loading metadata: ", metadata_file)

metadata <- read.csv(
  metadata_file,
  stringsAsFactors = FALSE
)

expression <- as.matrix(expression)
storage.mode(expression) <- "numeric"

message(
  "Expression matrix: ",
  nrow(expression),
  " probes x ",
  ncol(expression),
  " samples"
)

message(
  "Metadata: ",
  nrow(metadata),
  " samples x ",
  ncol(metadata),
  " fields"
)

# Validate and align samples ---------------------------------------------

if (anyDuplicated(colnames(expression))) {
  stop("Duplicate sample IDs found in the expression matrix.")
}

if (anyDuplicated(metadata$geo_accession)) {
  stop("Duplicate GEO accession IDs found in the metadata.")
}

sample_order <- match(
  colnames(expression),
  metadata$geo_accession
)

if (any(is.na(sample_order))) {
  missing_samples <- colnames(expression)[is.na(sample_order)]

  stop(
    "Expression samples missing from metadata: ",
    paste(missing_samples, collapse = ", ")
  )
}

metadata <- metadata[sample_order, , drop = FALSE]
rownames(metadata) <- metadata$geo_accession

stopifnot(
  identical(colnames(expression), rownames(metadata))
)

message("Expression and metadata samples match.")

# Prepare model variables ------------------------------------------------

metadata$age <- as.numeric(metadata$age)
metadata$sex <- factor(metadata$sex)
metadata$batch <- factor(metadata$batch)
metadata$individual_id <- factor(metadata$individual_id)
metadata$cell_type <- factor(metadata$cell_type)

model_variables <- c(
  "age",
  "sex",
  "batch",
  "individual_id",
  "cell_type"
)

if (anyNA(metadata[, model_variables])) {
  stop("Missing values found in one or more model variables.")
}

message("Model variables are complete.")

message(
  "Cell types: ",
  paste(levels(metadata$cell_type), collapse = ", ")
)

message(
  "Batches: ",
  paste(levels(metadata$batch), collapse = ", ")
)

message(
  "Individuals: ",
  nlevels(metadata$individual_id)
)

# Fit combined variance-partition model ----------------------------------

form <- ~ age +
  (1 | sex) +
  (1 | batch) +
  (1 | individual_id) +
  (1 | cell_type)

message("Fitting model:")
print(form)

variance_results <- fitExtractVarPartModel(
  expression,
  form,
  metadata,
  BPPARAM = SerialParam()
)


# Save results ------------------------------------------------------------

output_file <- file.path(
  results_dir,
  "immvar_variance_partition.csv"
)

write.csv(
  variance_results,
  output_file,
  row.names = TRUE
)

message("Results saved to: ", output_file)

