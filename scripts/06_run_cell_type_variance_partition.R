# Run variance partition analysis for ImmVar expression data

library(variancePartition)
library(BiocParallel)


# -------------------------------------------------------------------------
# File paths
# -------------------------------------------------------------------------

data_dir <- "data/processed"
results_dir <- "results/variance_partition"

dir.create(
  results_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# -------------------------------------------------------------------------
# Function to run variance partition for one cell type
# -------------------------------------------------------------------------

run_variance_partition <- function(expression_file,
                                   metadata_file,
                                   target_cell_type,
                                   output_prefix) {

  message("Loading expression data: ", expression_file)
  message("Loading metadata: ", metadata_file)

  # Load expression matrix
  expression <- read.csv(
    expression_file,
    row.names = 1,
    check.names = FALSE
  )

  # Load metadata
  metadata <- read.csv(
    metadata_file,
    stringsAsFactors = FALSE
  )

  # Subset metadata and expression to the requested cell type
  metadata <- metadata[
    metadata$cell_type == target_cell_type,
    ,
    drop = FALSE
  ]

  sample_ids <- metadata$geo_accession

  if (length(sample_ids) == 0) {
    stop("No samples found for cell type: ", target_cell_type)
  }

  expression <- expression[
    ,
    colnames(expression) %in% sample_ids,
    drop = FALSE
  ]

  message(
    target_cell_type,
    ": ",
    nrow(expression),
    " probes x ",
    ncol(expression),
    " samples"
  )


  # Convert expression data to a numeric matrix
  expression <- as.matrix(expression)
  storage.mode(expression) <- "numeric"


  # -----------------------------------------------------------------------
  # Validate sample identifiers
  # -----------------------------------------------------------------------

  if (anyDuplicated(colnames(expression))) {
    stop("Duplicate sample IDs found in the expression matrix.")
  }

  if (anyDuplicated(metadata$geo_accession)) {
    stop("Duplicate GEO accession IDs found in the metadata.")
  }


  # -----------------------------------------------------------------------
  # Align metadata rows with expression columns
  # -----------------------------------------------------------------------

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



  # -----------------------------------------------------------------------
  # Prepare model variables
  # -----------------------------------------------------------------------

  metadata$age <- as.numeric(metadata$age)
  metadata$sex <- factor(metadata$sex)
  metadata$batch <- factor(metadata$batch)
  metadata$individual_id <- factor(metadata$individual_id)

  model_variables <- c(
    "age",
    "sex",
    "batch",
    "individual_id"
  )

  if (anyNA(metadata[, model_variables])) {
    stop("Missing values found in one or more model variables.")
  }


  # -----------------------------------------------------------------------
  # Define and fit variance partition model
  # -----------------------------------------------------------------------

  form <- ~ age + (1 | sex) + (1 | batch) + (1 | individual_id)

  message("Fitting model: ")
  print(form)

  variance_results <- fitExtractVarPartModel(
    expression,
    form,
    metadata,
    BPPARAM = SerialParam()
  )


  # -----------------------------------------------------------------------
  # Save results
  # -----------------------------------------------------------------------

  output_file <- paste0(
    output_prefix,
    "_variance_partition.csv"
  )

  write.csv(
    variance_results,
    output_file,
    row.names = TRUE
  )

  message("Results saved to: ", output_file)

  return(variance_results)
}


# -------------------------------------------------------------------------
# Run CD4 analysis
# -------------------------------------------------------------------------

cd4_results <- run_variance_partition(
  expression_file = file.path(
    data_dir,
    "immvar_expression_log2.csv"
  ),
  metadata_file = file.path(
    data_dir,
    "immvar_metadata.csv"
  ),
  target_cell_type = "CD4",
  output_prefix = file.path(
    results_dir,
    "cd4"
  )
)


# -------------------------------------------------------------------------
# Run CD14 analysis
# -------------------------------------------------------------------------

cd14_results <- run_variance_partition(
  expression_file = file.path(
    data_dir,
    "immvar_expression_log2.csv"
  ),
  metadata_file = file.path(
    data_dir,
    "immvar_metadata.csv"
  ),
  target_cell_type = "CD14",
  output_prefix = file.path(
    results_dir,
    "cd14"
  )
)


# -------------------------------------------------------------------------
# Session information
# -------------------------------------------------------------------------

sessionInfo()
