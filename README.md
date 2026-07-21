# ImmVar variancePartition analysis

This project applies variancePartition package to the ImmVar CD4 T-cell
and CD14 monocyte expression datasets. The goal is to
quantify how much gene-expression variation is associated with cell type,
individual donor, age, sex, batch, and residual variation.

## GEO Data

- GSE56033 - CD4 T cells
- GSE56034 - CD14 monocytes

The downloaded GEO Series Matrix files are not committed to Git. Place them in
`data/raw/` before running the workflow:
data/raw/GSE56033_series_matrix.txt.gz
data/raw/GSE56034_series_matrix.txt.gz


## Project structure

immvar-variancepartition/
├── data/                     # downloaded and generated data (Git-ignored)
├── notebooks/                # metadata, expression, and combined-data QC
├── results/
│   └── variance_partition/   # model outputs, summaries, and figures
├── scripts/                  # numbered preprocessing and analysis scripts
├── environment.yml
└── README.md


## Environment

Create the Python environment and launch JupyterLab:

```bash
conda env create -f environment.yml
conda activate immvar-vp
jupyter lab
```

The statistical analysis also requires R and the following packages:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c("variancePartition", "BiocParallel"))
install.packages("ggplot2")
```

## Workflow

Run all commands from the project root.

### 1. Parse the GEO Series Matrix files

```
python scripts/01_parse_geo.py \
  data/raw/GSE56033_series_matrix.txt.gz cd4

python scripts/01_parse_geo.py \
  data/raw/GSE56034_series_matrix.txt.gz cd14
```

### 2. Clean the sample metadata

```
python scripts/02_clean_metadata.py \
  data/processed/cd4_metadata_raw.csv \
  data/processed/cd4_metadata.csv

python scripts/02_clean_metadata.py \
  data/processed/cd14_metadata_raw.csv \
  data/processed/cd14_metadata.csv
```

### 3. Merge and transform the expression data

```
python scripts/03_merge_datasets.py \
  data/processed/cd4_expression.csv \
  data/processed/cd14_expression.csv \
  data/processed/cd4_metadata.csv \
  data/processed/cd14_metadata.csv \
  data/processed/immvar_expression.csv \
  data/processed/immvar_expression_log2.csv \
  data/processed/immvar_metadata.csv
```

The Series Matrix values include a small number of non-positive measurements.
The merge script removes the 265 probes containing any non-positive value and
then applies a log2 transformation, leaving 21,462 probes for modeling.

### 4. Run variance-partition models

The primary model includes all samples:

```
Rscript scripts/04_run_combined_variance_partition.R
```

The paired-donor sensitivity analysis retains donors represented in both cell
types:

```
Rscript scripts/05_run_paired_variance_partition.R
```

The within-cell-type analyses fit CD4 and CD14 separately:

```
Rscript scripts/06_run_cell_type_variance_partition.R
```

The combined and paired models use:

```r
~ age + (1 | sex) + (1 | batch) +
  (1 | individual_id) + (1 | cell_type)
```

The cell-type-specific models omit `cell_type`


### 5. Generate summaries and figures

The summary script accepts a result CSV, a plot label, and an output prefix:

```
Rscript scripts/07_summarize_variance_partition.R \
  results/variance_partition/immvar_variance_partition.csv \
  "Combined ImmVar Dataset" \
  immvar

Rscript scripts/07_summarize_variance_partition.R \
  results/variance_partition/immvar_paired_variance_partition.csv \
  "Paired ImmVar Donors" \
  immvar_paired

Rscript scripts/07_summarize_variance_partition.R \
  results/variance_partition/cd4_variance_partition.csv \
  "CD4 T Cells" \
  cd4

Rscript scripts/07_summarize_variance_partition.R \
  results/variance_partition/cd14_variance_partition.csv \
  "CD14 Monocytes" \
  cd14
```

## Main results

Cell type was the dominant source of gene-expression variation in the combined
ImmVar dataset. The paired-donor sensitivity analysis produced nearly identical
results. Individual and residual variation were more dominant when CD4 and CD14
cells were analyzed separately. Figures and summaries are available in
results/variance_partition.




## References

- Hoffman GE, Schadt EE. variancePartition: interpreting drivers of variation
  in complex gene expression studies. BMC Bioinformatics. 2016.
- Raj T, et al. Polarization of the effects of autoimmune and neurodegenerative
  risk alleles in leukocytes. Science. 2014.
