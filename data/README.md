# Data directory

Downloaded and generated datasets are excluded from Git because the expression
matrices are large and can be recreated with the project scripts.

Download the GEO Series Matrix files for GSE56033 and GSE56034 and place them
in `data/raw/`:

```text
data/raw/GSE56033_series_matrix.txt.gz
data/raw/GSE56034_series_matrix.txt.gz
```

Running scripts `01` through `03` creates the analysis-ready files under
`data/processed/`.
