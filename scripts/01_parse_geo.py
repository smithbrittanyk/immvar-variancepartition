from io import StringIO
from pathlib import Path
import argparse
import csv
import gzip

import pandas as pd


def read_series_matrix(input_file: Path) -> list[str]:
    """Read a gzipped GEO Series Matrix file."""
    with gzip.open(input_file, "rt", encoding="utf-8") as handle:
        return handle.readlines()


def parse_expression(lines: list[str]) -> pd.DataFrame:
    """Extract the expression matrix from a GEO Series Matrix file."""
    table_start = lines.index("!series_matrix_table_begin\n")
    table_end = lines.index("!series_matrix_table_end\n")

    expression_text = "".join(lines[table_start + 1 : table_end])

    expression = pd.read_csv(
        StringIO(expression_text),
        sep="\t",
        index_col=0,
    )

    expression.index.name = "ID_REF"
    return expression


def parse_metadata(lines: list[str]) -> pd.DataFrame:
    """Extract sample-level metadata rows."""
    metadata_rows = []

    for line in lines:
        if line.startswith("!Sample_"):
            fields = next(csv.reader([line], delimiter="\t"))
            metadata_rows.append(fields)

    metadata = pd.DataFrame(metadata_rows).set_index(0).T
    metadata.index.name = "sample_number"

    return metadata


def parse_geo(input_file: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Parse expression and metadata from one GEO Series Matrix file."""
    lines = read_series_matrix(input_file)

    expression = parse_expression(lines)
    metadata = parse_metadata(lines)

    return expression, metadata


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Parse a GEO Series Matrix file into CSV files."
    )
    parser.add_argument(
        "input_file",
        type=Path,
        help="Path to the .txt.gz GEO Series Matrix file.",
    )
    parser.add_argument(
        "sample_type",
        choices=["cd4", "cd14"],
        help="Sample type used to name the output files.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("data/processed"),
        help="Directory for parsed CSV files.",
    )

    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    expression, metadata = parse_geo(args.input_file)


    """ Parsing Checks """
    assert expression.shape[1] == metadata.shape[0], (
        "Expression and metadata sample counts do not match."
    )

    assert expression.index.is_unique, (
        "Expression feature IDs are not unique."
    )

    assert not expression.columns.duplicated().any(), (
        "Expression sample names are duplicated."
    )


    """ Output processed data """ 
    expression_file = args.output_dir / f"{args.sample_type}_expression.csv"
    metadata_file = args.output_dir / f"{args.sample_type}_metadata_raw.csv"

    expression.to_csv(expression_file)
    metadata.to_csv(metadata_file)

    print(
        f"{args.sample_type.upper()} expression: "
        f"{expression.shape[0]:,} features x {expression.shape[1]:,} samples"
    )
    print(
        f"{args.sample_type.upper()} metadata: "
        f"{metadata.shape[0]:,} samples x {metadata.shape[1]:,} fields"
    )
    print(f"Wrote {expression_file}")
    print(f"Wrote {metadata_file}")

    print("\nFirst metadata fields:")
    print(metadata.columns[:10].tolist())
    

if __name__ == "__main__":
    main()
