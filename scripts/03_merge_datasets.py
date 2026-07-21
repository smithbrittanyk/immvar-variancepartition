from pathlib import Path
import argparse
import numpy as np
import pandas as pd


def load_datasets(
    cd4_expression_file: Path,
    cd14_expression_file: Path,
    cd4_metadata_file: Path,
    cd14_metadata_file: Path,
):
    """Load the processed expression and metadata tables."""
    cd4_expression = pd.read_csv(cd4_expression_file, index_col=0)
    cd14_expression = pd.read_csv(cd14_expression_file, index_col=0)

    cd4_metadata = pd.read_csv(cd4_metadata_file)
    cd14_metadata = pd.read_csv(cd14_metadata_file)

    return (
        cd4_expression,
        cd14_expression,
        cd4_metadata,
        cd14_metadata,
    )


def validate_probes(
    cd4_expression: pd.DataFrame,
    cd14_expression: pd.DataFrame,
) -> None:
    """Confirm that both expression matrices contain the same probes."""
    same_probes = cd4_expression.index.equals(cd14_expression.index)

    if not same_probes:
        raise ValueError(
            "CD4 and CD14 probe identifiers or their order do not match."
        )

    print("Probe identifiers and order match.")

def merge_datasets(
    cd4_expression: pd.DataFrame,
    cd14_expression: pd.DataFrame,
    cd4_metadata: pd.DataFrame,
    cd14_metadata: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Combine CD4 and CD14 expression and metadata tables."""
    if not cd4_metadata.columns.equals(cd14_metadata.columns):
        raise ValueError("CD4 and CD14 metadata columns do not match.")

    combined_expression = pd.concat(
        [cd4_expression, cd14_expression],
        axis=1,
    )

    combined_metadata = pd.concat(
        [cd4_metadata, cd14_metadata],
        axis=0,
        ignore_index=True,
    )

    return combined_expression, combined_metadata


def prepare_log2_expression(
    expression: pd.DataFrame,
) -> pd.DataFrame:
    """Remove probes with non-positive values and log2-transform expression."""
    positive_probes = (expression > 0).all(axis=1)
    removed_probes = (~positive_probes).sum()

    log2_expression = np.log2(expression.loc[positive_probes])

    print(f"Removed {removed_probes} probes with non-positive values.")
    print("Log2 expression:", log2_expression.shape)

    return log2_expression


def validate_samples(
    expression: pd.DataFrame,
    metadata: pd.DataFrame,
) -> None:
    """Confirm that expression and metadata samples match."""
    if expression.columns.duplicated().any():
        raise ValueError("Duplicate sample IDs found in expression data.")

    if metadata["geo_accession"].duplicated().any():
        raise ValueError("Duplicate GEO accessions found in metadata.")

    expression_samples = expression.columns
    metadata_samples = pd.Index(metadata["geo_accession"])

    if set(expression_samples) != set(metadata_samples):
        raise ValueError(
            "Expression and metadata contain different sample IDs."
        )

    if not expression_samples.equals(metadata_samples):
        raise ValueError(
            "Expression and metadata samples are not in the same order."
        )

    print("Expression and metadata samples match.")


def save_datasets(
    expression: pd.DataFrame,
    log2_expression: pd.DataFrame,
    metadata: pd.DataFrame,
    expression_file: Path,
    log2_expression_file: Path,
    metadata_file: Path,
) -> None:
    """Save the combined expression and metadata tables."""
    expression_file.parent.mkdir(parents=True, exist_ok=True)
    log2_expression_file.parent.mkdir(parents=True, exist_ok=True)
    metadata_file.parent.mkdir(parents=True, exist_ok=True)

    expression.to_csv(expression_file)
    log2_expression.to_csv(log2_expression_file)
    metadata.to_csv(metadata_file, index=False)

    print(f"Wrote {expression_file}")
    print(f"Wrote {log2_expression_file}")
    print(f"Wrote {metadata_file}")

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Merge CD4 and CD14 expression and metadata tables."
    )
    parser.add_argument(
        "cd4_expression",
        type=Path,
        help="Path to the processed CD4 expression CSV.",
    )
    parser.add_argument(
        "cd14_expression",
        type=Path,
        help="Path to the processed CD14 expression CSV.",
    )
    parser.add_argument(
        "cd4_metadata",
        type=Path,
        help="Path to the cleaned CD4 metadata CSV.",
    )
    parser.add_argument(
        "cd14_metadata",
        type=Path,
        help="Path to the cleaned CD14 metadata CSV.",
    )
    parser.add_argument(
        "output_expression",
        type=Path,
        help="Path for the combined expression CSV.",
    )

    parser.add_argument(
        "output_log2_expression",
        type=Path,
        help="Path for the filtered log2-transformed expression CSV.",
    )

    parser.add_argument(
        "output_metadata",
        type=Path,
        help="Path for the combined metadata CSV.",
    )

    args = parser.parse_args()

    (
        cd4_expression,
        cd14_expression,
        cd4_metadata,
        cd14_metadata,
    ) = load_datasets(
        args.cd4_expression,
        args.cd14_expression,
        args.cd4_metadata,
        args.cd14_metadata,
    )

    validate_probes(cd4_expression, cd14_expression)

    print("CD4 expression:", cd4_expression.shape)
    print("CD14 expression:", cd14_expression.shape)
    print("CD4 metadata:", cd4_metadata.shape)
    print("CD14 metadata:", cd14_metadata.shape)

    combined_expression, combined_metadata = merge_datasets(
        cd4_expression,
        cd14_expression,
        cd4_metadata,
        cd14_metadata,
    )

    print("Combined expression:", combined_expression.shape)
    print("Combined metadata:", combined_metadata.shape)

    validate_samples(combined_expression, combined_metadata)

    log2_expression = prepare_log2_expression(combined_expression)

    save_datasets(
        combined_expression,
        log2_expression,
        combined_metadata,
        args.output_expression,
        args.output_log2_expression,
        args.output_metadata,
    )

if __name__ == "__main__":
    main()
