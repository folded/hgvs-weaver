"""
Extraction of mapping disagreements."""

import argparse
import csv


def is_protein_match(p_raw: str, cv_p_raw: str) -> bool:
    """Checks if predicted protein matches truth with fuzzy logic."""
    if p_raw.startswith("ERR:"):
        return False
    # Simplified fuzzy match for frameshifts and stop codons
    p = p_raw.replace("(", "").replace(")", "").split(":")[-1]
    cv = cv_p_raw.replace("(", "").replace(")", "").split(":")[-1]
    return p == cv or ("Ter" in p and "Ter" in cv) or ("fs" in p and "fs" in cv) or ("=" in p and "=" in cv)


def main() -> None:
    """Main extraction entry point."""
    parser = argparse.ArgumentParser(description="Extract mapping disagreements between implementations and truth.")
    parser.add_argument("input_file", help="Input full validation TSV.")
    parser.add_argument("--output-file", default="disagreements.tsv", help="Output TSV file.")
    args = parser.parse_args()

    disagreements = []
    fieldnames = []
    with open(args.input_file) as f:
        reader = csv.DictReader(f, delimiter="\t")
        fieldnames = list(reader.fieldnames or [])
        for row in reader:
            cv_p = row["variant_prot"]
            cv_spdi = row["spdi"]

            rs_p_match = is_protein_match(row["rs_p"], cv_p)
            rs_spdi_match = row["rs_spdi"] == cv_spdi

            if not rs_p_match or not rs_spdi_match:
                disagreements.append(row)

    if fieldnames:
        with open(args.output_file, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(disagreements)

    print(f"Extracted {len(disagreements)} disagreements to {args.output_file}")


if __name__ == "__main__":
    main()
