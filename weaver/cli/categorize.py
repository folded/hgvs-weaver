"""
Categorization of mapping failures.

"""

import argparse
import csv


def is_p_match(pred: str, truth: str) -> bool:
    """Checks if predicted protein matches truth with fuzzy logic."""
    if pred.startswith("ERR:"):
        return False
    p = pred.replace("(", "").replace(")", "").split(":")[-1]
    t = truth.replace("(", "").replace(")", "").split(":")[-1]
    return p == t or ("Ter" in p and "Ter" in t) or ("fs" in p and "fs" in t) or ("=" in p and "=" in t)


def categorize_rs_failure(row: dict[str, str]) -> str:
    """Categorizes the reason for a weaver mapping failure."""
    rs_p = row["rs_p"]
    rs_spdi = row["rs_spdi"]
    cv_p = row["variant_prot"]
    cv_spdi = row["spdi"]

    if rs_p.startswith("ERR:RefMismatch"):
        return "RefSeq Data Mismatch"
    if rs_p.startswith("ERR:Validate"):
        return "Coordinate Out of Bounds"
    if rs_p.startswith("ERR:Parse"):
        return "Parsing Error"

    p_match = is_p_match(rs_p, cv_p)
    spdi_match = rs_spdi == cv_spdi

    if not p_match and not spdi_match:
        if "Ter" in cv_p and "Ter" not in rs_p:
            return "Protein: Stop Codon Logic"
        if "fs" in cv_p and "fs" not in rs_p:
            return "Protein: Frameshift Logic"
        if "del" in cv_spdi and "ins" in cv_spdi:
            return "Complex Indel"
        return "General Mapping Failure"

    if not p_match:
        if "Ter?" in rs_p:
            return "Protein: Uncertain Stop"
        return "Protein: Translation/Formatting"

    if not spdi_match:
        if "dup" in cv_spdi or "dup" in row["variant_nuc"]:
            return "Genomic: Duplication Shifting"
        return "Genomic: Normalization/Projection"

    return "Unknown"


def main() -> None:
    """Main entry point for categorization."""
    parser = argparse.ArgumentParser(description="Categorize weaver mapping failures.")
    parser.add_argument("input_file", help="Input disagreements TSV file.")
    args = parser.parse_args()

    categories: dict[str, int] = {}
    total = 0

    with open(args.input_file) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            cat = categorize_rs_failure(row)
            categories[cat] = categories.get(cat, 0) + 1
            total += 1

    print(f"Failure Categorization ({total} variants)")
    print("-" * 50)
    for cat, count in sorted(categories.items(), key=lambda x: x[1], reverse=True):
        print(f"{cat:35} | {count:6} ({count / total * 100:4.1f}%)")


if __name__ == "__main__":
    main()
