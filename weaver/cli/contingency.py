"""Generation of contingency tables for implementation agreement."""

import argparse
import csv


def clean_hgvs(s):
    if not s:
        return ""
    # Remove accession prefix
    if ":" in s:
        s = s.split(":")[-1]
    # Remove parentheses
    s = s.replace("(", "").replace(")", "")
    # Standardize Ter/*
    s = s.replace("Ter", "*")
    return s


def is_p_match(pred: str, truth: str) -> bool:
    """Checks if predicted protein matches truth with fuzzy logic."""
    if pred.startswith("ERR:"):
        return False
    p = clean_hgvs(pred)
    t = clean_hgvs(truth)
    return p == t or ("*" in p and "*" in t) or ("fs" in p and "fs" in t) or ("=" in p and "=" in t)


def print_table(title: str, stats: dict[str, int], _total: int) -> None:
    """Prints a formatted contingency table."""
    print(f"\n{title}")
    print("-" * 60)
    print(f"{'':20} | {'ref Match':15} | {'ref Mismatch':15}")
    print("-" * 60)
    print(f"{'weaver Match':20} | {stats['both']:15,d} | {stats['rs_only']:15,d}")
    print(f"{'weaver Mismatch':20} | {stats['ref_only']:15,d} | {stats['neither']:15,d}")
    print("-" * 60)


def main() -> None:
    """Main contingency entry point."""
    parser = argparse.ArgumentParser(description="Generate contingency tables for implementation agreement.")
    parser.add_argument("input_file", help="Input full validation TSV.")
    args = parser.parse_args()

    p_stats = {"both": 0, "rs_only": 0, "ref_only": 0, "neither": 0}
    spdi_stats = {"both": 0, "rs_only": 0, "ref_only": 0, "neither": 0}
    total = 0

    with open(args.input_file) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            total += 1
            cv_p = row["variant_prot"]
            cv_spdi = row["spdi"]

            rs_p_ok = is_p_match(row["rs_p"], cv_p)
            ref_p_ok = is_p_match(row["ref_p"], cv_p)

            if rs_p_ok and ref_p_ok:
                p_stats["both"] += 1
            elif rs_p_ok:
                p_stats["rs_only"] += 1
            elif ref_p_ok:
                p_stats["ref_only"] += 1
            else:
                p_stats["neither"] += 1

            rs_spdi_ok = row["rs_spdi"] == cv_spdi
            ref_spdi_ok = row["ref_spdi"] == cv_spdi

            if rs_spdi_ok and ref_spdi_ok:
                spdi_stats["both"] += 1
            elif rs_spdi_ok:
                spdi_stats["rs_only"] += 1
            elif ref_spdi_ok:
                spdi_stats["ref_only"] += 1
            else:
                spdi_stats["neither"] += 1

    print_table("Protein Translation Agreement", p_stats, total)
    print_table("SPDI Mapping Agreement", spdi_stats, total)


if __name__ == "__main__":
    main()
