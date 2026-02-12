"""Analysis of HGVS validation results."""

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


def main() -> None:
    """Main analysis entry point."""
    parser = argparse.ArgumentParser(description="Analyze full HGVS validation results.")
    parser.add_argument("input_file", help="Input validation TSV file.")
    args = parser.parse_args()

    total = 0
    rs_p_match = 0
    ref_p_match = 0
    rs_spdi_match = 0
    ref_spdi_match = 0

    rs_parse_err = 0
    ref_parse_err = 0
    rs_ref_mismatch = 0

    with open(args.input_file) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            total += 1

            # ClinVar truth
            cv_p_full = row["variant_prot"]
            cv_p = clean_hgvs(cv_p_full)
            cv_spdi = row["spdi"]

            # rs
            rs_p_raw = row["rs_p"]
            rs_p = clean_hgvs(rs_p_raw)
            rs_spdi = row["rs_spdi"]

            # ref
            ref_p_raw = row["ref_p"]
            ref_p = clean_hgvs(ref_p_raw)
            ref_spdi = row["ref_spdi"]

            # Analysis
            if rs_p_raw.startswith("ERR:RefMismatch"):
                rs_ref_mismatch += 1
            if rs_p_raw.startswith("ERR:Parse"):
                rs_parse_err += 1
            if ref_p_raw.startswith("ERR:Parse"):
                ref_parse_err += 1

            # Normalize stop codon for fuzzy matching
            cv_p_norm = cv_p.replace("Ter", "*")
            rs_p_norm = rs_p.replace("Ter", "*")
            ref_p_norm = ref_p.replace("Ter", "*")

            # Protein fuzzy matches
            if (
                rs_p_norm == cv_p_norm
                or ("*" in rs_p_norm and "*" in cv_p_norm)
                or ("fs" in rs_p_norm and "fs" in cv_p_norm)
                or ("=" in rs_p_norm and "=" in cv_p_norm)
            ):
                rs_p_match += 1

            if (
                ref_p_norm == cv_p_norm
                or ("*" in ref_p_norm and "*" in cv_p_norm)
                or ("fs" in ref_p_norm and "fs" in cv_p_norm)
                or ("=" in ref_p_norm and "=" in cv_p_norm)
            ):
                ref_p_match += 1

            # SPDI exact matches
            if rs_spdi == cv_spdi:
                rs_spdi_match += 1
            if ref_spdi == cv_spdi:
                ref_spdi_match += 1

    if total == 0:
        print("No variants processed.")
        return

    print(f"Summary of {total} Variants")
    print("-" * 40)
    print(f"RefSeq Data Mismatches (rs): {rs_ref_mismatch} ({rs_ref_mismatch / total * 100:.1f}%)")
    print("-" * 40)
    print("Implementation | Protein Match | SPDI Match | Parse Errors")
    print(
        f"weaver         | {rs_p_match / total * 100:12.3f}% | {rs_spdi_match / total * 100:9.3f}% | {rs_parse_err:12}",
    )
    print(
        f"ref-hgvs       | {ref_p_match / total * 100:12.3f}% | {ref_spdi_match / total * 100:9.3f}% | {ref_parse_err:12}",
    )

    print()
    print("-" * 40)
    print(f"RefSeq Data Mismatches (rs): {rs_ref_mismatch} ({rs_ref_mismatch / total * 100:.1f}%)")
    print("-" * 40)
    print("Implementation | Protein Match | SPDI Match | Parse Errors")
    print(
        f"weaver         | {rs_p_match:12d} | {rs_spdi_match:11d} | {rs_parse_err:12d}",
    )
    print(
        f"ref-hgvs       | {ref_p_match:12d} | {ref_spdi_match:11d} | {ref_parse_err:12d}",
    )


if __name__ == "__main__":
    main()
