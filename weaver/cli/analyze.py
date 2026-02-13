"""Analysis of HGVS validation results with integrated contingency reporting."""

import argparse
import csv
import re
from pathlib import Path


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
    if not pred or pred.startswith("ERR:"):
        return False
    p = clean_hgvs(pred)
    t = clean_hgvs(truth)
    return p == t or ("*" in p and "*" in t) or ("fs" in p and "fs" in t) or ("=" in p and "=" in t)


def main() -> None:
    """Main analysis entry point."""
    parser = argparse.ArgumentParser(description="Analyze full HGVS validation results.")
    parser.add_argument("input_file", help="Input validation TSV file.")
    parser.add_argument(
        "--update-readme",
        action="store_true",
        help="Update the project README.md with the latest results.",
    )
    args = parser.parse_args()

    total = 0
    rs_p_match = 0
    ref_p_match = 0
    rs_spdi_match = 0
    ref_spdi_match = 0

    rs_parse_err = 0
    ref_parse_err = 0
    rs_ref_mismatch = 0

    p_stats = {"both": 0, "rs_only": 0, "ref_only": 0, "neither": 0}
    spdi_stats = {"both": 0, "rs_only": 0, "ref_only": 0, "neither": 0}

    with open(args.input_file) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            total += 1

            # ClinVar truth
            cv_p = row["variant_prot"]
            cv_spdi = row["spdi"]

            # Analysis
            rs_p_raw = row["rs_p"]
            ref_p_raw = row["ref_p"]

            if rs_p_raw.startswith("ERR:RefMismatch"):
                rs_ref_mismatch += 1
            if rs_p_raw.startswith("ERR:Parse"):
                rs_parse_err += 1
            if ref_p_raw.startswith("ERR:Parse"):
                ref_parse_err += 1

            # Protein matches
            rs_p_ok = is_p_match(rs_p_raw, cv_p)
            ref_p_ok = is_p_match(ref_p_raw, cv_p)

            if rs_p_ok:
                rs_p_match += 1
            if ref_p_ok:
                ref_p_match += 1

            if rs_p_ok and ref_p_ok:
                p_stats["both"] += 1
            elif rs_p_ok:
                p_stats["rs_only"] += 1
            elif ref_p_ok:
                p_stats["ref_only"] += 1
            else:
                p_stats["neither"] += 1

            # SPDI matches
            rs_spdi_ok = row["rs_spdi"] == cv_spdi
            ref_spdi_ok = row["ref_spdi"] == cv_spdi

            if rs_spdi_ok:
                rs_spdi_match += 1
            if ref_spdi_ok:
                ref_spdi_match += 1

            if rs_spdi_ok and ref_spdi_ok:
                spdi_stats["both"] += 1
            elif rs_spdi_ok:
                spdi_stats["rs_only"] += 1
            elif ref_spdi_ok:
                spdi_stats["ref_only"] += 1
            else:
                spdi_stats["neither"] += 1

    if total == 0:
        print("No variants processed.")
        return

    # Calculate percentages
    rs_p_pct = rs_p_match / total * 100
    ref_p_pct = ref_p_match / total * 100
    rs_spdi_pct = rs_spdi_match / total * 100
    ref_spdi_pct = ref_spdi_match / total * 100

    # Determine bolds
    rs_p_str = f"{rs_p_pct:.3f}%"
    ref_p_str = f"{ref_p_pct:.3f}%"
    if rs_p_pct > ref_p_pct:
        rs_p_str = f"**{rs_p_str}**"
    elif ref_p_pct > rs_p_pct:
        ref_p_str = f"**{ref_p_str}**"

    rs_spdi_str = f"{rs_spdi_pct:.3f}%"
    ref_spdi_str = f"{ref_spdi_pct:.3f}%"
    if rs_spdi_pct > ref_spdi_pct:
        rs_spdi_str = f"**{rs_spdi_str}**"
    elif ref_spdi_pct > rs_spdi_pct:
        ref_spdi_str = f"**{ref_spdi_str}**"

    rs_err_str = f"{rs_parse_err:,}"
    ref_err_str = f"{ref_parse_err:,}"
    if rs_parse_err < ref_parse_err:
        rs_err_str = f"**{rs_err_str}**"
    elif ref_parse_err < rs_parse_err:
        ref_err_str = f"**{ref_err_str}**"

    report = []
    report.append(f"### Validation Results ({total:,} variants)")
    report.append("")
    report.append("Summary of results comparing `weaver` and `ref-hgvs` against ClinVar ground truth:")
    report.append("")
    report.append("| Implementation | Protein Match | SPDI Match  | Parse Errors |")
    report.append("| :------------- | :-----------: | :---------: | :----------: |")
    report.append(f"| weaver         |  {rs_p_str}  | {rs_spdi_str} | {rs_err_str} |")
    report.append(f"| ref-hgvs       |  {ref_p_str}  | {ref_spdi_str} | {ref_err_str} |")
    report.append("")
    report.append(f"RefSeq Data Mismatches: {rs_ref_mismatch:,} ({rs_ref_mismatch / total * 100:.1f}%)")
    report.append("")
    report.append("#### Protein Translation Agreement")
    report.append("")
    report.append("|                     | ref-hgvs Match | ref-hgvs Mismatch |")
    report.append("| :------------------ | :------------: | :---------------: |")
    report.append(f"| **weaver Match**    |     {p_stats['both']:,}     |     {p_stats['rs_only']:,}     |")
    report.append(f"| **weaver Mismatch** |     {p_stats['ref_only']:,}     |     {p_stats['neither']:,}     |")
    report.append("")
    report.append("#### SPDI Mapping Agreement")
    report.append("")
    report.append("|                     | ref-hgvs Match | ref-hgvs Mismatch |")
    report.append("| :------------------ | :------------: | :---------------: |")
    report.append(f"| **weaver Match**    |     {spdi_stats['both']:,}     |     {spdi_stats['rs_only']:,}     |")
    report.append(f"| **weaver Mismatch** |     {spdi_stats['ref_only']:,}     |     {spdi_stats['neither']:,}     |")

    out_text = "\n".join(report)
    print(out_text)

    if args.update_readme:
        readme_path = Path("README.md")
        if not readme_path.exists():
            # Try parent if running from a subdir
            readme_path = Path("../README.md")

        if readme_path.exists():
            content = readme_path.read_text()
            # Regex to find the section from ### Validation Results until the next section (starting with - **Variant Equivalence**)
            # or end of file.
            pattern = re.compile(r"### Validation Results \(.*?\).*?(?=\n- \*\*Variant Equivalence\*\*)", re.DOTALL)
            if pattern.search(content):
                new_content = pattern.sub(out_text.replace("\\", "\\\\"), content)
                readme_path.write_text(new_content)
                print(f"\n[Updated {readme_path}]")
            else:
                print("\n[Error: Could not find Validation Results section in README.md]")
        else:
            print("\n[Error: README.md not found]")


if __name__ == "__main__":
    main()
