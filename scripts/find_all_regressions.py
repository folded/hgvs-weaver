import csv
import os
import glob
from pathlib import Path

RESULTS_DIR = Path("benchmark_results")
TIP_FILE = RESULTS_DIR / "results_a551cff.tsv"


def normalize_p(p):
    if not p:
        return ""
    p = p.replace("(", "").replace(")", "")
    if ":" in p:
        p = p.split(":")[-1]
    return p


def get_key(row):
    return row.get("variant_nuc", "")


def is_prot_match(row):
    v_prot = normalize_p(row.get("variant_prot", ""))
    rs_p = normalize_p(row.get("rs_p", ""))
    return v_prot and v_prot == rs_p


def read_results(file_path):
    rows = []
    with open(file_path, "r", encoding="utf-8") as f:
        # Need to handle duplicate headers manually because DictReader will overwrite them
        line = f.readline()
        if not line:
            return []
        headers = line.strip().split("\t")

        # We only care about the FIRST occurrence of these headers
        # but DictReader uses the LAST.
        # Actually, let's just use the first few we care about.

        for line in f:
            values = line.strip().split("\t")
            row = dict(zip(headers, values))
            rows.append(row)
    return rows


def main():
    # 1. Load Tip failures
    tip_failures = {}
    tip_rows = read_results(TIP_FILE)
    for row in tip_rows:
        if not is_prot_match(row):
            key = get_key(row)
            if key:
                tip_failures[key] = {"row": row, "succeeded_in": []}

    print(f"Loaded {len(tip_failures)} protein failures from tip.")

    # 2. Check other files for successes
    other_files = glob.glob(str(RESULTS_DIR / "results_*.tsv"))
    for file_path in other_files:
        if "a551cff" in file_path:
            continue

        commit = os.path.basename(file_path).replace("results_", "").replace(".tsv", "")
        print(f"Checking {commit}...")

        other_rows = read_results(file_path)
        for row in other_rows:
            key = get_key(row)
            if key in tip_failures and is_prot_match(row):
                if commit not in tip_failures[key]["succeeded_in"]:
                    tip_failures[key]["succeeded_in"].append(commit)

    # 3. Filter for those that succeeded at least once
    regressions = {k: v for k, v in tip_failures.items() if v["succeeded_in"]}

    print(f"\nFound {len(regressions)} variants that succeeded at some point but fail at tip.")

    # 4. Save results
    output_path = "protein_regressions_all_history.txt"
    with open(output_path, "w") as f:
        f.write(f"Total History Regressions (Protein): {len(regressions)}\n")
        f.write("-" * 80 + "\n")
        # Sort by variant_nuc
        sorted_keys = sorted(regressions.keys())
        for key in sorted_keys:
            info = regressions[key]
            row = info["row"]
            f.write(f"Variant: {key}\n")
            f.write(f"  Truth: {row.get('variant_prot')}\n")
            f.write(f"  Tip Result: {row.get('rs_p')}\n")
            f.write(f"  Succeeded in: {', '.join(info['succeeded_in'])}\n")
            f.write("\n")

    print(f"Report saved to {output_path}")


if __name__ == "__main__":
    main()
