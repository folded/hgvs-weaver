import collections
import csv
import sys


def clean_hgvs(s):
    if not s:
        return ""
    # Remove accession prefix
    if ":" in s:
        s = s.split(":")[-1]
    # Remove parentheses
    s = s.replace("(", "").replace(")", "")
    # Standardize Ter/*
    return s.replace("Ter", "*")


def classify(row) -> str:
    rs_p = row["rs_p"]
    ref_p = row["ref_p"]
    gt_p = row["variant_prot"]

    c_rs = clean_hgvs(rs_p)
    c_gt = clean_hgvs(gt_p)
    c_ref = clean_hgvs(ref_p)

    if c_rs == c_ref and c_rs != "":
        return "Parity Match"

    if c_rs == c_gt and c_rs != "":
        return "ClinVar Match"

    if rs_p == "ERR" or rs_p.startswith("ERR:"):
        # Categorize known errors
        if "Transcript" in rs_p and "not found" in rs_p:
            return "Provider Error (Missing Transcript)"

        if "TranscriptMismatch" in rs_p:
            return "Weaver Error: Reference sequence mismatch"

        err_str = rs_p.split(":")[-1] if ":" in rs_p else "Generic"
        err_str = err_str.lstrip()
        return f"Weaver Error: {err_str}"

    # Biological/Nomenclature differences
    if "fs" in rs_p or "fs" in gt_p:
        return "Frameshift Difference"

    if "delins" in rs_p or "delins" in gt_p:
        return "AA Mismatch (delins)"

    if "[" in row["variant_nuc"] or "dup" in row["variant_nuc"]:
        return "Repeat/Duplication Mismatch"

    if "*" in rs_p or "*" in gt_p:
        return "Nonsense Notation Mismatch"

    return "Other Mismatch"


if len(sys.argv) < 2:
    print("Usage: classify_failures.py <results_tsv>")
    sys.exit(1)

stats = {}
total = 0
mismatches = collections.defaultdict(list)
success_count = 0

with open(sys.argv[1]) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        total += 1
        cat = classify(row)
        stats[cat] = stats.get(cat, 0) + 1
        if cat in ["Parity Match", "ClinVar Match"]:
            success_count += 1
        elif len(mismatches[cat]) < 100:
            mismatches[cat].append(row)

print(f"Total variants: {total}")
print(f"Total Successes: {success_count} ({(success_count / total) * 100:.2f}%)")
print(f"Total Failures:  {total - success_count} ({((total - success_count) / total) * 100:.2f}%)")

print("\nFailure breakdown (excluding successful matches):")
for cat, count in sorted(stats.items(), key=lambda x: x[1], reverse=True):
    print(f"  {cat}: {count} ({(count / (total - success_count)) * 100:.2f}% of failures)")

if mismatches:
    print("\nSample Mismatches:")
    # Sort samples to put delins first if we want to focus on them
    for cat, rows in sorted(mismatches.items(), key=lambda x: -len(x[1])):
        if cat == "Provider Error (Missing Transcript)":
            continue
        print(f"CAT: {cat}")
        for row in rows:
            print("-" * 10)
            print(f"NUC: {row['variant_nuc']}")
            print(f"GT:  {row['variant_prot']}")
            print(f"W:   {row['rs_p']}")
            print(f"R:   {row['ref_p']}")
        print("-" * 40)
