import csv
import collections


def analyze_spdi_mismatches(input_file, limit_examples=10):
    mismatches = 0
    categories = collections.Counter()
    ref_hgvs_behavior = collections.Counter()
    total = 0
    examples = []

    with open(input_file) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            total += 1
            cv_spdi = row["spdi"]
            rs_spdi = row["rs_spdi"]
            ref_spdi = row["ref_spdi"]

            if rs_spdi != cv_spdi:
                mismatches += 1

                cat = "Other"
                if rs_spdi.startswith("ERR:"):
                    err_type = rs_spdi.split(":")[1]
                    cat = f"Error: {err_type}"
                elif rs_spdi == "PANIC":
                    cat = "Panic"
                else:
                    cv_parts = cv_spdi.split(":")
                    rs_parts = rs_spdi.split(":")
                    if len(cv_parts) != 4 or len(rs_parts) != 4:
                        cat = "Malformed SPDI"
                    elif cv_parts[0] != rs_parts[0]:
                        cat = "Accession Mismatch"
                    elif cv_parts[1] != rs_parts[1]:
                        cat = "Position Mismatch"
                    elif cv_parts[2] != rs_parts[2]:
                        cat = "Reference Seq Mismatch"
                    elif cv_parts[3] != rs_parts[3]:
                        cat = "Alternate Seq Mismatch"

                categories[cat] += 1

                if ref_spdi == cv_spdi:
                    ref_hgvs_behavior["Matched ClinVar"] += 1
                elif ref_spdi == rs_spdi:
                    ref_hgvs_behavior["Matched Weaver (Shared Mismatch)"] += 1
                elif ref_spdi.startswith("ERR:"):
                    ref_hgvs_behavior["Error"] += 1
                else:
                    ref_hgvs_behavior["Unique Mismatch"] += 1

                if cat == "Position Mismatch" and len(examples) < limit_examples:
                    examples.append({"variant": row["variant_nuc"], "cv": cv_spdi, "rs": rs_spdi, "ref": ref_spdi})

    print(f"Total Variants: {total:,}")
    print(f"Total SPDI Mismatches: {mismatches:,} ({mismatches / total * 100:.2f}%)")
    print("\nBreakdown of Weaver Mismatches:")
    for cat, count in categories.most_common():
        print(f"  {cat:30}: {count:6,} ({count / mismatches * 100:5.2f}%)")

    print("\nRef-hgvs behavior when Weaver mismatches:")
    for cat, count in ref_hgvs_behavior.most_common():
        print(f"  {cat:30}: {count:6,} ({count / mismatches * 100:5.2f}%)")

    if examples:
        print("\nExamples of Position Mismatch (Weaver vs Truth vs Ref):")
        for ex in examples:
            print(f"  Variant: {ex['variant']}")
            print(f"    Truth:  {ex['cv']}")
            print(f"    Weaver: {ex['rs']}")
            print(f"    Ref:    {ex['ref']}")
            print("-" * 20)


if __name__ == "__main__":
    analyze_spdi_mismatches("results_100k_final.tsv")
