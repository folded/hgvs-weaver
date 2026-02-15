import collections
import csv


def get_mismatch_category(cv_spdi: str, rs_spdi: str) -> str:
    category = "Other"
    if rs_spdi.startswith("ERR:"):
        err_type = rs_spdi.split(":")[1]
        category = f"Error: {err_type}"
    elif rs_spdi == "PANIC":
        category = "Panic"
    else:
        cv_parts = cv_spdi.split(":")
        rs_parts = rs_spdi.split(":")
        expected_parts = 4
        if len(cv_parts) != expected_parts or len(rs_parts) != expected_parts:
            category = "Malformed SPDI"
        elif cv_parts[0] != rs_parts[0]:
            category = "Accession Mismatch"
        elif cv_parts[1] != rs_parts[1]:
            category = "Position Mismatch"
        elif cv_parts[2] != rs_parts[2]:
            category = "Reference Seq Mismatch"
        elif cv_parts[3] != rs_parts[3]:
            category = "Alternate Seq Mismatch"

    return category


def collect_mismatch_stats(
    input_file: str,
    limit_examples: int,
) -> tuple[int, int, collections.Counter[str], collections.Counter[str], list[dict[str, str]]]:
    mismatches = 0
    categories: collections.Counter[str] = collections.Counter()
    ref_hgvs_behavior: collections.Counter[str] = collections.Counter()
    total = 0
    examples: list[dict[str, str]] = []

    with open(input_file) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            total += 1
            cv_spdi = row["spdi"]
            rs_spdi = row["rs_spdi"]
            ref_spdi = row["ref_spdi"]

            if rs_spdi != cv_spdi:
                mismatches += 1
                cat = get_mismatch_category(cv_spdi, rs_spdi)
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
    return total, mismatches, categories, ref_hgvs_behavior, examples


def analyze_spdi_mismatches(input_file: str, limit_examples: int = 10) -> None:
    total, mismatches, categories, ref_hgvs_behavior, examples = collect_mismatch_stats(input_file, limit_examples)

    if total == 0:
        print("No variants found.")
        return

    print(f"Total Variants: {total:,}")
    if mismatches == 0:
        print("No mismatches found!")
        return

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
