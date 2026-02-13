import csv
import collections


def analyze_unsupported_vs_ref(input_file):
    total_unsupported = 0
    ref_hgvs_performance = collections.Counter()
    total = 0
    examples = []

    with open(input_file) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            total += 1
            cv_spdi = row["spdi"]
            rs_spdi = row["rs_spdi"]
            ref_spdi = row["ref_spdi"]

            if rs_spdi.startswith("ERR:Unsupported operation"):
                total_unsupported += 1

                if ref_spdi == cv_spdi:
                    ref_hgvs_performance["Matched ClinVar"] += 1
                    if len(examples) < 10:
                        examples.append({"variant": row["variant_nuc"], "cv": cv_spdi, "ref": ref_spdi})
                elif ref_spdi.startswith("ERR:"):
                    ref_hgvs_performance["Ref Error"] += 1
                else:
                    ref_hgvs_performance["Ref Unique Mismatch"] += 1

    print(f"Total Weaver 'Unsupported Operation' cases: {total_unsupported}")
    print("\nRef-hgvs performance in these cases:")
    for cat, count in ref_hgvs_performance.most_common():
        print(f"  {cat:30}: {count:6,} ({count / total_unsupported * 100:5.2f}%)")

    if examples:
        print("\nExamples where Weaver is Unsupported but Ref-HGVS matches ClinVar:")
        for ex in examples:
            print(f"  Variant: {ex['variant']}")
            print(f"    Truth: {ex['cv']}")
            print(f"    Ref:   {ex['ref']}")
            print("-" * 20)


if __name__ == "__main__":
    analyze_unsupported_vs_ref("results_100k_final.tsv")
