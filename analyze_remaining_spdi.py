import csv


def get_variant_type(row: dict[str, str]) -> str:
    nuc = row["variant_nuc"]
    v_type = "Unknown"
    if "dup" in nuc:
        v_type = "dup"
    elif "delins" in nuc:
        v_type = "delins"
    elif "del" in nuc:
        v_type = "del"
    elif "ins" in nuc:
        v_type = "ins"
    elif "[" in nuc:
        v_type = "repeat"
    elif ">" in nuc:
        v_type = "subst"
    return v_type


def analyze_spdi_mismatches(results_file: str) -> None:
    mismatches: list[dict[str, str]] = []
    unsupported: dict[str, int] = {}

    with open(results_file) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            truth = row["spdi"]
            weaver = row["rs_spdi"]

            if truth != weaver:
                mismatches.append(row)
                if weaver.startswith("ERR:Unsupported operation"):
                    reason = weaver.split("ERR:Unsupported operation: ")[1]
                    unsupported[reason] = unsupported.get(reason, 0) + 1
                elif weaver.startswith("ERR:"):
                    prefix = weaver.split(":")[1]
                    unsupported[f"Other ERR: {prefix}"] = unsupported.get(f"Other ERR: {prefix}", 0) + 1
                else:
                    unsupported["Real Mismatch"] = unsupported.get("Real Mismatch", 0) + 1
                    v_type = get_variant_type(row)
                    unsupported[f"Mismatch: {v_type}"] = unsupported.get(f"Mismatch: {v_type}", 0) + 1

    print(f"Total Mismatches: {len(mismatches)}")
    print("\nCategorized failures:")
    for reason, count in sorted(unsupported.items(), key=lambda x: x[1], reverse=True):
        print(f"  {count:5} - {reason}")

    print("\nTop 10 Real Mismatches:")
    real_mismatches = [m for m in mismatches if not m["rs_spdi"].startswith("ERR:")]
    max_to_show = 10
    for m in real_mismatches[:max_to_show]:
        print(f"  {m['variant_nuc']:30} Truth: {m['spdi']:50} Weaver: {m['rs_spdi']}")


if __name__ == "__main__":
    analyze_spdi_mismatches("results_spdi_fix.tsv")
