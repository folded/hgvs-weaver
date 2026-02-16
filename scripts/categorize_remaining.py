import csv


def normalize_p(p):
    if not p:
        return ""
    p = p.replace("(", "").replace(")", "")
    if ":" in p:
        p = p.split(":")[-1]
    return p


def main():
    mismatches = []
    with open("results_remaining.tsv", "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rs_p = normalize_p(row.get("rs_p", ""))
            truth_p = normalize_p(row.get("variant_prot", ""))

            if rs_p and truth_p and rs_p != truth_p:
                mismatches.append(row)

    print(f"Total protein mismatches found: {len(mismatches)}")

    # Categorization
    categories = {"synonymous": [], "premature_stop": [], "complex_delins": [], "other": []}

    for row in mismatches:
        v_p = normalize_p(row["rs_p"])
        h_p = normalize_p(row["variant_prot"])

        if "=" in h_p or "Ter=" in h_p:
            categories["synonymous"].append(row)
        elif "*" in h_p or "Ter" in h_p or "fs" in h_p:
            categories["premature_stop"].append(row)
        elif "delins" in h_p:
            categories["complex_delins"].append(row)
        else:
            categories["other"].append(row)

    for cat, list_rows in categories.items():
        print(f"Category {cat}: {len(list_rows)}")
        for r in list_rows[:3]:
            print(f"  {r['variant_nuc']} -> Result: {r['rs_p']} / Truth: {r['variant_prot']}")


if __name__ == "__main__":
    main()
