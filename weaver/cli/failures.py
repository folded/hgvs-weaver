"""Extraction of weaver mapping failures."""

import csv


def is_p_match(pred: str, truth: str) -> bool:
    """Checks if predicted protein matches truth with fuzzy logic."""
    if pred.startswith("ERR:"):
        return False
    p = pred.replace("(", "").replace(")", "").split(":")[-1]
    t = truth.replace("(", "").replace(")", "").split(":")[-1]
    return p == t or ("Ter" in p and "Ter" in t) or ("fs" in p and "fs" in t) or ("=" in p and "=" in t)


def main() -> None:
    """Main failures entry point."""
    input_file = "clinvar_full_validation_100k.tsv"
    count = 0
    max_count = 20
    print(f"{'Variant':30} | {'ClinVar':30} | {'weaver':30} | {'ref-hgvs':30}")
    print("-" * 130)

    try:
        with open(input_file) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                rs_ok = is_p_match(row["rs_p"], row["variant_prot"])
                ref_ok = is_p_match(row["ref_p"], row["variant_prot"])

                if ref_ok and not rs_ok:
                    v = row["variant_nuc"]
                    cv = row["variant_prot"].split(":")[-1]
                    rs = row["rs_p"]
                    ref = row["ref_p"]
                    print(f"{v:30} | {cv:30} | {rs:30} | {ref:30}")
                    count += 1
                    if count >= max_count:
                        break
    except FileNotFoundError:
        print(f"Error: {input_file} not found.")


if __name__ == "__main__":
    main()
