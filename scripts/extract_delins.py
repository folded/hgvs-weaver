import csv
import sys
from scripts.classify_failures import classify


def extract_delins(input_file, limit=50):
    with open(input_file, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        count = 0
        for row in reader:
            cat = classify(row)
            if cat == "AA Mismatch (delins)":
                print(f"CAT: {cat}")
                print(f"NUC: {row['variant_nuc']}")
                print(f"GT:  {row['variant_prot']}")
                print(f"W:   {row['rs_p']}")
                print(f"R:   {row['ref_p']}")
                print("-" * 40)
                count += 1
                if count >= limit:
                    break


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python extract_delins.py <results_tsv>")
    else:
        extract_delins(sys.argv[1])
