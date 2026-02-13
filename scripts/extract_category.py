import csv
import sys

from scripts.classify_failures import classify


def extract_category(input_file, category, limit=50) -> None:
    with open(input_file) as f:
        reader = csv.DictReader(f, delimiter="\t")
        count = 0
        for row in reader:
            cat = classify(row)
            if cat == category:
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
    if len(sys.argv) < 3:
        print("Usage: python extract_category.py <results_tsv> <category>")
    else:
        extract_category(sys.argv[1], sys.argv[2])
