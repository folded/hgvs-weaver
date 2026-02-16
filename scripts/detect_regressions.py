import csv
import sys
from pathlib import Path


def get_key(row):
    return f"{row['variant_nuc']}|{row['spdi']}"


def analyze_regressions(base_file, tip_file):
    base_results = {}
    with open(base_file) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            # Normalize protein for comparison
            p = row["rs_p"].replace("(", "").replace(")", "")
            if ":" in p:
                p = p.split(":")[-1]
            # Consider it "passed" if rs_p matches variant_prot OR rs_spdi matches spdi
            # (Adjust logic based on what we define as 'pass')
            gt_p = row["variant_prot"].replace("(", "").replace(")", "")
            if ":" in gt_p:
                gt_p = gt_p.split(":")[-1]

            p_pass = p == gt_p
            s_pass = row["rs_spdi"] == row["spdi"]
            base_results[get_key(row)] = (p_pass, s_pass, p, row["rs_spdi"])

    regressions = []
    with open(tip_file) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            key = get_key(row)
            if key in base_results:
                b_p_pass, b_s_pass, b_p, b_s = base_results[key]

                # Check Protein Regression
                curr_p = row["rs_p"].replace("(", "").replace(")", "")
                if ":" in curr_p:
                    curr_p = curr_p.split(":")[-1]
                gt_p = row["variant_prot"].replace("(", "").replace(")", "")
                if ":" in gt_p:
                    gt_p = gt_p.split(":")[-1]
                curr_p_pass = curr_p == gt_p

                # Check SPDI Regression
                curr_s_pass = row["rs_spdi"] == row["spdi"]

                if (b_p_pass and not curr_p_pass) or (b_s_pass and not curr_s_pass):
                    regressions.append(
                        {
                            "variant": row["variant_nuc"],
                            "truth_spdi": row["spdi"],
                            "regression_type": "Protein" if (b_p_pass and not curr_p_pass) else "SPDI",
                            "base_val": b_p if (b_p_pass and not curr_p_pass) else b_s,
                            "tip_val": curr_p if (b_p_pass and not curr_p_pass) else row["rs_spdi"],
                            "truth_val": gt_p if (b_p_pass and not curr_p_pass) else row["spdi"],
                        }
                    )

    return regressions


def main():
    if len(sys.argv) < 3:
        print("Usage: detect_regressions.py <base_tsv> <tip_tsv>")
        sys.exit(1)

    regs = analyze_regressions(sys.argv[1], sys.argv[2])
    print(f"Total Regressions Found: {len(regs)}")
    for r in regs[:20]:  # Show top 20
        print(
            f"[{r['regression_type']}] {r['variant']} (Base: {r['base_val']} -> Tip: {r['tip_val']}, Truth: {r['truth_val']})"
        )


if __name__ == "__main__":
    main()
