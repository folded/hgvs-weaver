import csv
import sys


def normalize_p(p_str: str | None) -> str:
    if not p_str or p_str == "ERR" or p_str.startswith("ERR:"):
        return p_str if p_str is not None else ""
    # Strip accession
    if ":" in p_str:
        p_str = p_str.split(":")[-1]
    # Strip parentheses
    return p_str.replace("(", "").replace(")", "")


def analyze(filepath: str) -> None:
    total = 0
    weaver_p_matches = 0
    weaver_spdi_matches = 0
    ref_p_matches = 0
    ref_spdi_matches = 0
    weaver_mismatch_ref_match = 0
    ref_mismatch_weaver_match = 0

    with open(filepath, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            total += 1
            gt_p = normalize_p(row["variant_prot"])
            rs_p = normalize_p(row["rs_p"])
            ref_p = normalize_p(row["ref_p"])

            rs_match = rs_p == gt_p
            ref_match = ref_p == gt_p

            if rs_match:
                weaver_p_matches += 1
            if ref_match:
                ref_p_matches += 1

            if ref_match and not rs_match:
                weaver_mismatch_ref_match += 1
            if rs_match and not ref_match:
                ref_mismatch_weaver_match += 1

            # SPDI match
            gt_spdi = row["spdi"]
            rs_spdi = row["rs_spdi"]
            ref_spdi = row["ref_spdi"]

            if rs_spdi == gt_spdi:
                weaver_spdi_matches += 1
            if ref_spdi == gt_spdi:
                ref_spdi_matches += 1

    if total == 0:
        print("No variants processed.")
        return

    print(f"Total variants: {total}")
    print(f"weaver Protein Match: {weaver_p_matches / total:.2%}")
    print(f"weaver SPDI Match: {weaver_spdi_matches / total:.2%}")
    print(f"ref-hgvs Protein Match: {ref_p_matches / total:.2%}")
    print(f"ref-hgvs SPDI Match: {ref_spdi_matches / total:.2%}")
    print(f"Weaver Mismatch but Ref Match Count: {weaver_mismatch_ref_match}")
    print(f"Ref Mismatch but Weaver Match Count: {ref_mismatch_weaver_match}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python analyze_validation.py <results_file>")
        sys.exit(1)
    analyze(sys.argv[1])
