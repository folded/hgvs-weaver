import csv
import json
import sys
from pathlib import Path

# Add repo root to sys.path to find weaver
REPO_ROOT = Path(__file__).parent.parent.resolve()
sys.path.append(str(REPO_ROOT))

import weaver
from weaver.cli import provider


def normalize_p(p_str: str | None) -> str:
    if not p_str or p_str == "ERR" or p_str.startswith("ERR:"):
        return p_str if p_str is not None else ""
    # Strip accession
    if ":" in p_str:
        p_str = p_str.split(":")[-1]
    # Strip parentheses
    p_str = p_str.replace("(", "").replace(")", "")
    # Standardize Ter/*
    return p_str.replace("Ter", "*")


def analyze(filepath: str, gff: str = "GRCh38_latest_genomic.gff.gz", fasta: str = "GRCh38_latest_genomic.fna") -> None:
    print(f"Loading provider from {gff} and {fasta}...")
    rp = provider.RefSeqDataProvider(gff, fasta)
    mapper = weaver.VariantMapper(rp)

    total = 0
    # Weaver counters
    w_identity = 0
    w_analogous = 0

    # Ref counters
    ref_identity = 0
    ref_analogous = 0

    with open(filepath, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            total += 1
            gt_p_str = row["variant_prot"]
            rs_p_str = row["rs_p"]
            ref_p_str = row["ref_p"]

            gt_p_norm = normalize_p(gt_p_str)
            rs_p_norm = normalize_p(rs_p_str)
            ref_p_norm = normalize_p(ref_p_str)

            # Weaver Check
            if rs_p_norm == gt_p_norm:
                w_identity += 1
            elif rs_p_str and not rs_p_str.startswith("ERR") and gt_p_str and not gt_p_str.startswith("ERR"):
                try:
                    ac_p = gt_p_str.split(":")[0] if ":" in gt_p_str else row["variant_nuc"].split(":")[0]
                    v_rs_str = rs_p_str if ":" in rs_p_str else f"{ac_p}:{rs_p_str}"
                    v_rs = weaver.parse(v_rs_str)
                    v_gt = weaver.parse(gt_p_str)
                    if mapper.equivalent(v_rs, v_gt, rp):
                        w_analogous += 1
                except Exception:
                    pass

            # Ref Check
            if ref_p_norm == gt_p_norm:
                ref_identity += 1
            elif ref_p_str and not ref_p_str.startswith("ERR") and gt_p_str and not gt_p_str.startswith("ERR"):
                try:
                    ac_p = gt_p_str.split(":")[0] if ":" in gt_p_str else row["variant_nuc"].split(":")[0]
                    v_ref_str = ref_p_str if ":" in ref_p_str else f"{ac_p}:{ref_p_str}"
                    v_ref = weaver.parse(v_ref_str)
                    v_gt = weaver.parse(gt_p_str)
                    if mapper.equivalent(v_ref, v_gt, rp):
                        ref_analogous += 1
                except Exception:
                    pass

    if total == 0:
        print("No variants processed.")
        return

    print(f"Total variants: {total}")
    print("-" * 30)
    print("Weaver Performance:")
    print(f"  Identity:  {w_identity} ({w_identity / total:.2%})")
    print(f"  Analogous: {w_analogous} ({w_analogous / total:.2%})")
    print(f"  Total:     {w_identity + w_analogous} ({(w_identity + w_analogous) / total:.2%})")

    print("-" * 30)
    print("Ref-HGVS Performance:")
    print(f"  Identity:  {ref_identity} ({ref_identity / total:.2%})")
    print(f"  Analogous: {ref_analogous} ({ref_analogous / total:.2%})")
    print(f"  Total:     {ref_identity + ref_analogous} ({(ref_identity + ref_analogous) / total:.2%})")

    # Output JSON for history update
    stats = {
        "total": total,
        "w_identity": w_identity,
        "w_analogous": w_analogous,
        "ref_identity": ref_identity,
        "ref_analogous": ref_analogous,
    }

    with open("current_stats.json", "w", encoding="utf-8") as f:
        json.dump(stats, f, indent=2)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python analyze_validation.py <results_file> [gff] [fasta]")
        sys.exit(1)

    results_file = sys.argv[1]
    gff = sys.argv[2] if len(sys.argv) > 2 else "GRCh38_latest_genomic.gff.gz"
    fasta = sys.argv[3] if len(sys.argv) > 3 else "GRCh38_latest_genomic.fna"

    analyze(results_file, gff, fasta)
