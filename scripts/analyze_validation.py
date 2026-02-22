import csv
import json
import logging
import sys
from pathlib import Path

# Add repo root to sys.path to find weaver
REPO_ROOT = Path(__file__).parent.parent.resolve()
sys.path.append(str(REPO_ROOT))

import weaver  # noqa: E402
from weaver.cli import provider  # noqa: E402

MIN_ARGS = 2
ARG_RESULTS_FILE_INDEX = 1
ARG_GFF_INDEX = 2
ARG_FASTA_INDEX = 3

DEFAULT_GFF = "GRCh38_latest_genomic.gff.gz"
DEFAULT_FASTA = "GRCh38_latest_genomic.fna"

logger = logging.getLogger(__name__)


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


def _check_equivalence(
    mapper: weaver.VariantMapper,
    rp: provider.RefSeqDataProvider,
    p_str: str,
    gt_p_str: str,
    row: dict[str, str],
) -> bool:
    try:
        ac_p = gt_p_str.split(":")[0] if ":" in gt_p_str else row["variant_nuc"].split(":")[0]
        v_str = p_str if ":" in p_str else f"{ac_p}:{p_str}"
        v = weaver.parse(v_str)
        v_gt = weaver.parse(gt_p_str)
        return mapper.equivalent(v, v_gt, rp)
    except (ValueError, RuntimeError, AttributeError):
        logger.debug("Equivalence check failed", exc_info=True)
        return False


def analyze(
    filepath: str,
    gff: str = DEFAULT_GFF,
    fasta: str = DEFAULT_FASTA,
) -> None:
    print(f"Loading provider from {gff} and {fasta}...")
    rp = provider.RefSeqDataProvider(gff, fasta)
    mapper = weaver.VariantMapper(rp)

    total = 0
    w_identity = 0
    w_analogous = 0
    ref_identity = 0
    ref_analogous = 0

    with open(filepath, newline="", encoding="utf-8") as f:
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
            elif (
                rs_p_str
                and not rs_p_str.startswith("ERR")
                and gt_p_str
                and not gt_p_str.startswith("ERR")
                and _check_equivalence(mapper, rp, rs_p_str, gt_p_str, row)
            ):
                w_analogous += 1

            # Ref Check
            if ref_p_norm == gt_p_norm:
                ref_identity += 1
            elif (
                ref_p_str
                and not ref_p_str.startswith("ERR")
                and gt_p_str
                and not gt_p_str.startswith("ERR")
                and _check_equivalence(mapper, rp, ref_p_str, gt_p_str, row)
            ):
                ref_analogous += 1

    if total == 0:
        print("No variants processed.")
        return

    _print_results(total, w_identity, w_analogous, ref_identity, ref_analogous)


def _print_results(total: int, w_id: int, w_an: int, r_id: int, r_an: int) -> None:
    print(f"Total variants: {total}")
    print("-" * 30)
    print("Weaver Performance:")
    print(f"  Identity:  {w_id} ({w_id / total:.2%})")
    print(f"  Analogous: {w_an} ({w_an / total:.2%})")
    print(f"  Total:     {w_id + w_an} ({(w_id + w_an) / total:.2%})")
    print("-" * 30)
    print("Ref-HGVS Performance:")
    print(f"  Identity:  {r_id} ({r_id / total:.2%})")
    print(f"  Analogous: {r_an} ({r_an / total:.2%})")
    print(f"  Total:     {r_id + r_an} ({(r_id + r_an) / total:.2%})")

    stats = {
        "total": total,
        "w_identity": w_id,
        "w_analogous": w_an,
        "ref_identity": r_id,
        "ref_analogous": r_an,
    }
    with open("current_stats.json", "w", encoding="utf-8") as f:
        json.dump(stats, f, indent=2)


if __name__ == "__main__":
    if len(sys.argv) < MIN_ARGS:
        print("Usage: python analyze_validation.py <results_file> [gff] [fasta]")
        sys.exit(1)

    results_file = sys.argv[1]
    gff_arg = sys.argv[2] if len(sys.argv) > MIN_ARGS else DEFAULT_GFF
    fasta_arg = sys.argv[3] if len(sys.argv) > (MIN_ARGS + 1) else DEFAULT_FASTA

    analyze(results_file, gff_arg, fasta_arg)
