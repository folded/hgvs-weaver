import csv
import sys

from weaver.cli.provider import RefSeqDataProvider


def verify_equivalence(results_file) -> None:
    provider = RefSeqDataProvider(
        gff_path="GRCh38_latest_genomic.gff.gz",
        fasta_path="GCF_000001405.40_GRCh38.p14_genomic.fna",
    )
    mismatches = 0
    equivalent = 0
    real_diff = 0

    with open(results_file) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for i, row in enumerate(reader):
            truth = row["spdi"]
            weaver_spdi = row["rs_spdi"]

            if not truth or not weaver_spdi or weaver_spdi.startswith("ERR:"):
                continue

            if truth == weaver_spdi:
                continue

            mismatches += 1

            # Verify equivalence
            # SPDI: ac:pos:del:ins
            try:
                t_parts = truth.split(":")
                w_parts = weaver_spdi.split(":")

                if t_parts[0] != w_parts[0]:
                    real_diff += 1
                    continue

                ac = t_parts[0]
                t_pos, t_del, t_ins = int(t_parts[1]), t_parts[2], t_parts[3]
                w_pos, w_del, w_ins = int(w_parts[1]), w_parts[2], w_parts[3]

                # Check if applying both to reference yields same result
                # We need a bit of context sequence
                start = min(t_pos, w_pos)
                end = max(t_pos + len(t_del), w_pos + len(w_del))
                padding = 50

                context_start = max(0, start - padding)
                context_end = end + padding

                ref_seq = provider.get_seq(ac, context_start, context_end, "g")

                # Apply T
                t_rel_pos = t_pos - context_start
                t_seq = ref_seq[:t_rel_pos] + t_ins + ref_seq[t_rel_pos + len(t_del) :]

                # Apply W
                w_rel_pos = w_pos - context_start
                w_seq = ref_seq[:w_rel_pos] + w_ins + ref_seq[w_rel_pos + len(w_del) :]

                if t_seq == w_seq:
                    equivalent += 1
                else:
                    real_diff += 1
                    if real_diff <= 5:
                        print(f"Mismatch {i}: {row['variant_nuc']}")
                        print(f"  Truth:  {truth}")
                        print(f"  Weaver: {weaver_spdi}")
                        # print(f"  Ref:    {ref_seq[t_rel_pos-5:t_rel_pos+len(t_del)+5]}")
            except Exception:
                # print(f"Error checking {row['variant_nuc']}: {e}")
                real_diff += 1

    print(f"\nSummary of {mismatches} lexical mismatches:")
    print(f"  Equivalent (biologically same): {equivalent} ({equivalent / mismatches:.1%})")
    print(f"  Real differences:               {real_diff}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python verify_spdi_equiv.py <results_file>")
        sys.exit(1)
    verify_equivalence(sys.argv[1])
