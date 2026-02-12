# /// script
# dependencies = [
#   "pysam",
#   "tqdm",
#   "parsley",
#   "bioutils",
# ]
# ///

"""Full validation script against ClinVar variants."""

import argparse
import concurrent.futures
import csv
import sys
import typing

try:
    import hgvs.parser
    import hgvs.variantmapper
except ImportError:
    print(
        "Error: 'hgvs' package not found. Please install it manually (e.g. without dependencies to avoid psycopg2) with:",
    )
    print("  pip install hgvs --no-deps")
    sys.exit(1)

try:
    import tqdm
except ImportError:
    print("Error: 'tqdm' package not found. Please install it manually with: pip install tqdm")
    sys.exit(1)

import weaver

from . import provider

_rp: provider.RefSeqDataProvider | None = None
_rs_mapper: weaver.VariantMapper | None = None
_ref_vm: hgvs.variantmapper.VariantMapper | None = None
_ref_hp: hgvs.parser.Parser | None = None


def init_worker(gff: str, fasta: str) -> None:
    """Initializes global mappers for worker processes."""
    global _rp, _rs_mapper, _ref_vm, _ref_hp
    _rp = provider.RefSeqDataProvider(gff, fasta)
    _rs_mapper = weaver.VariantMapper(_rp)
    _ref_hdp = provider.ReferenceHgvsDataProvider(_rp)
    _ref_vm = hgvs.variantmapper.VariantMapper(_ref_hdp)
    _ref_hp = hgvs.parser.Parser()


def hgvs_lib_to_spdi(v: typing.Any, data_provider: typing.Any) -> str | None:
    """Converts a standard hgvs library Variant object to SPDI string format."""
    try:
        ac = v.ac
        start_1 = v.posedit.pos.start.base
        end_1 = v.posedit.pos.end.base if v.posedit.pos.end else start_1
        if hasattr(v.posedit.edit, "ref"):
            ref = v.posedit.edit.ref or ""
            alt = v.posedit.edit.alt or ""
            if not ref or ref.isdigit():
                ref = data_provider.get_seq(ac, start_1 - 1, end_1, "g")
            if v.posedit.edit.type == "ins":
                return f"{ac}:{start_1}:{ref}:{alt}"
            return f"{ac}:{start_1 - 1}:{ref}:{alt}"
        return "UnsupportedType"
    except Exception:
        return "ERR:SPDI"


def process_variant(row: dict[str, str]) -> dict[str, str]:
    """Maps a single variant using both weaver and ref-hgvs for comparison."""
    nuc_hgvs = row["variant_nuc"]
    spdi_ac = row["spdi"].split(":")[0]

    rs_p = "ERR"
    rs_spdi = "ERR"
    ref_p = "ERR"
    ref_spdi = "ERR"

    # weaver block
    try:
        v_rs_raw = weaver.parse(nuc_hgvs)
        if not _rs_mapper:
            row["rs_p"] = row["rs_spdi"] = "ERR:MapperNotInit"
            return row

        v_rs = _rs_mapper.normalize_variant(v_rs_raw)
        try:
            if v_rs.coordinate_type == "c":
                v_p = _rs_mapper.c_to_p(v_rs)
                rs_p = v_p.format().split(":")[-1]
        except Exception as e:
            rs_p = f"ERR:{e!s}"

        try:
            if v_rs.coordinate_type != "g":
                vg_rs = _rs_mapper.c_to_g(v_rs, spdi_ac)
                vg_rs = _rs_mapper.normalize_variant(vg_rs)  # Normalize in genomic space for SPDI
            else:
                vg_rs = v_rs
            rs_spdi = vg_rs.to_spdi(_rp)
        except Exception as e:
            rs_spdi = f"ERR:{e!s}"
    except Exception as e:
        rs_p = rs_spdi = f"ERR:{e!s}"
    except BaseException:  # Catch absolutely everything including panics
        rs_p = rs_spdi = "PANIC"

    # ref-hgvs block
    try:
        if not _ref_hp or not _ref_vm:
            row["ref_p"] = row["ref_spdi"] = "ERR:RefMapperNotInit"
            return row

        v_ref = _ref_hp.parse_hgvs_variant(nuc_hgvs)
        try:
            if v_ref.type == "c":
                v_p_ref = _ref_vm.c_to_p(v_ref)
                ref_p = str(v_p_ref).split(":")[-1]
        except Exception as e:
            ref_p = f"ERR:{e!s}"

        try:
            vg_ref = _ref_vm.c_to_g(v_ref, spdi_ac) if v_ref.type != "g" else v_ref
            # Normalize in genomic space for SPDI using weaver normalizer
            try:
                vg_ref_rs = weaver.parse(str(vg_ref))
                vg_ref_rs = _rs_mapper.normalize_variant(vg_ref_rs)
                ref_spdi = vg_ref_rs.to_spdi(_rp)
            except Exception:
                # Fallback to unnormalized if weaver parsing fails
                ref_spdi = hgvs_lib_to_spdi(vg_ref, _rp)
        except Exception as e:
            ref_spdi = f"ERR:{e!s}"
    except Exception:
        ref_p = ref_spdi = "ERR:Parse"
    except BaseException:
        ref_p = ref_spdi = "PANIC"

    row["rs_p"] = rs_p
    row["rs_spdi"] = rs_spdi
    row["ref_p"] = ref_p
    row["ref_spdi"] = ref_spdi
    return row


def main() -> None:
    """Main entry point for validation."""
    parser = argparse.ArgumentParser(description="Full validation against ClinVar variants.")
    parser.add_argument("input_file", help="Input ClinVar TSV file.")
    parser.add_argument("--max-variants", type=int, default=None, help="Maximum variants to process.")
    parser.add_argument("--output-file", default="clinvar_full_validation.tsv", help="Output validation TSV.")
    parser.add_argument("--gff", default="GRCh38_latest_genomic.gff.gz", help="Reference GFF file.")
    parser.add_argument("--fasta", default="GCF_000001405.40_GRCh38.p14_genomic.fna", help="Reference FASTA file.")
    parser.add_argument("--workers", type=int, default=4, help="Number of worker processes.")
    args = parser.parse_args()

    with open(args.input_file) as f_in:
        reader = csv.DictReader(f_in, delimiter="\t")
        fieldnames = [*(reader.fieldnames or []), "rs_p", "rs_spdi", "ref_p", "ref_spdi"]
        rows: list[dict[str, str]] = (
            [next(reader) for _ in range(args.max_variants)] if args.max_variants else list(reader)
        )

    print(f"Processing {len(rows)} variants with ProcessPool...")

    with open(args.output_file, "w", newline="") as f_out:
        writer = csv.DictWriter(f_out, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()

        with concurrent.futures.ProcessPoolExecutor(
            max_workers=args.workers,
            initializer=init_worker,
            initargs=(args.gff, args.fasta),
        ) as executor:
            # map instead of executor.map to catch task-level errors
            results_iter = executor.map(process_variant, rows)

            pbar = tqdm.tqdm(total=len(rows))
            while True:
                try:
                    row_res = next(results_iter)
                    writer.writerow(row_res)
                    pbar.update(1)
                except StopIteration:
                    break
                except Exception as e:
                    print(f"\nWorker crashed: {e}")
                    pbar.update(1)
            pbar.close()


if __name__ == "__main__":
    main()
