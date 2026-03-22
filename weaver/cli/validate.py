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
import logging
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
# Pre-computed ferro normalize results: nuc_hgvs → normalized_string | "ERR:..."
_fh_results: dict[str, str] = {}


def init_worker(gff: str, fasta: str, fh_results: dict[str, str] | None = None) -> None:
    """Initializes global mappers for worker processes."""
    global _rp, _rs_mapper, _ref_vm, _ref_hp, _fh_results
    _rp = provider.RefSeqDataProvider(gff, fasta)
    _rs_mapper = weaver.VariantMapper(_rp)
    _ref_hdp = provider.ReferenceHgvsDataProvider(_rp)
    _ref_vm = hgvs.variantmapper.VariantMapper(_ref_hdp)
    _ref_hp = hgvs.parser.Parser()
    if fh_results:
        _fh_results = fh_results


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
    v_p = None
    rs_p = "ERR"
    rs_spdi = "ERR"
    try:
        v_rs_raw = weaver.parse(nuc_hgvs)
        if not _rs_mapper:
            res_row = row.copy()
            res_row.update(
                {
                    "rs_p": "ERR:MapperNotInit",
                    "rs_spdi": "ERR:MapperNotInit",
                    "rs_equiv": "ERR:MapperNotInit",
                    "ref_equiv": "ERR:MapperNotInit",
                },
            )
            return res_row

        v_rs = _rs_mapper.normalize_variant(v_rs_raw)
        try:
            if v_rs.coordinate_type == "c":
                v_p = _rs_mapper.c_to_p(v_rs)
                rs_p = v_p.format().split(":")[-1]
        except Exception as e:
            rs_p = f"ERR:{e!s}"

        try:
            rs_spdi = _rs_mapper.to_spdi(v_rs_raw, unambiguous=True)
        except Exception as e:
            rs_spdi = f"ERR:{e!s}"
    except Exception as e:
        rs_p = rs_spdi = f"ERR:{e!s}"
    except BaseException:
        rs_p = rs_spdi = "PANIC"

    # ref-hgvs block
    ref_p = "ERR"
    ref_spdi = "ERR"
    try:
        if not _ref_hp or not _ref_vm:
            ref_p = ref_spdi = "ERR:RefMapperNotInit"
        else:
            v_ref = _ref_hp.parse_hgvs_variant(nuc_hgvs)
            try:
                if v_ref.type == "c":
                    v_p_ref = _ref_vm.c_to_p(v_ref)
                    ref_p = str(v_p_ref).split(":")[-1]
            except Exception as e:
                ref_p = f"ERR:{e!s}"

            try:
                vg_ref = _ref_vm.c_to_g(v_ref, spdi_ac) if v_ref.type != "g" else v_ref
                ref_spdi = hgvs_lib_to_spdi(vg_ref, _rp)
            except Exception as e:
                ref_spdi = f"ERR:{e!s}"
    except Exception:
        ref_p = ref_spdi = "ERR:Parse"
    except BaseException:
        ref_p = ref_spdi = "PANIC"

    # ferro-hgvs block: look up pre-computed normalize result
    fh_parse = _fh_results.get(nuc_hgvs, "SKIP")

    # Equivalence Checks (Using weaver to judge both)
    rs_equiv = "Unknown"
    ref_equiv = "Unknown"
    gt_p_str = row.get("variant_prot", "")

    if gt_p_str and gt_p_str != "ERR" and not gt_p_str.startswith("ERR"):
        try:
            v_gt = weaver.parse(gt_p_str)
            # RS Equivalence
            if v_p:
                rs_equiv = str(_rs_mapper.equivalent_level(v_p, v_gt, _rp)).split(".")[-1]

            # REF Equivalence (judged by weaver)
            if ref_p and not ref_p.startswith("ERR"):
                try:
                    # Construct full protein string for weaver parsing
                    ref_p_val = f"{v_gt.ac}:{ref_p}" if ":" not in ref_p else ref_p
                    v_ref_p = weaver.parse(ref_p_val)
                    ref_equiv = str(_rs_mapper.equivalent_level(v_ref_p, v_gt, _rp)).split(".")[-1]
                except Exception:
                    try:
                        v_ref_p = weaver.parse(ref_p)
                        ref_equiv = str(_rs_mapper.equivalent_level(v_ref_p, v_gt, _rp)).split(".")[-1]
                    except Exception:
                        logging.exception("Failed to judge equivalence")
        except Exception:
            logging.exception("Failed to judge equivalence")

    res_row = row.copy()
    res_row.update(
        {
            "rs_p": rs_p or "",
            "rs_spdi": rs_spdi or "",
            "ref_p": ref_p or "",
            "ref_spdi": ref_spdi or "",
            "rs_equiv": rs_equiv,
            "ref_equiv": ref_equiv,
            "fh_parse": fh_parse,
        },
    )

    return res_row


def run_ferro_normalize(variants: list[str], reference_dir: str) -> dict[str, str]:
    """Batch-normalizes variants via the ferro CLI; returns nuc_hgvs → result mapping."""
    import shutil  # noqa: PLC0415
    import subprocess  # noqa: PLC0415
    import tempfile  # noqa: PLC0415

    ferro_bin = shutil.which("ferro")
    if not ferro_bin:
        print("Warning: 'ferro' binary not found in PATH; skipping ferro normalization.")
        return {}

    print(f"Running ferro normalize on {len(variants):,} variants (reference: {reference_dir})...")
    with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as tmp_in:
        tmp_in.write("\n".join(variants))
        tmp_in_path = tmp_in.name

    results: dict[str, str] = {}
    try:
        proc = subprocess.run(  # noqa: S603
            [ferro_bin, "normalize", "--reference", reference_dir, "-i", tmp_in_path, "-f", "text"],
            capture_output=True,
            text=True,
            check=False,
        )
        # ferro outputs one result per line; failures are printed to stderr with the input echoed
        # In text mode: success → normalized string; errors go to stderr as "ERR: <msg>"
        out_lines = proc.stdout.splitlines()
        for variant, out_line in zip(variants, out_lines):
            stripped = out_line.strip()
            if " -> " in stripped:
                # Changed: "NM_x:c.1A>G -> NM_x:c.1A>G"
                results[variant] = stripped.split(" -> ", 1)[1]
            elif stripped:
                results[variant] = stripped
            else:
                results[variant] = "ERR:EmptyOutput"
        # Any variant without a result line (e.g. crashed) gets marked as error
        for variant in variants:
            if variant not in results:
                results[variant] = "ERR:NoOutput"
    except Exception as e:
        print(f"Warning: ferro normalize failed: {e}")
    finally:
        import os  # noqa: PLC0415
        os.unlink(tmp_in_path)

    ok_count = sum(1 for v in results.values() if not v.startswith("ERR:"))
    print(f"ferro normalize: {ok_count:,}/{len(variants):,} succeeded.")
    return results


def main() -> None:
    """Main entry point for validation."""
    parser = argparse.ArgumentParser(description="Full validation against ClinVar variants.")
    parser.add_argument("input_file", help="Input ClinVar TSV file.")
    parser.add_argument("--max-variants", type=int, default=None, help="Maximum variants to process.")
    parser.add_argument("--output-file", default="clinvar_full_validation.tsv", help="Output validation TSV.")
    parser.add_argument("--gff", default="GRCh38_latest_genomic.gff.gz", help="Reference GFF file.")
    parser.add_argument("--fasta", default="GCF_000001405.40_GRCh38.p14_genomic.fna", help="Reference FASTA file.")
    parser.add_argument("--workers", type=int, default=4, help="Number of worker processes.")
    parser.add_argument(
        "--ferro-reference",
        default=None,
        help="Path to ferro reference directory (produced by 'ferro prepare'). "
             "When provided, ferro normalize is run in batch before validation.",
    )
    args = parser.parse_args()

    with open(args.input_file) as f_in:
        reader = csv.DictReader(f_in, delimiter="\t")
        base_fields = [
            f
            for f in (reader.fieldnames or [])
            if f not in {"rs_p", "rs_spdi", "ref_p", "ref_spdi", "rs_equiv", "ref_equiv", "equivalence_level", "fh_parse"}
        ]
        fieldnames = [*base_fields, "rs_p", "rs_spdi", "ref_p", "ref_spdi", "rs_equiv", "ref_equiv", "fh_parse"]
        rows: list[dict[str, str]] = (
            [next(reader) for _ in range(args.max_variants)] if args.max_variants else list(reader)
        )

    print(f"Processing {len(rows)} variants with ProcessPool...")

    # Pre-run ferro normalize in batch if reference directory supplied
    fh_results: dict[str, str] = {}
    if args.ferro_reference:
        all_nuc = [row["variant_nuc"] for row in rows]
        fh_results = run_ferro_normalize(all_nuc, args.ferro_reference)

    with open(args.output_file, "w", newline="") as f_out:
        writer = csv.DictWriter(f_out, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()

        with concurrent.futures.ProcessPoolExecutor(
            max_workers=args.workers,
            initializer=init_worker,
            initargs=(args.gff, args.fasta, fh_results),
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
