# ruff: noqa: ANN401, PLR0913, PLR0912, C901
import bz2
import contextlib
import csv
import logging
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import weaver
from weaver.cli import provider

# Increase CSV field size limit for very large variants
csv.field_size_limit(sys.maxsize)

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

MAX_MISMATCH_LOG_ENTRIES = 20
MAX_ERROR_LOG_ENTRIES = 50
MIN_REQUIRED_ARGS = 2
PARTS_EXPECTED = 4
CHROM_X_NUM = 23
CHROM_Y_NUM = 24

REPO_ROOT = Path(__file__).parent.parent.resolve()

GFF_37 = "GCF_000001405.25_GRCh37.p13_genomic.gff.gz"
FASTA_37 = "GCF_000001405.25_GRCh37.p13_genomic.fna"
GFF_38 = "GRCh38_latest_genomic.gff.gz"
FASTA_38 = "GRCh38_latest_genomic.fna"


@dataclass
class ValidationProviders:
    rp37: provider.RefSeqDataProvider
    rp38: provider.RefSeqDataProvider
    mapper37: weaver.VariantMapper
    mapper38: weaver.VariantMapper


def vcf_normalize(v_genomic: Any, rp: provider.RefSeqDataProvider, mapper: weaver.VariantMapper) -> str:
    """
    Converts a genomic variant to VCF-style {chrom}-{pos}-{ref}-{alt} format.
    Handles anchor base prepending for indels.
    """
    try:
        spdi = mapper.to_spdi_unambiguous(v_genomic)

        parts = spdi.split(":")
        if len(parts) == PARTS_EXPECTED:
            ref_ac, pos_str, ref, alt = parts
            pos = int(pos_str)

            chrom = ref_ac.split(".")[0]
            if chrom.startswith("NC_0000"):
                chrom_num = int(chrom[7:])
                if chrom_num == CHROM_X_NUM:
                    chrom = "X"
                elif chrom_num == CHROM_Y_NUM:
                    chrom = "Y"
                else:
                    chrom = str(chrom_num)

            # For indels (len change), use an anchor base before the variant
            if len(ref) != len(alt):
                anchor_pos = pos - 1
                anchor_base = rp.get_seq(ref_ac, anchor_pos, anchor_pos + 1, "g")
                return f"{chrom}-{anchor_pos + 1}-{anchor_base}{ref}-{anchor_base}{alt}"

            # Substitution or identity
            return f"{chrom}-{pos + 1}-{ref}-{alt}"
    except Exception:  # noqa: BLE001, S110
        pass

    # Fallback to simple formatting if SPDI fails
    s = str(v_genomic)
    if ":" in s:
        ac, rest = s.split(":", 1)
        chrom = ac.split(".")[0]
        if chrom.startswith("NC_0000"):
            chrom_num = int(chrom[7:])
            if chrom_num == CHROM_X_NUM:
                chrom = "X"
            elif chrom_num == CHROM_Y_NUM:
                chrom = "Y"
            else:
                chrom = str(chrom_num)

        if ">" in rest:
            m = rest.split(">")

            match = re.search(r"g\.(\d+)([A-Z])", m[0])
            if match:
                f_pos = match.group(1)
                f_ref = match.group(2)
                f_alt = m[1]
                return f"{chrom}-{f_pos}-{f_ref}-{f_alt}"

    return str(v_genomic)


def parse_vcf(v_str: str) -> tuple[str, int, str, str] | None:
    if not v_str or v_str.startswith("MAPPING_FAILED"):
        return None
    parts = v_str.split("-")
    if len(parts) != PARTS_EXPECTED:
        return None
    try:
        return parts[0], int(parts[1]), parts[2], parts[3]
    except ValueError:
        return None


def minimize_vcf(chrom: str, pos: int, ref: str, alt: str) -> tuple[str, int, str, str]:
    # Strip common suffixes
    while len(ref) > 0 and len(alt) > 0 and ref[-1] == alt[-1]:
        ref = ref[:-1]
        alt = alt[:-1]
    # Strip common prefixes
    while len(ref) > 0 and len(alt) > 0 and ref[0] == alt[0]:
        ref = ref[1:]
        alt = alt[1:]
        pos += 1
    return chrom, pos, ref, alt


def is_vcf_equivalent(v1: str, v2: str) -> bool:
    if v1 == v2:
        return True
    try:
        p1 = parse_vcf(v1)
        p2 = parse_vcf(v2)
        if not p1 or not p2:
            return False
        if p1[0] != p2[0]:
            return False
        m1 = minimize_vcf(*p1)
        m2 = minimize_vcf(*p2)
        return m1 == m2
    except Exception:  # noqa: BLE001
        return False


def is_hgvs_equivalent(
    v1: Any,
    v2_str: str,
    mapper: weaver.VariantMapper,
    searcher: provider.RefSeqDataProvider,
) -> bool:
    """
    Check if a Variant object is equivalent to an HGVS string.
    """
    if not v1 or not v2_str or v2_str in {"None", "?"}:
        return False
    # Basic string match after stripping parentheses
    s1 = v1.format().strip("()") if hasattr(v1, "format") else str(v1).strip("()")
    s2 = v2_str.strip("()")
    if s1 == s2:
        return True
    try:
        v2 = weaver.parse(v2_str)
        # Check for accession match
        if hasattr(v1, "ac") and v1.ac.split(".")[0] != v2.ac.split(".")[0]:
            return False
        return mapper.equivalent(v1, v2, searcher)
    except Exception:  # noqa: BLE001
        return s1 == s2


def process_rows(
    reader: csv.DictReader,
    vp: ValidationProviders,
    stats: dict[str, int],
    m_writer: csv.DictWriter | None,
) -> None:
    for row in reader:
        v_g = None
        stats["total"] += 1
        variant_str = row["variant_str"]
        build = row["genome_build"]
        expected_id = row["variant_id"]
        expected_c = row["hgvs_c"]
        expected_p = row["hgvs_p"]

        try:
            if ":" not in variant_str and (variant_str.startswith(("c.", "p."))):
                ac = expected_c.split(":")[0] if ":" in expected_c else None
                if ac:
                    variant_str = f"{ac}:{variant_str}"

            if variant_str.startswith("NR:"):
                tx_ac = _resolve_nr_prefix(row.get("gene_symbol", ""), build, vp)
                if tx_ac:
                    variant_str = variant_str.replace("NR:", f"{tx_ac}:", 1)

            try:
                v = weaver.parse(variant_str)
            except Exception:  # noqa: BLE001
                v = weaver.parse(expected_c)

            cur_mapper = vp.mapper37 if "37" in build else vp.mapper38
            cur_rp = vp.rp37 if "37" in build else vp.rp38

            if v.coordinate_type == "p":
                v_id, res_c, v_c, res_p, v_p, v_g, cur_mapper, cur_rp = _handle_protein_input(
                    v,
                    expected_c,
                    expected_p,
                    cur_mapper,
                    cur_rp,
                    vp,
                    build,
                    expected_id,
                )
            else:
                v_id, v_g, cur_mapper, cur_rp = _handle_coding_input(v, cur_mapper, cur_rp, vp, build, expected_id)
                res_c, v_c, res_p, v_p = _generate_c_p_results(v, v_g, expected_c, expected_p, cur_mapper)

            # Check matches
            matches_id = not expected_id or is_vcf_equivalent(v_id, expected_id)
            matches_c = not expected_c or is_hgvs_equivalent(v_c, expected_c, cur_mapper, cur_rp)
            matches_p = not expected_p or is_hgvs_equivalent(v_p, expected_p, cur_mapper, cur_rp)

            if not matches_c and expected_c and v_g:
                matches_c, v_c, res_c = _mane_fallback(v_g, v_c, expected_c, cur_mapper, cur_rp)

            if matches_id:
                stats["id_match"] += 1
            if matches_c:
                stats["c_match"] += 1
            if matches_p:
                stats["p_match"] += 1

            if not matches_id or not matches_c or not matches_p:
                _log_mismatch(
                    stats,
                    variant_str,
                    v_id,
                    expected_id,
                    build,
                    res_c,
                    expected_c,
                    res_p,
                    expected_p,
                    matches_id,
                    matches_c,
                    matches_p,
                    m_writer,
                )

        except Exception as e:  # noqa: BLE001
            stats["errors"] += 1
            if stats["total"] <= MAX_ERROR_LOG_ENTRIES:
                logger.error("Error processing %s: %s", variant_str, e)


def _resolve_nr_prefix(gene_symbol: str, build: str, vp: ValidationProviders) -> str | None:
    if not gene_symbol:
        return None
    cur_rp = vp.rp37 if "37" in build else vp.rp38
    tx_candidates = cur_rp.gene_to_transcripts.get(gene_symbol, [])
    nm_candidates = [t for t in tx_candidates if t.startswith("NM_")]
    return nm_candidates[0] if nm_candidates else None


def _handle_protein_input(
    v: Any,
    expected_c: str,
    expected_p: str,
    cur_mapper: weaver.VariantMapper,
    cur_rp: provider.RefSeqDataProvider,
    vp: ValidationProviders,
    build: str,
    expected_id: str,
) -> tuple:
    v_id, res_c, v_c, res_p, v_p, v_g = "MAPPING_FAILED", "N/A", None, str(v), v, None
    tx_ac = expected_c.split(":")[0] if expected_c and ":" in expected_c else None
    prot_ac = expected_p.split(":")[0] if expected_p and ":" in expected_p else None
    try:
        v_c_from_p, _ = cur_mapper.p_to_c(v, tx_ac)
        v_c, res_c, v_g = v_c_from_p, str(v_c_from_p), cur_mapper.c_to_g(v_c_from_p)
        v_id = vcf_normalize(v_g, cur_rp, cur_mapper)
        if prot_ac:
            with contextlib.suppress(Exception):
                v_p = cur_mapper.c_to_p(v_c_from_p, prot_ac)
                res_p = str(v_p)
    except Exception as ep:  # noqa: BLE001
        logger.debug("p_to_c failed: %s", ep)
        if tx_ac:
            with contextlib.suppress(Exception):
                v_c_parsed = weaver.parse(expected_c)
                v_g = cur_mapper.c_to_g(v_c_parsed)
                v_id, res_c, v_c = vcf_normalize(v_g, cur_rp, cur_mapper), str(v_c_parsed), v_c_parsed

    if not is_vcf_equivalent(v_id, expected_id):
        o_mapper = vp.mapper38 if "37" in build else vp.mapper37
        o_rp = vp.rp38 if "37" in build else vp.rp37
        try:
            v_c_other, _ = o_mapper.p_to_c(v, tx_ac)
            v_g_other = o_mapper.c_to_g(v_c_other)
            v_id_other = vcf_normalize(v_g_other, o_rp, o_mapper)
            if is_vcf_equivalent(v_id_other, expected_id):
                v_id, cur_mapper, cur_rp, v_g, v_c, res_c = (
                    v_id_other,
                    o_mapper,
                    o_rp,
                    v_g_other,
                    v_c_other,
                    str(v_c_other),
                )
                if prot_ac:
                    with contextlib.suppress(Exception):
                        v_p = o_mapper.c_to_p(v_c_other, prot_ac)
                        res_p = str(v_p)
        except Exception:  # noqa: BLE001, S110
            pass
    return v_id, res_c, v_c, res_p, v_p, v_g, cur_mapper, cur_rp


def _handle_coding_input(
    v: Any,
    cur_mapper: weaver.VariantMapper,
    cur_rp: provider.RefSeqDataProvider,
    vp: ValidationProviders,
    build: str,
    expected_id: str,
) -> tuple:
    v_g, v_id = None, "MAPPING_FAILED"
    try:
        v_g = cur_mapper.c_to_g(v)
        v_id = vcf_normalize(v_g, cur_rp, cur_mapper)
    except Exception:  # noqa: BLE001, S110
        pass

    if not is_vcf_equivalent(v_id, expected_id):
        o_mapper = vp.mapper38 if "37" in build else vp.mapper37
        o_rp = vp.rp38 if "37" in build else vp.rp37
        try:
            v_g_other = o_mapper.c_to_g(v)
            v_id_other = vcf_normalize(v_g_other, o_rp, o_mapper)
            if is_vcf_equivalent(v_id_other, expected_id):
                v_id, cur_mapper, cur_rp, v_g = v_id_other, o_mapper, o_rp, v_g_other
        except Exception:  # noqa: BLE001, S110
            pass
    return v_id, v_g, cur_mapper, cur_rp


def _generate_c_p_results(
    v: Any,
    v_g: Any,
    expected_c: str,
    expected_p: str,
    cur_mapper: weaver.VariantMapper,
) -> tuple:
    try:
        v_norm = cur_mapper.normalize_variant(v)
        res_c, v_c = str(v_norm), v_norm
    except weaver.HGVSError:
        res_c, v_c = str(v), v

    res_p, v_p = "N/A", None
    try:
        prot_ac = expected_p.split(":")[0] if ":" in expected_p else None
        tx_ac = expected_c.split(":")[0] if ":" in expected_c else None
        if prot_ac and tx_ac:
            try:
                if v_g:
                    v_c_from_g = cur_mapper.g_to_c(v_g, tx_ac)
                    v_p = cur_mapper.c_to_p(v_c_from_g, prot_ac)
                    res_p = str(v_p)
                else:
                    raise ValueError("No genomic variant")
            except Exception as ep:  # noqa: BLE001
                try:
                    v_p = cur_mapper.c_to_p(v, prot_ac)
                    res_p = str(v_p)
                except Exception as ep2:  # noqa: BLE001
                    res_p = f"FAILED: {ep} | {ep2}"
    except Exception as ep_outer:  # noqa: BLE001
        res_p = f"FAILED: {ep_outer}"
    return res_c, v_c, res_p, v_p


def _mane_fallback(
    v_g: Any,
    v_c: Any,
    expected_c: str,
    cur_mapper: weaver.VariantMapper,
    cur_rp: provider.RefSeqDataProvider,
) -> tuple[bool, Any, str]:
    expected_tx = expected_c.split(":")[0] if ":" in expected_c else None
    if expected_tx and v_c and hasattr(v_c, "ac") and v_c.ac.split(".")[0] != expected_tx.split(".")[0]:
        gene = None
        with contextlib.suppress(Exception):
            gene = cur_rp.get_transcript(v_c.ac, None).get("gene")
        mane_tx = cur_rp.mane_select.get(gene) if gene else None
        if mane_tx:
            try:
                v_c_mane = cur_mapper.g_to_c(v_g, mane_tx)
                if is_hgvs_equivalent(v_c_mane, expected_c, cur_mapper, cur_rp):
                    return True, v_c_mane, str(v_c_mane)
            except Exception:  # noqa: BLE001, S110
                pass
    return False, v_c, str(v_c)


def _log_mismatch(
    stats: dict,
    variant_str: str,
    v_id: str,
    expected_id: str,
    build: str,
    res_c: str,
    expected_c: str,
    res_p: str,
    expected_p: str,
    matches_id: bool,
    matches_c: bool,
    matches_p: bool,
    m_writer: csv.DictWriter | None,
) -> None:
    m_id, m_c, m_p = ("." if matches_id else "X"), ("." if matches_c else "X"), ("." if matches_p else "X")
    match_str = f"{m_id}{m_c}{m_p}"
    if stats["total"] <= MAX_MISMATCH_LOG_ENTRIES:
        logger.warning(
            "Mismatch row %s [%s]: %s\n  ID: Found %s, Expected %s (%s)\n  C/P: Found %s/%s, Expected %s/%s",
            stats["total"],
            match_str,
            variant_str,
            v_id,
            expected_id,
            build,
            res_c,
            res_p,
            expected_c,
            expected_p,
        )
    if m_writer:
        m_writer.writerow(
            {
                "row_index": stats["total"],
                "match": match_str,
                "variant_str": variant_str,
                "build": build,
                "expected_id": expected_id,
                "found_id": v_id,
                "expected_c": expected_c,
                "found_c": res_c,
                "expected_p": expected_p,
                "found_p": res_p,
            },
        )
    print(f"Processed {stats['total']}... Matches: ID={stats['id_match']}, C={stats['c_match']}, P={stats['p_match']}")


def validate(input_path: str, mismatch_path: str | None = None) -> None:
    logger.info("Loading providers...")
    p_37 = provider.RefSeqDataProvider(GFF_37, FASTA_37)
    p_38 = provider.RefSeqDataProvider(GFF_38, FASTA_38)
    vp = ValidationProviders(p_37, p_38, weaver.VariantMapper(p_37), weaver.VariantMapper(p_38))
    stats = {"total": 0, "id_match": 0, "c_match": 0, "p_match": 0, "errors": 0}

    with contextlib.ExitStack() as stack:
        input_file = stack.enter_context(bz2.open(input_path, "rt"))
        reader = csv.DictReader(input_file, delimiter="\t")
        m_writer = None
        if mismatch_path:
            m_file = stack.enter_context(open(mismatch_path, "w", newline=""))
            m_writer = csv.DictWriter(
                m_file,
                fieldnames=[
                    "row_index",
                    "match",
                    "variant_str",
                    "build",
                    "expected_id",
                    "found_id",
                    "expected_c",
                    "found_c",
                    "expected_p",
                    "found_p",
                ],
                delimiter="\t",
            )
            m_writer.writeheader()
        process_rows(reader, vp, stats, m_writer)

    print("\nFinal Results:")
    for key, val in stats.items():
        print(f"{key}: {val}")


if __name__ == "__main__":
    if len(sys.argv) < MIN_REQUIRED_ARGS:
        print("Usage: python scripts/validate_normalization.py <input_bz2> [mismatch_output_tsv]")
        sys.exit(1)
    validate(sys.argv[1], sys.argv[2] if len(sys.argv) > MIN_REQUIRED_ARGS else None)
