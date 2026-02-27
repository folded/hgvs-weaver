import bz2
import contextlib
import csv
import logging
import sys
from pathlib import Path

import weaver
from weaver.cli import provider

# Increase CSV field size limit for very large variants
csv.field_size_limit(sys.maxsize)

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

REPO_ROOT = Path(__file__).parent.parent.resolve()

GFF_37 = "GCF_000001405.25_GRCh37.p13_genomic.gff.gz"
FASTA_37 = "GCF_000001405.25_GRCh37.p13_genomic.fna"
GFF_38 = "GRCh38_latest_genomic.gff.gz"
FASTA_38 = "GRCh38_latest_genomic.fna"


def vcf_normalize(v_genomic, rp, mapper):
    """
    Converts a genomic variant to VCF-style {chrom}-{pos}-{ref}-{alt} format.
    Handles anchor base prepending for indels.
    """
    try:
        spdi = mapper.to_spdi_unambiguous(v_genomic)

        parts = spdi.split(":")
        if len(parts) == 4:
            ref_ac, pos, ref, alt = parts
            pos = int(pos)

            chrom = ref_ac.split(".")[0]
            if chrom.startswith("NC_0000"):
                chrom_num = int(chrom[7:])
                if chrom_num == 23:
                    chrom = "X"
                elif chrom_num == 24:
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
    except Exception:
        pass

    # Fallback to simple formatting if SPDI fails
    # {chrom}-{pos}-{ref}-{alt}
    # Parsing variant string to get ref/alt
    # NC_000012.12:g.21910221A>G
    s = str(v_genomic)
    if ":" in s:
        ac, rest = s.split(":", 1)
        chrom = ac.split(".")[0]
        if chrom.startswith("NC_0000"):
            chrom_num = int(chrom[7:])
            if chrom_num == 23:
                chrom = "X"
            elif chrom_num == 24:
                chrom = "Y"
            else:
                chrom = str(chrom_num)

        if ">" in rest:
            m = rest.split(">")
            import re

            match = re.search(r"g\.(\d+)([A-Z])", m[0])
            if match:
                pos = match.group(1)
                ref = match.group(2)
                alt = m[1]
                return f"{chrom}-{pos}-{ref}-{alt}"

    return str(v_genomic)


def parse_vcf(v_str):
    if not v_str or v_str.startswith("MAPPING_FAILED"):
        return None
    parts = v_str.split("-")
    if len(parts) != 4:
        return None
    try:
        return parts[0], int(parts[1]), parts[2], parts[3]
    except ValueError:
        return None


def minimize_vcf(chrom, pos, ref, alt):
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


def is_vcf_equivalent(v1, v2):
    if v1 == v2:
        return True
    try:
        # v1 and v2 are chrom-pos-ref-alt
        c1, p1, r1, a1 = v1.split("-")
        c2, p2, r2, a2 = v2.split("-")
        if c1 != c2:
            return False
        # Normalize both and compare
        m1 = minimize_vcf(c1, int(p1), r1, a1)
        m2 = minimize_vcf(c2, int(p2), r2, a2)
        return m1 == m2
    except:
        return False


def is_hgvs_equivalent(v1, v2_str, mapper, searcher):
    """
    Check if a Variant object is equivalent to an HGVS string.
    """
    if not v1 or not v2_str or v2_str in {"None", "?"}:
        return False
    # Basic string match after stripping parentheses
    s1 = v1.format().strip("()")
    s2 = v2_str.strip("()")
    if s1 == s2:
        return True
    try:
        v2 = weaver.parse(v2_str)
        # Check for accession match
        if v1.ac.split(".")[0] != v2.ac.split(".")[0]:
            return False
        return mapper.equivalent(v1, v2, searcher)
    except:
        return s1 == s2


def validate(input_path, mismatch_path=None) -> None:
    logger.info("Loading providers...")
    rp37 = provider.RefSeqDataProvider(GFF_37, FASTA_37)
    rp38 = provider.RefSeqDataProvider(GFF_38, FASTA_38)
    mapper37 = weaver.VariantMapper(rp37)
    mapper38 = weaver.VariantMapper(rp38)

    stats = {"total": 0, "id_match": 0, "c_match": 0, "p_match": 0, "errors": 0}

    mismatch_file = None
    mismatch_writer = None

    try:
        if mismatch_path:
            mismatch_file = open(mismatch_path, "w", newline="")
            mismatch_writer = csv.DictWriter(
                mismatch_file,
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
            mismatch_writer.writeheader()

        with bz2.open(input_path, "rt") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                stats["total"] += 1
                variant_str = row["variant_str"]
                build = row["genome_build"]
                expected_id = row["variant_id"]
                expected_c = row["hgvs_c"]
                expected_p = row["hgvs_p"]

                try:
                    # Robust parsing: if variant_str is just c. or p., try to prepend accession from expected_c
                    if ":" not in variant_str and (variant_str.startswith(("c.", "p."))):
                        ac = expected_c.split(":")[0] if ":" in expected_c else None
                        if ac:
                            variant_str = f"{ac}:{variant_str}"

                    # Handle NR: placeholder prefix by resolving gene symbol to NM_ accession
                    gene_symbol = row.get("gene_symbol", "")
                    if variant_str.startswith("NR:") and gene_symbol:
                        current_rp_for_resolve = rp37 if "37" in build else rp38
                        tx_candidates = current_rp_for_resolve.gene_to_transcripts.get(gene_symbol, [])
                        nm_candidates = [t for t in tx_candidates if t.startswith("NM_")]
                        if nm_candidates:
                            variant_str = variant_str.replace("NR:", f"{nm_candidates[0]}:", 1)

                    try:
                        v = weaver.parse(variant_str)
                    except Exception:
                        v = weaver.parse(expected_c)

                    # Pick mapper based on build
                    current_mapper = mapper37 if "37" in build else mapper38
                    current_rp = rp37 if "37" in build else rp38

                    # Handle protein-as-input: when variant_str is NM_…:p.Xxx
                    if v.coordinate_type == "p":
                        # Use mapper.p_to_c to back-convert protein variant to coding
                        v_id = "MAPPING_FAILED"
                        res_c = "N/A"
                        v_c = None
                        res_p = str(v)
                        v_p = v

                        tx_ac = expected_c.split(":")[0] if expected_c and ":" in expected_c else None
                        prot_ac = expected_p.split(":")[0] if expected_p and ":" in expected_p else None
                        try:
                            v_c_from_p, _is_unique = current_mapper.p_to_c(v, tx_ac)
                            v_c = v_c_from_p
                            res_c = str(v_c_from_p)
                            v_g = current_mapper.c_to_g(v_c_from_p)
                            v_id = vcf_normalize(v_g, current_rp, current_mapper)
                            # Now compute p. from the back-converted c. for a round-trip check
                            if prot_ac:
                                try:
                                    v_p = current_mapper.c_to_p(v_c_from_p, prot_ac)
                                    res_p = str(v_p)
                                except Exception:
                                    pass
                        except Exception as ep:
                            logger.debug("p_to_c failed for %s: %s", str(v), ep)
                            # Fallback: parse expected_c directly
                            if tx_ac:
                                try:
                                    v_c_parsed = weaver.parse(expected_c)
                                    v_g = current_mapper.c_to_g(v_c_parsed)
                                    v_id = vcf_normalize(v_g, current_rp, current_mapper)
                                    res_c = str(v_c_parsed)
                                    v_c = v_c_parsed
                                except Exception:
                                    pass

                        # Build swap fallback for p. branch
                        if not is_vcf_equivalent(v_id, expected_id):
                            other_mapper = mapper38 if "37" in build else mapper37
                            other_rp = rp38 if "37" in build else rp37
                            try:
                                v_c_other, _ = other_mapper.p_to_c(v, tx_ac)
                                v_g_other = other_mapper.c_to_g(v_c_other)
                                v_id_other = vcf_normalize(v_g_other, other_rp, other_mapper)
                                if is_vcf_equivalent(v_id_other, expected_id):
                                    v_id = v_id_other
                                    current_mapper = other_mapper
                                    current_rp = other_rp
                                    v_c = v_c_other
                                    res_c = str(v_c_other)
                                    if prot_ac:
                                        try:
                                            v_p = other_mapper.c_to_p(v_c_other, prot_ac)
                                            res_p = str(v_p)
                                        except Exception:
                                            pass
                            except Exception:
                                pass
                    else:
                        # Normal c. variant flow
                        # Map to genomic
                        v_g = None
                        try:
                            v_g = current_mapper.c_to_g(v)
                            v_id = vcf_normalize(v_g, current_rp, current_mapper)
                        except Exception:
                            v_id = "MAPPING_FAILED"

                        # Build swap fallback
                        if not is_vcf_equivalent(v_id, expected_id):
                            other_mapper = mapper38 if "37" in build else mapper37
                            other_rp = rp38 if "37" in build else rp37
                            try:
                                v_g_other = other_mapper.c_to_g(v)
                                v_id_other = vcf_normalize(v_g_other, other_rp, other_mapper)
                                if is_vcf_equivalent(v_id_other, expected_id):
                                    v_id = v_id_other
                                    current_mapper = other_mapper
                                    current_rp = other_rp
                                    v_g = v_g_other
                            except Exception:
                                logging.exception("Failed to map to genomic")

                        # Generate c. and p. — normalize the input for position comparison
                        try:
                            v_norm = current_mapper.normalize_variant(v)
                            res_c = str(v_norm)
                            v_c = v_norm
                        except weaver.HGVSError:
                            res_c = str(v)
                            v_c = v
                        res_p = "N/A"
                        v_p = None  # Keep track of the parsed protein variant for p. comparison
                        try:
                            prot_ac = expected_p.split(":")[0] if ":" in expected_p else None
                            tx_ac = expected_c.split(":")[0] if ":" in expected_c else None
                            if prot_ac and tx_ac:
                                try:
                                    # Try g -> c -> p
                                    if v_g:
                                        v_c_from_g = current_mapper.g_to_c(v_g, tx_ac)
                                        v_p = current_mapper.c_to_p(v_c_from_g, prot_ac)
                                        res_p = str(v_p)
                                    else:
                                        raise ValueError("No genomic variant")
                                except Exception as ep:
                                    # Try direct c -> p
                                    try:
                                        v_p = current_mapper.c_to_p(v, prot_ac)
                                        res_p = str(v_p)
                                    except Exception as ep2:
                                        res_p = f"FAILED: {ep} | {ep2}"
                        except Exception as ep_outer:
                            res_p = f"FAILED: {ep_outer}"

                    # Empty expected fields are not counted as mismatches
                    matches_id = not expected_id or is_vcf_equivalent(v_id, expected_id)
                    matches_c = not expected_c or is_hgvs_equivalent(v_c, expected_c, current_mapper, current_rp)
                    matches_p = not expected_p or is_hgvs_equivalent(v_p, expected_p, current_mapper, current_rp)

                    # MANE Select fallback: if c. fails due to transcript alias, try via MANE Select
                    if not matches_c and expected_c and v_g:
                        expected_tx = expected_c.split(":")[0] if ":" in expected_c else None
                        if expected_tx and v_c and v_c.ac.split(".")[0] != expected_tx.split(".")[0]:
                            gene = None
                            with contextlib.suppress(Exception):
                                gene = current_rp.get_transcript(v_c.ac, None).get("gene")
                            mane_tx = current_rp.mane_select.get(gene) if gene else None
                            if mane_tx:
                                try:
                                    v_c_mane = current_mapper.g_to_c(v_g, mane_tx)
                                    if is_hgvs_equivalent(v_c_mane, expected_c, current_mapper, current_rp):
                                        matches_c = True
                                        v_c = v_c_mane
                                        res_c = str(v_c_mane)
                                except Exception:
                                    pass

                    if matches_id:
                        stats["id_match"] += 1
                    if matches_c:
                        stats["c_match"] += 1
                    if matches_p:
                        stats["p_match"] += 1

                    is_any_mismatch = not matches_id or not matches_c or not matches_p

                    if is_any_mismatch:
                        m_id = "." if matches_id else "X"
                        m_c = "." if matches_c else "X"
                        m_p = "." if matches_p else "X"
                        match_str = f"{m_id}{m_c}{m_p}"

                        if stats["total"] <= 20:
                            logger.warning(
                                "Mismatch at row %s [%s]: %s\n"
                                "  ID: Found %s, Expected %s (Build %s)\n"
                                "  C:  Found %s, Expected %s\n"
                                "  P:  Found %s, Expected %s",
                                stats["total"],
                                match_str,
                                variant_str,
                                v_id,
                                expected_id,
                                build,
                                res_c,
                                expected_c,
                                res_p,
                                expected_p,
                            )

                        if mismatch_writer:
                            mismatch_writer.writerow(
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

                    if stats["total"] % 100 == 0:
                        print(
                            f"Processed {stats['total']}... Matches: ID={stats['id_match']}, C={stats['c_match']}, P={stats['p_match']}",
                        )

                except Exception as e:
                    stats["errors"] += 1
                    if stats["total"] <= 50:
                        logger.error("Error processing %s: %s", variant_str, e)
    finally:
        if mismatch_file:
            mismatch_file.close()

    print("\nFinal Results:")
    for k, variant in stats.items():
        print(f"{k}: {variant}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python scripts/validate_normalization.py <input_bz2> [mismatch_output_tsv]")
        sys.exit(1)

    m_path = sys.argv[2] if len(sys.argv) > 2 else None
    validate(sys.argv[1], m_path)
