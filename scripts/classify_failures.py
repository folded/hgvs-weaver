import collections
import csv
import os
import re
import sys
import typing

import weaver

# Try to initialize real provider if files exist
gff_path = "GRCh38_latest_genomic.gff.gz"
fasta_path = "GRCh38_latest_genomic.fna"

provider: typing.Any = None
searcher: typing.Any = None
mapper: weaver.VariantMapper = typing.cast("weaver.VariantMapper", None)

protein_gpff_path = "GRCh38_latest_protein.gpff.gz"
import gzip


class CombinedProvider:
    def __init__(self, gff_path: str, fasta_path: str, protein_gpff_path: str) -> None:
        from weaver.cli.provider import RefSeqDataProvider

        print(f"Loading RefSeq provider from {gff_path}...", file=sys.stderr)
        self.refseq = RefSeqDataProvider(gff_path, fasta_path)
        self.protein_seqs: dict[str, str] = {}
        if os.path.exists(protein_gpff_path):
            print(f"Loading protein sequences from {protein_gpff_path}...", file=sys.stderr)
            self._load_gpff(protein_gpff_path)

    def _load_gpff(self, path: str) -> None:
        current_ac: str | None = None
        in_origin = False
        seq_parts: list[str] = []

        with gzip.open(path, "rt", encoding="latin-1") as f:
            for line in f:
                if line.startswith("VERSION"):
                    current_ac = line.split()[1]
                elif line.startswith("ORIGIN"):
                    in_origin = True
                elif line.startswith("//"):
                    if current_ac:
                        # Remove whitespace and digits from sequence lines
                        clean_seq = "".join(seq_parts).upper()
                        self.protein_seqs[current_ac] = clean_seq
                        # Also store versionless key
                        base_ac = current_ac.split(".")[0]
                        self.protein_seqs[base_ac] = clean_seq
                    current_ac = None
                    in_origin = False
                    seq_parts = []
                elif in_origin:
                    parts = line.split()
                    if parts:
                        seq_parts.extend(parts[1:])
        print(f"Loaded {len(self.protein_seqs)} protein sequences.", file=sys.stderr)

    def get_identifier_type(self, identifier: str) -> "weaver.IdentifierType":
        res = self.refseq.get_identifier_type(identifier)
        if isinstance(res, str):
            # Map string to enum if necessary
            mapping = {
                "genomic_accession": weaver.IdentifierType.GenomicAccession,
                "transcript_accession": weaver.IdentifierType.TranscriptAccession,
                "protein_accession": weaver.IdentifierType.ProteinAccession,
                "gene_symbol": weaver.IdentifierType.GeneSymbol,
            }
            return mapping.get(res, weaver.IdentifierType.Unknown)
        return typing.cast("weaver.IdentifierType", res)

    def get_transcript(self, transcript_ac: str, reference_ac: str | None) -> typing.Any:
        return self.refseq.get_transcript(transcript_ac, reference_ac)

    def get_transcripts_for_region(self, chrom: str, start: int, end: int) -> list[str]:
        return self.refseq.get_transcripts_for_region(chrom, start, end)

    def get_seq(self, ac: str, start: int, end: int, kind: "weaver.IdentifierType") -> str | None:
        # Try exact match
        if ac == "AC":
            return None
        seq = self.protein_seqs.get(ac)

        # Try versionless match
        if seq is None and "." in ac:
            base_ac = ac.split(".")[0]
            seq = self.protein_seqs.get(base_ac)

        if seq is not None:
            s_idx = start if start >= 0 else 0
            e_idx = end if end >= 0 else len(seq)
            s = max(0, s_idx)
            e = min(len(seq), e_idx)
            return seq[int(s) : int(e)]

        # RefSeqDataProvider.get_seq expects a string for 'kind' in its implementation
        # but the protocol/stub use IdentifierType.
        kind_str = str(kind)
        return self.refseq.get_seq(ac, start, end, kind_str)

    def get_symbol_accessions(self, symbol: str, source_kind: str, target_kind: str) -> list[tuple[str, str]]:
        return self.refseq.get_symbol_accessions(symbol, source_kind, target_kind)


if os.path.exists(gff_path) and os.path.exists(fasta_path):
    print("Initializing CombinedProvider...", file=sys.stderr)
    provider = CombinedProvider(gff_path, fasta_path, protein_gpff_path)
    searcher = provider
    mapper = weaver.VariantMapper(provider=provider)
    print("Provider loaded.", file=sys.stderr)
else:
    print("Real data not found, utilizing Mock Provider.", file=sys.stderr)

    class HeuristicMockProvider:
        def __init__(self) -> None:
            self.seq_cache: dict[str, str] = {}

        def get_identifier_type(self, identifier: str) -> "weaver.IdentifierType":
            if ":" in identifier:
                ac = identifier.split(":")[0]
                if ac.startswith(("NM_", "XM_")):
                    return weaver.IdentifierType.TranscriptAccession
                if ac.startswith(("NP_", "XP_")):
                    return weaver.IdentifierType.ProteinAccession
                if ac.startswith("NC_"):
                    return weaver.IdentifierType.GenomicAccession
            return weaver.IdentifierType.ProteinAccession

        def get_transcript(self, _transcript_ac: str, _reference_ac: str | None) -> typing.Any:
            return {}

        def get_seq(self, _ac: str, _start: int, _end: int, _kind: "weaver.IdentifierType") -> str | None:
            return None

        def get_symbol_accessions(self, _symbol: str, _source_kind: str, _target_kind: str) -> list[tuple[str, str]]:
            return []

    class MockSearcher:
        def get_transcripts_for_region(self, _chrom: str, _start: int, _end: int) -> list[str]:
            return []

    provider = HeuristicMockProvider()
    searcher = MockSearcher()
    mapper = weaver.VariantMapper(provider=provider)


def clean_hgvs(s_raw: str) -> str:
    """Cleans and standardizes HGVS protein strings for simple comparison."""
    if not s_raw:
        return ""
    # Remove accession prefix
    s = s_raw
    if ":" in s:
        s = s.split(":")[-1]
    # Remove parentheses
    s = s.replace("(", "").replace(")", "")
    # Standardize Ter/*
    s = s.replace("Ter", "*")
    # Standardize silent variants: p.Ala465= -> p.=
    s = re.sub(r"p\.[A-Z][a-z][a-z]\d+=", "p.=", s)
    return s


def get_equivalence_level(v1_str: str, v2_str: str) -> typing.Optional["weaver.EquivalenceLevel"]:
    """Checks for biological equivalence using weaver's sparse reference projection."""
    if not weaver or not v1_str or not v2_str or v1_str.startswith("ERR") or v2_str.startswith("ERR"):
        return None
    try:
        # Seed heuristic provider removed

        # Robust parsing: ensure both have prefixes if one does
        if ":" in v1_str and ":" not in v2_str:
            ac = v1_str.split(":")[0]
            v2_str = f"{ac}:{v2_str}"
        elif ":" in v2_str and ":" not in v1_str:
            ac = v2_str.split(":")[0]
            v1_str = f"{ac}:{v1_str}"
        elif ":" not in v1_str and ":" not in v2_str:
            # Dummy prefix if none provided
            v1_str = f"AC:p.{v1_str.replace('p.', '')}"
            v2_str = f"AC:p.{v2_str.replace('p.', '')}"

        v1 = weaver.parse(v1_str)
        v2 = weaver.parse(v2_str)

        if el := mapper.equivalent_level(v1, v2, searcher):
            return el

    except Exception:
        return None
    return None


def check_consistency(v_nuc_str: str, v_prot_str: str, ref_prot_str: str) -> bool:
    """Checks if the nucleotide variant is biologically equivalent to the protein variant."""
    if not weaver or not v_nuc_str or not v_prot_str or v_prot_str.startswith("ERR"):
        return False

    try:
        # Normalize inputs
        v_prot_str = v_prot_str.replace("(", "").replace(")", "")

        prot_ac = None
        if ":" in ref_prot_str:
            prot_ac = ref_prot_str.split(":")[0]

        v_nuc = weaver.parse(v_nuc_str)

        # 1. Try recalculating p from c
        p_calc = None
        try:
            p_calc = mapper.c_to_p(v_nuc, prot_ac)
        except Exception:
            # Failed to calculate
            return False

        # 2. Check equivalence or string match
        try:
            parse_str = v_prot_str
            if ":" not in parse_str and prot_ac:
                parse_str = f"{prot_ac}:{parse_str}"

            v_prot = weaver.parse(parse_str)
            if mapper.equivalent(p_calc, v_prot, searcher):
                return True
        except (ValueError, RuntimeError, AttributeError):
            pass

        # Fallback to loose string comparison
        if mapper is not None:
            p_calc = mapper.c_to_p(v_nuc)
            return clean_hgvs(str(p_calc)) == clean_hgvs(v_prot_str)
        return False

    except (ValueError, RuntimeError, AttributeError):
        return False


def classify(row: dict[str, str]) -> str:  # noqa: C901, PLR0912
    rs_p = row["rs_p"]
    ref_p = row["ref_p"]
    gt_p = row["variant_prot"]
    var_nuc = row["variant_nuc"]

    c_rs = clean_hgvs(rs_p)
    c_gt = clean_hgvs(gt_p)
    c_ref = clean_hgvs(ref_p)

    category = "Other Mismatch"

    if c_rs != "" and c_rs == c_gt:
        is_consistent = check_consistency(var_nuc, rs_p, gt_p)
        category = f"ClinVar Match ({'Consistent' if is_consistent else 'Inconsistent'})"
    elif (el := get_equivalence_level(rs_p, gt_p)) and el in [
        weaver.EquivalenceLevel.Identity,
        weaver.EquivalenceLevel.Analogous,
    ]:
        el_name = repr(el).split(".")[-1]
        is_consistent = check_consistency(var_nuc, rs_p, gt_p)
        category = f"Biological Equivalence ({el_name}) ({'Consistent' if is_consistent else 'Inconsistent'}) - ClinVar"
    elif c_rs != "" and c_rs == c_ref:
        # It matches reference string, so it's a Parity Match "on paper"
        # Check if it's biologically Analogous to the GT (ClinVar)
        is_analogous_gt = False
        try:
            if mapper is not None:
                # Need to parse rs_p and gt_p to pass to equivalent_level
                # This block seems to be a copy-paste error from check_consistency
                # Reverting to original logic but ensuring mapper is not None if get_equivalence_level is called
                el = get_equivalence_level(rs_p, gt_p)
                if el == weaver.EquivalenceLevel.Analogous:
                    is_analogous_gt = True
        except (ValueError, RuntimeError, AttributeError):
            pass

        if is_analogous_gt:
            # Check consistency with NUC to confirm
            is_consistent = check_consistency(var_nuc, rs_p, gt_p)
            category = f"Parity Match (Analogous to GT) ({'Consistent' if is_consistent else 'Inconsistent'})"
        else:
            category = "Parity Match (Not Analogous to GT)"
    elif (el := get_equivalence_level(rs_p, ref_p)) and el in [
        weaver.EquivalenceLevel.Identity,
        weaver.EquivalenceLevel.Analogous,
    ]:
        el_name = repr(el).split(".")[-1]
        category = f"Biological Equivalence ({el_name}) - Parity"
    elif rs_p == "ERR" or rs_p.startswith("ERR:"):
        # Categorize known errors
        if "Transcript" in rs_p and "not found" in rs_p:
            category = "Provider Error (Missing Transcript)"
        elif "TranscriptMismatch" in rs_p:
            category = "Weaver Error: Reference sequence mismatch"
        else:
            err_str = rs_p.split(":")[-1] if ":" in rs_p else "Generic"
            err_str = err_str.lstrip()
            category = f"Weaver Error: {err_str}"
    elif "fs" in rs_p or "fs" in gt_p:
        category = "Frameshift Difference"
    elif "delins" in rs_p or "delins" in gt_p:
        category = "AA Mismatch (delins)"
    elif "[" in row["variant_nuc"] or "dup" in row["variant_nuc"]:
        category = "Repeat/Duplication Mismatch"
    elif "*" in rs_p or "*" in gt_p:
        category = "Nonsense Notation Mismatch"

    return category


def main() -> None:  # noqa: C901
    min_args = 2
    if len(sys.argv) < min_args:
        print("Usage: classify_failures.py <results_tsv>")
        sys.exit(1)

    stats: dict[str, int] = {}
    total = 0
    mismatches: dict[str, list[dict[str, str]]] = collections.defaultdict(list)
    success_count = 0

    max_samples = 100

    with open(sys.argv[1]) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            total += 1
            cat = classify(row)
            stats[cat] = stats.get(cat, 0) + 1
            if ("Match" in cat and "ClinVar" in cat) or ("Biological Equivalence" in cat and "ClinVar" in cat):
                success_count += 1
            elif len(mismatches[cat]) < max_samples:
                mismatches[cat].append(row)

    print(f"Total variants: {total}")
    print(f"Total Successes (Exact + Biological): {success_count} ({(success_count / total) * 100:.2f}%)")
    print(f"Total Failures:  {total - success_count} ({((total - success_count) / total) * 100:.2f}%)")

    print("\nFailure breakdown (excluding successful matches):")
    for cat, count in sorted(stats.items(), key=lambda x: x[1], reverse=True):
        if ("Match" in cat and "ClinVar" in cat) or ("Biological Equivalence" in cat and "ClinVar" in cat):
            continue
        print(f"  {cat}: {count} ({(count / (total - success_count)) * 100:.2f}% of failures)")

    print("\nSuccess breakdown:")
    for cat, count in sorted(stats.items(), key=lambda x: x[1], reverse=True):
        if ("Match" in cat and "ClinVar" in cat) or ("Biological Equivalence" in cat and "ClinVar" in cat):
            print(f"  {cat}: {count} ({(count / total) * 100:.2f}%)")

    if mismatches:
        print("\nSample Mismatches:")
        for cat, rows in sorted(mismatches.items(), key=lambda x: -len(x[1])):
            if (
                ("Match" in cat and "ClinVar" in cat)
                or ("Biological Equivalence" in cat and "ClinVar" in cat)
                or cat == "Provider Error (Missing Transcript)"
            ):
                continue
            print(f"CAT: {cat}")
            for row in rows:
                print("-" * 10)
                print(f"NUC: {row['variant_nuc']}")
                print(f"GT:  {row['variant_prot']}")
                print(f"W:   {row['rs_p']}")
                print(f"R:   {row['ref_p']}")
            print("-" * 40)


if __name__ == "__main__":
    main()
