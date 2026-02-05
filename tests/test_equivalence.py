"""Tests for the variant equivalence functionality."""

import typing

import weaver


class MockProvider:
    """Mock DataProvider for testing."""

    def get_transcript(self, ac: str, _ref: str | None) -> dict[str, typing.Any]:
        """Returns a mock transcript model."""
        return {
            "ac": ac,
            "gene": "TEST",
            "cds_start_index": 10,
            "cds_end_index": 20,
            "strand": 1,
            "reference_accession": "NC_TEST.1",
            "reference_alignment_method": "splign",
            "exons": [
                {
                    "transcript_start": 0,
                    "transcript_end": 100,
                    "reference_start": 1000,
                    "reference_end": 1100,
                    "alt_strand": 1,
                    "cigar": "100=",
                },
            ],
        }

    def get_seq(self, _ac: str, start: int, end: int, _kind: str | weaver.IdentifierType) -> str:
        """Returns a mock sequence."""
        # Index 10 is c.1.
        full_seq = "A" * 10 + "ATGGGGCCCAAA" + "A" * 2000
        return full_seq[start:end]

    def get_symbol_accessions(self, symbol: str, _s: str, t: str) -> list[tuple[typing.Any, str]]:
        """Maps mock symbols."""
        if symbol == "BRAF":
            if t == "p":
                # Returns as string
                return [("protein_accession", "NP_BRAF.1")]
            if t == "c":
                # Returns as enum
                return [(weaver.IdentifierType.TranscriptAccession, "NM_BRAF.1")]
            if t == "g":
                return [(weaver.IdentifierType.GenomicAccession, "NC_TEST.1")]
        if symbol == "NM_TEST" and t == "p":
            return [("protein_accession", "NP_TEST.1")]
        return [("gene_symbol", symbol)]

    def get_transcripts_for_region(self, _chrom: str, _start: int, _end: int) -> list[str]:
        """Returns transcripts for a region."""
        return ["NM_TEST"]

    def get_identifier_type(self, identifier: str) -> str | weaver.IdentifierType:
        """Identifies mock identifiers."""
        if identifier.startswith("NC_"):
            return weaver.IdentifierType.GenomicAccession
        if identifier.startswith("NM_"):
            return "transcript_accession"
        if identifier.startswith("NP_"):
            return "protein_accession"
        if identifier in ("BRAF", "NM_TEST") or "." not in identifier:
            return "gene_symbol"
        return "unknown"


def test_equivalence_g_vs_g() -> None:
    """Tests g. vs g. equivalence (with normalization)."""
    provider = MockProvider()
    mapper = weaver.VariantMapper(provider)

    # NC_TEST.1:g.1005_1006insA vs NC_TEST.1:g.1006_1007insA
    # Sequence at 1001..1010 is AAAAAAAAAA
    v1 = weaver.parse("NC_TEST.1:g.1005_1006insA")
    v2 = weaver.parse("NC_TEST.1:g.1006_1007insA")

    # These should be equivalent because of 3' shifting in A homopolymer
    assert mapper.equivalent(v1, v2, provider)  # type: ignore[attr-defined]


def test_equivalence_g_vs_c() -> None:
    """Tests g. vs c. equivalence."""
    provider = MockProvider()
    mapper = weaver.VariantMapper(provider)

    # c.1A>G maps to g.1011A>G
    vg = weaver.parse("NC_TEST.1:g.1011A>G")
    vc = weaver.parse("NM_TEST:c.1A>G")

    assert mapper.equivalent(vg, vc, provider)  # type: ignore[attr-defined]


def test_equivalence_c_vs_p() -> None:
    """Tests c. vs p. equivalence."""
    provider = MockProvider()
    mapper = weaver.VariantMapper(provider)

    # c.1A>G maps to p.Met1Val
    vc = weaver.parse("NM_TEST:c.1A>G")
    vp = weaver.parse("NP_TEST.1:p.Met1Val")

    assert mapper.equivalent(vc, vp, provider)  # type: ignore[attr-defined]


def test_equivalence_symbol_g() -> None:
    """Tests symbol-based genomic equivalence."""
    provider = MockProvider()
    mapper = weaver.VariantMapper(provider)

    # BRAF:g.1011A>G resolves to NC_TEST.1:g.1011A>G
    vg_symbol = weaver.parse("BRAF:g.1011A>G")
    vg_explicit = weaver.parse("NC_TEST.1:g.1011A>G")

    assert mapper.equivalent(vg_symbol, vg_explicit, provider)  # type: ignore[attr-defined]


def test_equivalence_symbol_c() -> None:
    """Tests symbol-based coding equivalence."""
    provider = MockProvider()
    mapper = weaver.VariantMapper(provider)

    # BRAF:c.1A>G resolves to NM_BRAF.1:c.1A>G
    vc_symbol = weaver.parse("BRAF:c.1A>G")
    vc_explicit = weaver.parse("NM_BRAF.1:c.1A>G")

    assert mapper.equivalent(vc_symbol, vc_explicit, provider)  # type: ignore[attr-defined]


def test_equivalence_g_vs_p() -> None:
    """Tests g. vs p. equivalence."""
    provider = MockProvider()
    mapper = weaver.VariantMapper(provider)

    # g.1011A>G -> c.1A>G -> p.Met1Val
    vg = weaver.parse("NC_TEST.1:g.1011A>G")
    vp = weaver.parse("NP_TEST.1:p.Met1Val")

    assert mapper.equivalent(vg, vp, provider)  # type: ignore[attr-defined]


def test_equivalence_gene_symbol() -> None:
    """Tests equivalence with gene symbols."""
    provider = MockProvider()
    mapper = weaver.VariantMapper(provider)

    # BRAF:p.V600E vs NP_BRAF.1:p.V600E
    v1 = weaver.parse("BRAF:p.Val600Glu")
    v2 = weaver.parse("NP_BRAF.1:p.Val600Glu")

    assert mapper.equivalent(v1, v2, provider)  # type: ignore[attr-defined]
