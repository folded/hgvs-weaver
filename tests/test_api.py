"""Tests for the weaver Python API."""

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

    def get_seq(self, _ac: str, start: int, end: int, _kind: str) -> str:
        """Returns a mock sequence."""
        # Mock sequence: just enough to handle translation or normalization
        # Index 10 is c.1.
        # ATG GGG CCC AAA ...
        return ("A" * 10 + "ATGGGGCCCAAA" + "A" * 100)[start:end]

    def get_symbol_accessions(self, symbol: str, _s: str, t: str) -> list[tuple[weaver.IdentifierType, str]]:
        """Maps mock symbols."""
        if t == "p":
            return [(weaver.IdentifierType.ProteinAccession, "NP_TEST.1")]
        return [(weaver.IdentifierType.GeneSymbol, symbol)]


def test_parsing() -> None:
    """Tests variant parsing."""
    v = weaver.parse("NM_0001.1:c.123A>G")
    assert v.ac == "NM_0001.1"
    assert v.coordinate_type == "c"
    assert v.format() == "NM_0001.1:c.123A>G"


def test_normalization() -> None:
    """Tests variant normalization (3' shifting)."""
    provider = MockProvider()
    mapper = weaver.VariantMapper(provider)

    # c.4_5del (GG) in GGG
    v = weaver.parse("NM_TEST:c.4_5del")
    v_norm = mapper.normalize_variant(v)
    # Should shift to 5_6del
    assert v_norm.format() == "NM_TEST:c.5_6delGG"


def test_c_to_p() -> None:
    """Tests c. to p. projection."""
    provider = MockProvider()
    mapper = weaver.VariantMapper(provider)

    # c.1A>G -> p.Met1Val
    v_c = weaver.parse("NM_TEST:c.1A>G")
    v_p = mapper.c_to_p(v_c)
    assert "p.(Met1Val)" in v_p.format()


def test_c_to_g() -> None:
    """Tests c. to g. mapping."""
    provider = MockProvider()
    mapper = weaver.VariantMapper(provider)

    # c.1A>G. cds_start=10. Index 10.
    # reference_start=1000 (index 0).
    # Index 10 -> g.1011
    v_c = weaver.parse("NM_TEST:c.1A>G")
    v_g = mapper.c_to_g(v_c, "NC_TEST.1")
    assert v_g.format() == "NC_TEST.1:g.1011A>G"


class MockMinusProvider:
    """Mock DataProvider for minus strand testing."""

    def get_transcript(self, ac: str, _ref: str | None) -> dict[str, typing.Any]:
        """Returns a mock transcript model on the minus strand."""
        return {
            "ac": ac,
            "gene": "TEST",
            "cds_start_index": 10,
            "cds_end_index": 20,
            "strand": -1,
            "reference_accession": "NC_TEST.1",
            "reference_alignment_method": "splign",
            "exons": [
                {
                    "transcript_start": 0,
                    "transcript_end": 100,
                    "reference_start": 1000,
                    "reference_end": 1100,
                    "alt_strand": -1,
                    "cigar": "100=",
                },
            ],
        }

    def get_seq(self, _ac: str, start: int, end: int, _kind: str) -> str:
        """Returns a mock sequence."""
        return ("A" * 1000)[start:end]

    def get_symbol_accessions(self, symbol: str, _s: str, _t: str) -> list[tuple[weaver.IdentifierType, str]]:
        """Maps mock symbols."""
        return [(weaver.IdentifierType.GeneSymbol, symbol)]


def test_minus_strand() -> None:
    """Tests mapping on the minus strand."""
    provider = MockMinusProvider()
    mapper = weaver.VariantMapper(provider)

    v_c = weaver.parse("NM_TEST:c.1A>G")
    v_g = mapper.c_to_g(v_c, "NC_TEST.1")
    # Coordinate 1090 is HgvsGenomicPos(1091).
    assert v_g.format() == "NC_TEST.1:g.1091T>C"


def test_intronic() -> None:
    """Tests mapping of intronic variants."""
    provider = MockProvider()
    mapper = weaver.VariantMapper(provider)

    # c.10+1A>G.
    v_c = weaver.parse("NM_TEST:c.10+1A>G")
    v_g = mapper.c_to_g(v_c, "NC_TEST.1")
    assert v_g.format() == "NC_TEST.1:g.1021A>G"
