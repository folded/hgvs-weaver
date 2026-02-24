"""Tests for the variant transform functionality."""

import weaver


def test_transform_met1_to_question() -> None:
    """Tests that p.(Met1Val) transforms to p.Met1? with HgvsQuestion convention."""
    settings = weaver.VariantTransformSettings(
        start_codon=weaver.StartCodonConvention.HgvsQuestion
    )

    v = weaver.parse("NP_000051.2:p.(Met1Val)")
    transformed = v.transform(settings)

    assert transformed.format() == "NP_000051.2:p.Met1?"


def test_transform_specific_unchanged() -> None:
    """Tests that Specific convention leaves p.(Met1Val) unchanged."""
    settings = weaver.VariantTransformSettings(
        start_codon=weaver.StartCodonConvention.Specific
    )

    v = weaver.parse("NP_000051.2:p.(Met1Val)")
    transformed = v.transform(settings)

    assert transformed.format() == "NP_000051.2:p.(Met1Val)"


def test_transform_non_met1_unchanged() -> None:
    """Tests that HgvsQuestion does not affect non-start-codon variants."""
    settings = weaver.VariantTransformSettings(
        start_codon=weaver.StartCodonConvention.HgvsQuestion
    )

    v = weaver.parse("NP_000051.2:p.(Gly2Arg)")
    transformed = v.transform(settings)

    # Position 2, not Met1 — should not be affected
    assert transformed.format() == "NP_000051.2:p.(Gly2Arg)"


def test_transform_silent_met1_unchanged() -> None:
    """Tests that silent Met1 variants are not converted to p.Met1?."""
    settings = weaver.VariantTransformSettings(
        start_codon=weaver.StartCodonConvention.HgvsQuestion
    )

    v = weaver.parse("NP_000051.2:p.(Met1=)")
    transformed = v.transform(settings)

    # Silent change — should not be converted to p.Met1?
    assert transformed.format() == "NP_000051.2:p.(Met1=)"


def test_transform_non_protein_unchanged() -> None:
    """Tests that non-protein variants pass through transform unchanged."""
    settings = weaver.VariantTransformSettings(
        start_codon=weaver.StartCodonConvention.HgvsQuestion
    )

    v = weaver.parse("NM_000051.3:c.1A>T")
    transformed = v.transform(settings)

    assert transformed.format() == "NM_000051.3:c.1A>T"


def test_transform_default_settings() -> None:
    """Tests that default settings (Specific) leave all variants unchanged."""
    settings = weaver.VariantTransformSettings()

    v = weaver.parse("NP_000051.2:p.(Met1Val)")
    transformed = v.transform(settings)

    assert transformed.format() == "NP_000051.2:p.(Met1Val)"


def test_transform_met1_leu_to_question() -> None:
    """Tests p.(Met1Leu) transforms to p.Met1? under HgvsQuestion."""
    settings = weaver.VariantTransformSettings(
        start_codon=weaver.StartCodonConvention.HgvsQuestion
    )

    v = weaver.parse("NP_000051.2:p.(Met1Leu)")
    transformed = v.transform(settings)

    assert transformed.format() == "NP_000051.2:p.Met1?"
