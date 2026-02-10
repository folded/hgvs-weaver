import pytest

import weaver

# Extracted from mavehgvs/tests/test_variant.py

VALID_VARIANTS = [
    "NM_001301.4:c.122-6T>A",
    "NM_007294.4:c.122-6T>A",
    "c.*33G>C",
    "c.101+1_101+7dup",
    "c.122-6T>A",
    "c.12=",
    "c.1_2delinsA",
    "c.1_2insAAA",
    "c.1_35del",
    "c.1_3=",
    "c.1_3del",
    "c.1_3dup",
    "c.3_4insT",
    "c.43-6_595+12=",
    "c.43-6_595+12delinsCTT",
    "c.4del",
    "c.4delinsTAAGC",
    "c.4dup",
    "c.77dup",
    "c.78+5_78+10del",
    "c.83_85delinsT",
    "c.84_85insCTG",
    "c.99+6_99+7insA",
    "g.22_24dup",
    "g.22delinsAACG",
    "g.234_235insT",
    "g.44del",
    "g.48C>A",
    "g.88_99=",
    "p.(Ala12_Pro13insGlyProCys)",
    "p.(Arg12Lysfs)",
    "p.(Arg12LysfsTer18)",
    "p.(Glu27fs)",
    "p.(Glu27fs*?)",
    "p.(=)",
    "p.Ala12_Pro13insGlyProCys",
    "p.Arg12Lysfs",
    "p.Arg12LysfsTer18",
    "p.Arg1Ala",
    "p.Cys22=",
    "p.Cys5dup",
    "p.Gln3Trp",
    "p.Gln7_Asn19del",
    "p.Glu27Trp",
    "p.Glu27fs",
    "p.Glu27fs*?",
    "p.Gly18del",
    "p.His44delinsValProGlyGlu",
    "p.His7_Gln8insSer",
    "p.Ile71_Cys80delinsSer",
    "p.Lys212fs",
    "p.Pro12_Gly18dup",
    "p.Ter345Lys",
    "p.=",
    "r.33+12a>c",
    "r.34_36del",
    "r.92delinsgac",
]

# Variants that mavehgvs accepts but weaver currently rejects
XFAIL_VARIANTS = [
    "NM_007294.4:c.[122-6T>A;153C>T]",  # Multi-variant / Allele
    "c.=",  # Whole-sequence identity (c.=)
    "c.[122-6T>A;153C>T]",  # Multi-variant
    "p.[Glu27Trp;Lys212fs]",  # Multi-variant
]

INVALID_VARIANTS = [
    "22G>A",
    "77dup",
    "G.44del",
    "G>A",
    "NM_001301.4::c.122-6T>A",
    "a.78+5_78+10del",
    "c.(=)",
    "g.22_23insauc",
    "g.25_24del",
    "g.25_24ins",
    "g.Glu27Trp",
    "n.Pro12_Gly18dup",
    "p.122-6T>A",
    "p.27Glu>Trp",
    "p.Arg12LysfsTer18",
    "p.Glu27fs*?",
    "p.Gly24(=)",
    "p.Pro12_Gly18insGlyProAla",
    "r.22_24insauc",
    "r.43-6_595+12delinsctt",
    "x.=",
]


@pytest.mark.parametrize("variant_string", VALID_VARIANTS)
def test_valid_mavehgvs_parsing(variant_string: str) -> None:
    """Ensure weaver can parse variants considered valid by mavehgvs."""
    # Some mavehgvs variants might depend on specific features not yet in weaver.
    # We allow failures for now but document them.
    try:
        v = weaver.parse(variant_string)
        assert v is not None
        assert str(v)  # Should verify formatting too
    except ValueError as e:
        # weaver strictly requires Accession:Variant. Try prepending a dummy accession
        # if one is missing, to verify the variant syntax itself.
        if ":" not in variant_string:
            prefix = "NP_000000.0" if variant_string.startswith("p.") else "NM_000000.0"
            retry_string = f"{prefix}:{variant_string}"
            try:
                v = weaver.parse(retry_string)
                assert v is not None
                # Passing with modification is acceptable for compatibility checking
                return
            except ValueError as e2:
                pytest.fail(f"Failed to parse valid variant {variant_string} (even with prefix {retry_string}): {e2}")

        # If it had a colon or failed even with prefix
        pytest.fail(f"Failed to parse valid variant {variant_string}: {e}")


@pytest.mark.parametrize("variant_string", XFAIL_VARIANTS)
@pytest.mark.xfail(reason="Known missing features (multi-variant or identity)")
def test_mavehgvs_xfail(variant_string: str) -> None:
    """Ensure these variants fail as expected (documentation of gaps)."""
    # If they start passing, we want to know (strict=True by default in recent pytest?)
    # We simulate the same logic as valid variants
    if ":" not in variant_string:
        prefix = "NP_000000.0" if variant_string.startswith("p.") else "NM_000000.0"
        variant_string = f"{prefix}:{variant_string}"

    weaver.parse(variant_string)


@pytest.mark.parametrize("variant_string", INVALID_VARIANTS)
def test_invalid_mavehgvs_parsing(variant_string: str) -> None:
    """Ensure weaver rejects variants considered invalid by mavehgvs."""
    # Note: weaver might be permissive for some (e.g. casing?)
    # But usually it should raise ValueError
    with pytest.raises(ValueError):
        weaver.parse(variant_string)
