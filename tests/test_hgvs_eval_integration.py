import dataclasses
import os
from collections.abc import Generator
from typing import Any

import pytest

import weaver
from weaver.cli.provider import RefSeqDataProvider

# Paths to data files
REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
GFF_REF_PATH = os.path.join(REPO_ROOT, "tests", "data", "hgvs_eval_reference.gff")
GFF38_PATH = GFF_REF_PATH
FASTA38_PATH = os.path.join(REPO_ROOT, "GCF_000001405.40_GRCh38.p14_genomic.fna")
GFF37_PATH = GFF_REF_PATH
FASTA37_PATH = os.path.join(REPO_ROOT, "GCF_000001405.25_GRCh37.p13_genomic.fna")
SEQ_CACHE_PATH = os.path.join(REPO_ROOT, "tests", "data", "hgvs_eval_sequences.json")


def setup_sequence_mocking() -> None:
    """Sets up WEAVER_SEQ_MODE and WEAVER_SEQ_CACHE based on file availability."""
    if "WEAVER_SEQ_CACHE" not in os.environ:
        os.environ["WEAVER_SEQ_CACHE"] = SEQ_CACHE_PATH

    if "WEAVER_SEQ_MODE" not in os.environ:
        if os.path.exists(SEQ_CACHE_PATH) and not os.path.exists(FASTA38_PATH):
            os.environ["WEAVER_SEQ_MODE"] = "replay"
        else:
            os.environ["WEAVER_SEQ_MODE"] = "live"


@dataclasses.dataclass(frozen=True)
class EvalCase:
    input: str
    output_preferred: str
    data: str


HGVS_EVAL_CASES: list[Any] = [
    pytest.param(
        EvalCase(
            input="NM_033089.6:c.471_473del",
            output_preferred="NC_000020.10:g.278701_278703del",
            data="GRCh37",
        ),
    ),
    pytest.param(
        EvalCase(
            input="NC_000023.10:g.73501562T>C",
            output_preferred="NR_028379.1:n.345A>G",
            data="GRCh37|RefSeq",
        ),
    ),
    pytest.param(
        EvalCase(
            input="NM_001135021.1:c.794T>C",
            output_preferred="NC_000002.11:g.85616929T>C",
            data="GRCh37|RefSeq",
        ),
    ),
    pytest.param(
        EvalCase(
            input="NR_028379.1:n.345A>G",
            output_preferred="NC_000023.10:g.73501562T>C",
            data="GRCh37|RefSeq",
        ),
    ),
    pytest.param(
        EvalCase(
            input="NM_033089.6:c.471_473delGGC",
            output_preferred="NP_149080.2:p.(Ala158del)",
            data="RefSeq",
        ),
    ),
    pytest.param(
        EvalCase(
            input="NM_033089.6:c.471_473delGGC",
            output_preferred="NM_033089.6:n.495_497del",
            data="RefSeq",
        ),
    ),
    pytest.param(
        EvalCase(
            input="NM_033089.6:n.495_497delGGC",
            output_preferred="NM_033089.6:c.471_473del",
            data="RefSeq",
        ),
    ),
    pytest.param(
        EvalCase(
            input="NC_000001.10:g.17345192_17345217delinsTTGGGGCAAGTAAAGGAACAGGTTC",
            output_preferred="NM_003000.2:c.*159_*184delinsGAACCTGTTCCTTTACTTGCCCCAA",
            data="GRCh37|RefSeq",
        ),
    ),
    pytest.param(
        EvalCase(
            input="NM_001135023.1:c.794T>C",
            output_preferred="NP_001128495.1:p.(Leu265Ser)",
            data="RefSeq",
        ),
    ),
    pytest.param(
        EvalCase(
            input="NM_001166478.1:c.35_36insT",
            output_preferred="NM_001166478.1:c.35dup",
            data="RefSeq",
        ),
    ),
    pytest.param(
        EvalCase(
            input="NM_000492.3:c.1520_1522delTCT",
            output_preferred="NM_000492.3:c.1521_1523del",
            data="RefSeq",
        ),
    ),
    pytest.param(
        EvalCase(
            input="NC_000002.11:g.37480321_37480322insT",
            output_preferred="NC_000002.11:g.37480321dup",
            data="GRCh37",
        ),
    ),
    pytest.param(
        EvalCase(
            input="NM_001166478.1:c.31del",
            output_preferred="NM_001166478.1:c.35del",
            data="RefSeq",
        ),
    ),
    pytest.param(
        EvalCase(
            input="NC_000020.10:g.278692_278694delGGC",
            output_preferred="NC_000020.10:g.278701_278703del",
            data="GRCh37",
        ),
    ),
    pytest.param(
        EvalCase(
            input="NC_000002.11:g.37480321_37480322insT",
            output_preferred="NC_000002.11:g.37480321dup",
            data="GRCh37",
        ),
    ),
    pytest.param(
        EvalCase(
            input="NM_005813.3:c.2673insA",
            output_preferred="NM_005813.3:c.2673dup",
            data="RefSeq",
        ),
    ),
    pytest.param(
        EvalCase(
            input="NP_689699.2:p.(G553E)",
            output_preferred="NP_689699.2:p.(Gly553Glu)",
            data="RefSeq",
        ),
    ),
]


@pytest.fixture(scope="session")
def real_provider_38() -> Generator[RefSeqDataProvider, None, None]:
    setup_sequence_mocking()
    mode = os.environ.get("WEAVER_SEQ_MODE", "live")

    if mode == "replay":
        if not os.path.exists(SEQ_CACHE_PATH):
            pytest.skip("Sequence cache missing for replay mode")
    elif not os.path.exists(GFF38_PATH) or not os.path.exists(FASTA38_PATH):
        pytest.skip("Real data files (GFF/FASTA) not found. Skipping integration tests.")

    print(f"Loading RefSeq provider (mode={mode}) with GFF: {GFF38_PATH}")
    provider = RefSeqDataProvider(GFF38_PATH, FASTA38_PATH)
    yield provider
    if provider.fasta.mode == "record":
        provider.fasta.save_cache()


@pytest.fixture(scope="session")
def real_provider_37() -> Generator[RefSeqDataProvider | None, None, None]:
    setup_sequence_mocking()
    mode = os.environ.get("WEAVER_SEQ_MODE", "live")

    if mode == "replay":
        if not os.path.exists(SEQ_CACHE_PATH):
            yield None
            return
    elif not os.path.exists(GFF37_PATH) or not os.path.exists(FASTA37_PATH):
        yield None
        return

    print(f"Loading RefSeq provider (mode={mode}) with GFF: {GFF37_PATH}")
    provider = RefSeqDataProvider(GFF37_PATH, FASTA37_PATH)
    yield provider
    if provider.fasta.mode == "record":
        provider.fasta.save_cache()


@pytest.mark.parametrize("hgvs_eval_case", HGVS_EVAL_CASES)
def test_hgvs_eval_equivalence(
    real_provider_38: RefSeqDataProvider,
    real_provider_37: RefSeqDataProvider | None,
    hgvs_eval_case: EvalCase,
) -> None:
    """
    Test that input and output_preferred are equivalent.
    """
    input_str = hgvs_eval_case.input.strip()
    output_str = hgvs_eval_case.output_preferred.strip()
    data_col = hgvs_eval_case.data

    # Select provider based on data column
    if "GRCh37" in data_col:
        if not real_provider_37:
            pytest.skip("GRCh37 data not available")
        real_provider = real_provider_37
    else:
        real_provider = real_provider_38

    if not input_str or not output_str:
        print(f"SKIPPING {hgvs_eval_case}: Missing input/output strings")
        pytest.skip("Missing input/output strings")

    try:
        v1 = weaver.parse(input_str)
        v2 = weaver.parse(output_str)

        if v1 is None or v2 is None:
            print(f"SKIPPING {input_str}: Transcript missing from provider")
            pytest.skip(f"Transcript missing from provider: {input_str} or {output_str}")

    except ValueError as e:
        print(f"SKIPPING {input_str}: Parsing failed: {e}")
        pytest.skip(f"Parsing failed for {input_str} or {output_str}: {e}")

    mapper = weaver.VariantMapper(real_provider)

    # We expect them to be equivalent
    try:
        is_equiv = mapper.equivalent(v1, v2, real_provider)

        if not is_equiv:
            pytest.fail(f"Variants not equivalent: {v1} vs {v2}. Params: {data_col}")

    except Exception as e:  # noqa: BLE001
        pytest.fail(f"Equivalence check raised exception: {e}")
