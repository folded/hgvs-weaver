import weaver
from weaver.cli.provider import RefSeqDataProvider


def test_intronic_cp() -> None:
    gff_path = "GCF_000001405.25_GRCh37.p13_genomic.gff.gz"
    fasta_path = "GCF_000001405.25_GRCh37.p13_genomic.fna"

    hdp = RefSeqDataProvider(gff_path, fasta_path)
    mapper = weaver.VariantMapper(hdp)

    # c.6833-1G>A is intronic
    c_str = "NM_013227.3:c.6833-1G>A"

    v_c = weaver.parse(c_str)
    assert v_c is not None

    # This should now return p.? instead of failing
    v_p = mapper.c_to_p(v_c)

    expected = "NP_037359.3:p.?"
    assert str(v_p) == expected
