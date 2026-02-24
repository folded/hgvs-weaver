import weaver
from weaver.cli.provider import RefSeqDataProvider


def test():
    gff_path = "GCF_000001405.25_GRCh37.p13_genomic.gff.gz"
    fasta_path = "GCF_000001405.25_GRCh37.p13_genomic.fna"

    hdp = RefSeqDataProvider(gff_path, fasta_path)
    mapper = weaver.VariantMapper(hdp)

    # c.6833-1G>A is intronic
    c_str = "NM_013227.3:c.6833-1G>A"

    try:
        v_c = weaver.parse(c_str)
        print(f"Parsed c: {v_c}")

        # This should now return p.? instead of failing
        v_p = mapper.c_to_p(v_c)
        print(f"Result p: {v_p}")

        expected = "NP_037359.3:p.?"
        if str(v_p) == expected:
            print("SUCCESS: Correctly returned p.?")
        else:
            print(f"FAILURE: Expected {expected}, got {v_p}")

    except Exception as e:
        print(f"ERROR: {e}")


if __name__ == "__main__":
    test()
