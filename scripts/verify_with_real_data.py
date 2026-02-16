import weaver
import weaver.cli.provider as p


def main() -> None:
    dp = p.RefSeqDataProvider(
        gff_path="GRCh38_latest_genomic.gff.gz",
        fasta_path="GCF_000001405.40_GRCh38.p14_genomic.fna",
    )
    mapper = weaver.VariantMapper(dp)

    variants = [
        ("NM_000038.6:c.1972_1973delinsAT", "p.Glu658Met"),
        ("NM_000527.5:c.669_683delinsAACTGCGGTAAACTGCGGTAAACT", "p.Asp224_Glu228delinsThrAlaValAsnCysGlyLysLeu"),
        ("NM_000478.6:c.876_882delinsT", "p.Gly293_Asp294del"),
        ("NM_001122606.1:c.763_768delinsTGAAGT", "p.Asn255_Pro256delinsTer"),
        ("NM_000465.4:c.66_70delinsTGCGT", "p.Pro24Ser"),
    ]

    print(f"{'Variant':<60} | {'Truth':<50} | {'Weaver (No Norm)':<50} | {'Weaver (With Norm)'}")
    print("-" * 220)

    for hgvs_c, truth_p in variants:
        try:
            var = weaver.parse(hgvs_c)
            p_no_norm = mapper.c_to_p(var)

            norm_var = mapper.normalize_variant(var)
            p_with_norm = mapper.c_to_p(norm_var)

            print(f"{hgvs_c:<60} | {truth_p:<50} | {p_no_norm!s:<50} | {p_with_norm}")
        except OSError as e:
            print(f"{hgvs_c:<60} | {truth_p:<50} | ERROR: {e}")


if __name__ == "__main__":
    main()
