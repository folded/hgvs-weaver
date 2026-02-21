use hgvs_weaver::*;

#[test]
fn test_spec_summary_variants() {
    let variants = vec![
        // DNA
        "NC_000001.11:g.1234=",
        "NC_000001.11:g.1234_1235insACGT",
        "NC_000001.11:g.1234_2345=",
        "NC_000001.11:g.1234_2345del",
        "NC_000001.11:g.1234_2345dup",
        "NC_000001.11:g.1234_2345inv",
        "NC_000001.11:g.1234del",
        "NC_000001.11:g.1234dup",
        "NC_000001.11:g.123_129delinsAC",
        "NC_000001.11:g.123delinsAC",
        "NC_000023.10:g.33038255C>A",
        "NG_012232.1(NM_004006.2):c.93+1G>T",
        // RNA
        "NM_004006.3:r.123_124insauc",
        "NM_004006.3:r.123_127del",
        "NM_004006.3:r.123_127delinsag",
        "NM_004006.3:r.123_345dup",
        "NM_004006.3:r.123_345inv",
        "NM_004006.3:r.123c>g",
        // Protein
        "NP_003070.3:p.Glu125_Ala132delinsGlyLeuHisArgPheIleValLeu",
        "NP_003997.1:p.Trp24Cys",
        "NP_003997.1:p.Trp24Ter",
        "NP_003997.2:p.Lys23_Val25del",
        "NP_003997.2:p.Lys23_Val25dup",
        "NP_003997.2:p.Val7del",
        "NP_003997.2:p.Val7dup",
        "NP_004371.2:p.Asn47delinsSerSerTer",
        "NP_0123456.1:p.Arg97fs",
    ];

    for v in variants {
        let result = parse_hgvs_variant(v);
        assert!(
            result.is_ok(),
            "Failed to parse spec example: {} - {:?}",
            v,
            result.err()
        );
        let parsed = result.unwrap();

        let rt = parsed.to_string();
        assert_eq!(rt, v, "Round trip failed for {}", v);
    }
}

#[test]
fn test_spec_summary_variants_normalized() {
    let variants = vec![
        ("NP_003997.1:p.W24*", "NP_003997.1:p.Trp24Ter"),
        ("NP_003997.1:p.(Trp24Cys)", "NP_003997.1:p.(Trp24Cys)"),
        (
            "NP_0123456.1:p.Arg97ProfsTer23",
            "NP_0123456.1:p.Arg97ProfsTer23",
        ),
    ];

    for (v, expected) in variants {
        let result = parse_hgvs_variant(v);
        assert!(result.is_ok(), "Failed to parse: {}", v);
        assert_eq!(result.unwrap().to_string(), expected);
    }
}

#[test]
fn test_gene_formatting() {
    let input = "NM_01234.5(BOGUS):c.65A>C";
    let v = parse_hgvs_variant(input).unwrap();
    assert_eq!(v.to_string(), input);

    // Test removal of gene name (structural check)
    if let SequenceVariant::Coding(mut cv) = v {
        cv.gene = None;
        assert_eq!(cv.to_string(), "NM_01234.5:c.65A>C");
    } else {
        panic!("Expected coding variant");
    }
}

#[test]
fn test_uncertain_intervals() {
    let variants = vec![
        "NC_000005.9:g.(90136803_90144453)_(90159675_90261231)dup",
        "NC_000009.11:g.(0_108337304)_(108337428_0)del",
        "NC_000005.9:g.(90136803_90159675)dup",
        "NC_000005.9:g.(90136803)_(90159675)dup",
    ];

    for v in variants {
        let parsed = parse_hgvs_variant(v)
            .unwrap_or_else(|_| panic!("Failed to parse uncertain interval: {}", v));
        assert_eq!(parsed.to_string(), v);
    }
}
