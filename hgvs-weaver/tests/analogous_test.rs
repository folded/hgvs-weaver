use hgvs_weaver::analogous_edit::{
    apply_aa_edit_to_sparse, apply_na_edit_to_sparse, project_aa_variant, reconcile_projections,
    ResidueToken, SparseReference,
};
use hgvs_weaver::structs::{AaEdit, NaEdit};

#[test]
fn test_residue_token_unification_known() {
    let t1 = vec![ResidueToken::Known("A".into())];
    let t2 = vec![ResidueToken::Known("A".into())];
    assert!(reconcile_projections(&t1, &t2));

    let t3 = vec![ResidueToken::Known("A".into())];
    let t4 = vec![ResidueToken::Known("C".into())];
    assert!(!reconcile_projections(&t3, &t4));
}

#[test]
fn test_residue_token_unification_unknown() {
    let t1 = vec![ResidueToken::Unknown(10), ResidueToken::Unknown(11)];
    let t2 = vec![ResidueToken::Unknown(10), ResidueToken::Unknown(11)];
    assert!(reconcile_projections(&t1, &t2));

    // Cross-mapping unknown to known
    let t3 = vec![ResidueToken::Unknown(10)];
    let t4 = vec![ResidueToken::Known("G".into())];
    assert!(reconcile_projections(&t3, &t4));

    // Consistency check: unknown 10 mapped to G and then compared to A
    let t5 = vec![ResidueToken::Unknown(10), ResidueToken::Unknown(10)];
    let t6 = vec![
        ResidueToken::Known("G".into()),
        ResidueToken::Known("A".into()),
    ];
    assert!(!reconcile_projections(&t5, &t6));
}

#[test]
fn test_analogous_duplication_shift() {
    let mut sref = SparseReference::new();
    // ABCABC -> ABCABCABC
    // dup of 1-3 (ABC) -> ABC ABC ABC
    // dup of 4-6 (ABC) -> ABC ABC ABC

    sref.set(1, "A".into()).unwrap();
    sref.set(2, "B".into()).unwrap();
    sref.set(3, "C".into()).unwrap();

    let edit = NaEdit::Dup {
        ref_: None,
        uncertain: false,
    };

    // v1: dup of 1-3. Projection returns [1, 2, 3, 1, 2, 3]
    let v1_seq = apply_na_edit_to_sparse(&edit, 1, 3, &sref);

    // v2: dup of 4-6. Projection returns [4, 5, 6, 4, 5, 6]
    let v2_seq = apply_na_edit_to_sparse(&edit, 4, 6, &sref);

    // They match if 1=4, 2=5, 3=6.
    assert!(v1_seq.is_analogous_to(&v2_seq));
}

#[test]
fn test_analogous_protein_repeats_shift() {
    let mut sref = SparseReference::new();
    // TrpTrpTrp -> TrpTrpTrpTrpTrp (Repeat TrpTrp twice)

    sref.set(50, "Trp".into()).unwrap();
    sref.set(51, "Trp".into()).unwrap();

    // v1: Repeat 50-51 twice. Projection: [50, 51, 50, 51]
    let repeat_edit = AaEdit::Repeat {
        ref_: None,
        min: 2,
        max: 2,
        uncertain: false,
    };
    let v1_seq = apply_aa_edit_to_sparse(&repeat_edit, 50, 51, &sref);

    // v2: Repeat 52-53 twice. Projection: [52, 53, 52, 53]
    let v2_seq = apply_aa_edit_to_sparse(&repeat_edit, 52, 53, &sref);

    assert!(v1_seq.is_analogous_to(&v2_seq));
}

#[test]
fn test_complex_unification_aliases() {
    let v1 = vec![ResidueToken::Unknown(10), ResidueToken::Known("A".into())];
    let v2 = vec![ResidueToken::Known("G".into()), ResidueToken::Unknown(20)];

    assert!(reconcile_projections(&v1, &v2));

    let v3 = vec![ResidueToken::Unknown(10), ResidueToken::Unknown(20)];
    let v4 = vec![
        ResidueToken::Known("G".into()),
        ResidueToken::Known("A".into()),
    ];

    assert!(reconcile_projections(&v3, &v4));
}

#[test]
fn test_inconsistent_aliases() {
    let v1 = vec![ResidueToken::Unknown(10), ResidueToken::Unknown(10)];
    let v2 = vec![
        ResidueToken::Known("A".into()),
        ResidueToken::Known("G".into()),
    ];
    assert!(!reconcile_projections(&v1, &v2));
}

#[test]
fn test_complex_dup_equivalence() {
    let mut sref = SparseReference::new();
    // GT: NP_001365188.1:p.Ala201_Val202insGlyProGlyAla
    // W:  p.(Gly198_Ala201dup)

    // Model the reference sequence for Weaver's sparse context.
    sref.set(198, "Gly".into()).unwrap();
    // 199 is unknown
    // 200 is unknown
    sref.set(201, "Ala".into()).unwrap();

    // GT: ins GlyProGlyAla after 201
    let ins_edit = AaEdit::Ins {
        alt: "GlyProGlyAla".into(),
        uncertain: false,
    };
    // Use project_aa_variant to compare the sequence over a window (e.g. 198..205).
    // GT: Ins at 201_202 (between 201 and 202).
    let v1_proj = project_aa_variant(&ins_edit, 201, 202, 198, 205, &sref);

    // W: dup 198-201
    // Dup range 198..201.
    let dup_edit = AaEdit::Dup {
        ref_: None,
        uncertain: false,
    };
    // Dup range 198..201.
    let v2_proj = project_aa_variant(&dup_edit, 198, 201, 198, 205, &sref);

    // Assertion: They should be analogous due to wildcard matching
    assert!(v1_proj.is_analogous_to(&v2_proj));
}

#[test]
fn test_cascading_unification() {
    // 1(unknown) matched to A(known)
    // 2(unknown) matched to 1(unknown) -> implies 2 matched to A
    // Then check if 2 matches A
    // v1: 1, 2, A
    // v2: A, 1, 2
    // Pair 1: (1, A) -> 1=A
    // Pair 2: (2, 1) -> 2=1 -> 2=A
    // Pair 3: (A, 2) -> A=2 -> A=A -> OK
    let v1 = vec![
        ResidueToken::Unknown(1),
        ResidueToken::Unknown(2),
        ResidueToken::Known("A".into()),
    ];
    let v2 = vec![
        ResidueToken::Known("A".into()),
        ResidueToken::Unknown(1),
        ResidueToken::Unknown(2),
    ];
    assert!(reconcile_projections(&v1, &v2));

    // Negative case:
    // 1 -> A
    // 2 -> 1 (implies 2 -> A)
    // Check if 2 matches B (Should Fail)
    // v3: 1, 2, B
    // v4: A, 1, 2
    // Pair 1: (1, A) -> 1=A
    // Pair 2: (2, 1) -> 2=1 -> 2=A
    // Pair 3: (B, 2) -> B=2 -> B=A -> FAIL
    let v3 = vec![
        ResidueToken::Unknown(1),
        ResidueToken::Unknown(2),
        ResidueToken::Known("B".into()),
    ];
    let v4 = vec![
        ResidueToken::Known("A".into()),
        ResidueToken::Unknown(1),
        ResidueToken::Unknown(2),
    ];
    assert!(!reconcile_projections(&v3, &v4));
}

#[test]
fn test_clinvar_regression_tyr165ter() -> Result<(), hgvs_weaver::error::HgvsError> {
    use hgvs_weaver::data::{DataProvider, IdentifierKind, IdentifierType, Transcript};
    use hgvs_weaver::equivalence::{EquivalenceLevel, VariantEquivalence};
    use hgvs_weaver::structs::TranscriptPos;

    // Regression test for NM_001350334.2:c.495_498del (Frameshift at codon 165)
    // Weaver: p.(Tyr165Ter) -> p.Tyr165Ter
    // ClinVar: p.Ala164_Tyr165insTer
    // These should now be Analogous thanks to offset fix.

    struct LocalMockProvider;
    impl DataProvider for LocalMockProvider {
        fn get_transcript(
            &self,
            _ac: &str,
            _ref: Option<&str>,
        ) -> Result<Box<dyn Transcript>, hgvs_weaver::error::HgvsError> {
            Err(hgvs_weaver::error::HgvsError::ValidationError(
                "Not implemented".into(),
            ))
        }
        fn get_seq(
            &self,
            _ac: &str,
            start: i32,
            end: i32,
            _kind: IdentifierType,
        ) -> Result<String, hgvs_weaver::error::HgvsError> {
            let mut seq = String::new();
            // 163: Leu, 164: Ala, 165: Tyr, 166: Arg
            for i in start..end {
                match i {
                    163 => seq.push('L'),
                    164 => seq.push('A'),
                    165 => seq.push('Y'),
                    166 => seq.push('R'),
                    _ => seq.push('X'),
                }
            }
            Ok(seq)
        }
        fn get_symbol_accessions(
            &self,
            _: &str,
            _: IdentifierKind,
            _: IdentifierKind,
        ) -> Result<Vec<(IdentifierType, String)>, hgvs_weaver::error::HgvsError> {
            Ok(vec![])
        }
        fn get_identifier_type(
            &self,
            _: &str,
        ) -> Result<IdentifierType, hgvs_weaver::error::HgvsError> {
            Ok(IdentifierType::ProteinAccession)
        }
        fn c_to_g(
            &self,
            _: &str,
            _: TranscriptPos,
            _: hgvs_weaver::structs::IntronicOffset,
        ) -> Result<(String, hgvs_weaver::structs::GenomicPos), hgvs_weaver::error::HgvsError>
        {
            Err(hgvs_weaver::error::HgvsError::ValidationError(
                "Not implemented".into(),
            ))
        }
    }

    struct LocalMockSearch;
    impl hgvs_weaver::data::TranscriptSearch for LocalMockSearch {
        fn get_transcripts_for_region(
            &self,
            _: &str,
            _: i32,
            _: i32,
        ) -> Result<Vec<String>, hgvs_weaver::error::HgvsError> {
            Ok(vec![])
        }
    }

    let hdp = LocalMockProvider;
    let search = LocalMockSearch;
    let eq = VariantEquivalence::new(&hdp, &search);

    let v1 = hgvs_weaver::parse_hgvs_variant("NP_001337263.1:p.Tyr165Ter")?;
    let v2 = hgvs_weaver::parse_hgvs_variant("NP_001337263.1:p.Ala164_Tyr165insTer")?;

    let level = eq.equivalent_level(&v1, &v2)?;
    assert_eq!(level, EquivalenceLevel::Analogous);
    Ok(())
}

#[test]
fn test_analogous_protein_truncation() -> Result<(), hgvs_weaver::error::HgvsError> {
    use hgvs_weaver::data::{DataProvider, IdentifierKind, IdentifierType, Transcript};
    use hgvs_weaver::equivalence::{EquivalenceLevel, VariantEquivalence};
    use hgvs_weaver::structs::{GenomicPos, IntronicOffset, TranscriptPos};

    // User Case: p.Tyr1433_Lys1434delinsTer vs p.(Tyr1433_Val3056del)
    // Both result in truncation at 1433.
    // delinsTer -> ...Tyr1433*
    // del -> ...Tyr1433 (end of sequence)

    struct TruncationMockProvider;
    impl DataProvider for TruncationMockProvider {
        fn get_transcript(
            &self,
            _ac: &str,
            _ref: Option<&str>,
        ) -> Result<Box<dyn Transcript>, hgvs_weaver::error::HgvsError> {
            Err(hgvs_weaver::error::HgvsError::UnsupportedOperation(
                "Not needed".into(),
            ))
        }
        fn get_seq(
            &self,
            _ac: &str,
            _start: i32,
            _end: i32,
            kind: IdentifierType,
        ) -> Result<String, hgvs_weaver::error::HgvsError> {
            if kind == IdentifierType::ProteinAccession {
                // Mock a long protein sequence
                return Ok("M".repeat(3500));
            }
            Ok("".to_string())
        }
        fn get_symbol_accessions(
            &self,
            _: &str,
            _: IdentifierKind,
            _: IdentifierKind,
        ) -> Result<Vec<(IdentifierType, String)>, hgvs_weaver::error::HgvsError> {
            Ok(vec![])
        }
        fn get_identifier_type(
            &self,
            _: &str,
        ) -> Result<IdentifierType, hgvs_weaver::error::HgvsError> {
            Ok(IdentifierType::ProteinAccession)
        }
        fn c_to_g(
            &self,
            _: &str,
            _: TranscriptPos,
            _: IntronicOffset,
        ) -> Result<(String, GenomicPos), hgvs_weaver::error::HgvsError> {
            Ok(("".into(), GenomicPos(0)))
        }
    }

    struct TruncationMockSearch;
    impl hgvs_weaver::data::TranscriptSearch for TruncationMockSearch {
        fn get_transcripts_for_region(
            &self,
            _: &str,
            _: i32,
            _: i32,
        ) -> Result<Vec<String>, hgvs_weaver::error::HgvsError> {
            Ok(vec![])
        }
    }

    let hdp = TruncationMockProvider;
    let search = TruncationMockSearch;
    let eq = VariantEquivalence::new(&hdp, &search);

    // GT: NP_000042.3:p.Tyr1433_Lys1434delinsTer
    let v_gt = hgvs_weaver::parse_hgvs_variant("NP_000042.3:p.Tyr1433_Lys1434delinsTer")?;
    // W: p.(Tyr1433_Val3056del)
    let v_w = hgvs_weaver::parse_hgvs_variant("NP_000042.3:p.(Tyr1433_Val3056del)")?;

    let lvl = eq.equivalent_level(&v_gt, &v_w)?;
    assert!(matches!(
        lvl,
        EquivalenceLevel::Analogous | EquivalenceLevel::Identity
    ));
    Ok(())
}

#[test]
fn test_analogous_clinvar_tyr165ter_mismatch() -> Result<(), hgvs_weaver::error::HgvsError> {
    use hgvs_weaver::data::{DataProvider, IdentifierKind, IdentifierType, Transcript};
    use hgvs_weaver::equivalence::{EquivalenceLevel, VariantEquivalence};
    use hgvs_weaver::structs::TranscriptPos;

    struct MockReproProvider;
    impl DataProvider for MockReproProvider {
        fn get_transcript(
            &self,
            _ac: &str,
            _ref: Option<&str>,
        ) -> Result<Box<dyn Transcript>, hgvs_weaver::error::HgvsError> {
            Err(hgvs_weaver::error::HgvsError::ValidationError(
                "Not implemented".into(),
            ))
        }
        fn get_seq(
            &self,
            _ac: &str,
            start: i32,
            end: i32,
            _kind: IdentifierType,
        ) -> Result<String, hgvs_weaver::error::HgvsError> {
            let mut seq = String::new();
            for i in start..end {
                match i {
                    163 => seq.push('L'),
                    164 => seq.push('A'),
                    165 => seq.push('Y'),
                    166 => seq.push('R'),
                    _ => seq.push('X'),
                }
            }
            Ok(seq)
        }
        fn get_symbol_accessions(
            &self,
            _: &str,
            _: IdentifierKind,
            _: IdentifierKind,
        ) -> Result<Vec<(IdentifierType, String)>, hgvs_weaver::error::HgvsError> {
            Ok(vec![])
        }
        fn get_identifier_type(
            &self,
            _: &str,
        ) -> Result<IdentifierType, hgvs_weaver::error::HgvsError> {
            Ok(IdentifierType::ProteinAccession)
        }
        fn c_to_g(
            &self,
            _: &str,
            _: TranscriptPos,
            _: hgvs_weaver::structs::IntronicOffset,
        ) -> Result<(String, hgvs_weaver::structs::GenomicPos), hgvs_weaver::error::HgvsError>
        {
            Err(hgvs_weaver::error::HgvsError::ValidationError(
                "Not implemented".into(),
            ))
        }
    }

    struct MockSearch;
    impl hgvs_weaver::data::TranscriptSearch for MockSearch {
        fn get_transcripts_for_region(
            &self,
            _: &str,
            _: i32,
            _: i32,
        ) -> Result<Vec<String>, hgvs_weaver::error::HgvsError> {
            Ok(vec![])
        }
    }

    let hdp = MockReproProvider;
    let search = MockSearch;
    let eq = VariantEquivalence::new(&hdp, &search);

    let v1 = hgvs_weaver::parse_hgvs_variant("NP_001337263.1:p.Tyr165Ter")?;
    let v2 = hgvs_weaver::parse_hgvs_variant("NP_001337263.1:p.Ala164_Tyr165insTer")?;

    let level = eq.equivalent_level(&v1, &v2)?;
    assert_eq!(level, EquivalenceLevel::Analogous);
    Ok(())
}

#[test]
fn test_analogous_repeat_equivalence() -> Result<(), hgvs_weaver::error::HgvsError> {
    use hgvs_weaver::data::Transcript;
    use hgvs_weaver::equivalence::VariantEquivalence;
    use hgvs_weaver::structs::{GenomicPos, IntronicOffset, TranscriptPos};
    use hgvs_weaver::{parse_hgvs_variant, DataProvider, IdentifierKind};

    struct MockRepeatEqProvider;
    impl DataProvider for MockRepeatEqProvider {
        fn get_seq(
            &self,
            _ac: &str,
            _start: i32,
            _end: i32,
            _kind: hgvs_weaver::data::IdentifierType,
        ) -> Result<String, hgvs_weaver::error::HgvsError> {
            Ok("X".repeat(489) + "PRS" + &"X".repeat(100))
        }
        fn get_transcript(
            &self,
            _ac: &str,
            _gene: Option<&str>,
        ) -> Result<Box<dyn Transcript>, hgvs_weaver::error::HgvsError> {
            panic!("Not implemented")
        }
        fn get_identifier_type(
            &self,
            _ac: &str,
        ) -> Result<hgvs_weaver::data::IdentifierType, hgvs_weaver::error::HgvsError> {
            Ok(hgvs_weaver::data::IdentifierType::ProteinAccession)
        }
        fn get_symbol_accessions(
            &self,
            _symbol: &str,
            _source: IdentifierKind,
            _target: IdentifierKind,
        ) -> Result<Vec<(hgvs_weaver::data::IdentifierType, String)>, hgvs_weaver::error::HgvsError>
        {
            Ok(vec![])
        }
        fn c_to_g(
            &self,
            _ac: &str,
            _c_pos: TranscriptPos,
            _offset: IntronicOffset,
        ) -> Result<(String, GenomicPos), hgvs_weaver::error::HgvsError> {
            panic!("Not implemented")
        }
    }

    struct MockSearcher;
    impl hgvs_weaver::data::TranscriptSearch for MockSearcher {
        fn get_transcripts_for_region(
            &self,
            _chrom: &str,
            _start: i32,
            _end: i32,
        ) -> Result<Vec<String>, hgvs_weaver::error::HgvsError> {
            Ok(vec![])
        }
    }

    let provider = MockRepeatEqProvider;
    let searcher = MockSearcher;
    let eq = VariantEquivalence::new(&provider, &searcher);

    let v1 = parse_hgvs_variant("NP_001365049.1:p.490PRS[1]")?;
    let v2 = parse_hgvs_variant("NP_001365049.1:p.=")?;
    let res = eq.equivalent_level(&v1, &v2)?;
    assert!(res.is_equivalent());

    let v3 = parse_hgvs_variant("NP_001365049.1:p.490PRS[2]")?;
    let v4 = parse_hgvs_variant("NP_001365049.1:p.490_492dup")?;
    let res2 = eq.equivalent_level(&v3, &v4)?;
    assert!(res2.is_equivalent());
    Ok(())
}

#[test]
fn test_analogous_fs_wildcard_unification() -> Result<(), hgvs_weaver::error::HgvsError> {
    use hgvs_weaver::data::Transcript;
    use hgvs_weaver::equivalence::VariantEquivalence;
    use hgvs_weaver::structs::{GenomicPos, IntronicOffset, TranscriptPos};
    use hgvs_weaver::{parse_hgvs_variant, DataProvider, IdentifierKind};

    struct MockWildcardProvider;
    impl DataProvider for MockWildcardProvider {
        fn get_seq(
            &self,
            _ac: &str,
            _start: i32,
            _end: i32,
            _kind: hgvs_weaver::data::IdentifierType,
        ) -> Result<String, hgvs_weaver::error::HgvsError> {
            Ok("X".repeat(96) + "R" + &"X".repeat(100))
        }
        fn get_transcript(
            &self,
            _ac: &str,
            _gene: Option<&str>,
        ) -> Result<Box<dyn Transcript>, hgvs_weaver::error::HgvsError> {
            panic!("Not implemented")
        }
        fn get_identifier_type(
            &self,
            _ac: &str,
        ) -> Result<hgvs_weaver::data::IdentifierType, hgvs_weaver::error::HgvsError> {
            Ok(hgvs_weaver::data::IdentifierType::ProteinAccession)
        }
        fn get_symbol_accessions(
            &self,
            _symbol: &str,
            _source: IdentifierKind,
            _target: IdentifierKind,
        ) -> Result<Vec<(hgvs_weaver::data::IdentifierType, String)>, hgvs_weaver::error::HgvsError>
        {
            Ok(vec![])
        }
        fn c_to_g(
            &self,
            _ac: &str,
            _c_pos: TranscriptPos,
            _offset: IntronicOffset,
        ) -> Result<(String, GenomicPos), hgvs_weaver::error::HgvsError> {
            panic!("Not implemented")
        }
    }

    struct MockSearcher;
    impl hgvs_weaver::data::TranscriptSearch for MockSearcher {
        fn get_transcripts_for_region(
            &self,
            _chrom: &str,
            _start: i32,
            _end: i32,
        ) -> Result<Vec<String>, hgvs_weaver::error::HgvsError> {
            Ok(vec![])
        }
    }

    let provider = MockWildcardProvider;
    let searcher = MockSearcher;
    let eq = VariantEquivalence::new(&provider, &searcher);

    let v1 = parse_hgvs_variant("NP_001.1:p.Arg97ProfsTer4")?;
    let v2 = parse_hgvs_variant("NP_001.1:p.Arg97_Arg97delinsProAlaValTer")?;
    let res = eq.equivalent_level(&v1, &v2)?;
    assert!(res.is_equivalent());

    let v3 = parse_hgvs_variant("NP_001.1:p.Arg97_Arg97delinsProAlaValLeuTer")?;
    let res3 = eq.equivalent_level(&v1, &v3)?;
    assert!(!res3.is_equivalent());

    let v4 = parse_hgvs_variant("NP_001.1:p.Arg97_Arg97delinsLeuAlaValTer")?;
    let res4 = eq.equivalent_level(&v1, &v4)?;
    assert!(!res4.is_equivalent());

    let v5 = parse_hgvs_variant("NP_001.1:p.Arg97_Arg97delinsProAlaValLys")?;
    let res5 = eq.equivalent_level(&v1, &v5)?;
    assert!(!res5.is_equivalent());
    Ok(())
}

#[test]
fn test_multi_unit_repeat_equivalence() -> Result<(), hgvs_weaver::error::HgvsError> {
    use hgvs_weaver::data::Transcript;
    use hgvs_weaver::equivalence::VariantEquivalence;
    use hgvs_weaver::structs::{GenomicPos, IntronicOffset, TranscriptPos};
    use hgvs_weaver::{parse_hgvs_variant, DataProvider, IdentifierKind};

    struct MockMultiRepeatProvider;
    impl DataProvider for MockMultiRepeatProvider {
        fn get_seq(
            &self,
            ac: &str,
            _start: i32,
            _end: i32,
            _kind: hgvs_weaver::data::IdentifierType,
        ) -> Result<String, hgvs_weaver::error::HgvsError> {
            if ac == "NP_001122316.1" {
                // ...pppsvsatg pgpgpgpgpg pgpgpappny s...
                // Residue 229 starts the GP repeat.
                // 228 X's + 8 GP units (16 chars) + 100 X's
                Ok("X".repeat(228) + "GPGPGPGPGPGPGPGP" + &"X".repeat(100))
            } else if ac == "NP_000067.1" {
                // ...lap apapapap apapvaapap apapapapap apapapdaap...
                // Residue 179 starts the AP repeat.
                // 178 X's + 12 AP units (24 chars) + 100 X's
                Ok("X".repeat(178) + "APAPAPAPAPAPAPAPAPAPAPAP" + &"X".repeat(100))
            } else {
                Ok("".to_string())
            }
        }
        fn get_transcript(
            &self,
            _ac: &str,
            _gene: Option<&str>,
        ) -> Result<Box<dyn Transcript>, hgvs_weaver::error::HgvsError> {
            panic!("Not implemented")
        }
        fn get_identifier_type(
            &self,
            _ac: &str,
        ) -> Result<hgvs_weaver::data::IdentifierType, hgvs_weaver::error::HgvsError> {
            Ok(hgvs_weaver::data::IdentifierType::ProteinAccession)
        }
        fn get_symbol_accessions(
            &self,
            _symbol: &str,
            _source: IdentifierKind,
            _target: IdentifierKind,
        ) -> Result<Vec<(hgvs_weaver::data::IdentifierType, String)>, hgvs_weaver::error::HgvsError>
        {
            Ok(vec![])
        }
        fn c_to_g(
            &self,
            _ac: &str,
            _c_pos: TranscriptPos,
            _offset: IntronicOffset,
        ) -> Result<(String, GenomicPos), hgvs_weaver::error::HgvsError> {
            panic!("Not implemented")
        }
    }

    struct MockSearcher;
    impl hgvs_weaver::data::TranscriptSearch for MockSearcher {
        fn get_transcripts_for_region(
            &self,
            _chrom: &str,
            _start: i32,
            _end: i32,
        ) -> Result<Vec<String>, hgvs_weaver::error::HgvsError> {
            Ok(vec![])
        }
    }

    let provider = MockMultiRepeatProvider;
    let searcher = MockSearcher;
    let eq = VariantEquivalence::new(&provider, &searcher);

    // Case 1: p.229GP[2] vs p.Gly233_Pro244del
    // 8 units -> 2 units = delete 6 units (12 residues)
    // 229GP refers to 229-230.
    // Gly233_Pro244del is 233 to 244 (12 residues).
    let v1 = parse_hgvs_variant("NP_001122316.1:p.229GP[2]")?;
    let v2 = parse_hgvs_variant("NP_001122316.1:p.Gly233_Pro244del")?;
    let res1 = eq.equivalent_level(&v1, &v2)?;
    assert!(
        res1.is_equivalent(),
        "GP repeat: {} vs {} -> {:?}",
        v1,
        v2,
        res1
    );

    // Case 2: p.179_180AP[5] vs p.Ala189_Pro194del
    // If W is 6 residues (3 units), and GT is AP[5], then initial was 8 units.
    // Ala189_Pro194del (residues 189, 190, 191, 192, 193, 194) is 6 residues.
    let v3 = parse_hgvs_variant("NP_000067.1:p.179_180AP[5]")?;
    let v4 = parse_hgvs_variant("NP_000067.1:p.Ala189_Pro194del")?;
    let res2 = eq.equivalent_level(&v3, &v4)?;
    assert!(
        res2.is_equivalent(),
        "AP repeat: {} vs {} -> {:?}",
        v3,
        v4,
        res2
    );

    Ok(())
}
