use hgvs_weaver::coords::{GenomicPos, IntronicOffset, SequenceVariant, TranscriptPos};
use hgvs_weaver::data::{
    DataProvider, ExonData, IdentifierKind, IdentifierType, Transcript, TranscriptData,
};
use hgvs_weaver::error::HgvsError;
use hgvs_weaver::mapper::VariantMapper;

struct RegressionProvider;
impl DataProvider for RegressionProvider {
    fn get_transcript(
        &self,
        ac: &str,
        _ref_ac: Option<&str>,
    ) -> Result<Box<dyn Transcript>, HgvsError> {
        let (cds_start, _protein_ac) = match ac {
            "NM_153046.3" => (0, "NP_694591.2"),
            "NM_058216.3" => (0, "NP_478123.1"),
            _ => (0, "NP_UNKNOWN"),
        };

        Ok(Box::new(TranscriptData {
            ac: ac.to_string(),
            gene: "TEST".to_string(),
            cds_start_index: Some(TranscriptPos(cds_start)),
            cds_end_index: Some(TranscriptPos(2000)),
            strand: 1,
            reference_accession: "NC_TEST".to_string(),
            exons: vec![ExonData {
                transcript_start: TranscriptPos(0),
                transcript_end: TranscriptPos(2000),
                reference_start: GenomicPos(1000),
                reference_end: GenomicPos(3000),
                alt_strand: 1,
                cigar: "2000M".to_string(),
            }],
        }))
    }
    fn get_seq(
        &self,
        _ac: &str,
        start: i32,
        end: i32,
        _kind: IdentifierType,
    ) -> Result<String, HgvsError> {
        let mut seq = "ACGC".repeat(1000).into_bytes();
        // For NM_058216.3:c.692_694delinsAA
        // Ser231: 691, 692, 693
        // If 691 is T, 692-693 replaced by AA -> TAA (Stop)
        if seq.len() > 690 {
            seq[690] = b'T';
        }
        let s = start as usize;
        let e = if end < 0 { seq.len() } else { end as usize };
        if s > seq.len() {
            return Ok("".to_string());
        }
        let actual_e = e.min(seq.len());
        Ok(String::from_utf8_lossy(&seq[s..actual_e]).to_string())
    }
    fn get_symbol_accessions(
        &self,
        ac: &str,
        _f: IdentifierKind,
        t: IdentifierKind,
    ) -> Result<Vec<(IdentifierType, String)>, HgvsError> {
        if t == IdentifierKind::Protein {
            match ac {
                "NM_153046.3" => Ok(vec![(
                    IdentifierType::ProteinAccession,
                    "NP_694591.2".to_string(),
                )]),
                "NM_058216.3" => Ok(vec![(
                    IdentifierType::ProteinAccession,
                    "NP_478123.1".to_string(),
                )]),
                _ => Ok(vec![]),
            }
        } else {
            Ok(vec![])
        }
    }
    fn get_identifier_type(&self, _id: &str) -> Result<IdentifierType, HgvsError> {
        Ok(IdentifierType::TranscriptAccession)
    }
    fn c_to_g(
        &self,
        transcript_ac: &str,
        pos: TranscriptPos,
        offset: IntronicOffset,
    ) -> Result<(String, GenomicPos), HgvsError> {
        let tx = self.get_transcript(transcript_ac, None)?;
        Ok((
            tx.reference_accession().to_string(),
            GenomicPos(pos.0 + offset.0),
        ))
    }
}

#[test]
fn test_regression_c_360_eq() -> Result<(), HgvsError> {
    let hdp = RegressionProvider;
    let mapper = VariantMapper::new(&hdp);

    // NM_153046.3:c.360=
    let v_raw = hgvs_weaver::parse_hgvs_variant("NM_153046.3:c.360=")?;
    let SequenceVariant::Coding(var_c) = v_raw else {
        panic!()
    };
    let p_var = mapper.c_to_p(&var_c, None)?;

    // Should be synonymous (Thr120=)
    // Current bug makes it a frameshift deletion of c.360
    assert!(
        !p_var.to_string().contains("fs"),
        "Should not be a frameshift, got {}",
        p_var
    );
    Ok(())
}

#[test]
fn test_regression_delins_stop() -> Result<(), HgvsError> {
    let hdp = RegressionProvider;
    let mapper = VariantMapper::new(&hdp);

    // NM_058216.3:c.692_694delinsAA
    let v_raw = hgvs_weaver::parse_hgvs_variant("NM_058216.3:c.692_694delinsAA")?;
    let SequenceVariant::Coding(var_c) = v_raw else {
        panic!()
    };
    let p_var = mapper.c_to_p(&var_c, None)?;

    // We want it to be a stop codon if it created one, not a frameshift.
    // In our mock: 691=T. Insert AA -> TAA (Stop) at 231.
    assert!(
        p_var.to_string().contains("Ter") || p_var.to_string().contains("*"),
        "Should contain stop codon, got {}",
        p_var
    );
    assert!(
        !p_var.to_string().contains("fs"),
        "Should not be a frameshift if stop is earlier, got {}",
        p_var
    );

    Ok(())
}

struct RepeatProvider;
// UTR... ATG (1-3) CAG (4-6) CAG (7-9) CAG (10-12) TAG (13-15)
// M Q Q Q *
impl DataProvider for RepeatProvider {
    fn get_transcript(
        &self,
        _ac: &str,
        _ref: Option<&str>,
    ) -> Result<Box<dyn Transcript>, HgvsError> {
        Ok(Box::new(TranscriptData {
            ac: "NM_001.1".to_string(),
            gene: "TEST".to_string(),
            cds_start_index: Some(TranscriptPos(0)),
            cds_end_index: Some(TranscriptPos(14)),
            strand: 1,
            reference_accession: "NC_001.1".to_string(),
            exons: vec![],
        }))
    }
    fn get_seq(
        &self,
        _ac: &str,
        start: i32,
        end: i32,
        _kind: IdentifierType,
    ) -> Result<String, HgvsError> {
        let full_seq = "ATGCAGCAGCAGTAG";
        let s = start as usize;
        let e = if end == -1 {
            full_seq.len()
        } else {
            end as usize
        };
        if s < full_seq.len() && e <= full_seq.len() {
            Ok(full_seq[s..e].to_string())
        } else {
            Ok("".to_string())
        }
    }
    fn get_symbol_accessions(
        &self,
        _: &str,
        _: IdentifierKind,
        _: IdentifierKind,
    ) -> Result<Vec<(IdentifierType, String)>, HgvsError> {
        Ok(vec![(
            IdentifierType::ProteinAccession,
            "NP_001.1".to_string(),
        )])
    }
    fn get_identifier_type(&self, _: &str) -> Result<IdentifierType, HgvsError> {
        Ok(IdentifierType::TranscriptAccession)
    }
    fn c_to_g(
        &self,
        _: &str,
        pos: TranscriptPos,
        _: hgvs_weaver::structs::IntronicOffset,
    ) -> Result<(String, GenomicPos), HgvsError> {
        Ok(("NC_001.1".to_string(), GenomicPos(pos.0)))
    }
}

#[test]
fn test_regression_gln4del_vs_ter() -> Result<(), HgvsError> {
    let hdp = RepeatProvider;
    let mapper = VariantMapper::new(&hdp);

    // c.4_6del
    let v_c = hgvs_weaver::parse_hgvs_variant("NM_001.1:c.4_6del")?;

    if let SequenceVariant::Coding(vc) = v_c {
        // Generate p.
        let vp = mapper.c_to_p(&vc, Some("NP_001.1"))?;
        let vp_str = vp.to_string();
        println!("Generated: {}", vp_str);

        // Assert we get p.Gln4del (or equivalent del) AND NOT p.Gln4Ter
        assert!(
            vp_str.contains("del"),
            "Expected 'del' in {}, got {}",
            vp_str,
            vp_str
        );
        assert!(
            !vp_str.contains("Ter"),
            "Did not expect 'Ter' in {}, got {}",
            vp_str,
            vp_str
        );
    } else {
        panic!("Parsed wrong type");
    }

    Ok(())
}

struct DelinsMismatchProvider;
impl DataProvider for DelinsMismatchProvider {
    fn get_transcript(
        &self,
        _ac: &str,
        _ref: Option<&str>,
    ) -> Result<Box<dyn Transcript>, HgvsError> {
        #[derive(Clone)]
        struct MockTranscript {
            exons: Vec<ExonData>,
        }
        impl Transcript for MockTranscript {
            fn ac(&self) -> &str {
                "NM_001008844.3"
            }
            fn gene(&self) -> &str {
                "TEST"
            }
            fn strand(&self) -> i32 {
                1
            }
            fn cds_start_index(&self) -> Option<TranscriptPos> {
                Some(TranscriptPos(0))
            }
            fn cds_end_index(&self) -> Option<TranscriptPos> {
                Some(TranscriptPos(4500))
            }
            fn reference_accession(&self) -> &str {
                "NC_000001.11"
            }
            fn exons(&self) -> &[ExonData] {
                &self.exons
            }
        }

        let exons = vec![ExonData {
            transcript_start: TranscriptPos(0),
            transcript_end: TranscriptPos(5000),
            reference_start: GenomicPos(0),
            reference_end: GenomicPos(5000),
            alt_strand: 1,
            cigar: "5000M".to_string(),
        }];

        Ok(Box::new(MockTranscript { exons }))
    }

    fn get_seq(
        &self,
        _ac: &str,
        start: i32,
        end: i32,
        _kind: IdentifierType,
    ) -> Result<String, HgvsError> {
        let effective_end = if end == -1 { 5000 } else { end };
        let mut seq = String::with_capacity((effective_end - start) as usize);
        for i in start..effective_end {
            if i >= 4497 && i <= 4499 {
                if i == 4497 {
                    seq.push('C');
                }
                if i == 4498 {
                    seq.push('C');
                }
                if i == 4499 {
                    seq.push('A');
                }
            } else {
                seq.push('G');
            }
        }
        Ok(seq)
    }

    fn get_symbol_accessions(
        &self,
        symbol: &str,
        source_kind: IdentifierKind,
        target_kind: IdentifierKind,
    ) -> Result<Vec<(IdentifierType, String)>, HgvsError> {
        if symbol == "NM_001008844.3"
            && source_kind == IdentifierKind::Transcript
            && target_kind == IdentifierKind::Protein
        {
            return Ok(vec![(
                IdentifierType::ProteinAccession,
                "NP_001008844.1".to_string(),
            )]);
        }
        Ok(vec![])
    }
    fn get_identifier_type(&self, _: &str) -> Result<IdentifierType, HgvsError> {
        Ok(IdentifierType::TranscriptAccession)
    }
    fn c_to_g(
        &self,
        _: &str,
        _: TranscriptPos,
        _: IntronicOffset,
    ) -> Result<(String, GenomicPos), HgvsError> {
        Ok(("".into(), GenomicPos(0)))
    }
}

#[test]
fn test_regression_pro_ile_mismatch() -> Result<(), HgvsError> {
    let provider = DelinsMismatchProvider;
    let mapper = VariantMapper::new(&provider);

    let v_nuc = hgvs_weaver::parse_hgvs_variant("NM_001008844.3:c.4498_4499delinsAT")?;

    if let SequenceVariant::Coding(c_var) = v_nuc {
        let v_p = mapper.c_to_p(&c_var, None)?;
        assert!(v_p.to_string().contains("Pro1500Ile"));
    } else {
        panic!("Parsed variant is not coding");
    }
    Ok(())
}

#[test]
fn test_regression_parse_clinvar_repeat() {
    use hgvs_weaver::parse_hgvs_variant;
    let v = parse_hgvs_variant("NP_001365049.1:p.490PRS[1]");
    assert!(v.is_ok(), "Failed to parse p.490PRS[1]: {:?}", v.err());
    let v_inner = v.unwrap();

    if let SequenceVariant::Protein(p) = v_inner {
        if let hgvs_weaver::edits::AaEdit::Repeat { ref_, min, max, .. } = p.posedit.edit {
            assert_eq!(ref_, Some("PRS".to_string()));
            assert_eq!(min, 1);
            assert_eq!(max, 1);
        } else {
            panic!(
                "Expected Repeat edit, got something else: {:?}",
                p.posedit.edit
            );
        }
    }
}
