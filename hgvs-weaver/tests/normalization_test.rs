use hgvs_weaver::data::{ExonData, TranscriptData};
use hgvs_weaver::structs::{GenomicPos, IntronicOffset, TranscriptPos};
use hgvs_weaver::*;

struct NormMockDataProvider;

impl DataProvider for NormMockDataProvider {
    fn get_seq(
        &self,
        ac: &str,
        start: i32,
        end: i32,
        _kind: hgvs_weaver::data::IdentifierType,
    ) -> Result<String, HgvsError> {
        let base_seq = if ac == "NM_SHIFT_BUG" {
            "CCATTTTTTT".to_string()
        } else if ac == "NM_PREMATURE_STOP" {
            "ATGCAACAAGATGATTAA".to_string()
        } else if ac == "NM_INFRAME_DEL" {
            "ATGGCTGCATGCGATTAA".to_string()
        } else if ac == "NM_CTERM_SUBST" {
            "ATGGCTGCATGCGATTAA".to_string()
        } else if ac == "NM_REPEAT_EXP" {
            "ATGGCTGCTGCTTTTTAA".to_string()
        } else if ac == "NM_REPEAT_CON" {
            "ATGGCAGCAGCAGCATTTTAA".to_string()
        } else {
            let mut s = String::new();
            s.push_str("AAAAAAAAAA"); // 10 A's
            s.push_str("ATGAAATAG");
            for _ in 0..20 {
                s.push_str("ATGC");
            }
            s
        };

        let s = start as usize;
        let e = if end == -1 {
            base_seq.len()
        } else {
            end as usize
        };
        if s > base_seq.len() {
            return Ok("".into());
        }
        let e = e.min(base_seq.len());
        Ok(base_seq[s..e].to_string())
    }

    fn get_transcript(
        &self,
        transcript_ac: &str,
        _reference_ac: Option<&str>,
    ) -> Result<Box<dyn Transcript>, HgvsError> {
        let exons = vec![ExonData {
            transcript_start: TranscriptPos(0),
            transcript_end: TranscriptPos(100),
            reference_start: GenomicPos(1000),
            reference_end: GenomicPos(1100),
            alt_strand: 1,
            cigar: "100M".to_string(),
        }];

        let (cds_start, cds_end) = match transcript_ac {
            "NM_0001.1" => (10, 19), // Met Lys *
            "NM_SHIFT_BUG" => (0, 30),
            "NM_PREMATURE_STOP" => (0, 18), // M Q Q D D * (18 bases)
            "NM_INFRAME_DEL" => (0, 18),    // M A B C D * (18 bases)
            "NM_CTERM_SUBST" => (0, 18),    // M A B C D *
            "NM_REPEAT_EXP" => (0, 18),     // M A A A F * (ATG GCT GCT GCT TTT TAA)
            "NM_REPEAT_CON" => (0, 21),     // M A A A A F * (ATG GCA GCA GCA GCA TTT TAA)
            _ => {
                return Err(HgvsError::DataProviderError(
                    "Transcript not found".to_string(),
                ))
            }
        };

        let td = TranscriptData {
            ac: transcript_ac.to_string(),
            gene: "NORM".to_string(),
            cds_start_index: Some(TranscriptPos(cds_start)),
            cds_end_index: Some(TranscriptPos(cds_end)),
            strand: 1,
            reference_accession: "NC_0001.10".to_string(),
            exons,
        };
        Ok(Box::new(td))
    }

    fn get_symbol_accessions(
        &self,
        symbol: &str,
        _sk: hgvs_weaver::data::IdentifierKind,
        tk: hgvs_weaver::data::IdentifierKind,
    ) -> Result<Vec<(hgvs_weaver::data::IdentifierType, String)>, HgvsError> {
        if tk == hgvs_weaver::data::IdentifierKind::Protein && symbol == "NM_0001.1" {
            return Ok(vec![(
                hgvs_weaver::data::IdentifierType::ProteinAccession,
                "NP_0001.1".to_string(),
            )]);
        }
        Ok(vec![(
            hgvs_weaver::data::IdentifierType::Unknown,
            symbol.to_string(),
        )])
    }

    fn get_identifier_type(
        &self,
        _identifier: &str,
    ) -> Result<hgvs_weaver::data::IdentifierType, HgvsError> {
        Ok(hgvs_weaver::data::IdentifierType::Unknown)
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
fn test_nonsense_normalization() {
    let hdp = NormMockDataProvider;
    let mapper = VariantMapper::new(&hdp);

    // c.4A>T changes AAA (Lys) to TAA (Stop).
    let var_c = parse_hgvs_variant("NM_0001.1:c.4A>T").unwrap();
    if let SequenceVariant::Coding(v) = var_c {
        let var_p = mapper.c_to_p(&v, Some("NP_0001.1")).unwrap();
        // Should be Lys2Ter
        assert_eq!(var_p.to_string(), "NP_0001.1:p.(Lys2Ter)");
    }
}

#[test]
fn test_extension_normalization() {
    let hdp = NormMockDataProvider;
    let mapper = VariantMapper::new(&hdp);

    // c.7T>G changes TAG (Stop) to GAG (Glu).
    let var_c = parse_hgvs_variant("NM_0001.1:c.7T>G").unwrap();
    if let SequenceVariant::Coding(v) = var_c {
        let var_p = mapper.c_to_p(&v, Some("NP_0001.1")).unwrap();
        // Should be {alt}extTer{length}
        assert!(var_p.to_string().contains("extTer"));
    }
}

#[test]
fn test_normalization_shift_bug() {
    let hdp = NormMockDataProvider;
    let mapper = VariantMapper::new(&hdp);

    // NM_SHIFT_BUG: Sequence "CCAT...". UTR=0.
    // c.1_2delinsAT. Ref=CC. Alt=AT.
    // Result Sequence: ATAT...
    // If shifted -> c.3_4delinsAT. Ref=AT. Alt=AT.
    // Result Sequence: CCAT...
    // The outputs are DIFFERENT. So it MUST NOT shift.

    let var_c = parse_hgvs_variant("NM_SHIFT_BUG:c.1_2delinsAT").unwrap();
    if let SequenceVariant::Coding(v) = var_c {
        let var_norm = mapper
            .normalize_variant(SequenceVariant::Coding(v))
            .unwrap();
        if let SequenceVariant::Coding(v_norm) = var_norm {
            // Should remain c.1_2
            // Because if it shifts to 3_4, it implies AT -> AT, which is Identity.
            let pos = v_norm.posedit.pos.unwrap();
            let start = pos.start.base.to_index().0;
            // 0-based index: 0
            assert_eq!(start, 0, "Variant should not have shifted!");
        }
    }
}

#[test]
fn test_premature_stop_formatting() {
    let hdp = NormMockDataProvider;
    let mapper = VariantMapper::new(&hdp);

    // NM_PREMATURE_STOP: M Q Q D D *.
    // c.4_9delinsCATTAA.
    // Replace CAA CAA (Q Q) with CAT TAA (H *).
    // Result: M H *.
    // Without fix: p.Gln2_Asp5delinsHis. (Deletion of QQD...).
    // With fix: p.Gln2_Gln3delinsHisTer.

    let var_c = parse_hgvs_variant("NM_PREMATURE_STOP:c.4_9delinsCATTAA").unwrap();
    if let SequenceVariant::Coding(v) = var_c {
        let var_p = mapper.c_to_p(&v, Some("NP_MOCK")).unwrap();
        assert_eq!(var_p.to_string(), "NP_MOCK:p.(Gln2_Gln3delinsHisTer)");
    }
}

#[test]
fn test_inframe_deletion_tail() {
    let hdp = NormMockDataProvider;
    let mapper = VariantMapper::new(&hdp);

    // NM_INFRAME_DEL: M A A C D *.
    // c.4_6del. Del A (second A).
    // Result: M A C D *.
    // Tail match: C D * (3 chars).
    // This is > 1 char. So it should be detected as Original Stop.
    // Result should be p.Ala2del. NOT delins...Ter.

    let var_c = parse_hgvs_variant("NM_INFRAME_DEL:c.4_6del").unwrap();
    if let SequenceVariant::Coding(v) = var_c {
        let var_p = mapper.c_to_p(&v, Some("NP_MOCK")).unwrap();
        assert_eq!(var_p.to_string(), "NP_MOCK:p.(Ala3del)");
    }
}

#[test]
fn test_cterm_substitution() {
    let hdp = NormMockDataProvider;
    let mapper = VariantMapper::new(&hdp);

    // NM_CTERM_SUBST: M A B C D *.
    // c.15T>G. D (GAT) -> E (GAG).
    // Tail match: *. (1 char).
    // Ref mismatch length: 1 (D).
    // Alt mismatch length: 1 (E).
    // Should be p.Asp5Glu. NOT delins...Ter.

    let var_c = parse_hgvs_variant("NM_CTERM_SUBST:c.15T>G").unwrap();
    if let SequenceVariant::Coding(v) = var_c {
        let var_p = mapper.c_to_p(&v, Some("NP_MOCK")).unwrap();
        assert_eq!(var_p.to_string(), "NP_MOCK:p.(Asp5Glu)");
    }
}

#[test]
fn test_repeat_expansion() {
    let hdp = NormMockDataProvider;
    let mapper = VariantMapper::new(&hdp);

    // NM_REPEAT_EXP: M A A A *. (ATG GCT GCT GCT TAA).
    // c.4GCT[5] -> p.Ala3_Ala4dup? Or p.Ala2_Ala3dup?
    // Ref has 3 copies. Var has 5. Net +2 copies (+6 bases).
    // Expecting duplication notation.

    let var_c = parse_hgvs_variant("NM_REPEAT_EXP:c.4GCT[5]").unwrap();
    if let SequenceVariant::Coding(v) = var_c {
        let var_p = mapper.c_to_p(&v, Some("NP_MOCK")).unwrap();
        assert!(
            var_p.to_string().contains("dup"),
            "Expected dup, got {}",
            var_p
        );
    }
}

#[test]
fn test_repeat_contraction() {
    let hdp = NormMockDataProvider;
    let mapper = VariantMapper::new(&hdp);

    // NM_REPEAT_CON: M A A A A *. (ATG GCA GCA GCA GCA TAA).
    // c.4GCA[2].
    // Ref has 4 copies. Var has 2. Net -2 copies (-6 bases).
    // Expecting deletion notation.

    let var_c = parse_hgvs_variant("NM_REPEAT_CON:c.4GCA[2]").unwrap();
    if let SequenceVariant::Coding(v) = var_c {
        let var_p = mapper.c_to_p(&v, Some("NP_MOCK")).unwrap();
        assert!(
            var_p.to_string().contains("del"),
            "Expected del, got {}",
            var_p
        );
    }
}
