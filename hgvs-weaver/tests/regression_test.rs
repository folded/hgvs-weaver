use hgvs_weaver::coords::{GenomicPos, IntronicOffset, TranscriptPos};
use hgvs_weaver::data::{
    DataProvider, IdentifierKind, IdentifierType, Transcript, TranscriptData, TranscriptSearch,
};
use hgvs_weaver::equivalence::{EquivalenceLevel, VariantEquivalence};
use hgvs_weaver::error::HgvsError;
use hgvs_weaver::parse_hgvs_variant;

struct MockDataProvider;
impl DataProvider for MockDataProvider {
    fn get_transcript(&self, ac: &str, _: Option<&str>) -> Result<Box<dyn Transcript>, HgvsError> {
        if ac == "NM_001166478.1" || ac == "NM_005813.3" {
            Ok(Box::new(TranscriptData {
                ac: ac.to_string(),
                gene: "TEST".to_string(),
                cds_start_index: Some(TranscriptPos(0)),
                cds_end_index: Some(TranscriptPos(3000)),
                strand: -1,
                reference_accession: "NC_000001.1".to_string(),
                exons: vec![],
            }))
        } else if ac == "NM_BRAF" {
            Ok(Box::new(TranscriptData {
                ac: ac.to_string(),
                gene: "BRAF".to_string(),
                cds_start_index: Some(TranscriptPos(0)),
                cds_end_index: Some(TranscriptPos(3000)),
                strand: 1, // Plus strand
                reference_accession: "NC_BRAF".to_string(),
                exons: vec![],
            }))
        } else {
            Err(HgvsError::ValidationError("Not found".into()))
        }
    }
    fn get_seq(&self, ac: &str, s: i32, e: i32, _k: IdentifierType) -> Result<String, HgvsError> {
        let effective_e = if e == -1 { 4000 } else { e };
        let len = (effective_e - s).max(0) as usize;
        let mut seq = vec!['N'; len];

        if ac == "NM_BRAF" || ac == "NC_BRAF" {
            // BRAF Val600 is GTG (1798-1800)
            // 1798 is 'G', 1799 is 'T', 1800 is 'G'
            for pos in s..effective_e {
                let char_idx = (pos - s) as usize;
                if pos == 1797 {
                    seq[char_idx] = 'G';
                } else if pos == 1798 {
                    seq[char_idx] = 'T';
                } else if pos == 1799 {
                    seq[char_idx] = 'G';
                }
            }
            return Ok(seq.into_iter().collect());
        }

        if s == 3966 && len == 1 {
            // Case 9: c.35 matches G3966.
            Ok("A".to_string())
        } else if s == 1328 && len == 1 {
            // Case 15: c.2673 matches G1328.
            Ok("T".to_string())
        } else {
            Ok(seq.into_iter().collect())
        }
    }
    fn get_symbol_accessions(
        &self,
        _s: &str,
        _f: IdentifierKind,
        _t: IdentifierKind,
    ) -> Result<Vec<(IdentifierType, String)>, HgvsError> {
        Ok(vec![])
    }
    fn get_identifier_type(&self, _id: &str) -> Result<IdentifierType, HgvsError> {
        Ok(IdentifierType::TranscriptAccession)
    }
    fn c_to_g(
        &self,
        _t: &str,
        pos: TranscriptPos,
        offset: IntronicOffset,
    ) -> Result<(String, GenomicPos), HgvsError> {
        // Reverse mapping for minus strand
        Ok((
            "NC_000001.1".to_string(),
            GenomicPos(4000 - pos.0 - offset.0),
        ))
    }
}

struct MockSearch;
impl TranscriptSearch for MockSearch {
    fn get_transcripts_for_region(
        &self,
        _: &str,
        _: i32,
        _: i32,
    ) -> Result<Vec<String>, HgvsError> {
        Ok(vec![])
    }
}

#[test]
fn test_repro_case9() -> Result<(), HgvsError> {
    let hdp = MockDataProvider;
    let search = MockSearch;
    let eq = VariantEquivalence::new(&hdp, &search);

    let v1 = parse_hgvs_variant("NM_001166478.1:c.35_36insT")?;
    let v2 = parse_hgvs_variant("NM_001166478.1:c.35dup")?;

    assert_eq!(eq.equivalent_level(&v1, &v2)?, EquivalenceLevel::Analogous);
    Ok(())
}

#[test]
fn test_repro_case15() -> Result<(), HgvsError> {
    let hdp = MockDataProvider;
    let search = MockSearch;
    let eq = VariantEquivalence::new(&hdp, &search);

    let v1 = parse_hgvs_variant("NM_005813.3:c.2673insA")?;
    let v2 = parse_hgvs_variant("NM_005813.3:c.2673dup")?;

    assert_eq!(eq.equivalent_level(&v1, &v2)?, EquivalenceLevel::Analogous);
    Ok(())
}

#[test]
fn test_braf_identity() -> Result<(), HgvsError> {
    let hdp = MockDataProvider;
    let search = MockSearch;
    let eq = VariantEquivalence::new(&hdp, &search);

    // c.1799T>A -> p.Val600Glu (predicted)
    let v1 = parse_hgvs_variant("NM_BRAF:c.1799T>A")?;

    // Target 1: p.Val600Glu (Observed/Experimental)
    // Comparisons of c. to p. (observed) are Analogous because c. implies p.(predicted).
    // The notation differs (parens vs no parens), so it's not Strict Identity.
    let v2_observed = parse_hgvs_variant("NP_BRAF:p.Val600Glu")?;
    assert_eq!(
        eq.equivalent_level(&v1, &v2_observed)?,
        EquivalenceLevel::Analogous
    );

    // Target 2: p.(Val600Glu) (Predicted)
    // Comparisons of c. to p.(predicted) should be Identity because c. implies p.(predicted) exactly.
    let v2_predicted = parse_hgvs_variant("NP_BRAF:p.(Val600Glu)")?;

    assert_eq!(
        eq.equivalent_level(&v1, &v2_predicted)?,
        EquivalenceLevel::Identity
    );

    Ok(())
}
