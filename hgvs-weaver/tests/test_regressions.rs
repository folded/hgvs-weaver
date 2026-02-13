use hgvs_weaver::data::{DataProvider, IdentifierType, IdentifierKind, TranscriptData, ExonData, Transcript};
use hgvs_weaver::coords::{GenomicPos, TranscriptPos, IntronicOffset, SequenceVariant};
use hgvs_weaver::mapper::VariantMapper;
use hgvs_weaver::error::HgvsError;

struct RegressionProvider;
impl DataProvider for RegressionProvider {
    fn get_transcript(&self, ac: &str, _ref_ac: Option<&str>) -> Result<Box<dyn Transcript>, HgvsError> {
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
    fn get_seq(&self, _ac: &str, start: i32, end: i32, _kind: IdentifierType) -> Result<String, HgvsError> {
        let mut seq = "ACGC".repeat(1000).into_bytes();
        // For NM_058216.3:c.692_694delinsAA
        // Ser231: 691, 692, 693
        // If 691 is T, 692-693 replaced by AA -> TAA (Stop)
        if seq.len() > 690 {
            seq[690] = b'T';
        }
        let s = start as usize;
        let e = if end < 0 { seq.len() } else { end as usize };
        if s > seq.len() { return Ok("".to_string()); }
        let actual_e = e.min(seq.len());
        Ok(String::from_utf8_lossy(&seq[s..actual_e]).to_string())
    }
    fn get_symbol_accessions(&self, ac: &str, _f: IdentifierKind, t: IdentifierKind) -> Result<Vec<(IdentifierType, String)>, HgvsError> { 
        if t == IdentifierKind::Protein {
             match ac {
                 "NM_153046.3" => Ok(vec![(IdentifierType::ProteinAccession, "NP_694591.2".to_string())]),
                 "NM_058216.3" => Ok(vec![(IdentifierType::ProteinAccession, "NP_478123.1".to_string())]),
                 _ => Ok(vec![]),
             }
        } else {
            Ok(vec![])
        }
    }
    fn get_identifier_type(&self, _id: &str) -> Result<IdentifierType, HgvsError> { Ok(IdentifierType::TranscriptAccession) }
    fn c_to_g(&self, transcript_ac: &str, pos: TranscriptPos, offset: IntronicOffset) -> Result<(String, GenomicPos), HgvsError> {
        let tx = self.get_transcript(transcript_ac, None)?;
        Ok((tx.reference_accession().to_string(), GenomicPos(pos.0 + offset.0)))
    }
}

#[test]
fn test_regression_c_360_eq() -> Result<(), HgvsError> {
    let hdp = RegressionProvider;
    let mapper = VariantMapper::new(&hdp);

    // NM_153046.3:c.360= 
    let v_raw = hgvs_weaver::parse_hgvs_variant("NM_153046.3:c.360=")?;
    let SequenceVariant::Coding(var_c) = v_raw else { panic!() };
    let p_var = mapper.c_to_p(&var_c, None)?;
    
    // Should be synonymous (Thr120=)
    // Current bug makes it a frameshift deletion of c.360
    assert!(!p_var.to_string().contains("fs"), "Should not be a frameshift, got {}", p_var);
    Ok(())
}

#[test]
fn test_regression_delins_stop() -> Result<(), HgvsError> {
    let hdp = RegressionProvider;
    let mapper = VariantMapper::new(&hdp);

    // NM_058216.3:c.692_694delinsAA
    let v_raw = hgvs_weaver::parse_hgvs_variant("NM_058216.3:c.692_694delinsAA")?;
    let SequenceVariant::Coding(var_c) = v_raw else { panic!() };
    let p_var = mapper.c_to_p(&var_c, None)?;
    
    // We want it to be a stop codon if it created one, not a frameshift.
    // In our mock: 691=T. Insert AA -> TAA (Stop) at 231.
    assert!(p_var.to_string().contains("Ter") || p_var.to_string().contains("*"), "Should contain stop codon, got {}", p_var);
    assert!(!p_var.to_string().contains("fs"), "Should not be a frameshift if stop is earlier, got {}", p_var);
    
    Ok(())
}
