use hgvs_weaver::coords::{GenomicPos, IntronicOffset, TranscriptPos};
use hgvs_weaver::data::{
    DataProvider, ExonData, IdentifierKind, IdentifierType, Transcript, TranscriptData,
    TranscriptSearch,
};
use hgvs_weaver::equivalence::{EquivalenceLevel, VariantEquivalence};
use hgvs_weaver::error::HgvsError;

struct SimpleProvider;
impl DataProvider for SimpleProvider {
    fn get_transcript(
        &self,
        ac: &str,
        _ref_ac: Option<&str>,
    ) -> Result<Box<dyn Transcript>, HgvsError> {
        Ok(Box::new(TranscriptData {
            ac: ac.to_string(),
            gene: "TEST".to_string(),
            cds_start_index: Some(TranscriptPos(0)),
            cds_end_index: Some(TranscriptPos(100)),
            strand: 1,
            reference_accession: "NC_TEST.1".to_string(),
            exons: vec![ExonData {
                transcript_start: TranscriptPos(0),
                transcript_end: TranscriptPos(100),
                reference_start: GenomicPos(1000),
                reference_end: GenomicPos(1100),
                alt_strand: 1,
                cigar: "100M".to_string(),
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
        let seq = "ACGT".repeat(1000);
        let s = start as usize;
        let e = end as usize;
        if s < seq.len() {
            let actual_e = e.min(seq.len());
            Ok(seq[s..actual_e].to_string())
        } else {
            Ok("".to_string())
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
        Ok(IdentifierType::GenomicAccession)
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

impl TranscriptSearch for SimpleProvider {
    fn get_transcripts_for_region(
        &self,
        _chrom: &str,
        _start: i32,
        _end: i32,
    ) -> Result<Vec<String>, HgvsError> {
        Ok(vec![])
    }
}

#[test]
fn test_equivalence_levels() -> Result<(), HgvsError> {
    let hdp = SimpleProvider;
    let eq = VariantEquivalence::new(&hdp, &hdp);

    // 1. Identity
    let v1 = hgvs_weaver::parse_hgvs_variant("NC_TEST.1:g.1001A>C")?;
    let v2 = hgvs_weaver::parse_hgvs_variant("NC_TEST.1:g.1001A>C")?;
    assert_eq!(eq.equivalent_level(&v1, &v2)?, EquivalenceLevel::Identity);

    // 2. Analogous (Normalization) - ins vs dup
    let v3 = hgvs_weaver::parse_hgvs_variant("NC_TEST.1:g.1005dupC")?;
    let v4 = hgvs_weaver::parse_hgvs_variant("NC_TEST.1:g.1005_1006insC")?;
    let lvl = eq.equivalent_level(&v3, &v4)?;
    assert!(lvl == EquivalenceLevel::Identity || lvl == EquivalenceLevel::Analogous);

    // 3. Parity (Functional - same protein)
    let p1 = hgvs_weaver::parse_hgvs_variant("NP_0001.1:p.Trp2Ter")?;
    let p2 = hgvs_weaver::parse_hgvs_variant("NP_0001.1:p.Trp2*")?;
    assert_eq!(eq.equivalent_level(&p1, &p2)?, EquivalenceLevel::Identity);

    Ok(())
}

#[test]
fn test_parity_match_unification() -> Result<(), HgvsError> {
    struct MissingSeqProvider;
    impl DataProvider for MissingSeqProvider {
        fn get_transcript(
            &self,
            ac: &str,
            _: Option<&str>,
        ) -> Result<Box<dyn Transcript>, HgvsError> {
            Ok(Box::new(TranscriptData {
                ac: ac.to_string(),
                gene: "TEST".to_string(),
                cds_start_index: Some(TranscriptPos(0)),
                cds_end_index: Some(TranscriptPos(1000)),
                strand: 1,
                reference_accession: "NC_TEST.1".to_string(),
                exons: vec![],
            }))
        }
        fn get_seq(&self, _: &str, _: i32, _: i32, _: IdentifierType) -> Result<String, HgvsError> {
            Err(HgvsError::Other("Missing sequence".into()))
        }
        fn get_symbol_accessions(
            &self,
            _: &str,
            _: IdentifierKind,
            _: IdentifierKind,
        ) -> Result<Vec<(IdentifierType, String)>, HgvsError> {
            Ok(vec![])
        }
        fn get_identifier_type(&self, _: &str) -> Result<IdentifierType, HgvsError> {
            Ok(IdentifierType::TranscriptAccession)
        }
        fn c_to_g(
            &self,
            _: &str,
            p: TranscriptPos,
            _: IntronicOffset,
        ) -> Result<(String, GenomicPos), HgvsError> {
            Ok(("NC_TEST.1".to_string(), GenomicPos(p.0)))
        }
    }
    impl TranscriptSearch for MissingSeqProvider {
        fn get_transcripts_for_region(
            &self,
            _: &str,
            _: i32,
            _: i32,
        ) -> Result<Vec<String>, HgvsError> {
            Ok(vec![])
        }
    }

    let hdp = MissingSeqProvider;
    let eq = VariantEquivalence::new(&hdp, &hdp);

    // Case: p.Ala201_Val202insGlyProGlyAla vs p.(Gly198_Ala201dup)
    // The insertion variant seeds 201=Ala, 202=Val.
    // The duplication variant seeds 198=Gly, 201=Ala.
    // The gaps 199, 200 remain Unknown.
    let v1 = hgvs_weaver::parse_hgvs_variant("NP_0001.1:p.Ala201_Val202insGlyProGlyAla")?;
    let v2 = hgvs_weaver::parse_hgvs_variant("NP_0001.1:p.Gly198_Ala201dup")?;

    // This should now return Weak equivalence via unification!
    assert_eq!(eq.equivalent_level(&v1, &v2)?, EquivalenceLevel::Analogous);

    Ok(())
}
