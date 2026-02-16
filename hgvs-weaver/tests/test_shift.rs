use hgvs_weaver::coords::{GenomicPos, IntronicOffset, TranscriptPos};
use hgvs_weaver::data::{
    DataProvider, ExonData, IdentifierKind, IdentifierType, Transcript, TranscriptData,
};
use hgvs_weaver::error::HgvsError;
use hgvs_weaver::mapper::VariantMapper;

struct HomopolymerProvider;
impl DataProvider for HomopolymerProvider {
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
        // Return 100 'A's
        let seq = "A".repeat(2000);
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

#[test]
fn test_ins_3_prime_shifting() -> Result<(), HgvsError> {
    let hdp = HomopolymerProvider;
    let mapper = VariantMapper::new(&hdp);

    // NC_TEST.1:g.1005_1006insA
    let v1 = hgvs_weaver::parse_hgvs_variant("NC_TEST.1:g.1005_1006insA")?;
    let SequenceVariant::Genomic(v1_g) = v1 else {
        panic!()
    };
    let nv1 = mapper.normalize_variant(SequenceVariant::Genomic(v1_g))?;

    // NC_TEST.1:g.1006_1007insA
    let v2 = hgvs_weaver::parse_hgvs_variant("NC_TEST.1:g.1006_1007insA")?;
    let SequenceVariant::Genomic(v2_g) = v2 else {
        panic!()
    };
    let nv2 = mapper.normalize_variant(SequenceVariant::Genomic(v2_g))?;

    assert_eq!(nv1.to_string(), nv2.to_string());
    Ok(())
}

use hgvs_weaver::SequenceVariant;
