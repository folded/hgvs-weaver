use hgvs_weaver::data::{ExonData, TranscriptData};
use hgvs_weaver::structs::{GenomicPos, IntronicOffset, TranscriptPos};
use hgvs_weaver::*;

struct MockDataProvider;

impl DataProvider for MockDataProvider {
    fn get_seq(
        &self,
        _ac: &str,
        start: i32,
        end: i32,
        _kind: hgvs_weaver::data::IdentifierType,
    ) -> Result<String, HgvsError> {
        let mut s = String::new();
        s.push_str("AAAAAAAAAA"); // 10 A's
        s.push_str("ATG"); // n.11 is c.1
        for _ in 0..25 {
            s.push_str("ATGC");
        }

        let start = start as usize;
        let end = if end == -1 { s.len() } else { end as usize };
        if start > s.len() {
            return Ok("".into());
        }
        let end = end.min(s.len());
        Ok(s[start..end].to_string())
    }

    fn get_transcript(
        &self,
        transcript_ac: &str,
        _reference_ac: Option<&str>,
    ) -> Result<Box<dyn Transcript>, HgvsError> {
        if transcript_ac == "NM_0001.3" {
            let exons = vec![ExonData {
                transcript_start: TranscriptPos(0),
                transcript_end: TranscriptPos(100),
                reference_start: GenomicPos(1000),
                reference_end: GenomicPos(1100),
                alt_strand: 1,
                cigar: "100M".to_string(),
            }];
            let td = TranscriptData {
                ac: "NM_0001.3".to_string(),
                gene: "MOCK".to_string(),
                cds_start_index: Some(TranscriptPos(10)), // n.11 is c.1
                cds_end_index: Some(TranscriptPos(50)),
                strand: 1,
                reference_accession: "NC_0001.10".to_string(),
                exons,
            };
            return Ok(Box::new(td));
        }
        Err(HgvsError::DataProviderError(
            "Transcript not found".to_string(),
        ))
    }

    fn get_symbol_accessions(
        &self,
        symbol: &str,
        _sk: hgvs_weaver::data::IdentifierKind,
        tk: hgvs_weaver::data::IdentifierKind,
    ) -> Result<Vec<(hgvs_weaver::data::IdentifierType, String)>, HgvsError> {
        if tk == hgvs_weaver::data::IdentifierKind::Protein && symbol == "NM_0001.3" {
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
fn test_mapper_c_to_p_subst() {
    let hdp = MockDataProvider;
    let mapper = VariantMapper::new(&hdp);

    let var_c = parse_hgvs_variant("NM_0001.3:c.1A>T").unwrap();
    if let SequenceVariant::Coding(v) = var_c {
        let var_p = mapper.c_to_p(&v, Some("NP_0001.1")).unwrap();
        assert_eq!(var_p.to_string(), "NP_0001.1:p.(Met1Leu)");
    }
}

#[test]
fn test_mapper_c_to_p_fs() {
    let hdp = MockDataProvider;
    let mapper = VariantMapper::new(&hdp);

    // c.2del deletes 'T' from 'ATG'
    let var_c = parse_hgvs_variant("NM_0001.3:c.2del").unwrap();
    if let SequenceVariant::Coding(v) = var_c {
        let var_p = mapper.c_to_p(&v, Some("NP_0001.1")).unwrap();
        assert!(var_p.to_string().contains("fsTer"));
    }
}

#[test]
fn test_mapper_g_to_c_3utr() {
    let hdp = MockDataProvider;
    let mapper = VariantMapper::new(&hdp);

    // Genomic 1052 (index 1051) -> n.52 -> c.*1
    let var_g = parse_hgvs_variant("NC_0001.10:g.1052A>T").unwrap();
    if let SequenceVariant::Genomic(v) = var_g {
        let var_c = mapper.g_to_c(&v, "NM_0001.3").unwrap();
        assert_eq!(var_c.to_string(), "NM_0001.3:c.*1A>T");
    }
}

#[test]
fn test_mapper_c_to_g_3utr() {
    let hdp = MockDataProvider;
    let mapper = VariantMapper::new(&hdp);

    let var_c = parse_hgvs_variant("NM_0001.3:c.*1A>T").unwrap();
    if let SequenceVariant::Coding(v) = var_c {
        let var_g = mapper.c_to_g(&v, Some("NC_0001.10")).unwrap();
        assert_eq!(var_g.to_string(), "NC_0001.10:g.1052A>T");
    }
}
