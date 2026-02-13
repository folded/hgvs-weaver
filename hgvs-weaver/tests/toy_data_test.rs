use hgvs_weaver::*;
use hgvs_weaver::data::TranscriptData;
use hgvs_weaver::structs::{TranscriptPos, GenomicPos, IntronicOffset};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;

#[derive(Serialize, Deserialize)]
struct ToyData {
    sequences: HashMap<String, String>,
    transcripts: HashMap<String, TranscriptData>,
}

struct JsonDataProvider {
    data: ToyData,
}

impl JsonDataProvider {
    fn new(path: &str) -> Self {
        let file = File::open(path).expect("Failed to open toy data file");
        let reader = BufReader::new(file);
        let data: ToyData = serde_json::from_reader(reader).expect("Failed to parse toy data");
        JsonDataProvider { data }
    }
}

impl DataProvider for JsonDataProvider {
    fn get_seq(&self, ac: &str, start: i32, end: i32, _kind: hgvs_weaver::data::IdentifierType) -> Result<String, HgvsError> {
        let seq = self.data.sequences.get(ac).ok_or_else(|| HgvsError::DataProviderError(format!("Sequence {} not found", ac)))?;
        let len = seq.len() as i32;
        let actual_end = if end == -1 { len } else { end };
        if start < 0 || actual_end > len || start > actual_end {
            return Err(HgvsError::DataProviderError("Sequence range out of bounds".into()));
        }
        Ok(seq[(start as usize)..(actual_end as usize)].to_string())
    }

    fn get_transcript(&self, transcript_ac: &str, _reference_accession: Option<&str>) -> Result<Box<dyn Transcript>, HgvsError> {
        let td = self.data.transcripts.get(transcript_ac)
            .ok_or_else(|| HgvsError::DataProviderError(format!("Transcript {} not found", transcript_ac)))?;
        Ok(Box::new(td.clone()))
    }

    fn get_symbol_accessions(&self, symbol: &str, _sk: hgvs_weaver::data::IdentifierKind, tk: hgvs_weaver::data::IdentifierKind) -> Result<Vec<(hgvs_weaver::data::IdentifierType, String)>, HgvsError> {
        if tk == hgvs_weaver::data::IdentifierKind::Protein {
            if symbol == "NM_PLUS.1" { return Ok(vec![(hgvs_weaver::data::IdentifierType::ProteinAccession, "NP_PLUS.1".to_string())]); }
            if symbol == "NM_MINUS.1" { return Ok(vec![(hgvs_weaver::data::IdentifierType::ProteinAccession, "NP_MINUS.1".to_string())]); }
        }
        Ok(vec![(hgvs_weaver::data::IdentifierType::Unknown, symbol.to_string())])
    }

    fn get_identifier_type(&self, _identifier: &str) -> Result<hgvs_weaver::data::IdentifierType, HgvsError> {
        Ok(hgvs_weaver::data::IdentifierType::Unknown)
    }

    fn c_to_g(&self, transcript_ac: &str, pos: TranscriptPos, offset: IntronicOffset) -> Result<(String, GenomicPos), HgvsError> {
        let tx = self.get_transcript(transcript_ac, None)?;
        Ok((tx.reference_accession().to_string(), GenomicPos(pos.0 + offset.0)))
    }
}

#[test]
fn test_toy_plus_strand_mapping() {
    let hdp = JsonDataProvider::new("../tests/data/toy_data.json");
    let mapper = VariantMapper::new(&hdp);

    // Genomic variant NC_TOY.1:g.25A>T
    let var_g = parse_hgvs_variant("NC_TOY.1:g.25A>T").unwrap();
    if let SequenceVariant::Genomic(v) = var_g {
        let var_c = mapper.g_to_c(&v, "NM_PLUS.1").unwrap();
        assert_eq!(var_c.to_string(), "NM_PLUS.1:c.1A>T");
    }

    let var_c = parse_hgvs_variant("NM_PLUS.1:c.1A>T").unwrap();
    if let SequenceVariant::Coding(v) = var_c {
        let var_p = mapper.c_to_p(&v, Some("NP_PLUS.1")).unwrap();
        assert_eq!(var_p.to_string(), "NP_PLUS.1:p.(Met1Leu)");
    }
}

#[test]
fn test_toy_minus_strand_mapping() {
    let hdp = JsonDataProvider::new("../tests/data/toy_data.json");
    let mapper = VariantMapper::new(&hdp);

    let var_g = parse_hgvs_variant("NC_TOY.1:g.236T>G").unwrap();
    if let SequenceVariant::Genomic(v) = var_g {
        let var_c = mapper.g_to_c(&v, "NM_MINUS.1").unwrap();
        assert_eq!(var_c.to_string(), "NM_MINUS.1:c.32A>C");
    }

    // Test c. to g. on minus strand
    let var_c = parse_hgvs_variant("NM_MINUS.1:c.32A>C").unwrap();
    if let SequenceVariant::Coding(v) = var_c {
        let var_g = mapper.c_to_g(&v, "NC_TOY.1").unwrap();
        assert_eq!(var_g.to_string(), "NC_TOY.1:g.236T>G");
    }
}
