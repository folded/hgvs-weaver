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
            return Ok(vec![(hgvs_weaver::data::IdentifierType::ProteinAccession, "MOCK".to_string())]);
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

fn run_c_to_p(hgvsc: &str, expected_p: &str) {
    let hdp = JsonDataProvider::new("../tests/data/toy_data.json");
    let mapper = VariantMapper::new(&hdp);
    let var_c = parse_hgvs_variant(hgvsc).unwrap();
    if let SequenceVariant::Coding(v) = var_c {
        let var_p = mapper.c_to_p(&v, Some("MOCK")).unwrap();
        assert_eq!(var_p.to_string(), expected_p, "Conversion failed for {}", hgvsc);
    } else {
        panic!("Expected coding variant for {}", hgvsc);
    }
}

#[test]
fn test_parity_substitutions() {
    run_c_to_p("NM_999999.1:c.6A>G", "MOCK:p.(Lys2=)");
    run_c_to_p("NM_999999.1:c.6A>T", "MOCK:p.(Lys2Asn)");
    run_c_to_p("NM_999996.1:c.8C>A", "MOCK:p.(Ser3Ter)");
}

#[test]
fn test_parity_insertions() {
    run_c_to_p("NM_999999.1:c.6_7insGGG", "MOCK:p.(Lys2_Ala3insGly)");
    // Rust currently produces p.(Ala8ValfsTer?) for frameshifts with unknown stop.
    // The Python test says p.(Ala8ValfsTer?).
    run_c_to_p("NM_999999.1:c.22_23insT", "MOCK:p.(Ala8ValfsTer?)");
    run_c_to_p("NM_999999.1:c.8_9insTT", "MOCK:p.(Lys4Ter)");
}

#[test]
fn test_parity_deletions() {
    run_c_to_p("NM_999999.1:c.10_12del", "MOCK:p.(Lys4del)");
    run_c_to_p("NM_999999.1:c.4_15del", "MOCK:p.(Lys2_Ala5del)");
    run_c_to_p("NM_999995.1:c.4_6del", "MOCK:p.(Lys3del)");
    run_c_to_p("NM_999994.1:c.4_9del", "MOCK:p.(Lys3_Lys4del)");
    run_c_to_p("NM_999999.1:c.5_7del", "MOCK:p.(Lys2_Ala3delinsThr)");
    run_c_to_p("NM_999993.1:c.13_24del", "MOCK:p.(Arg5_Ala8del)");
}

#[test]
fn test_parity_frameshifts() {
    run_c_to_p("NM_999999.1:c.11_12del", "MOCK:p.(Lys4SerfsTer?)");
    run_c_to_p("NM_999997.1:c.7del", "MOCK:p.(Ala3ArgfsTer6)");
}

#[test]
fn test_parity_indels() {
    run_c_to_p("NM_999999.1:c.11_12delinsTCCCA", "MOCK:p.(Lys4delinsIlePro)");
    run_c_to_p("NM_999999.1:c.11_18delinsTCCCA", "MOCK:p.(Lys4_Phe6delinsIlePro)");
}
