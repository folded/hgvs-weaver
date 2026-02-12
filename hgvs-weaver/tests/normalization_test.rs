use hgvs_weaver::*;
use hgvs_weaver::structs::{TranscriptPos, GenomicPos};
use hgvs_weaver::data::{ExonData, TranscriptData};

struct NormMockDataProvider;

impl DataProvider for NormMockDataProvider {
    fn get_seq(&self, _ac: &str, _start: i32, _end: i32, _kind: hgvs_weaver::data::IdentifierType) -> Result<String, HgvsError> {
        let mut s = String::new();
        s.push_str("AAAAAAAAAA"); // 10 A's
        // NM_0001.1: ATG AAA TAG (9 bases)
        s.push_str("ATGAAATAG");
        for _ in 0..20 {
            s.push_str("ATGC");
        }
        Ok(s)
    }

    fn get_transcript(&self, transcript_ac: &str, _reference_ac: Option<&str>) -> Result<Box<dyn Transcript>, HgvsError> {
        let exons = vec![
            ExonData {
                transcript_start: TranscriptPos(0),
                transcript_end: TranscriptPos(100),
                reference_start: GenomicPos(1000),
                reference_end: GenomicPos(1100),
                alt_strand: 1,
                cigar: "100M".to_string(),
            }
        ];
        
        let (cds_start, cds_end) = match transcript_ac {
            "NM_0001.1" => (10, 19), // Met Lys *
            _ => return Err(HgvsError::DataProviderError("Transcript not found".to_string())),
        };

        let td = TranscriptData {
            ac: transcript_ac.to_string(),
            gene: "NORM".to_string(),
            cds_start_index: Some(TranscriptPos(cds_start)),
            cds_end_index: Some(TranscriptPos(cds_end)),
            strand: 1,
            reference_accession: "NC_0001.10".to_string(),
            exons
        };
        Ok(Box::new(td))
    }

    fn get_symbol_accessions(&self, symbol: &str, _sk: hgvs_weaver::data::IdentifierKind, tk: hgvs_weaver::data::IdentifierKind) -> Result<Vec<(hgvs_weaver::data::IdentifierType, String)>, HgvsError> {
        if tk == hgvs_weaver::data::IdentifierKind::Protein && symbol == "NM_0001.1" {
            return Ok(vec![(hgvs_weaver::data::IdentifierType::ProteinAccession, "NP_0001.1".to_string())]);
        }
        Ok(vec![(hgvs_weaver::data::IdentifierType::Unknown, symbol.to_string())])
    }

    fn get_identifier_type(&self, _identifier: &str) -> Result<hgvs_weaver::data::IdentifierType, HgvsError> {
        Ok(hgvs_weaver::data::IdentifierType::Unknown)
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
        // Should be ext*X
        assert!(var_p.to_string().contains("ext*"));
    }
}
