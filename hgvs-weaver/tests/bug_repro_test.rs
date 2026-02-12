use hgvs_weaver::*;
use hgvs_weaver::structs::{TranscriptPos, GenomicPos};
use hgvs_weaver::data::{ExonData, TranscriptData};

struct BugMockDataProvider;

impl DataProvider for BugMockDataProvider {
    fn get_seq(&self, _ac: &str, _start: i32, _end: i32, _kind: hgvs_weaver::data::IdentifierType) -> Result<String, HgvsError> {
        let mut s = String::new();
        for _ in 0..10000 { s.push_str("AAA"); } // 30000 A's
        let mut full_seq = "A".repeat(30000);
        
        // Codon 1575: 4722, 4723, 4724 (0-based)
        full_seq.replace_range(4722..4725, "AGA"); 
        
        // Codon 310: 927..930
        full_seq.replace_range(927..930, "TCA");
        
        Ok(full_seq)
    }

    fn get_transcript(&self, transcript_ac: &str, _reference_ac: Option<&str>) -> Result<Box<dyn Transcript>, HgvsError> {
        let exons = vec![
            ExonData {
                transcript_start: TranscriptPos(0),
                transcript_end: TranscriptPos(30000),
                reference_start: GenomicPos(1000),
                reference_end: GenomicPos(31000),
                alt_strand: 1,
                cigar: "30000M".to_string(),
            }
        ];
        
        let td = TranscriptData {
            ac: transcript_ac.to_string(),
            gene: "BUG".to_string(),
            cds_start_index: Some(TranscriptPos(0)),
            cds_end_index: Some(TranscriptPos(20000)),
            strand: 1,
            reference_accession: "NC_0001.10".to_string(),
            exons
        };
        Ok(Box::new(td))
    }

    fn get_symbol_accessions(&self, _symbol: &str, _sk: hgvs_weaver::data::IdentifierKind, _tk: hgvs_weaver::data::IdentifierKind) -> Result<Vec<(hgvs_weaver::data::IdentifierType, String)>, HgvsError> {
        Ok(vec![(hgvs_weaver::data::IdentifierType::ProteinAccession, "NP001".to_string())])
    }

    fn get_identifier_type(&self, _identifier: &str) -> Result<hgvs_weaver::data::IdentifierType, HgvsError> {
        Ok(hgvs_weaver::data::IdentifierType::Unknown)
    }
}

#[test]
fn test_repro_redundant_identity_interval() {
    let hdp = BugMockDataProvider;
    let mapper = VariantMapper::new(&hdp);
    
    let var_c = parse_hgvs_variant("NM_001039476.3:c.929_930=").unwrap();
    if let SequenceVariant::Coding(v) = var_c {
        let var_p = mapper.c_to_p(&v, None).unwrap();
        // Should be p.(Ser310=) NOT p.(Ser310_Ser310=)
        assert_eq!(var_p.to_string(), "NP001:p.(Ser310=)");
    }
}

#[test]
fn test_repro_delins_identity_mismatch() {
    let hdp = BugMockDataProvider;
    let mapper = VariantMapper::new(&hdp);

    let var_c = parse_hgvs_variant("NM_000051.4:c.4723_4725delinsAGA").unwrap();
    if let SequenceVariant::Coding(v) = var_c {
        let var_p = mapper.c_to_p(&v, None).unwrap();
        assert_eq!(var_p.to_string(), "NP001:p.(Arg1575=)");
    }
}

#[test]
fn test_repro_redundant_interval_multiple_codons() {
    let hdp = BugMockDataProvider;
    let mapper = VariantMapper::new(&hdp);
    
    let var_c = parse_hgvs_variant("NM_001008844.3:c.4498_4499delinsAT").unwrap();
    if let SequenceVariant::Coding(v) = var_c {
        let var_p = mapper.c_to_p(&v, None).unwrap();
        println!("DEBUG: NM_001008844.3:c.4498_4499delinsAT -> {}", var_p.to_string());
        // Should be Lys1500Ile (since AT is substation)
        assert!(!var_p.to_string().contains("_"));
    }
}

#[test]
fn test_repro_aa_mismatch_delins() {
    let hdp = BugMockDataProvider;
    let mapper = VariantMapper::new(&hdp);

    // NM_031935.3:c.15010_15011delinsTT
    // Codon 5004 is 15010, 15011, 15012.
    // Let's ensure 3rd base (15012) in BugMockDataProvider results in Phe (TTT).
    // Indices: 15009, 15010, 15011.
    // My BugMockDataProvider has all A's by default. AAA -> Lys.
    // If we replace 15009, 15010 with TT, and 15011 is A, we get TTA -> Leu.
    // GT says Phe (TTT).
    let var_c = parse_hgvs_variant("NM_031935.3:c.15010_15011delinsTT").unwrap();
    if let SequenceVariant::Coding(v) = var_c {
        let var_p = mapper.c_to_p(&v, None).unwrap();
        // Just print it for now to see what Weaver produces.
        println!("DEBUG: NM_031935.3:c.15010_15011delinsTT -> {}", var_p.to_string());
    }
}
