use hgvs_weaver::data::{DataProvider, ExonData, IdentifierKind, IdentifierType, TranscriptData};
use hgvs_weaver::structs::{GenomicPos, IntronicOffset, TranscriptPos};
use hgvs_weaver::{parse_hgvs_variant, HgvsError, SequenceVariant, Transcript, VariantMapper};

struct MockDataProvider {
    transcripts: std::collections::HashMap<String, TranscriptData>,
    sequences: std::collections::HashMap<String, String>,
}

impl DataProvider for MockDataProvider {
    fn get_transcript(
        &self,
        ac: &str,
        _ref_ac: Option<&str>,
    ) -> Result<Box<dyn Transcript>, HgvsError> {
        self.transcripts
            .get(ac)
            .cloned()
            .map(|tx| Box::new(tx) as Box<dyn Transcript>)
            .ok_or_else(|| HgvsError::ValidationError(format!("Transcript {} not found", ac)))
    }

    fn get_seq(
        &self,
        ac: &str,
        start: i32,
        end: i32,
        _kind: IdentifierType,
    ) -> Result<String, HgvsError> {
        let seq = self
            .sequences
            .get(ac)
            .ok_or_else(|| HgvsError::ValidationError(format!("Sequence {} not found", ac)))?;
        let start = start as usize;
        let end = if end == -1 { seq.len() } else { end as usize };
        if start > seq.len() || end > seq.len() || start > end {
            return Err(HgvsError::ValidationError(format!(
                "Invalid seq range: {}-{} for len {}",
                start,
                end,
                seq.len()
            )));
        }
        Ok(seq[start..end].to_string())
    }

    fn get_symbol_accessions(
        &self,
        _symbol: &str,
        _src: IdentifierKind,
        _target: IdentifierKind,
    ) -> Result<Vec<(IdentifierType, String)>, HgvsError> {
        Ok(vec![])
    }
    fn get_identifier_type(&self, _identifier: &str) -> Result<IdentifierType, HgvsError> {
        Ok(IdentifierType::TranscriptAccession)
    }
    fn c_to_g(
        &self,
        _ac: &str,
        _pos: TranscriptPos,
        _offset: IntronicOffset,
    ) -> Result<(String, GenomicPos), HgvsError> {
        Ok(("NC_000001.1".to_string(), GenomicPos(0)))
    }
}

fn run_regression_test(
    dp: MockDataProvider,
    variant_str: &str,
    protein_ac: &str,
    expected_p: &str,
) {
    let mapper = VariantMapper::new(&dp);
    let variant = parse_hgvs_variant(variant_str).unwrap();
    if let SequenceVariant::Coding(cv) = variant {
        // We use the variant as-is to match ClinVar representation exactly.
        // Normalization can sometimes shift complex variants in ways that change the translation representation.
        let p = mapper.c_to_p(&cv, Some(protein_ac)).expect("c_to_p failed");
        let expected = format!("{}:{}", protein_ac, expected_p);
        assert_eq!(p.to_string(), expected, "Variant: {}", variant_str);
    } else {
        panic!("Parsed variant is not Coding");
    }
}

#[test]
fn test_regression_nm_000038_6() {
    // APC c.1972_1973delinsAT -> p.(Glu658Met)
    let seq = include_str!("data/NM_000038.6.seq").trim();
    let ac = "NM_000038.6";
    let tx = TranscriptData {
        ac: ac.to_string(),
        gene: "APC".to_string(),
        cds_start_index: Some(TranscriptPos(59)),
        cds_end_index: Some(TranscriptPos(8590)),
        strand: 1,
        reference_accession: "NC_000017.11".to_string(),
        exons: vec![
            ExonData {
                transcript_start: TranscriptPos(0),
                transcript_end: TranscriptPos(41),
                reference_start: GenomicPos(112737884),
                reference_end: GenomicPos(112737924),
                alt_strand: 1,
                cigar: "41=".to_string(),
            },
            ExonData {
                transcript_start: TranscriptPos(41),
                transcript_end: TranscriptPos(194),
                reference_start: GenomicPos(112754872),
                reference_end: GenomicPos(112755024),
                alt_strand: 1,
                cigar: "153=".to_string(),
            },
            ExonData {
                transcript_start: TranscriptPos(194),
                transcript_end: TranscriptPos(279),
                reference_start: GenomicPos(112766325),
                reference_end: GenomicPos(112766409),
                alt_strand: 1,
                cigar: "85=".to_string(),
            },
            ExonData {
                transcript_start: TranscriptPos(279),
                transcript_end: TranscriptPos(481),
                reference_start: GenomicPos(112767188),
                reference_end: GenomicPos(112767389),
                alt_strand: 1,
                cigar: "202=".to_string(),
            },
            ExonData {
                transcript_start: TranscriptPos(481),
                transcript_end: TranscriptPos(590),
                reference_start: GenomicPos(112775628),
                reference_end: GenomicPos(112775736),
                alt_strand: 1,
                cigar: "109=".to_string(),
            },
            ExonData {
                transcript_start: TranscriptPos(590),
                transcript_end: TranscriptPos(704),
                reference_start: GenomicPos(112780789),
                reference_end: GenomicPos(112780902),
                alt_strand: 1,
                cigar: "114=".to_string(),
            },
            ExonData {
                transcript_start: TranscriptPos(704),
                transcript_end: TranscriptPos(788),
                reference_start: GenomicPos(112792445),
                reference_end: GenomicPos(112792528),
                alt_strand: 1,
                cigar: "84=".to_string(),
            },
            ExonData {
                transcript_start: TranscriptPos(788),
                transcript_end: TranscriptPos(893),
                reference_start: GenomicPos(112801278),
                reference_end: GenomicPos(112801382),
                alt_strand: 1,
                cigar: "105=".to_string(),
            },
            ExonData {
                transcript_start: TranscriptPos(893),
                transcript_end: TranscriptPos(992),
                reference_start: GenomicPos(112815494),
                reference_end: GenomicPos(112815592),
                alt_strand: 1,
                cigar: "99=".to_string(),
            },
            ExonData {
                transcript_start: TranscriptPos(992),
                transcript_end: TranscriptPos(1371),
                reference_start: GenomicPos(112818965),
                reference_end: GenomicPos(112819343),
                alt_strand: 1,
                cigar: "379=".to_string(),
            },
            ExonData {
                transcript_start: TranscriptPos(1371),
                transcript_end: TranscriptPos(1467),
                reference_start: GenomicPos(112821895),
                reference_end: GenomicPos(112821990),
                alt_strand: 1,
                cigar: "96=".to_string(),
            },
            ExonData {
                transcript_start: TranscriptPos(1467),
                transcript_end: TranscriptPos(1607),
                reference_start: GenomicPos(112827107),
                reference_end: GenomicPos(112827246),
                alt_strand: 1,
                cigar: "140=".to_string(),
            },
            ExonData {
                transcript_start: TranscriptPos(1607),
                transcript_end: TranscriptPos(1685),
                reference_start: GenomicPos(112827928),
                reference_end: GenomicPos(112828005),
                alt_strand: 1,
                cigar: "78=".to_string(),
            },
            ExonData {
                transcript_start: TranscriptPos(1685),
                transcript_end: TranscriptPos(1802),
                reference_start: GenomicPos(112828855),
                reference_end: GenomicPos(112828971),
                alt_strand: 1,
                cigar: "117=".to_string(),
            },
            ExonData {
                transcript_start: TranscriptPos(1802),
                transcript_end: TranscriptPos(2017),
                reference_start: GenomicPos(112834950),
                reference_end: GenomicPos(112835164),
                alt_strand: 1,
                cigar: "215=".to_string(),
            },
            ExonData {
                transcript_start: TranscriptPos(2017),
                transcript_end: TranscriptPos(10704),
                reference_start: GenomicPos(112837552),
                reference_end: GenomicPos(112846238),
                alt_strand: 1,
                cigar: "8687=".to_string(),
            },
        ],
    };
    let mut transcripts = std::collections::HashMap::new();
    transcripts.insert(ac.to_string(), tx);
    let mut sequences = std::collections::HashMap::new();
    sequences.insert(ac.to_string(), seq.to_string());

    let dp = MockDataProvider {
        transcripts,
        sequences,
    };
    run_regression_test(
        dp,
        "NM_000038.6:c.1972_1973delinsAT",
        "NP_000029.2",
        "p.(Glu658Met)",
    );
}

#[test]
fn test_regression_nm_000527_5() {
    // LDLR c.669_683delinsAACTGCGGTAAACTGCGGTAAACT -> p.(Asp224_Glu228delinsThrAlaValAsnCysGlyLysLeu)
    let seq = include_str!("data/NM_000527.5.seq").trim();
    let ac = "NM_000527.5";
    let tx = TranscriptData {
        ac: ac.to_string(),
        gene: "LDLR".to_string(),
        cds_start_index: Some(TranscriptPos(86)),
        cds_end_index: Some(TranscriptPos(2668)),
        strand: 1,
        reference_accession: "NC_000019.10".to_string(),
        exons: vec![
            ExonData {
                transcript_start: TranscriptPos(0),
                transcript_end: TranscriptPos(153),
                reference_start: GenomicPos(11089462),
                reference_end: GenomicPos(11089614),
                alt_strand: 1,
                cigar: "153=".to_string(),
            },
            ExonData {
                transcript_start: TranscriptPos(153),
                transcript_end: TranscriptPos(276),
                reference_start: GenomicPos(11100222),
                reference_end: GenomicPos(11100344),
                alt_strand: 1,
                cigar: "123=".to_string(),
            },
            ExonData {
                transcript_start: TranscriptPos(276),
                transcript_end: TranscriptPos(399),
                reference_start: GenomicPos(11102663),
                reference_end: GenomicPos(11102785),
                alt_strand: 1,
                cigar: "123=".to_string(),
            },
            ExonData {
                transcript_start: TranscriptPos(399),
                transcript_end: TranscriptPos(780),
                reference_start: GenomicPos(11105219),
                reference_end: GenomicPos(11105599),
                alt_strand: 1,
                cigar: "381=".to_string(),
            },
            ExonData {
                transcript_start: TranscriptPos(780),
                transcript_end: TranscriptPos(903),
                reference_start: GenomicPos(11106564),
                reference_end: GenomicPos(11106686),
                alt_strand: 1,
                cigar: "123=".to_string(),
            },
            ExonData {
                transcript_start: TranscriptPos(903),
                transcript_end: TranscriptPos(1026),
                reference_start: GenomicPos(11107391),
                reference_end: GenomicPos(11107513),
                alt_strand: 1,
                cigar: "123=".to_string(),
            },
            ExonData {
                transcript_start: TranscriptPos(1026),
                transcript_end: TranscriptPos(1146),
                reference_start: GenomicPos(11110651),
                reference_end: GenomicPos(11110770),
                alt_strand: 1,
                cigar: "120=".to_string(),
            },
            ExonData {
                transcript_start: TranscriptPos(1146),
                transcript_end: TranscriptPos(1272),
                reference_start: GenomicPos(11111513),
                reference_end: GenomicPos(11111638),
                alt_strand: 1,
                cigar: "126=".to_string(),
            },
            ExonData {
                transcript_start: TranscriptPos(1272),
                transcript_end: TranscriptPos(1444),
                reference_start: GenomicPos(11113277),
                reference_end: GenomicPos(11113448),
                alt_strand: 1,
                cigar: "172=".to_string(),
            },
            ExonData {
                transcript_start: TranscriptPos(1444),
                transcript_end: TranscriptPos(1672),
                reference_start: GenomicPos(11113534),
                reference_end: GenomicPos(11113761),
                alt_strand: 1,
                cigar: "228=".to_string(),
            },
            ExonData {
                transcript_start: TranscriptPos(1672),
                transcript_end: TranscriptPos(1791),
                reference_start: GenomicPos(11116093),
                reference_end: GenomicPos(11116211),
                alt_strand: 1,
                cigar: "119=".to_string(),
            },
            ExonData {
                transcript_start: TranscriptPos(1791),
                transcript_end: TranscriptPos(1931),
                reference_start: GenomicPos(11116858),
                reference_end: GenomicPos(11116997),
                alt_strand: 1,
                cigar: "140=".to_string(),
            },
            ExonData {
                transcript_start: TranscriptPos(1931),
                transcript_end: TranscriptPos(2073),
                reference_start: GenomicPos(11120091),
                reference_end: GenomicPos(11120232),
                alt_strand: 1,
                cigar: "142=".to_string(),
            },
            ExonData {
                transcript_start: TranscriptPos(2073),
                transcript_end: TranscriptPos(2226),
                reference_start: GenomicPos(11120369),
                reference_end: GenomicPos(11120521),
                alt_strand: 1,
                cigar: "153=".to_string(),
            },
            ExonData {
                transcript_start: TranscriptPos(2226),
                transcript_end: TranscriptPos(2397),
                reference_start: GenomicPos(11123173),
                reference_end: GenomicPos(11123343),
                alt_strand: 1,
                cigar: "171=".to_string(),
            },
            ExonData {
                transcript_start: TranscriptPos(2397),
                transcript_end: TranscriptPos(2475),
                reference_start: GenomicPos(11128007),
                reference_end: GenomicPos(11128084),
                alt_strand: 1,
                cigar: "78=".to_string(),
            },
            ExonData {
                transcript_start: TranscriptPos(2475),
                transcript_end: TranscriptPos(2633),
                reference_start: GenomicPos(11129512),
                reference_end: GenomicPos(11129669),
                alt_strand: 1,
                cigar: "158=".to_string(),
            },
            ExonData {
                transcript_start: TranscriptPos(2633),
                transcript_end: TranscriptPos(5173),
                reference_start: GenomicPos(11131280),
                reference_end: GenomicPos(11133819),
                alt_strand: 1,
                cigar: "2540=".to_string(),
            },
        ],
    };
    let mut transcripts = std::collections::HashMap::new();
    transcripts.insert(ac.to_string(), tx);
    let mut sequences = std::collections::HashMap::new();
    sequences.insert(ac.to_string(), seq.to_string());

    let dp = MockDataProvider {
        transcripts,
        sequences,
    };
    run_regression_test(
        dp,
        "NM_000527.5:c.669_683delinsAACTGCGGTAAACTGCGGTAAACT",
        "NP_000518.1",
        "p.(Asp224_Glu228delinsThrAlaValAsnCysGlyLysLeu)",
    );
}

#[test]
fn test_regression_nm_000478_6() {
    // ALPL c.876_882delinsT -> p.(Gly293_Asp294del)
    let seq = include_str!("data/NM_000478.6.seq").trim();
    let ac = "NM_000478.6";
    let tx = TranscriptData {
        ac: ac.to_string(),
        gene: "ALPL".to_string(),
        cds_start_index: Some(TranscriptPos(199)),
        cds_end_index: Some(TranscriptPos(1773)),
        strand: 1,
        reference_accession: "NC_000001.11".to_string(),
        exons: vec![
            ExonData {
                transcript_start: TranscriptPos(0),
                transcript_end: TranscriptPos(95),
                reference_start: GenomicPos(21509422),
                reference_end: GenomicPos(21509516),
                alt_strand: 1,
                cigar: "95=".to_string(),
            },
            ExonData {
                transcript_start: TranscriptPos(95),
                transcript_end: TranscriptPos(260),
                reference_start: GenomicPos(21553977),
                reference_end: GenomicPos(21554141),
                alt_strand: 1,
                cigar: "165=".to_string(),
            },
            ExonData {
                transcript_start: TranscriptPos(260),
                transcript_end: TranscriptPos(380),
                reference_start: GenomicPos(21560625),
                reference_end: GenomicPos(21560744),
                alt_strand: 1,
                cigar: "120=".to_string(),
            },
            ExonData {
                transcript_start: TranscriptPos(380),
                transcript_end: TranscriptPos(496),
                reference_start: GenomicPos(21561096),
                reference_end: GenomicPos(21561211),
                alt_strand: 1,
                cigar: "116=".to_string(),
            },
            ExonData {
                transcript_start: TranscriptPos(496),
                transcript_end: TranscriptPos(671),
                reference_start: GenomicPos(21563109),
                reference_end: GenomicPos(21563283),
                alt_strand: 1,
                cigar: "175=".to_string(),
            },
            ExonData {
                transcript_start: TranscriptPos(671),
                transcript_end: TranscriptPos(847),
                reference_start: GenomicPos(21564040),
                reference_end: GenomicPos(21564215),
                alt_strand: 1,
                cigar: "176=".to_string(),
            },
            ExonData {
                transcript_start: TranscriptPos(847),
                transcript_end: TranscriptPos(991),
                reference_start: GenomicPos(21568103),
                reference_end: GenomicPos(21568246),
                alt_strand: 1,
                cigar: "144=".to_string(),
            },
            ExonData {
                transcript_start: TranscriptPos(991),
                transcript_end: TranscriptPos(1061),
                reference_start: GenomicPos(21570304),
                reference_end: GenomicPos(21570373),
                alt_strand: 1,
                cigar: "70=".to_string(),
            },
            ExonData {
                transcript_start: TranscriptPos(1061),
                transcript_end: TranscriptPos(1196),
                reference_start: GenomicPos(21573664),
                reference_end: GenomicPos(21573798),
                alt_strand: 1,
                cigar: "135=".to_string(),
            },
            ExonData {
                transcript_start: TranscriptPos(1196),
                transcript_end: TranscriptPos(1388),
                reference_start: GenomicPos(21575732),
                reference_end: GenomicPos(21575923),
                alt_strand: 1,
                cigar: "192=".to_string(),
            },
            ExonData {
                transcript_start: TranscriptPos(1388),
                transcript_end: TranscriptPos(1508),
                reference_start: GenomicPos(21576521),
                reference_end: GenomicPos(21576640),
                alt_strand: 1,
                cigar: "120=".to_string(),
            },
            ExonData {
                transcript_start: TranscriptPos(1508),
                transcript_end: TranscriptPos(2536),
                reference_start: GenomicPos(21577382),
                reference_end: GenomicPos(21578409),
                alt_strand: 1,
                cigar: "1028=".to_string(),
            },
        ],
    };
    let mut transcripts = std::collections::HashMap::new();
    transcripts.insert(ac.to_string(), tx);
    let mut sequences = std::collections::HashMap::new();
    sequences.insert(ac.to_string(), seq.to_string());

    let dp = MockDataProvider {
        transcripts,
        sequences,
    };
    run_regression_test(
        dp,
        "NM_000478.6:c.876_882delinsT",
        "NP_000469.3",
        "p.(Gly293_Asp294del)",
    );
}

#[test]
fn test_regression_nm_001122606_1() {
    // CASP8 c.763_768delinsTGAAGT -> p.(Asn255Ter)
    // Note: ClinVar expects p.(Asn255_Pro256delinsTer), but Weaver simplifies this to p.(Asn255Ter)
    // which is biologically equivalent.
    let seq = include_str!("data/NM_001122606.1.seq").trim();
    let ac = "NM_001122606.1";
    let tx = TranscriptData {
        ac: ac.to_string(),
        gene: "CASP8".to_string(),
        cds_start_index: Some(TranscriptPos(180)),
        cds_end_index: Some(TranscriptPos(1415)),
        strand: 1,
        reference_accession: "NC_000002.12".to_string(),
        exons: vec![
            ExonData {
                transcript_start: TranscriptPos(0),
                transcript_end: TranscriptPos(244),
                reference_start: GenomicPos(120469105),
                reference_end: GenomicPos(120469348),
                alt_strand: -1,
                cigar: "244=".to_string(),
            },
            ExonData {
                transcript_start: TranscriptPos(244),
                transcript_end: TranscriptPos(363),
                reference_start: GenomicPos(120456650),
                reference_end: GenomicPos(120456768),
                alt_strand: -1,
                cigar: "119=".to_string(),
            },
            ExonData {
                transcript_start: TranscriptPos(363),
                transcript_end: TranscriptPos(577),
                reference_start: GenomicPos(120455356),
                reference_end: GenomicPos(120455569),
                alt_strand: -1,
                cigar: "214=".to_string(),
            },
            ExonData {
                transcript_start: TranscriptPos(577),
                transcript_end: TranscriptPos(736),
                reference_start: GenomicPos(120448969),
                reference_end: GenomicPos(120449127),
                alt_strand: -1,
                cigar: "159=".to_string(),
            },
            ExonData {
                transcript_start: TranscriptPos(736),
                transcript_end: TranscriptPos(921),
                reference_start: GenomicPos(120447840),
                reference_end: GenomicPos(120448024),
                alt_strand: -1,
                cigar: "185=".to_string(),
            },
            ExonData {
                transcript_start: TranscriptPos(921),
                transcript_end: TranscriptPos(1044),
                reference_start: GenomicPos(120446304),
                reference_end: GenomicPos(120446426),
                alt_strand: -1,
                cigar: "123=".to_string(),
            },
            ExonData {
                transcript_start: TranscriptPos(1044),
                transcript_end: TranscriptPos(1108),
                reference_start: GenomicPos(120442598),
                reference_end: GenomicPos(120442661),
                alt_strand: -1,
                cigar: "64=".to_string(),
            },
            ExonData {
                transcript_start: TranscriptPos(1108),
                transcript_end: TranscriptPos(1273),
                reference_start: GenomicPos(120441729),
                reference_end: GenomicPos(120441893),
                alt_strand: -1,
                cigar: "165=".to_string(),
            },
            ExonData {
                transcript_start: TranscriptPos(1273),
                transcript_end: TranscriptPos(3752),
                reference_start: GenomicPos(120426147),
                reference_end: GenomicPos(120428625),
                alt_strand: -1,
                cigar: "2479=".to_string(),
            },
        ],
    };
    let mut transcripts = std::collections::HashMap::new();
    transcripts.insert(ac.to_string(), tx);
    let mut sequences = std::collections::HashMap::new();
    sequences.insert(ac.to_string(), seq.to_string());

    let dp = MockDataProvider {
        transcripts,
        sequences,
    };
    run_regression_test(
        dp,
        "NM_001122606.1:c.763_768delinsTGAAGT",
        "NP_001116078.1",
        "p.(Asn255Ter)",
    );
}

#[test]
fn test_regression_nm_000465_4() {
    // ABCA4 c.66_70delinsTGCGT -> p.(Pro24Ser)
    let seq = include_str!("data/NM_000465.4.seq").trim();
    let ac = "NM_000465.4";
    let tx = TranscriptData {
        ac: ac.to_string(),
        gene: "ABCA4".to_string(),
        cds_start_index: Some(TranscriptPos(114)),
        cds_end_index: Some(TranscriptPos(2447)),
        strand: 1,
        reference_accession: "NC_000001.11".to_string(),
        exons: vec![
            ExonData {
                transcript_start: TranscriptPos(0),
                transcript_end: TranscriptPos(272),
                reference_start: GenomicPos(214809411),
                reference_end: GenomicPos(214809682),
                alt_strand: -1,
                cigar: "272=".to_string(),
            },
            ExonData {
                transcript_start: TranscriptPos(272),
                transcript_end: TranscriptPos(329),
                reference_start: GenomicPos(214797060),
                reference_end: GenomicPos(214797116),
                alt_strand: -1,
                cigar: "57=".to_string(),
            },
            ExonData {
                transcript_start: TranscriptPos(329),
                transcript_end: TranscriptPos(478),
                reference_start: GenomicPos(214792296),
                reference_end: GenomicPos(214792444),
                alt_strand: -1,
                cigar: "149=".to_string(),
            },
            ExonData {
                transcript_start: TranscriptPos(478),
                transcript_end: TranscriptPos(1428),
                reference_start: GenomicPos(214780559),
                reference_end: GenomicPos(214781508),
                alt_strand: -1,
                cigar: "950=".to_string(),
            },
            ExonData {
                transcript_start: TranscriptPos(1428),
                transcript_end: TranscriptPos(1509),
                reference_start: GenomicPos(214769231),
                reference_end: GenomicPos(214769311),
                alt_strand: -1,
                cigar: "81=".to_string(),
            },
            ExonData {
                transcript_start: TranscriptPos(1509),
                transcript_end: TranscriptPos(1682),
                reference_start: GenomicPos(214767481),
                reference_end: GenomicPos(214767653),
                alt_strand: -1,
                cigar: "173=".to_string(),
            },
            ExonData {
                transcript_start: TranscriptPos(1682),
                transcript_end: TranscriptPos(1791),
                reference_start: GenomicPos(214752446),
                reference_end: GenomicPos(214752554),
                alt_strand: -1,
                cigar: "109=".to_string(),
            },
            ExonData {
                transcript_start: TranscriptPos(1791),
                transcript_end: TranscriptPos(1924),
                reference_start: GenomicPos(214745721),
                reference_end: GenomicPos(214745853),
                alt_strand: -1,
                cigar: "133=".to_string(),
            },
            ExonData {
                transcript_start: TranscriptPos(1924),
                transcript_end: TranscriptPos(2017),
                reference_start: GenomicPos(214745066),
                reference_end: GenomicPos(214745158),
                alt_strand: -1,
                cigar: "93=".to_string(),
            },
            ExonData {
                transcript_start: TranscriptPos(2017),
                transcript_end: TranscriptPos(2115),
                reference_start: GenomicPos(214730410),
                reference_end: GenomicPos(214730507),
                alt_strand: -1,
                cigar: "98=".to_string(),
            },
            ExonData {
                transcript_start: TranscriptPos(2115),
                transcript_end: TranscriptPos(5478),
                reference_start: GenomicPos(214725645),
                reference_end: GenomicPos(214729007),
                alt_strand: -1,
                cigar: "3363=".to_string(),
            },
        ],
    };
    let mut transcripts = std::collections::HashMap::new();
    transcripts.insert(ac.to_string(), tx);
    let mut sequences = std::collections::HashMap::new();
    sequences.insert(ac.to_string(), seq.to_string());

    let dp = MockDataProvider {
        transcripts,
        sequences,
    };
    run_regression_test(
        dp,
        "NM_000465.4:c.66_70delinsTGCGT",
        "NP_000456.2",
        "p.(Pro24Ser)",
    );
}
