use hgvs_weaver::utils::{decompose_aa, Residue};

#[test]
fn test_strict_decompose() {
    // 3-letter codes
    assert_eq!(decompose_aa("AlaCysAsp").unwrap(), vec![Residue::Ala, Residue::Cys, Residue::Asp]);
    
    // 1-letter codes
    assert_eq!(decompose_aa("ACD").unwrap(), vec![Residue::Ala, Residue::Cys, Residue::Asp]);
    
    // Mixed - should fail now
    assert!(decompose_aa("AlaCD").is_err());
    assert!(decompose_aa("ACysD").is_err());
    
    // Ambiguous
    // "ADCYSK" -> 3-letter attempt finds "CYS" but fails on others?
    // Wait, my logic: 
    // Attempt 1 (3-letter): "ADC" invalid 3-letter -> fails all_3.
    // Attempt 2 (1-letter): "A", "D", "C", "Y", "S", "K" all valid 1-letter -> succeeds.
    // So "ADCYSK" should parse as 1-letter codes now!
    // Before it greedily took "CYS". Now strict 3-letter loop will see "ADC" is not 3-letter and bail.
    assert_eq!(decompose_aa("ADCYSK").unwrap(), vec![Residue::Ala, Residue::Asp, Residue::Cys, Residue::Tyr, Residue::Ser, Residue::Lys]);
    
    // "AlaCys" -> 3-letter: Ala, Cys -> OK
    assert_eq!(decompose_aa("AlaCys").unwrap(), vec![Residue::Ala, Residue::Cys]);
    
    // "AC" -> 1-letter: A, C -> OK
     assert_eq!(decompose_aa("AC").unwrap(), vec![Residue::Ala, Residue::Cys]);
}
