use hgvs_weaver::{parse_hgvs_variant, HgvsParser, Rule};
use pest::Parser;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

#[test]
fn test_gauntlet_variants() {
    let gauntlet_path = Path::new("../tests/data/gauntlet");
    let file = File::open(gauntlet_path).expect("Failed to open gauntlet file");
    let reader = BufReader::new(file);

    // List of variants that are syntactically valid (parser accepts) but semantically invalid (validator rejects)
    // and thus fail our strict parsing.
    let known_strict_failures = [
        "AC_01234.5:g.1_22A>T", // Deletion sequence length 1 != interval length 22
    ];

    for (line_num, line) in reader.lines().enumerate() {
        let line = line.expect("Failed to read line");
        let trimmed = line.trim();

        if trimmed.is_empty() || trimmed.starts_with("#") || trimmed.starts_with("!") {
            continue;
        }

        if known_strict_failures.contains(&trimmed) {
            println!(
                "Skipping known strict failure line {}: {}",
                line_num + 1,
                trimmed
            );
            continue;
        }

        match parse_hgvs_variant(trimmed) {
            Ok(v) => {
                let rt = v.to_string();
                // Some variants might normalize, so we check if the round-tripped string parses to the same variant.
                let v2 = parse_hgvs_variant(&rt)
                    .unwrap_or_else(|_| panic!("Round-tripped string failed to parse: {}", rt));
                assert_eq!(v, v2, "Round trip structure mismatch for {}", trimmed);
            }
            Err(e) => match HgvsParser::parse(Rule::hgvs_variant, trimmed) {
                Ok(_) => panic!(
                    "Line {}: '{}' - Raw parser succeeded, AST conversion failed: {}",
                    line_num + 1,
                    trimmed,
                    e
                ),
                Err(pe) => panic!(
                    "Line {}: '{}' - Parsing failed: {}",
                    line_num + 1,
                    trimmed,
                    pe
                ),
            },
        }
    }
}

#[test]
fn test_reject_variants() {
    let reject_path = Path::new("../tests/data/reject");
    if !reject_path.exists() {
        return;
    }
    let file = File::open(reject_path).expect("Failed to open reject file");
    let reader = BufReader::new(file);

    for (line_num, line) in reader.lines().enumerate() {
        let line = line.expect("Failed to read line");
        let trimmed = line.trim();

        if trimmed.is_empty() || trimmed.starts_with("#") {
            continue;
        }

        let variant = trimmed.split('\t').next().unwrap();

        let result = parse_hgvs_variant(variant);
        assert!(
            result.is_err(),
            "Line {}: Variant '{}' should have been rejected",
            line_num + 1,
            variant
        );
    }
}

#[test]
fn test_gene_names() {
    let variants = vec![
        "NM_01234.5(BOGUS):c.22+1A>T",
        "NM_01234.5(BOGUS-EXCELLENT):c.22+1A>T",
        "NM_01234.5(BOGUS_EXCELLENT):c.22+1A>T",
        "BOGUS:c.22+1A>T",
        "BOGUS-EXCELLENT:c.22+1A>T",
    ];

    for v in variants {
        assert!(
            parse_hgvs_variant(v).is_ok(),
            "Failed to parse variant with gene name: {}",
            v
        );
    }
}

#[test]
fn test_invalid_gene_names() {
    let variants = vec![
        "NM_01234.5(1BOGUS):c.22+1A>T",
        "NM_01234.5(-BOGUS):c.22+1A>T",
        "NM_01234.5(BOGUS-):c.22+1A>T",
    ];

    for v in variants {
        assert!(
            parse_hgvs_variant(v).is_err(),
            "Variant with invalid gene name should have failed: {}",
            v
        );
    }
}
