use crate::error::HgvsError;
use std::fmt::{self, Display, Formatter};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Residue {
    Ala,
    Arg,
    Asn,
    Asp,
    Cys,
    Gln,
    Glu,
    Gly,
    His,
    Ile,
    Leu,
    Lys,
    Met,
    Phe,
    Pro,
    Ser,
    Thr,
    Trp,
    Tyr,
    Val,
    // Ambiguous / Special
    Asx, // B
    Glx, // Z
    Xaa, // X
    Sec, // U
    Ter, // *
}

impl Display for Residue {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        let s = match self {
            Residue::Ala => "A",
            Residue::Arg => "R",
            Residue::Asn => "N",
            Residue::Asp => "D",
            Residue::Cys => "C",
            Residue::Gln => "Q",
            Residue::Glu => "E",
            Residue::Gly => "G",
            Residue::His => "H",
            Residue::Ile => "I",
            Residue::Leu => "L",
            Residue::Lys => "K",
            Residue::Met => "M",
            Residue::Phe => "F",
            Residue::Pro => "P",
            Residue::Ser => "S",
            Residue::Thr => "T",
            Residue::Trp => "W",
            Residue::Tyr => "Y",
            Residue::Val => "V",
            Residue::Asx => "B",
            Residue::Glx => "Z",
            Residue::Xaa => "X",
            Residue::Sec => "U",
            Residue::Ter => "*",
        };
        write!(f, "{}", s)
    }
}

pub fn reverse_complement(seq: &str) -> String {
    seq.chars().rev().map(complement_dna_char).collect()
}

pub fn complement_dna_char(c: char) -> char {
    match c {
        'A' => 'T',
        'T' => 'A',
        'C' => 'G',
        'G' => 'C',
        'N' => 'N',
        'a' => 't',
        't' => 'a',
        'c' => 'g',
        'g' => 'c',
        'n' => 'n',
        'U' => 'A',
        'u' => 'a',
        _ => c,
    }
}

pub fn aa1_to_aa3(aa1: char) -> &'static str {
    match aa1.to_ascii_uppercase() {
        'A' => "Ala",
        'R' => "Arg",
        'N' => "Asn",
        'D' => "Asp",
        'C' => "Cys",
        'E' => "Glu",
        'Q' => "Gln",
        'G' => "Gly",
        'H' => "His",
        'I' => "Ile",
        'L' => "Leu",
        'K' => "Lys",
        'M' => "Met",
        'F' => "Phe",
        'P' => "Pro",
        'S' => "Ser",
        'T' => "Thr",
        'W' => "Trp",
        'Y' => "Tyr",
        'V' => "Val",
        '*' => "Ter",
        'X' => "Xaa",
        _ => "Xaa",
    }
}

pub fn seq1_to_aa3(seq1: &str) -> String {
    seq1.chars().map(aa1_to_aa3).collect()
}

pub fn translate_cds(cds: &str) -> String {
    let mut aa = String::new();
    for i in (0..cds.len()).step_by(3) {
        if i + 3 > cds.len() {
            break;
        }
        let codon = &cds[i..i + 3];
        let res = translate_codon(codon).unwrap_or('X');
        aa.push(res);
        if res == '*' {
            break;
        }
    }
    aa
}

/// Translates a single 3-base codon to a 1-letter amino acid.
pub fn translate_codon(codon: &str) -> Option<char> {
    match codon.to_uppercase().as_str() {
        "TTT" | "TTC" => Some('F'),
        "TTA" | "TTG" => Some('L'),
        "CTT" | "CTC" | "CTA" | "CTG" => Some('L'),
        "ATT" | "ATC" | "ATA" => Some('I'),
        "ATG" => Some('M'),
        "GTT" | "GTC" | "GTA" | "GTG" => Some('V'),
        "TCT" | "TCC" | "TCA" | "TCG" => Some('S'),
        "CCT" | "CCC" | "CCA" | "CCG" => Some('P'),
        "ACT" | "ACC" | "ACA" | "ACG" => Some('T'),
        "GCT" | "GCC" | "GCA" | "GCG" => Some('A'),
        "TAT" | "TAC" => Some('Y'),
        "TAA" | "TAG" | "TGA" => Some('*'),
        "CAT" | "CAC" => Some('H'),
        "CAA" | "CAG" => Some('Q'),
        "AAT" | "AAC" => Some('N'),
        "AAA" | "AAG" => Some('K'),
        "GAT" | "GAC" => Some('D'),
        "GAA" | "GAG" => Some('E'),
        "TGT" | "TGC" => Some('C'),
        "TGG" => Some('W'),
        "CGT" | "CGC" | "CGA" | "CGG" => Some('R'),
        "AGT" | "AGC" => Some('S'),
        "AGA" | "AGG" => Some('R'),
        "GGT" | "GGC" | "GGA" | "GGG" => Some('G'),
        _ => None,
    }
}

/// Returns all DNA codons that encode the given 1-letter amino acid.
pub fn codons_for_aa(aa: char) -> Vec<&'static str> {
    match aa {
        'F' => vec!["TTT", "TTC"],
        'L' => vec!["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
        'I' => vec!["ATT", "ATC", "ATA"],
        'M' => vec!["ATG"],
        'V' => vec!["GTT", "GTC", "GTA", "GTG"],
        'S' => vec!["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
        'P' => vec!["CCT", "CCC", "CCA", "CCG"],
        'T' => vec!["ACT", "ACC", "ACA", "ACG"],
        'A' => vec!["GCT", "GCC", "GCA", "GCG"],
        'Y' => vec!["TAT", "TAC"],
        '*' => vec!["TAA", "TAG", "TGA"],
        'H' => vec!["CAT", "CAC"],
        'Q' => vec!["CAA", "CAG"],
        'N' => vec!["AAT", "AAC"],
        'K' => vec!["AAA", "AAG"],
        'D' => vec!["GAT", "GAC"],
        'E' => vec!["GAA", "GAG"],
        'C' => vec!["TGT", "TGC"],
        'W' => vec!["TGG"],
        'R' => vec!["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
        'G' => vec!["GGT", "GGC", "GGA", "GGG"],
        _ => vec![],
    }
}

pub fn aa3_to_aa1(aa3: &str) -> String {
    match aa3.to_lowercase().as_str() {
        "ala" => "A".to_string(),
        "arg" => "R".to_string(),
        "asn" => "N".to_string(),
        "asp" => "D".to_string(),
        "cys" => "C".to_string(),
        "gln" => "Q".to_string(),
        "glu" => "E".to_string(),
        "gly" => "G".to_string(),
        "his" => "H".to_string(),
        "ile" => "I".to_string(),
        "leu" => "L".to_string(),
        "lys" => "K".to_string(),
        "met" => "M".to_string(),
        "phe" => "F".to_string(),
        "pro" => "P".to_string(),
        "ser" => "S".to_string(),
        "thr" => "T".to_string(),
        "trp" => "W".to_string(),
        "tyr" => "Y".to_string(),
        "val" => "V".to_string(),
        "asx" => "B".to_string(),
        "glx" => "Z".to_string(),
        "xaa" => "X".to_string(),
        "ter" | "stop" | "*" => "*".to_string(),
        // 1-letter codes are returned as-is (uppercase)
        s if s.len() == 1 => s.to_uppercase(),
        _ => "X".to_string(),
    }
}

pub fn normalize_aa(s: &str) -> String {
    aa3_to_aa1(s)
}

fn aa3_chunk_to_residue(s: &str) -> Option<Residue> {
    match s.to_lowercase().as_str() {
        "ala" => Some(Residue::Ala),
        "arg" => Some(Residue::Arg),
        "asn" => Some(Residue::Asn),
        "asp" => Some(Residue::Asp),
        "cys" => Some(Residue::Cys),
        "gln" => Some(Residue::Gln),
        "glu" => Some(Residue::Glu),
        "gly" => Some(Residue::Gly),
        "his" => Some(Residue::His),
        "ile" => Some(Residue::Ile),
        "leu" => Some(Residue::Leu),
        "lys" => Some(Residue::Lys),
        "met" => Some(Residue::Met),
        "phe" => Some(Residue::Phe),
        "pro" => Some(Residue::Pro),
        "ser" => Some(Residue::Ser),
        "thr" => Some(Residue::Thr),
        "trp" => Some(Residue::Trp),
        "tyr" => Some(Residue::Tyr),
        "val" => Some(Residue::Val),
        "asx" => Some(Residue::Asx),
        "glx" => Some(Residue::Glx),
        "xaa" => Some(Residue::Xaa),
        "ter" | "stop" => Some(Residue::Ter),
        "sec" => Some(Residue::Sec),
        _ => None,
    }
}

pub fn decompose_aa(s: &str) -> Result<Vec<Residue>, HgvsError> {
    if s.is_empty() {
        return Ok(Vec::new());
    }

    // Attempt 1: Try to parse as only 3-letter codes
    let mut res_3 = Vec::new();
    let chars: Vec<char> = s.chars().collect();
    let mut i = 0;
    let mut all_3 = true;

    while i < chars.len() {
        if i + 3 <= chars.len() {
            let chunk: String = chars[i..i + 3].iter().collect();
            let res = aa3_chunk_to_residue(&chunk);

            if let Some(r) = res {
                res_3.push(r);
                i += 3;
                continue;
            }
        }
        all_3 = false;
        break;
    }

    if all_3 && i == chars.len() {
        return Ok(res_3);
    }

    // Attempt 2: Try to parse as only 1-letter codes
    let mut res_1 = Vec::new();
    let mut all_1 = true;
    for c in s.chars() {
        let res = match c {
            'A' => Some(Residue::Ala),
            'R' => Some(Residue::Arg),
            'N' => Some(Residue::Asn),
            'D' => Some(Residue::Asp),
            'C' => Some(Residue::Cys),
            'Q' => Some(Residue::Gln),
            'E' => Some(Residue::Glu),
            'G' => Some(Residue::Gly),
            'H' => Some(Residue::His),
            'I' => Some(Residue::Ile),
            'L' => Some(Residue::Leu),
            'K' => Some(Residue::Lys),
            'M' => Some(Residue::Met),
            'F' => Some(Residue::Phe),
            'P' => Some(Residue::Pro),
            'S' => Some(Residue::Ser),
            'T' => Some(Residue::Thr),
            'W' => Some(Residue::Trp),
            'Y' => Some(Residue::Tyr),
            'V' => Some(Residue::Val),
            'B' => Some(Residue::Asx),
            'Z' => Some(Residue::Glx),
            'X' => Some(Residue::Xaa),
            'U' => Some(Residue::Sec),
            '*' => Some(Residue::Ter),
            _ => None,
        };

        if let Some(r) = res {
            res_1.push(r);
        } else {
            all_1 = false;
            break;
        }
    }

    if all_1 {
        return Ok(res_1);
    }

    Err(HgvsError::Other(format!(
        "Could not decompose protein sequence '{}' into exclusively 3-letter or 1-letter codes",
        s
    )))
}
