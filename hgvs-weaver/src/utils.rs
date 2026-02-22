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
        let res = match codon.to_uppercase().as_str() {
            "TTT" | "TTC" => 'F',
            "TTA" | "TTG" => 'L',
            "CTT" | "CTC" | "CTA" | "CTG" => 'L',
            "ATT" | "ATC" | "ATA" => 'I',
            "ATG" => 'M',
            "GTT" | "GTC" | "GTA" | "GTG" => 'V',
            "TCT" | "TCC" | "TCA" | "TCG" => 'S',
            "CCT" | "CCC" | "CCA" | "CCG" => 'P',
            "ACT" | "ACC" | "ACA" | "ACG" => 'T',
            "GCT" | "GCC" | "GCA" | "GCG" => 'A',
            "TAT" | "TAC" => 'Y',
            "TAA" | "TAG" | "TGA" => '*',
            "CAT" | "CAC" => 'H',
            "CAA" | "CAG" => 'Q',
            "AAT" | "AAC" => 'N',
            "AAA" | "AAG" => 'K',
            "GAT" | "GAC" => 'D',
            "GAA" | "GAG" => 'E',
            "TGT" | "TGC" => 'C',
            "TGG" => 'W',
            "CGT" | "CGC" | "CGA" | "CGG" => 'R',
            "AGT" | "AGC" => 'S',
            "AGA" | "AGG" => 'R',
            "GGT" | "GGC" | "GGA" | "GGG" => 'G',
            _ => 'X',
        };
        aa.push(res);
        if res == '*' {
            break;
        }
    }
    aa
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
            let res = match chunk.to_lowercase().as_str() {
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
            };

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
