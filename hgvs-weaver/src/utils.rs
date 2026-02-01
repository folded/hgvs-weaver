pub fn reverse_complement(seq: &str) -> String {
    seq.chars().rev().map(complement_dna_char).collect()
}

pub fn complement_dna_char(c: char) -> char {
    match c {
        'A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C', 'N' => 'N',
        'a' => 't', 't' => 'a', 'c' => 'g', 'g' => 'c', 'n' => 'n',
        'U' => 'A', 'u' => 'a',
        _ => c
    }
}

pub fn aa1_to_aa3(aa1: char) -> &'static str {
    match aa1.to_ascii_uppercase() {
        'A' => "Ala", 'R' => "Arg", 'N' => "Asn", 'D' => "Asp", 'C' => "Cys",
        'E' => "Glu", 'Q' => "Gln", 'G' => "Gly", 'H' => "His", 'I' => "Ile",
        'L' => "Leu", 'K' => "Lys", 'M' => "Met", 'F' => "Phe", 'P' => "Pro",
        'S' => "Ser", 'T' => "Thr", 'W' => "Trp", 'Y' => "Tyr", 'V' => "Val",
        '*' => "Ter", 'X' => "Xaa", _ => "Xaa"
    }
}

pub fn seq1_to_aa3(seq1: &str) -> String {
    seq1.chars().map(aa1_to_aa3).collect()
}

pub fn translate_cds(cds: &str) -> String {
    let mut aa = String::new();
    for i in (0..cds.len()).step_by(3) {
        if i + 3 > cds.len() { break; }
        let codon = &cds[i..i+3];
        let res = match codon.to_uppercase().as_str() {
            "TTT" | "TTC" => 'F', "TTA" | "TTG" => 'L',
            "CTT" | "CTC" | "CTA" | "CTG" => 'L',
            "ATT" | "ATC" | "ATA" => 'I', "ATG" => 'M',
            "GTT" | "GTC" | "GTA" | "GTG" => 'V',
            "TCT" | "TCC" | "TCA" | "TCG" => 'S',
            "CCT" | "CCC" | "CCA" | "CCG" => 'P',
            "ACT" | "ACC" | "ACA" | "ACG" => 'T',
            "GCT" | "GCC" | "GCA" | "GCG" => 'A',
            "TAT" | "TAC" => 'Y', "TAA" | "TAG" | "TGA" => '*',
            "CAT" | "CAC" => 'H', "CAA" | "CAG" => 'Q',
            "AAT" | "AAC" => 'N', "AAA" | "AAG" => 'K',
            "GAT" | "GAC" => 'D', "GAA" | "GAG" => 'E',
            "TGT" | "TGC" => 'C', "TGG" => 'W',
            "CGT" | "CGC" | "CGA" | "CGG" => 'R',
            "AGT" | "AGC" => 'S', "AGA" | "AGG" => 'R',
            "GGT" | "GGC" | "GGA" | "GGG" => 'G',
            _ => 'X',
        };
        aa.push(res);
        if res == '*' { break; }
    }
    aa
}
