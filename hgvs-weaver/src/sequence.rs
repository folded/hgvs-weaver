use crate::error::HgvsError;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Sequence(String);

impl Sequence {
    pub fn new(seq: String) -> Self {
        Sequence(seq)
    }

    pub fn as_str(&self) -> &str {
        &self.0
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    pub fn subseq(&self, start: usize, end: usize) -> Result<Self, HgvsError> {
        if start > self.len() || end > self.len() || start > end {
             // Clamping or error? Let's error for safety.
            return Err(HgvsError::ValidationError(format!("Invalid subsequence range {}..{} for length {}", start, end, self.len())));
        }
        Ok(Sequence(self.0[start..end].to_string()))
    }

    pub fn reverse_complement(&self) -> Self {
        let rev_comp: String = self.0.chars()
            .rev()
            .map(|c| match c {
                'A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C', 'N' => 'N',
                'a' => 't', 't' => 'a', 'c' => 'g', 'g' => 'c', 'n' => 'n',
                'U' => 'A', 'u' => 'a',
                _ => c
            })
            .collect();
        Sequence(rev_comp)
    }

    pub fn transcribe(&self) -> Self {
        // DNA to RNA (T -> U)
        let rna: String = self.0.chars()
            .map(|c| match c {
                'T' => 'U', 't' => 'u',
                _ => c
            })
            .collect();
        Sequence(rna)
    }

    pub fn translate(&self, cds_start: usize) -> Result<String, HgvsError> {
        // Simple translation. Assume self is DNA or RNA (we'll treat T and U equivalent for codon table).
        // cds_start is 0-based index relative to the start of this sequence.
        
        if cds_start >= self.len() {
             return Err(HgvsError::ValidationError(format!("CDS start {} out of bounds for length {}", cds_start, self.len())));
        }

        let seq_slice = &self.0[cds_start..];
        let mut aa = String::new();
        
        for i in (0..seq_slice.len()).step_by(3) {
            if i + 3 > seq_slice.len() { break; }
            let codon = &seq_slice[i..i+3];
            let res = match codon.to_uppercase().as_str() {
                "TTT" | "TTC" | "UUU" | "UUC" => 'F',
                "TTA" | "TTG" | "UUA" | "UUG" | "CTT" | "CTC" | "CTA" | "CTG" | "CUU" | "CUC" | "CUA" | "CUG" => 'L',
                "ATT" | "ATC" | "ATA" | "AUU" | "AUC" | "AUA" => 'I',
                "ATG" | "AUG" => 'M',
                "GTT" | "GTC" | "GTA" | "GTG" | "GUU" | "GUC" | "GUA" | "GUG" => 'V',
                "TCT" | "TCC" | "TCA" | "TCG" | "UCU" | "UCC" | "UCA" | "UCG" | "AGT" | "AGC" | "AGU" => 'S',
                "CCT" | "CCC" | "CCA" | "CCG" | "CCU" => 'P', 
                "ACT" | "ACC" | "ACA" | "ACG" | "ACU" => 'T',
                "GCT" | "GCC" | "GCA" | "GCG" | "GCU" => 'A',
                "TAT" | "TAC" | "UAU" | "UAC" => 'Y',
                "TAA" | "TAG" | "TGA" | "UAA" | "UAG" | "UGA" => '*',
                "CAT" | "CAC" | "CAU" => 'H', 
                "CAA" | "CAG" => 'Q',
                "AAT" | "AAC" | "AAU" => 'N', 
                "AAA" | "AAG" => 'K',
                "GAT" | "GAC" | "GAU" => 'D', 
                "GAA" | "GAG" => 'E',
                "TGT" | "TGC" | "UGU" | "UGC" => 'C', 
                "TGG" | "UGG" => 'W',
                "CGT" | "CGC" | "CGA" | "CGG" | "CGU" | "AGA" | "AGG" => 'R', 
                "GGT" | "GGC" | "GGA" | "GGG" | "GGU" => 'G',
                _ => 'X',
            };
            aa.push(res);
            if res == '*' { break; }
        }
        Ok(aa)
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reverse_complement() {
        let s = Sequence::new("ATGC".to_string());
        assert_eq!(s.reverse_complement().as_str(), "GCAT");
        
        let s2 = Sequence::new("atgc".to_string());
        assert_eq!(s2.reverse_complement().as_str(), "gcat");
    }

    #[test]
    fn test_translate() {
        // Met Ala Stop
        let s = Sequence::new("ATGGCTTAA".to_string());
        assert_eq!(s.translate(0).unwrap(), "MA*");
        
        // Offset
        let s2 = Sequence::new("CCATGGCTTAA".to_string());
        assert_eq!(s2.translate(2).unwrap(), "MA*");
    }
}
