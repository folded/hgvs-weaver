use crate::error::HgvsError;
use crate::data::{DataProvider, IdentifierType};

/// Trait for a sequence of bases (DNA, RNA, or Protein).
pub trait Sequence {
    /// Returns an iterator over the bases.
    fn iter(&self) -> Box<dyn Iterator<Item = char> + '_>;

    /// Returns the length of the sequence.
    fn len(&self) -> usize;

    /// Returns true if the sequence is empty.
    fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Fetches a subsequence as a MemSequence.
    fn subseq(&self, start: usize, end: usize) -> Result<MemSequence, HgvsError> {
        if start > self.len() || end > self.len() || start > end {
            return Err(HgvsError::ValidationError(format!("Invalid subsequence range {}..{} for length {}", start, end, self.len())));
        }
        let seq: String = self.iter().skip(start).take(end - start).collect();
        Ok(MemSequence(seq))
    }

    /// Returns a lazy slice of the sequence.
    fn slice(&self, start: usize, end: usize) -> SliceSequence<'_> where Self: Sized {
        SliceSequence { inner: self, start, end }
    }

    /// Converts the sequence back to a String.
    fn to_string(&self) -> String {
        self.iter().collect()
    }
}

/// An in-memory sequence.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MemSequence(pub String);

impl MemSequence {
    pub fn new(seq: String) -> Self {
        MemSequence(seq)
    }

    pub fn as_str(&self) -> &str {
        &self.0
    }
}

impl Sequence for MemSequence {
    fn iter(&self) -> Box<dyn Iterator<Item = char> + '_> {
        Box::new(self.0.chars())
    }

    fn len(&self) -> usize {
        self.0.len()
    }
}

/// A sequence that is fetched lazily from a data provider.
pub struct LazySequence<'a> {
    pub hdp: &'a dyn DataProvider,
    pub ac: String,
    pub start: usize, // 0-based
    pub end: usize,   // 0-based, exclusive
    pub kind: IdentifierType,
}

impl<'a> Sequence for LazySequence<'a> {
    fn iter(&self) -> Box<dyn Iterator<Item = char> + '_> {
        let seq = self.hdp.get_seq(&self.ac, self.start as i32, self.end as i32, self.kind).unwrap_or_default();
        Box::new(seq.chars().collect::<Vec<_>>().into_iter()) // Still collecting to string/vec for now to avoid borrow checker hell with provider
    }

    fn len(&self) -> usize {
        self.end - self.start
    }
}

/// Adapter for reverse-complementation.
pub struct RevCompSequence<'a> {
    pub inner: &'a dyn Sequence,
}

impl<'a> Sequence for RevCompSequence<'a> {
    fn iter(&self) -> Box<dyn Iterator<Item = char> + '_> {
        Box::new(self.inner.iter().collect::<Vec<_>>().into_iter().rev().map(|c| match c {
            'A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C', 'N' => 'N',
            'a' => 't', 't' => 'a', 'c' => 'g', 'g' => 'c', 'n' => 'n',
            'U' => 'A', 'u' => 'a',
            _ => c
        }))
    }

    fn len(&self) -> usize {
        self.inner.len()
    }
}

/// Adapter for transcription (T -> U).
pub struct TranscribedSequence<'a> {
    pub inner: &'a dyn Sequence,
}

impl<'a> Sequence for TranscribedSequence<'a> {
    fn iter(&self) -> Box<dyn Iterator<Item = char> + '_> {
        Box::new(self.inner.iter().map(|c| match c {
            'T' => 'U', 't' => 'u',
            _ => c
        }))
    }

    fn len(&self) -> usize {
        self.inner.len()
    }
}

/// Adapter for splicing multiple sequences together.
pub struct SplicedSequence<'a> {
    pub pieces: Vec<&'a dyn Sequence>,
}

impl<'a> Sequence for SplicedSequence<'a> {
    fn iter(&self) -> Box<dyn Iterator<Item = char> + '_> {
        Box::new(self.pieces.iter().flat_map(|p| p.iter()))
    }

    fn len(&self) -> usize {
        self.pieces.iter().map(|p| p.len()).sum()
    }
}

/// Adapter for translation (Nucleotides -> Amino Acids).
pub struct TranslatedSequence<'a> {
    pub inner: &'a dyn Sequence,
}

impl<'a> Sequence for TranslatedSequence<'a> {
    fn iter(&self) -> Box<dyn Iterator<Item = char> + '_> {
        Box::new(TranslateIterator::new(self.inner.iter()))
    }

    fn len(&self) -> usize {
        self.inner.len() / 3
    }
}

struct TranslateIterator<'a> {
    inner: Box<dyn Iterator<Item = char> + 'a>,
    done: bool,
}

impl<'a> TranslateIterator<'a> {
    fn new(inner: Box<dyn Iterator<Item = char> + 'a>) -> Self {
        TranslateIterator { inner, done: false }
    }
}

impl<'a> Iterator for TranslateIterator<'a> {
    type Item = char;

    fn next(&mut self) -> Option<Self::Item> {
        if self.done { return None; }
        
        let c1 = self.inner.next()?;
        let c2 = self.inner.next()?;
        let c3 = self.inner.next()?;
        
        let codon = [c1, c2, c3];
        let aa = match codon {
            ['T', 'T', 'T'] | ['T', 'T', 'C'] | ['U', 'U', 'U'] | ['U', 'U', 'C'] => 'F',
            ['T', 'T', 'A'] | ['T', 'T', 'G'] | ['U', 'U', 'A'] | ['U', 'U', 'G'] | 
            ['C', 'T', 'T'] | ['C', 'T', 'C'] | ['C', 'T', 'A'] | ['C', 'T', 'G'] | 
            ['C', 'U', 'T'] | ['C', 'U', 'C'] | ['C', 'U', 'A'] | ['C', 'U', 'G'] => 'L',
            ['A', 'T', 'T'] | ['A', 'T', 'C'] | ['A', 'T', 'A'] | ['A', 'U', 'T'] | ['A', 'U', 'C'] | ['A', 'U', 'A'] => 'I',
            ['A', 'T', 'G'] | ['A', 'U', 'G'] => 'M',
            ['G', 'T', 'T'] | ['G', 'T', 'C'] | ['G', 'T', 'A'] | ['G', 'T', 'G'] | 
            ['G', 'U', 'T'] | ['G', 'U', 'C'] | ['G', 'U', 'A'] | ['G', 'U', 'G'] => 'V',
            ['T', 'C', 'T'] | ['T', 'C', 'C'] | ['T', 'C', 'A'] | ['T', 'C', 'G'] | 
            ['U', 'C', 'T'] | ['U', 'C', 'C'] | ['U', 'C', 'A'] | ['U', 'C', 'G'] | 
            ['A', 'G', 'T'] | ['A', 'G', 'C'] | ['A', 'G', 'U'] => 'S',
            ['C', 'C', 'T'] | ['C', 'C', 'C'] | ['C', 'C', 'A'] | ['C', 'C', 'G'] | ['C', 'C', 'U'] => 'P', 
            ['A', 'C', 'T'] | ['A', 'C', 'C'] | ['A', 'C', 'A'] | ['A', 'C', 'G'] | ['A', 'C', 'U'] => 'T',
            ['G', 'C', 'T'] | ['G', 'C', 'C'] | ['G', 'C', 'A'] | ['G', 'C', 'G'] | ['G', 'C', 'U'] => 'A',
            ['T', 'A', 'T'] | ['T', 'A', 'C'] | ['U', 'A', 'U'] | ['U', 'A', 'C'] => 'Y',
            ['T', 'A', 'A'] | ['T', 'A', 'G'] | ['T', 'G', 'A'] | ['U', 'A', 'A'] | ['U', 'A', 'G'] | ['U', 'G', 'A'] => '*',
            ['C', 'A', 'T'] | ['C', 'A', 'C'] | ['C', 'A', 'U'] => 'H', 
            ['C', 'A', 'A'] | ['C', 'A', 'G'] => 'Q',
            ['A', 'A', 'T'] | ['A', 'A', 'C'] | ['A', 'A', 'U'] => 'N', 
            ['A', 'A', 'A'] | ['A', 'A', 'G'] => 'K',
            ['G', 'A', 'T'] | ['G', 'A', 'C'] | ['G', 'A', 'U'] => 'D', 
            ['G', 'A', 'A'] | ['G', 'A', 'G'] => 'E',
            ['T', 'G', 'T'] | ['T', 'G', 'C'] | ['U', 'G', 'T'] | ['U', 'G', 'C'] => 'C', 
            ['T', 'G', 'G'] | ['U', 'G', 'G'] => 'W',
            ['C', 'G', 'T'] | ['C', 'G', 'C'] | ['C', 'G', 'A'] | ['C', 'G', 'G'] | ['C', 'G', 'U'] | 
            ['A', 'G', 'A'] | ['A', 'G', 'G'] => 'R', 
            ['G', 'G', 'T'] | ['G', 'G', 'C'] | ['G', 'G', 'A'] | ['G', 'G', 'G'] | ['G', 'G', 'U'] => 'G',
            _ => 'X',
        };
        
        if aa == '*' { self.done = true; }
        Some(aa)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mem_sequence() {
        let s = MemSequence::new("ATGC".to_string());
        assert_eq!(s.len(), 4);
        assert_eq!(s.to_string(), "ATGC");
    }

    #[test]
    fn test_rev_comp() {
        let s = MemSequence::new("ATGC".to_string());
        let rc = RevCompSequence { inner: &s };
        assert_eq!(rc.to_string(), "GCAT");
    }

    #[test]
    fn test_translate() {
        let s = MemSequence::new("ATGGCTTAA".to_string());
        let t = TranslatedSequence { inner: &s };
        assert_eq!(t.to_string(), "MA*");
    }

    #[test]
    fn test_spliced_sequence() {
        let s1 = MemSequence::new("ATG".to_string());
        let s2 = MemSequence::new("GCT".to_string());
        let s3 = MemSequence::new("TAA".to_string());
        let spliced = SplicedSequence { pieces: vec![&s1, &s2, &s3] };
        assert_eq!(spliced.len(), 9);
        assert_eq!(spliced.to_string(), "ATGGCTTAA");
        
        let trans = TranslatedSequence { inner: &spliced };
        assert_eq!(trans.to_string(), "MA*");
    }

    #[test]
    fn test_slice_sequence() {
        let s = MemSequence::new("ATGGCTTAA".to_string());
        let slice = s.slice(3, 6);
        assert_eq!(slice.len(), 3);
        assert_eq!(slice.to_string(), "GCT");
    }
}

/// Adapter for slicing a sequence.
pub struct SliceSequence<'a> {
    pub inner: &'a dyn Sequence,
    pub start: usize,
    pub end: usize,
}

impl<'a> Sequence for SliceSequence<'a> {
    fn iter(&self) -> Box<dyn Iterator<Item = char> + '_> {
        Box::new(self.inner.iter().skip(self.start).take(self.end - self.start))
    }

    fn len(&self) -> usize {
        self.end - self.start
    }
}
