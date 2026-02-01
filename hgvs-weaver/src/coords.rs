use std::ops::{Add, Sub};
use serde::{Serialize, Deserialize};

/// HGVS coordinate system anchor point.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum Anchor {
    /// Start of the transcript sequence.
    TranscriptStart,
    /// First nucleotide of the start codon (A of ATG).
    CdsStart,
    /// Last nucleotide of the stop codon.
    CdsEnd,
}

/// Internal 0-based genomic position.
#[derive(Debug, PartialEq, Eq, PartialOrd, Ord, Clone, Copy, Serialize, Deserialize)]
pub struct GenomicPos(pub i32);

impl Add<i32> for GenomicPos {
    type Output = Self;
    fn add(self, rhs: i32) -> Self { GenomicPos(self.0 + rhs) }
}
impl Sub<i32> for GenomicPos {
    type Output = Self;
    fn sub(self, rhs: i32) -> Self { GenomicPos(self.0 - rhs) }
}
impl Sub<GenomicPos> for GenomicPos {
    type Output = i32;
    fn sub(self, rhs: GenomicPos) -> i32 { self.0 - rhs.0 }
}

impl GenomicPos {
    /// Converts the internal 0-based index to a 1-based `HgvsGenomicPos`.
    pub fn to_hgvs(&self) -> HgvsGenomicPos {
        (*self).into()
    }
}

/// Internal 0-based transcript position relative to sequence start.
#[derive(Debug, PartialEq, Eq, PartialOrd, Ord, Clone, Copy, Serialize, Deserialize)]
pub struct TranscriptPos(pub i32);

impl Add<i32> for TranscriptPos {
    type Output = Self;
    fn add(self, rhs: i32) -> Self { TranscriptPos(self.0 + rhs) }
}
impl Sub<i32> for TranscriptPos {
    type Output = Self;
    fn sub(self, rhs: i32) -> Self { TranscriptPos(self.0 - rhs) }
}
impl Sub<TranscriptPos> for TranscriptPos {
    type Output = i32;
    fn sub(self, rhs: TranscriptPos) -> i32 { self.0 - rhs.0 }
}

impl TranscriptPos {
    /// Converts the internal 0-based index to a 1-based `HgvsTranscriptPos`.
    pub fn to_hgvs(&self) -> HgvsTranscriptPos {
        (*self).into()
    }
}

/// Internal 0-based amino acid position.
#[derive(Debug, PartialEq, Eq, PartialOrd, Ord, Clone, Copy, Serialize, Deserialize)]
pub struct ProteinPos(pub i32);

impl Add<i32> for ProteinPos {
    type Output = Self;
    fn add(self, rhs: i32) -> Self { ProteinPos(self.0 + rhs) }
}
impl Sub<i32> for ProteinPos {
    type Output = Self;
    fn sub(self, rhs: i32) -> Self { ProteinPos(self.0 - rhs) }
}
impl Sub<ProteinPos> for ProteinPos {
    type Output = i32;
    fn sub(self, rhs: ProteinPos) -> i32 { self.0 - rhs.0 }
}

impl ProteinPos {
    /// Converts the internal 0-based index to a 1-based `HgvsProteinPos`.
    pub fn to_hgvs(&self) -> HgvsProteinPos {
        (*self).into()
    }
}

impl From<ProteinPos> for HgvsProteinPos {
    fn from(pos: ProteinPos) -> Self {
        HgvsProteinPos(pos.0 + 1)
    }
}

/// Distance from an exonic anchor into an intron.
#[derive(Debug, PartialEq, Eq, PartialOrd, Ord, Clone, Copy, Serialize, Deserialize)]
pub struct IntronicOffset(pub i32);

impl Add<i32> for IntronicOffset {
    type Output = Self;
    fn add(self, rhs: i32) -> Self { IntronicOffset(self.0 + rhs) }
}
impl Sub<i32> for IntronicOffset {
    type Output = Self;
    fn sub(self, rhs: i32) -> Self { IntronicOffset(self.0 - rhs) }
}
impl Sub<IntronicOffset> for IntronicOffset {
    type Output = i32;
    fn sub(self, rhs: IntronicOffset) -> i32 { self.0 - rhs.0 }
}

// --- 1-based HGVS tagged types ---

/// 1-based HGVS genomic coordinate.
#[derive(Debug, PartialEq, Eq, Clone, Copy, Serialize, Deserialize)]
pub struct HgvsGenomicPos(pub i32);
impl HgvsGenomicPos {
    /// Converts a 1-based HGVS genomic coordinate to an internal 0-based `GenomicPos`.
    pub fn to_index(&self) -> GenomicPos {
        (*self).into()
    }
}

impl From<HgvsGenomicPos> for GenomicPos {
    fn from(pos: HgvsGenomicPos) -> Self {
        GenomicPos(pos.0 - 1)
    }
}

/// 1-based HGVS transcript coordinate (c. or n.).
#[derive(Debug, PartialEq, Eq, Clone, Copy, Serialize, Deserialize)]
pub struct HgvsTranscriptPos(pub i32);
impl HgvsTranscriptPos {
    /// Converts a 1-based HGVS transcript coordinate to a 0-based internal sequence index.
    /// Correctly handles the non-existent position 0 in HGVS cDNA coordinates.
    pub fn to_index(&self) -> TranscriptPos {
        if self.0 > 0 { TranscriptPos(self.0 - 1) } else { TranscriptPos(self.0) }
    }
}

impl From<HgvsTranscriptPos> for TranscriptPos {
    fn from(pos: HgvsTranscriptPos) -> Self {
        pos.to_index()
    }
}

/// 1-based HGVS protein coordinate.
#[derive(Debug, PartialEq, Eq, Clone, Copy, Serialize, Deserialize)]
pub struct HgvsProteinPos(pub i32);
impl HgvsProteinPos {
    /// Converts a 1-based HGVS protein coordinate to an internal 0-based `ProteinPos`.
    pub fn to_index(&self) -> ProteinPos {
        (*self).into()
    }
}

impl From<HgvsProteinPos> for ProteinPos {
    fn from(pos: HgvsProteinPos) -> Self {
        ProteinPos(pos.0 - 1)
    }
}

use crate::structs::{GVariant, CVariant, PVariant, MVariant, NVariant, RVariant, Variant};

/// A complete HGVS variant spanning any coordinate system.
#[derive(Debug, PartialEq, Clone, Serialize, Deserialize)]
#[serde(tag = "variant_type")]
pub enum SequenceVariant {
    Genomic(GVariant),
    Coding(CVariant),
    Protein(PVariant),
    Mitochondrial(MVariant),
    NonCoding(NVariant),
    Rna(RVariant),
}

impl Variant for SequenceVariant {
    fn ac(&self) -> &str {
        match self {
            SequenceVariant::Genomic(v) => v.ac(),
            SequenceVariant::Coding(v) => v.ac(),
            SequenceVariant::Protein(v) => v.ac(),
            SequenceVariant::Mitochondrial(v) => v.ac(),
            SequenceVariant::NonCoding(v) => v.ac(),
            SequenceVariant::Rna(v) => v.ac(),
        }
    }
    fn gene(&self) -> Option<&str> {
        match self {
            SequenceVariant::Genomic(v) => v.gene(),
            SequenceVariant::Coding(v) => v.gene(),
            SequenceVariant::Protein(v) => v.gene(),
            SequenceVariant::Mitochondrial(v) => v.gene(),
            SequenceVariant::NonCoding(v) => v.gene(),
            SequenceVariant::Rna(v) => v.gene(),
        }
    }
    fn coordinate_type(&self) -> &str {
        match self {
            SequenceVariant::Genomic(v) => v.coordinate_type(),
            SequenceVariant::Coding(v) => v.coordinate_type(),
            SequenceVariant::Protein(v) => v.coordinate_type(),
            SequenceVariant::Mitochondrial(v) => v.coordinate_type(),
            SequenceVariant::NonCoding(v) => v.coordinate_type(),
            SequenceVariant::Rna(v) => v.coordinate_type(),
        }
    }
}

impl From<TranscriptPos> for HgvsTranscriptPos {
    fn from(pos: TranscriptPos) -> Self {
        if pos.0 >= 0 { HgvsTranscriptPos(pos.0 + 1) } else { HgvsTranscriptPos(pos.0) }
    }
}

impl From<GenomicPos> for HgvsGenomicPos {
    fn from(pos: GenomicPos) -> Self {
        HgvsGenomicPos(pos.0 + 1)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_genomic_conversions() {
        // 0-based -> 1-based HGVS
        assert_eq!(GenomicPos(0).to_hgvs(), HgvsGenomicPos(1));
        assert_eq!(GenomicPos(999).to_hgvs(), HgvsGenomicPos(1000));
        assert_eq!(HgvsGenomicPos::from(GenomicPos(42)), HgvsGenomicPos(43));

        // 1-based HGVS -> 0-based
        assert_eq!(HgvsGenomicPos(1).to_index(), GenomicPos(0));
        assert_eq!(HgvsGenomicPos(1000).to_index(), GenomicPos(999));
        assert_eq!(GenomicPos::from(HgvsGenomicPos(43)), GenomicPos(42));
    }

    #[test]
    fn test_transcript_conversions_jump() {
        // HGVS cDNA coordinates skip 0.
        // ..., c.-2, c.-1, c.1, c.2, ...

        // 0-based -> 1-based HGVS
        assert_eq!(TranscriptPos(-2).to_hgvs(), HgvsTranscriptPos(-2));
        assert_eq!(TranscriptPos(-1).to_hgvs(), HgvsTranscriptPos(-1));
        assert_eq!(TranscriptPos(0).to_hgvs(), HgvsTranscriptPos(1));
        assert_eq!(TranscriptPos(1).to_hgvs(), HgvsTranscriptPos(2));
        assert_eq!(HgvsTranscriptPos::from(TranscriptPos(0)), HgvsTranscriptPos(1));
        assert_eq!(HgvsTranscriptPos::from(TranscriptPos(-1)), HgvsTranscriptPos(-1));

        // 1-based HGVS -> 0-based
        assert_eq!(HgvsTranscriptPos(-2).to_index(), TranscriptPos(-2));
        assert_eq!(HgvsTranscriptPos(-1).to_index(), TranscriptPos(-1));
        assert_eq!(HgvsTranscriptPos(1).to_index(), TranscriptPos(0));
        assert_eq!(HgvsTranscriptPos(2).to_index(), TranscriptPos(1));
        assert_eq!(TranscriptPos::from(HgvsTranscriptPos(1)), TranscriptPos(0));
        assert_eq!(TranscriptPos::from(HgvsTranscriptPos(-1)), TranscriptPos(-1));
    }

    #[test]
    fn test_protein_conversions() {
        // 0-based -> 1-based HGVS
        assert_eq!(ProteinPos(0).to_hgvs(), HgvsProteinPos(1));
        assert_eq!(ProteinPos(10).to_hgvs(), HgvsProteinPos(11));
        assert_eq!(HgvsProteinPos::from(ProteinPos(5)), HgvsProteinPos(6));

        // 1-based HGVS -> 0-based
        assert_eq!(HgvsProteinPos(1).to_index(), ProteinPos(0));
        assert_eq!(HgvsProteinPos(11).to_index(), ProteinPos(10));
        assert_eq!(ProteinPos::from(HgvsProteinPos(6)), ProteinPos(5));
    }

    #[test]
    fn test_invalid_edge_cases_no_panic() {
        // Test position 0 in systems that shouldn't have it according to docs
        // Even if invalid in HGVS, we want to ensure the code doesn't panic.

        // Transcript systems
        assert_eq!(HgvsTranscriptPos(0).to_index(), TranscriptPos(0));
        assert_eq!(TranscriptPos(0).to_hgvs(), HgvsTranscriptPos(1));

        // Protein system
        assert_eq!(HgvsProteinPos(0).to_index(), ProteinPos(-1));
        assert_eq!(ProteinPos(-1).to_hgvs(), HgvsProteinPos(0));

        // Extremely large/small values
        assert_eq!(HgvsGenomicPos(i32::MAX).to_index(), GenomicPos(i32::MAX - 1));
        assert_eq!(HgvsTranscriptPos(i32::MIN).to_index(), TranscriptPos(i32::MIN));
    }

    #[test]
    fn test_arithmetic() {
        // GenomicPos
        assert_eq!(GenomicPos(100) + 10, GenomicPos(110));
        assert_eq!(GenomicPos(100) - 10, GenomicPos(90));
        assert_eq!(GenomicPos(100) - GenomicPos(90), 10);

        // TranscriptPos
        assert_eq!(TranscriptPos(100) + 10, TranscriptPos(110));
        assert_eq!(TranscriptPos(100) - 10, TranscriptPos(90));
        assert_eq!(TranscriptPos(100) - TranscriptPos(90), 10);

        // ProteinPos
        assert_eq!(ProteinPos(100) + 10, ProteinPos(110));
        assert_eq!(ProteinPos(100) - 10, ProteinPos(90));
        assert_eq!(ProteinPos(100) - ProteinPos(90), 10);

        // IntronicOffset
        assert_eq!(IntronicOffset(10) + 5, IntronicOffset(15));
        assert_eq!(IntronicOffset(10) - 5, IntronicOffset(5));
        assert_eq!(IntronicOffset(10) - IntronicOffset(5), 5);
    }
}
