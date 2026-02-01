use crate::error::HgvsError;
pub use crate::coords::{Anchor, GenomicPos, TranscriptPos, ProteinPos, IntronicOffset, HgvsGenomicPos, HgvsTranscriptPos, HgvsProteinPos, SequenceVariant};
pub use crate::edits::{NaEdit, AaEdit};
pub use crate::data::IdentifierKind;
use serde::{Serialize, Deserialize};

/// Common trait for all HGVS variants.
pub trait Variant {
    /// Returns the primary accession (e.g., "NM_000051.3").
    fn ac(&self) -> &str;
    /// Returns the optional gene symbol (e.g., "ATM").
    fn gene(&self) -> Option<&str>;
    /// Returns the coordinate type code ("g", "c", "p", etc.).
    fn coordinate_type(&self) -> &str;
}

macro_rules! impl_variant {
    ($struct_name:ident, $type_code:expr) => {
        impl Variant for $struct_name {
            fn ac(&self) -> &str { &self.ac }
            fn gene(&self) -> Option<&str> { self.gene.as_deref() }
            fn coordinate_type(&self) -> &str { $type_code }
        }
    };
}

/// Represents a genomic variant (g.).
#[derive(Debug, PartialEq, Clone, Serialize, Deserialize)]
pub struct GVariant {
    pub ac: String,
    pub gene: Option<String>,
    pub posedit: PosEdit<SimpleInterval, NaEdit>,
}
impl_variant!(GVariant, "g");

/// Represents a coding cDNA variant (c.).
#[derive(Debug, PartialEq, Clone, Serialize, Deserialize)]
pub struct CVariant {
    pub ac: String,
    pub gene: Option<String>,
    pub posedit: PosEdit<BaseOffsetInterval, NaEdit>,
}
impl_variant!(CVariant, "c");

/// Represents a protein variant (p.).
#[derive(Debug, PartialEq, Clone, Serialize, Deserialize)]
pub struct PVariant {
    pub ac: String,
    pub gene: Option<String>,
    pub posedit: PosEdit<AaInterval, AaEdit>,
}
impl_variant!(PVariant, "p");

/// Represents a mitochondrial variant (m.).
#[derive(Debug, PartialEq, Clone, Serialize, Deserialize)]
pub struct MVariant {
    pub ac: String,
    pub gene: Option<String>,
    pub posedit: PosEdit<SimpleInterval, NaEdit>,
}
impl_variant!(MVariant, "m");

/// Represents a non-coding transcript variant (n.).
#[derive(Debug, PartialEq, Clone, Serialize, Deserialize)]
pub struct NVariant {
    pub ac: String,
    pub gene: Option<String>,
    pub posedit: PosEdit<BaseOffsetInterval, NaEdit>,
}
impl_variant!(NVariant, "n");

/// Represents an RNA variant (r.).
#[derive(Debug, PartialEq, Clone, Serialize, Deserialize)]
pub struct RVariant {
    pub ac: String,
    pub gene: Option<String>,
    pub posedit: PosEdit<BaseOffsetInterval, NaEdit>,
}
impl_variant!(RVariant, "r");

/// Combines an interval and an edit (e.g., `123A>G`).
#[derive(Debug, PartialEq, Clone, Serialize, Deserialize)]
pub struct PosEdit<I, E> {
    /// The location of the variant.
    pub pos: Option<I>,
    /// The type of change (substitution, deletion, etc.).
    pub edit: E,
    /// Whether the variant is uncertain (indicated by `?`).
    pub uncertain: bool,
    /// Whether the variant is predicted (indicated by `()`).
    pub predicted: bool,
}

/// An interval spanning simple genomic or mitochondrial coordinates.
#[derive(Debug, PartialEq, Clone, Serialize, Deserialize)]
pub struct SimpleInterval {
    pub start: SimplePosition,
    pub end: Option<SimplePosition>,
    pub uncertain: bool,
}

impl SimpleInterval {
    pub fn length(&self) -> Result<i32, HgvsError> {
        match &self.end {
            Some(end) => Ok(end.base.0 - self.start.base.0 + 1),
            None => Ok(1),
        }
    }
}

/// A simple position in genomic or mitochondrial coordinates.
#[derive(Debug, PartialEq, Clone, Copy, Serialize, Deserialize)]
pub struct SimplePosition {
    pub base: HgvsGenomicPos,
    pub end: Option<HgvsGenomicPos>,
    pub uncertain: bool,
}

/// An interval spanning cDNA, n. or r. coordinates.
#[derive(Debug, PartialEq, Clone, Serialize, Deserialize)]
pub struct BaseOffsetInterval {
    pub start: BaseOffsetPosition,
    pub end: Option<BaseOffsetPosition>,
    pub uncertain: bool,
}

impl BaseOffsetInterval {
    pub fn length(&self) -> Result<i32, HgvsError> {
        match &self.end {
            Some(end) => {
                if self.start.anchor != end.anchor || self.start.offset.is_some() || end.offset.is_some() {
                     return Err(HgvsError::UnsupportedOperation("Complex interval length calculation not implemented".into()));
                }
                Ok(end.base.0 - self.start.base.0 + 1)
            }
            None => Ok(1),
        }
    }
}

/// A position in cDNA, n. or r. coordinates.
#[derive(Debug, PartialEq, Clone, Copy, Serialize, Deserialize)]
pub struct BaseOffsetPosition {
    pub base: HgvsTranscriptPos,
    pub offset: Option<IntronicOffset>,
    pub anchor: Anchor,
    pub uncertain: bool,
}

/// An interval spanning amino acid positions.
#[derive(Debug, PartialEq, Clone, Serialize, Deserialize)]
pub struct AaInterval {
    pub start: AAPosition,
    pub end: Option<AAPosition>,
    pub uncertain: bool,
}

impl AaInterval {
    pub fn length(&self) -> Result<i32, HgvsError> {
        match &self.end {
            Some(end) => Ok(end.base.0 - self.start.base.0 + 1),
            None => Ok(1),
        }
    }
}

/// A position in a protein sequence.
#[derive(Debug, PartialEq, Clone, Serialize, Deserialize)]
pub struct AAPosition {
    pub base: HgvsProteinPos,
    pub aa: String,
    pub uncertain: bool,
}
