use crate::error::HgvsError;
pub use crate::coords::{Anchor, GenomicPos, TranscriptPos, ProteinPos, IntronicOffset, HgvsGenomicPos, HgvsTranscriptPos, HgvsProteinPos, SequenceVariant};
pub use crate::edits::{NaEdit, AaEdit};
pub use crate::data::{IdentifierKind, IdentifierType};
use serde::{Serialize, Deserialize};

/// Common trait for all HGVS variants.
pub trait Variant {
    /// Returns the primary accession (e.g., "NM_000051.3").
    fn ac(&self) -> &str;
    /// Returns the optional gene symbol (e.g., "ATM").
    fn gene(&self) -> Option<&str>;
    /// Returns the coordinate type code ("g", "c", "p", etc.).
    fn coordinate_type(&self) -> &str;
    /// Converts the variant to an SPDI string representation.
    fn to_spdi(&self, data_provider: &dyn crate::data::DataProvider) -> Result<String, HgvsError>;
}

macro_rules! impl_variant {
    ($struct_name:ident, $type_code:expr) => {
        impl Variant for $struct_name {
            fn ac(&self) -> &str { &self.ac }
            fn gene(&self) -> Option<&str> { self.gene.as_deref() }
            fn coordinate_type(&self) -> &str { $type_code }
            fn to_spdi(&self, data_provider: &dyn crate::data::DataProvider) -> Result<String, HgvsError> {
                self.posedit.to_spdi(&self.ac, data_provider)
            }
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

impl<I, E> PosEdit<I, E>
where
    I: IntervalSpdi,
    E: EditSpdi,
{
    pub fn to_spdi(&self, ac: &str, data_provider: &dyn crate::data::DataProvider) -> Result<String, HgvsError> {
        let (start, end, spdi_ac) = if let Some(pos) = &self.pos {
            pos.spdi_interval(ac, data_provider)?
        } else {
            return Err(HgvsError::ValidationError("SPDI requires a position".into()));
        };

        // SPDI is 0-based. HGVS is 1-based (mostly).
        // IntervalSpdi trait handles the coordinate conversion to 0-based.

        self.edit.to_spdi(&spdi_ac, start, end, data_provider)
    }
}

pub trait IntervalSpdi {
    /// Returns (start, end, ac) as 0-based integer coordinates and the accession to use for SPDI.
    fn spdi_interval(&self, ac: &str, data_provider: &dyn crate::data::DataProvider) -> Result<(i32, i32, String), HgvsError>;
}

impl IntervalSpdi for SimpleInterval {
    fn spdi_interval(&self, ac: &str, _data_provider: &dyn crate::data::DataProvider) -> Result<(i32, i32, String), HgvsError> {
        let start = self.start.base.to_index().0;
        let end = self.end.as_ref().map_or(start + 1, |e| e.base.to_index().0 + 1);
        Ok((start, end, ac.to_string()))
    }
}

impl IntervalSpdi for BaseOffsetInterval {
    fn spdi_interval(&self, ac: &str, data_provider: &dyn crate::data::DataProvider) -> Result<(i32, i32, String), HgvsError> {
        // Resolving transcript/non-coding positions to genomic coordinates for SPDI.
        // This ensures the resulting SPDI uses a chromosomal accession.
        let g_start = data_provider.c_to_g(ac, self.start.base.to_index(), self.start.offset.unwrap_or(IntronicOffset(0)))?;
        let g_end = if let Some(e) = &self.end {
            data_provider.c_to_g(ac, e.base.to_index(), e.offset.unwrap_or(IntronicOffset(0)))?
        } else {
            (g_start.0.clone(), GenomicPos(g_start.1.0 + 1))
        };

        if g_start.0 != g_end.0 {
            return Err(HgvsError::UnsupportedOperation("Interval spans multiple genomic accessions".into()));
        }

        Ok((g_start.1.0, g_end.1.0, g_start.0))
    }
}


impl IntervalSpdi for AaInterval {
    fn spdi_interval(&self, _ac: &str, _data_provider: &dyn crate::data::DataProvider) -> Result<(i32, i32, String), HgvsError> {
         Err(HgvsError::UnsupportedOperation("SPDI not supported for protein variants (yet)".into()))
    }
}

pub trait EditSpdi {
    fn to_spdi(&self, ac: &str, start: i32, end: i32, data_provider: &dyn crate::data::DataProvider) -> Result<String, HgvsError>;
}

impl EditSpdi for NaEdit {
    fn to_spdi(&self, ac: &str, start: i32, end: i32, data_provider: &dyn crate::data::DataProvider) -> Result<String, HgvsError> {
        match self {
            NaEdit::RefAlt { ref_, alt, .. } => {
                let r_seq = if let Some(r) = ref_ {
                    if r.chars().all(|c| c.is_ascii_digit()) {
                        data_provider.get_seq(ac, start, end, IdentifierType::Unknown)?
                    } else {
                        r.clone()
                    }
                } else {
                    data_provider.get_seq(ac, start, end, IdentifierType::Unknown)?
                };
                
                let a_seq = if ref_.is_none() && alt.is_none() {
                    r_seq.clone()
                } else {
                    alt.as_deref().unwrap_or("").to_string()
                };
                
                let (p_start, r_strip, a_strip) = strip_common_prefix_suffix(start, &r_seq, &a_seq);
                Ok(format!("{}:{}:{}:{}", ac, p_start, r_strip, a_strip))
            }
            NaEdit::Del { ref_, .. } => {
                let r_seq = if let Some(r) = ref_ {
                    if r.chars().all(|c| c.is_ascii_digit()) {
                        data_provider.get_seq(ac, start, end, IdentifierType::Unknown)?
                    } else {
                        r.clone()
                    }
                } else {
                    data_provider.get_seq(ac, start, end, IdentifierType::Unknown)?
                };
                Ok(format!("{}:{}:{}:", ac, start, r_seq))
            }
            NaEdit::Ins { alt, .. } => {
                let a_seq = alt.as_deref().unwrap_or("");
                Ok(format!("{}:{}:{}:{}", ac, start, "", a_seq))
            }
            NaEdit::Dup { ref_, .. } => {
                let r_seq = if let Some(r) = ref_ {
                    r.clone()
                } else {
                     data_provider.get_seq(ac, start, end, IdentifierType::Unknown)?
                };
                Ok(format!("{}:{}:{}:{}", ac, end, "", r_seq))
            }
            NaEdit::Repeat { ref_, max, .. } => {
                let unit = if let Some(r) = ref_ {
                    r.clone()
                } else {
                    data_provider.get_seq(ac, start, end, IdentifierType::Unknown)?
                };
                let ins_seq = unit.repeat(*max as usize);
                let r_seq = data_provider.get_seq(ac, start, end, IdentifierType::Unknown)?;
                
                let (p_start, r_strip, a_strip) = strip_common_prefix_suffix(start, &r_seq, &ins_seq);
                Ok(format!("{}:{}:{}:{}", ac, p_start, r_strip, a_strip))
            }
            NaEdit::Inv { .. } => {
                let r_seq = data_provider.get_seq(ac, start, end, IdentifierType::Unknown)?;
                let a_seq = crate::sequence::rev_comp(&r_seq);
                
                let (p_start, r_strip, a_strip) = strip_common_prefix_suffix(start, &r_seq, &a_seq);
                Ok(format!("{}:{}:{}:{}", ac, p_start, r_strip, a_strip))
            }
            _ => Err(HgvsError::UnsupportedOperation(format!("Edit type {:?} not yet supported for SPDI", self)))
        }
    }
}

pub fn strip_common_prefix_suffix(start: i32, ref_seq: &str, alt_seq: &str) -> (i32, String, String) {
    let mut r_bytes = ref_seq.as_bytes();
    let mut a_bytes = alt_seq.as_bytes();
    let mut p_start = start;

    // Strip prefix
    let mut prefix_len = 0;
    while prefix_len < r_bytes.len() && prefix_len < a_bytes.len() && r_bytes[prefix_len] == a_bytes[prefix_len] {
        prefix_len += 1;
    }
    r_bytes = &r_bytes[prefix_len..];
    a_bytes = &a_bytes[prefix_len..];
    p_start += prefix_len as i32;

    // Strip suffix
    let mut suffix_len = 0;
    while suffix_len < r_bytes.len() && suffix_len < a_bytes.len() {
        let r_idx = r_bytes.len() - 1 - suffix_len;
        let a_idx = a_bytes.len() - 1 - suffix_len;
        if r_bytes[r_idx] == a_bytes[a_idx] {
            suffix_len += 1;
        } else {
            break;
        }
    }
    r_bytes = &r_bytes[..r_bytes.len() - suffix_len];
    a_bytes = &a_bytes[..a_bytes.len() - suffix_len];

    (
        p_start,
        String::from_utf8_lossy(r_bytes).to_string(),
        String::from_utf8_lossy(a_bytes).to_string(),
    )
}

impl EditSpdi for AaEdit {
    fn to_spdi(&self, _ac: &str, _start: i32, _end: i32, _data_provider: &dyn crate::data::DataProvider) -> Result<String, HgvsError> {
          Err(HgvsError::UnsupportedOperation("SPDI not supported for protein variants (yet)".into()))
    }
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
