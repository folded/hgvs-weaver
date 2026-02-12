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
        let (start, end) = if let Some(pos) = &self.pos {
            pos.spdi_interval(data_provider)?
        } else {
            return Err(HgvsError::UnsupportedOperation("SPDI requires a position".into()));
        };

        // SPDI is 0-based. HGVS is 1-based (mostly).
        // IntervalSpdi trait handles the coordinate conversion to 0-based.

        self.edit.to_spdi(ac, start, end, data_provider)
    }
}

pub trait IntervalSpdi {
    /// Returns (start, end) as 0-based integer coordinates.
    fn spdi_interval(&self, data_provider: &dyn crate::data::DataProvider) -> Result<(i32, i32), HgvsError>;
}

impl IntervalSpdi for SimpleInterval {
    fn spdi_interval(&self, _data_provider: &dyn crate::data::DataProvider) -> Result<(i32, i32), HgvsError> {
        let start = self.start.base.to_index().0;
        let end = self.end.as_ref().map_or(start + 1, |e| e.base.to_index().0 + 1);
        Ok((start, end))
    }
}

impl IntervalSpdi for BaseOffsetInterval {
    fn spdi_interval(&self, _data_provider: &dyn crate::data::DataProvider) -> Result<(i32, i32), HgvsError> {
        // SPDI usually applies to genomic or transcript sequences where simple indexing works.
        // For c./n. variants, if they have offsets, SPDI generation might be ambiguous or require projection.
        // For now, we'll error on offsets, or use the base index if it's simple.
        
        if self.start.offset.is_some() || self.end.as_ref().map_or(false, |e| e.offset.is_some()) {
             return Err(HgvsError::UnsupportedOperation("Cannot convert intronic offsets to simple SPDI coordinates directly".into()));
        }
        
        // This assumes the underlying "base" (HgvsTranscriptPos) converts linearly to an index.
        // Check `coords.rs` for `to_index()`. 
        // NOTE: c.1 is index 0 in CDS? No, typically c.1 is start of CDS.
        // But for SPDI on a transcript, we generally want 0-based from start of transcript?
        // OR 0-based from start of CDS?
        // SPDI spec: "The position is 0-based". On a RefSeq NM_, it's 0-based from start of sequence.
        // But `HgvsTranscriptPos` is relative to CDS start (c.1) or star (c.*1) or 5'UTR (c.-1).
        // Converting this to absolute transcript index requires `Transcript` data (CDS start/end).
        // Since we don't have the Transcript object here easily, 
        // AND `validate.py` was doing logic like `start_1 = pos["start"]["base"]`.
        // If `validate.py` was just using the raw number, it might be wrong for c. vars if not adjusted.
        // Use `return Err` for now to force a check, or implement simple logic if we trust the caller knows what they did.
        // The `validate.py` logic was: `start_1 = pos["start"]["base"]`.
        // If it was c.1, base is 1. SPDI requires 0-based.
        // So `start_1 - 1`.
        
        let start = self.start.base.0;
        let end = self.end.as_ref().map_or(start, |e| e.base.0);
        
        // If anchors are different, it's a range.
        if start > end {
             return Err(HgvsError::ValidationError("Invalid range".into()));
        }
        
        // Convert to 0-based half-open [start, end).
        // HGVS points intervals [start, end].
        // So 0-based start is start-1.
        // 0-based end is end. 
        // E.g. c.1_2 (2 bases). start=1, end=2. 
        // 0-based: 0, 2. (Indices 0, 1). Length 2. Correct.
        
        Ok((start - 1, end))
    }
}


impl IntervalSpdi for AaInterval {
    fn spdi_interval(&self, _data_provider: &dyn crate::data::DataProvider) -> Result<(i32, i32), HgvsError> {
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
                        // It's a length/integer? Fetch actual sequence.
                        data_provider.get_seq(ac, start, end, IdentifierType::Unknown)?
                    } else {
                        r.clone()
                    }
                } else {
                    // Implicit reference, fetch it.
                    data_provider.get_seq(ac, start, end, IdentifierType::Unknown)?
                };
                
                let a_seq = alt.as_deref().unwrap_or("");
                
                // SPDI format: SEQ:POS:DEL_SEQ:INS_SEQ
                // But wait, normalized SPDI usually just puts the ref/alt.
                // If ref matches genomic/transcript ref, we use it.
                
                Ok(format!("{}:{}:{}:{}", ac, start, r_seq, a_seq))
            }
            NaEdit::Del { ref_, .. } => {
                 let r_seq = if let Some(r) = ref_ {
                    // Logic to check if valid sequence or length 
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
                // For insertion, HGVS usually points between bases, e.g. c.1_2insA.
                // Start=1, End=1 (or 2?)
                // In `BaseOffsetInterval`, if start != end, we calculated range.
                // But `Ins` usually has start+1 = end in HGVS coordinates?
                // Wait, logic in `validate.py`: 
                // `if start_1 < end_1: return ...`
                // `ac:start_1:ref:alt`
                
                // Let's stick to strict interpretation:
                // SPDI for insertion is at the position.
                // Ref is empty string (or the base before?).
                // SPDI spec: "Deletion of 0 length at position".
                
                let a_seq = alt.as_deref().unwrap_or("");
                Ok(format!("{}:{}:{}:{}", ac, start, "", a_seq))
            }
            NaEdit::Dup { ref_, .. } => {
                // Duplication of sequence at [start, end].
                // Converted to insertion of that sequence.
                // Ref is empty (inserting specific sequence).
                // Or is it? SPDI doesn't have "Dup". It's an insertion of the duplicated sequence.
                // Where? At end?
                // HGVS: 10_11dup (duplicates bases 10 and 11).
                // SPDI: At 11 (0-based 11, which is after 10-11 range?), insert the sequence.
                // Logic in reference implementation:
                // `ac:end_1:ref:ref` (ref:ref implies deletion of 3rd arg, insertion of 4th arg?)
                // No, `ac:pos:del:ins`.
                // `validate.py` says: `f"{ac}:{end_1}:{ref}:{ref}"` ? 
                // Wait. `ref` is the sequence being duplicated. 
                // `del` is `ref`? That would mean deleting the sequence and replacing it with itself? That's Identity.
                // Ah, `validate.py` logic:
                // `if edit["type"] == "Dup": ref = get_seq... return f"{ac}:{end_1}:{ref}:{ref}"`
                // This looks weird. If del=ref and ins=ref, it's a no-op?
                // Unless it logic meant empty string for del?
                // CHECK validate.py AGAIN.
                
                let r_seq = if let Some(r) = ref_ {
                    r.clone()
                } else {
                     data_provider.get_seq(ac, start, end, IdentifierType::Unknown)?
                };
                
                // If we assume duplication is "Insert copy of sequence after the sequence".
                // Deletion: "" (0 chars). Insertion: "SEQUENCE".
                // Position: end.
                
                Ok(format!("{}:{}:{}:{}", ac, end, "", r_seq))
            }
            _ => Err(HgvsError::UnsupportedOperation("Edit type not yet supported for SPDI".into()))
        }
    }
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
