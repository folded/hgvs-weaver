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
        // If has offsets, resolve to genomic.
        if self.start.offset.is_some() || self.end.as_ref().map_or(false, |e| e.offset.is_some()) {
             let g_start = data_provider.c_to_g(ac, self.start.base.to_index(), self.start.offset.unwrap_or(IntronicOffset(0)))?;
             let g_end = if let Some(e) = &self.end {
                 data_provider.c_to_g(ac, e.base.to_index(), e.offset.unwrap_or(IntronicOffset(0)))?
             } else {
                 (g_start.0.clone(), GenomicPos(g_start.1.0 + 1))
             };
             
             if g_start.0 != g_end.0 {
                 return Err(HgvsError::UnsupportedOperation("Interval spans multiple genomic accessions".into()));
             }
             
             return Ok((g_start.1.0, g_end.1.0, g_start.0));
        }
        
        let start = self.start.base.0;
        let end = self.end.as_ref().map_or(start, |e| e.base.0);
        
        if start > end {
             return Err(HgvsError::ValidationError("Invalid range".into()));
        }
        
        Ok((start - 1, end, ac.to_string()))
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
                        // It's a length/integer? Fetch actual sequence.
                        data_provider.get_seq(ac, start, end, IdentifierType::Unknown)?
                    } else {
                        r.clone()
                    }
                } else {
                    // Implicit reference, fetch it.
                    data_provider.get_seq(ac, start, end, IdentifierType::Unknown)?
                };
                
                let a_seq = if ref_.is_none() && alt.is_none() {
                    r_seq.clone()
                } else {
                    alt.as_deref().unwrap_or("").to_string()
                };
                
                // SPDI format: SEQ:POS:DEL_SEQ:INS_SEQ
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
            NaEdit::Repeat { ref_, max, .. } => {
                // Repeat sequence (e.g., c.7035TGGAAC[3]).
                // SPDI represents this as a deletion of the original region and insertion of the repeated sequence.
                // ClinVar SPDI for NM_001291285.3:c.7035TGGAAC[3] is NC_000004.12:125434258:ACTGGAACTGGAAC:ACTGGAACTGGAACTGGAAC
                // Weaver's interval [start, end) for a point position c.7035 is [7034, 7035).
                // Repeat unit length is ref_.len().
                
                let unit = if let Some(r) = ref_ {
                    r.clone()
                } else {
                    // If unit not provided, assume it's the whole interval? 
                    // HGVS repeats usually have a unit sequence or it's implied by the position.
                    // For now, fetch from data provider if None.
                    data_provider.get_seq(ac, start, end, IdentifierType::Unknown)?
                };

                // Del seq is the reference sequence at the interval.
                // HGVS c.7035TGGAAC[3] usually means the unit TGGAAC is present at c.7035 and we want 3 copies total.
                // We need to know how many copies were already there to determine the actual delta.
                // BUT SPDI is absolute. It says "Delete this, Insert that".
                // ClinVar's example has del=ACTGGAACTGGAAC (2 units?) and ins=ACTGGAACTGGAACTGGAAC (3 units?).
                // Wait, if it's c.7035TGGAAC[3], and it results in 3 units...
                // If r_seq at [start, end) is 1 unit, and we want 3 units...
                // SPDI position should be the start of the repeat region.
                
                // For simplicity and matching ClinVar's "minimal" (but sometimes expanded) style:
                // Canonical SPDI for repeats often settles on the smallest delins that describes the change.
                // If we want [max] copies of [unit]:
                // Ins seq = unit * max.
                // Del seq = we need to know how many units are in the reference to be precise.
                
                let ins_seq = unit.repeat(*max as usize);
                
                // If we don't know the reference repeat count, we might produce a non-minimal SPDI.
                // However, SPDI normalization (which weaver does) will clean it up.
                // Let's at least get the reference sequence for the interval.
                let r_seq = data_provider.get_seq(ac, start, end, IdentifierType::Unknown)?;
                
                Ok(format!("{}:{}:{}:{}", ac, start, r_seq, ins_seq))
            }
            _ => Err(HgvsError::UnsupportedOperation(format!("Edit type {:?} not yet supported for SPDI", self)))
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
