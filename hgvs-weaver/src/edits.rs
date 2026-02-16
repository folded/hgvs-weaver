use crate::coords::SequenceVariant;
use crate::error::HgvsError;
use serde::{Deserialize, Serialize};

/// Nucleic acid edits (substitutions, deletions, insertions, etc.).
#[derive(Debug, PartialEq, Clone, Serialize, Deserialize)]
#[serde(tag = "type")]
pub enum NaEdit {
    /// Substitution, deletion, or insertion represented by reference and alternate sequences.
    RefAlt {
        ref_: Option<String>,
        alt: Option<String>,
        uncertain: bool,
    },
    /// Deletion of a sequence.
    Del {
        ref_: Option<String>,
        uncertain: bool,
    },
    /// Insertion of a sequence.
    Ins {
        alt: Option<String>,
        uncertain: bool,
    },
    /// Duplication of a sequence.
    Dup {
        ref_: Option<String>,
        uncertain: bool,
    },
    /// Inversion of a sequence.
    Inv {
        ref_: Option<String>,
        uncertain: bool,
    },
    /// Conversion to another variant sequence.
    Con {
        con: Box<SequenceVariant>,
        uncertain: bool,
    },
    /// Repeat sequence (e.g., `[10]`).
    Repeat {
        ref_: Option<String>,
        min: i32,
        max: i32,
        uncertain: bool,
    },
    /// Copy number change.
    NACopy { copy: i32, uncertain: bool },
    /// No change (identity).
    None,
}

/// Amino acid edits (substitutions, frameshifts, extensions, etc.).
#[derive(Debug, PartialEq, Clone, Serialize, Deserialize)]
#[serde(tag = "type")]
pub enum AaEdit {
    /// Simple substitution.
    Subst {
        ref_: String,
        alt: String,
        uncertain: bool,
    },
    /// Deletion of amino acids.
    Del { ref_: String, uncertain: bool },
    /// Insertion of amino acids.
    Ins { alt: String, uncertain: bool },
    /// Deletion-insertion.
    DelIns {
        ref_: String,
        alt: String,
        uncertain: bool,
    },
    /// Substitution, deletion, or insertion represented by reference and alternate sequences.
    RefAlt {
        ref_: Option<String>,
        alt: Option<String>,
        uncertain: bool,
    },
    /// Frameshift.
    Fs {
        ref_: String,
        alt: String,
        term: Option<String>,
        length: Option<String>,
        uncertain: bool,
    },
    /// Extension of the protein (stop codon loss).
    Ext {
        ref_: String,
        alt: String,
        aaterm: Option<String>,
        length: Option<String>,
        uncertain: bool,
    },
    /// Repeat sequence.
    Repeat {
        ref_: Option<String>,
        min: i32,
        max: i32,
        uncertain: bool,
    },
    /// Duplication.
    Dup {
        ref_: Option<String>,
        uncertain: bool,
    },
    /// Silent variant (no change).
    Identity { uncertain: bool },
    /// Special cases (e.g., `p.0`, `p.?`).
    Special { value: String, uncertain: bool },
    /// Placeholder for no edit.
    None,
}

impl NaEdit {
    /// Calculates deletion and insertion lengths for the edit.
    pub fn del_ins_lengths(&self, ilen: i32) -> Result<(i32, i32), HgvsError> {
        match self {
            NaEdit::RefAlt { ref_, alt, .. } => {
                let del_len = if let Some(r) = ref_ {
                    if r.is_empty() {
                        0
                    } else if r.chars().all(|c| c.is_ascii_digit()) {
                        r.parse::<i32>().unwrap_or(ilen)
                    } else {
                        r.len() as i32
                    }
                } else {
                    0
                };
                let ins_len = alt.as_ref().map_or(0, |a| {
                    if a.chars().all(|c| c.is_ascii_digit()) {
                        a.parse::<i32>().unwrap_or(0)
                    } else {
                        a.len() as i32
                    }
                });
                Ok((del_len, ins_len))
            }
            NaEdit::Del { ref_, .. } => {
                let del_len = if let Some(r) = ref_ {
                    if r.chars().all(|c| c.is_ascii_digit()) {
                        r.parse::<i32>().unwrap_or(ilen)
                    } else {
                        r.len() as i32
                    }
                } else {
                    ilen
                };
                Ok((del_len, 0))
            }
            NaEdit::Ins { alt, .. } => {
                let ins_len = alt.as_ref().map_or(0, |a| {
                    if a.chars().all(|c| c.is_ascii_digit()) {
                        a.parse::<i32>().unwrap_or(0)
                    } else {
                        a.len() as i32
                    }
                });
                Ok((0, ins_len))
            }
            NaEdit::Dup { .. } => Ok((0, ilen)),
            NaEdit::Inv { .. } => Ok((ilen, ilen)),
            _ => Err(HgvsError::UnsupportedOperation(
                "Not implemented for this edit type".into(),
            )),
        }
    }

    /// Returns the reverse complement of the edit (used for minus-strand mapping).
    pub fn reverse_complement(&self) -> NaEdit {
        match self {
            NaEdit::RefAlt {
                ref_,
                alt,
                uncertain,
            } => {
                let r = ref_.as_ref().map(|s| {
                    if s.chars().all(|c| c.is_ascii_digit()) {
                        s.clone()
                    } else {
                        crate::utils::reverse_complement(s)
                    }
                });
                let a = alt.as_ref().map(|s| {
                    if s.chars().all(|c| c.is_ascii_digit()) {
                        s.clone()
                    } else {
                        crate::utils::reverse_complement(s)
                    }
                });
                NaEdit::RefAlt {
                    ref_: r,
                    alt: a,
                    uncertain: *uncertain,
                }
            }
            NaEdit::Del { ref_, uncertain } => {
                let r = ref_.as_ref().map(|s| {
                    if s.chars().all(|c| c.is_ascii_digit()) {
                        s.clone()
                    } else {
                        crate::utils::reverse_complement(s)
                    }
                });
                NaEdit::Del {
                    ref_: r,
                    uncertain: *uncertain,
                }
            }
            NaEdit::Ins { alt, uncertain } => {
                let a = alt.as_ref().map(|s| {
                    if s.chars().all(|c| c.is_ascii_digit()) {
                        s.clone()
                    } else {
                        crate::utils::reverse_complement(s)
                    }
                });
                NaEdit::Ins {
                    alt: a,
                    uncertain: *uncertain,
                }
            }
            NaEdit::Dup { ref_, uncertain } => {
                let r = ref_.as_ref().map(|s| crate::utils::reverse_complement(s));
                NaEdit::Dup {
                    ref_: r,
                    uncertain: *uncertain,
                }
            }
            NaEdit::Inv { ref_, uncertain } => {
                let r = ref_.as_ref().map(|s| crate::utils::reverse_complement(s));
                NaEdit::Inv {
                    ref_: r,
                    uncertain: *uncertain,
                }
            }
            NaEdit::Repeat {
                ref_,
                min,
                max,
                uncertain,
            } => {
                let r = ref_.as_ref().map(|s| crate::utils::reverse_complement(s));
                NaEdit::Repeat {
                    ref_: r,
                    min: *min,
                    max: *max,
                    uncertain: *uncertain,
                }
            }
            _ => self.clone(),
        }
    }

    /// Applies a function to all sequence strings within the edit.
    pub fn map_sequence<F>(self, f: F) -> NaEdit
    where
        F: Fn(&str) -> String + Copy,
    {
        match self {
            NaEdit::RefAlt {
                ref_,
                alt,
                uncertain,
            } => NaEdit::RefAlt {
                ref_: ref_.map(|s| f(&s)),
                alt: alt.map(|s| f(&s)),
                uncertain,
            },
            NaEdit::Del { ref_, uncertain } => NaEdit::Del {
                ref_: ref_.map(|s| f(&s)),
                uncertain,
            },
            NaEdit::Ins { alt, uncertain } => NaEdit::Ins {
                alt: alt.map(|s| f(&s)),
                uncertain,
            },
            NaEdit::Dup { ref_, uncertain } => NaEdit::Dup {
                ref_: ref_.map(|s| f(&s)),
                uncertain,
            },
            NaEdit::Inv { ref_, uncertain } => NaEdit::Inv {
                ref_: ref_.map(|s| f(&s)),
                uncertain,
            },
            NaEdit::Repeat {
                ref_,
                min,
                max,
                uncertain,
            } => NaEdit::Repeat {
                ref_: ref_.map(|s| f(&s)),
                min,
                max,
                uncertain,
            },
            _ => self,
        }
    }
}
