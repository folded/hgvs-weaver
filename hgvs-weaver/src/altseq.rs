use crate::error::HgvsError;
use crate::sequence::{MemSequence, Sequence, SliceSequence, SplicedSequence, TranslatedSequence};
use crate::structs::{Anchor, CVariant, NaEdit, ProteinPos, TranscriptPos};

/// Represents the data for a transcript with a variant applied.
pub struct AltTranscriptData {
    pub transcript_sequence: String,
    pub aa_sequence: String,
    pub cds_start_index: TranscriptPos,
    pub cds_end_index: TranscriptPos,
    pub protein_accession: String,
    pub is_frameshift: bool,
    /// The index of the first affected amino acid.
    pub variant_start_aa: Option<ProteinPos>,
    pub frameshift_start: Option<ProteinPos>,
    pub is_substitution: bool,
    pub is_ambiguous: bool,
    /// The original cDNA variant.
    pub c_variant: CVariant,
}

pub struct AltSeqBuilder<'a> {
    pub var_c: &'a CVariant,
    pub transcript_sequence: &'a dyn Sequence,
    pub cds_start_index: TranscriptPos,
    pub cds_end_index: TranscriptPos,
    pub protein_accession: String,
}

impl<'a> AltSeqBuilder<'a> {
    pub fn build_altseq(&self) -> Result<AltTranscriptData, HgvsError> {
        let (start_idx, end_idx) = self.get_variant_indices()?;

        // --- Validate reference sequence ---
        match &self.var_c.posedit.edit {
            NaEdit::RefAlt { ref_: Some(r), .. }
            | NaEdit::Del { ref_: Some(r), .. }
            | NaEdit::Dup { ref_: Some(r), .. } => {
                if !r.is_empty() && !r.chars().all(|c| c.is_ascii_digit()) {
                    let actual_ref = SliceSequence {
                        inner: self.transcript_sequence,
                        start: start_idx,
                        end: end_idx,
                    }
                    .to_string();
                    if actual_ref != *r {
                        return Err(HgvsError::TranscriptMismatch {
                            expected: r.to_string(),
                            found: actual_ref,
                            start: start_idx,
                            end: end_idx,
                        });
                    }
                }
            }
            _ => {}
        }
        // --- End validation ---

        let pos = self
            .var_c
            .posedit
            .pos
            .as_ref()
            .ok_or_else(|| HgvsError::ValidationError("Missing position".into()))?;
        let pos_start_c_0 = pos.start.base.to_index();
        let variant_start_aa = Some(ProteinPos(pos_start_c_0.0.max(0) / 3));

        let (is_substitution, is_frameshift, alt_transcript) = match &self.var_c.posedit.edit {
            NaEdit::RefAlt { ref_, alt, .. } => {
                let is_identity = ref_.is_none() && alt.is_none();
                let is_ins = ref_.is_none()
                    && !is_identity
                    && self
                        .var_c
                        .posedit
                        .pos
                        .as_ref()
                        .is_some_and(|p| p.end.is_some());
                let alt_str = alt.as_deref().unwrap_or("");

                let (is_subst, is_fs, res) = if is_identity {
                    (false, false, self.transcript_sequence.to_string())
                } else if is_ins {
                    let alt_seq = MemSequence(alt_str.to_string());
                    let ins_pos = (start_idx + 1).min(self.transcript_sequence.len());
                    let p1 = SliceSequence {
                        inner: self.transcript_sequence,
                        start: 0,
                        end: ins_pos,
                    };
                    let p3 = SliceSequence {
                        inner: self.transcript_sequence,
                        start: ins_pos,
                        end: self.transcript_sequence.len(),
                    };
                    (
                        false,
                        alt_str.len() % 3 != 0,
                        SplicedSequence {
                            pieces: vec![
                                &p1 as &dyn Sequence,
                                &alt_seq as &dyn Sequence,
                                &p3 as &dyn Sequence,
                            ],
                        }
                        .to_string(),
                    )
                } else {
                    let alt_seq = MemSequence(alt_str.to_string());
                    let r_len = if let Some(r) = ref_ {
                        if r.chars().all(|c| c.is_ascii_digit()) {
                            end_idx - start_idx
                        } else {
                            r.len()
                        }
                    } else {
                        end_idx - start_idx
                    };

                    let p1 = SliceSequence {
                        inner: self.transcript_sequence,
                        start: 0,
                        end: start_idx,
                    };
                    let p3 = SliceSequence {
                        inner: self.transcript_sequence,
                        start: end_idx,
                        end: self.transcript_sequence.len(),
                    };
                    let is_subst =
                        ref_.is_some() && alt.is_some() && r_len == 1 && alt_str.len() == 1;
                    let is_fs = (alt_str.len() as i32 - r_len as i32) % 3 != 0;
                    (
                        is_subst,
                        is_fs,
                        SplicedSequence {
                            pieces: vec![
                                &p1 as &dyn Sequence,
                                &alt_seq as &dyn Sequence,
                                &p3 as &dyn Sequence,
                            ],
                        }
                        .to_string(),
                    )
                };

                (is_subst, is_fs, res)
            }
            NaEdit::Del { ref_, .. } => {
                let r_len = if let Some(r) = ref_ {
                    if r.chars().all(|c| c.is_ascii_digit()) {
                        end_idx - start_idx
                    } else {
                        r.len()
                    }
                } else {
                    end_idx - start_idx
                };

                let p1 = SliceSequence {
                    inner: self.transcript_sequence,
                    start: 0,
                    end: start_idx,
                };
                let p3 = SliceSequence {
                    inner: self.transcript_sequence,
                    start: end_idx,
                    end: self.transcript_sequence.len(),
                };
                let res = SplicedSequence {
                    pieces: vec![&p1 as &dyn Sequence, &p3 as &dyn Sequence],
                }
                .to_string();
                (false, (r_len as i32) % 3 != 0, res)
            }
            NaEdit::Ins { alt: Some(alt), .. } => {
                let alt_seq = MemSequence(alt.clone());
                let ins_pos = (start_idx + 1).min(self.transcript_sequence.len());
                let p1 = SliceSequence {
                    inner: self.transcript_sequence,
                    start: 0,
                    end: ins_pos,
                };
                let p3 = SliceSequence {
                    inner: self.transcript_sequence,
                    start: ins_pos,
                    end: self.transcript_sequence.len(),
                };
                let res = SplicedSequence {
                    pieces: vec![
                        &p1 as &dyn Sequence,
                        &alt_seq as &dyn Sequence,
                        &p3 as &dyn Sequence,
                    ],
                }
                .to_string();
                (false, (alt.len() as i32) % 3 != 0, res)
            }
            NaEdit::Dup { ref_, .. } => {
                let dup_str = if let Some(r) = ref_ {
                    if r.chars().all(|c| c.is_ascii_digit()) {
                        SliceSequence {
                            inner: self.transcript_sequence,
                            start: start_idx,
                            end: end_idx,
                        }
                        .to_string()
                    } else {
                        r.clone()
                    }
                } else {
                    SliceSequence {
                        inner: self.transcript_sequence,
                        start: start_idx,
                        end: end_idx,
                    }
                    .to_string()
                };

                let dup_seq = MemSequence(dup_str.clone());
                let ins_pos = end_idx.min(self.transcript_sequence.len());
                let p1 = SliceSequence {
                    inner: self.transcript_sequence,
                    start: 0,
                    end: ins_pos,
                };
                let p3 = SliceSequence {
                    inner: self.transcript_sequence,
                    start: ins_pos,
                    end: self.transcript_sequence.len(),
                };
                let res = SplicedSequence {
                    pieces: vec![
                        &p1 as &dyn Sequence,
                        &dup_seq as &dyn Sequence,
                        &p3 as &dyn Sequence,
                    ],
                }
                .to_string();
                (false, (dup_str.len() as i32) % 3 != 0, res)
            }
            NaEdit::Inv { .. } => {
                let sub = SliceSequence {
                    inner: self.transcript_sequence,
                    start: start_idx,
                    end: end_idx,
                }
                .to_string();
                let inv_str = crate::utils::reverse_complement(&sub);
                let inv_seq = MemSequence(inv_str);
                let p1 = SliceSequence {
                    inner: self.transcript_sequence,
                    start: 0,
                    end: start_idx,
                };
                let p3 = SliceSequence {
                    inner: self.transcript_sequence,
                    start: end_idx,
                    end: self.transcript_sequence.len(),
                };
                let res = SplicedSequence {
                    pieces: vec![
                        &p1 as &dyn Sequence,
                        &inv_seq as &dyn Sequence,
                        &p3 as &dyn Sequence,
                    ],
                }
                .to_string();
                (false, false, res)
            }
            NaEdit::Repeat { min, ref_, .. } => {
                let unit = if let Some(r) = ref_ {
                    r.clone()
                } else {
                    SliceSequence {
                        inner: self.transcript_sequence,
                        start: start_idx,
                        end: end_idx,
                    }
                    .to_string()
                };

                let mut current_idx = start_idx;
                loop {
                    if current_idx + unit.len() > self.transcript_sequence.len() {
                        break;
                    }
                    let next_unit = SliceSequence {
                        inner: self.transcript_sequence,
                        start: current_idx,
                        end: current_idx + unit.len(),
                    }
                    .to_string();
                    if next_unit == unit {
                        current_idx += unit.len();
                    } else {
                        break;
                    }
                }

                let mut total_str = String::new();
                for _ in 0..*min {
                    total_str.push_str(&unit);
                }
                let alt_seq = MemSequence(total_str);
                let p1 = SliceSequence {
                    inner: self.transcript_sequence,
                    start: 0,
                    end: start_idx,
                };
                let p3 = SliceSequence {
                    inner: self.transcript_sequence,
                    start: current_idx,
                    end: self.transcript_sequence.len(),
                };
                let res = SplicedSequence {
                    pieces: vec![
                        &p1 as &dyn Sequence,
                        &alt_seq as &dyn Sequence,
                        &p3 as &dyn Sequence,
                    ],
                }
                .to_string();

                let net_change =
                    (unit.len() as i32 * (*min as i32)) - (current_idx as i32 - start_idx as i32);
                (false, net_change % 3 != 0, res)
            }
            NaEdit::None => (false, false, self.transcript_sequence.to_string()),
            _ => {
                return Err(HgvsError::UnsupportedOperation(
                    "Unsupported edit for altseq".into(),
                ))
            }
        };

        let cds_start = self.cds_start_index.0 as usize;
        let alt_transcript_seq = MemSequence(alt_transcript);
        let aa_sequence = if cds_start < alt_transcript_seq.len() {
            let slice = SliceSequence {
                inner: &alt_transcript_seq,
                start: cds_start,
                end: alt_transcript_seq.len(),
            };
            TranslatedSequence { inner: &slice }.to_string()
        } else {
            "".to_string()
        };

        // Find cds_end_i in the new sequence.
        // It's original_cds_end + net_change.
        // Or simply finding the length of the alt_transcript_seq.
        // Actually, we should be careful with 3' UTR.
        let net_change = alt_transcript_seq.len() as i32 - self.transcript_sequence.len() as i32;
        let cds_end_i = self.cds_end_index.0 + net_change;

        let mut is_fs = is_frameshift;
        if is_fs {
            if let Some(v_start_aa) = variant_start_aa {
                if let Some(c) = aa_sequence.chars().nth(v_start_aa.0 as usize) {
                    if c == '*' {
                        is_fs = false;
                    }
                }
            }
        }

        Ok(AltTranscriptData {
            transcript_sequence: alt_transcript_seq.0,
            aa_sequence,
            cds_start_index: self.cds_start_index,
            cds_end_index: TranscriptPos(cds_end_i),
            protein_accession: self.protein_accession.clone(),
            is_frameshift: is_fs,
            variant_start_aa,
            frameshift_start: if is_fs { variant_start_aa } else { None },
            is_substitution,
            is_ambiguous: false,
            c_variant: self.var_c.clone(),
        })
    }

    fn get_variant_indices(&self) -> Result<(usize, usize), HgvsError> {
        let pos = self
            .var_c
            .posedit
            .pos
            .as_ref()
            .ok_or_else(|| HgvsError::ValidationError("Missing position".into()))?;
        let start = self.pos_to_idx(&pos.start)?;
        let mut end = if let Some(e) = &pos.end {
            self.pos_to_idx(e)?
        } else {
            start
        };
        end += 1;
        Ok((start, end))
    }

    fn pos_to_idx(&self, pos: &crate::structs::BaseOffsetPosition) -> Result<usize, HgvsError> {
        let base_idx_0 = pos.base.to_index();

        if pos.offset.is_some() && pos.offset.unwrap().0 != 0 {
            return Err(HgvsError::UnsupportedOperation(
                "Intronic variants not yet supported in c_to_p".into(),
            ));
        }

        let idx = match pos.anchor {
            Anchor::TranscriptStart => {
                let i = base_idx_0.0;
                if i < 0 {
                    return Err(HgvsError::ValidationError(format!(
                        "Position {} before transcript start",
                        i
                    )));
                }
                i as usize
            }
            Anchor::CdsStart => {
                let i = (self.cds_start_index.0 + base_idx_0.0) as i32;
                if i < 0 {
                    return Err(HgvsError::ValidationError(format!(
                        "Position {} before transcript start",
                        i
                    )));
                }
                i as usize
            }
            Anchor::CdsEnd => {
                let i = (self.cds_end_index.0 + base_idx_0.0) as i32;
                if i < 0 {
                    return Err(HgvsError::ValidationError(format!(
                        "Position {} before transcript start",
                        i
                    )));
                }
                i as usize
            }
        };
        Ok(idx)
    }
}
