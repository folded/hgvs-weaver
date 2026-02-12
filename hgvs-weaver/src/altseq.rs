use crate::error::HgvsError;
use crate::structs::{CVariant, NaEdit, Anchor, TranscriptPos, ProteinPos};

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
    pub transcript_sequence: String,
    pub cds_start_index: TranscriptPos,
    pub cds_end_index: TranscriptPos,
    pub protein_accession: String,
}

impl<'a> AltSeqBuilder<'a> {
    pub fn build_altseq(&self) -> Result<AltTranscriptData, HgvsError> {
        let mut seq: Vec<char> = self.transcript_sequence.chars().collect();
        let mut cds_end_i = self.cds_end_index.0;

        let (start_idx, end_idx) = self.get_variant_indices()?;

        let check_bounds = |s: usize, e: usize, len: usize| -> Result<(), HgvsError> {
            if s > len || e > len || s > e {
                return Err(HgvsError::ValidationError(format!("Splicing bounds error: start={}, end={}, len={}", s, e, len)));
            }
            Ok(())
        };

        // --- Validate reference sequence ---
        match &self.var_c.posedit.edit {
            NaEdit::RefAlt { ref_: Some(r), .. } | NaEdit::Del { ref_: Some(r), .. } | NaEdit::Dup { ref_: Some(r), .. } => {
                if !r.is_empty() && !r.chars().all(|c| c.is_ascii_digit()) {
                    check_bounds(start_idx, end_idx, self.transcript_sequence.len())?;
                    let actual_ref = &self.transcript_sequence[start_idx..end_idx];
                    if actual_ref != r {
                        return Err(HgvsError::ValidationError(format!(
                            "Reference sequence mismatch: expected {}, found {} at transcript indices {}..{}",
                            r, actual_ref, start_idx, end_idx
                        )));
                    }
                }
            }
            _ => {}
        }
        // --- End validation ---

        let pos = self.var_c.posedit.pos.as_ref().ok_or_else(|| HgvsError::ValidationError("Missing position".into()))?;
        let pos_start_c_0 = pos.start.base.to_index();
        let variant_start_aa = Some(ProteinPos(pos_start_c_0.0.max(0) / 3));

        let (is_substitution, is_frameshift) = match &self.var_c.posedit.edit {
            NaEdit::RefAlt { ref_, alt, .. } => {
                let is_ins = ref_.is_none() && self.var_c.posedit.pos.as_ref().is_some_and(|p| p.end.is_some());

                let r_len = if is_ins {
                    0
                } else if let Some(r) = ref_ {
                    if r.is_empty() { (end_idx - start_idx) as i32 }
                    else if r.chars().all(|c| c.is_ascii_digit()) { r.parse::<i32>().unwrap_or((end_idx - start_idx) as i32) }
                    else { r.len() as i32 }
                } else {
                    (end_idx - start_idx) as i32
                };

                let a_str = alt.as_deref().unwrap_or("");
                let a_len = a_str.len() as i32;

                if ref_.is_none() && alt.is_none() {
                    return Ok(AltTranscriptData {
                        transcript_sequence: self.transcript_sequence.clone(),
                        aa_sequence: crate::utils::translate_cds(&self.transcript_sequence[self.cds_start_index.0 as usize..]),
                        cds_start_index: self.cds_start_index,
                        cds_end_index: self.cds_end_index,
                        protein_accession: self.protein_accession.clone(),
                        is_frameshift: false,
                        variant_start_aa,
                        frameshift_start: None,
                        is_substitution: false,
                        is_ambiguous: false,
                        c_variant: self.var_c.clone(),
                    });
                }

                if let Some(r) = ref_ {
                    if r == a_str {
                        return Ok(AltTranscriptData {
                            transcript_sequence: self.transcript_sequence.clone(),
                            aa_sequence: crate::utils::translate_cds(&self.transcript_sequence[self.cds_start_index.0 as usize..]),
                            cds_start_index: self.cds_start_index,
                            cds_end_index: self.cds_end_index,
                            protein_accession: self.protein_accession.clone(),
                            is_frameshift: false,
                            variant_start_aa,
                            frameshift_start: None,
                            is_substitution: false,
                            is_ambiguous: false,
                            c_variant: self.var_c.clone(),
                        });
                    }
                }

                let net_change = a_len - r_len;
                let is_fs = net_change % 3 != 0;
                cds_end_i += net_change;

                let mut is_subst = false;
                if let (Some(_r), Some(a)) = (ref_, alt) {
                    if r_len == 1 && a.len() == 1 && start_idx + 1 == end_idx { is_subst = true; }
                    let a_chars: Vec<char> = a.chars().collect();
                    check_bounds(start_idx, end_idx, seq.len())?;
                    // println!("DEBUG: splice #1 at {}..{} with {:?} (len {})", start_idx, end_idx, a_chars, a_chars.len());
                    seq.splice(start_idx..end_idx, a_chars);
                } else if let Some(a) = alt {
                    let a_chars: Vec<char> = a.chars().collect();
                    if is_ins {
                        let ins_pos = (start_idx + 1).min(seq.len());
                        check_bounds(ins_pos, ins_pos, seq.len())?;
                        seq.splice(ins_pos..ins_pos, a_chars);
                    } else {
                        check_bounds(start_idx, end_idx, seq.len())?;
                        // println!("DEBUG: splice #2 at {}..{} with {:?} (len {})", start_idx, end_idx, a_chars, a_chars.len());
                        seq.splice(start_idx..end_idx, a_chars);
                    }
                } else {
                    check_bounds(start_idx, end_idx, seq.len())?;
                    seq.splice(start_idx..end_idx, std::iter::empty());
                }
                (is_subst, is_fs)
            }
            NaEdit::Del { ref_, .. } => {
                let r_len = if let Some(r) = ref_ {
                    if r.chars().all(|c| c.is_ascii_digit()) { r.parse::<i32>().unwrap_or((end_idx - start_idx) as i32) }
                    else { r.len() as i32 }
                } else {
                    (end_idx - start_idx) as i32
                };
                cds_end_i -= r_len;
                let is_fs = r_len % 3 != 0;
                check_bounds(start_idx, end_idx, seq.len())?;
                seq.splice(start_idx..end_idx, std::iter::empty());
                (false, is_fs)
            }
            NaEdit::Ins { alt: Some(alt), .. } => {
                let a_len = alt.len() as i32;
                cds_end_i += a_len;
                let is_fs = a_len % 3 != 0;
                let a_chars: Vec<char> = alt.chars().collect();
                // For an insertion c.1_2insA, start_idx=0, end_idx=2 (exclusive).
                // We want to insert at index 1 (between base 1 and 2).
                let ins_pos = (start_idx + 1).min(seq.len());
                check_bounds(ins_pos, ins_pos, seq.len())?;
                seq.splice(ins_pos..ins_pos, a_chars);
                (false, is_fs)
            }
            NaEdit::Dup { ref_, .. } => {
                let dup_chars: Vec<char> = if let Some(r) = ref_ {
                    if r.chars().all(|c| c.is_ascii_digit()) {
                        check_bounds(start_idx, end_idx, self.transcript_sequence.len())?;
                        self.transcript_sequence[start_idx..end_idx].chars().collect()
                    } else {
                        r.chars().collect()
                    }
                } else {
                    check_bounds(start_idx, end_idx, self.transcript_sequence.len())?;
                    self.transcript_sequence[start_idx..end_idx].chars().collect()
                };
                let a_len = dup_chars.len() as i32;
                cds_end_i += a_len;
                let is_fs = a_len % 3 != 0;
                // For a duplication c.1_2dup, start_idx=0, end_idx=2 (exclusive).
                // We want to insert after base 2, which is index 2.
                let ins_pos = end_idx.min(seq.len());
                check_bounds(ins_pos, ins_pos, seq.len())?;
                seq.splice(ins_pos..ins_pos, dup_chars);
                (false, is_fs)
            }
            NaEdit::Inv { .. } => {
                check_bounds(start_idx, end_idx, self.transcript_sequence.len())?;
                let sub = &self.transcript_sequence[start_idx..end_idx];
                let inv_seq = crate::utils::reverse_complement(sub);
                let inv_chars: Vec<char> = inv_seq.chars().collect();
                seq.splice(start_idx..end_idx, inv_chars);
                (false, false)
            }
            NaEdit::Repeat { min, .. } => {
                check_bounds(start_idx, end_idx, self.transcript_sequence.len())?;
                let unit = &self.transcript_sequence[start_idx..end_idx];
                // In HGVS c.7035TGGAAC[3], min=max=3 usually.
                // We assume start_idx..end_idx is the unit.
                // The total sequence becomes unit repeated 'min' times.
                let mut total_seq = String::new();
                for _ in 0..*min {
                    total_seq.push_str(unit);
                }
                let total_chars: Vec<char> = total_seq.chars().collect();
                let net_change = (total_seq.len() as i32) - (unit.len() as i32);
                cds_end_i += net_change;
                let is_fs = net_change % 3 != 0;
                // println!("DEBUG: splice at {}..{} with {:?} (len {})", start_idx, end_idx, total_chars, total_chars.len());
                seq.splice(start_idx..end_idx, total_chars);
                (false, is_fs)
            }
            NaEdit::None => (false, false),
            _ => return Err(HgvsError::UnsupportedOperation("Unsupported edit for altseq".into())),
        };

        let transcript_sequence: String = seq.iter().collect();
        let cds_start = self.cds_start_index.0 as usize;
        let cds_seq = if cds_start < transcript_sequence.len() {
            &transcript_sequence[cds_start..]
        } else { "" };
        
        let aa_sequence = crate::utils::translate_cds(cds_seq);

        Ok(AltTranscriptData {
            transcript_sequence,
            aa_sequence,
            cds_start_index: self.cds_start_index,
            cds_end_index: TranscriptPos(cds_end_i),
            protein_accession: self.protein_accession.clone(),
            is_frameshift,
            variant_start_aa,
            frameshift_start: if is_frameshift { variant_start_aa } else { None },
            is_substitution,
            is_ambiguous: false,
            c_variant: self.var_c.clone(),
        })
    }

    fn get_variant_indices(&self) -> Result<(usize, usize), HgvsError> {
        let pos = self.var_c.posedit.pos.as_ref().ok_or_else(|| HgvsError::ValidationError("Missing position".into()))?;
        let start = self.pos_to_idx(&pos.start)?;
        let mut end = if let Some(e) = &pos.end { self.pos_to_idx(e)? } else { start };
        end += 1;
        Ok((start, end))
    }

    fn pos_to_idx(&self, pos: &crate::structs::BaseOffsetPosition) -> Result<usize, HgvsError> {
        let base_idx_0 = pos.base.to_index();

        if pos.offset.is_some() && pos.offset.unwrap().0 != 0 {
             return Err(HgvsError::UnsupportedOperation("Intronic variants not yet supported in c_to_p".into()));
        }

        let idx = match pos.anchor {
            Anchor::TranscriptStart => {
                let i = base_idx_0.0;
                if i < 0 { return Err(HgvsError::ValidationError(format!("Position {} before transcript start", i))); }
                i as usize
            }
            Anchor::CdsStart => {
                let i = (self.cds_start_index.0 + base_idx_0.0) as i32;
                if i < 0 { return Err(HgvsError::ValidationError(format!("Position {} before transcript start", i))); }
                i as usize
            }
            Anchor::CdsEnd => {
                let i = (self.cds_end_index.0 + base_idx_0.0) as i32;
                if i < 0 { return Err(HgvsError::ValidationError(format!("Position {} before transcript start", i))); }
                i as usize
            }
        };
        Ok(idx)
    }
}
