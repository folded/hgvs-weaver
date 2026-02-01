use crate::error::HgvsError;
use crate::structs::{PVariant, AAPosition, AaInterval, AaEdit, PosEdit, ProteinPos};
use crate::altseq::AltTranscriptData;
use crate::fmt::aa1_to_aa3;

pub struct AltSeqToHgvsp<'a> {
    pub ref_aa: String,
    pub alt_data: &'a AltTranscriptData,
}

impl<'a> AltSeqToHgvsp<'a> {
    pub fn build_hgvsp(&self) -> Result<PVariant, HgvsError> {
        let alt_aa = &self.alt_data.aa_sequence;
        let ref_chars: Vec<char> = self.ref_aa.chars().collect();
        let alt_chars: Vec<char> = alt_aa.chars().collect();

        if self.ref_aa == *alt_aa {
            let c_pos = self.alt_data.c_variant.posedit.pos.as_ref().ok_or_else(|| HgvsError::ValidationError("Missing position".into()))?;
            let start_0 = ProteinPos(c_pos.start.base.to_index().0 / 3);
            let end_0 = if let Some(e) = &c_pos.end {
                Some(ProteinPos(e.base.to_index().0 / 3))
            } else {
                None
            };
            return self.create_variant(start_0, end_0, None, None, None, None, false, true);
        }

        // Find first difference
        let mut start_idx = self.alt_data.variant_start_aa.unwrap_or(ProteinPos(0)).0 as usize;
        while start_idx < ref_chars.len() && start_idx < alt_chars.len() && ref_chars[start_idx] == alt_chars[start_idx] {
            start_idx += 1;
        }

        if self.alt_data.is_frameshift {
            let ref_curr = aa1_to_aa3(ref_chars.get(start_idx).cloned().unwrap_or('*')).to_string();
            let alt_curr = aa1_to_aa3(alt_chars.get(start_idx).cloned().unwrap_or('*')).to_string();

            // Find first stop in alt_aa starting from start_idx
            let mut stop_idx = None;
            for (i, &c) in alt_chars.iter().enumerate().skip(start_idx) {
                if c == '*' {
                    stop_idx = Some(i);
                    break;
                }
            }

            let term = Some("Ter".to_string());
            let (length, _uncertain) = if let Some(idx) = stop_idx {
                ((idx - start_idx + 1).to_string(), false)
            } else {
                ("?".to_string(), true)
            };

            return self.create_variant(
                ProteinPos(start_idx as i32),
                None,
                Some(ref_curr),
                Some(alt_curr),
                term,
                Some(length),
                true,
                false
            );
        }

        // Non-frameshift
        let mut ref_end = ref_chars.len();
        let mut alt_end = alt_chars.len();
        while ref_end > start_idx && alt_end > start_idx && ref_chars[ref_end-1] == alt_chars[alt_end-1] {
            ref_end -= 1;
            alt_end -= 1;
        }

        // 1. Check for Nonsense that might be an in-frame deletion-insertion
        if alt_chars.get(start_idx) == Some(&'*') {
             let c_len = self.alt_data.c_variant.posedit.pos.as_ref()
                .map(|p| p.length().unwrap_or(0))
                .unwrap_or(0) as usize;

             if c_len > 0 && c_len % 3 == 0 {
                 let aa_del_len = c_len / 3;
                 if aa_del_len >= 1 {
                     let ref_end_idx = start_idx + aa_del_len;
                     if ref_end_idx <= ref_chars.len() {
                         let start_pos_0 = ProteinPos(start_idx as i32);
                         let end_pos_0 = ProteinPos((ref_end_idx - 1) as i32);
                         let del_seq: String = ref_chars[start_idx..ref_end_idx].iter().map(|c| aa1_to_aa3(*c)).collect::<Vec<&str>>().join("");
                         return Ok(PVariant {
                            ac: self.alt_data.protein_accession.clone(),
                            gene: None,
                            posedit: PosEdit {
                                pos: Some(AaInterval {
                                    start: AAPosition { base: start_pos_0.to_hgvs(), aa: aa1_to_aa3(ref_chars.get(start_idx).cloned().unwrap_or('*')).to_string(), uncertain: false },
                                    end: Some(AAPosition { base: end_pos_0.to_hgvs(), aa: aa1_to_aa3(ref_chars.get(end_pos_0.0 as usize).cloned().unwrap_or('*')).to_string(), uncertain: false }),
                                    uncertain: false,
                                }),
                                edit: AaEdit::DelIns { ref_: del_seq, alt: "Ter".to_string(), uncertain: false },
                                uncertain: false,
                                predicted: false,
                            }
                        });
                     }
                 }
             }
        }

        let del_seq: String = ref_chars[start_idx..ref_end].iter().map(|c| aa1_to_aa3(*c)).collect::<Vec<&str>>().join("");
        let ins_seq: String = alt_chars[start_idx..alt_end].iter().map(|c| aa1_to_aa3(*c)).collect::<Vec<&str>>().join("");

        if del_seq.len() == 3 && ins_seq == "Ter" { // 3 letters for 1 AA
             return self.create_variant(
                ProteinPos(start_idx as i32),
                None,
                Some(del_seq),
                Some("Ter".to_string()),
                None,
                None,
                false,
                false
            );
        }

        if del_seq.len() == 3 && ins_seq.len() == 3 {
            return self.create_variant(
                ProteinPos(start_idx as i32),
                None,
                Some(del_seq),
                Some(ins_seq),
                None,
                None,
                false,
                false
            );
        }

        // Detect duplication
        if del_seq.is_empty() && !ins_seq.is_empty() {
            let aa_ins_len = ins_seq.len() / 3;
            if start_idx >= aa_ins_len {
                let prev_seq: String = ref_chars[start_idx - aa_ins_len..start_idx].iter().map(|c| aa1_to_aa3(*c)).collect::<Vec<&str>>().join("");
                if prev_seq == ins_seq {
                    let start_pos_0 = ProteinPos((start_idx - aa_ins_len) as i32);
                    let end_pos_0 = ProteinPos((start_idx - 1) as i32);
                    let aa_start = aa1_to_aa3(ref_chars[start_pos_0.0 as usize]).to_string();
                    let aa_end = aa1_to_aa3(ref_chars[end_pos_0.0 as usize]).to_string();

                    return Ok(PVariant {
                        ac: self.alt_data.protein_accession.clone(),
                        gene: None,
                        posedit: PosEdit {
                            pos: Some(AaInterval {
                                start: AAPosition { base: start_pos_0.to_hgvs(), aa: aa_start, uncertain: false },
                                end: if aa_ins_len > 1 {
                                    Some(AAPosition { base: end_pos_0.to_hgvs(), aa: aa_end, uncertain: false })
                                } else {
                                    None
                                },
                                uncertain: false,
                            }),
                            edit: AaEdit::Dup { ref_: Some(ins_seq), uncertain: false },
                            uncertain: false,
                            predicted: false,
                        }
                    });
                }
            }
        }

        // Detect pure insertion
        if del_seq.is_empty() && !ins_seq.is_empty() {
            let start_pos_0 = ProteinPos((start_idx as i32).saturating_sub(1));
            let end_pos_0 = ProteinPos(start_idx as i32);
            let aa_start = aa1_to_aa3(ref_chars.get(start_pos_0.0 as usize).cloned().unwrap_or('*')).to_string();
            let aa_end = aa1_to_aa3(ref_chars.get(end_pos_0.0 as usize).cloned().unwrap_or('*')).to_string();
            return Ok(PVariant {
                ac: self.alt_data.protein_accession.clone(),
                gene: None,
                posedit: PosEdit {
                    pos: Some(AaInterval {
                        start: AAPosition { base: start_pos_0.to_hgvs(), aa: aa_start, uncertain: false },
                        end: Some(AAPosition { base: end_pos_0.to_hgvs(), aa: aa_end, uncertain: false }),
                        uncertain: false,
                    }),
                    edit: AaEdit::Ins { alt: ins_seq, uncertain: false },
                    uncertain: false,
                    predicted: false,
                }
            });
        }

        // Del / DelIns
        let start_pos_0 = ProteinPos(start_idx as i32);
        let end_pos_0 = if ref_end > start_idx + 1 { Some(ProteinPos((ref_end - 1) as i32)) } else { None };
        let aa_start = aa1_to_aa3(ref_chars.get(start_idx).cloned().unwrap_or('*')).to_string();
        let aa_end = end_pos_0.map(|e| aa1_to_aa3(ref_chars.get(e.0 as usize).cloned().unwrap_or('*')).to_string());

        let edit = if ins_seq.is_empty() {
            AaEdit::Del { ref_: del_seq, uncertain: false }
        } else {
            AaEdit::DelIns { ref_: del_seq, alt: ins_seq, uncertain: false }
        };

        Ok(PVariant {
            ac: self.alt_data.protein_accession.clone(),
            gene: None,
            posedit: PosEdit {
                pos: Some(AaInterval {
                    start: AAPosition { base: start_pos_0.to_hgvs(), aa: aa_start, uncertain: false },
                    end: end_pos_0.map(|e| e.to_hgvs()).map(|base| AAPosition { base, aa: aa_end.unwrap(), uncertain: false }),
                    uncertain: false,
                }),
                edit,
                uncertain: false,
                predicted: false,
            }
        })
    }

    #[allow(clippy::too_many_arguments)]
    fn create_variant(
        &self,
        start_0: ProteinPos,
        end_0: Option<ProteinPos>,
        ref_aa: Option<String>,
        alt_aa: Option<String>,
        term: Option<String>,
        length: Option<String>,
        is_fs: bool,
        is_silent: bool,
    ) -> Result<PVariant, HgvsError> {
        let ref_chars: Vec<char> = self.ref_aa.chars().collect();
        let aa_start = aa1_to_aa3(ref_chars.get(start_0.0 as usize).cloned().unwrap_or('*')).to_string();

        let edit = if is_silent {
            AaEdit::Identity { uncertain: false }
        } else if alt_aa.as_ref().is_some_and(|a| a == "Ter") {
            // Nonsense
            AaEdit::Subst { ref_: ref_aa.unwrap_or_default(), alt: "Ter".to_string(), uncertain: false }
        } else if is_fs {
            let len_str = length.map(|l| l.replace("Ter", ""));
            AaEdit::Fs { ref_: "".into(), alt: alt_aa.unwrap_or_default(), term, length: len_str, uncertain: false }
        } else {
            AaEdit::Subst { ref_: ref_aa.unwrap_or_default(), alt: alt_aa.unwrap_or_default(), uncertain: false }
        };

        let aa_end = end_0.map(|e| aa1_to_aa3(ref_chars.get(e.0 as usize).cloned().unwrap_or('*')).to_string());

        Ok(PVariant {
            ac: self.alt_data.protein_accession.clone(),
            gene: None,
            posedit: PosEdit {
                pos: Some(AaInterval {
                    start: AAPosition { base: start_0.to_hgvs(), aa: aa_start, uncertain: false },
                    end: end_0.map(|e| e.to_hgvs()).map(|base| AAPosition { base, aa: aa_end.unwrap(), uncertain: false }),
                    uncertain: false,
                }),
                edit,
                uncertain: false,
                predicted: false,
            }
        })
    }
}
