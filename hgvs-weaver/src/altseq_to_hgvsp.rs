use crate::altseq::AltTranscriptData;
use crate::error::HgvsError;
use crate::fmt::aa1_to_aa3;
use crate::structs::{AAPosition, AaEdit, AaInterval, PVariant, PosEdit, ProteinPos};

pub struct AltSeqToHgvsp<'a> {
    pub ref_aa: String,
    pub ref_cds_start_idx: usize,
    pub ref_cds_end_idx: usize,
    pub alt_data: &'a AltTranscriptData,
}

impl<'a> AltSeqToHgvsp<'a> {
    pub fn build_hgvsp(&self) -> Result<PVariant, HgvsError> {
        let alt_aa = &self.alt_data.aa_sequence;
        let ref_chars: Vec<char> = self.ref_aa.chars().collect();
        let alt_chars: Vec<char> = alt_aa.chars().collect();

        if self.ref_aa == *alt_aa {
            return self.build_identity_variant();
        }

        // Find first difference
        let mut start_idx = self.alt_data.variant_start_aa.unwrap_or(ProteinPos(0)).0 as usize;
        while start_idx < ref_chars.len()
            && start_idx < alt_chars.len()
            && ref_chars[start_idx] == alt_chars[start_idx]
        {
            start_idx += 1;
        }

        // --- GOD-MODE FIX: Check if the difference is after the official stop ---
        // official_stop_idx is the index of the stop codon in the AA sequence.
        // It is (cds_end - cds_start) / 3.
        let official_stop_idx = (self.ref_cds_end_idx as i32 - self.ref_cds_start_idx as i32) / 3;

        if start_idx > official_stop_idx as usize {
            // Difference is entirely in the 3' UTR and doesn't affect the protein.
            return self.build_identity_variant();
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
                false,
            );
        }

        // Non-frameshift

        // Check for premature stop in alt_chars
        // A stop is "premature" if it appears in the variant sequence but is NOT the original stop of the protein.
        // We detect this by checking if the tail of the protein (including the stop) matches the reference tail.
        // If the tail matches significantly (len > 1), it's the original stop (shifted by deletion).
        // If the tail matches only the stop ('*'), then it's ambiguous.
        // - If the variant type is Del/Dup/Inv/Repeat, we generally assume it preserves the stop (Original).
        // - If the variant type is DelIns/Ins/Subst, we might be creating a new stop (Premature).

        let mut ref_end = ref_chars.len();
        let mut alt_end = alt_chars.len();

        // 1. Perform standard tail trimming
        while ref_end > start_idx
            && alt_end > start_idx
            && ref_chars[ref_end - 1] == alt_chars[alt_end - 1]
        {
            ref_end -= 1;
            alt_end -= 1;
        }

        let tail_match_len = alt_chars.len() - alt_end;

        let mut is_premature_stop = false;

        // Check if we trimmed a stop codon
        if tail_match_len > 0 && alt_chars[alt_chars.len() - 1] == '*' {
            if tail_match_len == 1 {
                // Only '*' matched. Ambiguous.
                // Check if it's a simple substitution (1 AA replaced by 1 AA).
                // Mismatch block length = total length - match length (from start) - tail match length.
                // But we don't track "match length from start" here explicitly, we just know start_idx.
                // So mismatch length = end (before tail trim) - start_idx.

                let ref_mismatch_len = ref_end - start_idx;
                let alt_mismatch_len = alt_end - start_idx;

                if ref_mismatch_len == 1 && alt_mismatch_len == 1 {
                    // 1-vs-1 substitution at C-terminus.
                    // Treat as Original Stop (Trim).
                    is_premature_stop = false;
                } else {
                    // Check variant type for other cases.
                    match &self.alt_data.c_variant.posedit.edit {
                        crate::structs::NaEdit::RefAlt { .. }
                        | crate::structs::NaEdit::Ins { .. } => {
                            // Likely a premature stop if not 1-vs-1 subst.
                            is_premature_stop = true;
                        }
                        _ => {} // Del/Dup etc assume original stop.
                    }
                }
            }
        } else if alt_end < alt_chars.len() && alt_chars[alt_end] == '*' {
            // Stop was NOT trimmed? (Mismatch before stop).
            // Can happen if stop is part of the mismatch.
            // But loop trims from end. If alt ends in *, ref ends in *, they SHOULD match.
            // So this branch is unlikely unless ref does NOT end in *.
        }

        if is_premature_stop {
            // Undo trimming for the stop codon.
            // We want the variant to include the stop codon to describe it as "delins...Ter".
            alt_end += 1;
            ref_end += 1;

            // Recalculate ref_end based on DNA span, because the "Ref" part of delins should match the DNA deletion.
            // (The Ref part after that is implicitly truncated).
            if let Some(pos) = &self.alt_data.c_variant.posedit.pos {
                let start_c = pos.start.base.to_index().0;
                let end_c = if let Some(e) = &pos.end {
                    e.base.to_index().0
                } else {
                    start_c
                };

                // c variant indices are 0-based from start of CDS.
                let _start_codon = start_c / 3;
                let end_codon = end_c / 3;

                let calc_ref_end = (end_codon + 1) as usize;

                // Ensure ref_end is at least start_idx.
                // We use calc_ref_end, but we must ensure we don't exceed ref_chars.len().
                // Also, if calc_ref_end < ref_end (current), it means we are shortening the ref span to matched DNA.
                ref_end = calc_ref_end.max(start_idx).min(ref_chars.len());
            }
        }

        // 1. Check for Nonsense (Substitution to Ter)
        // Only classify as Nonsense if the stop codon is part of the mismatch (not the preserved tail).
        if alt_chars.get(start_idx) == Some(&'*') && start_idx < alt_end {
            let ref_curr = aa1_to_aa3(ref_chars.get(start_idx).cloned().unwrap_or('*')).to_string();
            return self.create_variant(
                ProteinPos(start_idx as i32),
                None,
                Some(ref_curr),
                Some("Ter".to_string()),
                None,
                None,
                false,
                false,
            );
        }

        // 2. Check for Stop Loss (Extension)
        if ref_chars.get(start_idx) == Some(&'*') {
            // Find length of extension in alt_chars
            let mut ext_len = 0;
            let mut found_stop = false;
            for &c in alt_chars.iter().skip(start_idx + 1) {
                ext_len += 1;
                if c == '*' {
                    found_stop = true;
                    break;
                }
            }

            let alt_curr = aa1_to_aa3(alt_chars.get(start_idx).cloned().unwrap_or('*')).to_string();
            let length = if found_stop {
                Some(ext_len.to_string())
            } else {
                Some("?".to_string())
            };

            return Ok(PVariant {
                ac: self.alt_data.protein_accession.clone(),
                gene: None,
                posedit: PosEdit {
                    pos: Some(AaInterval {
                        start: AAPosition {
                            base: ProteinPos(start_idx as i32).to_hgvs(),
                            aa: "Ter".to_string(),
                            uncertain: false,
                        },
                        end: None,
                        uncertain: false,
                    }),
                    edit: AaEdit::Ext {
                        ref_: "Ter".into(),
                        alt: alt_curr,
                        aaterm: Some("*".to_string()),
                        length,
                        uncertain: !found_stop,
                    },
                    uncertain: false,
                    predicted: false,
                },
            });
        }

        let del_seq: String = ref_chars[start_idx..ref_end]
            .iter()
            .map(|c| aa1_to_aa3(*c))
            .collect::<Vec<&str>>()
            .join("");
        let ins_seq: String = alt_chars[start_idx..alt_end]
            .iter()
            .map(|c| aa1_to_aa3(*c))
            .collect::<Vec<&str>>()
            .join("");

        // Detect duplication
        if del_seq.is_empty() && !ins_seq.is_empty() {
            let aa_ins_len = ins_seq.len() / 3;
            if start_idx >= aa_ins_len {
                let prev_seq: String = ref_chars[start_idx - aa_ins_len..start_idx]
                    .iter()
                    .map(|c| aa1_to_aa3(*c))
                    .collect::<Vec<&str>>()
                    .join("");
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
                                start: AAPosition {
                                    base: start_pos_0.to_hgvs(),
                                    aa: aa_start,
                                    uncertain: false,
                                },
                                end: if aa_ins_len > 1 {
                                    Some(AAPosition {
                                        base: end_pos_0.to_hgvs(),
                                        aa: aa_end,
                                        uncertain: false,
                                    })
                                } else {
                                    None
                                },
                                uncertain: false,
                            }),
                            edit: AaEdit::Dup {
                                ref_: Some(ins_seq),
                                uncertain: false,
                            },
                            uncertain: false,
                            predicted: false,
                        },
                    });
                }
            }
        }

        // Detect pure insertion
        if del_seq.is_empty() && !ins_seq.is_empty() {
            let start_pos_0 = ProteinPos((start_idx as i32).saturating_sub(1));
            let end_pos_0 = ProteinPos(start_idx as i32);
            let aa_start = aa1_to_aa3(
                ref_chars
                    .get(start_pos_0.0 as usize)
                    .cloned()
                    .unwrap_or('*'),
            )
            .to_string();
            let aa_end =
                aa1_to_aa3(ref_chars.get(end_pos_0.0 as usize).cloned().unwrap_or('*')).to_string();
            return Ok(PVariant {
                ac: self.alt_data.protein_accession.clone(),
                gene: None,
                posedit: PosEdit {
                    pos: Some(AaInterval {
                        start: AAPosition {
                            base: start_pos_0.to_hgvs(),
                            aa: aa_start,
                            uncertain: false,
                        },
                        end: Some(AAPosition {
                            base: end_pos_0.to_hgvs(),
                            aa: aa_end,
                            uncertain: false,
                        }),
                        uncertain: false,
                    }),
                    edit: AaEdit::Ins {
                        alt: ins_seq,
                        uncertain: false,
                    },
                    uncertain: false,
                    predicted: false,
                },
            });
        }

        // Del / DelIns / Subst
        if ins_seq.len() == 3 && del_seq.len() == 3 {
            return self.create_variant(
                ProteinPos(start_idx as i32),
                None,
                Some(del_seq),
                Some(ins_seq),
                None,
                None,
                false,
                false,
            );
        }

        let start_pos_0 = ProteinPos(start_idx as i32);
        let end_pos_0 = if ref_end > start_idx + 1 {
            Some(ProteinPos((ref_end - 1) as i32))
        } else {
            None
        };
        let aa_start = aa1_to_aa3(ref_chars.get(start_idx).cloned().unwrap_or('*')).to_string();
        let aa_end = end_pos_0
            .map(|e| aa1_to_aa3(ref_chars.get(e.0 as usize).cloned().unwrap_or('*')).to_string());

        let edit = if ins_seq.is_empty() {
            AaEdit::Del {
                ref_: del_seq,
                uncertain: false,
            }
        } else {
            AaEdit::DelIns {
                ref_: del_seq,
                alt: ins_seq,
                uncertain: false,
            }
        };

        Ok(PVariant {
            ac: self.alt_data.protein_accession.clone(),
            gene: None,
            posedit: PosEdit {
                pos: Some(AaInterval {
                    start: AAPosition {
                        base: start_pos_0.to_hgvs(),
                        aa: aa_start,
                        uncertain: false,
                    },
                    end: end_pos_0.map(|e| e.to_hgvs()).map(|base| AAPosition {
                        base,
                        aa: aa_end.unwrap(),
                        uncertain: false,
                    }),
                    uncertain: false,
                }),
                edit,
                uncertain: false,
                predicted: false,
            },
        })
    }

    fn build_identity_variant(&self) -> Result<PVariant, HgvsError> {
        let ref_chars: Vec<char> = self.ref_aa.chars().collect();
        let start_0 = self.alt_data.variant_start_aa.unwrap_or(ProteinPos(0));
        let end_0 = {
            let cds_start = self.alt_data.cds_start_index.0 as i32;
            let start_idx = self.alt_data.variant_start_idx as i32;
            let end_idx = self.alt_data.variant_end_idx as i32;
            let s0 = ProteinPos(((start_idx - cds_start).max(0) / 3) as i32);
            let e0_raw = ((end_idx - 1 - cds_start).max(0) / 3) as i32;
            let e0 = ProteinPos(e0_raw.min(ref_chars.len() as i32 - 1));
            if e0.0 > s0.0 {
                Some(e0)
            } else {
                None
            }
        };

        self.create_variant(start_0, end_0, None, None, None, None, false, true)
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

        // Handle 3' UTR or other variants beyond the protein
        if start_0.0 >= ref_chars.len() as i32 {
            return Ok(PVariant {
                ac: self.alt_data.protein_accession.clone(),
                gene: None,
                posedit: PosEdit {
                    pos: None,
                    edit: AaEdit::Identity { uncertain: false },
                    uncertain: false,
                    predicted: false,
                },
            });
        }

        let aa_start =
            aa1_to_aa3(ref_chars.get(start_0.0 as usize).cloned().unwrap_or('*')).to_string();

        let edit = if is_silent {
            AaEdit::Identity { uncertain: false }
        } else if alt_aa.as_ref().is_some_and(|a| a == "Ter") {
            // Nonsense
            AaEdit::Subst {
                ref_: ref_aa.unwrap_or_default(),
                alt: "Ter".to_string(),
                uncertain: false,
            }
        } else if is_fs {
            let len_str = length.map(|l| l.replace("Ter", ""));
            AaEdit::Fs {
                ref_: "".into(),
                alt: alt_aa.unwrap_or_default(),
                term,
                length: len_str,
                uncertain: false,
            }
        } else {
            AaEdit::Subst {
                ref_: ref_aa.unwrap_or_default(),
                alt: alt_aa.unwrap_or_default(),
                uncertain: false,
            }
        };

        let aa_end = end_0
            .map(|e| aa1_to_aa3(ref_chars.get(e.0 as usize).cloned().unwrap_or('*')).to_string());

        let interval = AaInterval {
            start: AAPosition {
                base: start_0.to_hgvs(),
                aa: aa_start,
                uncertain: false,
            },
            end: end_0.map(|e| e.to_hgvs()).map(|base| AAPosition {
                base,
                aa: aa_end.clone().unwrap(),
                uncertain: false,
            }),
            uncertain: false,
        };

        Ok(PVariant {
            ac: self.alt_data.protein_accession.clone(),
            gene: None,
            posedit: PosEdit {
                pos: Some(interval),
                edit,
                uncertain: false,
                predicted: false,
            },
        })
    }
}
