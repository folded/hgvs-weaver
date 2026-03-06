use crate::altseq::AltSeqBuilder;
use crate::altseq_to_hgvsp::AltSeqToHgvsp;
use crate::data::{DataProvider, IdentifierKind, IdentifierType, Transcript, TranscriptSearch};
use crate::error::HgvsError;
use crate::sequence::{MemSequence, RevCompSequence, Sequence, TranslatedSequence};
use crate::structs::{
    BaseOffsetInterval, BaseOffsetPosition, CVariant, GVariant, NVariant, PVariant,
};
use crate::transcript_mapper::TranscriptMapper;

/// Converts a 0-based transcript index to a fresh `BaseOffsetPosition` via `n_to_c`.
///
/// The returned position has `uncertain = false`; callers that need to propagate
/// an existing uncertainty flag should overwrite that field after calling this.
fn n_to_c_position(am: &TranscriptMapper, n: i32) -> Result<BaseOffsetPosition, HgvsError> {
    let (c_pos, offset, anchor) = am.n_to_c(crate::coords::TranscriptPos(n))?;
    Ok(BaseOffsetPosition {
        base: c_pos.to_hgvs(),
        offset: if offset.0 != 0 { Some(offset) } else { None },
        anchor,
        uncertain: false,
    })
}

fn ins_anchor_and_end(start: usize, end: usize, is_ins: bool) -> (usize, usize) {
    let actual_end = if is_ins { end - 1 } else { end };
    let anchor = if is_ins { actual_end } else { start };
    (anchor, actual_end)
}

fn make_base_offset_position(
    base: crate::coords::HgvsTranscriptPos,
    combined_offset: i32,
    anchor: crate::coords::Anchor,
) -> BaseOffsetPosition {
    BaseOffsetPosition {
        base,
        offset: if combined_offset != 0 {
            Some(crate::structs::IntronicOffset(combined_offset))
        } else {
            None
        },
        anchor,
        uncertain: false,
    }
}

fn make_simple_position(base: crate::coords::HgvsGenomicPos) -> crate::structs::SimplePosition {
    crate::structs::SimplePosition {
        base,
        end: None,
        uncertain: false,
    }
}

fn apply_strand_complement(
    edit: crate::edits::NaEdit,
    strand: crate::data::Strand,
) -> crate::edits::NaEdit {
    if strand == crate::data::Strand::Minus {
        edit.map_sequence(|s| {
            let seq = MemSequence(s.to_string());
            RevCompSequence { inner: &seq }.to_string()
        })
    } else {
        edit
    }
}

fn extract_edit_sequences(
    hdp: &dyn DataProvider,
    ac: &str,
    start: usize,
    end: usize,
    edit: &crate::edits::NaEdit,
) -> Result<Option<(String, String)>, HgvsError> {
    let result = match edit {
        crate::edits::NaEdit::RefAlt { ref_, alt, .. } => Some((
            ref_.clone().unwrap_or_default(),
            alt.clone().unwrap_or_default(),
        )),
        crate::edits::NaEdit::Del { ref_: Some(s), .. } => Some((s.clone(), String::new())),
        crate::edits::NaEdit::Del { ref_: None, .. } => Some((String::new(), String::new())),
        crate::edits::NaEdit::Ins { alt: Some(s), .. } => Some((String::new(), s.clone())),
        crate::edits::NaEdit::Ins { alt: None, .. } => Some((String::new(), String::new())),
        crate::edits::NaEdit::Dup { ref_: Some(s), .. } => Some((s.clone(), String::new())),
        crate::edits::NaEdit::Dup { ref_: None, .. } => Some((String::new(), String::new())),
        crate::edits::NaEdit::Repeat { ref_, max, .. } => {
            let r = if let Some(r) = ref_ {
                r.clone()
            } else {
                hdp.get_seq(
                    ac,
                    start as i32,
                    end as i32,
                    IdentifierType::GenomicAccession,
                )?
            };
            let a = r.repeat(*max as usize);
            Some((r, a))
        }
        crate::edits::NaEdit::Inv { .. } => {
            let r = hdp.get_seq(
                ac,
                start as i32,
                end as i32,
                IdentifierType::GenomicAccession,
            )?;
            let a = crate::sequence::rev_comp(&r);
            Some((r, a))
        }
        _ => None,
    };
    Ok(result)
}

fn checked_usize(val: i32, context: &str) -> Result<usize, HgvsError> {
    if val < 0 {
        Err(HgvsError::ValidationError(format!(
            "Invalid negative value for {}: {}",
            context, val
        )))
    } else {
        Ok(val as usize)
    }
}

/// High-level mapper for transforming variants between coordinate systems.
pub struct VariantMapper<'a> {
    /// Data provider used to retrieve transcript and sequence information.
    pub hdp: &'a dyn DataProvider,
}

impl<'a> VariantMapper<'a> {
    /// Creates a new `VariantMapper` with the given data provider.
    pub fn new(hdp: &'a dyn DataProvider) -> Self {
        VariantMapper { hdp }
    }

    /// Transforms a genomic variant (`g.`) to a coding cDNA variant (`c.`).
    pub fn g_to_c(&self, var_g: &GVariant, transcript_ac: &str) -> Result<CVariant, HgvsError> {
        let transcript = self.hdp.get_transcript(transcript_ac, Some(&var_g.ac))?;
        let am = TranscriptMapper::new(transcript)?;

        let pos = var_g
            .posedit
            .pos
            .as_ref()
            .ok_or_else(|| HgvsError::ValidationError("Missing genomic position".into()))?;
        let g_start_0 = pos.start.base.to_index();
        let (mut n_pos, mut offset) = am.g_to_n(g_start_0)?;

        if let Some(end_g_simple) = &pos.end {
            let g_end_0 = end_g_simple.base.to_index();
            let (mut n_pos_e, mut offset_e) = am.g_to_n(g_end_0)?;

            if n_pos.0 > n_pos_e.0 {
                std::mem::swap(&mut n_pos, &mut n_pos_e);
                std::mem::swap(&mut offset, &mut offset_e);
            }

            let (c_pos_index, c_offset, anchor) = am.n_to_c(n_pos)?;
            let (c_pos_e_index, c_offset_e, anchor_e) = am.n_to_c(n_pos_e)?;

            let pos_c =
                make_base_offset_position(c_pos_index.to_hgvs(), c_offset.0 + offset.0, anchor);

            let pos_c_e = make_base_offset_position(
                c_pos_e_index.to_hgvs(),
                c_offset_e.0 + offset_e.0,
                anchor_e,
            );

            let edit = apply_strand_complement(var_g.posedit.edit.clone(), am.transcript.strand());

            return Ok(CVariant {
                ac: transcript_ac.to_string(),
                gene: var_g.gene.clone(),
                posedit: crate::structs::PosEdit {
                    pos: Some(crate::structs::BaseOffsetInterval {
                        start: pos_c,
                        end: Some(pos_c_e),
                        uncertain: false,
                    }),
                    edit,
                    uncertain: var_g.posedit.uncertain,
                    predicted: var_g.posedit.predicted,
                },
            });
        }

        let (c_pos_index, c_offset, anchor) = am.n_to_c(n_pos)?;
        let pos_c = make_base_offset_position(c_pos_index.to_hgvs(), c_offset.0 + offset.0, anchor);

        let edit = apply_strand_complement(var_g.posedit.edit.clone(), am.transcript.strand());

        Ok(CVariant {
            ac: transcript_ac.to_string(),
            gene: var_g.gene.clone(),
            posedit: crate::structs::PosEdit {
                pos: Some(crate::structs::BaseOffsetInterval {
                    start: pos_c,
                    end: None,
                    uncertain: false,
                }),
                edit,
                uncertain: var_g.posedit.uncertain,
                predicted: var_g.posedit.predicted,
            },
        })
    }

    /// Transforms a coding cDNA variant (`c.`) to a genomic variant (`g.`).
    pub fn c_to_g(
        &self,
        var_c: &CVariant,
        reference_ac: Option<&str>,
    ) -> Result<GVariant, HgvsError> {
        let transcript = self.hdp.get_transcript(&var_c.ac, reference_ac)?;
        let am = TranscriptMapper::new(transcript)?;

        let pos = var_c
            .posedit
            .pos
            .as_ref()
            .ok_or_else(|| HgvsError::ValidationError("Missing cDNA position".into()))?;
        let n_pos = am.c_to_n(pos.start.base.to_index(), pos.start.anchor)?;
        let g_pos = am.n_to_g(
            n_pos,
            pos.start
                .offset
                .unwrap_or(crate::structs::IntronicOffset(0)),
        )?;

        if let Some(end_c) = &pos.end {
            let n_pos_e = am.c_to_n(end_c.base.to_index(), end_c.anchor)?;
            let g_pos_e = am.n_to_g(
                n_pos_e,
                end_c.offset.unwrap_or(crate::structs::IntronicOffset(0)),
            )?;
            let mut pos_g = make_simple_position(g_pos.to_hgvs());
            let mut pos_g_e = make_simple_position(g_pos_e.to_hgvs());

            if pos_g.base.0 > pos_g_e.base.0 {
                std::mem::swap(&mut pos_g, &mut pos_g_e);
            }

            let edit = apply_strand_complement(var_c.posedit.edit.clone(), am.transcript.strand());

            return Ok(GVariant {
                ac: reference_ac
                    .unwrap_or_else(|| am.transcript.reference_accession())
                    .to_string(),
                gene: var_c.gene.clone(),
                posedit: crate::structs::PosEdit {
                    pos: Some(crate::structs::SimpleInterval {
                        start: pos_g,
                        end: Some(pos_g_e),
                        uncertain: false,
                    }),
                    edit,
                    uncertain: var_c.posedit.uncertain,
                    predicted: var_c.posedit.predicted,
                },
            });
        }

        let pos_g = make_simple_position(g_pos.to_hgvs());
        let edit = apply_strand_complement(var_c.posedit.edit.clone(), am.transcript.strand());

        Ok(GVariant {
            ac: reference_ac
                .unwrap_or_else(|| am.transcript.reference_accession())
                .to_string(),
            gene: var_c.gene.clone(),
            posedit: crate::structs::PosEdit {
                pos: Some(crate::structs::SimpleInterval {
                    start: pos_g,
                    end: None,
                    uncertain: false,
                }),
                edit,
                uncertain: var_c.posedit.uncertain,
                predicted: var_c.posedit.predicted,
            },
        })
    }

    /// Transforms a non-coding cDNA variant (`n.`) to a genomic variant (`g.`).
    pub fn n_to_g(
        &self,
        var_n: &crate::structs::NVariant,
        reference_ac: Option<&str>,
    ) -> Result<GVariant, HgvsError> {
        let transcript = self.hdp.get_transcript(&var_n.ac, reference_ac)?;
        let am = TranscriptMapper::new(transcript)?;

        let pos = var_n
            .posedit
            .pos
            .as_ref()
            .ok_or_else(|| HgvsError::ValidationError("Missing cDNA position".into()))?;
        let n_pos = am.c_to_n(pos.start.base.to_index(), pos.start.anchor)?;
        let g_pos = am.n_to_g(
            n_pos,
            pos.start
                .offset
                .unwrap_or(crate::structs::IntronicOffset(0)),
        )?;

        if let Some(end_n) = &pos.end {
            let n_pos_e = am.c_to_n(end_n.base.to_index(), end_n.anchor)?;
            let g_pos_e = am.n_to_g(
                n_pos_e,
                end_n.offset.unwrap_or(crate::structs::IntronicOffset(0)),
            )?;
            let mut pos_g = make_simple_position(g_pos.to_hgvs());
            let mut pos_g_e = make_simple_position(g_pos_e.to_hgvs());

            if pos_g.base.0 > pos_g_e.base.0 {
                std::mem::swap(&mut pos_g, &mut pos_g_e);
            }

            let edit = apply_strand_complement(var_n.posedit.edit.clone(), am.transcript.strand());

            return Ok(GVariant {
                ac: reference_ac
                    .unwrap_or_else(|| am.transcript.reference_accession())
                    .to_string(),
                gene: var_n.gene.clone(),
                posedit: crate::structs::PosEdit {
                    pos: Some(crate::structs::SimpleInterval {
                        start: pos_g,
                        end: Some(pos_g_e),
                        uncertain: false,
                    }),
                    edit,
                    uncertain: var_n.posedit.uncertain,
                    predicted: var_n.posedit.predicted,
                },
            });
        }

        let pos_g = make_simple_position(g_pos.to_hgvs());
        let edit = apply_strand_complement(var_n.posedit.edit.clone(), am.transcript.strand());

        Ok(GVariant {
            ac: reference_ac
                .unwrap_or_else(|| am.transcript.reference_accession())
                .to_string(),
            gene: var_n.gene.clone(),
            posedit: crate::structs::PosEdit {
                pos: Some(crate::structs::SimpleInterval {
                    start: pos_g,
                    end: None,
                    uncertain: false,
                }),
                edit,
                uncertain: var_n.posedit.uncertain,
                predicted: var_n.posedit.predicted,
            },
        })
    }

    /// Discovers all possible cDNA consequences for a genomic variant.
    pub fn g_to_c_all(
        &self,
        var_g: &GVariant,
        searcher: &dyn TranscriptSearch,
    ) -> Result<Vec<CVariant>, HgvsError> {
        let pos = var_g
            .posedit
            .pos
            .as_ref()
            .ok_or_else(|| HgvsError::ValidationError("Missing position".into()))?;
        let start_0 = pos.start.base.to_index().0;
        let end_0 = pos
            .end
            .as_ref()
            .map_or(start_0 + 1, |e| e.base.to_index().0 + 1);

        let transcripts = searcher.get_transcripts_for_region(&var_g.ac, start_0, end_0)?;
        if transcripts.is_empty() {
            return Err(HgvsError::ValidationError(format!(
                "No transcripts found for region {}:{}-{}",
                var_g.ac, start_0, end_0
            )));
        }

        let mut results = Vec::new();
        let mut errors = Vec::new();

        for tx_ac in transcripts {
            match self.g_to_c(var_g, &tx_ac) {
                Ok(vc) => results.push(vc),
                Err(e) => errors.push(format!("{}: {}", tx_ac, e)),
            }
        }

        if results.is_empty() && !errors.is_empty() {
            return Err(HgvsError::ValidationError(format!(
                "All mapping attempts failed: {}",
                errors.join("; ")
            )));
        }

        Ok(results)
    }

    /// Transforms a coding cDNA variant (`c.`) to a protein variant (`p.`).
    pub fn c_to_p(
        &self,
        var_c: &CVariant,
        protein_ac: Option<&str>,
    ) -> Result<PVariant, HgvsError> {
        let transcript_ac = &var_c.ac;
        let pro_ac_str = if let Some(ac) = protein_ac {
            ac.to_string()
        } else {
            self.hdp
                .get_symbol_accessions(
                    transcript_ac,
                    IdentifierKind::Transcript,
                    IdentifierKind::Protein,
                )?
                .first()
                .ok_or_else(|| {
                    HgvsError::ValidationError(format!(
                        "No protein accession found for {}",
                        transcript_ac
                    ))
                })?
                .1
                .clone()
        };

        // Handle intronic variants by returning p.?
        if let Some(pos) = &var_c.posedit.pos {
            let has_offset = pos.start.offset.is_some_and(|o| o.0 != 0)
                || pos
                    .end
                    .as_ref()
                    .is_some_and(|e| e.offset.is_some_and(|o| o.0 != 0));

            if has_offset {
                return Ok(PVariant {
                    ac: pro_ac_str,
                    gene: var_c.gene.clone(),
                    posedit: crate::structs::PosEdit {
                        pos: None,
                        edit: crate::edits::AaEdit::Special {
                            value: "?".to_string(),
                            uncertain: false,
                        },
                        uncertain: false,
                        predicted: false,
                    },
                });
            }

            // 5'UTR variants (negative c. position) → p.?
            let start_base = pos.start.base.0;
            if pos.start.anchor == crate::coords::Anchor::CdsStart && start_base < 0 {
                return Ok(PVariant {
                    ac: pro_ac_str,
                    gene: var_c.gene.clone(),
                    posedit: crate::structs::PosEdit {
                        pos: None,
                        edit: crate::edits::AaEdit::Special {
                            value: "?".to_string(),
                            uncertain: false,
                        },
                        uncertain: false,
                        predicted: false,
                    },
                });
            }
        }

        let transcript = self.hdp.get_transcript(transcript_ac, None)?;
        let ref_seq = self.hdp.get_seq(
            transcript_ac,
            0,
            -1,
            IdentifierKind::Transcript.into_identifier_type(),
        )?;

        let cds_start_tx = transcript
            .cds_start_index()
            .ok_or_else(|| HgvsError::ValidationError("Missing CDS start".into()))?;
        let cds_end_tx = transcript
            .cds_end_index()
            .ok_or_else(|| HgvsError::ValidationError("Missing CDS end".into()))?;
        let cds_start_idx = checked_usize(cds_start_tx.0, "CDS start")?;
        let cds_end_idx = checked_usize(cds_end_tx.0, "CDS end")?;

        let ref_seq_obj = MemSequence(ref_seq);

        if ref_seq_obj.len() < cds_end_idx {
            return Err(HgvsError::ValidationError(format!(
                "Transcript sequence too short (len={}, expected at least {})",
                ref_seq_obj.len(),
                cds_end_idx
            )));
        }

        if cds_start_idx > ref_seq_obj.len() {
            return Err(HgvsError::ValidationError(format!(
                "CDS start {} out of sequence bounds {}",
                cds_start_idx,
                ref_seq_obj.len()
            )));
        }

        // Translate from the already-fetched transcript sequence (avoids a second provider call).
        let cds_slice = ref_seq_obj.slice(cds_start_idx, ref_seq_obj.len());
        let ref_aa = TranslatedSequence { inner: &cds_slice }.to_string();

        let builder = AltSeqBuilder {
            var_c,
            transcript_sequence: &ref_seq_obj,
            cds_start_index: cds_start_tx,
            cds_end_index: cds_end_tx,
            protein_accession: pro_ac_str,
        };
        let alt_data = builder.build_altseq()?;

        let hgvsp_builder = AltSeqToHgvsp {
            ref_aa,
            ref_cds_start_idx: cds_start_idx,
            ref_cds_end_idx: cds_end_idx,
            alt_data: &alt_data,
        };
        let mut var_p = hgvsp_builder.build_hgvsp()?;
        var_p.posedit.predicted = true;
        Ok(var_p)
    }

    /// Back-converts a single amino acid substitution (p.Xxx###Yyy) to a coding variant (c.).
    ///
    /// Only handles `AaEdit::Subst` currently. For ambiguous back-conversions (multiple
    /// codons with the same minimum nucleotide distance), the first alphabetically is chosen
    /// and `is_unique` is set to false.
    ///
    /// Returns `(CVariant, is_unique)`.
    pub fn p_to_c(
        &self,
        var_p: &PVariant,
        transcript_ac: Option<&str>,
    ) -> Result<(CVariant, bool), HgvsError> {
        let ac = &var_p.ac;

        // Extract position and edit
        let pos = var_p
            .posedit
            .pos
            .as_ref()
            .ok_or_else(|| HgvsError::ValidationError("Missing protein position".into()))?;

        let edit = &var_p.posedit.edit;

        // Only handle Subst for now
        let alt_aa_str = match edit {
            crate::edits::AaEdit::Subst { alt, .. } => alt.clone(),
            _ => {
                return Err(HgvsError::UnsupportedOperation(
                    "p_to_c only supports single amino acid substitutions".into(),
                ));
            }
        };

        let raw_pos = pos.start.base.0;
        if raw_pos <= 0 {
            return Err(HgvsError::ValidationError(format!(
                "Protein position {raw_pos} is not valid (must be >= 1)"
            )));
        }
        let aa_pos = checked_usize(raw_pos, "protein position")?; // 1-based
        let ref_aa_str = &pos.start.aa;

        // Convert 3-letter or 1-letter to single char
        let ref_aa_1 = crate::utils::aa3_to_aa1(ref_aa_str);
        let alt_aa_1 = crate::utils::aa3_to_aa1(&alt_aa_str);

        let _ref_aa = ref_aa_1
            .chars()
            .next()
            .ok_or_else(|| HgvsError::ValidationError("Invalid reference AA".into()))?;
        let alt_aa = alt_aa_1
            .chars()
            .next()
            .ok_or_else(|| HgvsError::ValidationError("Invalid alternate AA".into()))?;

        // Resolve transcript accession
        let tx_ac = if let Some(ta) = transcript_ac {
            ta.to_string()
        } else {
            self.hdp
                .get_symbol_accessions(ac, IdentifierKind::Protein, IdentifierKind::Transcript)?
                .first()
                .ok_or_else(|| {
                    HgvsError::ValidationError(format!("No transcript accession found for {}", ac))
                })?
                .1
                .clone()
        };

        // Get transcript
        let transcript = self.hdp.get_transcript(&tx_ac, None)?;
        let cds_start = transcript
            .cds_start_index()
            .ok_or_else(|| HgvsError::ValidationError("No CDS start for transcript".into()))?
            .0 as usize;

        // Get transcript sequence
        let tx_seq_str = self.hdp.get_seq(
            &tx_ac,
            0,
            -1,
            IdentifierKind::Transcript.into_identifier_type(),
        )?;

        // Calculate codon position (0-based in transcript)
        let codon_start = cds_start + (aa_pos - 1) * 3;
        let codon_end = codon_start + 3;
        if codon_end > tx_seq_str.len() {
            return Err(HgvsError::ValidationError(format!(
                "Codon position {}-{} out of range for sequence length {}",
                codon_start,
                codon_end,
                tx_seq_str.len()
            )));
        }

        let ref_codon = tx_seq_str[codon_start..codon_end].to_uppercase();
        let ref_codon_bytes = ref_codon.as_bytes();

        // Get all codons for the target AA
        let alt_codons = crate::utils::codons_for_aa(alt_aa);
        if alt_codons.is_empty() {
            return Err(HgvsError::ValidationError(format!(
                "No codons found for amino acid '{}'",
                alt_aa
            )));
        }

        // Score each candidate by nucleotide differences, pick minimum
        let mut scored: Vec<(usize, &str)> = alt_codons
            .iter()
            .map(|codon| {
                let diffs = codon
                    .as_bytes()
                    .iter()
                    .zip(ref_codon_bytes.iter())
                    .filter(|(a, b)| a != b)
                    .count();
                (diffs, *codon)
            })
            .collect();
        scored.sort();

        let best_diffs = scored[0].0;
        let best_codons: Vec<&str> = scored
            .iter()
            .filter(|(d, _)| *d == best_diffs)
            .map(|(_, c)| *c)
            .collect();
        let is_unique = best_codons.len() == 1;
        let best_codon = best_codons[0];

        // Build the c. variant for the changed nucleotides
        let mut changes: Vec<(usize, u8, u8)> = Vec::new();
        for i in 0..3 {
            if ref_codon_bytes[i] != best_codon.as_bytes()[i] {
                // c. position is 1-based from CDS start
                let c_pos = (aa_pos - 1) * 3 + i + 1;
                changes.push((c_pos, ref_codon_bytes[i], best_codon.as_bytes()[i]));
            }
        }

        let na_edit = if changes.len() == 1 {
            crate::edits::NaEdit::RefAlt {
                ref_: Some(String::from(changes[0].1 as char)),
                alt: Some(String::from(changes[0].2 as char)),
                uncertain: false,
            }
        } else {
            // Multiple nucleotide changes → delins
            let ref_nts: String = changes.iter().map(|(_, r, _)| *r as char).collect();
            let alt_nts: String = changes.iter().map(|(_, _, a)| *a as char).collect();
            crate::edits::NaEdit::RefAlt {
                ref_: Some(ref_nts),
                alt: Some(alt_nts),
                uncertain: false,
            }
        };

        let start_c_pos = changes.first().map(|(p, _, _)| *p).unwrap_or(1);
        let end_c_pos = changes.last().map(|(p, _, _)| *p).unwrap_or(start_c_pos);

        use crate::coords::{Anchor, HgvsTranscriptPos};
        use crate::structs::BaseOffsetPosition;

        let c_pos = BaseOffsetInterval {
            start: BaseOffsetPosition {
                base: HgvsTranscriptPos(start_c_pos as i32),
                offset: None,
                anchor: Anchor::CdsStart,
                uncertain: false,
            },
            end: if start_c_pos == end_c_pos {
                None
            } else {
                Some(BaseOffsetPosition {
                    base: HgvsTranscriptPos(end_c_pos as i32),
                    offset: None,
                    anchor: Anchor::CdsStart,
                    uncertain: false,
                })
            },
            uncertain: false,
        };

        let c_variant = CVariant {
            ac: tx_ac.clone(),
            gene: var_p.gene.clone(),
            posedit: crate::structs::PosEdit {
                pos: Some(c_pos),
                edit: na_edit,
                uncertain: false,
                predicted: false,
            },
        };

        Ok((c_variant, is_unique))
    }

    /// Normalizes a variant to its 3' most position.
    pub fn normalize_variant(
        &self,
        var: crate::SequenceVariant,
    ) -> Result<crate::SequenceVariant, HgvsError> {
        match var {
            crate::SequenceVariant::Coding(v_c) => Ok(crate::SequenceVariant::Coding(
                self.normalize_coding_variant(v_c)?,
            )),
            crate::SequenceVariant::Genomic(v_g) => Ok(crate::SequenceVariant::Genomic(
                self.normalize_genomic_variant(v_g)?,
            )),
            crate::SequenceVariant::NonCoding(v_n) => Ok(crate::SequenceVariant::NonCoding(
                self.normalize_noncoding_variant(v_n)?,
            )),
            _ => Ok(var),
        }
    }

    fn normalize_coding_variant(&self, mut v_c: CVariant) -> Result<CVariant, HgvsError> {
        let transcript = self.hdp.get_transcript(&v_c.ac, None)?;
        if let Some(pos) = &mut v_c.posedit.pos {
            let (start_idx, end_idx) = self.get_c_indices(pos, &transcript)?;
            let is_ins = matches!(&v_c.posedit.edit, crate::edits::NaEdit::Ins { .. });
            let (ins_anchor, actual_end) = ins_anchor_and_end(start_idx, end_idx, is_ins);
            let (new_start, _new_end) = self.shift_3_prime(
                &v_c.ac,
                IdentifierKind::Transcript,
                ins_anchor,
                actual_end,
                &v_c.posedit.edit,
            )?;

            if new_start != ins_anchor {
                // Re-derive HGVS positions via n_to_c so the c.0 gap is handled
                // correctly (shifting past the 5'UTR/CDS boundary must not produce
                // the invalid HgvsTranscriptPos(0)).
                let am = TranscriptMapper::new(dyn_clone::clone_box(&*transcript))?;
                if is_ins {
                    // Insertion anchor is the "after" base. HGVS spans [before, after].
                    pos.start = n_to_c_position(&am, (new_start as i32) - 1)?;
                    if pos.end.is_some() {
                        pos.end = Some(n_to_c_position(&am, new_start as i32)?);
                    }
                } else {
                    let width = end_idx - start_idx;
                    pos.start = n_to_c_position(&am, new_start as i32)?;
                    if let Some(e) = &mut pos.end {
                        *e = n_to_c_position(&am, (new_start + width - 1) as i32)?;
                    }
                }
            }

            self.update_del_dup_ref(
                &mut v_c.posedit.edit,
                &v_c.ac,
                IdentifierKind::Transcript,
                new_start,
                end_idx - start_idx,
            )?;

            // After 3'-shifting, convert Ins → Dup when the inserted sequence
            // exactly matches the reference at the insertion point.
            if is_ins {
                if let crate::edits::NaEdit::Ins {
                    alt: Some(ins_seq),
                    uncertain,
                } = v_c.posedit.edit.clone()
                {
                    let pos_after = v_c.posedit.pos.as_ref().unwrap();
                    // Only attempt for exonic positions (no intronic offsets).
                    if pos_after.start.offset.is_none()
                        && pos_after.end.as_ref().map_or(true, |e| e.offset.is_none())
                    {
                        if let Ok((cur_n_start, _)) = self.get_c_indices(pos_after, &transcript) {
                            let n = ins_seq.len() as i32;
                            let check_start = cur_n_start as i32 - n + 1;
                            if check_start >= 0 {
                                if let Ok(ref_seq) = self.hdp.get_seq(
                                    &v_c.ac,
                                    check_start,
                                    cur_n_start as i32 + 1,
                                    IdentifierKind::Transcript.into_identifier_type(),
                                ) {
                                    if ref_seq == ins_seq {
                                        let am2 = TranscriptMapper::new(dyn_clone::clone_box(
                                            &*transcript,
                                        ))?;
                                        if let Some(pos_mut) = &mut v_c.posedit.pos {
                                            pos_mut.start = n_to_c_position(&am2, check_start)?;
                                            pos_mut.end = if check_start != cur_n_start as i32 {
                                                Some(n_to_c_position(&am2, cur_n_start as i32)?)
                                            } else {
                                                None
                                            };
                                        }
                                        v_c.posedit.edit = crate::edits::NaEdit::Dup {
                                            ref_: Some(ins_seq),
                                            uncertain,
                                        };
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        Ok(v_c)
    }

    fn normalize_genomic_variant(&self, mut v_g: GVariant) -> Result<GVariant, HgvsError> {
        if let Some(pos) = &mut v_g.posedit.pos {
            let start_i = pos.start.base.to_index().0;
            if start_i < 0 {
                return Err(HgvsError::ValidationError(format!(
                    "Genomic start position {} is negative or uncertain",
                    start_i
                )));
            }
            let mut start_idx = checked_usize(start_i, "genomic start index")?;
            let is_ins = matches!(&v_g.posedit.edit, crate::edits::NaEdit::Ins { .. });
            let end_idx = if let Some(e) = &pos.end {
                let end_i = e.base.to_index().0;
                if end_i < 0 {
                    return Err(HgvsError::ValidationError(format!(
                        "Genomic end position {} is negative or uncertain",
                        end_i
                    )));
                }
                let idx = checked_usize(end_i, "genomic end index")?;
                if is_ins {
                    start_idx = idx;
                    idx
                } else {
                    idx.checked_add(1).ok_or_else(|| {
                        HgvsError::ValidationError("Genomic end position overflow".into())
                    })?
                }
            } else {
                start_idx.checked_add(1).ok_or_else(|| {
                    HgvsError::ValidationError("Genomic start position overflow".into())
                })?
            };

            let (new_start, new_end) = self.shift_3_prime(
                &v_g.ac,
                IdentifierKind::Genomic,
                start_idx,
                end_idx,
                &v_g.posedit.edit,
            )?;
            if new_start != start_idx {
                let shift = (new_start as i32) - (start_idx as i32);
                pos.start.base.0 += shift;
                if let Some(e) = &mut pos.end {
                    e.base.0 += shift;
                }
            }

            self.update_del_dup_ref(
                &mut v_g.posedit.edit,
                &v_g.ac,
                IdentifierKind::Genomic,
                new_start,
                new_end - new_start,
            )?;
        }
        Ok(v_g)
    }

    fn normalize_noncoding_variant(&self, mut v_n: NVariant) -> Result<NVariant, HgvsError> {
        let transcript = self.hdp.get_transcript(&v_n.ac, None)?;
        if let Some(pos) = &mut v_n.posedit.pos {
            let (start_idx, end_idx) = self.get_c_indices(pos, &transcript)?;
            let is_ins = matches!(&v_n.posedit.edit, crate::edits::NaEdit::Ins { .. });
            let (ins_anchor, actual_end) = ins_anchor_and_end(start_idx, end_idx, is_ins);
            let (new_start, new_end) = self.shift_3_prime(
                &v_n.ac,
                IdentifierKind::Transcript,
                ins_anchor,
                actual_end,
                &v_n.posedit.edit,
            )?;

            if new_start != ins_anchor {
                let shift = (new_start as i32) - (ins_anchor as i32);
                pos.start.base.0 += shift;
                if let Some(e) = &mut pos.end {
                    e.base.0 += shift;
                }
            }

            self.update_del_dup_ref(
                &mut v_n.posedit.edit,
                &v_n.ac,
                IdentifierKind::Transcript,
                new_start,
                new_end - new_start,
            )?;
        }
        Ok(v_n)
    }

    fn update_del_dup_ref(
        &self,
        edit: &mut crate::edits::NaEdit,
        ac: &str,
        kind: IdentifierKind,
        new_start: usize,
        span: usize,
    ) -> Result<(), HgvsError> {
        if let crate::edits::NaEdit::Del { ref_: r, .. }
        | crate::edits::NaEdit::Dup { ref_: r, .. } = edit
        {
            *r = Some(self.hdp.get_seq(
                ac,
                new_start as i32,
                (new_start + span) as i32,
                kind.into_identifier_type(),
            )?);
        }
        Ok(())
    }

    pub fn get_c_indices(
        &self,
        pos: &BaseOffsetInterval,
        transcript: &Box<dyn Transcript>,
    ) -> Result<(usize, usize), HgvsError> {
        let am = TranscriptMapper::new(dyn_clone::clone_box(&**transcript))?;
        let n_start = am.c_to_n(pos.start.base.to_index(), pos.start.anchor)?;
        let n_end = if let Some(e) = &pos.end {
            am.c_to_n(e.base.to_index(), e.anchor)?
        } else {
            n_start
        };
        if n_start.0 < 0 {
            return Err(HgvsError::ValidationError(format!(
                "Transcript start position {:?} maps to a negative index; \
                 normalization of 5'UTR-spanning variants is not yet supported",
                pos.start.base
            )));
        }
        let end_val = n_end.0.checked_add(1).filter(|&v| v >= 0).ok_or_else(|| {
            HgvsError::ValidationError(format!(
                "Transcript end position {:?} is out of sequence bounds",
                pos.end.as_ref().map(|e| e.base)
            ))
        })?;
        Ok((
            checked_usize(n_start.0, "transcript start index")?,
            checked_usize(end_val, "transcript end index")?,
        ))
    }

    fn shift_3_prime(
        &self,
        ac: &str,
        kind: IdentifierKind,
        start: usize,
        end: usize,
        edit: &crate::edits::NaEdit,
    ) -> Result<(usize, usize), HgvsError> {
        let (ref_owned, alt_owned) = match extract_edit_sequences(self.hdp, ac, start, end, edit)? {
            None => return Ok((start, end)),
            Some(seqs) => seqs,
        };
        let ref_str = ref_owned.as_str();
        let alt_str = alt_owned.as_str();

        if ref_str == alt_str && matches!(edit, crate::edits::NaEdit::RefAlt { .. }) {
            return Ok((start, end));
        }

        let mut curr_start = start;
        let mut curr_end = end;
        let mut chunk_size = 128;

        let mut chunk_start = end;
        let mut chunk = self.hdp.get_seq(
            ac,
            chunk_start as i32,
            (chunk_start + chunk_size) as i32,
            kind.into_identifier_type(),
        )?;
        let mut chunk_bytes = chunk.as_bytes();

        let is_del_or_dup = matches!(
            edit,
            crate::edits::NaEdit::Del { .. } | crate::edits::NaEdit::Dup { .. }
        );

        if is_del_or_dup
            || (!ref_str.is_empty() && alt_str.is_empty())
            || (matches!(edit, crate::edits::NaEdit::RefAlt { .. })
                && (end - start) != alt_str.len())
        {
            // Deletion, Duplication, or DelIns with a non-empty range
            let mut current_ref = if ref_str.is_empty() {
                self.hdp.get_seq(
                    ac,
                    curr_start as i32,
                    curr_end as i32,
                    kind.into_identifier_type(),
                )?
            } else {
                ref_str.to_string()
            };

            if current_ref.is_empty() {
                return Ok((curr_start, curr_end));
            }

            loop {
                if (curr_end - chunk_start) >= chunk_bytes.len() {
                    if chunk_bytes.len() < chunk_size {
                        break;
                    }
                    chunk_start += chunk_bytes.len();
                    chunk_size = std::cmp::min(chunk_size * 2, 4096);
                    chunk = self.hdp.get_seq(
                        ac,
                        chunk_start as i32,
                        (chunk_start + chunk_size) as i32,
                        kind.into_identifier_type(),
                    )?;
                    chunk_bytes = chunk.as_bytes();
                    if chunk_bytes.is_empty() {
                        break;
                    }
                }

                // To shift a delins/del/dup, the next base must match the first base of the range being shifted.
                // And the range must be "internally" repetitive or we must match the whole range?
                // Standard 3' shift: if seq[start] == seq[end], then [start, end) -> [start+1, end+1) is equivalent.
                let first_ref_byte = current_ref.as_bytes()[0];
                if first_ref_byte == chunk_bytes[curr_end - chunk_start] {
                    curr_start += 1;
                    curr_end += 1;
                    // Update current_ref for the next iteration (it's the sequence at the new [start, end))
                    current_ref = self.hdp.get_seq(
                        ac,
                        curr_start as i32,
                        curr_end as i32,
                        kind.into_identifier_type(),
                    )?;
                    if current_ref.is_empty() {
                        break;
                    }
                } else {
                    break;
                }
            }
        } else if start == end && !alt_str.is_empty() {
            // Pure Insertion (start == end)
            let alt_bytes = alt_str.as_bytes();
            let n = alt_bytes.len();
            if n == 0 {
                return Ok((curr_start, curr_end));
            }

            loop {
                if (curr_end - chunk_start) >= chunk_bytes.len() {
                    if chunk_bytes.len() < chunk_size {
                        break;
                    }
                    chunk_start += chunk_bytes.len();
                    chunk_size = std::cmp::min(chunk_size * 2, 4096);
                    chunk = self.hdp.get_seq(
                        ac,
                        chunk_start as i32,
                        (chunk_start + chunk_size) as i32,
                        kind.into_identifier_type(),
                    )?;
                    chunk_bytes = chunk.as_bytes();
                    if chunk_bytes.is_empty() {
                        break;
                    }
                }

                // For an insertion to shift right, the base we pass MUST match the base we are putting "behind" it.
                // If we insert 'ABC' at pos 1 in 'XABC', we can move it to pos 2 only if ref[1] == 'A'.
                // 'X | ABC' -> 'XA | BCA' -> 'XAB | CAB' -> 'XABC | ABC'.
                // So at each step k, we need ref[end+k] == alt[k % n].
                if chunk_bytes[curr_end - chunk_start] == alt_bytes[(curr_start - start) % n] {
                    curr_start += 1;
                    curr_end += 1;
                } else {
                    break;
                }
            }
        }
        Ok((curr_start, curr_end))
    }

    fn shift_5_prime(
        &self,
        ac: &str,
        kind: IdentifierKind,
        start: usize,
        end: usize,
        edit: &crate::edits::NaEdit,
    ) -> Result<(usize, usize), HgvsError> {
        let (ref_owned, alt_owned) = match extract_edit_sequences(self.hdp, ac, start, end, edit)? {
            None => return Ok((start, end)),
            Some(seqs) => seqs,
        };
        let ref_str = ref_owned.as_str();
        let alt_str = alt_owned.as_str();

        if ref_str == alt_str && matches!(edit, crate::edits::NaEdit::RefAlt { .. }) {
            return Ok((start, end));
        }

        let mut curr_start = start;
        let mut curr_end = end;

        let is_del_or_dup = matches!(
            edit,
            crate::edits::NaEdit::Del { .. } | crate::edits::NaEdit::Dup { .. }
        );

        if is_del_or_dup
            || (!ref_str.is_empty() && alt_str.is_empty())
            || (matches!(edit, crate::edits::NaEdit::RefAlt { .. })
                && (end - start) != alt_str.len())
        {
            // Deletion, Duplication, or DelIns with a non-empty range
            let mut current_ref = if ref_str.is_empty() {
                self.hdp.get_seq(
                    ac,
                    curr_start as i32,
                    curr_end as i32,
                    kind.into_identifier_type(),
                )?
            } else {
                ref_str.to_string()
            };

            if current_ref.is_empty() {
                return Ok((curr_start, curr_end));
            }

            loop {
                if curr_start == 0 {
                    break;
                }

                let prev_base_pos = curr_start - 1;
                let prev_base = self.hdp.get_seq(
                    ac,
                    prev_base_pos as i32,
                    curr_start as i32,
                    kind.into_identifier_type(),
                )?;
                if prev_base.is_empty() {
                    break;
                }

                let last_ref_byte = current_ref.as_bytes()[current_ref.len() - 1];
                if prev_base.as_bytes()[0] == last_ref_byte {
                    curr_start -= 1;
                    curr_end -= 1;
                    current_ref = self.hdp.get_seq(
                        ac,
                        curr_start as i32,
                        curr_end as i32,
                        kind.into_identifier_type(),
                    )?;
                } else {
                    break;
                }
            }
        } else if start == end && !alt_str.is_empty() {
            // Pure Insertion
            let alt_bytes = alt_str.as_bytes();
            let n = alt_bytes.len();

            loop {
                if curr_start == 0 {
                    break;
                }
                let prev_base_pos = curr_start - 1;
                let prev_base = self.hdp.get_seq(
                    ac,
                    prev_base_pos as i32,
                    curr_start as i32,
                    kind.into_identifier_type(),
                )?;
                if prev_base.is_empty() {
                    break;
                }

                let rel_pos =
                    (curr_start as i64 - start as i64 + n as i64 - 1).rem_euclid(n as i64) as usize;
                if prev_base.as_bytes()[0] == alt_bytes[rel_pos] {
                    curr_start -= 1;
                    curr_end -= 1;
                } else {
                    break;
                }
            }
        }
        Ok((curr_start, curr_end))
    }

    pub fn expand_unambiguous_range(
        &self,
        ac: &str,
        kind: IdentifierKind,
        start: usize,
        end: usize,
        edit: &crate::edits::NaEdit,
    ) -> Result<(usize, usize), HgvsError> {
        // Substitutions in homopolymers are NOT expanded in ClinVar/SPDI standard.
        // We only expand length-changing variants (Del, Ins, Dup, Repeat).
        let is_length_changing = match edit {
            crate::edits::NaEdit::RefAlt { alt, .. } => {
                let r_len = end - start;
                let a_len = alt.as_deref().unwrap_or("").len();
                r_len != a_len
            }
            crate::edits::NaEdit::Del { .. }
            | crate::edits::NaEdit::Ins { .. }
            | crate::edits::NaEdit::Dup { .. }
            | crate::edits::NaEdit::Repeat { .. } => true,
            _ => false,
        };

        if !is_length_changing {
            return Ok((start, end));
        }

        let (s_5, _) = self.shift_5_prime(ac, kind, start, end, edit)?;
        let (_, e_3) = self.shift_3_prime(ac, kind, start, end, edit)?;
        Ok((s_5, e_3))
    }

    pub fn to_spdi(
        &self,
        var: &crate::SequenceVariant,
        unambiguous: bool,
    ) -> Result<String, HgvsError> {
        if unambiguous {
            self.to_spdi_unambiguous(var)
        } else {
            // 1. Resolve to genomic if possible.
            let g_var_obj = match var {
                crate::SequenceVariant::Genomic(v) => v.clone(),
                crate::SequenceVariant::Coding(v) => self.c_to_g(v, None)?,
                crate::SequenceVariant::NonCoding(v) => self.n_to_g(v, None)?,
                _ => {
                    return Err(HgvsError::UnsupportedOperation(
                        "SPDI only for genomic/coding/non-coding".into(),
                    ))
                }
            };

            // 2. Normalize (3' shift, minimal delins)
            let g_norm_var = self.normalize_variant(crate::SequenceVariant::Genomic(g_var_obj))?;
            let g_norm = match g_norm_var {
                crate::SequenceVariant::Genomic(v) => v,
                _ => unreachable!(),
            };
            g_norm.posedit.to_spdi(&g_norm.ac, &*self.hdp)
        }
    }

    pub fn to_spdi_unambiguous(&self, var: &crate::SequenceVariant) -> Result<String, HgvsError> {
        // 1. Resolve to genomic if possible. Unambiguous SPDI is ideally on chromosomal coordinates.
        let g_var_obj = match var {
            crate::SequenceVariant::Genomic(v) => v.clone(),
            crate::SequenceVariant::Coding(v) => self.c_to_g(v, None)?,
            crate::SequenceVariant::NonCoding(v) => self.n_to_g(v, None)?,
            _ => {
                return Err(HgvsError::UnsupportedOperation(
                    "SPDI expansion only for genomic/coding/non-coding".into(),
                ))
            }
        };

        // 2. Normalize (3' shift, minimal delins)
        let g_norm_var = self.normalize_variant(crate::SequenceVariant::Genomic(g_var_obj))?;
        let g_norm = match g_norm_var {
            crate::SequenceVariant::Genomic(v) => v,
            _ => unreachable!(),
        };

        let ac = &g_norm.ac;
        if let Some(pos) = &g_norm.posedit.pos {
            let start_i = pos.start.base.to_index().0;
            if start_i < 0 {
                return Err(HgvsError::ValidationError(format!(
                    "Genomic start position {} is negative or uncertain",
                    start_i
                )));
            }
            let mut start_idx = checked_usize(start_i, "genomic start index")?;
            let is_ins = matches!(&g_norm.posedit.edit, crate::edits::NaEdit::Ins { .. });
            let end_idx = if let Some(e) = &pos.end {
                let end_i = e.base.to_index().0;
                if end_i < 0 {
                    return Err(HgvsError::ValidationError(format!(
                        "Genomic end position {} is negative or uncertain",
                        end_i
                    )));
                }
                let idx = checked_usize(end_i, "genomic end index")?;

                if is_ins {
                    // Mirror normalize_variant: use the end position as the insertion
                    // anchor so that expand_unambiguous_range sees start == end and
                    // correctly expands the ambiguous run.
                    start_idx = idx;
                    idx
                } else {
                    idx.checked_add(1).ok_or_else(|| {
                        HgvsError::ValidationError("Genomic end position overflow".into())
                    })?
                }
            } else {
                start_idx.checked_add(1).ok_or_else(|| {
                    HgvsError::ValidationError("Genomic start position overflow".into())
                })?
            };

            // 3. Expand range to cover ambiguity
            let (u_start, u_end) = self.expand_unambiguous_range(
                ac,
                IdentifierKind::Genomic,
                start_idx,
                end_idx,
                &g_norm.posedit.edit,
            )?;

            // 4. Construct expanded sequences
            let r_seq = self.hdp.get_seq(
                ac,
                u_start as i32,
                u_end as i32,
                IdentifierType::GenomicAccession,
            )?;

            let rel_start = start_idx - u_start;
            let rel_end = end_idx - u_start;

            let alt_storage;
            let alt_str = match &g_norm.posedit.edit {
                crate::edits::NaEdit::RefAlt { alt, .. } => alt.as_deref().unwrap_or(""),
                crate::edits::NaEdit::Ins { alt: Some(s), .. } => s.as_str(),
                crate::edits::NaEdit::Del { .. } => "",
                crate::edits::NaEdit::Dup { ref_: Some(s), .. } => {
                    alt_storage = format!("{}{}", s, s);
                    &alt_storage
                }
                crate::edits::NaEdit::Repeat { ref_, max, .. } => {
                    let unit = if let Some(u) = ref_ {
                        u.clone()
                    } else {
                        self.hdp.get_seq(
                            ac,
                            start_idx as i32,
                            end_idx as i32,
                            IdentifierType::GenomicAccession,
                        )?
                    };
                    alt_storage = unit.repeat(*max as usize);
                    &alt_storage
                }
                crate::edits::NaEdit::Inv { .. } => {
                    let s = self.hdp.get_seq(
                        ac,
                        start_idx as i32,
                        end_idx as i32,
                        IdentifierType::GenomicAccession,
                    )?;
                    alt_storage = crate::sequence::rev_comp(&s);
                    &alt_storage
                }
                _ => return g_norm.posedit.to_spdi(ac, &*self.hdp),
            };

            let a_seq = format!("{}{}{}", &r_seq[..rel_start], alt_str, &r_seq[rel_end..]);

            Ok(format!("{}:{}:{}:{}", ac, u_start, r_seq, a_seq))
        } else {
            g_norm.posedit.to_spdi(&g_norm.ac, &*self.hdp) // Fallback for identity?
        }
    }
}
