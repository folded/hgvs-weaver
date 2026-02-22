use crate::altseq::AltSeqBuilder;
use crate::altseq_to_hgvsp::AltSeqToHgvsp;
use crate::data::{DataProvider, IdentifierKind, IdentifierType, Transcript, TranscriptSearch};
use crate::error::HgvsError;
use crate::sequence::{LazySequence, MemSequence, RevCompSequence, Sequence, TranslatedSequence};
use crate::structs::{BaseOffsetInterval, CVariant, GVariant, PVariant};
use crate::transcript_mapper::TranscriptMapper;

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
        let (n_pos, offset) = am.g_to_n(g_start_0)?;
        let (c_pos_index, c_offset, anchor) = am.n_to_c(n_pos)?;

        let pos_c = crate::structs::BaseOffsetPosition {
            base: c_pos_index.to_hgvs(),
            offset: if c_offset.0 + offset.0 != 0 {
                Some(crate::structs::IntronicOffset(c_offset.0 + offset.0))
            } else {
                None
            },
            anchor,
            uncertain: false,
        };

        if let Some(end_g_simple) = &pos.end {
            let g_end_0 = end_g_simple.base.to_index();
            let (n_pos_e, offset_e) = am.g_to_n(g_end_0)?;
            let (c_pos_e_index, c_offset_e, anchor_e) = am.n_to_c(n_pos_e)?;
            let pos_c_e = crate::structs::BaseOffsetPosition {
                base: c_pos_e_index.to_hgvs(),
                offset: if c_offset_e.0 + offset_e.0 != 0 {
                    Some(crate::structs::IntronicOffset(c_offset_e.0 + offset_e.0))
                } else {
                    None
                },
                anchor: anchor_e,
                uncertain: false,
            };
            let mut edit = var_g.posedit.edit.clone();
            if am.transcript.strand() == -1 {
                edit = edit.map_sequence(|s| {
                    let seq = crate::sequence::MemSequence(s.to_string());
                    let rc = crate::sequence::RevCompSequence { inner: &seq };
                    rc.to_string()
                });
            }

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

        let mut edit = var_g.posedit.edit.clone();
        if am.transcript.strand() == -1 {
            edit = edit.map_sequence(|s| {
                let seq = crate::sequence::MemSequence(s.to_string());
                let rc = crate::sequence::RevCompSequence { inner: &seq };
                rc.to_string()
            });
        }

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
            let mut pos_g = crate::structs::SimplePosition {
                base: g_pos.to_hgvs(),
                end: None,
                uncertain: false,
            };
            let mut pos_g_e = crate::structs::SimplePosition {
                base: g_pos_e.to_hgvs(),
                end: None,
                uncertain: false,
            };

            if pos_g.base.0 > pos_g_e.base.0 {
                std::mem::swap(&mut pos_g, &mut pos_g_e);
            }

            let mut edit = var_c.posedit.edit.clone();
            if am.transcript.strand() == -1 {
                edit = edit.map_sequence(|s| {
                    let seq = MemSequence(s.to_string());
                    let rc = RevCompSequence { inner: &seq };
                    rc.to_string()
                });
            }

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

        let pos_g = crate::structs::SimplePosition {
            base: g_pos.to_hgvs(),
            end: None,
            uncertain: false,
        };
        let mut edit = var_c.posedit.edit.clone();
        if am.transcript.strand() == -1 {
            edit = edit.map_sequence(|s| {
                let seq = crate::sequence::MemSequence(s.to_string());
                let rc = crate::sequence::RevCompSequence { inner: &seq };
                rc.to_string()
            });
        }

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
            let mut pos_g = crate::structs::SimplePosition {
                base: g_pos.to_hgvs(),
                end: None,
                uncertain: false,
            };
            let mut pos_g_e = crate::structs::SimplePosition {
                base: g_pos_e.to_hgvs(),
                end: None,
                uncertain: false,
            };

            if pos_g.base.0 > pos_g_e.base.0 {
                std::mem::swap(&mut pos_g, &mut pos_g_e);
            }

            let mut edit = var_n.posedit.edit.clone();
            if am.transcript.strand() == -1 {
                edit = edit.map_sequence(|s| {
                    let seq = MemSequence(s.to_string());
                    let rc = RevCompSequence { inner: &seq };
                    rc.to_string()
                });
            }

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

        let pos_g = crate::structs::SimplePosition {
            base: g_pos.to_hgvs(),
            end: None,
            uncertain: false,
        };
        let mut edit = var_n.posedit.edit.clone();
        if am.transcript.strand() == -1 {
            edit = edit.map_sequence(|s| {
                let seq = crate::sequence::MemSequence(s.to_string());
                let rc = crate::sequence::RevCompSequence { inner: &seq };
                rc.to_string()
            });
        }

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
        let mut results = Vec::new();

        for tx_ac in transcripts {
            if let Ok(vc) = self.g_to_c(var_g, &tx_ac) {
                results.push(vc);
            }
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

        let transcript = self.hdp.get_transcript(transcript_ac, None)?;
        let ref_seq = self.hdp.get_seq(
            transcript_ac,
            0,
            -1,
            IdentifierKind::Transcript.into_identifier_type(),
        )?;

        let cds_start_idx = transcript
            .cds_start_index()
            .ok_or_else(|| HgvsError::ValidationError("Missing CDS start".into()))?
            .0 as usize;
        let cds_end_idx = transcript
            .cds_end_index()
            .ok_or_else(|| HgvsError::ValidationError("Missing CDS end".into()))?
            .0 as usize;

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

        // Use Sequence abstraction for translation
        let trans_obj = TranslatedSequence {
            inner: &LazySequence {
                hdp: self.hdp,
                ac: transcript_ac.to_string(),
                start: cds_start_idx,
                end: ref_seq_obj.len(),
                kind: IdentifierType::TranscriptAccession,
            },
        };
        let ref_aa = trans_obj.to_string();

        let builder = AltSeqBuilder {
            var_c,
            transcript_sequence: &ref_seq_obj,
            cds_start_index: transcript.cds_start_index().unwrap(),
            cds_end_index: transcript.cds_end_index().unwrap(),
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

    /// Normalizes a variant to its 3' most position.
    pub fn normalize_variant(
        &self,
        var: crate::SequenceVariant,
    ) -> Result<crate::SequenceVariant, HgvsError> {
        match var {
            crate::SequenceVariant::Coding(mut v_c) => {
                let transcript = self.hdp.get_transcript(&v_c.ac, None)?;
                if let Some(pos) = &mut v_c.posedit.pos {
                    let (start_idx, end_idx) = self.get_c_indices(pos, &transcript)?;
                    let is_ins = matches!(&v_c.posedit.edit, crate::edits::NaEdit::Ins { .. });
                    let actual_end = if is_ins { end_idx - 1 } else { end_idx };

                    let (new_start, _new_end) = self.shift_3_prime(
                        &v_c.ac,
                        IdentifierKind::Transcript,
                        start_idx,
                        actual_end,
                        &v_c.posedit.edit,
                    )?;

                    if new_start != start_idx {
                        let shift = (new_start as i32) - (start_idx as i32);
                        pos.start.base.0 += shift;
                        if let Some(e) = &mut pos.end {
                            e.base.0 += shift;
                        }
                    }

                    // Update sequences for Del/Dup to match reference at new position
                    match &mut v_c.posedit.edit {
                        crate::edits::NaEdit::Del { ref_: r, .. }
                        | crate::edits::NaEdit::Dup { ref_: r, .. } => {
                            let seq = self.hdp.get_seq(
                                &v_c.ac,
                                new_start as i32,
                                (new_start + (end_idx - start_idx)) as i32,
                                IdentifierKind::Transcript.into_identifier_type(),
                            )?;
                            *r = Some(seq);
                        }
                        _ => {}
                    }
                }
                Ok(crate::SequenceVariant::Coding(v_c))
            }
            crate::SequenceVariant::Genomic(mut v_g) => {
                if let Some(pos) = &mut v_g.posedit.pos {
                    let mut start_idx = pos.start.base.to_index().0 as usize;
                    let is_ins = matches!(&v_g.posedit.edit, crate::edits::NaEdit::Ins { .. });
                    let end_idx = pos.end.as_ref().map_or(start_idx + 1, |e| {
                        let idx = e.base.to_index().0 as usize;
                        if is_ins {
                            start_idx = idx;
                            idx
                        } else {
                            idx + 1
                        }
                    });

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

                    // Update sequences for Del/Dup to match reference at new position
                    match &mut v_g.posedit.edit {
                        crate::edits::NaEdit::Del { ref_: r, .. }
                        | crate::edits::NaEdit::Dup { ref_: r, .. } => {
                            let seq = self.hdp.get_seq(
                                &v_g.ac,
                                new_start as i32,
                                new_end as i32,
                                IdentifierKind::Genomic.into_identifier_type(),
                            )?;
                            *r = Some(seq);
                        }
                        _ => {}
                    }
                }
                Ok(crate::SequenceVariant::Genomic(v_g))
            }
            crate::SequenceVariant::NonCoding(mut v_n) => {
                let transcript = self.hdp.get_transcript(&v_n.ac, None)?;
                if let Some(pos) = &mut v_n.posedit.pos {
                    let (start_idx, end_idx) = self.get_n_indices(pos, &transcript)?;
                    let is_ins = matches!(&v_n.posedit.edit, crate::edits::NaEdit::Ins { .. });
                    let actual_end = if is_ins { end_idx - 1 } else { end_idx };

                    let (new_start, new_end) = self.shift_3_prime(
                        &v_n.ac,
                        IdentifierKind::Transcript,
                        start_idx,
                        actual_end,
                        &v_n.posedit.edit,
                    )?;

                    if new_start != start_idx {
                        let shift = (new_start as i32) - (start_idx as i32);
                        pos.start.base.0 += shift;
                        if let Some(e) = &mut pos.end {
                            e.base.0 += shift;
                        }
                    }

                    // Update sequences for Del/Dup to match reference at new position
                    match &mut v_n.posedit.edit {
                        crate::edits::NaEdit::Del { ref_: r, .. }
                        | crate::edits::NaEdit::Dup { ref_: r, .. } => {
                            let seq = self.hdp.get_seq(
                                &v_n.ac,
                                new_start as i32,
                                new_end as i32,
                                IdentifierKind::Transcript.into_identifier_type(),
                            )?;
                            *r = Some(seq);
                        }
                        _ => {}
                    }
                }
                Ok(crate::SequenceVariant::NonCoding(v_n))
            }
            _ => Ok(var),
        }
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
        Ok((n_start.0 as usize, (n_end.0 + 1) as usize))
    }

    pub fn get_n_indices(
        &self,
        pos: &BaseOffsetInterval,
        transcript: &Box<dyn Transcript>,
    ) -> Result<(usize, usize), HgvsError> {
        // For n. variants, we don't need CDS.
        // We use TranscriptMapper but we should be careful about anchors.
        // Actually TranscriptMapper::new might fail if CDS is missing?
        // Let's check TranscriptMapper::new.
        let am = TranscriptMapper::new(dyn_clone::clone_box(&**transcript))?;
        let n_start = am.c_to_n(pos.start.base.to_index(), pos.start.anchor)?;
        let n_end = if let Some(e) = &pos.end {
            am.c_to_n(e.base.to_index(), e.anchor)?
        } else {
            n_start
        };
        Ok((n_start.0 as usize, (n_end.0 + 1) as usize))
    }

    fn shift_3_prime(
        &self,
        ac: &str,
        kind: IdentifierKind,
        start: usize,
        end: usize,
        edit: &crate::edits::NaEdit,
    ) -> Result<(usize, usize), HgvsError> {
        let storage_r;
        let storage_a;
        let (ref_str, alt_str) = match edit {
            crate::edits::NaEdit::RefAlt { ref_, alt, .. } => {
                (ref_.as_deref().unwrap_or(""), alt.as_deref().unwrap_or(""))
            }
            crate::edits::NaEdit::Del { ref_: Some(s), .. } => (s.as_str(), ""),
            crate::edits::NaEdit::Del { ref_: None, .. } => ("", ""),
            crate::edits::NaEdit::Ins { alt: Some(s), .. } => ("", s.as_str()),
            crate::edits::NaEdit::Ins { alt: None, .. } => ("", ""),
            crate::edits::NaEdit::Dup { ref_: Some(s), .. } => (s.as_str(), ""),
            crate::edits::NaEdit::Dup { ref_: None, .. } => ("", ""),
            crate::edits::NaEdit::Repeat { ref_, max, .. } => {
                storage_r = if let Some(r) = ref_ {
                    r.clone()
                } else {
                    self.hdp.get_seq(
                        ac,
                        start as i32,
                        end as i32,
                        IdentifierType::GenomicAccession,
                    )?
                };
                storage_a = storage_r.repeat(*max as usize);
                (storage_r.as_str(), storage_a.as_str())
            }
            crate::edits::NaEdit::Inv { .. } => {
                storage_r = self.hdp.get_seq(
                    ac,
                    start as i32,
                    end as i32,
                    IdentifierType::GenomicAccession,
                )?;
                storage_a = crate::sequence::rev_comp(&storage_r);
                (storage_r.as_str(), storage_a.as_str())
            }
            _ => return Ok((start, end)),
        };

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
                // If we insert 'A' at pos 1 in 'XA', we can move it to pos 2 only if ref[1] == 'A'.
                // So 'X | A' -> 'X A |'. Both result in 'XAA'.
                if chunk_bytes[curr_end - chunk_start] == alt_bytes[0] {
                    // For 1-base insertions, shifting is simple.
                    // For multi-base, we'd need to "rotate" the alt string (TODO).
                    if n == 1 {
                        curr_start += 1;
                        curr_end += 1;
                    } else {
                        break;
                    }
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
        let storage_r;
        let storage_a;
        let (ref_str, alt_str) = match edit {
            crate::edits::NaEdit::RefAlt { ref_, alt, .. } => {
                (ref_.as_deref().unwrap_or(""), alt.as_deref().unwrap_or(""))
            }
            crate::edits::NaEdit::Del { ref_: Some(s), .. } => (s.as_str(), ""),
            crate::edits::NaEdit::Del { ref_: None, .. } => ("", ""),
            crate::edits::NaEdit::Ins { alt: Some(s), .. } => ("", s.as_str()),
            crate::edits::NaEdit::Ins { alt: None, .. } => ("", ""),
            crate::edits::NaEdit::Dup { ref_: Some(s), .. } => (s.as_str(), ""),
            crate::edits::NaEdit::Dup { ref_: None, .. } => ("", ""),
            crate::edits::NaEdit::Repeat { ref_, max, .. } => {
                storage_r = if let Some(r) = ref_ {
                    r.clone()
                } else {
                    self.hdp.get_seq(
                        ac,
                        start as i32,
                        end as i32,
                        IdentifierType::GenomicAccession,
                    )?
                };
                storage_a = storage_r.repeat(*max as usize);
                (storage_r.as_str(), storage_a.as_str())
            }
            crate::edits::NaEdit::Inv { .. } => {
                storage_r = self.hdp.get_seq(
                    ac,
                    start as i32,
                    end as i32,
                    IdentifierType::GenomicAccession,
                )?;
                storage_a = crate::sequence::rev_comp(&storage_r);
                (storage_r.as_str(), storage_a.as_str())
            }
            _ => return Ok((start, end)),
        };

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

                if prev_base.as_bytes()[0] == alt_bytes[n - 1] {
                    if n == 1 {
                        curr_start -= 1;
                        curr_end -= 1;
                    } else {
                        break;
                    }
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
            let start_idx = pos.start.base.to_index().0 as usize;
            let is_ins = matches!(&g_norm.posedit.edit, crate::edits::NaEdit::Ins { .. });
            let end_idx = pos.end.as_ref().map_or(start_idx + 1, |e| {
                let idx = e.base.to_index().0 as usize;
                if is_ins {
                    idx
                } else {
                    idx + 1
                }
            });

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
