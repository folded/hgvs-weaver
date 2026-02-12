use crate::error::HgvsError;
use crate::data::{DataProvider, TranscriptSearch, IdentifierKind, Transcript};
use crate::structs::{GVariant, CVariant, PVariant, BaseOffsetInterval};
use crate::transcript_mapper::TranscriptMapper;
use crate::altseq::AltSeqBuilder;
use crate::altseq_to_hgvsp::AltSeqToHgvsp;

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

        let pos = var_g.posedit.pos.as_ref().ok_or_else(|| HgvsError::ValidationError("Missing genomic position".into()))?;
        let g_start_0 = pos.start.base.to_index();
        let (n_pos, offset) = am.g_to_n(g_start_0)?;
        let (c_pos_index, c_offset, anchor) = am.n_to_c(n_pos)?;

        let pos_c = crate::structs::BaseOffsetPosition {
            base: c_pos_index.to_hgvs(),
            offset: if c_offset.0 + offset.0 != 0 { Some(crate::structs::IntronicOffset(c_offset.0 + offset.0)) } else { None },
            anchor,
            uncertain: false,
        };

        if let Some(end_g_simple) = &pos.end {
            let g_end_0 = end_g_simple.base.to_index();
            let (n_pos_e, offset_e) = am.g_to_n(g_end_0)?;
            let (c_pos_e_index, c_offset_e, anchor_e) = am.n_to_c(n_pos_e)?;
            let pos_c_e = crate::structs::BaseOffsetPosition {
                base: c_pos_e_index.to_hgvs(),
                offset: if c_offset_e.0 + offset_e.0 != 0 { Some(crate::structs::IntronicOffset(c_offset_e.0 + offset_e.0)) } else { None },
                anchor: anchor_e,
                uncertain: false,
            };
            let mut edit = var_g.posedit.edit.clone();
            if am.transcript.strand() == -1 { edit = edit.reverse_complement(); }

            return Ok(CVariant {
                ac: transcript_ac.to_string(),
                gene: var_g.gene.clone(),
                posedit: crate::structs::PosEdit {
                    pos: Some(crate::structs::BaseOffsetInterval { start: pos_c, end: Some(pos_c_e), uncertain: false }),
                    edit,
                    uncertain: var_g.posedit.uncertain,
                    predicted: var_g.posedit.predicted,
                }
            });
        }

        let mut edit = var_g.posedit.edit.clone();
        if am.transcript.strand() == -1 { edit = edit.reverse_complement(); }

        Ok(CVariant {
            ac: transcript_ac.to_string(),
            gene: var_g.gene.clone(),
            posedit: crate::structs::PosEdit {
                pos: Some(crate::structs::BaseOffsetInterval { start: pos_c, end: None, uncertain: false }),
                edit,
                uncertain: var_g.posedit.uncertain,
                predicted: var_g.posedit.predicted,
            }
        })
    }

    /// Transforms a coding cDNA variant (`c.`) to a genomic variant (`g.`).
    pub fn c_to_g(&self, var_c: &CVariant, reference_ac: &str) -> Result<GVariant, HgvsError> {
        let transcript = self.hdp.get_transcript(&var_c.ac, Some(reference_ac))?;
        let am = TranscriptMapper::new(transcript)?;

        let pos = var_c.posedit.pos.as_ref().ok_or_else(|| HgvsError::ValidationError("Missing cDNA position".into()))?;
        let n_pos = am.c_to_n(pos.start.base.to_index(), pos.start.anchor)?;
        let g_pos = am.n_to_g(n_pos, pos.start.offset.unwrap_or(crate::structs::IntronicOffset(0)))?;

        if let Some(end_c) = &pos.end {
            let n_pos_e = am.c_to_n(end_c.base.to_index(), end_c.anchor)?;
            let g_pos_e = am.n_to_g(n_pos_e, end_c.offset.unwrap_or(crate::structs::IntronicOffset(0)))?;
            let mut pos_g = crate::structs::SimplePosition { base: g_pos.to_hgvs(), end: None, uncertain: false };
            let mut pos_g_e = crate::structs::SimplePosition { base: g_pos_e.to_hgvs(), end: None, uncertain: false };

            if pos_g.base.0 > pos_g_e.base.0 {
                std::mem::swap(&mut pos_g, &mut pos_g_e);
            }

            let mut edit = var_c.posedit.edit.clone();
            if am.transcript.strand() == -1 { edit = edit.reverse_complement(); }

            return Ok(GVariant {
                ac: reference_ac.to_string(),
                gene: var_c.gene.clone(),
                posedit: crate::structs::PosEdit {
                    pos: Some(crate::structs::SimpleInterval { start: pos_g, end: Some(pos_g_e), uncertain: false }),
                    edit,
                    uncertain: var_c.posedit.uncertain,
                    predicted: var_c.posedit.predicted,
                }
            });
        }

        let pos_g = crate::structs::SimplePosition { base: g_pos.to_hgvs(), end: None, uncertain: false };
        let mut edit = var_c.posedit.edit.clone();
        if am.transcript.strand() == -1 { edit = edit.reverse_complement(); }

        Ok(GVariant {
            ac: reference_ac.to_string(),
            gene: var_c.gene.clone(),
            posedit: crate::structs::PosEdit {
                pos: Some(crate::structs::SimpleInterval { start: pos_g, end: None, uncertain: false }),
                edit,
                uncertain: var_c.posedit.uncertain,
                predicted: var_c.posedit.predicted,
            }
        })
    }

    /// Transforms a non-coding cDNA variant (`n.`) to a genomic variant (`g.`).
    pub fn n_to_g(&self, var_n: &crate::structs::NVariant, reference_ac: &str) -> Result<GVariant, HgvsError> {
        let transcript = self.hdp.get_transcript(&var_n.ac, Some(reference_ac))?;
        let am = TranscriptMapper::new(transcript)?;

        let pos = var_n.posedit.pos.as_ref().ok_or_else(|| HgvsError::ValidationError("Missing cDNA position".into()))?;
        let n_pos = am.c_to_n(pos.start.base.to_index(), pos.start.anchor)?;
        let g_pos = am.n_to_g(n_pos, pos.start.offset.unwrap_or(crate::structs::IntronicOffset(0)))?;

        if let Some(end_n) = &pos.end {
            let n_pos_e = am.c_to_n(end_n.base.to_index(), end_n.anchor)?;
            let g_pos_e = am.n_to_g(n_pos_e, end_n.offset.unwrap_or(crate::structs::IntronicOffset(0)))?;
            let mut pos_g = crate::structs::SimplePosition { base: g_pos.to_hgvs(), end: None, uncertain: false };
            let mut pos_g_e = crate::structs::SimplePosition { base: g_pos_e.to_hgvs(), end: None, uncertain: false };

            if pos_g.base.0 > pos_g_e.base.0 {
                std::mem::swap(&mut pos_g, &mut pos_g_e);
            }

            let mut edit = var_n.posedit.edit.clone();
            if am.transcript.strand() == -1 { edit = edit.reverse_complement(); }

            return Ok(GVariant {
                ac: reference_ac.to_string(),
                gene: var_n.gene.clone(),
                posedit: crate::structs::PosEdit {
                    pos: Some(crate::structs::SimpleInterval { start: pos_g, end: Some(pos_g_e), uncertain: false }),
                    edit,
                    uncertain: var_n.posedit.uncertain,
                    predicted: var_n.posedit.predicted,
                }
            });
        }

        let pos_g = crate::structs::SimplePosition { base: g_pos.to_hgvs(), end: None, uncertain: false };
        let mut edit = var_n.posedit.edit.clone();
        if am.transcript.strand() == -1 { edit = edit.reverse_complement(); }

        Ok(GVariant {
            ac: reference_ac.to_string(),
            gene: var_n.gene.clone(),
            posedit: crate::structs::PosEdit {
                pos: Some(crate::structs::SimpleInterval { start: pos_g, end: None, uncertain: false }),
                edit,
                uncertain: var_n.posedit.uncertain,
                predicted: var_n.posedit.predicted,
            }
        })
    }

    /// Discovers all possible cDNA consequences for a genomic variant.
    pub fn g_to_c_all(&self, var_g: &GVariant, searcher: &dyn TranscriptSearch) -> Result<Vec<CVariant>, HgvsError> {
        let pos = var_g.posedit.pos.as_ref().ok_or_else(|| HgvsError::ValidationError("Missing position".into()))?;
        let start_0 = pos.start.base.to_index().0;
        let end_0 = pos.end.as_ref().map_or(start_0 + 1, |e| e.base.to_index().0 + 1);

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
    pub fn c_to_p(&self, var_c: &CVariant, protein_ac: Option<&str>) -> Result<PVariant, HgvsError> {
        let transcript_ac = &var_c.ac;
        let pro_ac_str = if let Some(ac) = protein_ac {
            ac.to_string()
        } else {
            self.hdp.get_symbol_accessions(transcript_ac, IdentifierKind::Transcript, IdentifierKind::Protein)?
                .first()
                .ok_or_else(|| HgvsError::ValidationError(format!("No protein accession found for {}", transcript_ac)))?
                .1.clone()
        };

        let transcript = self.hdp.get_transcript(transcript_ac, None)?;
        let ref_seq = self.hdp.get_seq(transcript_ac, 0, -1, IdentifierKind::Transcript.into_identifier_type())?;

        let cds_start_idx = transcript.cds_start_index().ok_or_else(|| HgvsError::ValidationError("Missing CDS start".into()))?.0 as usize;
        let cds_end_idx = transcript.cds_end_index().ok_or_else(|| HgvsError::ValidationError("Missing CDS end".into()))?.0 as usize;

        if ref_seq.len() < cds_end_idx {
            return Err(HgvsError::ValidationError(format!("Transcript sequence too short (len={}, expected at least {})", ref_seq.len(), cds_end_idx)));
        }

        if cds_start_idx > ref_seq.len() {
            return Err(HgvsError::ValidationError(format!("CDS start {} out of sequence bounds {}", cds_start_idx, ref_seq.len())));
        }

        let ref_aa = crate::utils::translate_cds(&ref_seq[cds_start_idx..]);


        let builder = AltSeqBuilder {
            var_c,
            transcript_sequence: ref_seq,
            cds_start_index: transcript.cds_start_index().unwrap(),
            cds_end_index: transcript.cds_end_index().unwrap(),
            protein_accession: pro_ac_str,
        };
        let alt_data = builder.build_altseq()?;

        let hgvsp_builder = AltSeqToHgvsp {
            ref_aa,
            alt_data: &alt_data,
        };
        let mut var_p = hgvsp_builder.build_hgvsp()?;
        var_p.posedit.predicted = true;
        Ok(var_p)
    }

    /// Normalizes a variant to its 3' most position.
    pub fn normalize_variant(&self, var: crate::SequenceVariant) -> Result<crate::SequenceVariant, HgvsError> {
        match var {
            crate::SequenceVariant::Coding(mut v_c) => {
                let transcript = self.hdp.get_transcript(&v_c.ac, None)?;
                if let Some(pos) = &mut v_c.posedit.pos {
                    let (start_idx, end_idx) = self.get_c_indices(pos, &transcript)?;
                    let is_ins = matches!(&v_c.posedit.edit, crate::edits::NaEdit::Ins { .. });
                    let actual_end = if is_ins { end_idx - 1 } else { end_idx };

                    let (new_start, _new_end) = self.shift_3_prime(&v_c.ac, IdentifierKind::Transcript, start_idx, actual_end, &v_c.posedit.edit)?;

                    if new_start != start_idx {
                        let shift = (new_start as i32) - (start_idx as i32);
                        pos.start.base.0 += shift;
                        if let Some(e) = &mut pos.end {
                            e.base.0 += shift;
                        }
                    }

                    // Update sequences for Del/Dup to match reference at new position
                    match &mut v_c.posedit.edit {
                        crate::edits::NaEdit::Del { ref_: r, .. } | crate::edits::NaEdit::Dup { ref_: r, .. } => {
                            let seq = self.hdp.get_seq(&v_c.ac, new_start as i32, (new_start + (end_idx - start_idx)) as i32, IdentifierKind::Transcript.into_identifier_type())?;
                            *r = Some(seq);
                        }
                        _ => {}
                    }
                }
                Ok(crate::SequenceVariant::Coding(v_c))
            }
            crate::SequenceVariant::Genomic(mut v_g) => {
                if let Some(pos) = &mut v_g.posedit.pos {
                    let start_idx = pos.start.base.to_index().0 as usize;
                    let is_ins = matches!(&v_g.posedit.edit, crate::edits::NaEdit::Ins { .. });
                    let end_idx = pos.end.as_ref().map_or(start_idx + 1, |e| {
                        let idx = e.base.to_index().0 as usize;
                        if is_ins { idx } else { idx + 1 }
                    });

                    let (new_start, new_end) = self.shift_3_prime(&v_g.ac, IdentifierKind::Genomic, start_idx, end_idx, &v_g.posedit.edit)?;
                    if new_start != start_idx {
                        let shift = (new_start as i32) - (start_idx as i32);
                        pos.start.base.0 += shift;
                        if let Some(e) = &mut pos.end {
                            e.base.0 += shift;
                        }
                    }

                    // Update sequences for Del/Dup to match reference at new position
                    match &mut v_g.posedit.edit {
                        crate::edits::NaEdit::Del { ref_: r, .. } | crate::edits::NaEdit::Dup { ref_: r, .. } => {
                            let seq = self.hdp.get_seq(&v_g.ac, new_start as i32, new_end as i32, IdentifierKind::Genomic.into_identifier_type())?;
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

                    let (new_start, new_end) = self.shift_3_prime(&v_n.ac, IdentifierKind::Transcript, start_idx, actual_end, &v_n.posedit.edit)?;

                    if new_start != start_idx {
                        let shift = (new_start as i32) - (start_idx as i32);
                        pos.start.base.0 += shift;
                        if let Some(e) = &mut pos.end {
                            e.base.0 += shift;
                        }
                    }

                    // Update sequences for Del/Dup to match reference at new position
                    match &mut v_n.posedit.edit {
                        crate::edits::NaEdit::Del { ref_: r, .. } | crate::edits::NaEdit::Dup { ref_: r, .. } => {
                            let seq = self.hdp.get_seq(&v_n.ac, new_start as i32, new_end as i32, IdentifierKind::Transcript.into_identifier_type())?;
                            *r = Some(seq);
                        }
                        _ => {}
                    }
                }
                Ok(crate::SequenceVariant::NonCoding(v_n))
            }
            _ => Ok(var)
        }
    }

    pub fn get_c_indices(&self, pos: &BaseOffsetInterval, transcript: &Box<dyn Transcript>) -> Result<(usize, usize), HgvsError> {
        let am = TranscriptMapper::new(dyn_clone::clone_box(&**transcript))?;
        let n_start = am.c_to_n(pos.start.base.to_index(), pos.start.anchor)?;
        let n_end = if let Some(e) = &pos.end {
            am.c_to_n(e.base.to_index(), e.anchor)?
        } else {
            n_start
        };
        Ok((n_start.0 as usize, (n_end.0 + 1) as usize))
    }

    pub fn get_n_indices(&self, pos: &BaseOffsetInterval, transcript: &Box<dyn Transcript>) -> Result<(usize, usize), HgvsError> {
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

    fn shift_3_prime(&self, ac: &str, kind: IdentifierKind, start: usize, end: usize, edit: &crate::edits::NaEdit) -> Result<(usize, usize), HgvsError> {
        let (ref_str, alt_str) = match edit {
            crate::edits::NaEdit::RefAlt { ref_, alt, .. } => (ref_.as_deref().unwrap_or(""), alt.as_deref().unwrap_or("")),
            crate::edits::NaEdit::Del { ref_: Some(s), .. } => (s.as_str(), ""),
            crate::edits::NaEdit::Del { ref_: None, .. } => ("", ""),
            crate::edits::NaEdit::Ins { alt: Some(s), .. } => ("", s.as_str()),
            crate::edits::NaEdit::Ins { alt: None, .. } => ("", ""),
            crate::edits::NaEdit::Dup { ref_: Some(s), .. } => (s.as_str(), ""),
            crate::edits::NaEdit::Dup { ref_: None, .. } => ("", ""),
            _ => return Ok((start, end)),
        };

        if ref_str == alt_str && matches!(edit, crate::edits::NaEdit::RefAlt { .. }) { return Ok((start, end)); }

        let mut curr_start = start;
        let mut curr_end = end;
        let mut chunk_size = 128;

        let mut chunk_start = end;
        let mut chunk = self.hdp.get_seq(ac, chunk_start as i32, (chunk_start + chunk_size) as i32, kind.into_identifier_type())?;
        let mut chunk_bytes = chunk.as_bytes();

        let is_del_or_dup = matches!(edit, crate::edits::NaEdit::Del { .. } | crate::edits::NaEdit::Dup { .. });

        if is_del_or_dup || (!ref_str.is_empty() && alt_str.is_empty()) {
            // Deletion or Duplication (explicit or implicit)
            let mut ref_chunk = self.hdp.get_seq(ac, curr_start as i32, (curr_start + 1) as i32, kind.into_identifier_type())?;
            if ref_chunk.is_empty() { return Ok((curr_start, curr_end)); }
            let mut ref_byte = ref_chunk.as_bytes()[0];

            loop {
                if (curr_end - chunk_start) >= chunk_bytes.len() {
                    if chunk_bytes.len() < chunk_size { break; }
                    chunk_start += chunk_bytes.len();
                    chunk_size = std::cmp::min(chunk_size * 2, 4096);
                    chunk = self.hdp.get_seq(ac, chunk_start as i32, (chunk_start + chunk_size) as i32, kind.into_identifier_type())?;
                    chunk_bytes = chunk.as_bytes();
                    if chunk_bytes.is_empty() { break; }
                }

                if ref_byte == chunk_bytes[curr_end - chunk_start] {
                    curr_start += 1;
                    curr_end += 1;
                    ref_chunk = self.hdp.get_seq(ac, curr_start as i32, (curr_start + 1) as i32, kind.into_identifier_type())?;
                    if ref_chunk.is_empty() { break; }
                    ref_byte = ref_chunk.as_bytes()[0];
                } else {
                    break;
                }
            }
        } else if (ref_str.is_empty() && start == end) && !alt_str.is_empty() {
            // Insertion
            let alt_bytes = alt_str.as_bytes();
            let n = alt_bytes.len();
            if n == 0 { return Ok((curr_start, curr_end)); }
            let mut i = 0;
            loop {
                if (curr_end - chunk_start) >= chunk_bytes.len() {
                    if chunk_bytes.len() < chunk_size { break; }
                    chunk_start += chunk_bytes.len();
                    chunk_size = std::cmp::min(chunk_size * 2, 4096);
                    chunk = self.hdp.get_seq(ac, chunk_start as i32, (chunk_start + chunk_size) as i32, kind.into_identifier_type())?;
                    chunk_bytes = chunk.as_bytes();
                    if chunk_bytes.is_empty() { break; }
                }

                if chunk_bytes[curr_end - chunk_start] == alt_bytes[i % n] {
                    curr_start += 1;
                    curr_end += 1;
                    i += 1;
                } else {
                    break;
                }
            }
        }
        Ok((curr_start, curr_end))
    }
}
