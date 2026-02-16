use crate::data::Transcript;
use crate::error::HgvsError;
use crate::structs::{Anchor, GenomicPos, IntronicOffset, TranscriptPos};

/// Handles coordinate transformations within a single transcript.
pub struct TranscriptMapper {
    /// The transcript model providing exon and CDS information.
    pub transcript: Box<dyn Transcript>,
}

impl TranscriptMapper {
    /// Creates a new `TranscriptMapper` for the given transcript.
    pub fn new(transcript: Box<dyn Transcript>) -> Result<Self, HgvsError> {
        Ok(TranscriptMapper { transcript })
    }

    /// Maps a 0-based genomic position to a 0-based transcript position and intronic offset.
    pub fn g_to_n(&self, g_pos: GenomicPos) -> Result<(TranscriptPos, IntronicOffset), HgvsError> {
        let mut n_pos = 0;
        for exon in self.transcript.exons() {
            let (e_start, e_end) = (exon.reference_start, exon.reference_end);
            // e_start and e_end are 0-based inclusive
            if g_pos.0 >= e_start.0 && g_pos.0 <= e_end.0 {
                let offset_in_exon = if exon.alt_strand == 1 {
                    g_pos.0 - e_start.0
                } else {
                    e_end.0 - g_pos.0
                };
                return Ok((TranscriptPos(n_pos + offset_in_exon), IntronicOffset(0)));
            }
            n_pos += (e_end.0 - e_start.0) + 1;
        }

        // Handle intronic positions by finding the nearest exon
        let mut best_n = TranscriptPos(0);
        let mut best_offset = i32::MAX;
        let mut curr_n = 0;
        for exon in self.transcript.exons() {
            let (e_start, e_end) = (exon.reference_start, exon.reference_end);
            let d_start = (g_pos.0 - e_start.0).abs();
            let d_end = (g_pos.0 - e_end.0).abs();
            let d = d_start.min(d_end);
            if d < best_offset.abs() {
                best_offset = if g_pos.0 < e_start.0 {
                    if exon.alt_strand == 1 {
                        g_pos.0 - e_start.0
                    } else {
                        e_start.0 - g_pos.0
                    }
                } else {
                    if exon.alt_strand == 1 {
                        g_pos.0 - e_end.0
                    } else {
                        e_end.0 - g_pos.0
                    }
                };
                best_n = if g_pos.0 < e_start.0 {
                    TranscriptPos(curr_n)
                } else {
                    TranscriptPos(curr_n + (e_end.0 - e_start.0))
                };
            }
            curr_n += (e_end.0 - e_start.0) + 1;
        }
        Ok((best_n, IntronicOffset(best_offset)))
    }

    /// Maps a 0-based transcript position to a 0-based cDNA position and anchor.
    pub fn n_to_c(
        &self,
        n_pos: TranscriptPos,
    ) -> Result<(TranscriptPos, IntronicOffset, Anchor), HgvsError> {
        if let (Some(cds_start), Some(cds_end)) = (
            self.transcript.cds_start_index(),
            self.transcript.cds_end_index(),
        ) {
            if n_pos < cds_start {
                Ok((
                    TranscriptPos(n_pos.0 - cds_start.0),
                    IntronicOffset(0),
                    Anchor::CdsStart,
                ))
            } else if n_pos > cds_end {
                Ok((
                    TranscriptPos(n_pos.0 - cds_end.0 - 1),
                    IntronicOffset(0),
                    Anchor::CdsEnd,
                ))
            } else {
                Ok((
                    TranscriptPos(n_pos.0 - cds_start.0),
                    IntronicOffset(0),
                    Anchor::CdsStart,
                ))
            }
        } else {
            Ok((n_pos, IntronicOffset(0), Anchor::TranscriptStart))
        }
    }

    /// Maps a cDNA position and anchor to a 0-based transcript position.
    pub fn c_to_n(&self, c_pos: TranscriptPos, anchor: Anchor) -> Result<TranscriptPos, HgvsError> {
        match anchor {
            Anchor::TranscriptStart => Ok(c_pos),
            Anchor::CdsStart => {
                let cds_start = self
                    .transcript
                    .cds_start_index()
                    .ok_or_else(|| HgvsError::ValidationError("Missing CDS start".into()))?;
                Ok(TranscriptPos(cds_start.0 + c_pos.0))
            }
            Anchor::CdsEnd => {
                let cds_end = self
                    .transcript
                    .cds_end_index()
                    .ok_or_else(|| HgvsError::ValidationError("Missing CDS end".into()))?;
                Ok(TranscriptPos(cds_end.0 + 1 + c_pos.0))
            }
        }
    }

    /// Maps a 0-based transcript position and offset to a 0-based genomic position.
    pub fn n_to_g(
        &self,
        n_pos: TranscriptPos,
        offset: IntronicOffset,
    ) -> Result<GenomicPos, HgvsError> {
        let mut curr_n = 0;
        for exon in self.transcript.exons() {
            let (e_start, e_end) = (exon.reference_start, exon.reference_end);
            let e_len = (e_end.0 - e_start.0) + 1;
            if n_pos.0 >= curr_n && n_pos.0 < curr_n + e_len {
                let offset_in_exon = n_pos.0 - curr_n;
                let g_base = if exon.alt_strand == 1 {
                    e_start.0 + offset_in_exon
                } else {
                    e_end.0 - offset_in_exon
                };
                return Ok(GenomicPos(
                    g_base
                        + if exon.alt_strand == 1 {
                            offset.0
                        } else {
                            -offset.0
                        },
                ));
            }
            curr_n += e_len;
        }
        Err(HgvsError::ValidationError(
            "Transcript position out of exon bounds".into(),
        ))
    }
}
