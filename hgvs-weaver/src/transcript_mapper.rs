use crate::data::{ExonData, Transcript};
use crate::error::HgvsError;
use crate::structs::{Anchor, GenomicPos, IntronicOffset, TranscriptPos};

/// Handles coordinate transformations within a single transcript.
pub struct TranscriptMapper {
    /// The transcript model providing exon and CDS information.
    pub transcript: Box<dyn Transcript>,
    /// Sorted exons (transcript order).
    pub exons: Vec<ExonData>,
    /// CIGAR mappers for exons with non-trivial alignments.
    pub cigar_mappers: Vec<Option<crate::cigar::CigarMapper>>,
}

impl TranscriptMapper {
    /// Creates a new `TranscriptMapper` for the given transcript.
    pub fn new(transcript: Box<dyn Transcript>) -> Result<Self, HgvsError> {
        let mut exons = transcript.exons().to_vec();
        if transcript.strand() == crate::data::Strand::Plus {
            exons.sort_by_key(|e| e.reference_start.0);
        } else {
            exons.sort_by_key(|e| std::cmp::Reverse(e.reference_start.0));
        }
        let mut cigar_mappers = Vec::with_capacity(exons.len());
        for exon in &exons {
            if !exon.cigar.is_empty()
                && exon.cigar != format!("{}M", (exon.reference_end.0 - exon.reference_start.0) + 1)
                && exon.cigar != format!("{}=", (exon.reference_end.0 - exon.reference_start.0) + 1)
            {
                cigar_mappers.push(Some(crate::cigar::CigarMapper::new(&exon.cigar)?));
            } else {
                cigar_mappers.push(None);
            }
        }
        Ok(TranscriptMapper {
            transcript,
            exons,
            cigar_mappers,
        })
    }

    /// Maps a 0-based genomic position to a 0-based transcript position and intronic offset.
    pub fn g_to_n(&self, g_pos: GenomicPos) -> Result<(TranscriptPos, IntronicOffset), HgvsError> {
        let mut n_pos = 0;
        for (i, exon) in self.exons.iter().enumerate() {
            let (e_start, e_end) = (exon.reference_start, exon.reference_end);
            // e_start and e_end are 0-based inclusive
            if g_pos.0 >= e_start.0 && g_pos.0 <= e_end.0 {
                let offset_in_exon = if let Some(cm) = &self.cigar_mappers[i] {
                    let g_offset = if exon.alt_strand == crate::data::Strand::Plus {
                        g_pos.0 - e_start.0
                    } else {
                        e_end.0 - g_pos.0
                    };
                    let (t_offset, _intronic, _op) = cm.map_ref_to_tgt(g_offset, "start", true)?;
                    t_offset
                } else {
                    if exon.alt_strand == crate::data::Strand::Plus {
                        g_pos.0 - e_start.0
                    } else {
                        e_end.0 - g_pos.0
                    }
                };
                return Ok((TranscriptPos(n_pos + offset_in_exon), IntronicOffset(0)));
            }
            n_pos += if let Some(cm) = &self.cigar_mappers[i] {
                cm.tgt_len()
            } else {
                (e_end.0 - e_start.0) + 1
            };
        }

        // Handle intronic positions by finding the nearest exon
        let mut best_n = TranscriptPos(0);
        let mut best_dist = i32::MAX;
        let mut best_offset = 0;
        let mut curr_n = 0;
        for (i, exon) in self.exons.iter().enumerate() {
            let (e_start, e_end) = (exon.reference_start, exon.reference_end);
            let d_start = (g_pos.0 - e_start.0).abs();
            let d_end = (g_pos.0 - e_end.0).abs();
            let d = d_start.min(d_end);
            if d < best_dist {
                best_dist = d;
                best_offset = if g_pos.0 < e_start.0 {
                    if exon.alt_strand == crate::data::Strand::Plus {
                        g_pos.0 - e_start.0
                    } else {
                        e_start.0 - g_pos.0
                    }
                } else {
                    if exon.alt_strand == crate::data::Strand::Plus {
                        g_pos.0 - e_end.0
                    } else {
                        e_end.0 - g_pos.0
                    }
                };
                best_n = if g_pos.0 < e_start.0 {
                    TranscriptPos(curr_n)
                } else {
                    let e_tgt_len = if let Some(cm) = &self.cigar_mappers[i] {
                        cm.tgt_len()
                    } else {
                        (e_end.0 - e_start.0) + 1
                    };
                    TranscriptPos(curr_n + e_tgt_len - 1)
                };
            }
            curr_n += if let Some(cm) = &self.cigar_mappers[i] {
                cm.tgt_len()
            } else {
                (e_end.0 - e_start.0) + 1
            };
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
        for (i, exon) in self.exons.iter().enumerate() {
            let (e_start, e_end) = (exon.reference_start, exon.reference_end);
            let e_tgt_len = if let Some(cm) = &self.cigar_mappers[i] {
                cm.tgt_len()
            } else {
                (e_end.0 - e_start.0) + 1
            };
            if n_pos.0 >= curr_n && n_pos.0 < curr_n + e_tgt_len {
                let offset_in_exon_tgt = n_pos.0 - curr_n;
                let g_offset = if let Some(cm) = &self.cigar_mappers[i] {
                    let (g_off, _intronic, _op) =
                        cm.map_tgt_to_ref(offset_in_exon_tgt, "start", true)?;
                    g_off
                } else {
                    offset_in_exon_tgt
                };
                let g_base = if exon.alt_strand == crate::data::Strand::Plus {
                    e_start.0 + g_offset
                } else {
                    e_end.0 - g_offset
                };
                return Ok(GenomicPos(
                    g_base
                        + if exon.alt_strand == crate::data::Strand::Plus {
                            offset.0
                        } else {
                            -offset.0
                        },
                ));
            }
            curr_n += e_tgt_len;
        }
        Err(HgvsError::ValidationError(
            "Transcript position out of exon bounds".into(),
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data::{ExonData, TranscriptData};

    fn create_mock_transcript(
        strand: crate::data::Strand,
        exons: Vec<ExonData>,
    ) -> Box<dyn Transcript> {
        Box::new(TranscriptData {
            ac: "NM_0001.1".to_string(),
            gene: "TEST".to_string(),
            cds_start_index: None,
            cds_end_index: None,
            strand,
            reference_accession: "NC_000001.1".to_string(),
            exons,
        })
    }

    #[test]
    fn test_g_to_n_minus_strand_order() {
        // Exons in genomic ascending order: [1000, 1100], [2000, 2100]
        // For minus strand, transcript order should be [2000, 2100] then [1000, 1100]
        let exons = vec![
            ExonData {
                transcript_start: TranscriptPos(0), // Placeholder
                transcript_end: TranscriptPos(100),
                reference_start: GenomicPos(1000),
                reference_end: GenomicPos(1100),
                alt_strand: crate::data::Strand::Minus,
                cigar: "101M".to_string(),
            },
            ExonData {
                transcript_start: TranscriptPos(101),
                transcript_end: TranscriptPos(201),
                reference_start: GenomicPos(2000),
                reference_end: GenomicPos(2100),
                alt_strand: crate::data::Strand::Minus,
                cigar: "101M".to_string(),
            },
        ];
        let tx = create_mock_transcript(crate::data::Strand::Minus, exons);
        let mapper = TranscriptMapper::new(tx).unwrap();

        // Genomic 2100 should be n.0
        let (n_pos, offset) = mapper.g_to_n(GenomicPos(2100)).unwrap();
        // CURRENT BEHAVIOR (BUGGY):
        // It visits [1000, 1100] first, n_pos starts at 0.
        // Then it visits [2000, 2100], n_pos = 101.
        // offset_in_exon = 2100 - 2100 = 0 (for alt_strand = -1)
        // Result: n.101.
        // EXPECTED: n.0
        assert_eq!(n_pos.0, 0, "Minus strand genomic 2100 should be n.0");
        assert_eq!(offset.0, 0);
    }

    #[test]
    fn test_g_to_n_intronic_offset() {
        let exons = vec![ExonData {
            transcript_start: TranscriptPos(0),
            transcript_end: TranscriptPos(10),
            reference_start: GenomicPos(1000),
            reference_end: GenomicPos(1010),
            alt_strand: crate::data::Strand::Plus,
            cigar: "11M".to_string(),
        }];
        let tx = create_mock_transcript(crate::data::Strand::Plus, exons);
        let mapper = TranscriptMapper::new(tx).unwrap();

        // Genomic 999 is 1bp upstream of exon start (1000)
        let (n_pos, offset) = mapper.g_to_n(GenomicPos(999)).unwrap();
        assert_eq!(n_pos.0, 0);
        assert_eq!(offset.0, -1);

        // Genomic 1011 is 1bp downstream of exon end (1010)
        let (n_pos, offset) = mapper.g_to_n(GenomicPos(1011)).unwrap();
        assert_eq!(n_pos.0, 10);
        assert_eq!(offset.0, 1);
    }

    #[test]
    fn test_g_to_n_cigar() {
        // Exon: Genomic [1000, 1010] (11bp)
        // CIGAR: 5=1D5= (Ref 11bp -> Tgt 10bp)
        // Ref index: 0 1 2 3 4 5 6 7 8 9 10
        // Tgt index: 0 1 2 3 4 _ 5 6 7 8 9
        let exons = vec![ExonData {
            transcript_start: TranscriptPos(0),
            transcript_end: TranscriptPos(9),
            reference_start: GenomicPos(1000),
            reference_end: GenomicPos(1010),
            alt_strand: crate::data::Strand::Plus,
            cigar: "5=1D5=".to_string(),
        }];
        let tx = create_mock_transcript(crate::data::Strand::Plus, exons);
        let mapper = TranscriptMapper::new(tx).unwrap();

        // g.1000 -> n.0
        assert_eq!(mapper.g_to_n(GenomicPos(1000)).unwrap().0 .0, 0);
        // g.1004 -> n.4
        assert_eq!(mapper.g_to_n(GenomicPos(1004)).unwrap().0 .0, 4);
        // g.1005 is deleted in transcript (1D)
        // CigarMapper::map_ref_to_tgt(5, "start", true) for "5=1D5="
        // ref_pos: [0, 5, 6, 11]
        // tgt_pos: [0, 5, 5, 10]
        // pos=5 is in op 1 (1D). end_strategy="start" -> mapped_pos = tgt_pos[1] - 1 = 4.
        assert_eq!(mapper.g_to_n(GenomicPos(1005)).unwrap().0 .0, 4);
        // g.1006 -> n.5
        assert_eq!(mapper.g_to_n(GenomicPos(1006)).unwrap().0 .0, 5);
    }
}
