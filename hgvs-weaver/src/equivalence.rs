use crate::error::HgvsError;
use crate::data::{DataProvider, TranscriptSearch, IdentifierKind};
use crate::structs::{
    SequenceVariant, GVariant, PVariant, Variant,
    SimpleInterval, SimplePosition, BaseOffsetInterval, BaseOffsetPosition,
    GenomicPos, TranscriptPos, NaEdit
};
use crate::mapper::VariantMapper;

pub struct VariantEquivalence<'a> {
    pub hdp: &'a dyn DataProvider,
    pub searcher: &'a dyn TranscriptSearch,
    pub mapper: VariantMapper<'a>,
}

impl<'a> VariantEquivalence<'a> {
    pub fn new(hdp: &'a dyn DataProvider, searcher: &'a dyn TranscriptSearch) -> Self {
        VariantEquivalence {
            hdp,
            searcher,
            mapper: VariantMapper::new(hdp),
        }
    }

    pub fn equivalent(&self, var1: &SequenceVariant, var2: &SequenceVariant) -> Result<bool, HgvsError> {
        // Expand gene symbols if present
        let vars1 = self.expand_if_gene_symbol(var1)?;
        let vars2 = self.expand_if_gene_symbol(var2)?;

        for v1 in &vars1 {
            for v2 in &vars2 {
                if self.are_equivalent_single(v1, v2)? {
                    return Ok(true);
                }
            }
        }
        Ok(false)
    }

    /// Fills in missing sequence information for deletions and duplications.
    fn fill_implicit_sequence(&self, var: &SequenceVariant) -> Result<SequenceVariant, HgvsError> {
        match var {
            SequenceVariant::Genomic(v) => {
                if let Some(pos) = &v.posedit.pos {
                    let start = pos.start.base.to_index().0 as usize;
                    let end = pos.end.as_ref().map_or(start + 1, |e| e.base.to_index().0 as usize + 1);

                    let mut new_v = v.clone();
                    new_v.posedit.edit = self.fill_na_edit(&v.ac, IdentifierKind::Genomic, start, end, v.posedit.edit.clone())?;
                    Ok(SequenceVariant::Genomic(new_v))
                } else {
                    Ok(var.clone())
                }
            },
            SequenceVariant::Coding(v) => {
                let transcript = self.hdp.get_transcript(&v.ac, None)?;
                if let Some(pos) = &v.posedit.pos {
                     let (start, end) = self.mapper.get_c_indices(pos, &transcript)?;
                     let mut new_v = v.clone();
                     new_v.posedit.edit = self.fill_na_edit(&v.ac, IdentifierKind::Transcript, start, end, v.posedit.edit.clone())?;
                     Ok(SequenceVariant::Coding(new_v))
                } else {
                    Ok(var.clone())
                }
            },
            SequenceVariant::NonCoding(v) => {
                 let transcript = self.hdp.get_transcript(&v.ac, None)?;
                if let Some(pos) = &v.posedit.pos {
                     let (start, end) = self.mapper.get_n_indices(pos, &transcript)?;
                     let mut new_v = v.clone();
                     new_v.posedit.edit = self.fill_na_edit(&v.ac, IdentifierKind::Transcript, start, end, v.posedit.edit.clone())?;
                     Ok(SequenceVariant::NonCoding(new_v))
                } else {
                    Ok(var.clone())
                }
            },
            SequenceVariant::Rna(v) => {
                 let transcript = self.hdp.get_transcript(&v.ac, None)?;
                if let Some(pos) = &v.posedit.pos {
                     let (start, end) = self.mapper.get_n_indices(pos, &transcript)?;
                     let mut new_v = v.clone();
                     new_v.posedit.edit = self.fill_na_edit(&v.ac, IdentifierKind::Transcript, start, end, v.posedit.edit.clone())?;
                     Ok(SequenceVariant::Rna(new_v))
                } else {
                    Ok(var.clone())
                }
            },
            SequenceVariant::Mitochondrial(v) => {
                if let Some(pos) = &v.posedit.pos {
                    let start = pos.start.base.to_index().0 as usize;
                    let end = pos.end.as_ref().map_or(start + 1, |e| e.base.to_index().0 as usize + 1);

                    let mut new_v = v.clone();
                    new_v.posedit.edit = self.fill_na_edit(&v.ac, IdentifierKind::Genomic, start, end, v.posedit.edit.clone())?;
                    Ok(SequenceVariant::Mitochondrial(new_v))
                } else {
                    Ok(var.clone())
                }
            },
            _ => Ok(var.clone()),
        }
    }

    fn fill_na_edit(&self, ac: &str, kind: IdentifierKind, start: usize, end: usize, edit: crate::edits::NaEdit) -> Result<crate::edits::NaEdit, HgvsError> {
         match edit {
            crate::edits::NaEdit::Del { ref_: None, uncertain } => {
                let seq = self.hdp.get_seq(ac, start as i32, end as i32, kind.into_identifier_type())?;
                Ok(crate::edits::NaEdit::Del { ref_: Some(seq), uncertain })
            },
             crate::edits::NaEdit::Dup { ref_: None, uncertain } => {
                let seq = self.hdp.get_seq(ac, start as i32, end as i32, kind.into_identifier_type())?;
                Ok(crate::edits::NaEdit::Dup { ref_: Some(seq), uncertain })
            },
            _ => Ok(edit)
         }
    }


    fn expand_if_gene_symbol(&self, var: &SequenceVariant) -> Result<Vec<SequenceVariant>, HgvsError> {
        let ac = var.ac();

        // Use DataProvider to determine if this is a symbol or an accession.
        let id_type = self.hdp.get_identifier_type(ac)?;

        if id_type == crate::data::IdentifierType::GeneSymbol {
            let target_kind = match var {
                SequenceVariant::Protein(_) => IdentifierKind::Protein,
                SequenceVariant::Coding(_) | SequenceVariant::NonCoding(_) | SequenceVariant::Rna(_) => IdentifierKind::Transcript,
                _ => IdentifierKind::Genomic,
            };

            // Try symbol expansion.
            let accessions = self.hdp.get_symbol_accessions(ac, IdentifierKind::Genomic, target_kind)?;

            if !accessions.is_empty() {
                let mut expanded = Vec::new();
                for (ac_type, new_ac) in accessions {
                    // Only expand into compatible types.
                    let is_compatible = match (var, ac_type) {
                        (SequenceVariant::Protein(_), crate::data::IdentifierType::ProteinAccession) => true,
                        (SequenceVariant::Coding(_) | SequenceVariant::NonCoding(_) | SequenceVariant::Rna(_), crate::data::IdentifierType::TranscriptAccession) => true,
                        (SequenceVariant::Genomic(_) | SequenceVariant::Mitochondrial(_), crate::data::IdentifierType::GenomicAccession) => true,
                        // Allow g. on transcripts if specifically provided
                        (SequenceVariant::Genomic(_), crate::data::IdentifierType::TranscriptAccession) => true,
                        _ => false,
                    };

                    if is_compatible {
                        let mut v = var.clone();
                        match &mut v {
                            SequenceVariant::Genomic(v_g) => v_g.ac = new_ac,
                            SequenceVariant::Coding(v_c) => v_c.ac = new_ac,
                            SequenceVariant::Protein(v_p) => v_p.ac = new_ac,
                            SequenceVariant::Mitochondrial(v_m) => v_m.ac = new_ac,
                            SequenceVariant::NonCoding(v_n) => v_n.ac = new_ac,
                            SequenceVariant::Rna(v_r) => v_r.ac = new_ac,
                        }
                        expanded.push(v);
                    }
                }
                if !expanded.is_empty() {
                    return Ok(expanded);
                }
            }
        }

        // Default: return as-is.
        Ok(vec![var.clone()])
    }

    fn are_equivalent_single(&self, var1: &SequenceVariant, var2: &SequenceVariant) -> Result<bool, HgvsError> {
        match (var1, var2) {
            // Nucleotide vs Nucleotide
            (SequenceVariant::Genomic(v1), SequenceVariant::Genomic(v2)) => self.n_vs_n_equivalent(v1, v2),
            (SequenceVariant::Coding(v1), SequenceVariant::Coding(v2)) => self.n_vs_n_equivalent_c(v1, v2),
            (SequenceVariant::NonCoding(v1), SequenceVariant::NonCoding(v2)) => self.n_vs_n_equivalent_n(v1, v2),

            (SequenceVariant::Genomic(v1), SequenceVariant::Coding(v2)) => self.g_vs_c_equivalent(v1, v2),
            (SequenceVariant::Coding(v1), SequenceVariant::Genomic(v2)) => self.g_vs_c_equivalent(v2, v1),

            (SequenceVariant::Genomic(v1), SequenceVariant::NonCoding(v2)) => self.g_vs_n_equivalent(v1, v2),
            (SequenceVariant::NonCoding(v1), SequenceVariant::Genomic(v2)) => self.g_vs_n_equivalent(v2, v1),

            (SequenceVariant::Coding(v1), SequenceVariant::NonCoding(v2)) => self.c_vs_n_equivalent(v1, v2),
            (SequenceVariant::NonCoding(v1), SequenceVariant::Coding(v2)) => self.c_vs_n_equivalent(v2, v1),

            // Nucleotide vs Protein
            (SequenceVariant::Genomic(v1), SequenceVariant::Protein(v2)) => self.g_vs_p_equivalent(v1, v2),
            (SequenceVariant::Protein(v1), SequenceVariant::Genomic(v2)) => self.g_vs_p_equivalent(v2, v1),
            (SequenceVariant::Coding(v1), SequenceVariant::Protein(v2)) => self.c_vs_p_equivalent(v1, v2),
            (SequenceVariant::Protein(v1), SequenceVariant::Coding(v2)) => self.c_vs_p_equivalent(v2, v1),
            (SequenceVariant::NonCoding(v1), SequenceVariant::Protein(v2)) => self.n_vs_p_equivalent(v1, v2),
            (SequenceVariant::Protein(v1), SequenceVariant::NonCoding(v2)) => self.n_vs_p_equivalent(v2, v1),

            // Protein vs Protein
            (SequenceVariant::Protein(v1), SequenceVariant::Protein(v2)) => self.p_vs_p_equivalent(v1, v2),

            // Fallback
            _ => {
                if var1.ac() == var2.ac() && var1.coordinate_type() == var2.coordinate_type() {
                    return Ok(self.normalize_format(&var1.to_string()) == self.normalize_format(&var2.to_string()));
                }
                Ok(false)
            }
        }
    }

    fn normalize_format(&self, s: &str) -> String {
        let mut s = s.replace(['(', ')', '?'], "");
        // Normalize 3-letter AA codes to 1-letter
        let replacements = [
            ("Ala", "A"), ("Arg", "R"), ("Asn", "N"), ("Asp", "D"), ("Cys", "C"),
            ("Gln", "Q"), ("Glu", "E"), ("Gly", "G"), ("His", "H"), ("Ile", "I"),
            ("Leu", "L"), ("Lys", "K"), ("Met", "M"), ("Phe", "F"), ("Pro", "P"),
            ("Ser", "S"), ("Thr", "T"), ("Trp", "W"), ("Tyr", "Y"), ("Val", "V"),
            ("Asx", "B"), ("Glx", "Z"), ("Xaa", "X"), ("Ter", "*"),
        ];
        for (from, to) in replacements {
            s = s.replace(from, to);
        }
        s
    }

    fn n_vs_n_equivalent(&self, v1: &GVariant, v2: &GVariant) -> Result<bool, HgvsError> {
        let nv1 = self.mapper.normalize_variant(SequenceVariant::Genomic(v1.clone()))?;
        let nv2 = self.mapper.normalize_variant(SequenceVariant::Genomic(v2.clone()))?;

        // Fill implicit sequences
        let nv1_filled = self.fill_implicit_sequence(&nv1)?;
        let nv2_filled = self.fill_implicit_sequence(&nv2)?;

        // Normalize Ins to Dup
        let nv1_dup = self.normalize_ins_to_dup(&nv1_filled)?;
        let nv2_dup = self.normalize_ins_to_dup(&nv2_filled)?;

        let s1 = self.normalize_format(&nv1_dup.to_string());
        let s2 = self.normalize_format(&nv2_dup.to_string());

        Ok(s1 == s2)
    }

    fn normalize_ins_to_dup(&self, var: &SequenceVariant) -> Result<SequenceVariant, HgvsError> {
        match var {
            SequenceVariant::Genomic(v) => {
                if let Some(pos) = &v.posedit.pos {
                    if let NaEdit::Ins { alt: Some(seq), uncertain } = &v.posedit.edit {
                        let start_0 = pos.start.base.to_index();
                        if let Some((check_start, start_idx, edit)) = self.try_normalize_to_dup(&v.ac, IdentifierKind::Genomic, start_0.0, seq, *uncertain)? {
                             let mut new_v = v.clone();
                             new_v.posedit.pos = Some(SimpleInterval {
                                 start: SimplePosition {
                                     base: GenomicPos(check_start).to_hgvs(),
                                     end: None,
                                     uncertain: false,
                                 },
                                 end: if check_start != start_idx {
                                     Some(SimplePosition {
                                         base: GenomicPos(start_idx).to_hgvs(),
                                         end: None,
                                         uncertain: false,
                                     })
                                 } else {
                                     None
                                 },
                                 uncertain: false,
                             });
                             new_v.posedit.edit = edit;
                             return Ok(SequenceVariant::Genomic(new_v));
                        }
                    }
                }
                Ok(var.clone())
            },
            SequenceVariant::Coding(v) => {
                if let Some(pos) = &v.posedit.pos {
                    if let NaEdit::Ins { alt: Some(seq), uncertain } = &v.posedit.edit {
                         if pos.start.offset.is_some() || pos.end.as_ref().map_or(false, |e| e.offset.is_some()) {
                             return Ok(var.clone());
                         }
                         let transcript = self.hdp.get_transcript(&v.ac, None)?;
                         let (start_idx_usize, _) = self.mapper.get_c_indices(pos, &transcript)?;
                         let start_idx = start_idx_usize as i32;

                         if let Some((check_start, last_idx, edit)) = self.try_normalize_to_dup(&v.ac, IdentifierKind::Transcript, start_idx, seq, *uncertain)? {
                              let mut new_v = v.clone();
                              let am = crate::transcript_mapper::TranscriptMapper::new(transcript)?;
                              let (c_pos_index, _, anchor) = am.n_to_c(TranscriptPos(check_start))?;
                              new_v.posedit.pos = Some(BaseOffsetInterval {
                                  start: BaseOffsetPosition {
                                      base: c_pos_index.to_hgvs(),
                                      offset: None,
                                      anchor,
                                      uncertain: false,
                                  },
                                  end: if check_start != last_idx {
                                      let (c_pos_e_index, _, anchor_e) = am.n_to_c(TranscriptPos(last_idx))?;
                                      Some(BaseOffsetPosition {
                                          base: c_pos_e_index.to_hgvs(),
                                          offset: None,
                                          anchor: anchor_e,
                                          uncertain: false,
                                      })
                                  } else {
                                      None
                                  },
                                  uncertain: false,
                              });
                              new_v.posedit.edit = edit;
                              return Ok(SequenceVariant::Coding(new_v));
                         }
                    }
                }
                Ok(var.clone())
            },
            SequenceVariant::NonCoding(v) => {
                if let Some(pos) = &v.posedit.pos {
                    if let NaEdit::Ins { alt: Some(seq), uncertain } = &v.posedit.edit {
                         if pos.start.offset.is_some() || pos.end.as_ref().map_or(false, |e| e.offset.is_some()) {
                             return Ok(var.clone());
                         }
                         let transcript = self.hdp.get_transcript(&v.ac, None)?;
                         let (start_idx_usize, _) = self.mapper.get_n_indices(pos, &transcript)?;
                         let start_idx = start_idx_usize as i32;

                         if let Some((check_start, last_idx, edit)) = self.try_normalize_to_dup(&v.ac, IdentifierKind::Transcript, start_idx, seq, *uncertain)? {
                              let mut new_v = v.clone();
                              let am = crate::transcript_mapper::TranscriptMapper::new(transcript)?;
                              let (c_pos_index, _, anchor) = am.n_to_c(TranscriptPos(check_start))?;
                              new_v.posedit.pos = Some(BaseOffsetInterval {
                                  start: BaseOffsetPosition {
                                      base: c_pos_index.to_hgvs(),
                                      offset: None,
                                      anchor,
                                      uncertain: false,
                                  },
                                  end: if check_start != last_idx {
                                      let (c_pos_e_index, _, anchor_e) = am.n_to_c(TranscriptPos(last_idx))?;
                                      Some(BaseOffsetPosition {
                                          base: c_pos_e_index.to_hgvs(),
                                          offset: None,
                                          anchor: anchor_e,
                                          uncertain: false,
                                      })
                                  } else {
                                      None
                                  },
                                  uncertain: false,
                              });
                              new_v.posedit.edit = edit;
                              return Ok(SequenceVariant::NonCoding(new_v));
                         }
                    }
                }
                Ok(var.clone())
            },
            _ => Ok(var.clone()),
        }
    }

    fn try_normalize_to_dup(&self, ac: &str, kind: IdentifierKind, start_idx: i32, seq: &str, uncertain: bool) -> Result<Option<(i32, i32, NaEdit)>, HgvsError> {
        let len = seq.len() as i32;
        let check_start = start_idx - len + 1;
        if check_start < 0 {
            return Ok(None);
        }
        let ref_seq = self.hdp.get_seq(ac, check_start, start_idx + 1, kind.into_identifier_type())?;
        if ref_seq == *seq {
            Ok(Some((check_start, start_idx, NaEdit::Dup { ref_: Some(seq.to_string()), uncertain })))
        } else {
            Ok(None)
        }
    }

    fn n_vs_n_equivalent_c(&self, v1: &crate::structs::CVariant, v2: &crate::structs::CVariant) -> Result<bool, HgvsError> {
        let tx1 = self.hdp.get_transcript(&v1.ac, None)?;
        let g1 = self.mapper.c_to_g(v1, tx1.reference_accession())?;
        let tx2 = self.hdp.get_transcript(&v2.ac, None)?;
        let g2 = self.mapper.c_to_g(v2, tx2.reference_accession())?;
        self.n_vs_n_equivalent(&g1, &g2)
    }

    fn n_vs_n_equivalent_n(&self, v1: &crate::structs::NVariant, v2: &crate::structs::NVariant) -> Result<bool, HgvsError> {
        let tx1 = self.hdp.get_transcript(&v1.ac, None)?;
        let g1 = self.mapper.n_to_g(v1, tx1.reference_accession())?;
        let tx2 = self.hdp.get_transcript(&v2.ac, None)?;
        let g2 = self.mapper.n_to_g(v2, tx2.reference_accession())?;
        self.n_vs_n_equivalent(&g1, &g2)
    }

    fn g_vs_c_equivalent(&self, vg: &crate::structs::GVariant, vc: &crate::structs::CVariant) -> Result<bool, HgvsError> {
        let g2 = self.mapper.c_to_g(vc, &vg.ac)?;
        self.n_vs_n_equivalent(vg, &g2)
    }

    fn g_vs_n_equivalent(&self, vg: &crate::structs::GVariant, vn: &crate::structs::NVariant) -> Result<bool, HgvsError> {
        let g2 = self.mapper.n_to_g(vn, &vg.ac)?;
        self.n_vs_n_equivalent(vg, &g2)
    }

    fn c_vs_n_equivalent(&self, vc: &crate::structs::CVariant, vn: &crate::structs::NVariant) -> Result<bool, HgvsError> {
        let tx = self.hdp.get_transcript(&vc.ac, None)?;
        let ref_ac = tx.reference_accession();
        let g1 = self.mapper.c_to_g(vc, ref_ac)?;
        let g2 = self.mapper.n_to_g(vn, ref_ac)?;
        self.n_vs_n_equivalent(&g1, &g2)
    }

    fn n_vs_p_equivalent(&self, vn: &crate::structs::NVariant, vp: &crate::structs::PVariant) -> Result<bool, HgvsError> {
        let tx = self.hdp.get_transcript(&vn.ac, None)?;
        let ref_ac = tx.reference_accession();
        let vg = self.mapper.n_to_g(vn, ref_ac)?;
        self.g_vs_p_equivalent(&vg, vp)
    }

    fn g_vs_p_equivalent(&self, vg: &crate::structs::GVariant, vp: &crate::structs::PVariant) -> Result<bool, HgvsError> {
        let c_variants = self.mapper.g_to_c_all(vg, self.searcher)?;
        for vc in c_variants {
            if self.c_vs_p_equivalent(&vc, vp)? {
                return Ok(true);
            }
        }
        Ok(false)
    }

    fn c_vs_p_equivalent(&self, vc: &crate::structs::CVariant, vp: &crate::structs::PVariant) -> Result<bool, HgvsError> {
        let vp_generated = self.mapper.c_to_p(vc, Some(&vp.ac))?;
        Ok(self.normalize_format(&vp_generated.to_string()) == self.normalize_format(&vp.to_string()))
    }

    fn p_vs_p_equivalent(&self, v1: &PVariant, v2: &PVariant) -> Result<bool, HgvsError> {
        if v1.ac == v2.ac {
            return Ok(self.normalize_format(&v1.to_string()) == self.normalize_format(&v2.to_string()));
        }
        Ok(false)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data::{IdentifierType, IdentifierKind, TranscriptData, ExonData, Transcript};
    use crate::coords::{GenomicPos, TranscriptPos};

    struct MockDataProvider;
    impl DataProvider for MockDataProvider {
        fn get_transcript(&self, ac: &str, _ref_ac: Option<&str>) -> Result<Box<dyn Transcript>, HgvsError> {
            if ac == "NM_000123.4" {
                Ok(Box::new(TranscriptData {
                    ac: "NM_000123.4".to_string(),
                    gene: "ABC".to_string(),
                    cds_start_index: Some(TranscriptPos(0)),
                    cds_end_index: Some(TranscriptPos(19)),
                    strand: 1,
                    reference_accession: "NC_000001.11".to_string(),
                    exons: vec![ExonData {
                        transcript_start: TranscriptPos(0),
                        transcript_end: TranscriptPos(19),
                        reference_start: GenomicPos(0),
                        reference_end: GenomicPos(19),
                        alt_strand: 1,
                        cigar: "20M".to_string(),
                    }],
                }))
            } else {
                Err(HgvsError::ValidationError("Not found".into()))
            }
        }
        fn get_seq(&self, _ac: &str, start: i32, end: i32, _kind: IdentifierType) -> Result<String, HgvsError> {
            let seq = "ACGTACGTACGTACGTACGT"; // A=0, C=1, G=2, T=3, A=4, ...
            let s = start as usize;
            let e = end as usize;
            if s < seq.len() && e <= seq.len() {
                Ok(seq[s..e].to_string())
            } else {
                Ok("".to_string())
            }
        }
        fn get_symbol_accessions(&self, _s: &str, _f: IdentifierKind, _t: IdentifierKind) -> Result<Vec<(IdentifierType, String)>, HgvsError> { Ok(vec![]) }
        fn get_identifier_type(&self, _id: &str) -> Result<IdentifierType, HgvsError> { Ok(IdentifierType::GenomicAccession) }
    }

    struct MockSearch;
    impl TranscriptSearch for MockSearch {
        fn get_transcripts_for_region(&self, _ac: &str, _s: i32, _e: i32) -> Result<Vec<String>, HgvsError> { Ok(vec![]) }
    }

    #[test]
    fn test_normalize_ins_to_dup() -> Result<(), HgvsError> {
        let hdp = MockDataProvider;
        let search = MockSearch;
        let eq = VariantEquivalence::new(&hdp, &search);

        // Genomic: NC_000001.11:g.2_3insC (base 2 index 1 is C)
        let var_g = crate::parse_hgvs_variant("NC_000001.11:g.2_3insC")?;
        let norm_g = eq.normalize_ins_to_dup(&var_g)?;
        assert_eq!(norm_g.to_string(), "NC_000001.11:g.2dupC");

        // Coding: NM_000123.4:c.2_3insC
        let var_c = crate::parse_hgvs_variant("NM_000123.4:c.2_3insC")?;
        let norm_c = eq.normalize_ins_to_dup(&var_c)?;
        assert_eq!(norm_c.to_string(), "NM_000123.4:c.2dupC");

        // NonCoding: NM_000123.4:n.2_3insC
        let var_n = crate::parse_hgvs_variant("NM_000123.4:n.2_3insC")?;
        let norm_n = eq.normalize_ins_to_dup(&var_n)?;
        assert_eq!(norm_n.to_string(), "NM_000123.4:n.2dupC");

        // Multi-base Genomic: NC_000001.11:g.4_5insACGT (bases 1-4 are ACGT)
        let var_gm = crate::parse_hgvs_variant("NC_000001.11:g.4_5insACGT")?;
        let norm_gm = eq.normalize_ins_to_dup(&var_gm)?;
        assert_eq!(norm_gm.to_string(), "NC_000001.11:g.1_4dupACGT");

        Ok(())
    }
}
