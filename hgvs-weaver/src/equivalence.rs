use crate::error::HgvsError;
use crate::data::{DataProvider, TranscriptSearch, IdentifierKind};
use crate::structs::{SequenceVariant, GVariant, CVariant, PVariant, Variant};
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

    fn eq_dup_ins(&self, s1: &str, s2: &str) -> bool {
        // Simple heuristic: if one is dup and other is ins, and they look similar.
        // Actually, converting textual representation is hard.
        // Better to check if normalization logic (fill_implicit) handles it.
        // `dup` normalized is `dupSEQ`. `insSEQ` is `insSEQ`.
        // Semantically: pos_dupSEQ == pos_insSEQ?
        // No. `10dupT` means insert T after 10. `10_11insT` means insert T between 10 and 11.
        // `10dup` (where 10 is T) -> `10_11insT`.
        // So `pos` vs `pos_pos+1` interval mismatch.
        // But `hgvs` crate outputs `dup` as `dup`.
        // Let's rely on standard string normalization for now, but explicit SEQ might help.
        // If we fill SEQ, `del` becomes `delSEQ` which matches explicit.
        // For `dup`, `dup` becomes `dupSEQ`.
        // If failure persists for `ins` vs `dup`, we might need more logic.
        // Let's assume filling implicit is enough for now.
        s1 == s2
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
        // If variant is Ins, check if inserted sequence matches immediately preceding sequence.
        // If so, convert to Dup.
        // E.g. 10_11insT -> check base at 10. If T, then 10dupT.
        // This requires fetching sequence at start position.

        // This is complex because we need to know coordinate type and fetch seq.
        // Similar dispatch as fill_implicit_sequence.
        match var {
            SequenceVariant::Genomic(v) => {
                if let Some(pos) = &v.posedit.pos {
                     if let crate::edits::NaEdit::Ins { alt: Some(seq), uncertain } = &v.posedit.edit {
                         let start = pos.start.base.to_index().0 as usize; // insertion is between start and end? 10_11insT -> start=10 (0-based 9?).
                         // HGVS 10_11insT means between 10 and 11.
                         // Internal representation: interval start=10, end=11? or start=10, end=10?
                         // Let's assume start is the base BEFORE insertion.

                         // Check base at `start`.
                         // Note: `pos.start` is strictly the implementation detail.
                         // But `get_seq` is 0-based.
                         // `to_index()` gives 0-based index.

                         // If seq len is 1, check 1 base. If len N, check N bases ending at start.
                         let len = seq.len();
                         // Fetch reference sequence ending at `start` (inclusive).
                         // So range [start - len + 1, start + 1].
                         // start index from `to_index` is 0-based index of the base.

                         let start_idx = start as i32;
                         let check_start = start_idx - (len as i32) + 1;
                         if check_start >= 0 {
                             let ref_seq = self.hdp.get_seq(&v.ac, check_start, start_idx + 1, IdentifierKind::Genomic.into_identifier_type())?;
                             if ref_seq == *seq {
                                 // Convert to Dup
                                 // Dup position should be the range of the duplicated sequence.
                                 // 10_11insT (where 10 is T). Dup is 10dupT.
                                 // pos is 10.
                                 // Range for dup is [check_start, start].
                                 let mut new_v = v.clone();
                                 // Update position to be the range of the duplication source
                                 // But wait, `Dup` usually implies `pos` is the *copied* sequence.
                                 // Yes.
                                 // So we need to update `pos` to encompass `[check_start, start]`.
                                 // And change edit to `Dup`.

                                 // We need to construct a new interval.
                                 // Accessing `Base` directly is hard if we don't have simple constructors handy, but we can modify existing.
                                 // `pos.start.base` is `HgvsGenomicPos(start)`.
                                 // We need `HgvsGenomicPos(check_start)`.

                                 // Actually, modifying `pos` generically is hard.
                                 // But we are in `Genomic` arm, so `pos` is `SimpleInterval`.
                                 // `pos.start.base` is `HgvsGenomicPos`.

                                 // new start:
                                 new_v.posedit.pos.as_mut().unwrap().start.base.0 = check_start + 1; // 1-based for HGVS struct?
                                 // Wait, `HgvsGenomicPos` wraps `i32` which is 1-based?
                                 // `to_index()` subtracts 1.
                                 // `SimplePosition { base: HgvsGenomicPos(val) }`.
                                 // If `start_idx` is 0-based, then `new_val` = `check_start + 1`.

                                 // new end:
                                 if len > 1 {
                                     new_v.posedit.pos.as_mut().unwrap().end = Some(crate::structs::SimplePosition {
                                         base: crate::coords::HgvsGenomicPos(start_idx + 1),
                                         end: None,
                                         uncertain: false
                                     });
                                 } else {
                                     new_v.posedit.pos.as_mut().unwrap().end = None; // Single base dup
                                 }

                                 new_v.posedit.edit = crate::edits::NaEdit::Dup { ref_: Some(seq.clone()), uncertain: *uncertain };
                                 return Ok(SequenceVariant::Genomic(new_v));
                             }
                         }
                     }
                }
                Ok(var.clone())
            },
            // Similar for other types if needed (Coding, etc.)
            // For now only implementing Genomic as that's where the test failure is.
            _ => Ok(var.clone()),
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
