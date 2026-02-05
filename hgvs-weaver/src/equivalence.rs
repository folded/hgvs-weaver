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

    pub fn are_equivalent(&self, var1: &SequenceVariant, var2: &SequenceVariant) -> Result<bool, HgvsError> {
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

    fn expand_if_gene_symbol(&self, var: &SequenceVariant) -> Result<Vec<SequenceVariant>, HgvsError> {
        let ac = var.ac();
        if !ac.contains('.') && !ac.starts_with("NC_") && !ac.starts_with("NM_") && !ac.starts_with("NP_") && !ac.starts_with("NR_") {
            // Likely a gene symbol
            let target_kind = match var {
                SequenceVariant::Protein(_) => IdentifierKind::Protein,
                _ => IdentifierKind::Transcript,
            };
            let accessions = self.hdp.get_symbol_accessions(ac, IdentifierKind::Genomic, target_kind).unwrap_or_default();
            if accessions.is_empty() {
                return Ok(vec![var.clone()]);
            }
            
            let mut expanded = Vec::new();
            for new_ac in accessions {
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
            Ok(expanded)
        } else {
            Ok(vec![var.clone()])
        }
    }

    fn are_equivalent_single(&self, var1: &SequenceVariant, var2: &SequenceVariant) -> Result<bool, HgvsError> {
        match (var1, var2) {
            // Nucleotide vs Nucleotide
            (SequenceVariant::Genomic(v1), SequenceVariant::Genomic(v2)) => self.n_vs_n_equivalent(v1, v2),
            (SequenceVariant::Coding(v1), SequenceVariant::Coding(v2)) => self.n_vs_n_equivalent_c(v1, v2),
            (SequenceVariant::Genomic(v1), SequenceVariant::Coding(v2)) => self.g_vs_c_equivalent(v1, v2),
            (SequenceVariant::Coding(v1), SequenceVariant::Genomic(v2)) => self.g_vs_c_equivalent(v2, v1),

            // Nucleotide vs Protein
            (SequenceVariant::Genomic(v1), SequenceVariant::Protein(v2)) => self.g_vs_p_equivalent(v1, v2),
            (SequenceVariant::Protein(v1), SequenceVariant::Genomic(v2)) => self.g_vs_p_equivalent(v2, v1),
            (SequenceVariant::Coding(v1), SequenceVariant::Protein(v2)) => self.c_vs_p_equivalent(v1, v2),
            (SequenceVariant::Protein(v1), SequenceVariant::Coding(v2)) => self.c_vs_p_equivalent(v2, v1),

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
        s.replace(['(', ')', '?'], "")
    }

    fn n_vs_n_equivalent(&self, v1: &GVariant, v2: &GVariant) -> Result<bool, HgvsError> {
        let nv1 = self.mapper.normalize_variant(SequenceVariant::Genomic(v1.clone()))?;
        let nv2 = self.mapper.normalize_variant(SequenceVariant::Genomic(v2.clone()))?;
        Ok(self.normalize_format(&nv1.to_string()) == self.normalize_format(&nv2.to_string()))
    }

    fn n_vs_n_equivalent_c(&self, v1: &CVariant, v2: &CVariant) -> Result<bool, HgvsError> {
        let g1 = self.mapper.c_to_g(v1, &self.hdp.get_transcript(&v1.ac, None)?.reference_accession())?;
        let g2 = self.mapper.c_to_g(v2, &self.hdp.get_transcript(&v2.ac, None)?.reference_accession())?;
        self.n_vs_n_equivalent(&g1, &g2)
    }

    fn g_vs_c_equivalent(&self, vg: &GVariant, vc: &CVariant) -> Result<bool, HgvsError> {
        let g2 = self.mapper.c_to_g(vc, &vg.ac)?;
        self.n_vs_n_equivalent(vg, &g2)
    }

    fn g_vs_p_equivalent(&self, vg: &GVariant, vp: &PVariant) -> Result<bool, HgvsError> {
        let c_variants = self.mapper.g_to_c_all(vg, self.searcher)?;
        for vc in c_variants {
            if self.c_vs_p_equivalent(&vc, vp)? {
                return Ok(true);
            }
        }
        Ok(false)
    }

    fn c_vs_p_equivalent(&self, vc: &CVariant, vp: &PVariant) -> Result<bool, HgvsError> {
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
