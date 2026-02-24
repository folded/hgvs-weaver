use crate::coords::SequenceVariant;
use crate::edits::AaEdit;
use crate::structs::{AaInterval, AAPosition, PVariant, PosEdit};

/// Controls how start-codon variants are represented in protein notation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StartCodonConvention {
    /// Keep the specific predicted amino acid change (e.g., `p.(Met1Val)`). Default.
    Specific,
    /// Use the HGVS `p.Met1?` notation for any non-silent change at the first position.
    HgvsQuestion,
}

/// Settings that control how a variant is transformed before formatting or comparison.
#[derive(Debug, Clone)]
pub struct VariantTransformSettings {
    /// Convention to use for start-codon protein variants.
    pub start_codon: StartCodonConvention,
}

impl Default for VariantTransformSettings {
    fn default() -> Self {
        VariantTransformSettings {
            start_codon: StartCodonConvention::Specific,
        }
    }
}

impl VariantTransformSettings {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn with_start_codon(mut self, convention: StartCodonConvention) -> Self {
        self.start_codon = convention;
        self
    }
}

/// Applies transform settings to a `SequenceVariant`, returning a new variant.
///
/// Currently only protein variants are transformed (e.g., start-codon convention).
/// All other variant types are returned unchanged.
pub fn transform_variant(
    var: &SequenceVariant,
    settings: &VariantTransformSettings,
) -> SequenceVariant {
    match var {
        SequenceVariant::Protein(vp) => {
            SequenceVariant::Protein(transform_protein(vp, settings))
        }
        _ => var.clone(),
    }
}

fn transform_protein(vp: &PVariant, settings: &VariantTransformSettings) -> PVariant {
    match settings.start_codon {
        StartCodonConvention::Specific => vp.clone(),
        StartCodonConvention::HgvsQuestion => apply_hgvs_question(vp),
    }
}

/// Converts a start-codon protein variant to the `p.Met1?` notation if applicable.
///
/// Applies when:
/// - The variant is at position 1 (single amino acid, no range end)
/// - The reference amino acid is Met (any case, 1- or 3-letter)
/// - The edit is not silent (identity)
fn apply_hgvs_question(vp: &PVariant) -> PVariant {
    if let Some(pos) = &vp.posedit.pos {
        let is_pos1 = pos.start.base.0 == 1 && pos.end.is_none();
        let aa = pos.start.aa.to_lowercase();
        let is_met = aa == "met" || aa == "m";
        let is_silent = vp.posedit.edit.is_identity();

        if is_pos1 && is_met && !is_silent {
            return PVariant {
                ac: vp.ac.clone(),
                gene: vp.gene.clone(),
                posedit: PosEdit {
                    pos: Some(AaInterval {
                        start: AAPosition {
                            base: pos.start.base,
                            aa: pos.start.aa.clone(),
                            uncertain: false,
                        },
                        end: None,
                        uncertain: false,
                    }),
                    edit: AaEdit::Special {
                        value: "?".to_string(),
                        uncertain: false,
                    },
                    uncertain: false,
                    predicted: false,
                },
            };
        }
    }
    vp.clone()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parse_hgvs_variant;

    #[test]
    fn test_transform_met1_to_question() {
        let settings = VariantTransformSettings {
            start_codon: StartCodonConvention::HgvsQuestion,
        };

        // p.(Met1Val) → p.Met1?
        let var = parse_hgvs_variant("NP_000051.2:p.(Met1Val)").unwrap();
        let transformed = transform_variant(&var, &settings);
        assert_eq!(transformed.to_string(), "NP_000051.2:p.Met1?");
    }

    #[test]
    fn test_transform_specific_unchanged() {
        let settings = VariantTransformSettings {
            start_codon: StartCodonConvention::Specific,
        };

        // Specific: p.(Met1Val) stays as is
        let var = parse_hgvs_variant("NP_000051.2:p.(Met1Val)").unwrap();
        let transformed = transform_variant(&var, &settings);
        assert_eq!(transformed.to_string(), "NP_000051.2:p.(Met1Val)");
    }

    #[test]
    fn test_transform_non_met1_unchanged() {
        let settings = VariantTransformSettings {
            start_codon: StartCodonConvention::HgvsQuestion,
        };

        // Not position 1: p.(Gly2Arg) → unchanged
        let var = parse_hgvs_variant("NP_000051.2:p.(Gly2Arg)").unwrap();
        let transformed = transform_variant(&var, &settings);
        assert_eq!(transformed.to_string(), "NP_000051.2:p.(Gly2Arg)");
    }

    #[test]
    fn test_transform_silent_met1_unchanged() {
        let settings = VariantTransformSettings {
            start_codon: StartCodonConvention::HgvsQuestion,
        };

        // Silent: p.(Met1=) → unchanged (identity)
        let var = parse_hgvs_variant("NP_000051.2:p.(Met1=)").unwrap();
        let transformed = transform_variant(&var, &settings);
        assert_eq!(transformed.to_string(), "NP_000051.2:p.(Met1=)");
    }

    #[test]
    fn test_transform_non_protein_unchanged() {
        let settings = VariantTransformSettings {
            start_codon: StartCodonConvention::HgvsQuestion,
        };

        // Non-protein variant: unchanged
        let var = parse_hgvs_variant("NM_000051.3:c.1A>T").unwrap();
        let transformed = transform_variant(&var, &settings);
        assert_eq!(transformed.to_string(), "NM_000051.3:c.1A>T");
    }
}
