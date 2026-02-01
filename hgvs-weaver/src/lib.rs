use pest_derive::Parser as PestParser;
use pest::Parser;

#[derive(PestParser)]
#[grammar = "grammar.pest"]
pub struct HgvsParser;

/// Parses an HGVS string into a `SequenceVariant`.
pub fn parse_hgvs_variant(hgvs_str: &str) -> Result<SequenceVariant, HgvsError> {
    let mut pairs = HgvsParser::parse(Rule::hgvs_variant, hgvs_str)
        .map_err(|e| HgvsError::PestError(e.to_string()))?;

    let pair = pairs.next().ok_or_else(|| HgvsError::PestError("Empty input".into()))?;
    let inner = pair.into_inner().next().ok_or_else(|| HgvsError::PestError("Missing inner variant".into()))?;

    match inner.as_rule() {
        Rule::g_variant => {
            let mut inner = inner.into_inner();
            let ac = inner.next().ok_or_else(|| HgvsError::PestError("Missing accession".into()))?.as_str().to_string();
            let gene_expr_pair = inner.next().ok_or_else(|| HgvsError::PestError("Missing gene expr".into()))?;
            let gene = parse_gene_expr(gene_expr_pair);
            let posedit = parser::parse_g_posedit(inner.next().ok_or_else(|| HgvsError::PestError("Missing posedit".into()))?)?;
            Ok(SequenceVariant::Genomic(GVariant { ac, gene, posedit }))
        }
        Rule::c_variant => {
            let mut inner = inner.into_inner();
            let ac = inner.next().ok_or_else(|| HgvsError::PestError("Missing accession".into()))?.as_str().to_string();
            let gene_expr_pair = inner.next().ok_or_else(|| HgvsError::PestError("Missing gene expr".into()))?;
            let gene = parse_gene_expr(gene_expr_pair);
            let posedit = parser::parse_c_posedit(inner.next().ok_or_else(|| HgvsError::PestError("Missing posedit".into()))?)?;
            Ok(SequenceVariant::Coding(CVariant { ac, gene, posedit }))
        }
        Rule::p_variant => {
            let mut inner = inner.into_inner();
            let ac = inner.next().ok_or_else(|| HgvsError::PestError("Missing accession".into()))?.as_str().to_string();
            let gene_expr_pair = inner.next().ok_or_else(|| HgvsError::PestError("Missing gene expr".into()))?;
            let gene = parse_gene_expr(gene_expr_pair);
            let posedit = parser::parse_p_posedit(inner.next().ok_or_else(|| HgvsError::PestError("Missing posedit".into()))?)?;
            Ok(SequenceVariant::Protein(PVariant { ac, gene, posedit }))
        }
        Rule::m_variant => {
            let mut inner = inner.into_inner();
            let ac = inner.next().ok_or_else(|| HgvsError::PestError("Missing accession".into()))?.as_str().to_string();
            let gene_expr_pair = inner.next().ok_or_else(|| HgvsError::PestError("Missing gene expr".into()))?;
            let gene = parse_gene_expr(gene_expr_pair);
            let posedit = parser::parse_g_posedit(inner.next().ok_or_else(|| HgvsError::PestError("Missing posedit".into()))?)?;
            Ok(SequenceVariant::Mitochondrial(MVariant { ac, gene, posedit }))
        }
        Rule::n_variant => {
            let mut inner = inner.into_inner();
            let ac = inner.next().ok_or_else(|| HgvsError::PestError("Missing accession".into()))?.as_str().to_string();
            let gene_expr_pair = inner.next().ok_or_else(|| HgvsError::PestError("Missing gene expr".into()))?;
            let gene = parse_gene_expr(gene_expr_pair);
            let posedit = parser::parse_c_posedit(inner.next().ok_or_else(|| HgvsError::PestError("Missing posedit".into()))?)?;
            Ok(SequenceVariant::NonCoding(NVariant { ac, gene, posedit }))
        }
        Rule::r_variant => {
            let mut inner = inner.into_inner();
            let ac = inner.next().ok_or_else(|| HgvsError::PestError("Missing accession".into()))?.as_str().to_string();
            let gene_expr_pair = inner.next().ok_or_else(|| HgvsError::PestError("Missing gene expr".into()))?;
            let gene = parse_gene_expr(gene_expr_pair);
            let posedit = parser::parse_c_posedit(inner.next().ok_or_else(|| HgvsError::PestError("Missing posedit".into()))?)?;
            Ok(SequenceVariant::Rna(RVariant { ac, gene, posedit }))
        }
        _ => Err(HgvsError::PestError("Unsupported variant type".into())),
    }
}

fn parse_gene_expr(pair: pest::iterators::Pair<Rule>) -> Option<String> {
    let s = pair.as_str();
    if s.is_empty() { return None; }
    Some(s.replace(['(', ')'], ""))
}

pub mod structs;
pub mod coords;
pub mod parser;
pub mod mapper;
pub mod error;
pub mod transcript_mapper;
pub mod cigar;
pub mod data;
pub mod utils;
pub mod edits;
pub mod fmt;
pub mod altseq;
pub mod altseq_to_hgvsp;

// Re-exports for public usage
pub use structs::{GVariant, CVariant, PVariant, MVariant, NVariant, RVariant, Variant};
pub use coords::{SequenceVariant};
pub use data::{DataProvider, TranscriptSearch, Transcript, IdentifierKind};
pub use mapper::VariantMapper;
pub use error::HgvsError;
