use dyn_clone::DynClone;
use crate::error::HgvsError;
use crate::structs::{GenomicPos, TranscriptPos};
use serde::{Serialize, Deserialize};

pub trait Exon: DynClone {
    fn transcript_start(&self) -> TranscriptPos;
    fn transcript_end(&self) -> TranscriptPos;
    fn reference_start(&self) -> GenomicPos;
    fn reference_end(&self) -> GenomicPos;
    fn alt_strand(&self) -> i32;
    fn cigar(&self) -> &str;
}
dyn_clone::clone_trait_object!(Exon);

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ExonData {
    pub transcript_start: TranscriptPos,
    pub transcript_end: TranscriptPos,
    pub reference_start: GenomicPos,
    pub reference_end: GenomicPos,
    pub alt_strand: i32,
    pub cigar: String,
}

impl Exon for ExonData {
    fn transcript_start(&self) -> TranscriptPos { self.transcript_start }
    fn transcript_end(&self) -> TranscriptPos { self.transcript_end }
    fn reference_start(&self) -> GenomicPos { self.reference_start }
    fn reference_end(&self) -> GenomicPos { self.reference_end }
    fn alt_strand(&self) -> i32 { self.alt_strand }
    fn cigar(&self) -> &str { &self.cigar }
}

pub trait Transcript: DynClone {
    fn ac(&self) -> &str;
    fn gene(&self) -> &str;
    fn cds_start_index(&self) -> Option<TranscriptPos>;
    fn cds_end_index(&self) -> Option<TranscriptPos>;
    fn strand(&self) -> i32;
    fn reference_accession(&self) -> &str;
    fn exons(&self) -> &[ExonData];
}
dyn_clone::clone_trait_object!(Transcript);

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TranscriptData {
    pub ac: String,
    pub gene: String,
    pub cds_start_index: Option<TranscriptPos>,
    pub cds_end_index: Option<TranscriptPos>,
    pub strand: i32,
    pub reference_accession: String,
    pub exons: Vec<ExonData>,
}

impl Transcript for TranscriptData {
    fn ac(&self) -> &str { &self.ac }
    fn gene(&self) -> &str { &self.gene }
    fn cds_start_index(&self) -> Option<TranscriptPos> { self.cds_start_index }
    fn cds_end_index(&self) -> Option<TranscriptPos> { self.cds_end_index }
    fn strand(&self) -> i32 { self.strand }
    fn reference_accession(&self) -> &str { &self.reference_accession }
    fn exons(&self) -> &[ExonData] { &self.exons }
}

/// Interface for retrieving transcript and sequence data.
pub trait DataProvider {
    fn get_transcript(&self, transcript_ac: &str, reference_ac: Option<&str>) -> Result<Box<dyn Transcript>, HgvsError>;
    fn get_seq(&self, ac: &str, start: i32, end: i32, kind: IdentifierKind) -> Result<String, HgvsError>;
    fn get_symbol_accessions(&self, symbol: &str, source_kind: IdentifierKind, target_kind: IdentifierKind) -> Result<Vec<String>, HgvsError>;
}

/// Interface for discovering transcripts by region.
pub trait TranscriptSearch {
    fn get_transcripts_for_region(&self, chrom: &str, start: i32, end: i32) -> Result<Vec<String>, HgvsError>;
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum IdentifierKind {
    Genomic,
    Transcript,
    Protein,
}
