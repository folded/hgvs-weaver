/// Error types for HGVS operations.
#[derive(thiserror::Error, Debug)]
pub enum HgvsError {
    /// Failure during parsing of an HGVS string or CIGAR string.
    #[error("Parse error: {0}")]
    PestError(String),
    /// Failure during validation of transcript or exon metadata.
    #[error("Validation error: {0}")]
    ValidationError(String),
    /// Failure during data retrieval from a `DataProvider`.
    #[error("Data provider error: {0}")]
    DataProviderError(String),
    /// Attempted an operation that is not yet supported.
    #[error("Unsupported operation: {0}")]
    UnsupportedOperation(String),
    /// Error specifically related to CIGAR string mapping.
    #[error("CIGAR error: {0}")]
    CigarError(String),
    /// Catch-all for other error types.
    #[error("Other error: {0}")]
    Other(String),
}
