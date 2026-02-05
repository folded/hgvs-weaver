"""weaver: High-performance HGVS variant mapping and validation engine."""

from typing import Protocol, TypedDict

from ._weaver import (
    Variant as Variant,
)
from ._weaver import (
    VariantMapper as VariantMapper,
)
from ._weaver import (
    parse as parse,
)


class ExonData(TypedDict):
    """Represents an exon's coordinates and alignment.

    Coordinates are 0-based:
    - transcript_start: inclusive start index in transcript.
    - transcript_end: exclusive end index in transcript.
    - reference_start: inclusive start index on genomic reference.
    - reference_end: inclusive end index on genomic reference.
    """

    transcript_start: int
    transcript_end: int
    reference_start: int
    reference_end: int
    alt_strand: int  # 1 for plus, -1 for minus
    cigar: str  # Extended CIGAR string (e.g., "100=")


class TranscriptData(TypedDict):
    """Represents a full transcript model.

    Coordinates are 0-based:
    - cds_start_index: inclusive index of the first base of the start codon.
    - cds_end_index: inclusive index of the last base of the stop codon.
    """

    ac: str
    gene: str
    cds_start_index: int | None
    cds_end_index: int | None
    strand: int  # 1 or -1
    reference_accession: str  # Genomic accession (e.g., NC_000001.11)
    exons: list[ExonData]


class TranscriptSearch(Protocol):
    """Optional interface for regional discovery."""

    def get_transcripts_for_region(self, chrom: str, start: int, end: int) -> list[str]:
        """Return list of transcript accessions overlapping the given genomic region."""
        ...


class DataProvider(Protocol):
    """Required interface for the object passed to VariantMapper."""

    def get_transcript(self, transcript_ac: str, reference_ac: str | None) -> TranscriptData:
        """Retrieve transcript model for the given accession.

        If reference_ac is provided, returns the alignment for that specific reference.
        """
        ...

    def get_seq(self, ac: str, start: int, end: int, kind: str) -> str:
        """Fetch sequence for an accession.

        - kind: one of 'g' (genomic), 'c' (transcript), 'p' (protein).
        - start/end: 0-based half-open (interbase) coordinates.
        - end=-1: fetch to the end of the sequence.
        """
        ...

    def get_symbol_accessions(self, symbol: str, source_kind: str, target_kind: str) -> list[str]:
        """Map identifiers between different namespaces.

        Usually used to map gene symbols to transcript accessions.
        """
        ...

    def get_identifier_type(self, identifier: str) -> str:
        """Identify what type of identifier a string is.

        Should return one of:
        - 'genomic_accession'
        - 'transcript_accession'
        - 'protein_accession'
        - 'gene_symbol'
        - 'unknown'
        """
        ...
