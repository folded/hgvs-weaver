# weaver

<img src="https://raw.githubusercontent.com/folded/hgvs-weaver/main/docs/source/_static/weaver.svg" alt="weaver" width=200>

High-performance HGVS variant mapping and validation engine.

Registered on PyPI as `hgvs-weaver`.
Registered on Crates.io as `hgvs-weaver`.

## Overview

`weaver` is a high-performance engine for parsing, validating, and mapping HGVS variants. It provides a robust Python interface backed by a core implementation in Rust, designed for high-throughput variant interpretation pipelines.

### Correctness through Type Safety

A key feature of `weaver` is its use of Rust's type system to ensure coordinate system integrity. Internally, the library employs "tagged integers" to represent positions in different coordinate spaces:

- **`GenomicPos`**: 0-based inclusive genomic coordinates.
- **`TranscriptPos`**: 0-based inclusive transcript coordinates (distance from the transcription start site).
- **`ProteinPos`**: 0-based inclusive amino acid positions.

The library also uses explicit types for HGVS-style 1-based coordinates:

- **`HgvsGenomicPos`**: 1-based genomic coordinates.
- **`HgvsTranscriptPos`**: 1-based cDNA/non-coding coordinates, which correctly skip the non-existent position 0 (e.g., jumps from `-1` to `1`).
- **`HgvsProteinPos`**: 1-based amino acid positions.

By enforcing these types at compile time in Rust, `weaver` prevents common off-by-one errors and accidental mixing of coordinate systems during complex mapping operations (e.g., from `g.` to `c.` to `p.`).

### Supported HGVS Features

`weaver` supports a wide range of HGVS variant types and operations:

- **Parsing**: robust parsing of `g.`, `m.`, `c.`, `n.`, `r.`, and `p.` variants.
- **Mapping**:
    - Genomic to Coding (`g.` to `c.`).
    - Coding to Genomic (`c.` to `g.`).
    - Coding to Protein (`c.` to `p.`) with full translation.
- **Normalization**: Automatic 3' shifting of variants in repetitive regions.
- **Complex Edits**: Support for deletions (`del`), insertions (`ins`), duplications (`dup`), inversions (`inv`), and repeats (`[n]`).

#### Examples

```python
import weaver

# Parsing and Formatting
v = weaver.parse("NM_000051.3:c.123A>G")
print(v.format()) # "NM_000051.3:c.123A>G"

# Mapping c. to p.
# (Requires a DataProvider, see below)
v_p = mapper.c_to_p(v)
print(v_p.format()) # "NP_000042.3:p.(Lys41Arg)"

# Normalization (3' shifting)
v_raw = weaver.parse("NM_000051.3:c.4_5del")
v_norm = mapper.normalize_variant(v_raw)
print(v_norm.format()) # e.g., "NM_000051.3:c.5_6del"
```

## Data Provider Implementation

To perform mapping operations, `weaver` requires an object that implements the `DataProvider` protocol. This object is responsible for providing transcript models and reference sequences.

### Coordinate Expectations

When implementing a `DataProvider`, you must provide coordinates in the following formats:

- **Transcript Models**:
    - `cds_start_index`: The 0-based inclusive index of the first base of the start codon (A of ATG) relative to the transcript start.
    - `cds_end_index`: The 0-based inclusive index of the last base of the stop codon relative to the transcript start.
    - **Exons**:
        - `transcript_start`: 0-based inclusive start index in the transcript.
        - `transcript_end`: 0-based **exclusive** end index in the transcript.
        - `reference_start`: 0-based inclusive start index on the genomic reference.
        - `reference_end`: 0-based inclusive end index on the genomic reference.

- **Sequence Retrieval**:
    - `get_seq(ac, start, end, kind)`: Should return the sequence for accession `ac`. `start` and `end` are 0-based half-open (interbase) coordinates.

### Python Protocol

```python
class DataProvider(Protocol):
    def get_transcript(self, transcript_ac: str, reference_ac: str | None) -> TranscriptData:
        """Return a dictionary matching the TranscriptData structure."""
        ...

    def get_seq(self, ac: str, start: int, end: int, kind: str | IdentifierType) -> str:
        """Fetch sequence for an accession. kind is an IdentifierType."""
        ...

    def get_symbol_accessions(self, symbol: str, source_kind: str, target_kind: str) -> list[tuple[str, str]] | list[tuple[IdentifierType, str]]:
        """Map gene symbols to accessions (e.g., 'ATM' -> [('transcript_accession', 'NM_000051.3')])."""
        ...

    def get_identifier_type(self, identifier: str) -> str | IdentifierType:
        """Identify what type of identifier a string is (e.g., 'genomic_accession', 'gene_symbol')."""
        ...
```

## Dataset

This repository includes a dataset of 100,000 variants sampled from ClinVar (August 2025 release) for validation purposes, located in `data/clinvar_variants_100k.tsv`.

**ClinVar License & Terms**:
ClinVar data is public domain and available for use under the terms of the [National Library of Medicine (NLM)](https://www.ncbi.nlm.nih.gov/home/about/policies/). Use of ClinVar data must adhere to their [citation and data use policies](https://www.ncbi.nlm.nih.gov/clinvar/docs/maintenance_use/).

## Installation

```sh
pip install hgvs-weaver
```

## Usage

```python
import weaver

# Parse a variant
var = weaver.parse("NM_000051.3:c.123A>G")
print(var.ac)  # NM_000051.3
print(var.format())  # NM_000051.3:c.123A>G
```

## Validation

`weaver` has been extensively validated against ClinVar data to ensure accuracy and parity with the standard Python HGVS implementation.

### Running Validation

To rerun the validation, you need the RefSeq annotation and genomic sequence files:

1. **Download Required Files**:

   ```sh
   # Download RefSeq GFF
   curl -O https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.gff.gz

   # Download RefSeq FASTA and decompress
   curl -O https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz
   gunzip GRCh38_latest_genomic.fna.gz
   ```

2. **Install Validation Dependencies**:

   ```sh
   pip install pysam tqdm bioutils parsley
   pip install hgvs --no-deps  # Avoids psycopg2 build requirement
   ```

3. **Run Validation**:
   You can run the validation using the installed entry point (if you installed with `[validation]` extra):

   ```sh
   weaver-validate data/clinvar_variants_100k.tsv \
       --output-file results.tsv \
       --gff GRCh38_latest_genomic.gff.gz \
       --fasta GRCh38_latest_genomic.fna
   ```

   Alternatively, if you use `uv`, you can run the script directly from the source tree without manually installing dependencies (it will use the PEP 723 metadata to auto-install them):

   ```sh
   uv run weaver/cli/validate.py data/clinvar_variants_100k.tsv ...
   ```

### Validation Results (100,000 variants)

Summary of results comparing `weaver` and `ref-hgvs` against ClinVar ground truth:

| Implementation | Protein Match | SPDI Match  | Parse Errors |
| :------------- | :-----------: | :---------: | :----------: |
| weaver         |  **98.590%**  | **97.968%** |    **1**     |
| ref-hgvs       |    98.116%    |   94.022%   |     394      |

RefSeq Data Mismatches: 0 (0.0%)

#### Protein Translation Agreement

|                     | ref-hgvs Match | ref-hgvs Mismatch |
| :------------------ | :------------: | :---------------: |
| **weaver Match**    |     98,103     |        487        |
| **weaver Mismatch** |       13       |       1,397       |

#### SPDI Mapping Agreement

|                     | ref-hgvs Match | ref-hgvs Mismatch |
| :------------------ | :------------: | :---------------: |
| **weaver Match**    |     93,919     |       4,049       |
| **weaver Mismatch** |      103       |       1,929       |

- **Variant Equivalence**: Check if two variants are biologically equivalent using advanced cross-coordinate mapping and normalization. [See Algorithm](docs/source/equivalence_logic.md).
