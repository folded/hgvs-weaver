# HGVS Coordinate Systems

This document summarizes the coordinate numbering rules for various reference sequences as defined by the HGVS nomenclature.

## General Principles

- **1-based Indexing**: All coordinate systems are 1-based. The first nucleotide or amino acid is 1.
- **No Position 0**: There is no position 0. Numbering jumps from -1 to 1.
- **Sequential Numbering**: Genomic (`g.`), Mitochondrial (`m.`), and Circular (`o.`) DNA use simple sequential numbering starting from the first nucleotide of the reference sequence.

---

## Coding DNA (`c.`)

The coding DNA coordinate system is relative to the transcript's coding sequence (CDS).

- **Start Codon**: `c.1` is the **A** of the `ATG` translation initiation codon.
- **5' UTR**: Nucleotides upstream of the start codon are numbered `c.-1`, `c.-2`, etc. (moving further upstream).
- **Stop Codon**: The last nucleotide of the translation termination codon (stop codon) is the end of the `c.` range.
- **3' UTR**: Nucleotides downstream of the stop codon are numbered `c.*1`, `c.*2`, etc. (moving further downstream). `c.*1` is the base immediately following the stop codon.
- **Introns**:
    - **5' end of intron**: Numbered relative to the last nucleotide of the preceding exon, followed by a `+` (e.g., `c.87+1`).
    - **3' end of intron**: Numbered relative to the first nucleotide of the following exon, followed by a `-` (e.g., `c.88-1`).
    - **UTR Introns**: Introns within UTRs follow the same rule, using the UTR coordinate as the anchor (e.g., `c.-85+1` or `c.*37-1`).

---

## Non-coding Transcript (`n.`)

Used for transcripts that do not code for a protein (e.g., lncRNA).

- **Sequential Numbering**: Numbered `n.1`, `n.2`, etc., from the first to the last nucleotide of the transcript.
- **No Asterisk Notation**: Unlike `c.`, `n.` does not use the `*` notation for downstream regions; it uses simple sequential numbering for the entire transcript.
- **Introns**: Numbered similarly to `c.` introns but relative to `n.` positions (e.g., `n.45+1`).

---

## RNA (`r.`)

RNA coordinates follow the associated DNA reference sequence.

- **Coding RNA**: Follows `c.` numbering (e.g., `r.123` corresponds to `c.123`).
- **Non-coding RNA**: Follows `n.` numbering (e.g., `r.45` corresponds to `n.45`).
- **Exonic Only**: RNA reference sequences do not contain introns.

---

## Protein (`p.`)

- **Initiator Methionine**: The first amino acid (usually Methionine) is `p.1`.
- **Sequential Numbering**: Amino acids are numbered sequentially `p.1`, `p.2`, etc.
- **Stop Codon**: Represented as `Ter` or `*`.

---

## Implementation Notes for `weaver`

- **Internal 0-based representation**: The library uses `GenomicPos`, `TranscriptPos`, and `ProteinPos` as internal 0-based indices.
- **Coordinate Mapping**:
    - `c.1` (1-based) <-> Index `0` (internal).
    - `c.-1` (1-based) <-> Index `-1` (internal).
    - `c.*1` (1-based) <-> Index `cds_end + 1` (internal).
- **Skip-Zero Behavior**: The conversion logic must ensure that `TranscriptPos(0)` maps to HGVS `1` and `TranscriptPos(-1)` maps to HGVS `-1` to correctly respect the "No Position 0" rule.
