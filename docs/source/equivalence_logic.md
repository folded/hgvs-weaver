# Variant Equivalence Algorithm

This document outlines the algorithm used by `hgvs-weaver` to determine if two variants are biologically equivalent. The `VariantEquivalence` struct manages this process, handling gene symbol expansion, normalization, and cross-coordinate mapping.

## Overview

The equivalence check proceeds in two main phases:

1. **Gene Symbol Expansion**: If a variant uses a gene symbol (e.g., `BRAF:c.1799T>A`) instead of an accession, it is resolved to all valid accessions for that gene.
2. **Comparison**: All combinations of the expanded variants are compared using strategies specific to their coordinate types.

## Algorithm Flowchart

```mermaid
flowchart TD
    Start([Start: Compare var1, var2]) --> Expansion1[Expand var1]
    Expansion1 --> Expansion2[Expand var2]
    Expansion2 --> LoopStart{For each v1 in vars1\nFor each v2 in vars2}
    
    LoopStart -->|Next Pair| CheckType{Check Types}
    
    CheckType -- Both Genomic (g.) --> NvsN[Normalize & String Compare]
    CheckType -- Both Coding (c.) --> CvsC[Map c. to g. \nNormalize & String Compare]
    CheckType -- g. vs c. --> GvsC[Map c. to g. \nNormalize & String Compare]
    
    CheckType -- g. vs p. --> GvsP[Map g. to all c. (transcripts)\nProject c. to p.\nCompare p. strings]
    CheckType -- c. vs p. --> CvsP[Project c. to p.\nCompare p. strings]
    CheckType -- Both Protein (p.) --> PvsP[Direct String Compare]
    
    CheckType -- Other/Mismatch --> Fallback[Simple String Compare]
    
    NvsN --> IsMatch{Match?}
    CvsC --> IsMatch
    GvsC --> IsMatch
    GvsP --> IsMatch
    CvsP --> IsMatch
    PvsP --> IsMatch
    Fallback --> IsMatch
    
    IsMatch -- Yes --> ReturnTrue([Return TRUE])
    IsMatch -- No --> LoopStart
    
    LoopStart -- No More Pairs --> ReturnFalse([Return FALSE])

    subgraph SymbolExpansion [Gene Symbol Expansion]
        direction TB
        Input[Variant] --> IsSymbol{Is Gene Symbol?}
        IsSymbol -- Yes --> GetKind{Target Kind?}
        GetKind -- Protein --> FetchP[Fetch Protein Accessions]
        GetKind -- Transcript --> FetchT[Fetch Transcript Accessions]
        GetKind -- Genomic --> FetchG[Fetch Genomic Accessions]
        
        FetchP --> Filter[Filter by Compatibility]
        FetchT --> Filter
        FetchG --> Filter
        
        Filter --> Output[List of Accession-based Variants]
        IsSymbol -- No --> AsIs[Return Original Variant]
    end
```

## Detailed Steps

### 1. Gene Symbol Expansion

Before comparison, each variant is checked to see if it uses a gene symbol (e.g., `BRAF`) via `DataProvider::get_identifier_type`.

- **Smart Selection**: The expansion targets a specific accession type based on the variant's coordinate system:
    - `p.` variants $\rightarrow$ `ProteinAccession`
    - `c.` / `n.` / `r.` variants $\rightarrow$ `TranscriptAccession`
    - `g.` / `m.` variants $\rightarrow$ `GenomicAccession`
- **Compatibility Filter**: Returned accessions are filtered to ensure they match the variant type (e.g., ensuring a `c.` variant doesn't get assigned a protein accession).

### 2. Comparison Strategies

Code: [hgvs-weaver/src/equivalence.rs](file:///Users/tjs/repo/hgvs-rs/hgvs-weaver/src/equivalence.rs)

#### Genomic vs. Genomic (`n_vs_n_equivalent`)

1. Both variants are normalized to their 3'-most representation using `VariantMapper::normalize_variant`.
2. The normalized strings are compared (ignoring parentheses/uncertainties).

#### Coding vs. Coding (`n_vs_n_equivalent_c`)

1. Both `c.` variants are mapped to genomic coordinates (`c_to_g`) using their respective transcript references.
2. The resulting `g.` variants are compared using the **Genomic vs. Genomic** strategy.

#### Genomic vs. Coding (`g_vs_c_equivalent`)

1. The `c.` variant is mapped to `g.` coordinates *on the reference sequence of the genomic variant*.
2. The resulting `g.` variants are compared.

#### Genomic vs. Protein (`g_vs_p_equivalent`)

1. The `g.` variant is mapped to **all** overlapping transcripts (`g_to_c_all` via `TranscriptSearch`).
2. For each discovered `c.` variant, it is projected to protein (`c_to_p`).
3. The projected `p.` variant is compared to the target `p.` variant.
4. If *any* path results in a match, the variants are equivalent.

#### Coding vs. Protein (`c_vs_p_equivalent`)

1. The `c.` variant is projected to protein (`c_to_p`) using the accession of the target `p.` variant.
2. The projected `p.` string is compared to the target `p.` string.
