# hgvs-weaver Code Review

> Prepared 2026-02-24. Covers all Rust source files in `hgvs-weaver/src/`.
> Updated 2026-02-25 to reflect fixes applied in subsequent commits.

---

## Open findings

| # | Severity | File | Description |
|---|----------|------|-------------|
| 5 | Major | `structs.rs` / `equivalence.rs` | `DataProvider::c_to_g` contract mismatch (CDS-relative vs transcript-relative) |
| 12 | Minor | `structs.rs` | `SimplePosition.end` for uncertain positions not yet parsed/formatted |

---

## Detailed findings

### 5. `DataProvider::c_to_g` interface contract mismatch (Major)

**File:** `hgvs-weaver/src/structs.rs` (impl `IntervalSpdi for BaseOffsetInterval`), `equivalence.rs`

`DataProvider::c_to_g` takes a `TranscriptPos` documented as "CDS-relative, 0-based", but the callers in `structs.rs` (line 162) and `equivalence.rs` pass the result of `HgvsTranscriptPos::to_index()` which is a CDS-relative index adjusted for the ±1 offset, but not shifted by `cds_start`.

This means the SPDI coordinates produced by `BaseOffsetInterval::spdi_interval` for c. variants will be wrong for any transcript where `cds_start ≠ 0`. The `equivalence.rs` comment at line 173 already hints at this: `"keep consistent with spdi_interval which is currently 1-based"`.

**Recommended fix:** Clarify the contract of `DataProvider::c_to_g`. If it expects transcript-relative coords, callers must add `cds_start_index()` before calling. If it expects CDS-relative coords, the doc comment and type should reflect that. The `TranscriptMapper` already has `c_to_n` which correctly handles the offset — consider routing through it.

---

### 12. `SimplePosition.end` not yet parsed/formatted (Minor)

**File:** `hgvs-weaver/src/structs.rs`

`SimplePosition.end` represents the uncertainty window for positions like `(1_3)_(7_10)`, where the exact breakpoint is unknown. The field is intentional and correctly models this HGVS notation.

**Remaining work:** The parser and formatter do not yet round-trip uncertain genomic positions. Until they do, `end` will remain `None` in practice.

---

## Type system opportunities

### A. Replace `strand: i32` with a `Strand` enum

`Exon::alt_strand` and `Transcript::strand()` return `i32` where only `1` and `-1` are meaningful. Every consumer checks `== 1` or `== -1` with no exhaustiveness guarantee. A `Strand { Plus, Minus }` enum would make invalid strand values unrepresentable and enable exhaustive matching.

### C. Anchor consistency enforcement

`BaseOffsetPosition` can be constructed with any `Anchor` regardless of the enclosing variant type (e.g., a `NVariant` with `Anchor::CdsStart`). Consider a builder or `From` impl that enforces the correct anchor per coordinate system.

---

## Test coverage gaps

The following areas have limited or no test coverage:

| File | What is still untested |
|------|------------------------|
| `transcript_mapper.rs` | `n_to_c`, `c_to_n` |
| `altseq.rs` | `AltSeqBuilder::build_altseq` for most edit types |
| `altseq_to_hgvsp.rs` | Most protein consequence cases (frameshift, delins, synonymous) |

`transcript_mapper.rs` now has unit tests for `g_to_n` (exonic, intronic, CIGAR, and minus-strand ordering). `mapper.rs` is covered by `mapping_test.rs` and `test_shift.rs`, including multi-base insertion shifting.
