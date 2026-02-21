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

<!-- PERFORMANCE_GRAPH_START -->

<p align="center">
<svg xmlns:xlink="http://www.w3.org/1999/xlink" width="836.594797pt" height="406.611094pt" viewBox="0 0 836.594797 406.611094" xmlns="http://www.w3.org/2000/svg" version="1.1">
 <metadata>
  <rdf:RDF xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:cc="http://creativecommons.org/ns#" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
   <cc:Work>
    <dc:type rdf:resource="http://purl.org/dc/dcmitype/StillImage"/>
    <dc:date>2026-02-21T23:52:18.536757</dc:date>
    <dc:format>image/svg+xml</dc:format>
    <dc:creator>
     <cc:Agent>
      <dc:title>Matplotlib v3.10.8, https://matplotlib.org/</dc:title>
     </cc:Agent>
    </dc:creator>
   </cc:Work>
  </rdf:RDF>
 </metadata>
 <defs>
  <style type="text/css">*{stroke-linejoin: round; stroke-linecap: butt}</style>
 </defs>
 <g id="figure_1">
  <g id="patch_1">
   <path d="M 0 406.611094 
L 836.594797 406.611094 
L 836.594797 0 
L 0 0 
z
" style="fill: #ffffff"/>
  </g>
  <g id="axes_1">
   <g id="patch_2">
    <path d="M 64.030078 354.623437 
L 599.698078 354.623437 
L 599.698078 40.303437 
L 64.030078 40.303437 
z
" style="fill: #ffffff"/>
   </g>
   <g id="matplotlib.axis_1">
    <g id="xtick_1">
     <g id="line2d_1">
      <path d="M 134.443439 354.623437 
L 134.443439 40.303437 
" clip-path="url(#pc5b34efcd2)" style="fill: none; stroke: #cccccc; stroke-width: 1.5; stroke-linecap: round"/>
     </g>
     <g id="text_1">
      <!-- v0.1.0 -->
      <g style="fill: #262626" transform="translate(111.971213 378.933828) scale(0.165 -0.165)">
       <defs>
        <path id="ArialMT-76" d="M 1344 0 
L 81 3319 
L 675 3319 
L 1388 1331 
Q 1503 1009 1600 663 
Q 1675 925 1809 1294 
L 2547 3319 
L 3125 3319 
L 1869 0 
L 1344 0 
z
" transform="scale(0.015625)"/>
        <path id="ArialMT-30" d="M 266 2259 
Q 266 3072 433 3567 
Q 600 4063 929 4331 
Q 1259 4600 1759 4600 
Q 2128 4600 2406 4451 
Q 2684 4303 2865 4023 
Q 3047 3744 3150 3342 
Q 3253 2941 3253 2259 
Q 3253 1453 3087 958 
Q 2922 463 2592 192 
Q 2263 -78 1759 -78 
Q 1097 -78 719 397 
Q 266 969 266 2259 
z
M 844 2259 
Q 844 1131 1108 757 
Q 1372 384 1759 384 
Q 2147 384 2411 759 
Q 2675 1134 2675 2259 
Q 2675 3391 2411 3762 
Q 2147 4134 1753 4134 
Q 1366 4134 1134 3806 
Q 844 3388 844 2259 
z
" transform="scale(0.015625)"/>
        <path id="ArialMT-2e" d="M 581 0 
L 581 641 
L 1222 641 
L 1222 0 
L 581 0 
z
" transform="scale(0.015625)"/>
        <path id="ArialMT-31" d="M 2384 0 
L 1822 0 
L 1822 3584 
Q 1619 3391 1289 3197 
Q 959 3003 697 2906 
L 697 3450 
Q 1169 3672 1522 3987 
Q 1875 4303 2022 4600 
L 2384 4600 
L 2384 0 
z
" transform="scale(0.015625)"/>
       </defs>
       <use xlink:href="#ArialMT-76"/>
       <use xlink:href="#ArialMT-30" transform="translate(50 0)"/>
       <use xlink:href="#ArialMT-2e" transform="translate(105.615234 0)"/>
       <use xlink:href="#ArialMT-31" transform="translate(133.398438 0)"/>
       <use xlink:href="#ArialMT-2e" transform="translate(189.013672 0)"/>
       <use xlink:href="#ArialMT-30" transform="translate(216.796875 0)"/>
      </g>
     </g>
    </g>
    <g id="xtick_2">
     <g id="line2d_2">
      <path d="M 266.057199 354.623437 
L 266.057199 40.303437 
" clip-path="url(#pc5b34efcd2)" style="fill: none; stroke: #cccccc; stroke-width: 1.5; stroke-linecap: round"/>
     </g>
     <g id="text_2">
      <!-- v0.1.1 -->
      <g style="fill: #262626" transform="translate(243.584972 378.933828) scale(0.165 -0.165)">
       <use xlink:href="#ArialMT-76"/>
       <use xlink:href="#ArialMT-30" transform="translate(50 0)"/>
       <use xlink:href="#ArialMT-2e" transform="translate(105.615234 0)"/>
       <use xlink:href="#ArialMT-31" transform="translate(133.398438 0)"/>
       <use xlink:href="#ArialMT-2e" transform="translate(189.013672 0)"/>
       <use xlink:href="#ArialMT-31" transform="translate(216.796875 0)"/>
      </g>
     </g>
    </g>
    <g id="xtick_3">
     <g id="line2d_3">
      <path d="M 397.670958 354.623437 
L 397.670958 40.303437 
" clip-path="url(#pc5b34efcd2)" style="fill: none; stroke: #cccccc; stroke-width: 1.5; stroke-linecap: round"/>
     </g>
     <g id="text_3">
      <!-- v0.2.0 -->
      <g style="fill: #262626" transform="translate(375.198731 378.933828) scale(0.165 -0.165)">
       <defs>
        <path id="ArialMT-32" d="M 3222 541 
L 3222 0 
L 194 0 
Q 188 203 259 391 
Q 375 700 629 1000 
Q 884 1300 1366 1694 
Q 2113 2306 2375 2664 
Q 2638 3022 2638 3341 
Q 2638 3675 2398 3904 
Q 2159 4134 1775 4134 
Q 1369 4134 1125 3890 
Q 881 3647 878 3216 
L 300 3275 
Q 359 3922 746 4261 
Q 1134 4600 1788 4600 
Q 2447 4600 2831 4234 
Q 3216 3869 3216 3328 
Q 3216 3053 3103 2787 
Q 2991 2522 2730 2228 
Q 2469 1934 1863 1422 
Q 1356 997 1212 845 
Q 1069 694 975 541 
L 3222 541 
z
" transform="scale(0.015625)"/>
       </defs>
       <use xlink:href="#ArialMT-76"/>
       <use xlink:href="#ArialMT-30" transform="translate(50 0)"/>
       <use xlink:href="#ArialMT-2e" transform="translate(105.615234 0)"/>
       <use xlink:href="#ArialMT-32" transform="translate(133.398438 0)"/>
       <use xlink:href="#ArialMT-2e" transform="translate(189.013672 0)"/>
       <use xlink:href="#ArialMT-30" transform="translate(216.796875 0)"/>
      </g>
     </g>
    </g>
    <g id="xtick_4">
     <g id="line2d_4">
      <path d="M 529.284717 354.623437 
L 529.284717 40.303437 
" clip-path="url(#pc5b34efcd2)" style="fill: none; stroke: #cccccc; stroke-width: 1.5; stroke-linecap: round"/>
     </g>
     <g id="text_4">
      <!-- 0.3.0 (dev) -->
      <g style="fill: #262626" transform="translate(489.851006 379.134922) scale(0.165 -0.165)">
       <defs>
        <path id="ArialMT-33" d="M 269 1209 
L 831 1284 
Q 928 806 1161 595 
Q 1394 384 1728 384 
Q 2125 384 2398 659 
Q 2672 934 2672 1341 
Q 2672 1728 2419 1979 
Q 2166 2231 1775 2231 
Q 1616 2231 1378 2169 
L 1441 2663 
Q 1497 2656 1531 2656 
Q 1891 2656 2178 2843 
Q 2466 3031 2466 3422 
Q 2466 3731 2256 3934 
Q 2047 4138 1716 4138 
Q 1388 4138 1169 3931 
Q 950 3725 888 3313 
L 325 3413 
Q 428 3978 793 4289 
Q 1159 4600 1703 4600 
Q 2078 4600 2393 4439 
Q 2709 4278 2876 4000 
Q 3044 3722 3044 3409 
Q 3044 3113 2884 2869 
Q 2725 2625 2413 2481 
Q 2819 2388 3044 2092 
Q 3269 1797 3269 1353 
Q 3269 753 2831 336 
Q 2394 -81 1725 -81 
Q 1122 -81 723 278 
Q 325 638 269 1209 
z
" transform="scale(0.015625)"/>
        <path id="ArialMT-20" transform="scale(0.015625)"/>
        <path id="ArialMT-28" d="M 1497 -1347 
Q 1031 -759 709 28 
Q 388 816 388 1659 
Q 388 2403 628 3084 
Q 909 3875 1497 4659 
L 1900 4659 
Q 1522 4009 1400 3731 
Q 1209 3300 1100 2831 
Q 966 2247 966 1656 
Q 966 153 1900 -1347 
L 1497 -1347 
z
" transform="scale(0.015625)"/>
        <path id="ArialMT-64" d="M 2575 0 
L 2575 419 
Q 2259 -75 1647 -75 
Q 1250 -75 917 144 
Q 584 363 401 755 
Q 219 1147 219 1656 
Q 219 2153 384 2558 
Q 550 2963 881 3178 
Q 1213 3394 1622 3394 
Q 1922 3394 2156 3267 
Q 2391 3141 2538 2938 
L 2538 4581 
L 3097 4581 
L 3097 0 
L 2575 0 
z
M 797 1656 
Q 797 1019 1065 703 
Q 1334 388 1700 388 
Q 2069 388 2326 689 
Q 2584 991 2584 1609 
Q 2584 2291 2321 2609 
Q 2059 2928 1675 2928 
Q 1300 2928 1048 2622 
Q 797 2316 797 1656 
z
" transform="scale(0.015625)"/>
        <path id="ArialMT-65" d="M 2694 1069 
L 3275 997 
Q 3138 488 2766 206 
Q 2394 -75 1816 -75 
Q 1088 -75 661 373 
Q 234 822 234 1631 
Q 234 2469 665 2931 
Q 1097 3394 1784 3394 
Q 2450 3394 2872 2941 
Q 3294 2488 3294 1666 
Q 3294 1616 3291 1516 
L 816 1516 
Q 847 969 1125 678 
Q 1403 388 1819 388 
Q 2128 388 2347 550 
Q 2566 713 2694 1069 
z
M 847 1978 
L 2700 1978 
Q 2663 2397 2488 2606 
Q 2219 2931 1791 2931 
Q 1403 2931 1139 2672 
Q 875 2413 847 1978 
z
" transform="scale(0.015625)"/>
        <path id="ArialMT-29" d="M 791 -1347 
L 388 -1347 
Q 1322 153 1322 1656 
Q 1322 2244 1188 2822 
Q 1081 3291 891 3722 
Q 769 4003 388 4659 
L 791 4659 
Q 1378 3875 1659 3084 
Q 1900 2403 1900 1659 
Q 1900 816 1576 28 
Q 1253 -759 791 -1347 
z
" transform="scale(0.015625)"/>
       </defs>
       <use xlink:href="#ArialMT-30"/>
       <use xlink:href="#ArialMT-2e" transform="translate(55.615234 0)"/>
       <use xlink:href="#ArialMT-33" transform="translate(83.398438 0)"/>
       <use xlink:href="#ArialMT-2e" transform="translate(139.013672 0)"/>
       <use xlink:href="#ArialMT-30" transform="translate(166.796875 0)"/>
       <use xlink:href="#ArialMT-20" transform="translate(222.412109 0)"/>
       <use xlink:href="#ArialMT-28" transform="translate(250.195312 0)"/>
       <use xlink:href="#ArialMT-64" transform="translate(283.496094 0)"/>
       <use xlink:href="#ArialMT-65" transform="translate(339.111328 0)"/>
       <use xlink:href="#ArialMT-76" transform="translate(394.726562 0)"/>
       <use xlink:href="#ArialMT-29" transform="translate(444.726562 0)"/>
      </g>
     </g>
    </g>
    <g id="text_5">
     <!-- Release -->
     <g style="fill: #262626" transform="translate(306.182828 396.628594) scale(0.14 -0.14)">
      <defs>
       <path id="ArialMT-52" d="M 503 0 
L 503 4581 
L 2534 4581 
Q 3147 4581 3465 4457 
Q 3784 4334 3975 4021 
Q 4166 3709 4166 3331 
Q 4166 2844 3850 2509 
Q 3534 2175 2875 2084 
Q 3116 1969 3241 1856 
Q 3506 1613 3744 1247 
L 4541 0 
L 3778 0 
L 3172 953 
Q 2906 1366 2734 1584 
Q 2563 1803 2427 1890 
Q 2291 1978 2150 2013 
Q 2047 2034 1813 2034 
L 1109 2034 
L 1109 0 
L 503 0 
z
M 1109 2559 
L 2413 2559 
Q 2828 2559 3062 2645 
Q 3297 2731 3419 2920 
Q 3541 3109 3541 3331 
Q 3541 3656 3305 3865 
Q 3069 4075 2559 4075 
L 1109 4075 
L 1109 2559 
z
" transform="scale(0.015625)"/>
       <path id="ArialMT-6c" d="M 409 0 
L 409 4581 
L 972 4581 
L 972 0 
L 409 0 
z
" transform="scale(0.015625)"/>
       <path id="ArialMT-61" d="M 2588 409 
Q 2275 144 1986 34 
Q 1697 -75 1366 -75 
Q 819 -75 525 192 
Q 231 459 231 875 
Q 231 1119 342 1320 
Q 453 1522 633 1644 
Q 813 1766 1038 1828 
Q 1203 1872 1538 1913 
Q 2219 1994 2541 2106 
Q 2544 2222 2544 2253 
Q 2544 2597 2384 2738 
Q 2169 2928 1744 2928 
Q 1347 2928 1158 2789 
Q 969 2650 878 2297 
L 328 2372 
Q 403 2725 575 2942 
Q 747 3159 1072 3276 
Q 1397 3394 1825 3394 
Q 2250 3394 2515 3294 
Q 2781 3194 2906 3042 
Q 3031 2891 3081 2659 
Q 3109 2516 3109 2141 
L 3109 1391 
Q 3109 606 3145 398 
Q 3181 191 3288 0 
L 2700 0 
Q 2613 175 2588 409 
z
M 2541 1666 
Q 2234 1541 1622 1453 
Q 1275 1403 1131 1340 
Q 988 1278 909 1158 
Q 831 1038 831 891 
Q 831 666 1001 516 
Q 1172 366 1500 366 
Q 1825 366 2078 508 
Q 2331 650 2450 897 
Q 2541 1088 2541 1459 
L 2541 1666 
z
" transform="scale(0.015625)"/>
       <path id="ArialMT-73" d="M 197 991 
L 753 1078 
Q 800 744 1014 566 
Q 1228 388 1613 388 
Q 2000 388 2187 545 
Q 2375 703 2375 916 
Q 2375 1106 2209 1216 
Q 2094 1291 1634 1406 
Q 1016 1563 777 1677 
Q 538 1791 414 1992 
Q 291 2194 291 2438 
Q 291 2659 392 2848 
Q 494 3038 669 3163 
Q 800 3259 1026 3326 
Q 1253 3394 1513 3394 
Q 1903 3394 2198 3281 
Q 2494 3169 2634 2976 
Q 2775 2784 2828 2463 
L 2278 2388 
Q 2241 2644 2061 2787 
Q 1881 2931 1553 2931 
Q 1166 2931 1000 2803 
Q 834 2675 834 2503 
Q 834 2394 903 2306 
Q 972 2216 1119 2156 
Q 1203 2125 1616 2013 
Q 2213 1853 2448 1751 
Q 2684 1650 2818 1456 
Q 2953 1263 2953 975 
Q 2953 694 2789 445 
Q 2625 197 2315 61 
Q 2006 -75 1616 -75 
Q 969 -75 630 194 
Q 291 463 197 991 
z
" transform="scale(0.015625)"/>
      </defs>
      <use xlink:href="#ArialMT-52"/>
      <use xlink:href="#ArialMT-65" transform="translate(72.216797 0)"/>
      <use xlink:href="#ArialMT-6c" transform="translate(127.832031 0)"/>
      <use xlink:href="#ArialMT-65" transform="translate(150.048828 0)"/>
      <use xlink:href="#ArialMT-61" transform="translate(205.664062 0)"/>
      <use xlink:href="#ArialMT-73" transform="translate(261.279297 0)"/>
      <use xlink:href="#ArialMT-65" transform="translate(311.279297 0)"/>
     </g>
    </g>
   </g>
   <g id="matplotlib.axis_2">
    <g id="ytick_1">
     <g id="line2d_5">
      <path d="M 64.030078 333.668771 
L 599.698078 333.668771 
" clip-path="url(#pc5b34efcd2)" style="fill: none; stroke: #cccccc; stroke-width: 1.5; stroke-linecap: round"/>
     </g>
     <g id="text_6">
      <!-- 86 -->
      <g style="fill: #262626" transform="translate(33.178984 339.573966) scale(0.165 -0.165)">
       <defs>
        <path id="ArialMT-38" d="M 1131 2484 
Q 781 2613 612 2850 
Q 444 3088 444 3419 
Q 444 3919 803 4259 
Q 1163 4600 1759 4600 
Q 2359 4600 2725 4251 
Q 3091 3903 3091 3403 
Q 3091 3084 2923 2848 
Q 2756 2613 2416 2484 
Q 2838 2347 3058 2040 
Q 3278 1734 3278 1309 
Q 3278 722 2862 322 
Q 2447 -78 1769 -78 
Q 1091 -78 675 323 
Q 259 725 259 1325 
Q 259 1772 486 2073 
Q 713 2375 1131 2484 
z
M 1019 3438 
Q 1019 3113 1228 2906 
Q 1438 2700 1772 2700 
Q 2097 2700 2305 2904 
Q 2513 3109 2513 3406 
Q 2513 3716 2298 3927 
Q 2084 4138 1766 4138 
Q 1444 4138 1231 3931 
Q 1019 3725 1019 3438 
z
M 838 1322 
Q 838 1081 952 856 
Q 1066 631 1291 507 
Q 1516 384 1775 384 
Q 2178 384 2440 643 
Q 2703 903 2703 1303 
Q 2703 1709 2433 1975 
Q 2163 2241 1756 2241 
Q 1359 2241 1098 1978 
Q 838 1716 838 1322 
z
" transform="scale(0.015625)"/>
        <path id="ArialMT-36" d="M 3184 3459 
L 2625 3416 
Q 2550 3747 2413 3897 
Q 2184 4138 1850 4138 
Q 1581 4138 1378 3988 
Q 1113 3794 959 3422 
Q 806 3050 800 2363 
Q 1003 2672 1297 2822 
Q 1591 2972 1913 2972 
Q 2475 2972 2870 2558 
Q 3266 2144 3266 1488 
Q 3266 1056 3080 686 
Q 2894 316 2569 119 
Q 2244 -78 1831 -78 
Q 1128 -78 684 439 
Q 241 956 241 2144 
Q 241 3472 731 4075 
Q 1159 4600 1884 4600 
Q 2425 4600 2770 4297 
Q 3116 3994 3184 3459 
z
M 888 1484 
Q 888 1194 1011 928 
Q 1134 663 1356 523 
Q 1578 384 1822 384 
Q 2178 384 2434 671 
Q 2691 959 2691 1453 
Q 2691 1928 2437 2201 
Q 2184 2475 1800 2475 
Q 1419 2475 1153 2201 
Q 888 1928 888 1484 
z
" transform="scale(0.015625)"/>
       </defs>
       <use xlink:href="#ArialMT-38"/>
       <use xlink:href="#ArialMT-36" transform="translate(55.615234 0)"/>
      </g>
     </g>
    </g>
    <g id="ytick_2">
     <g id="line2d_6">
      <path d="M 64.030078 291.759438 
L 599.698078 291.759438 
" clip-path="url(#pc5b34efcd2)" style="fill: none; stroke: #cccccc; stroke-width: 1.5; stroke-linecap: round"/>
     </g>
     <g id="text_7">
      <!-- 88 -->
      <g style="fill: #262626" transform="translate(33.178984 297.664633) scale(0.165 -0.165)">
       <use xlink:href="#ArialMT-38"/>
       <use xlink:href="#ArialMT-38" transform="translate(55.615234 0)"/>
      </g>
     </g>
    </g>
    <g id="ytick_3">
     <g id="line2d_7">
      <path d="M 64.030078 249.850104 
L 599.698078 249.850104 
" clip-path="url(#pc5b34efcd2)" style="fill: none; stroke: #cccccc; stroke-width: 1.5; stroke-linecap: round"/>
     </g>
     <g id="text_8">
      <!-- 90 -->
      <g style="fill: #262626" transform="translate(33.178984 255.755299) scale(0.165 -0.165)">
       <defs>
        <path id="ArialMT-39" d="M 350 1059 
L 891 1109 
Q 959 728 1153 556 
Q 1347 384 1650 384 
Q 1909 384 2104 503 
Q 2300 622 2425 820 
Q 2550 1019 2634 1356 
Q 2719 1694 2719 2044 
Q 2719 2081 2716 2156 
Q 2547 1888 2255 1720 
Q 1963 1553 1622 1553 
Q 1053 1553 659 1965 
Q 266 2378 266 3053 
Q 266 3750 677 4175 
Q 1088 4600 1706 4600 
Q 2153 4600 2523 4359 
Q 2894 4119 3086 3673 
Q 3278 3228 3278 2384 
Q 3278 1506 3087 986 
Q 2897 466 2520 194 
Q 2144 -78 1638 -78 
Q 1100 -78 759 220 
Q 419 519 350 1059 
z
M 2653 3081 
Q 2653 3566 2395 3850 
Q 2138 4134 1775 4134 
Q 1400 4134 1122 3828 
Q 844 3522 844 3034 
Q 844 2597 1108 2323 
Q 1372 2050 1759 2050 
Q 2150 2050 2401 2323 
Q 2653 2597 2653 3081 
z
" transform="scale(0.015625)"/>
       </defs>
       <use xlink:href="#ArialMT-39"/>
       <use xlink:href="#ArialMT-30" transform="translate(55.615234 0)"/>
      </g>
     </g>
    </g>
    <g id="ytick_4">
     <g id="line2d_8">
      <path d="M 64.030078 207.940771 
L 599.698078 207.940771 
" clip-path="url(#pc5b34efcd2)" style="fill: none; stroke: #cccccc; stroke-width: 1.5; stroke-linecap: round"/>
     </g>
     <g id="text_9">
      <!-- 92 -->
      <g style="fill: #262626" transform="translate(33.178984 213.845966) scale(0.165 -0.165)">
       <use xlink:href="#ArialMT-39"/>
       <use xlink:href="#ArialMT-32" transform="translate(55.615234 0)"/>
      </g>
     </g>
    </g>
    <g id="ytick_5">
     <g id="line2d_9">
      <path d="M 64.030078 166.031438 
L 599.698078 166.031438 
" clip-path="url(#pc5b34efcd2)" style="fill: none; stroke: #cccccc; stroke-width: 1.5; stroke-linecap: round"/>
     </g>
     <g id="text_10">
      <!-- 94 -->
      <g style="fill: #262626" transform="translate(33.178984 171.936633) scale(0.165 -0.165)">
       <defs>
        <path id="ArialMT-34" d="M 2069 0 
L 2069 1097 
L 81 1097 
L 81 1613 
L 2172 4581 
L 2631 4581 
L 2631 1613 
L 3250 1613 
L 3250 1097 
L 2631 1097 
L 2631 0 
L 2069 0 
z
M 2069 1613 
L 2069 3678 
L 634 1613 
L 2069 1613 
z
" transform="scale(0.015625)"/>
       </defs>
       <use xlink:href="#ArialMT-39"/>
       <use xlink:href="#ArialMT-34" transform="translate(55.615234 0)"/>
      </g>
     </g>
    </g>
    <g id="ytick_6">
     <g id="line2d_10">
      <path d="M 64.030078 124.122104 
L 599.698078 124.122104 
" clip-path="url(#pc5b34efcd2)" style="fill: none; stroke: #cccccc; stroke-width: 1.5; stroke-linecap: round"/>
     </g>
     <g id="text_11">
      <!-- 96 -->
      <g style="fill: #262626" transform="translate(33.178984 130.027299) scale(0.165 -0.165)">
       <use xlink:href="#ArialMT-39"/>
       <use xlink:href="#ArialMT-36" transform="translate(55.615234 0)"/>
      </g>
     </g>
    </g>
    <g id="ytick_7">
     <g id="line2d_11">
      <path d="M 64.030078 82.212771 
L 599.698078 82.212771 
" clip-path="url(#pc5b34efcd2)" style="fill: none; stroke: #cccccc; stroke-width: 1.5; stroke-linecap: round"/>
     </g>
     <g id="text_12">
      <!-- 98 -->
      <g style="fill: #262626" transform="translate(33.178984 88.117966) scale(0.165 -0.165)">
       <use xlink:href="#ArialMT-39"/>
       <use xlink:href="#ArialMT-38" transform="translate(55.615234 0)"/>
      </g>
     </g>
    </g>
    <g id="ytick_8">
     <g id="line2d_12">
      <path d="M 64.030078 40.303438 
L 599.698078 40.303438 
" clip-path="url(#pc5b34efcd2)" style="fill: none; stroke: #cccccc; stroke-width: 1.5; stroke-linecap: round"/>
     </g>
     <g id="text_13">
      <!-- 100 -->
      <g style="fill: #262626" transform="translate(24.003438 46.208633) scale(0.165 -0.165)">
       <use xlink:href="#ArialMT-31"/>
       <use xlink:href="#ArialMT-30" transform="translate(55.615234 0)"/>
       <use xlink:href="#ArialMT-30" transform="translate(111.230469 0)"/>
      </g>
     </g>
    </g>
    <g id="text_14">
     <!-- Match % -->
     <g style="fill: #262626" transform="translate(17.220937 224.693437) rotate(-90) scale(0.14 -0.14)">
      <defs>
       <path id="ArialMT-4d" d="M 475 0 
L 475 4581 
L 1388 4581 
L 2472 1338 
Q 2622 884 2691 659 
Q 2769 909 2934 1394 
L 4031 4581 
L 4847 4581 
L 4847 0 
L 4263 0 
L 4263 3834 
L 2931 0 
L 2384 0 
L 1059 3900 
L 1059 0 
L 475 0 
z
" transform="scale(0.015625)"/>
       <path id="ArialMT-74" d="M 1650 503 
L 1731 6 
Q 1494 -44 1306 -44 
Q 1000 -44 831 53 
Q 663 150 594 308 
Q 525 466 525 972 
L 525 2881 
L 113 2881 
L 113 3319 
L 525 3319 
L 525 4141 
L 1084 4478 
L 1084 3319 
L 1650 3319 
L 1650 2881 
L 1084 2881 
L 1084 941 
Q 1084 700 1114 631 
Q 1144 563 1211 522 
Q 1278 481 1403 481 
Q 1497 481 1650 503 
z
" transform="scale(0.015625)"/>
       <path id="ArialMT-63" d="M 2588 1216 
L 3141 1144 
Q 3050 572 2676 248 
Q 2303 -75 1759 -75 
Q 1078 -75 664 370 
Q 250 816 250 1647 
Q 250 2184 428 2587 
Q 606 2991 970 3192 
Q 1334 3394 1763 3394 
Q 2303 3394 2647 3120 
Q 2991 2847 3088 2344 
L 2541 2259 
Q 2463 2594 2264 2762 
Q 2066 2931 1784 2931 
Q 1359 2931 1093 2626 
Q 828 2322 828 1663 
Q 828 994 1084 691 
Q 1341 388 1753 388 
Q 2084 388 2306 591 
Q 2528 794 2588 1216 
z
" transform="scale(0.015625)"/>
       <path id="ArialMT-68" d="M 422 0 
L 422 4581 
L 984 4581 
L 984 2938 
Q 1378 3394 1978 3394 
Q 2347 3394 2619 3248 
Q 2891 3103 3008 2847 
Q 3125 2591 3125 2103 
L 3125 0 
L 2563 0 
L 2563 2103 
Q 2563 2525 2380 2717 
Q 2197 2909 1863 2909 
Q 1613 2909 1392 2779 
Q 1172 2650 1078 2428 
Q 984 2206 984 1816 
L 984 0 
L 422 0 
z
" transform="scale(0.015625)"/>
       <path id="ArialMT-25" d="M 372 3481 
Q 372 3972 619 4315 
Q 866 4659 1334 4659 
Q 1766 4659 2048 4351 
Q 2331 4044 2331 3447 
Q 2331 2866 2045 2552 
Q 1759 2238 1341 2238 
Q 925 2238 648 2547 
Q 372 2856 372 3481 
z
M 1350 4272 
Q 1141 4272 1002 4090 
Q 863 3909 863 3425 
Q 863 2984 1003 2804 
Q 1144 2625 1350 2625 
Q 1563 2625 1702 2806 
Q 1841 2988 1841 3469 
Q 1841 3913 1700 4092 
Q 1559 4272 1350 4272 
z
M 1353 -169 
L 3859 4659 
L 4316 4659 
L 1819 -169 
L 1353 -169 
z
M 3334 1075 
Q 3334 1569 3581 1911 
Q 3828 2253 4300 2253 
Q 4731 2253 5014 1945 
Q 5297 1638 5297 1041 
Q 5297 459 5011 145 
Q 4725 -169 4303 -169 
Q 3888 -169 3611 142 
Q 3334 453 3334 1075 
z
M 4316 1866 
Q 4103 1866 3964 1684 
Q 3825 1503 3825 1019 
Q 3825 581 3965 400 
Q 4106 219 4313 219 
Q 4528 219 4667 400 
Q 4806 581 4806 1063 
Q 4806 1506 4665 1686 
Q 4525 1866 4316 1866 
z
" transform="scale(0.015625)"/>
      </defs>
      <use xlink:href="#ArialMT-4d"/>
      <use xlink:href="#ArialMT-61" transform="translate(83.300781 0)"/>
      <use xlink:href="#ArialMT-74" transform="translate(138.916016 0)"/>
      <use xlink:href="#ArialMT-63" transform="translate(166.699219 0)"/>
      <use xlink:href="#ArialMT-68" transform="translate(216.699219 0)"/>
      <use xlink:href="#ArialMT-20" transform="translate(272.314453 0)"/>
      <use xlink:href="#ArialMT-25" transform="translate(300.097656 0)"/>
     </g>
    </g>
   </g>
   <g id="patch_3">
    <path d="M 88.378624 2135.770104 
L 134.443439 2135.770104 
L 134.443439 318.497592 
L 88.378624 318.497592 
z
" clip-path="url(#pc5b34efcd2)" style="fill: #3498db; stroke: #ffffff; stroke-width: 1.5; stroke-linejoin: miter"/>
   </g>
   <g id="patch_4">
    <path d="M 219.992383 2135.770104 
L 266.057199 2135.770104 
L 266.057199 318.497592 
L 219.992383 318.497592 
z
" clip-path="url(#pc5b34efcd2)" style="fill: #3498db; stroke: #ffffff; stroke-width: 1.5; stroke-linejoin: miter"/>
   </g>
   <g id="patch_5">
    <path d="M 351.606142 2135.770104 
L 397.670958 2135.770104 
L 397.670958 169.67755 
L 351.606142 169.67755 
z
" clip-path="url(#pc5b34efcd2)" style="fill: #3498db; stroke: #ffffff; stroke-width: 1.5; stroke-linejoin: miter"/>
   </g>
   <g id="patch_6">
    <path d="M 483.219901 2135.770104 
L 529.284717 2135.770104 
L 529.284717 169.195592 
L 483.219901 169.195592 
z
" clip-path="url(#pc5b34efcd2)" style="fill: #3498db; stroke: #ffffff; stroke-width: 1.5; stroke-linejoin: miter"/>
   </g>
   <g id="patch_7">
    <path d="M 88.378624 318.497592 
L 134.443439 318.497592 
L 134.443439 146.417869 
L 88.378624 146.417869 
z
" clip-path="url(#pc5b34efcd2)" style="fill: #85c1e9; stroke: #ffffff; stroke-width: 1.5; stroke-linejoin: miter"/>
   </g>
   <g id="patch_8">
    <path d="M 219.992383 318.497592 
L 266.057199 318.497592 
L 266.057199 146.417869 
L 219.992383 146.417869 
z
" clip-path="url(#pc5b34efcd2)" style="fill: #85c1e9; stroke: #ffffff; stroke-width: 1.5; stroke-linejoin: miter"/>
   </g>
   <g id="patch_9">
    <path d="M 351.606142 169.67755 
L 397.670958 169.67755 
L 397.670958 67.376867 
L 351.606142 67.376867 
z
" clip-path="url(#pc5b34efcd2)" style="fill: #85c1e9; stroke: #ffffff; stroke-width: 1.5; stroke-linejoin: miter"/>
   </g>
   <g id="patch_10">
    <path d="M 483.219901 169.195592 
L 529.284717 169.195592 
L 529.284717 66.035768 
L 483.219901 66.035768 
z
" clip-path="url(#pc5b34efcd2)" style="fill: #85c1e9; stroke: #ffffff; stroke-width: 1.5; stroke-linejoin: miter"/>
   </g>
   <g id="patch_11">
    <path d="M 134.443439 2135.770104 
L 180.508255 2135.770104 
L 180.508255 363.906355 
L 134.443439 363.906355 
z
" clip-path="url(#pc5b34efcd2)" style="fill: #27ae60; stroke: #ffffff; stroke-width: 1.5; stroke-linejoin: miter"/>
   </g>
   <g id="patch_12">
    <path d="M 266.057199 2135.770104 
L 312.122014 2135.770104 
L 312.122014 363.906355 
L 266.057199 363.906355 
z
" clip-path="url(#pc5b34efcd2)" style="fill: #27ae60; stroke: #ffffff; stroke-width: 1.5; stroke-linejoin: miter"/>
   </g>
   <g id="patch_13">
    <path d="M 397.670958 2135.770104 
L 443.735773 2135.770104 
L 443.735773 179.924382 
L 397.670958 179.924382 
z
" clip-path="url(#pc5b34efcd2)" style="fill: #27ae60; stroke: #ffffff; stroke-width: 1.5; stroke-linejoin: miter"/>
   </g>
   <g id="patch_14">
    <path d="M 529.284717 2135.770104 
L 575.349533 2135.770104 
L 575.349533 179.924382 
L 529.284717 179.924382 
z
" clip-path="url(#pc5b34efcd2)" style="fill: #27ae60; stroke: #ffffff; stroke-width: 1.5; stroke-linejoin: miter"/>
   </g>
   <g id="patch_15">
    <path d="M 134.443439 363.906355 
L 180.508255 363.906355 
L 180.508255 265.545149 
L 134.443439 265.545149 
z
" clip-path="url(#pc5b34efcd2)" style="fill: #7dcea0; stroke: #ffffff; stroke-width: 1.5; stroke-linejoin: miter"/>
   </g>
   <g id="patch_16">
    <path d="M 266.057199 363.906355 
L 312.122014 363.906355 
L 312.122014 265.545149 
L 266.057199 265.545149 
z
" clip-path="url(#pc5b34efcd2)" style="fill: #7dcea0; stroke: #ffffff; stroke-width: 1.5; stroke-linejoin: miter"/>
   </g>
   <g id="patch_17">
    <path d="M 397.670958 179.924381 
L 443.735773 179.924381 
L 443.735773 77.770382 
L 397.670958 77.770382 
z
" clip-path="url(#pc5b34efcd2)" style="fill: #7dcea0; stroke: #ffffff; stroke-width: 1.5; stroke-linejoin: miter"/>
   </g>
   <g id="patch_18">
    <path d="M 529.284717 179.924381 
L 575.349533 179.924381 
L 575.349533 77.770382 
L 529.284717 77.770382 
z
" clip-path="url(#pc5b34efcd2)" style="fill: #7dcea0; stroke: #ffffff; stroke-width: 1.5; stroke-linejoin: miter"/>
   </g>
   <g id="patch_19">
    <path d="M 64.030078 354.623438 
L 64.030078 40.303438 
" style="fill: none; stroke: #cccccc; stroke-width: 1.875; stroke-linejoin: miter; stroke-linecap: square"/>
   </g>
   <g id="patch_20">
    <path d="M 599.698078 354.623438 
L 599.698078 40.303438 
" style="fill: none; stroke: #cccccc; stroke-width: 1.875; stroke-linejoin: miter; stroke-linecap: square"/>
   </g>
   <g id="patch_21">
    <path d="M 64.030078 354.623437 
L 599.698078 354.623437 
" style="fill: none; stroke: #cccccc; stroke-width: 1.875; stroke-linejoin: miter; stroke-linecap: square"/>
   </g>
   <g id="patch_22">
    <path d="M 64.030078 40.303437 
L 599.698078 40.303437 
" style="fill: none; stroke: #cccccc; stroke-width: 1.875; stroke-linejoin: miter; stroke-linecap: square"/>
   </g>
   <g id="text_15">
    <!-- Protein Projection Performance (100k ClinVar Variants) -->
    <g style="fill: #262626" transform="translate(112.607203 20.303437) scale(0.18 -0.18)">
     <defs>
      <path id="ArialMT-50" d="M 494 0 
L 494 4581 
L 2222 4581 
Q 2678 4581 2919 4538 
Q 3256 4481 3484 4323 
Q 3713 4166 3852 3881 
Q 3991 3597 3991 3256 
Q 3991 2672 3619 2267 
Q 3247 1863 2275 1863 
L 1100 1863 
L 1100 0 
L 494 0 
z
M 1100 2403 
L 2284 2403 
Q 2872 2403 3119 2622 
Q 3366 2841 3366 3238 
Q 3366 3525 3220 3729 
Q 3075 3934 2838 4000 
Q 2684 4041 2272 4041 
L 1100 4041 
L 1100 2403 
z
" transform="scale(0.015625)"/>
      <path id="ArialMT-72" d="M 416 0 
L 416 3319 
L 922 3319 
L 922 2816 
Q 1116 3169 1280 3281 
Q 1444 3394 1641 3394 
Q 1925 3394 2219 3213 
L 2025 2691 
Q 1819 2813 1613 2813 
Q 1428 2813 1281 2702 
Q 1134 2591 1072 2394 
Q 978 2094 978 1738 
L 978 0 
L 416 0 
z
" transform="scale(0.015625)"/>
      <path id="ArialMT-6f" d="M 213 1659 
Q 213 2581 725 3025 
Q 1153 3394 1769 3394 
Q 2453 3394 2887 2945 
Q 3322 2497 3322 1706 
Q 3322 1066 3130 698 
Q 2938 331 2570 128 
Q 2203 -75 1769 -75 
Q 1072 -75 642 372 
Q 213 819 213 1659 
z
M 791 1659 
Q 791 1022 1069 705 
Q 1347 388 1769 388 
Q 2188 388 2466 706 
Q 2744 1025 2744 1678 
Q 2744 2294 2464 2611 
Q 2184 2928 1769 2928 
Q 1347 2928 1069 2612 
Q 791 2297 791 1659 
z
" transform="scale(0.015625)"/>
      <path id="ArialMT-69" d="M 425 3934 
L 425 4581 
L 988 4581 
L 988 3934 
L 425 3934 
z
M 425 0 
L 425 3319 
L 988 3319 
L 988 0 
L 425 0 
z
" transform="scale(0.015625)"/>
      <path id="ArialMT-6e" d="M 422 0 
L 422 3319 
L 928 3319 
L 928 2847 
Q 1294 3394 1984 3394 
Q 2284 3394 2536 3286 
Q 2788 3178 2913 3003 
Q 3038 2828 3088 2588 
Q 3119 2431 3119 2041 
L 3119 0 
L 2556 0 
L 2556 2019 
Q 2556 2363 2490 2533 
Q 2425 2703 2258 2804 
Q 2091 2906 1866 2906 
Q 1506 2906 1245 2678 
Q 984 2450 984 1813 
L 984 0 
L 422 0 
z
" transform="scale(0.015625)"/>
      <path id="ArialMT-6a" d="M 419 3928 
L 419 4581 
L 981 4581 
L 981 3928 
L 419 3928 
z
M -294 -1288 
L -188 -809 
Q -19 -853 78 -853 
Q 250 -853 334 -739 
Q 419 -625 419 -169 
L 419 3319 
L 981 3319 
L 981 -181 
Q 981 -794 822 -1034 
Q 619 -1347 147 -1347 
Q -81 -1347 -294 -1288 
z
" transform="scale(0.015625)"/>
      <path id="ArialMT-66" d="M 556 0 
L 556 2881 
L 59 2881 
L 59 3319 
L 556 3319 
L 556 3672 
Q 556 4006 616 4169 
Q 697 4388 901 4523 
Q 1106 4659 1475 4659 
Q 1713 4659 2000 4603 
L 1916 4113 
Q 1741 4144 1584 4144 
Q 1328 4144 1222 4034 
Q 1116 3925 1116 3625 
L 1116 3319 
L 1763 3319 
L 1763 2881 
L 1116 2881 
L 1116 0 
L 556 0 
z
" transform="scale(0.015625)"/>
      <path id="ArialMT-6d" d="M 422 0 
L 422 3319 
L 925 3319 
L 925 2853 
Q 1081 3097 1340 3245 
Q 1600 3394 1931 3394 
Q 2300 3394 2536 3241 
Q 2772 3088 2869 2813 
Q 3263 3394 3894 3394 
Q 4388 3394 4653 3120 
Q 4919 2847 4919 2278 
L 4919 0 
L 4359 0 
L 4359 2091 
Q 4359 2428 4304 2576 
Q 4250 2725 4106 2815 
Q 3963 2906 3769 2906 
Q 3419 2906 3187 2673 
Q 2956 2441 2956 1928 
L 2956 0 
L 2394 0 
L 2394 2156 
Q 2394 2531 2256 2718 
Q 2119 2906 1806 2906 
Q 1569 2906 1367 2781 
Q 1166 2656 1075 2415 
Q 984 2175 984 1722 
L 984 0 
L 422 0 
z
" transform="scale(0.015625)"/>
      <path id="ArialMT-6b" d="M 425 0 
L 425 4581 
L 988 4581 
L 988 1969 
L 2319 3319 
L 3047 3319 
L 1778 2088 
L 3175 0 
L 2481 0 
L 1384 1697 
L 988 1316 
L 988 0 
L 425 0 
z
" transform="scale(0.015625)"/>
      <path id="ArialMT-43" d="M 3763 1606 
L 4369 1453 
Q 4178 706 3683 314 
Q 3188 -78 2472 -78 
Q 1731 -78 1267 223 
Q 803 525 561 1097 
Q 319 1669 319 2325 
Q 319 3041 592 3573 
Q 866 4106 1370 4382 
Q 1875 4659 2481 4659 
Q 3169 4659 3637 4309 
Q 4106 3959 4291 3325 
L 3694 3184 
Q 3534 3684 3231 3912 
Q 2928 4141 2469 4141 
Q 1941 4141 1586 3887 
Q 1231 3634 1087 3207 
Q 944 2781 944 2328 
Q 944 1744 1114 1308 
Q 1284 872 1643 656 
Q 2003 441 2422 441 
Q 2931 441 3284 734 
Q 3638 1028 3763 1606 
z
" transform="scale(0.015625)"/>
      <path id="ArialMT-56" d="M 1803 0 
L 28 4581 
L 684 4581 
L 1875 1253 
Q 2019 853 2116 503 
Q 2222 878 2363 1253 
L 3600 4581 
L 4219 4581 
L 2425 0 
L 1803 0 
z
" transform="scale(0.015625)"/>
     </defs>
     <use xlink:href="#ArialMT-50"/>
     <use xlink:href="#ArialMT-72" transform="translate(66.699219 0)"/>
     <use xlink:href="#ArialMT-6f" transform="translate(100 0)"/>
     <use xlink:href="#ArialMT-74" transform="translate(155.615234 0)"/>
     <use xlink:href="#ArialMT-65" transform="translate(183.398438 0)"/>
     <use xlink:href="#ArialMT-69" transform="translate(239.013672 0)"/>
     <use xlink:href="#ArialMT-6e" transform="translate(261.230469 0)"/>
     <use xlink:href="#ArialMT-20" transform="translate(316.845703 0)"/>
     <use xlink:href="#ArialMT-50" transform="translate(344.628906 0)"/>
     <use xlink:href="#ArialMT-72" transform="translate(411.328125 0)"/>
     <use xlink:href="#ArialMT-6f" transform="translate(444.628906 0)"/>
     <use xlink:href="#ArialMT-6a" transform="translate(500.244141 0)"/>
     <use xlink:href="#ArialMT-65" transform="translate(522.460938 0)"/>
     <use xlink:href="#ArialMT-63" transform="translate(578.076172 0)"/>
     <use xlink:href="#ArialMT-74" transform="translate(628.076172 0)"/>
     <use xlink:href="#ArialMT-69" transform="translate(655.859375 0)"/>
     <use xlink:href="#ArialMT-6f" transform="translate(678.076172 0)"/>
     <use xlink:href="#ArialMT-6e" transform="translate(733.691406 0)"/>
     <use xlink:href="#ArialMT-20" transform="translate(789.306641 0)"/>
     <use xlink:href="#ArialMT-50" transform="translate(817.089844 0)"/>
     <use xlink:href="#ArialMT-65" transform="translate(883.789062 0)"/>
     <use xlink:href="#ArialMT-72" transform="translate(939.404297 0)"/>
     <use xlink:href="#ArialMT-66" transform="translate(972.705078 0)"/>
     <use xlink:href="#ArialMT-6f" transform="translate(1000.488281 0)"/>
     <use xlink:href="#ArialMT-72" transform="translate(1056.103516 0)"/>
     <use xlink:href="#ArialMT-6d" transform="translate(1089.404297 0)"/>
     <use xlink:href="#ArialMT-61" transform="translate(1172.705078 0)"/>
     <use xlink:href="#ArialMT-6e" transform="translate(1228.320312 0)"/>
     <use xlink:href="#ArialMT-63" transform="translate(1283.935547 0)"/>
     <use xlink:href="#ArialMT-65" transform="translate(1333.935547 0)"/>
     <use xlink:href="#ArialMT-20" transform="translate(1389.550781 0)"/>
     <use xlink:href="#ArialMT-28" transform="translate(1417.333984 0)"/>
     <use xlink:href="#ArialMT-31" transform="translate(1450.634766 0)"/>
     <use xlink:href="#ArialMT-30" transform="translate(1506.25 0)"/>
     <use xlink:href="#ArialMT-30" transform="translate(1561.865234 0)"/>
     <use xlink:href="#ArialMT-6b" transform="translate(1617.480469 0)"/>
     <use xlink:href="#ArialMT-20" transform="translate(1667.480469 0)"/>
     <use xlink:href="#ArialMT-43" transform="translate(1695.263672 0)"/>
     <use xlink:href="#ArialMT-6c" transform="translate(1767.480469 0)"/>
     <use xlink:href="#ArialMT-69" transform="translate(1789.697266 0)"/>
     <use xlink:href="#ArialMT-6e" transform="translate(1811.914062 0)"/>
     <use xlink:href="#ArialMT-56" transform="translate(1867.529297 0)"/>
     <use xlink:href="#ArialMT-61" transform="translate(1926.853516 0)"/>
     <use xlink:href="#ArialMT-72" transform="translate(1982.46875 0)"/>
     <use xlink:href="#ArialMT-20" transform="translate(2015.769531 0)"/>
     <use xlink:href="#ArialMT-56" transform="translate(2043.552734 0)"/>
     <use xlink:href="#ArialMT-61" transform="translate(2102.876953 0)"/>
     <use xlink:href="#ArialMT-72" transform="translate(2158.492188 0)"/>
     <use xlink:href="#ArialMT-69" transform="translate(2191.792969 0)"/>
     <use xlink:href="#ArialMT-61" transform="translate(2214.009766 0)"/>
     <use xlink:href="#ArialMT-6e" transform="translate(2269.625 0)"/>
     <use xlink:href="#ArialMT-74" transform="translate(2325.240234 0)"/>
     <use xlink:href="#ArialMT-73" transform="translate(2353.023438 0)"/>
     <use xlink:href="#ArialMT-29" transform="translate(2403.023438 0)"/>
    </g>
   </g>
   <g id="legend_1">
    <g id="patch_23">
     <path d="M 613.711438 139.788125 
L 826.094797 139.788125 
Q 829.394797 139.788125 829.394797 136.488125 
L 829.394797 43.603437 
Q 829.394797 40.303437 826.094797 40.303437 
L 613.711438 40.303437 
Q 610.411438 40.303437 610.411438 43.603437 
L 610.411438 136.488125 
Q 610.411438 139.788125 613.711438 139.788125 
z
" style="fill: #ffffff; opacity: 0.8; stroke: #cccccc; stroke-width: 1.5; stroke-linejoin: miter"/>
    </g>
    <g id="patch_24">
     <path d="M 617.011438 58.713828 
L 650.011438 58.713828 
L 650.011438 47.163828 
L 617.011438 47.163828 
z
" style="fill: #3498db; stroke: #ffffff; stroke-width: 1.5; stroke-linejoin: miter"/>
    </g>
    <g id="text_16">
     <!-- Weaver Identity -->
     <g style="fill: #262626" transform="translate(663.211438 58.713828) scale(0.165 -0.165)">
      <defs>
       <path id="ArialMT-57" d="M 1294 0 
L 78 4581 
L 700 4581 
L 1397 1578 
Q 1509 1106 1591 641 
Q 1766 1375 1797 1488 
L 2669 4581 
L 3400 4581 
L 4056 2263 
Q 4303 1400 4413 641 
Q 4500 1075 4641 1638 
L 5359 4581 
L 5969 4581 
L 4713 0 
L 4128 0 
L 3163 3491 
Q 3041 3928 3019 4028 
Q 2947 3713 2884 3491 
L 1913 0 
L 1294 0 
z
" transform="scale(0.015625)"/>
       <path id="ArialMT-49" d="M 597 0 
L 597 4581 
L 1203 4581 
L 1203 0 
L 597 0 
z
" transform="scale(0.015625)"/>
       <path id="ArialMT-79" d="M 397 -1278 
L 334 -750 
Q 519 -800 656 -800 
Q 844 -800 956 -737 
Q 1069 -675 1141 -563 
Q 1194 -478 1313 -144 
Q 1328 -97 1363 -6 
L 103 3319 
L 709 3319 
L 1400 1397 
Q 1534 1031 1641 628 
Q 1738 1016 1872 1384 
L 2581 3319 
L 3144 3319 
L 1881 -56 
Q 1678 -603 1566 -809 
Q 1416 -1088 1222 -1217 
Q 1028 -1347 759 -1347 
Q 597 -1347 397 -1278 
z
" transform="scale(0.015625)"/>
      </defs>
      <use xlink:href="#ArialMT-57"/>
      <use xlink:href="#ArialMT-65" transform="translate(92.634766 0)"/>
      <use xlink:href="#ArialMT-61" transform="translate(148.25 0)"/>
      <use xlink:href="#ArialMT-76" transform="translate(203.865234 0)"/>
      <use xlink:href="#ArialMT-65" transform="translate(253.865234 0)"/>
      <use xlink:href="#ArialMT-72" transform="translate(309.480469 0)"/>
      <use xlink:href="#ArialMT-20" transform="translate(342.78125 0)"/>
      <use xlink:href="#ArialMT-49" transform="translate(370.564453 0)"/>
      <use xlink:href="#ArialMT-64" transform="translate(398.347656 0)"/>
      <use xlink:href="#ArialMT-65" transform="translate(453.962891 0)"/>
      <use xlink:href="#ArialMT-6e" transform="translate(509.578125 0)"/>
      <use xlink:href="#ArialMT-74" transform="translate(565.193359 0)"/>
      <use xlink:href="#ArialMT-69" transform="translate(592.976562 0)"/>
      <use xlink:href="#ArialMT-74" transform="translate(615.193359 0)"/>
      <use xlink:href="#ArialMT-79" transform="translate(642.976562 0)"/>
     </g>
    </g>
    <g id="patch_25">
     <path d="M 617.011438 82.246953 
L 650.011438 82.246953 
L 650.011438 70.696953 
L 617.011438 70.696953 
z
" style="fill: #85c1e9; stroke: #ffffff; stroke-width: 1.5; stroke-linejoin: miter"/>
    </g>
    <g id="text_17">
     <!-- Weaver Analogous -->
     <g style="fill: #262626" transform="translate(663.211438 82.246953) scale(0.165 -0.165)">
      <defs>
       <path id="ArialMT-41" d="M -9 0 
L 1750 4581 
L 2403 4581 
L 4278 0 
L 3588 0 
L 3053 1388 
L 1138 1388 
L 634 0 
L -9 0 
z
M 1313 1881 
L 2866 1881 
L 2388 3150 
Q 2169 3728 2063 4100 
Q 1975 3659 1816 3225 
L 1313 1881 
z
" transform="scale(0.015625)"/>
       <path id="ArialMT-67" d="M 319 -275 
L 866 -356 
Q 900 -609 1056 -725 
Q 1266 -881 1628 -881 
Q 2019 -881 2231 -725 
Q 2444 -569 2519 -288 
Q 2563 -116 2559 434 
Q 2191 0 1641 0 
Q 956 0 581 494 
Q 206 988 206 1678 
Q 206 2153 378 2554 
Q 550 2956 876 3175 
Q 1203 3394 1644 3394 
Q 2231 3394 2613 2919 
L 2613 3319 
L 3131 3319 
L 3131 450 
Q 3131 -325 2973 -648 
Q 2816 -972 2473 -1159 
Q 2131 -1347 1631 -1347 
Q 1038 -1347 672 -1080 
Q 306 -813 319 -275 
z
M 784 1719 
Q 784 1066 1043 766 
Q 1303 466 1694 466 
Q 2081 466 2343 764 
Q 2606 1063 2606 1700 
Q 2606 2309 2336 2618 
Q 2066 2928 1684 2928 
Q 1309 2928 1046 2623 
Q 784 2319 784 1719 
z
" transform="scale(0.015625)"/>
       <path id="ArialMT-75" d="M 2597 0 
L 2597 488 
Q 2209 -75 1544 -75 
Q 1250 -75 995 37 
Q 741 150 617 320 
Q 494 491 444 738 
Q 409 903 409 1263 
L 409 3319 
L 972 3319 
L 972 1478 
Q 972 1038 1006 884 
Q 1059 663 1231 536 
Q 1403 409 1656 409 
Q 1909 409 2131 539 
Q 2353 669 2445 892 
Q 2538 1116 2538 1541 
L 2538 3319 
L 3100 3319 
L 3100 0 
L 2597 0 
z
" transform="scale(0.015625)"/>
      </defs>
      <use xlink:href="#ArialMT-57"/>
      <use xlink:href="#ArialMT-65" transform="translate(92.634766 0)"/>
      <use xlink:href="#ArialMT-61" transform="translate(148.25 0)"/>
      <use xlink:href="#ArialMT-76" transform="translate(203.865234 0)"/>
      <use xlink:href="#ArialMT-65" transform="translate(253.865234 0)"/>
      <use xlink:href="#ArialMT-72" transform="translate(309.480469 0)"/>
      <use xlink:href="#ArialMT-20" transform="translate(342.78125 0)"/>
      <use xlink:href="#ArialMT-41" transform="translate(365.064453 0)"/>
      <use xlink:href="#ArialMT-6e" transform="translate(431.763672 0)"/>
      <use xlink:href="#ArialMT-61" transform="translate(487.378906 0)"/>
      <use xlink:href="#ArialMT-6c" transform="translate(542.994141 0)"/>
      <use xlink:href="#ArialMT-6f" transform="translate(565.210938 0)"/>
      <use xlink:href="#ArialMT-67" transform="translate(620.826172 0)"/>
      <use xlink:href="#ArialMT-6f" transform="translate(676.441406 0)"/>
      <use xlink:href="#ArialMT-75" transform="translate(732.056641 0)"/>
      <use xlink:href="#ArialMT-73" transform="translate(787.671875 0)"/>
     </g>
    </g>
    <g id="patch_26">
     <path d="M 617.011438 105.981172 
L 650.011438 105.981172 
L 650.011438 94.431172 
L 617.011438 94.431172 
z
" style="fill: #27ae60; stroke: #ffffff; stroke-width: 1.5; stroke-linejoin: miter"/>
    </g>
    <g id="text_18">
     <!-- Ref-HGVS Identity -->
     <g style="fill: #262626" transform="translate(663.211438 105.981172) scale(0.165 -0.165)">
      <defs>
       <path id="ArialMT-2d" d="M 203 1375 
L 203 1941 
L 1931 1941 
L 1931 1375 
L 203 1375 
z
" transform="scale(0.015625)"/>
       <path id="ArialMT-48" d="M 513 0 
L 513 4581 
L 1119 4581 
L 1119 2700 
L 3500 2700 
L 3500 4581 
L 4106 4581 
L 4106 0 
L 3500 0 
L 3500 2159 
L 1119 2159 
L 1119 0 
L 513 0 
z
" transform="scale(0.015625)"/>
       <path id="ArialMT-47" d="M 2638 1797 
L 2638 2334 
L 4578 2338 
L 4578 638 
Q 4131 281 3656 101 
Q 3181 -78 2681 -78 
Q 2006 -78 1454 211 
Q 903 500 622 1047 
Q 341 1594 341 2269 
Q 341 2938 620 3517 
Q 900 4097 1425 4378 
Q 1950 4659 2634 4659 
Q 3131 4659 3532 4498 
Q 3934 4338 4162 4050 
Q 4391 3763 4509 3300 
L 3963 3150 
Q 3859 3500 3706 3700 
Q 3553 3900 3268 4020 
Q 2984 4141 2638 4141 
Q 2222 4141 1919 4014 
Q 1616 3888 1430 3681 
Q 1244 3475 1141 3228 
Q 966 2803 966 2306 
Q 966 1694 1177 1281 
Q 1388 869 1791 669 
Q 2194 469 2647 469 
Q 3041 469 3416 620 
Q 3791 772 3984 944 
L 3984 1797 
L 2638 1797 
z
" transform="scale(0.015625)"/>
       <path id="ArialMT-53" d="M 288 1472 
L 859 1522 
Q 900 1178 1048 958 
Q 1197 738 1509 602 
Q 1822 466 2213 466 
Q 2559 466 2825 569 
Q 3091 672 3220 851 
Q 3350 1031 3350 1244 
Q 3350 1459 3225 1620 
Q 3100 1781 2813 1891 
Q 2628 1963 1997 2114 
Q 1366 2266 1113 2400 
Q 784 2572 623 2826 
Q 463 3081 463 3397 
Q 463 3744 659 4045 
Q 856 4347 1234 4503 
Q 1613 4659 2075 4659 
Q 2584 4659 2973 4495 
Q 3363 4331 3572 4012 
Q 3781 3694 3797 3291 
L 3216 3247 
Q 3169 3681 2898 3903 
Q 2628 4125 2100 4125 
Q 1550 4125 1298 3923 
Q 1047 3722 1047 3438 
Q 1047 3191 1225 3031 
Q 1400 2872 2139 2705 
Q 2878 2538 3153 2413 
Q 3553 2228 3743 1945 
Q 3934 1663 3934 1294 
Q 3934 928 3725 604 
Q 3516 281 3123 101 
Q 2731 -78 2241 -78 
Q 1619 -78 1198 103 
Q 778 284 539 648 
Q 300 1013 288 1472 
z
" transform="scale(0.015625)"/>
      </defs>
      <use xlink:href="#ArialMT-52"/>
      <use xlink:href="#ArialMT-65" transform="translate(72.216797 0)"/>
      <use xlink:href="#ArialMT-66" transform="translate(127.832031 0)"/>
      <use xlink:href="#ArialMT-2d" transform="translate(155.615234 0)"/>
      <use xlink:href="#ArialMT-48" transform="translate(188.916016 0)"/>
      <use xlink:href="#ArialMT-47" transform="translate(261.132812 0)"/>
      <use xlink:href="#ArialMT-56" transform="translate(338.916016 0)"/>
      <use xlink:href="#ArialMT-53" transform="translate(405.615234 0)"/>
      <use xlink:href="#ArialMT-20" transform="translate(472.314453 0)"/>
      <use xlink:href="#ArialMT-49" transform="translate(500.097656 0)"/>
      <use xlink:href="#ArialMT-64" transform="translate(527.880859 0)"/>
      <use xlink:href="#ArialMT-65" transform="translate(583.496094 0)"/>
      <use xlink:href="#ArialMT-6e" transform="translate(639.111328 0)"/>
      <use xlink:href="#ArialMT-74" transform="translate(694.726562 0)"/>
      <use xlink:href="#ArialMT-69" transform="translate(722.509766 0)"/>
      <use xlink:href="#ArialMT-74" transform="translate(744.726562 0)"/>
      <use xlink:href="#ArialMT-79" transform="translate(772.509766 0)"/>
     </g>
    </g>
    <g id="patch_27">
     <path d="M 617.011438 129.715391 
L 650.011438 129.715391 
L 650.011438 118.165391 
L 617.011438 118.165391 
z
" style="fill: #7dcea0; stroke: #ffffff; stroke-width: 1.5; stroke-linejoin: miter"/>
    </g>
    <g id="text_19">
     <!-- Ref-HGVS Analogous -->
     <g style="fill: #262626" transform="translate(663.211438 129.715391) scale(0.165 -0.165)">
      <use xlink:href="#ArialMT-52"/>
      <use xlink:href="#ArialMT-65" transform="translate(72.216797 0)"/>
      <use xlink:href="#ArialMT-66" transform="translate(127.832031 0)"/>
      <use xlink:href="#ArialMT-2d" transform="translate(155.615234 0)"/>
      <use xlink:href="#ArialMT-48" transform="translate(188.916016 0)"/>
      <use xlink:href="#ArialMT-47" transform="translate(261.132812 0)"/>
      <use xlink:href="#ArialMT-56" transform="translate(338.916016 0)"/>
      <use xlink:href="#ArialMT-53" transform="translate(405.615234 0)"/>
      <use xlink:href="#ArialMT-20" transform="translate(472.314453 0)"/>
      <use xlink:href="#ArialMT-41" transform="translate(494.597656 0)"/>
      <use xlink:href="#ArialMT-6e" transform="translate(561.296875 0)"/>
      <use xlink:href="#ArialMT-61" transform="translate(616.912109 0)"/>
      <use xlink:href="#ArialMT-6c" transform="translate(672.527344 0)"/>
      <use xlink:href="#ArialMT-6f" transform="translate(694.744141 0)"/>
      <use xlink:href="#ArialMT-67" transform="translate(750.359375 0)"/>
      <use xlink:href="#ArialMT-6f" transform="translate(805.974609 0)"/>
      <use xlink:href="#ArialMT-75" transform="translate(861.589844 0)"/>
      <use xlink:href="#ArialMT-73" transform="translate(917.205078 0)"/>
     </g>
    </g>
   </g>
  </g>
 </g>
 <defs>
  <clipPath id="pc5b34efcd2">
   <rect x="64.030078" y="40.303437" width="535.668" height="314.32"/>
  </clipPath>
 </defs>
</svg>

</p>

<!-- PERFORMANCE_GRAPH_END -->

### Validation Results (100,000 variants)

Summary of results comparing `weaver` and `ref-hgvs` against ClinVar ground truth:

| Implementation | Identity Match | Analogous Match | Total Success |
| :------------- | :------------: | :-------------: | :-----------: |
| **weaver**     |   **93.85%**   |    **4.92%**    |  **98.77%**   |
| ref-hgvs       |     93.34%     |      4.88%      |    98.22%     |

RefSeq Data Mismatches: 0 (0.0%)


#### Comparison with Reference Tool
As shown in the visualization above, `weaver` now consistently outperforms existing implementations in protein projection accuracy and biological equivalence.

- **Variant Equivalence**: Check if two variants are biologically equivalent using advanced cross-coordinate mapping and normalization. [See Algorithm](docs/source/equivalence_logic.md).

## Troubleshooting

### `stub_gen` Linking Error

If you encounter a linking error like `dyld: symbol not found in flat namespace '_PyBool_Type'` when running `stub_gen`, it is because the `extension-module` feature (enabled by default for Python bindings) conflicts with standalone binary execution on some platforms (like macOS).

To fix this, run `stub_gen` with default features disabled:

```bash
cargo run --bin stub_gen --no-default-features
```
