"RefSeq data provider implementation."

import collections
import gzip
import json
import logging
import os
import sys
import typing

try:
    import hgvs.dataproviders.interface
except ImportError:
    print(
        "Error: 'hgvs' package not found. Please install it manually (e.g. without dependencies to avoid psycopg2) with:",
    )
    print(
        "  pip install hgvs --no-deps",
    )
    sys.exit(1)

try:
    import pysam
except ImportError:
    print("Error: 'pysam' package not found. Please install it manually with: pip install pysam")
    sys.exit(1)


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class IdentifierKind:
    """Enum for sequence identifier types."""

    Genomic = "g"
    Transcript = "c"
    Protein = "p"


class SequenceProxy:
    """Proxy for accessing genomic sequences, supporting recording and replay."""

    def __init__(self, fasta_path: str, cache_path: str | None = None, mode: str | None = None) -> None:
        self.fasta_path = fasta_path
        self.cache_path = cache_path or os.environ.get("WEAVER_SEQ_CACHE")
        self.mode = mode or os.environ.get("WEAVER_SEQ_MODE", "live")
        self.cache: dict[str, str] = {}
        self.fasta: typing.Any = None
        self.references: list[str] = []

        if self.mode == "replay":
            if self.cache_path and os.path.exists(self.cache_path):
                try:
                    with open(self.cache_path) as f:
                        self.cache = json.load(f)

                    # Try to get references from manifest first
                    manifest_data = self.cache.get("_manifest")
                    if isinstance(manifest_data, dict):
                        manifest = typing.cast("dict[str, list[str]]", manifest_data)
                        fname = os.path.basename(self.fasta_path)
                        if fname in manifest:
                            self.references = manifest[fname]
                    else:
                        # Fallback: extract unique accessions from keys like "AC:start-end"
                        refs = set()
                        for key in self.cache:
                            if ":" in key and not key.startswith("_"):
                                refs.add(key.split(":")[0])
                        self.references = list(refs)
                except (OSError, json.JSONDecodeError) as e:
                    logger.error("Failed to load sequence cache from %s: %s", self.cache_path, e)
            else:
                logger.warning("Replay mode enabled but cache file not found: %s", self.cache_path)
        elif os.path.exists(fasta_path):
            try:
                self.fasta = pysam.FastaFile(fasta_path)
                self.references = list(self.fasta.references)
            except Exception as e:
                logger.error("Failed to load FASTA from %s: %s", fasta_path, e)
        elif self.mode != "replay":
            logger.warning("FASTA file not found: %s", fasta_path)

    def fetch(self, ac: str, start: int, end: int | None = None) -> str:
        # Normalize end for consistency in the key (None and -1 should be same)
        key_end = end
        if end == -1:
            key_end = None
        key = f"{ac}:{start}-{key_end}"

        if self.mode == "replay":
            if key in self.cache:
                return self.cache[key]
            logger.error("Missing sequence in cache for %s", key)
            return ""

        if not self.fasta:
            return ""

        try:
            # pysam handles end=-1 or None as end-of-ref
            seq = str(self.fasta.fetch(ac, start, end).upper())
            if self.mode == "record":
                self.cache[key] = seq
            return seq
        except Exception as e:
            logger.error("Error fetching from FASTA for %s: %s", key, e)
            return ""

    def save_cache(self) -> None:
        if self.mode == "record" and self.cache_path:
            try:
                # Merge with existing cache if it exists
                existing = {}
                if os.path.exists(self.cache_path):
                    with open(self.cache_path) as f:
                        existing = json.load(f)

                # Update sequences
                existing.update(self.cache)

                # Update manifest
                manifest = existing.setdefault("_manifest", {})
                fname = os.path.basename(self.fasta_path)
                manifest[fname] = self.references

                os.makedirs(os.path.dirname(self.cache_path), exist_ok=True)
                with open(self.cache_path, "w") as f:
                    json.dump(existing, f, indent=2)
                logger.info("Saved cache with %d sequences and manifest for %s", len(existing) - 1, fname)
            except Exception as e:
                logger.error("Failed to save sequence cache to %s: %s", self.cache_path, e)


class RefSeqDataProvider:
    """DataProvider implementation using RefSeq GFF and FASTA files."""

    def __init__(self, gff_path: str, fasta_path: str) -> None:
        """Initializes the provider and loads the GFF file.

        Args:
          gff_path: Path to the RefSeq GFF3 file (can be gzipped).
          fasta_path: Path to the indexed RefSeq genomic FASTA file.
        """
        self.gff_path = gff_path
        self.fasta_path = fasta_path
        # (tx_ac, chrom) -> TranscriptData
        self.transcripts: dict[tuple[str, str], typing.Any] = {}
        # tx_ac -> list of ref_ac
        self.tx_to_refs: dict[str, list[str]] = collections.defaultdict(list)
        self.gene_to_transcripts: dict[str, list[str]] = collections.defaultdict(list)
        self.chrom_to_transcripts: dict[str, list[typing.Any]] = collections.defaultdict(list)
        self.accession_map: dict[str, tuple[str, str]] = {}  # protein_id -> tx_id

        self._load_gff()
        self.fasta = SequenceProxy(fasta_path)

        self._transcript_cache: dict[tuple[str, str], str] = {}

    def _load_gff(self) -> None:
        """Parses the GFF file into internal transcript models."""
        logger.info("Loading RefSeq GFF from %s into memory...", self.gff_path)

        # (tx_ac, chrom) -> { 'exons': set(), 'cds': set(), 'info': {} }
        tx_data: dict[tuple[str, str], dict[str, typing.Any]] = collections.defaultdict(
            lambda: {"exons": set(), "cds": set(), "info": {}},
        )
        genes: dict[str, str] = {}
        id_to_tx_ac: dict[str, str] = {}

        opener = gzip.open if self.gff_path.endswith(".gz") else open

        with opener(self.gff_path, "rt") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 9:
                    continue

                chrom, source, feature_type, start, end, _score, strand, _frame, attr_str = parts
                attrs: dict[str, str] = {}
                for item in attr_str.split(";"):
                    if "=" in item:
                        k, v = item.split("=", 1)
                        attrs[k] = v

                feat_id = attrs.get("ID")
                parent = attrs.get("Parent")
                tx_ac = attrs.get("transcript_id")

                if not tx_ac:
                    if feat_id and (feat_id.startswith(("rna-NM_", "rna-NR_"))):
                        tx_ac = feat_id[4:]
                    elif parent and (parent.startswith(("rna-NM_", "rna-NR_"))):
                        tx_ac = parent[4:]
                    elif parent in id_to_tx_ac:
                        tx_ac = id_to_tx_ac[parent]

                if feature_type == "gene":
                    if feat_id:
                        genes[feat_id] = attrs.get("gene", attrs.get("Name", ""))

                elif feature_type in [
                    "mRNA",
                    "transcript",
                    "tRNA",
                    "ncRNA",
                    "lnc_RNA",
                    "rRNA",
                    "scRNA",
                    "snRNA",
                    "snoRNA",
                ]:
                    if tx_ac:
                        if feat_id:
                            id_to_tx_ac[feat_id] = tx_ac
                        key = (tx_ac, chrom)
                        if not tx_data[key]["info"] or source in ["RefSeq", "BestRefSeq"]:
                            tx_data[key]["info"] = {
                                "strand": strand,
                                "parent": parent,
                                "gene_name": attrs.get("gene", ""),
                            }

                elif feature_type == "exon":
                    if tx_ac:
                        tx_data[(tx_ac, chrom)]["exons"].add((int(start), int(end)))

                elif feature_type == "CDS" and tx_ac:
                    tx_data[(tx_ac, chrom)]["cds"].add((int(start), int(end)))
                    prot_id = attrs.get("protein_id")
                    if prot_id:
                        tx_data[(tx_ac, chrom)]["info"]["protein_id"] = prot_id

        logger.info("Finalizing %d transcript-reference pairs...", len(tx_data))
        for (tx_id, chrom), data in tx_data.items():
            info = data["info"]
            if not info or "strand" not in info:
                continue

            exons_raw: list[tuple[int, int]] = list(data["exons"])
            if not exons_raw:
                continue

            cds: list[tuple[int, int]] = list(data["cds"])
            strand = info["strand"]
            protein_id = info.get("protein_id", "")

            py_exons: list[dict[str, typing.Any]] = []
            current_tx_pos = 0
            exons_for_tx = sorted(exons_raw, key=lambda x: x[0], reverse=(strand == "-"))

            g_min = min(e[0] for e in exons_raw)
            g_max = max(e[1] for e in exons_raw)

            for i, (e_start, e_end) in enumerate(exons_for_tx):
                exon_len = e_end - e_start + 1
                py_exons.append(
                    {
                        "ord": i,
                        "transcript_start": current_tx_pos,
                        "transcript_end": current_tx_pos + exon_len,
                        "reference_start": e_start - 1,
                        "reference_end": e_end - 1,
                        "alt_strand": 1 if strand == "+" else -1,
                        "cigar": f"{exon_len}=",
                    },
                )
                current_tx_pos += exon_len

            cds_start_idx = None
            cds_end_idx = None
            if cds:
                # Biological start codon base
                g_cds_start = min(c[0] for c in cds) if strand == "+" else max(c[1] for c in cds)
                # Biological last base of stop codon
                g_cds_end = max(c[1] for c in cds) if strand == "+" else min(c[0] for c in cds)

                cds_start_idx = self._genomic_to_tx(g_cds_start, py_exons, "+" if strand == "+" else "-")
                cds_end_idx = self._genomic_to_tx(g_cds_end, py_exons, "+" if strand == "+" else "-")

            gene_name = info["gene_name"]
            if not gene_name and info["parent"]:
                gene_name = genes.get(info["parent"], "")

            record = {
                "ac": tx_id,
                "gene": gene_name,
                "cds_start_index": cds_start_idx,
                "cds_end_index": cds_end_idx,
                "strand": 1 if strand == "+" else -1,
                "reference_accession": chrom,
                "exons": py_exons,
                "protein_id": protein_id,
                "start": g_min,
                "end": g_max,
            }
            self.transcripts[(tx_id, chrom)] = record
            self.tx_to_refs[tx_id].append(chrom)
            self.chrom_to_transcripts[chrom].append(record)

            if protein_id:
                self.accession_map[protein_id] = (tx_id, chrom)
            if gene_name:
                self.gene_to_transcripts[gene_name].append(tx_id)

        logger.info("GFF loading complete.")

    def _genomic_to_tx(self, g_pos: int, exons: list[dict[str, typing.Any]], strand: str) -> int | None:
        """Maps genomic position to transcript index.

        Args:
          g_pos: 1-based genomic coordinate.
          exons: List of exon dictionaries.
          strand: "+" or "-.

        Returns:
          0-based transcript index or None if not in exons.
        """
        g_0 = g_pos - 1
        for exon in exons:
            if exon["reference_start"] <= g_0 <= exon["reference_end"]:
                if strand == "+":
                    return exon["transcript_start"] + (g_0 - exon["reference_start"])
                return exon["transcript_start"] + (exon["reference_end"] - g_0)
        return None

    def get_seq(self, ac: str, start: int, end: int, kind: str, force_plus: bool = False) -> str:
        """Retrieves a sequence from the provider.
        ...
        """
        kind = str(kind).lower()
        if "proteinaccession" in kind or kind in {"p", "protein_accession"}:
            res = self.accession_map.get(ac)
            if not res:
                return ""
            tx_ac, chrom = res
            tx_seq = self._get_tx_seq(tx_ac, chrom, 0, -1, force_plus=force_plus).upper()

            tx_info = self.transcripts.get((tx_ac, chrom))
            if not tx_info or tx_info.get("cds_start_index") is None:
                return ""

            cds_start = tx_info["cds_start_index"]
            cds_end = tx_info["cds_end_index"]

            if cds_start is not None and cds_end is not None:
                cds_seq = tx_seq[cds_start : cds_end + 1]
                translated = self._translate_cds(cds_seq)
                if end == -1 or end is None:
                    return translated[start:]
                return translated[start:end]
            return ""

        if "transcript" in kind.lower() or kind.lower() == "c":
            # Need to decide which reference to use if multiple exist
            refs = self.tx_to_refs.get(ac)
            if not refs:
                return ""
            # Prioritize standard NC chromosomes that exist in the current FASTA
            ref_ac = next((r for r in refs if r.startswith("NC_0000") and r in self.fasta.references), None)
            if not ref_ac:
                ref_ac = next((r for r in refs if r in self.fasta.references), refs[0])

            return self._get_tx_seq(ac, ref_ac, start, end, force_plus=force_plus).upper()

        if "genomic" in kind.lower() or kind == "g":
            pass  # Fall through to fasta.fetch

        try:
            if end == -1 or end is None:
                return str(self.fasta.fetch(ac, start).upper())
            return str(self.fasta.fetch(ac, start, end).upper())
        except Exception as e:
            # Try fuzzy lookup for genomic fetch too?
            # pysam.fetch usually needs exact match.
            logger.error("Error fetching genomic seq for %s: %s", ac, e)
            return ""

    def _get_full_tx_seq(self, tx_ac: str, ref_ac: str) -> str:
        """Builds and caches the full transcript sequence (respecting strand)."""
        tx = self.transcripts.get((tx_ac, ref_ac))
        if not tx:
            return ""
        seq_parts = []
        # Exons are already ordered 5' -> 3' in tx["exons"]
        for exon in tx["exons"]:
            s = self.fasta.fetch(tx["reference_accession"], exon["reference_start"], exon["reference_end"] + 1)
            if tx["strand"] == -1:
                s = self.reverse_complement(s)
            seq_parts.append(s)
        return "".join(seq_parts)

    def _get_full_tx_seq_cached(self, tx_ac: str, ref_ac: str) -> str:
        try:
            return self._transcript_cache[(tx_ac, ref_ac)]
        except KeyError:
            self._transcript_cache[(tx_ac, ref_ac)] = self._get_full_tx_seq(tx_ac, ref_ac)
            return self._transcript_cache[(tx_ac, ref_ac)]

    def _get_tx_seq(self, tx_ac: str, ref_ac: str, start: int, end: int, force_plus: bool = False) -> str:
        """Builds a transcript sequence from genomic exons with optional transformations."""
        tx = self.transcripts.get((tx_ac, ref_ac))
        if not tx:
            return ""

        full_seq = self._get_full_tx_seq_cached(tx_ac, ref_ac)

        if force_plus and tx["strand"] == -1:
            # reverse_complement(natural_seq) gives genomic-plus orientation
            full_seq = self.reverse_complement(full_seq)

        if end == -1 or end is None:
            return full_seq[start:]
        return full_seq[start:end]

    def get_transcript(self, transcript_ac: str, reference_ac: str | None) -> typing.Any:
        """Returns the transcript model for the given accession."""
        if reference_ac:
            tx = self.transcripts.get((transcript_ac, reference_ac))
            if tx:
                return tx

        refs = self.tx_to_refs.get(transcript_ac)

        if not refs:
            raise ValueError(f"Transcript {transcript_ac} not found")

        # Prioritize standard NC chromosomes that exist in the current FASTA
        ref_ac = next((r for r in refs if r.startswith("NC_0000") and r in self.fasta.references), None)
        if not ref_ac:
            ref_ac = next((r for r in refs if r in self.fasta.references), refs[0])

        return self.transcripts[(transcript_ac, ref_ac)]

    def get_transcripts_for_region(self, chrom: str, start: int, end: int) -> list[str]:
        """Finds transcripts overlapping a genomic region."""
        transcripts = self.chrom_to_transcripts.get(chrom, [])
        results = []
        for tx in transcripts:
            # Check overlap: (start1 <= end2) and (end1 >= start2)
            if start <= tx["end"] and end >= tx["start"]:
                results.append(tx["ac"])
        return list(set(results))

    def get_symbol_accessions(self, symbol: str, source_kind: str, target_kind: str) -> list[tuple[str, str]]:
        """Maps gene symbols to transcript accessions."""
        if source_kind == IdentifierKind.Transcript and target_kind == IdentifierKind.Protein:
            # Ambiguous if multiple refs, but usually protein is same
            refs = self.tx_to_refs.get(symbol)
            if refs:
                tx = self.transcripts.get((symbol, refs[0]))
                if tx and tx.get("protein_id"):
                    return [("protein_accession", tx["protein_id"])]
        if symbol in self.gene_to_transcripts:
            return [("transcript_accession", tx_ac) for tx_ac in self.gene_to_transcripts[symbol]]
        return [("gene_symbol", symbol)]

    def get_identifier_type(self, identifier: str) -> str:
        """Determines the type of the identifier."""
        if identifier.startswith(("NC_", "NW_", "NT_")):
            return "genomic_accession"
        if identifier.startswith(("NM_", "NR_", "XM_", "XR_")):
            return "transcript_accession"
        if identifier.startswith(("NP_", "XP_")):
            return "protein_accession"
        # If it's a known gene symbol
        if identifier in self.gene_to_transcripts:
            return "gene_symbol"
        # Fallback heuristic: simplistic assumption that anything else might be a gene symbol
        # if it doesn't contain a colon (which would imply pre-parsed variant)
        if ":" not in identifier:
            return "gene_symbol"
        return "unknown"

    def reverse_complement(self, seq: str) -> str:
        """Returns the reverse complement of a sequence."""
        complement = str.maketrans("ATCGNautcgn", "TAGCNtagcgn")
        return seq.translate(complement)[::-1]

    def _translate_cds(self, seq: str) -> str:
        """Translates a CDS nucleotide sequence into an amino acid sequence (1-letter codes)."""
        codon_table = {
            "ATA": "I",
            "ATC": "I",
            "ATT": "I",
            "ATG": "M",
            "ACA": "T",
            "ACC": "T",
            "ACG": "T",
            "ACT": "T",
            "AAC": "N",
            "AAT": "N",
            "AAA": "K",
            "AAG": "K",
            "AGC": "S",
            "AGT": "S",
            "AGA": "R",
            "AGG": "R",
            "CTA": "L",
            "CTC": "L",
            "CTG": "L",
            "CTT": "L",
            "CCA": "P",
            "CCC": "P",
            "CCG": "P",
            "CCT": "P",
            "CAC": "H",
            "CAT": "H",
            "CAA": "Q",
            "CAG": "Q",
            "CGA": "R",
            "CGC": "R",
            "CGG": "R",
            "CGT": "R",
            "GTA": "V",
            "GTC": "V",
            "GTG": "V",
            "GTT": "V",
            "GCA": "A",
            "GCC": "A",
            "GCG": "A",
            "GCT": "A",
            "GAC": "D",
            "GAT": "D",
            "GAA": "E",
            "GAG": "E",
            "GGA": "G",
            "GGC": "G",
            "GGG": "G",
            "GGT": "G",
            "TCA": "S",
            "TCC": "S",
            "TCG": "S",
            "TCT": "S",
            "TTC": "F",
            "TTT": "F",
            "TTA": "L",
            "TTG": "L",
            "TAC": "Y",
            "TAT": "Y",
            "TAA": "*",
            "TAG": "*",
            "TGC": "C",
            "TGT": "C",
            "TGA": "*",
            "TGG": "W",
        }
        protein = []
        seq = seq.upper()
        for i in range(0, len(seq) - len(seq) % 3, 3):
            codon = seq[i : i + 3]
            protein.append(codon_table.get(codon, "X"))
        return "".join(protein)

    def c_to_g(self, transcript_ac: str, pos: int, offset: int) -> tuple[str, int]:
        """Resolves a transcript position and offset to a genomic accession and position.

        Args:
          transcript_ac: The transcript accession.
          pos: 0-based transcript index.
          offset: intronic offset.

        Returns:
          A tuple of (genomic_accession, 0-based_genomic_position).
        """
        tx = self.get_transcript(transcript_ac, None)
        chrom = tx["reference_accession"]
        strand = tx["strand"]
        exons = tx["exons"]

        # This logic should match _genomic_to_tx but in reverse.
        # Transcript position 'pos' is relative to the start of the transcript sequence.
        # We need to find the exon containing this position.
        for exon in exons:
            exon["reference_end"] - exon["reference_start"] + 1
            if exon["transcript_start"] <= pos < exon["transcript_end"]:
                # Found the exon
                offset_in_exon = pos - exon["transcript_start"]
                if strand == 1:
                    g_base = exon["reference_start"] + offset_in_exon
                    g_pos = g_base + offset
                else:
                    g_base = exon["reference_end"] - offset_in_exon
                    g_pos = g_base - offset
                return chrom, g_pos

        # If not in exons (should not happen for valid cDNA variants without offset,
        # but for offset calculation we might be at exon boundary)
        # Handle cases where pos might be exactly at boundary for intronic mapping
        if pos < 0:
            # Before first exon
            exon = exons[0]
            if strand == 1:
                return chrom, exon["reference_start"] + pos + offset
            return chrom, exon["reference_end"] - pos - offset

        # After last exon or in intron
        # Just use the last exon for default projection if not found
        exon = exons[-1]
        if pos >= exon["transcript_end"]:
            delta = pos - (exon["transcript_end"] - 1)
            if strand == 1:
                return chrom, exon["reference_end"] + delta + offset
            return chrom, exon["reference_start"] - delta - offset

        # If it's a valid transcript position but 'exons' loop failed?
        return chrom, 0  # Fallback

    def to_json(self) -> None:
        return None


class ReferenceHgvsDataProvider(hgvs.dataproviders.interface.Interface):
    """Bridge between weaver DataProvider and hgvs library Interface."""

    def __init__(self, refseq_provider: RefSeqDataProvider) -> None:
        self.url = "local://refseq"
        self.required_version = "1.1"
        super().__init__()
        self.rp = refseq_provider

    def get_seq(self, ac: str, start: int | None = None, end: int | None = None) -> str:
        kind = "g" if ac.startswith("NC_") else "c"
        if ac.startswith("NP_"):
            kind = "p"
        return self.rp.get_seq(ac, start or 0, end or -1, kind)

    def get_tx_info(
        self,
        tx_ac: str,
        alt_ac: str | None = None,
        _alt_aln_method: str | None = None,
    ) -> dict[str, typing.Any] | None:
        try:
            tx = self.rp.get_transcript(tx_ac, alt_ac)
            return {
                "hgnc": tx["gene"],
                "cds_start_i": tx["cds_start_index"],
                "cds_end_i": tx["cds_end_index"] + 1 if tx["cds_end_index"] is not None else None,
                "strand": tx["strand"],
                "alt_ac": tx["reference_accession"],
                "alt_aln_method": "transcript",
            }
        except Exception:
            return None

    def get_tx_exons(
        self,
        tx_ac: str,
        alt_ac: str | None = None,
        _alt_aln_method: str | None = None,
    ) -> list[dict[str, typing.Any]] | None:
        try:
            tx = self.rp.get_transcript(tx_ac, alt_ac)
            res = []
            exons_transcript = sorted(tx["exons"], key=lambda x: x["transcript_start"])
            for i, e in enumerate(exons_transcript):
                res.append(
                    {
                        "tx_ac": tx_ac,
                        "alt_ac": tx["reference_accession"],
                        "alt_aln_method": "transcript",
                        "ord": i,
                        "tx_start_i": e["transcript_start"],
                        "tx_end_i": e["transcript_end"],
                        "alt_start_i": e["reference_start"],
                        "alt_end_i": e["reference_end"] + 1,
                        "alt_strand": tx["strand"],
                        "cigar": e["cigar"],
                    },
                )
            res.sort(key=lambda x: x["alt_start_i"])
            return res
        except Exception:
            return None

    def get_tx_identity_info(self, tx_ac: str) -> dict[str, typing.Any] | None:
        try:
            tx = self.rp.get_transcript(tx_ac, None)
            total_len = tx["exons"][-1]["transcript_end"]
            return {
                "hgnc": tx["gene"],
                "lengths": [total_len],
                "tx_ac": tx_ac,
                "alt_acs": [tx["reference_accession"]],
                "cds_start_i": tx["cds_start_index"],
                "cds_end_i": tx["cds_end_index"] + 1 if tx["cds_end_index"] is not None else None,
            }
        except Exception:
            return None

    def data_version(self) -> str:
        return "1.1"

    def schema_version(self) -> str:
        return "1.1"

    def get_assembly_map(self, _assembly_name: str) -> dict[str, typing.Any]:
        return {}

    def get_gene_info(self, _gene: str) -> dict[str, typing.Any]:
        return {}

    def get_pro_ac_for_tx_ac(self, tx_ac: str) -> str | None:
        try:
            tx = self.rp.get_transcript(tx_ac, None)
            return str(tx.get("protein_id"))
        except Exception:
            return None

    def get_tx_for_gene(self, _gene: str) -> list[str]:
        return []

    def get_tx_for_region(self, _alt_ac: str, _alt_aln_method: str, _start: int, _end: int) -> list[str]:
        return []

    def get_tx_mapping_options(self, tx_ac: str) -> list[dict[str, typing.Any]]:
        try:
            refs = self.rp.tx_to_refs.get(tx_ac, [])
            return [{"tx_ac": tx_ac, "alt_ac": r, "alt_aln_method": "transcript"} for r in refs]
        except Exception:
            return []

    def get_similar_transcripts(self, _tx_ac: str) -> list[str]:
        return []

    def get_acs_for_protein_seq(self, _seq: str) -> list[str]:
        return []
