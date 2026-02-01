"RefSeq data provider implementation."

import collections
import gzip
import logging
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
        self.fasta = pysam.FastaFile(fasta_path)

    def _load_gff(self) -> None:
        """Parses the GFF file into internal transcript models."""
        logger.info("Loading RefSeq GFF from %s into memory...", self.gff_path)

        # (tx_ac, chrom) -> { 'exons': set(), 'cds': set(), 'info': {} }
        tx_data: dict[tuple[str, str], dict[str, typing.Any]] = collections.defaultdict(
            lambda: {"exons": set(), "cds": set(), "info": {}},
        )
        genes: dict[str, str] = {}

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

        Args:
          ac: Accession to fetch.
          start: 0-based interbase start.
          end: 0-based interbase end (-1 for end of sequence).
          kind: "g", "c", or "p".
          force_plus: If True, always returns the plus-strand genomic sequence.

        Returns:
          The requested sequence string.
        """
        if kind == "p":
            res = self.accession_map.get(ac)
            if not res:
                return ""
            tx_ac, chrom = res
            return self._get_tx_seq(tx_ac, chrom, start, end, force_plus=force_plus).upper()

        if kind == "c":
            # Need to decide which reference to use if multiple exist
            refs = self.tx_to_refs.get(ac)
            if not refs:
                return ""
            # Prioritize standard NC chromosomes
            ref_ac = next((r for r in refs if r.startswith("NC_0000")), refs[0])
            return self._get_tx_seq(ac, ref_ac, start, end, force_plus=force_plus).upper()

        try:
            if end == -1 or end is None:
                return str(self.fasta.fetch(ac, start).upper())
            return str(self.fasta.fetch(ac, start, end).upper())
        except Exception as e:
            logger.error("Error fetching genomic seq for %s: %s", ac, e)
            return ""

    def _get_tx_seq(self, tx_ac: str, ref_ac: str, start: int, end: int, force_plus: bool = False) -> str:
        """Builds a transcript sequence from genomic exons."""
        tx = self.transcripts.get((tx_ac, ref_ac))
        if not tx:
            return ""
        seq_parts = []
        exons = tx["exons"]
        if force_plus:
            # Sort by genomic coordinate ascending
            exons = sorted(exons, key=lambda x: x["reference_start"])

        for exon in exons:
            s = self.fasta.fetch(tx["reference_accession"], exon["reference_start"], exon["reference_end"] + 1)
            if not force_plus and tx["strand"] == -1:
                s = self.reverse_complement(s)
            seq_parts.append(s)
        full_seq = "".join(seq_parts)
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

        ref_ac = next((r for r in refs if r.startswith("NC_0000")), refs[0])
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

    def get_symbol_accessions(self, symbol: str, source_kind: str, target_kind: str) -> list[str]:
        """Maps gene symbols to transcript accessions."""
        if source_kind == IdentifierKind.Transcript and target_kind == IdentifierKind.Protein:
            # Ambiguous if multiple refs, but usually protein is same
            refs = self.tx_to_refs.get(symbol)
            if refs:
                tx = self.transcripts.get((symbol, refs[0]))
                if tx and tx.get("protein_id"):
                    return [tx["protein_id"]]
        if symbol in self.gene_to_transcripts:
            return self.gene_to_transcripts[symbol]
        return [symbol]

    def reverse_complement(self, seq: str) -> str:
        """Returns the reverse complement of a sequence."""
        complement = str.maketrans("ATCGNautcgn", "TAGCNtagcgn")
        return seq.translate(complement)[::-1]

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
