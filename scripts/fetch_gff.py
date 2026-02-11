import argparse
import gzip
import io
import sys
from datetime import UTC, datetime
from xml.etree.ElementTree import Element

import requests
from defusedxml import ElementTree

BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
FTP_BASE = "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/annotation_releases/"

# Thresholds for NCBI Annotation Releases
CURRENT_YEAR_THRESHOLD = 2021
HISTORICAL_YEAR_THRESHOLD = 2014

# Mapping of Annotation Releases to FTP paths
RELEASE_MAP = {
    "109": "109/GCF_000001405.38_GRCh38.p12/GCF_000001405.38_GRCh38.p12_genomic.gff.gz",
    "105": "105.20220307/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_genomic.gff.gz",
}


def fetch_transcript_dates(accession: str) -> tuple[str | None, str | None]:
    """
    Fetches creation and update dates for a transcript accession.
    Returns (create_date_str, update_date_str) or (None, None).
    """
    print(f"Fetching dates for {accession}...", file=sys.stderr)
    try:
        r = requests.get(
            BASE_URL + "esummary.fcgi",
            params={"db": "nucleotide", "id": accession, "retmode": "json", "version": "2.0"},
            timeout=10,
        )
        r.raise_for_status()
        data = r.json()
        if "result" in data:
            for uid, info in data["result"].items():
                if uid == "uids":
                    continue
                return info.get("createdate"), info.get("updatedate")
    except requests.RequestException as e:
        print(f"Warning: Could not fetch dates for {accession}: {e}", file=sys.stderr)
    return None, None


def map_to_annotation_release(update_date_str: str | None) -> str | None:
    """
    Maps an update date to a likely NCBI Annotation Release.
    Date format: YYYY/MM/DD
    """
    if not update_date_str:
        return None

    try:
        dt = datetime.strptime(update_date_str, "%Y/%m/%d").replace(tzinfo=UTC)
        year = dt.year

        # Logic:
        # If updated recently (>2021), use current (return None -> fallback).
        # If updated 2014-2021, try Release 109 (GRCh38).
        # If updated < 2014, try Release 105 (GRCh37) - though user seems to want GRCh38 usually.
        # But 'hgvs_eval.gff' has NC_000001.11 (GRCh38).
        # So we really want "Historical GRCh38".

        if year > CURRENT_YEAR_THRESHOLD:
            return None  # Use current

        if year >= HISTORICAL_YEAR_THRESHOLD:
            return "109"

        return "109"  # Try 109 even for older ones if they existed on GRCh38?
        # Or maybe 105 if they want GRCh37.
        # Given the user's specific case (NM_152486.2, 2018 update), 109 is perfect.

    except ValueError:
        return None


def stream_gff_from_ftp_multi(release_key: str, transcript_ids: list[str]) -> dict[str, list[str]]:
    """
    Streams GFF from NCBI FTP and filters for multiple transcripts.
    Returns dict mapping tid -> matching lines.
    """
    if release_key not in RELEASE_MAP:
        return {}

    file_path = RELEASE_MAP[release_key]
    url = FTP_BASE + file_path

    print(
        f"Streaming GFF from Annotation Release {release_key} for {len(transcript_ids)} transcripts...",
        file=sys.stderr,
    )
    print(f"URL: {url}", file=sys.stderr)

    results: dict[str, list[str]] = {tid: [] for tid in transcript_ids}
    ids_set = set(transcript_ids)

    try:
        with requests.get(url, stream=True, timeout=60) as r:
            r.raise_for_status()
            d = gzip.GzipFile(fileobj=r.raw)
            text_stream = io.TextIOWrapper(d, encoding="utf-8")

            for line in text_stream:
                if line.startswith("#"):
                    continue

                # Substring check for efficiency
                for tid in ids_set:
                    if tid in line:
                        results[tid].append(line.strip())
    except requests.RequestException as e:
        print(f"Error streaming GFF from release {release_key}: {e}", file=sys.stderr)

    return results


# --- Existing Builder for Fallback ---
class RefSeqGFFBuilder:
    def __init__(self, transcript_id: str) -> None:
        self.transcript_id = transcript_id

    def fetch_and_build(self) -> tuple[list[str] | None, str | None]:
        # 1. Search
        gene_id, err = self._get_gene_id()
        if gene_id is None:
            return None, err
        root, err = self._fetch_gene_record(gene_id)
        if root is None:
            return None, err

        selected_node, err = self._select_transcript_node(root)
        if selected_node is None:
            return None, err

        # Extract Coords
        coords_node = selected_node.find("Gene-commentary_genomic-coords")
        loc_node = coords_node.find("Seq-loc") if coords_node is not None else None
        if loc_node is None:
            return None, "No genomic coords or Seq-loc node"
        intervals = self.parse_seq_loc(loc_node)
        if not intervals:
            return None, "Empty intervals"

        # Resolve Chromosome
        gi = intervals[0][3]
        chrom_acc = self._resolve_chromosome_accession(gi)

        # Build GFF Lines
        lines = ["##gff-version 3"]

        min_start = min(i[0] for i in intervals) + 1
        max_end = max(i[1] for i in intervals) + 1
        strand = intervals[0][2]

        attr_parts = [
            f"ID=rna-{self.transcript_id}",
            f"Parent=gene-{gene_id}",
            f"Dbxref=GeneID:{gene_id},Genbank:{self.transcript_id}",
            f"Name={self.transcript_id}",
            "gbkey=mRNA",
            f"gene={gene_id}",
            "product=transcript",
            f"transcript_id={self.transcript_id}",
        ]
        attr_base = ";".join(attr_parts)
        lines.append(f"{chrom_acc}\tBestRefSeq\tmRNA\t{min_start}\t{max_end}\t.\t{strand}\t.\t{attr_base}")

        intervals.sort(key=lambda x: x[0])
        for i, (start, end, st, _) in enumerate(intervals):
            e_attr_parts = [
                f"ID=exon-{self.transcript_id}-{i + 1}",
                f"Parent=rna-{self.transcript_id}",
                f"Dbxref=GeneID:{gene_id},Genbank:{self.transcript_id}",
                "gbkey=mRNA",
                f"gene={gene_id}",
                "product=transcript",
                f"transcript_id={self.transcript_id}",
            ]
            attr = ";".join(e_attr_parts)
            lines.append(f"{chrom_acc}\tBestRefSeq\texon\t{start + 1}\t{end + 1}\t.\t{st}\t.\t{attr}")

        # CDS Logic
        self._append_cds_lines(selected_node, lines, chrom_acc, strand, gene_id)

        return lines, None

    def _get_gene_id(self) -> tuple[str | None, str | None]:
        term = self.transcript_id.split(".")[0]
        print(f"Searching for gene associated with {term}...", file=sys.stderr)
        try:
            r = requests.get(
                BASE_URL + "esearch.fcgi",
                params={"db": "gene", "term": term, "retmode": "json"},
                timeout=10,
            )
            r.raise_for_status()
            data = r.json()
            if not data.get("esearchresult", {}).get("idlist"):
                return None, "Gene not found"
            return data["esearchresult"]["idlist"][0], None
        except requests.RequestException as e:
            return None, f"Search failed: {e}"

    def _fetch_gene_record(self, gene_id: str) -> tuple[Element | None, str | None]:
        print(f"Found Gene ID: {gene_id}. Fetching record...", file=sys.stderr)
        try:
            r = requests.get(
                BASE_URL + "efetch.fcgi",
                params={"db": "gene", "id": gene_id, "retmode": "xml"},
                timeout=20,
            )
            r.raise_for_status()
            return ElementTree.fromstring(r.text), None
        except requests.RequestException as e:
            return None, f"Fetch failed: {e}"
        except ElementTree.ParseError as e:
            return None, f"XML Parse failed: {e}"

    def _select_transcript_node(self, root: Element) -> tuple[Element | None, str | None]:
        target_acc = self.transcript_id.split(".")[0]
        target_ver = self.transcript_id.split(".")[1] if "." in self.transcript_id else None

        candidates = []
        for comm in root.findall(".//Gene-commentary"):
            type_node = comm.find("Gene-commentary_type")
            if type_node is not None and type_node.attrib.get("value") == "mRNA":
                acc_node = comm.find("Gene-commentary_accession")
                if acc_node is None or acc_node.text is None:
                    continue
                acc = acc_node.text
                ver_node = comm.find("Gene-commentary_version")
                ver = ver_node.text if ver_node is not None else None
                if acc == target_acc:
                    candidates.append((comm, ver))

        if not candidates:
            return None, "Transcript not found in Gene record"

        selected_node = None
        if target_ver:
            for node, v in candidates:
                if v == target_ver:
                    selected_node = node
                    break

        if not selected_node:
            print(f"Warning: Version {target_ver} not found. Available: {[v for _, v in candidates]}", file=sys.stderr)
            candidates.sort(key=lambda x: int(x[1]) if x[1] and x[1].isdigit() else 0, reverse=True)
            selected_node = candidates[0][0]
            print(f"Falling back to version {candidates[0][1]}", file=sys.stderr)

        return selected_node, None

    def _resolve_chromosome_accession(self, gi: str) -> str:
        print(f" resolving GI {gi}...", file=sys.stderr)
        chrom_acc = f"gi|{gi}"
        try:
            r = requests.get(
                BASE_URL + "esummary.fcgi",
                params={"db": "nucleotide", "id": gi, "retmode": "json"},
                timeout=10,
            )
            r.raise_for_status()
            res = r.json()
            if "result" in res and str(gi) in res["result"]:
                info = res["result"][str(gi)]
                chrom_acc = info.get("accessionversion", info.get("caption", chrom_acc))
            print(f"Resolved to {chrom_acc}", file=sys.stderr)
        except requests.RequestException as e:
            print(f"Error resolving GI {gi}: {e}", file=sys.stderr)
        return chrom_acc

    def _append_cds_lines(
        self,
        selected_node: Element,
        lines: list[str],
        chrom_acc: str,
        strand: str,
        gene_id: str,
    ) -> None:
        products = selected_node.find("Gene-commentary_products")
        if not products:
            return

        for p in products.findall("Gene-commentary"):
            p_type_node = p.find("Gene-commentary_type")
            if p_type_node is not None and p_type_node.attrib.get("value") == "peptide":
                p_acc_node = p.find("Gene-commentary_accession")
                if p_acc_node is None or p_acc_node.text is None:
                    continue
                prot_acc = p_acc_node.text
                p_ver_node = p.find("Gene-commentary_version")
                prot_ver = p_ver_node.text if p_ver_node is not None else ""
                prot_id = f"{prot_acc}.{prot_ver}" if prot_ver else prot_acc

                p_coords_node = p.find("Gene-commentary_genomic-coords")
                if p_coords_node is not None:
                    loc_node = p_coords_node.find("Seq-loc")
                    if loc_node is not None:
                        cds_ints = self.parse_seq_loc(loc_node)
                        self.append_cds(lines, cds_ints, chrom_acc, strand, prot_id, self.transcript_id, gene_id)

    def parse_seq_loc(self, seq_loc: Element) -> list[tuple[int, int, str, str]]:
        ints: list[tuple[int, int, str, str]] = []
        # Handle mix
        mix = seq_loc.find("Seq-loc_mix/Seq-loc-mix")
        if mix is not None:
            for sl in mix.findall("Seq-loc"):
                ints.extend(self.parse_seq_loc(sl))
            return ints

        # Handle int
        sl_int = seq_loc.find("Seq-loc_int/Seq-interval")
        if sl_int is not None:
            f_node = sl_int.find("Seq-interval_from")
            t_node = sl_int.find("Seq-interval_to")
            if f_node is not None and t_node is not None and f_node.text is not None and t_node.text is not None:
                start = int(f_node.text)
                end = int(t_node.text)
                strand_node = sl_int.find("Seq-interval_strand/Na-strand")
                strand = "-" if strand_node is not None and strand_node.attrib.get("value") == "minus" else "+"
                id_node = sl_int.find("Seq-interval_id/Seq-id/Seq-id_gi")
                gi = id_node.text if id_node is not None and id_node.text is not None else "0"
                return [(start, end, strand, gi)]
        return []

    def append_cds(  # noqa: PLR0913
        self,
        lines: list[str],
        intervals: list[tuple[int, int, str, str]],
        chrom: str,
        strand: str,
        prot_id: str,
        trans_id: str,
        gene_id: str,
    ) -> None:
        if not intervals:
            return
        # Sort
        ordered = sorted(intervals, key=lambda x: x[0], reverse=(strand == "-"))
        raw_intervals = sorted(intervals, key=lambda x: x[0])

        # Phase calc
        phases = {}
        current_phase = 0
        for s, e, _, _ in ordered:
            phases[s] = current_phase
            interval_len = e - s + 1
            rem = (interval_len - current_phase) % 3
            current_phase = (3 - rem) % 3 if rem != 0 else 0

        for s, e, st, _ in raw_intervals:
            ph = phases[s]
            attr = (
                f"ID=cds-{prot_id};Parent=rna-{trans_id};Dbxref=GeneID:{gene_id},Genbank:{prot_id};"
                f"Name={prot_id};gbkey=CDS;gene={gene_id};product=protein;protein_id={prot_id}"
            )
            lines.append(f"{chrom}\tBestRefSeq\tCDS\t{s + 1}\t{e + 1}\t.\t{st}\t{ph}\t{attr}")


def group_transcripts_by_release(
    tids: list[str],
    force_release: str | None = None,
) -> tuple[dict[str, list[str]], list[str]]:
    """Group transcripts by their likely NCBI Annotation Release."""
    release_groups: dict[str, list[str]] = {}
    fallback_ids = []

    for tid in tids:
        if force_release:
            release = force_release
        else:
            _, udate = fetch_transcript_dates(tid)
            release = map_to_annotation_release(udate)

        if release:
            release_groups.setdefault(release, []).append(tid)
        else:
            fallback_ids.append(tid)
    return release_groups, fallback_ids


def process_release_groups(release_groups: dict[str, list[str]], final_output: dict[str, list[str]]) -> list[str]:
    """Fetch GFF lines for each release group."""
    fallback_ids = []
    for release, ids in release_groups.items():
        results = stream_gff_from_ftp_multi(release, ids)
        for tid, lines in results.items():
            if lines and any("\tmRNA\t" in line for line in lines):
                print(f"Successfully retrieved historical GFF lines from Release {release} for {tid}.", file=sys.stderr)
                final_output[tid] = lines
            else:
                print(f"Historical fetch failed for {tid}, will attempt fallback.", file=sys.stderr)
                fallback_ids.append(tid)
    return fallback_ids


def process_fallbacks(fallback_ids: list[str], final_output: dict[str, list[str]]) -> None:
    """Attempt live RefSeq fallback for transcripts not found in historical releases."""
    for tid in fallback_ids:
        if not final_output.get(tid):
            print(f"Attempting live RefSeq fallback for {tid}...", file=sys.stderr)
            builder = RefSeqGFFBuilder(tid)
            lines, err = builder.fetch_and_build()
            if lines:
                # Filter out header
                final_output[tid] = [line for line in lines if not line.startswith("##gff-version")]
            else:
                print(f"Error for {tid}: {err}", file=sys.stderr)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("transcript_ids", nargs="+")
    parser.add_argument("--release", help="Force a specific NCBI Annotation Release (e.g., 105, 109)")
    args = parser.parse_args()

    tids = args.transcript_ids
    final_output: dict[str, list[str]] = {tid: [] for tid in tids}

    # 1. Group by release
    release_groups, initial_fallbacks = group_transcripts_by_release(tids, args.release)

    # 2. Fetch from releases
    more_fallbacks = process_release_groups(release_groups, final_output)

    # 3. Fallback for others
    process_fallbacks(initial_fallbacks + more_fallbacks, final_output)

    # 4. Print Output
    print("##gff-version 3")
    for tid in tids:
        if final_output.get(tid):
            print("\n".join(final_output[tid]))
        else:
            print(f"# Failed to fetch data for {tid}")


if __name__ == "__main__":
    main()
