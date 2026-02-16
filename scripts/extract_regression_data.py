import weaver
from weaver.cli.provider import RefSeqDataProvider


def extract_data(provider, transcript_ac):
    try:
        tx = provider.get_transcript(transcript_ac, None)
    except ValueError:
        print(f"Transcript {transcript_ac} not found")
        return None

    print(f"DEBUG: {transcript_ac} refs: {provider.tx_to_refs.get(transcript_ac)}")
    ref_ac = tx["reference_accession"]
    print(f"DEBUG: {transcript_ac} ref_ac: {ref_ac}")
    print(f"DEBUG: {ref_ac} in FASTA: {ref_ac in provider.fasta.references}")

    seq = provider.get_seq(transcript_ac, 0, -1, "c")
    return {
        "ac": transcript_ac,
        "seq": seq,
        "cds_start": tx["cds_start_index"],
        "cds_end": tx["cds_end_index"],
    }


def main():
    provider = RefSeqDataProvider(
        gff_path="GRCh38_latest_genomic.gff.gz", fasta_path="GCF_000001405.40_GRCh38.p14_genomic.fna"
    )

    variants = ["NM_000038.6", "NM_000527.5", "NM_000478.6", "NM_001122606.1", "NM_000465.4"]

    for ac in variants:
        data = extract_data(provider, ac)
        if data:
            print(f"Transcript: {ac}")
            print(f"CDS Start: {data['cds_start']}")
            print(f"CDS End: {data['cds_end']}")
            # Print a small region around the variant sites if possible,
            # but for tests, we might need more context.
            # I'll just print the full sequence for now as I'll need to
            # potentially trim it for the Rust test.
            print(f"Sequence Length: {len(data['seq'])}")
            print(f"Full Sequence: {data['seq']}")
            print("-" * 40)


if __name__ == "__main__":
    main()
