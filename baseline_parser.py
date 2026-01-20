"""
Baseline Python FASTQ Parser for RBPscan
- R1: Extract 11nt motif from forward reads
- R2: Extract hairpin region and count edits from reverse reads
"""
import os
import re
import time
import gzip
import pandas as pd


def read_fastq(fq_path):
    """Generator that yields (read_id, sequence) from FASTQ file"""
    opener = gzip.open if fq_path.endswith('.gz') else open
    with opener(fq_path, 'rt') as f:
        while True:
            header = f.readline().strip()
            if not header:
                break
            seq = f.readline().strip()
            f.readline()  # + line
            f.readline()  # quality
            read_id = header.split(' ')[0].lstrip('@')
            yield read_id, seq


def parse_fastq(fq_path, parse_motif=True, motif_regex=None, hairpin_regex=None):
    """
    parse_motif=True:  Processes R1 to find the 11nt motif.
    parse_motif=False: Processes R2 to find the Hairpin and Edit Counts.

    Configurable regex patterns:
    - motif_regex: Pattern with capture group for motif (default: TTCTGGCTGACATA(.{11})ATACAATCAGATATGCA)
    - hairpin_regex: Pattern with capture group for hairpin (default: (AA[TC]C[TC]){6}(?=AATTT))
    """
    # Default patterns
    if motif_regex is None:
        motif_regex = r'TTCTGGCTGACATA(.{11})ATACAATCAGATATGCA'
    if hairpin_regex is None:
        hairpin_regex = r"(AA[TC]C[TC]AA[TC]C[TC]AA[TC]C[TC]AA[TC]C[TC]AA[TC]C[TC]AA[TC]C[TC])(?=AATTT)"

    # Pre-compile regex for speed
    if parse_motif:
        pattern = re.compile(motif_regex)
    else:
        pattern = re.compile(hairpin_regex)

    records = []
    file_type = "R1 (Motif)" if parse_motif else "R2 (Hairpin)"

    print(f"Parsing {file_type}: {os.path.basename(fq_path)}")
    start = time.time()

    for read_id, seq in read_fastq(fq_path):
        if parse_motif:
            # R1 LOGIC: 11nt motif
            match = pattern.search(seq)
            motif = match.group(1) if match else None
            records.append({'read_id': read_id, 'motif': motif})
        else:
            # R2 LOGIC: Hairpin (Direct Strand)
            match = pattern.search(seq)
            if match:
                hp_region = match.group(1)
                edits_count = hp_region.count("AATCC")
            else:
                hp_region, edits_count = None, None

            records.append({
                'read_id': read_id,
                'hp': hp_region,
                'edits_count': edits_count
            })

    elapsed = time.time() - start
    df = pd.DataFrame(records)
    print(f"Parsed {len(df)} reads in {elapsed:.2f}s ({len(df)/elapsed:.0f} reads/sec)")
    return df


def merge_paired_reads(fwd_df, rev_df):
    """Merge FWD and REV reads by read_id"""
    merged = pd.merge(fwd_df, rev_df, on='read_id', how='inner')
    merged = merged.dropna(subset=['motif', 'hp', 'edits_count'])
    print(f"Merged: {len(merged)} paired reads with valid motif+hairpin")
    return merged


def benchmark_parser(fq_path, n_reads_list=[10000, 100000, 1000000], parse_motif=True):
    """Benchmark parsing at different read counts"""
    import gzip
    import tempfile

    results = []

    for n_reads in n_reads_list:
        # Create temp file with n_reads
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fq', delete=False) as tmp:
            tmp_path = tmp.name
            with gzip.open(fq_path, 'rt') as f:
                for i, line in enumerate(f):
                    if i >= n_reads * 4:
                        break
                    tmp.write(line)

        # Time the parsing
        start = time.time()
        df = parse_fastq(tmp_path, parse_motif=parse_motif)
        elapsed = time.time() - start

        results.append({
            'n_reads': n_reads,
            'time_sec': elapsed,
            'reads_per_sec': n_reads / elapsed
        })

        os.unlink(tmp_path)
        print(f"  {n_reads:>10,} reads: {elapsed:.3f}s ({n_reads/elapsed:,.0f} reads/sec)")

    return pd.DataFrame(results)


if __name__ == "__main__":
    import sys

    # Test with sample data
    r1_path = "raw_data/Lib2040_1866_2700_rep1_CKDL250033540-1A_2357TCLT4_L6_1.fq.gz"
    r2_path = "raw_data/Lib2040_1866_2700_rep1_CKDL250033540-1A_2357TCLT4_L6_2.fq.gz"

    if len(sys.argv) > 1 and sys.argv[1] == "benchmark":
        print("=== Benchmarking Python Parser ===\n")
        print("R1 (Motif) parsing:")
        r1_results = benchmark_parser(r1_path, parse_motif=True)
        print("\nR2 (Hairpin) parsing:")
        r2_results = benchmark_parser(r2_path, parse_motif=False)
    else:
        # Quick test with 100k reads
        print("=== Testing Python Parser (100k reads) ===\n")

        # Create test files
        import gzip
        import tempfile

        for label, path, is_motif in [("R1", r1_path, True), ("R2", r2_path, False)]:
            with tempfile.NamedTemporaryFile(mode='w', suffix='.fq', delete=False) as tmp:
                with gzip.open(path, 'rt') as f:
                    for i, line in enumerate(f):
                        if i >= 100000 * 4:
                            break
                        tmp.write(line)
                tmp_path = tmp.name

            df = parse_fastq(tmp_path, parse_motif=is_motif)
            print(f"  Non-null: {df.iloc[:, 1].notna().sum()}")
            print()
            os.unlink(tmp_path)
