#!/usr/bin/env python3
"""
Benchmark: Python vs Rust FASTQ Parser
Compares performance at 10k, 100k, and 1M reads
"""
import os
import subprocess
import tempfile
import time
import gzip
import re
import pandas as pd
import matplotlib.pyplot as plt

# Paths
R1_PATH = "raw_data/Lib2040_1866_2700_rep1_CKDL250033540-1A_2357TCLT4_L6_1.fq.gz"
R2_PATH = "raw_data/Lib2040_1866_2700_rep1_CKDL250033540-1A_2357TCLT4_L6_2.fq.gz"
RUST_BIN = "fastq_parser_rs/target/release/fastq_parser"

# Patterns
MOTIF_PATTERN = re.compile(r'TTCTGGCTGACATA(.{11})ATACAATCAGATATGCA')
HAIRPIN_PATTERN = re.compile(r"(AA[TC]C[TC]AA[TC]C[TC]AA[TC]C[TC]AA[TC]C[TC]AA[TC]C[TC]AA[TC]C[TC])AATTT")


def create_test_file(src_path, n_reads, suffix='.fq'):
    """Extract n_reads from gzipped FASTQ to temp file"""
    with tempfile.NamedTemporaryFile(mode='w', suffix=suffix, delete=False) as tmp:
        with gzip.open(src_path, 'rt') as f:
            for i, line in enumerate(f):
                if i >= n_reads * 4:
                    break
                tmp.write(line)
        return tmp.name


def python_parse_motif(fq_path):
    """Python baseline: parse motifs from FASTQ"""
    records = []
    with open(fq_path, 'r') as f:
        while True:
            header = f.readline().strip()
            if not header:
                break
            seq = f.readline().strip()
            f.readline()
            f.readline()
            read_id = header.split(' ')[0].lstrip('@')
            match = MOTIF_PATTERN.search(seq)
            motif = match.group(1) if match else None
            records.append((read_id, motif))
    return records


def python_parse_hairpin(fq_path):
    """Python baseline: parse hairpins from FASTQ"""
    records = []
    with open(fq_path, 'r') as f:
        while True:
            header = f.readline().strip()
            if not header:
                break
            seq = f.readline().strip()
            f.readline()
            f.readline()
            read_id = header.split(' ')[0].lstrip('@')
            match = HAIRPIN_PATTERN.search(seq)
            if match:
                hp = match.group(1)
                edits = hp.count("AATCC")
            else:
                hp, edits = None, None
            records.append((read_id, hp, edits))
    return records


def rust_parse(fq_path, mode):
    """Run Rust parser"""
    result = subprocess.run(
        [RUST_BIN, '--input', fq_path, '--mode', mode, '--output', '/dev/null'],
        capture_output=True, text=True
    )
    return result.stderr


def benchmark_single(n_reads, mode='motif'):
    """Benchmark Python vs Rust for a given read count and mode"""
    src = R1_PATH if mode == 'motif' else R2_PATH
    parse_fn = python_parse_motif if mode == 'motif' else python_parse_hairpin

    print(f"\n{'='*60}")
    print(f"Benchmarking {mode.upper()} mode with {n_reads:,} reads")
    print('='*60)

    # Create test file
    tmp_path = create_test_file(src, n_reads)

    try:
        # Python benchmark
        start = time.time()
        py_result = parse_fn(tmp_path)
        py_time = time.time() - start
        py_matches = sum(1 for r in py_result if r[1] is not None)

        # Rust benchmark
        start = time.time()
        rust_output = rust_parse(tmp_path, mode)
        rust_time = time.time() - start

        # Parse Rust output for match count
        rust_matches = 0
        if 'Matches:' in rust_output:
            rust_matches = int(rust_output.split('Matches:')[1].split()[0])

        # Results
        speedup = py_time / rust_time if rust_time > 0 else float('inf')

        print(f"\nPython:  {py_time:.3f}s ({n_reads/py_time:,.0f} reads/sec) - {py_matches:,} matches")
        print(f"Rust:    {rust_time:.3f}s ({n_reads/rust_time:,.0f} reads/sec) - {rust_matches:,} matches")
        print(f"Speedup: {speedup:.1f}x faster")

        return {
            'n_reads': n_reads,
            'mode': mode,
            'python_time': py_time,
            'python_reads_per_sec': n_reads / py_time,
            'rust_time': rust_time,
            'rust_reads_per_sec': n_reads / rust_time,
            'speedup': speedup,
            'python_matches': py_matches,
            'rust_matches': rust_matches
        }
    finally:
        os.unlink(tmp_path)


def run_all_benchmarks():
    """Run benchmarks for all configurations"""
    read_counts = [10_000, 100_000, 1_000_000]
    modes = ['motif', 'hairpin']

    results = []
    for mode in modes:
        for n_reads in read_counts:
            result = benchmark_single(n_reads, mode)
            results.append(result)

    return pd.DataFrame(results)


def plot_results(df):
    """Create benchmark visualization"""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    for idx, mode in enumerate(['motif', 'hairpin']):
        ax = axes[idx]
        mode_df = df[df['mode'] == mode]

        x = range(len(mode_df))
        width = 0.35

        # Bar chart
        bars1 = ax.bar([i - width/2 for i in x], mode_df['python_reads_per_sec'] / 1e6,
                       width, label='Python', color='#3498db', alpha=0.8)
        bars2 = ax.bar([i + width/2 for i in x], mode_df['rust_reads_per_sec'] / 1e6,
                       width, label='Rust', color='#e74c3c', alpha=0.8)

        ax.set_xlabel('Number of Reads')
        ax.set_ylabel('Throughput (Million reads/sec)')
        ax.set_title(f'{mode.upper()} Mode Performance')
        ax.set_xticks(x)
        ax.set_xticklabels([f'{n:,}' for n in mode_df['n_reads']])
        ax.legend()
        ax.grid(axis='y', alpha=0.3)

        # Add speedup labels
        for i, row in enumerate(mode_df.itertuples()):
            ax.annotate(f'{row.speedup:.1f}x',
                       xy=(i + width/2, row.rust_reads_per_sec / 1e6),
                       xytext=(0, 5), textcoords='offset points',
                       ha='center', fontsize=9, fontweight='bold')

    plt.tight_layout()
    plt.savefig('benchmark_results.png', dpi=150, bbox_inches='tight')
    print(f"\nPlot saved to: benchmark_results.png")
    plt.show()


if __name__ == '__main__':
    print("FASTQ Parser Benchmark: Python vs Rust")
    print("=" * 60)

    # Check Rust binary exists
    if not os.path.exists(RUST_BIN):
        print(f"Error: Rust binary not found at {RUST_BIN}")
        print("Run: cd fastq_parser_rs && cargo build --release")
        exit(1)

    # Run benchmarks
    results_df = run_all_benchmarks()

    # Summary table
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(results_df[['mode', 'n_reads', 'python_time', 'rust_time', 'speedup']].to_string(index=False))

    # Save results
    results_df.to_csv('benchmark_results.csv', index=False)
    print(f"\nResults saved to: benchmark_results.csv")

    # Plot if matplotlib available
    try:
        plot_results(results_df)
    except Exception as e:
        print(f"Could not create plot: {e}")
