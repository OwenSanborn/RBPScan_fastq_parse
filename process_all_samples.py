"""
Process all RBPscan samples.
1. Parse R1 (motif) and R2 (hairpin)
2. Merge paired reads
3. Compute editing scores per motif
4. Output summarized CSV
"""
import os
import sys
from pathlib import Path
import pandas as pd
import numpy as np
from fastq_parser import parse_motif, parse_hairpin

# ============ CONFIGURE THESE PATHS ============
INPUT_DIR = "/path/to/01.RawData/RBPscan"  # Directory containing sample folders
OUTPUT_DIR = "/path/to/out"                 # Where to save results
# ===============================================


def editing_counts(reads: pd.DataFrame, min_occurrence: int = 10, sites: int = 6) -> pd.DataFrame:
    """
    Compute editing scores and efficiencies per motif per sample.

    Returns DataFrame with columns:
        sample, motif, occurrence, reads_edited, total_edits, all_reads, score, efficiency_1
    """
    cols = ['sample', 'motif', 'edits_count']
    if reads is None or reads.empty or not set(cols).issubset(reads.columns):
        return pd.DataFrame(columns=[
            'sample', 'motif', 'occurrence', 'reads_edited',
            'total_edits', 'all_reads', 'score', 'efficiency_1'
        ])

    # Keep only needed columns, drop NaN
    df = reads[cols].dropna().copy()
    df['edits_count'] = pd.to_numeric(df['edits_count'], errors='coerce').fillna(0).astype(np.int16)
    df['sample'] = df['sample'].astype('category')
    df['motif'] = df['motif'].astype('category')

    # Per-sample total reads
    sample_sizes = df.groupby('sample', observed=True).size()

    # Group by sample + motif and aggregate
    g = df.groupby(['sample', 'motif'], observed=True)['edits_count']
    out = g.agg(
        occurrence='size',
        total_edits='sum',
        reads_edited=lambda s: np.count_nonzero(s.values > 0)
    ).reset_index()

    # Filter low-support motifs
    out = out[out['occurrence'] >= min_occurrence]

    # Add derived metrics
    out['all_reads'] = out['sample'].map(sample_sizes)
    out['score'] = out['total_edits'] / out['occurrence']
    out['efficiency_1'] = out['total_edits'] / (sites * out['occurrence'] + 1)

    # Sort
    out = out.sort_values(
        ['sample', 'occurrence', 'total_edits'],
        ascending=[True, False, False]
    ).reset_index(drop=True)

    return out


def find_fastq_files(sample_dir):
    """Find R1 and R2 files in a sample directory."""
    r1_files = list(Path(sample_dir).glob("*_1.fq.gz"))
    r2_files = list(Path(sample_dir).glob("*_2.fq.gz"))
    if not r1_files or not r2_files:
        return None, None
    return str(r1_files[0]), str(r2_files[0])


def process_sample(sample_dir, output_dir):
    """Process a single sample: parse, merge, summarize."""
    sample_name = Path(sample_dir).name
    r1_path, r2_path = find_fastq_files(sample_dir)

    if not r1_path or not r2_path:
        print(f"  Skipping {sample_name}: Missing R1 or R2 files")
        return None

    print(f"\n{'='*60}")
    print(f"Processing: {sample_name}")
    print('='*60)

    # 1. Parse R1 (motifs)
    print("\n[1/4] Parsing R1 (motifs)...")
    motif_df = parse_motif(r1_path)

    # 2. Parse R2 (hairpins)
    print("\n[2/4] Parsing R2 (hairpins)...")
    hairpin_df = parse_hairpin(r2_path)

    # 3. Merge by read_id
    print("\n[3/4] Merging paired reads...")
    merged = pd.merge(motif_df, hairpin_df, on='read_id', how='inner')
    merged = merged.dropna(subset=['motif', 'hp', 'edits_count'])
    merged['sample'] = sample_name

    print(f"  R1 reads: {len(motif_df):,}")
    print(f"  R2 reads: {len(hairpin_df):,}")
    print(f"  Merged (valid motif + hairpin): {len(merged):,}")

    if merged.empty:
        print(f"  WARNING: No valid paired reads for {sample_name}")
        return None

    # 4. Compute editing scores
    print("\n[4/4] Computing editing scores...")
    edits_summary = editing_counts(merged, min_occurrence=10)

    print(f"  Unique motifs (>= 10 reads): {len(edits_summary):,}")

    # Save individual sample summary
    out_path = Path(output_dir) / f"{sample_name}_edits.csv"
    edits_summary.to_csv(out_path, index=False)
    print(f"  Saved: {out_path}")

    return edits_summary


def process_all(input_dir, output_dir):
    """Process all sample directories."""
    input_path = Path(input_dir)
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Find all sample directories
    sample_dirs = sorted([d for d in input_path.iterdir() if d.is_dir()])
    print(f"Found {len(sample_dirs)} samples to process")

    all_results = []

    for sample_dir in sample_dirs:
        try:
            result = process_sample(sample_dir, output_path)
            if result is not None and not result.empty:
                all_results.append(result)
        except Exception as e:
            print(f"  ERROR processing {sample_dir.name}: {e}")
            import traceback
            traceback.print_exc()
            continue

    # Combine all results
    if all_results:
        print(f"\n{'='*60}")
        print("FINAL SUMMARY")
        print('='*60)

        combined = pd.concat(all_results, ignore_index=True)
        combined_path = output_path / "final_consolidated_edits.csv"
        combined.to_csv(combined_path, index=False)

        print(f"\nSamples processed: {len(all_results)}")
        print(f"Total motif-sample pairs: {len(combined):,}")
        print(f"\nSaved: {combined_path}")

        # Quick stats
        print("\nPer-sample summary:")
        summary = combined.groupby('sample').agg({
            'motif': 'count',
            'score': 'mean'
        }).rename(columns={'motif': 'n_motifs', 'score': 'avg_score'})
        print(summary.to_string())

        return combined

    return None


if __name__ == "__main__":
    if len(sys.argv) >= 3:
        INPUT_DIR = sys.argv[1]
        OUTPUT_DIR = sys.argv[2]

    if INPUT_DIR == "/path/to/01.RawData/RBPscan":
        print("Usage: python process_all_samples.py <input_dir> <output_dir>")
        print("\nExample:")
        print("  python process_all_samples.py /path/to/01.RawData/RBPscan /path/to/out")
        sys.exit(1)

    process_all(INPUT_DIR, OUTPUT_DIR)
