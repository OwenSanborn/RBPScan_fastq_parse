"""
Process all RBPscan samples.
1. Parse R1 (motif) and R2 (hairpin) with fast Rust parser
2. Merge paired reads
3. Summarize by motif with editing_counts()
4. Output: sample, motif, occurrence, total_edits, reads_edited, all_reads, score, efficiency_1
"""
import argparse
import pandas as pd
import numpy as np
from pathlib import Path
from fastq_parser import parse_paired


PAIR_PATTERNS = [
    ("*_R1_001.fastq.gz", "*_R2_001.fastq.gz"),
    ("*_R1.fastq.gz",     "*_R2.fastq.gz"),
    ("*_R1_001.fq.gz",    "*_R2_001.fq.gz"),
    ("*_R1.fq.gz",        "*_R2.fq.gz"),
    ("*_1.fastq.gz",      "*_2.fastq.gz"),
    ("*_1.fq.gz",         "*_2.fq.gz"),
    ("*_R1.fq",           "*_R2.fq"),
    ("*_1.fq",            "*_2.fq"),
]


def find_paired_files(folder):
    """
    Find R1/R2 paired files in a folder using multiple naming conventions.

    Returns:
        (Path, Path) if found, or (None, None) if no paired files found.
    """
    folder = Path(folder)
    for pat_r1, pat_r2 in PAIR_PATTERNS:
        r1_files = list(folder.glob(pat_r1))
        r2_files = list(folder.glob(pat_r2))
        if r1_files and r2_files:
            if len(r1_files) > 1:
                print(f"  Warning: multiple R1 files match '{pat_r1}' in {folder.name}, using first")
            if len(r2_files) > 1:
                print(f"  Warning: multiple R2 files match '{pat_r2}' in {folder.name}, using first")
            return r1_files[0], r2_files[0]
    return None, None


def editing_counts(reads, min_occurrence=10, sites=6):
    """
    Summarize edits by motif.

    Input: DataFrame with columns [sample, motif, edits_count]
    Output: DataFrame with columns [sample, motif, occurrence, total_edits,
            reads_edited, all_reads, score, efficiency_1]
    """
    df = reads[['sample', 'motif', 'edits_count']].dropna().copy()
    df['edits_count'] = pd.to_numeric(df['edits_count'], errors='coerce').fillna(0).astype(int)

    sample_sizes = df.groupby('sample').size()

    out = df.groupby(['sample', 'motif']).agg(
        occurrence=('edits_count', 'size'),
        total_edits=('edits_count', 'sum'),
        reads_edited=('edits_count', lambda x: (x > 0).sum())
    ).reset_index()

    out = out[out['occurrence'] >= min_occurrence]
    out['all_reads'] = out['sample'].map(sample_sizes)
    out['score'] = out['total_edits'] / out['occurrence']
    out['efficiency_1'] = out['total_edits'] / (sites * out['occurrence'] + 1)

    return out.sort_values(['sample', 'occurrence'], ascending=[True, False]).reset_index(drop=True)


def apply_library_key(
    summary_df: pd.DataFrame,
    library_key_path: str,
    sequence_col: str = "sequence",
    name_col: str = "name",
    sites: int = 6,
) -> pd.DataFrame:
    """
    Map observed motifs to known library variants by substring match.

    Args:
        summary_df: Output of editing_counts() with columns [sample, motif, ...]
        library_key_path: Path to CSV with known variants (columns: name, sequence)
        sequence_col: Column name in key CSV containing variant sequences
        name_col: Column name in key CSV containing variant names
        sites: Number of editing sites (used for efficiency_1 calculation)

    Returns:
        DataFrame with columns: sample, variant_name, variant_seq, occurrence,
        total_edits, reads_edited, all_reads, score, efficiency_1
        One row per (sample, variant), including 0-count rows for unobserved variants.
    """
    key_df = pd.read_csv(library_key_path)
    if sequence_col not in key_df.columns:
        raise ValueError(f"Column '{sequence_col}' not found in library key. Available: {list(key_df.columns)}")
    if name_col not in key_df.columns:
        raise ValueError(f"Column '{name_col}' not found in library key. Available: {list(key_df.columns)}")

    rows = []
    for sample in summary_df['sample'].unique():
        sample_df = summary_df[summary_df['sample'] == sample]
        all_reads = int(sample_df['all_reads'].iloc[0]) if len(sample_df) > 0 else 0

        for _, variant_row in key_df.iterrows():
            variant_seq = variant_row[sequence_col]
            variant_name = variant_row[name_col]

            matches = sample_df[sample_df['motif'].str.contains(variant_seq, regex=False, na=False)]

            if len(matches) > 0:
                occurrence = int(matches['occurrence'].sum())
                total_edits = int(matches['total_edits'].sum())
                reads_edited = int(matches['reads_edited'].sum())
                score = total_edits / occurrence if occurrence > 0 else 0.0
                efficiency_1 = total_edits / (sites * occurrence + 1) if occurrence > 0 else 0.0
            else:
                occurrence = 0
                total_edits = 0
                reads_edited = 0
                score = 0.0
                efficiency_1 = 0.0

            rows.append({
                'sample': sample,
                'variant_name': variant_name,
                'variant_seq': variant_seq,
                'occurrence': occurrence,
                'total_edits': total_edits,
                'reads_edited': reads_edited,
                'all_reads': all_reads,
                'score': score,
                'efficiency_1': efficiency_1,
            })

    return pd.DataFrame(rows)


def process_all(input_dir, output_dir, min_occurrence=10, library_key=None):
    """
    Process all sample directories.

    Args:
        input_dir: Directory containing sample folders (each with paired FASTQ files)
                   or a flat directory of paired FASTQ files
        output_dir: Where to save results
        min_occurrence: Minimum read count per motif to include in output
        library_key: Optional path to library key CSV for variant mapping

    Returns:
        DataFrame with all summarized results
    """
    input_path = Path(input_dir)
    if not input_path.exists():
        raise FileNotFoundError(f"Input directory not found: {input_dir}")

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Collect (name, r1, r2) tuples from subdirectories
    candidates = []
    lib_folders = sorted([f for f in input_path.iterdir() if f.is_dir()])
    for lib_folder in lib_folders:
        r1, r2 = find_paired_files(lib_folder)
        if r1 and r2:
            candidates.append((lib_folder.name, r1, r2))

    # Flat directory fallback: no subdir had paired files
    if not candidates:
        r1, r2 = find_paired_files(input_path)
        if r1 and r2:
            candidates.append((input_path.name, r1, r2))

    n = len(candidates)
    print(f"Total libraries found: {n}")

    if n == 0:
        print("No paired FASTQ files found. Check your input directory structure.")
        return pd.DataFrame()

    all_summaries = []
    failures = []

    for i, (lib_name, r1_file, r2_file) in enumerate(candidates):
        print(f"\n[{i+1}/{n}] Processing: {lib_name}")

        try:
            # 1. Parse paired reads (fast Rust parser)
            df = parse_paired(str(r1_file), str(r2_file))
            df = df.dropna().copy()
            df['sample'] = lib_name

            print(f"   Paired reads: {len(df):,}")

            # 2. Summarize by motif
            summary = editing_counts(df, min_occurrence=min_occurrence)

            # 3. Save summary
            sample_out = output_path / f"{lib_name}_edits.csv"
            summary.to_csv(sample_out, index=False)

            all_summaries.append(summary)
            print(f"   Unique motifs: {len(summary):,}")
            print(f"   Saved: {sample_out}")

        except Exception as e:
            print(f"   Error in {lib_name}: {e}")
            import traceback
            traceback.print_exc()
            failures.append((lib_name, str(e)))

    # 4. Consolidate all results
    if all_summaries:
        final_df = pd.concat(all_summaries, ignore_index=True)
        final_path = output_path / "final_consolidated_edits.csv"
        final_df.to_csv(final_path, index=False)

        print(f"\n{'='*60}")
        print(f"DONE!")
        print(f"Samples processed: {len(all_summaries)}")
        print(f"Total motif-sample pairs: {len(final_df):,}")
        print(f"Saved: {final_path}")

        # 5. Apply library key if provided
        if library_key:
            print(f"\nApplying library key: {library_key}")
            key_df = apply_library_key(final_df, library_key)
            key_path = output_path / "final_library_key_counts.csv"
            key_df.to_csv(key_path, index=False)
            print(f"Saved: {key_path}")

        if failures:
            print(f"\nFailed samples ({len(failures)}):")
            for name, err in failures:
                print(f"  {name}: {err}")

        return final_df

    if failures:
        print(f"\nFailed samples ({len(failures)}):")
        for name, err in failures:
            print(f"  {name}: {err}")

    return pd.DataFrame()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Process all RBPscan samples.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("input_dir", help="Directory containing sample folders or flat paired FASTQ files")
    parser.add_argument("output_dir", help="Where to save results")
    parser.add_argument("--min-occurrence", type=int, default=10,
                        help="Minimum reads per motif to include in output")
    parser.add_argument("--library-key", metavar="PATH",
                        help="Path to library key CSV for variant mapping (columns: name, sequence)")
    args = parser.parse_args()

    process_all(
        args.input_dir,
        args.output_dir,
        min_occurrence=args.min_occurrence,
        library_key=args.library_key,
    )
