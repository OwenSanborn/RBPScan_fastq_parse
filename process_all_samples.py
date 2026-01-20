"""
Process all RBPscan samples.
1. Parse R1 (motif) and R2 (hairpin) with fast Rust parser
2. Merge paired reads
3. Summarize by motif with editing_counts()
4. Output: sample, motif, occurrence, total_edits, reads_edited, all_reads, score, efficiency_1
"""
import pandas as pd
import numpy as np
from pathlib import Path
from fastq_parser import parse_paired


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


def process_all(input_dir, output_dir):
    """
    Process all sample directories.

    Args:
        input_dir: Directory containing sample folders (each with *_1.fq.gz and *_2.fq.gz)
        output_dir: Where to save results

    Returns:
        DataFrame with all summarized results
    """
    input_path = Path(input_dir)
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    all_summaries = []
    lib_folders = sorted([f for f in input_path.iterdir() if f.is_dir()])

    print(f"Total libraries found: {len(lib_folders)}")

    for lib_folder in lib_folders:
        r1_files = list(lib_folder.glob("*_1.fq.gz"))
        r2_files = list(lib_folder.glob("*_2.fq.gz"))

        if r1_files and r2_files:
            lib_name = lib_folder.name
            print(f"\nProcessing: {lib_name}")

            try:
                # 1. Parse paired reads (fast Rust parser)
                df = parse_paired(str(r1_files[0]), str(r2_files[0]))
                df = df.dropna().copy()
                df['sample'] = lib_name

                print(f"   Paired reads: {len(df):,}")

                # 2. Summarize by motif
                summary = editing_counts(df, min_occurrence=10)

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

        return final_df

    return pd.DataFrame()


if __name__ == "__main__":
    import sys

    if len(sys.argv) < 3:
        print("Usage: python process_all_samples.py <input_dir> <output_dir>")
        print("\nExample:")
        print("  python process_all_samples.py /path/to/01.RawData/RBPscan /path/to/out")
        sys.exit(1)

    INPUT_DIR = sys.argv[1]
    OUTPUT_DIR = sys.argv[2]

    process_all(INPUT_DIR, OUTPUT_DIR)
