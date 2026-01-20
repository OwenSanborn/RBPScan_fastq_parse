import os
import re
import pandas as pd
import numpy as np
import pyfastx

# Parse FASTQ file to extract motifs (FWD) or hairpins (REV)
def parse_fastq(fq_path, parse_motif=None):
    """
    Parses both the 11nt motif and the RC Hairpin from a single forward read.
    """
    records = []
    fq = pyfastx.Fastx(fq_path)
    
    print(f"📥 Single-pass parsing: {os.path.basename(fq_path)}")

    for name, seq, qual in fq:
        read_id = name.split(' ')[0]
        
        # 1. Search for the 11nt Motif (Forward)
        motif_match = re.search(r'TTCTGGCTGACATA(.{11})ATACAATCAGATATGCA', seq)
        motif_val = motif_match.group(1) if motif_match else None
        
        # 2. Search for the Hairpin RC (Forward)
        # Anchor AAATT followed by 6 repeats of [AG]G[AG]TT
        hp_match = re.search(r"AAATT(([AG]G[AG]TT){6})", seq)
        
        hp_region = None
        edits_count = None
        
        if hp_match:
            hp_region = hp_match.group(1)
            # Counting 'GGATT' as the forward-strand equivalent of 'AATCC'
            edits_count = hp_region.count("GGATT")

        records.append({
            'read_id': read_id,
            'motif': motif_val,
            'hp_rc': hp_region,
            'edits_count': edits_count
        })

    df = pd.DataFrame(records)
    
    # Optional: Filter for reads that actually matched something to save space
    # df = df.dropna(subset=['motif', 'hp_rc'], how='all')
    
    print(f"✅ Parsed {len(df)} reads.")
    return df
        
#Merge FWD and REV reads by read_id
def merge_paired_reads(fwd_df, rev_df):
    fwd_ids = set(fwd_df['read_id'])
    rev_ids = set(rev_df['read_id'])
    shared_ids = fwd_ids.intersection(rev_ids)

    print(f"🔍 FWD IDs: {len(fwd_ids)}")
    print(f"🔍 REV IDs: {len(rev_ids)}")
    print(f"🔗 Shared IDs: {len(shared_ids)}")

    if len(shared_ids) == 0:
        # Print some examples
        print("🧪 FWD example:", next(iter(fwd_ids)))
        print("🧪 REV example:", next(iter(rev_ids)))
        # Print 5 near matches (string similarity)
        import difflib
        example = next(iter(fwd_ids))
        print("🔍 Close matches in REV:")
        print(difflib.get_close_matches(example, rev_ids, n=5, cutoff=0.8))

    merged = pd.merge(fwd_df, rev_df, on='read_id', how='inner')
    print(f"✅ Merged {len(merged)} reads.")
    merged = merged.dropna(subset=['motif', 'hp', 'edits_count'])
    print(f"🧼 After dropping nulls: {len(merged)} reads remain.")
    return merged

# Compute editing scores and efficiencies per motif per sample
def editing_counts(reads: pd.DataFrame, min_occurrence: int = 10, sites: int = 6, verbose: bool = False) -> pd.DataFrame:
    """
    Vectorized version: computes motif editing metrics per sample with a single groupby.
    - min_occurrence: require at least this many reads for a (sample, motif) to be kept.
    - sites: denominator used in efficiency_1 (your code used 6; keep it configurable).
    """
    cols = ['sample','motif','edits_count']
    if reads is None or reads.empty or not set(cols).issubset(reads.columns):
        if verbose:
            print("⚠️ No reads or required columns missing.")
        return pd.DataFrame(columns=['sample','motif','occurrence','reads_edited','total_edits','all_reads','score','efficiency_1'])

    # Keep only what we need and ensure dtypes are fast
    df = reads.loc[:, cols].dropna()
    # Make sure edits_count is numeric (and small ints for speed/memory)
    df['edits_count'] = pd.to_numeric(df['edits_count'], errors='coerce').fillna(0).astype(np.int16)
    # Categorical speeds up groupby a lot when many repeats
    df['sample'] = df['sample'].astype('category')
    df['motif']  = df['motif'].astype('category')

    # Per-sample total reads (for all_reads)
    sample_sizes = df.groupby('sample', observed=True).size()

    # Group once and aggregate everything we need
    g = df.groupby(['sample','motif'], observed=True)['edits_count']
    out = g.agg(
        occurrence='size',
        total_edits='sum',
        reads_edited=lambda s: np.count_nonzero(s.values > 0)
    ).reset_index()

    # Filter low-support motifs
    out = out[out['occurrence'] > min_occurrence]

    # Add per-sample totals and derived metrics
    out['all_reads'] = out['sample'].map(sample_sizes)
    out['score'] = out['total_edits'] / out['occurrence']
    out['efficiency_1'] = out['total_edits'] / (sites * out['occurrence'] + 1)

    # Optional: tidy ordering
    out = out.sort_values(['sample','occurrence','total_edits'], ascending=[True, False, False]).reset_index(drop=True)

    if verbose:
        n_motifs = out.groupby('sample')['motif'].nunique()
        for s, k in n_motifs.items():
            print(f"✔️ Motifs kept in sample {s}: {k} (>= {min_occurrence} reads)")
        print(f"\n📈 Editing scores calculated for {len(out)} motif-sample pairs.")

    return out


import os
import re
import pandas as pd
from current_rbpscan_parse import parse_fastq, editing_counts
from pathlib import Path

def process_fwd_samples(fq_folder, base):
    """
    Process forward reads, annotate with metadata, and handle 'NA' miRNA.
    """
    
    # 1. SETUP PATHS
    reads_out_dir = os.path.join(base, "processed_data/RBPscan/")
    edits_out_dir = os.path.join(base, "out/")
    os.makedirs(reads_out_dir, exist_ok=True)
    os.makedirs(edits_out_dir, exist_ok=True)
    
    target_directory = Path(fq_folder)
    all_edits_list = [] # To store individual dataframes for final merge
    
    # 2. ITERATE THROUGH SAMPLES
    for sample_dir in sorted(target_directory.iterdir()):
        if not sample_dir.is_dir():
            continue
            
        sample_name = sample_dir.name
        fwd_files = list(sample_dir.glob("*_1.fq.gz"))
        
        if not fwd_files:
            print(f"⚠️ Skipping {sample_name}: No *_1.fq.gz found.")
            continue
        
        fwd_path = str(fwd_files[0])
        
        try:
            print(f"Processing: {sample_name}...")
            
            # 3. PARSE FORWARD READS
            fwd_df = parse_fastq(fwd_path, parse_motif=True)
            if fwd_df.empty:
                print(f"  ⚠️ {sample_name} is empty.")
                continue
            
            # 4. ANNOTATE METADATA (Adapted from your R logic)
            fwd_df['sample'] = sample_name
            
            # 5. SAVE PROCESSED READS
            reads_csv = os.path.join(reads_out_dir, f"{sample_name}_reads.csv")
            fwd_df.to_csv(reads_csv, index=False)
            
            # 6. CALCULATE EDITS
            edits = editing_counts(fwd_df)
            
            # Attach metadata to the summary table
            meta_cols = ['sample', 'library', 'protein', 'mirna', 'rep']
            for col in meta_cols:
                edits[col] = fwd_df[col].iloc[0]

            # Save individual edit file
            edits_csv = os.path.join(edits_out_dir, f"{sample_name}_edits.csv")
            edits.to_csv(edits_csv, index=False)
            
            all_edits_list.append(edits)
            print(f"  ✅ Done.")
            
        except Exception as e:
            print(f"  ❌ Error in {sample_name}: {e}")
            continue

    # 7. CONSOLIDATE ALL RESULTS
    if all_edits_list:
        final_df = pd.concat(all_edits_list, ignore_index=True)
        final_csv = os.path.join(base, "final_consolidated_edits.csv")
        final_df.to_csv(final_csv, index=False)
        print(f"\n✨ SUCCESS: Consolidated file saved to {final_csv}")

