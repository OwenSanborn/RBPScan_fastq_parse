"""
Python wrapper for the fast Rust FASTQ parser.
Just call these functions - no Rust knowledge needed!
"""
import subprocess
import pandas as pd
from pathlib import Path

# Path to the Rust binary (relative to this file)
RUST_BIN = Path(__file__).parent / "fastq_parser_rs/target/release/fastq_parser"

_COLUMNS = {
    "motif": ["read_id", "motif"],
    "hairpin": ["read_id", "hp", "edits_count"],
}


def parse_motif(fq_path, output_path=None, pattern=None):
    """
    Parse R1 reads to extract 11nt motifs.

    Args:
        fq_path: Path to R1 FASTQ file (gzipped or plain)
        output_path: Optional output TSV path (returns DataFrame if None)
        pattern: Optional custom regex pattern

    Returns:
        DataFrame with columns: read_id, motif

    Example:
        df = parse_motif("sample_R1.fq.gz")
        print(df.head())
    """
    return _run_parser(fq_path, "motif", output_path, pattern)


def parse_hairpin(fq_path, output_path=None, pattern=None):
    """
    Parse R2 reads to extract hairpin regions and edit counts.

    Args:
        fq_path: Path to R2 FASTQ file (gzipped or plain)
        output_path: Optional output TSV path (returns DataFrame if None)
        pattern: Optional custom regex pattern

    Returns:
        DataFrame with columns: read_id, hp, edits_count

    Example:
        df = parse_hairpin("sample_R2.fq.gz")
        print(f"Average edits: {df['edits_count'].mean():.2f}")
    """
    return _run_parser(fq_path, "hairpin", output_path, pattern)


def parse_paired(r1_path, r2_path):
    """
    Parse paired R1/R2 files and merge by read_id.

    Returns:
        DataFrame with columns: read_id, motif, hp, edits_count

    Example:
        df = parse_paired("sample_R1.fq.gz", "sample_R2.fq.gz")
        # Filter to reads with valid motif AND hairpin
        valid = df.dropna()
    """
    motif_df = parse_motif(r1_path)
    hairpin_df = parse_hairpin(r2_path)
    merged = pd.merge(motif_df, hairpin_df, on="read_id", how="inner")
    return merged


def _run_parser(fq_path, mode, output_path=None, pattern=None):
    """Internal: run the Rust parser and return results."""
    import tempfile

    if not RUST_BIN.exists():
        raise FileNotFoundError(
            f"Rust binary not found at {RUST_BIN}\n"
            "Build it with: cd fastq_parser_rs && cargo build --release"
        )

    # Use temp file if no output specified
    if output_path is None:
        tmp = tempfile.NamedTemporaryFile(suffix='.tsv', delete=False)
        output_path = tmp.name
        tmp.close()
        return_df = True
    else:
        return_df = False

    fq_path = Path(fq_path)
    if not fq_path.exists():
        raise FileNotFoundError(f"Input file not found: {fq_path}")

    # Build command
    cmd = [str(RUST_BIN), "-i", str(fq_path), "-m", mode, "-o", str(output_path)]
    if pattern:
        cmd.extend(["-p", pattern])

    # Run parser
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"Parser failed: {result.stderr}")

    # Print timing info
    print(result.stderr.strip())

    if Path(output_path).stat().st_size == 0:
        if return_df:
            Path(output_path).unlink()
        return pd.DataFrame(columns=_COLUMNS[mode])

    if return_df:
        df = pd.read_csv(output_path, sep='\t', na_values=[""], keep_default_na=True)
        Path(output_path).unlink()  # Clean up temp file
        return df
    else:
        return pd.read_csv(output_path, sep='\t', na_values=[""], keep_default_na=True)


# Quick test
if __name__ == "__main__":
    print("=== Testing Python Wrapper ===\n")

    r1 = "raw_data/Lib2040_1866_2700_rep1_CKDL250033540-1A_2357TCLT4_L6_1.fq.gz"
    r2 = "raw_data/Lib2040_1866_2700_rep1_CKDL250033540-1A_2357TCLT4_L6_2.fq.gz"

    print("Parsing R1 (motif)...")
    motif_df = parse_motif(r1)
    print(f"  Shape: {motif_df.shape}")
    print(f"  Motifs found: {motif_df['motif'].notna().sum():,}")
    print(motif_df.head())

    print("\nParsing R2 (hairpin)...")
    hp_df = parse_hairpin(r2)
    print(f"  Shape: {hp_df.shape}")
    print(f"  Hairpins found: {hp_df['hp'].notna().sum():,}")
    print(hp_df.head())

    print("\nEdit count distribution:")
    print(hp_df['edits_count'].value_counts().sort_index())
