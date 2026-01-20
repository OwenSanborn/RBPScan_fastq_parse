# RBPScan FASTQ Parser

Fast FASTQ parser for RBPscan motif and hairpin extraction. Built with Rust for speed, with a simple Python interface.

## Performance

| Dataset | Python | Rust | Speedup |
|---------|--------|------|---------|
| 9.2M reads (gzipped) | 11-13s | 7-8s | **1.5x** |

## Installation

### 1. Clone the repo

```bash
git clone https://github.com/YOUR_USERNAME/RBPScan_fastq_parse.git
cd RBPScan_fastq_parse
```

### 2. Build the Rust binary

Requires [Rust](https://rustup.rs/) to be installed.

```bash
cd fastq_parser_rs
cargo build --release
cd ..
```

### 3. (Optional) Install Python dependencies

```bash
pip install pandas
```

## Usage

### Python (Recommended)

```python
from fastq_parser import parse_motif, parse_hairpin, parse_paired

# Parse R1 motifs - returns pandas DataFrame
motif_df = parse_motif("sample_R1.fq.gz")

# Parse R2 hairpins
hairpin_df = parse_hairpin("sample_R2.fq.gz")

# Parse both and merge by read_id
merged_df = parse_paired("sample_R1.fq.gz", "sample_R2.fq.gz")

# Custom regex pattern
df = parse_motif("file.fq.gz", pattern=r"CUSTOM(.{11})PATTERN")
```

### Command Line

```bash
# Motif extraction (R1)
./fastq_parser_rs/target/release/fastq_parser \
    -i input_R1.fq.gz \
    -m motif \
    -o motifs.tsv

# Hairpin extraction (R2)
./fastq_parser_rs/target/release/fastq_parser \
    -i input_R2.fq.gz \
    -m hairpin \
    -o hairpins.tsv
```

#### Options

| Flag | Description |
|------|-------------|
| `-i, --input` | Input FASTQ file (gzipped or plain) |
| `-m, --mode` | `motif` (R1) or `hairpin` (R2) |
| `-o, --output` | Output TSV file (stdout if omitted) |
| `-p, --pattern` | Custom regex pattern (optional) |
| `-t, --threads` | Number of threads (0 = auto) |

## Output Format

### Motif (R1)

| Column | Description |
|--------|-------------|
| `read_id` | FASTQ read identifier |
| `motif` | 11bp motif sequence (empty if no match) |

### Hairpin (R2)

| Column | Description |
|--------|-------------|
| `read_id` | FASTQ read identifier |
| `hp` | 30bp hairpin region |
| `edits_count` | Number of AATCC (edited) sites (0-6) |

## Default Patterns

- **Motif (R1)**: `TTCTGGCTGACATA(.{11})ATACAATCAGATATGCA`
- **Hairpin (R2)**: `(AA[TC]C[TC]){6}AATTT` - captures 6 repeats before AATTT anchor

## Files

```
fastq_parser.py          # Python wrapper (import this)
fastq_parser_rs/         # Rust source code
  ├── Cargo.toml         # Rust dependencies
  └── src/main.rs        # Parser implementation
baseline_parser.py       # Pure Python baseline (slower)
benchmark.py             # Performance comparison script
```
