# RBPScan FASTQ Parser

Fast parser for RBPscan FASTQ data. Extracts motifs (R1) and hairpin edit counts (R2), merges by read, and summarizes editing rates per motif. Optionally maps results to a known variant library.

---

## Quick start

### 1. Get the code

```bash
git clone https://github.com/OwenSanborn/RBPScan_fastq_parse.git
cd RBPScan_fastq_parse
```

### 2. Run setup (one time only)

```bash
bash setup.sh
```

This installs Rust (if needed), builds the fast parser, and installs Python dependencies. Takes about 1 minute the first time.

### 3. Run on your data

```bash
python process_all_samples.py /path/to/raw_data/ /path/to/output/
```

Your `raw_data/` folder can contain either:
- **Subdirectories**, one per sample (e.g. `raw_data/Sample1/`, `raw_data/Sample2/`) — each containing paired FASTQ files
- **A flat folder** with a single pair of FASTQ files

Supported file naming conventions: `_R1_001.fastq.gz`, `_R1.fq.gz`, `_1.fq.gz`, and others (auto-detected).

---

## Output

Results are saved to your output folder:

| File | Contents |
|------|----------|
| `<sample>_edits.csv` | Per-sample motif summary |
| `final_consolidated_edits.csv` | All samples combined |
| `final_library_key_counts.csv` | Per-variant counts (only if `--library-key` is used) |

Each row in the summary:

| Column | Description |
|--------|-------------|
| `sample` | Sample name |
| `motif` | 11bp motif sequence |
| `occurrence` | Number of reads containing this motif |
| `total_edits` | Total edited bases across all reads |
| `reads_edited` | Reads with at least one edit |
| `all_reads` | Total paired reads for this sample |
| `score` | Mean edits per read (`total_edits / occurrence`) |
| `efficiency_1` | Editing efficiency (`total_edits / (6 × occurrence + 1)`) |

---

## Library key (variant mapping)

If you have a defined set of known variants, you can map observed motifs to them. Create a CSV with `name` and `sequence` columns:

```
name,sequence
Sox2_motif,AATCAATGG
Oct4_motif,TTTGCATA
Klf4_motif,CCGCCCGC
```

Then run:

```bash
python process_all_samples.py raw_data/ output/ --library-key variants.csv
```

This produces `final_library_key_counts.csv` with one row per `(sample, variant)`. Variants with no observed reads appear as 0-count rows.

---

## All options

```
python process_all_samples.py <input_dir> <output_dir> [options]

Options:
  --min-occurrence N    Minimum reads per motif to report (default: 10)
  --library-key PATH    Path to variant library CSV
```

---

## Requirements

- Python 3.8+
- pandas, numpy (installed by `setup.sh`)
- Rust (installed automatically by `setup.sh`)

---

## Files

```
setup.sh                 # One-time setup script
process_all_samples.py   # Main script — run this
fastq_parser.py          # Python wrapper for the Rust parser
fastq_parser_rs/         # Rust source (built by setup.sh)
requirements.txt         # Python dependencies
```
