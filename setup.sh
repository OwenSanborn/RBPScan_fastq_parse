#!/usr/bin/env bash
# RBPScan FASTQ Parser — one-time setup
set -e

echo "=== RBPScan FASTQ Parser Setup ==="
echo ""

# ── 1. Rust ──────────────────────────────────────────────────────────────────
if ! command -v cargo &>/dev/null; then
    echo "Rust not found. Installing via rustup..."
    curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y --quiet
    # Activate for this session
    source "$HOME/.cargo/env"
    echo "Rust installed."
else
    echo "Rust found: $(cargo --version)"
fi

# ── 2. Build the Rust binary ──────────────────────────────────────────────────
echo ""
echo "Building fast parser (this takes ~1 minute the first time)..."
cd fastq_parser_rs
cargo build --release --quiet
cd ..
echo "Build complete."

# ── 3. Python dependencies ────────────────────────────────────────────────────
echo ""
echo "Installing Python dependencies..."
pip install -r requirements.txt --quiet
echo "Done."

echo ""
echo "=== Setup complete! ==="
echo ""
echo "Run the parser:"
echo "  python process_all_samples.py <input_dir> <output_dir>"
echo ""
echo "With a library key:"
echo "  python process_all_samples.py <input_dir> <output_dir> --library-key variants.csv"
