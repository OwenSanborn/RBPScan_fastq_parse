use clap::Parser;
use flate2::read::MultiGzDecoder;
use memchr::memmem;
use rayon::prelude::*;
use regex::Regex;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::time::Instant;

#[derive(Parser, Debug)]
#[command(name = "fastq_parser")]
#[command(about = "Ultra-fast FASTQ parser for RBPscan motif/hairpin extraction")]
struct Args {
    /// Input FASTQ file (gzipped or plain)
    #[arg(short, long)]
    input: String,

    /// Parse mode: 'motif' for R1, 'hairpin' for R2
    #[arg(short, long, default_value = "motif")]
    mode: String,

    /// Custom regex pattern (optional). For motif, must have 1 capture group.
    #[arg(short, long)]
    pattern: Option<String>,

    /// Number of threads (0 = auto)
    #[arg(short, long, default_value = "0")]
    threads: usize,

    /// Output file (default: stdout)
    #[arg(short, long)]
    output: Option<String>,
}

/// Read FASTQ records from a file (gzipped or plain)
fn read_fastq_records(path: &str) -> Vec<(String, String)> {
    let file = File::open(path).expect("Cannot open file");
    let mut reader: Box<dyn BufRead> = if path.ends_with(".gz") {
        Box::new(BufReader::with_capacity(8 * 1024 * 1024, MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::with_capacity(8 * 1024 * 1024, file))
    };

    let mut records = Vec::with_capacity(10_000_000);
    let mut header = String::with_capacity(256);
    let mut seq = String::with_capacity(256);
    let mut plus = String::with_capacity(8);
    let mut qual = String::with_capacity(256);

    loop {
        header.clear();
        seq.clear();
        plus.clear();
        qual.clear();

        if reader.read_line(&mut header).unwrap() == 0 {
            break;
        }
        reader.read_line(&mut seq).unwrap();
        reader.read_line(&mut plus).unwrap();
        reader.read_line(&mut qual).unwrap();

        let read_id = header
            .trim()
            .split_whitespace()
            .next()
            .unwrap_or(&header)
            .trim_start_matches('@')
            .to_string();

        records.push((read_id, seq.trim().to_string()));
    }
    records
}

// Pre-filter anchors for fast rejection (SIMD-accelerated with memchr)
const MOTIF_LEFT: &[u8] = b"TTCTGGCTGACATA";
const MOTIF_RIGHT: &[u8] = b"ATACAATCAGATATGCA";
const HAIRPIN_ANCHOR: &[u8] = b"AATTT";
const EDIT_MARKER: &[u8] = b"AATCC";

/// Parse motif from R1 reads with SIMD pre-filter optimization
fn parse_motif(records: &[(String, String)], pattern: &Regex) -> Vec<(String, Option<String>)> {
    // Pre-build SIMD finders (amortize setup cost)
    let left_finder = memmem::Finder::new(MOTIF_LEFT);
    let right_finder = memmem::Finder::new(MOTIF_RIGHT);

    records
        .par_iter()
        .map(|(read_id, seq)| {
            let seq_bytes = seq.as_bytes();
            // SIMD-accelerated substring search
            let motif = if left_finder.find(seq_bytes).is_some()
                && right_finder.find(seq_bytes).is_some()
            {
                pattern.captures(seq).map(|caps| caps[1].to_string())
            } else {
                None
            };
            (read_id.clone(), motif)
        })
        .collect()
}

/// Parse hairpin from R2 reads with SIMD pre-filter optimization
fn parse_hairpin(
    records: &[(String, String)],
    pattern: &Regex,
) -> Vec<(String, Option<String>, Option<u8>)> {
    let anchor_finder = memmem::Finder::new(HAIRPIN_ANCHOR);
    let edit_finder = memmem::Finder::new(EDIT_MARKER);

    records
        .par_iter()
        .map(|(read_id, seq)| {
            let seq_bytes = seq.as_bytes();
            // SIMD-accelerated anchor check
            if anchor_finder.find(seq_bytes).is_none() {
                return (read_id.clone(), None, None);
            }
            if let Some(caps) = pattern.captures(seq) {
                let hp = caps[1].to_string();
                // SIMD-accelerated edit counting
                let edits = edit_finder.find_iter(hp.as_bytes()).count() as u8;
                (read_id.clone(), Some(hp), Some(edits))
            } else {
                (read_id.clone(), None, None)
            }
        })
        .collect()
}

fn main() {
    let args = Args::parse();

    if args.threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global()
            .unwrap();
    }

    let default_motif = r"TTCTGGCTGACATA(.{11})ATACAATCAGATATGCA";
    // Hairpin pattern: capture group excludes AATTT anchor (no lookahead needed)
    let default_hairpin =
        r"(AA[TC]C[TC]AA[TC]C[TC]AA[TC]C[TC]AA[TC]C[TC]AA[TC]C[TC]AA[TC]C[TC])AATTT";

    let pattern_str = args.pattern.as_deref().unwrap_or(match args.mode.as_str() {
        "motif" => default_motif,
        "hairpin" => default_hairpin,
        _ => panic!("Mode must be 'motif' or 'hairpin'"),
    });

    eprintln!("Reading FASTQ: {}", args.input);
    let start = Instant::now();

    let records = read_fastq_records(&args.input);
    let read_time = start.elapsed();
    eprintln!(
        "Read {} records in {:.2}s",
        records.len(),
        read_time.as_secs_f64()
    );

    let parse_start = Instant::now();

    // Setup buffered output
    let mut output: Box<dyn Write> = match &args.output {
        Some(path) => Box::new(BufWriter::with_capacity(
            8 * 1024 * 1024,
            File::create(path).expect("Cannot create output file"),
        )),
        None => Box::new(BufWriter::with_capacity(8 * 1024 * 1024, std::io::stdout())),
    };

    match args.mode.as_str() {
        "motif" => {
            let pattern = Regex::new(pattern_str).expect("Invalid regex pattern");
            writeln!(output, "read_id\tmotif").unwrap();
            let results = parse_motif(&records, &pattern);
            let parse_time = parse_start.elapsed();

            let mut matches = 0;
            for (read_id, motif) in &results {
                if motif.is_some() {
                    matches += 1;
                }
                writeln!(output, "{}\t{}", read_id, motif.as_deref().unwrap_or("")).unwrap();
            }

            let total_time = start.elapsed();
            eprintln!(
                "Parsed in {:.2}s | Total: {:.2}s | {:.0} reads/sec | Matches: {} ({:.1}%)",
                parse_time.as_secs_f64(),
                total_time.as_secs_f64(),
                records.len() as f64 / total_time.as_secs_f64(),
                matches,
                100.0 * matches as f64 / records.len() as f64
            );
        }
        "hairpin" => {
            let pattern = Regex::new(pattern_str).expect("Invalid regex pattern");
            writeln!(output, "read_id\thp\tedits_count").unwrap();
            let results = parse_hairpin(&records, &pattern);
            let parse_time = parse_start.elapsed();

            let mut matches = 0;
            for (read_id, hp, edits) in &results {
                if hp.is_some() {
                    matches += 1;
                }
                writeln!(
                    output,
                    "{}\t{}\t{}",
                    read_id,
                    hp.as_deref().unwrap_or(""),
                    edits.map(|e| e.to_string()).unwrap_or_default()
                )
                .unwrap();
            }

            let total_time = start.elapsed();
            eprintln!(
                "Parsed in {:.2}s | Total: {:.2}s | {:.0} reads/sec | Matches: {} ({:.1}%)",
                parse_time.as_secs_f64(),
                total_time.as_secs_f64(),
                records.len() as f64 / total_time.as_secs_f64(),
                matches,
                100.0 * matches as f64 / records.len() as f64
            );
        }
        _ => panic!("Mode must be 'motif' or 'hairpin'"),
    }
}
