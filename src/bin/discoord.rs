use clap::{CommandFactory, Parser};
use dfam_coord::{
    derive_assembly_name, detect_format_and_compression, find_reference_file, init_thread_pool,
    load_reference, output_results, parse_delimited_file, parse_fasta, parse_stockholm,
    process_sequences, write_delimited_output, write_fasta_output, write_stockholm_output,
    LogLevel, SequenceRecord,
};
use std::collections::HashMap;
use std::env;
use std::fs;
use std::path::{Path, PathBuf};
use std::process;
use std::time::Instant;

#[derive(Parser, Debug)]
#[command(name = "discoord", author, version, about, long_about = None)]
struct Args {
    /// Directory containing reference files
    #[arg(short = 'd', long)]
    reference_dir: Option<String>,

    /// Default reference file for sequences without an assembly_id
    #[arg(short = 'r', long)]
    reference_default: Option<String>,

    /// Enable mapping (Aho-Corasick) for invalid identifiers
    #[arg(short, long, default_value = "false")]
    map_sequences: bool,

    /// Use Boyer-Moore instead of Aho-Corasick for the mapping step.
    /// Only meaningful with -m.
    #[arg(long, default_value = "false")]
    boyer_moore: bool,

    /// Log level (Summary, PerRecord or Detailed)
    #[arg(short, long, default_value = "summary")]
    log_level: LogLevel,

    /// Input file paths (Fasta, Stockholm or Tab/Comma Delimited files)
    #[arg()]
    input: Vec<String>,

    /// Threads to use for parallel processing
    #[arg(short = 'x', long)]
    threads: Option<usize>,

    /// Generate a new output file per input file in this directory (or <input_file>.discoord if ".")
    #[arg(short = 'o', long)]
    output_dir: Option<String>,

    /// When mapping, remove a sequence if every genomic hit for it is already
    /// occupied by an earlier sequence in the file.  Sequences with at least one
    /// unoccupied hit are kept and assigned to that position as usual.
    /// Removed sequences are reported as "removed_remapped_duplicate" in the summary.
    #[arg(short = 'u', long, default_value = "false")]
    remove_duplicates: bool,
}

fn main() {
    let args = Args::parse();
    let t_start = Instant::now();
    let debug_mode = false;
    let mut validation_failed = false;

    if args.input.is_empty() {
        eprintln!("Error: No input files provided.");
        eprintln!("{}", Args::command().render_long_help());
        process::exit(1);
    }

    println!("##\n## DisCoord Version {}\n##", env!("CARGO_PKG_VERSION"));

    if let Some(n) = args.threads {
        init_thread_pool(n);
        println!("## Threads: {}", n);
    } else {
        println!("## Threads: all available");
    }

    if args.map_sequences {
        let method = if args.boyer_moore { "boyer-moore" } else { "aho-corasick" };
        println!("## Mapping method: {}", method);
    }

    let mut sequences_by_assembly: HashMap<Option<String>, Vec<SequenceRecord>> = HashMap::new();
    let mut all_metadata: Vec<String> = Vec::new();

    if args.output_dir.is_some() {
        let dir_path = Path::new(args.output_dir.as_ref().unwrap());
        if !dir_path.is_dir() {
            fs::create_dir_all(dir_path).expect("Failed to create output directory");
        }
    }

    for input_file in &args.input {
        let (is_gzip, format) = detect_format_and_compression(input_file)
            .expect("Failed to detect input format and compression");
        println!(
            "## File: {}, Format: {}, Compression: {}",
            input_file,
            format,
            if is_gzip { "Gzip" } else { "None" }
        );

        let (mut sequences, metadata) = match format {
            "Fasta" => parse_fasta(input_file, is_gzip),
            "Stockholm" => match parse_stockholm(input_file, is_gzip) {
                Ok(result) => result,
                Err(e) => {
                    eprintln!("Error parsing Stockholm file {}: {}", input_file, e);
                    process::exit(1);
                }
            },
            "TabDelimited" | "CommaDelimited" => parse_delimited_file(input_file),
            _ => {
                eprintln!("Unsupported format for file: {}", input_file);
                process::exit(1);
            }
        };
        println!("##       Sequences: {}", sequences.len());

        let metadata_start_idx = all_metadata.len();
        all_metadata.extend(metadata.into_iter());
        for sequence in &mut sequences {
            sequence.metadata_idx += metadata_start_idx;
        }

        for record in sequences {
            sequences_by_assembly
                .entry(record.assembly_id.clone())
                .or_default()
                .push(record);
        }
    }
    println!("## Total references assemblies: {}", sequences_by_assembly.len());

    let mut processed_sequences: Vec<SequenceRecord> = Vec::new();

    if let Some(ref_dir) = &args.reference_dir {
        println!("## Reference directory: {}", ref_dir);
    }
    if let Some(ref_def) = &args.reference_default {
        println!("## Reference default: {}", ref_def);
    }

    // When using -r with -m, derive the assembly name from the reference filename
    // so successfully remapped records get their assembly_id updated.
    let remapped_assembly: Option<String> = if args.map_sequences {
        args.reference_default.as_ref().map(|path| derive_assembly_name(path))
    } else {
        None
    };
    if let Some(ref name) = remapped_assembly {
        println!("## Remapped assembly name: {}", name);
    }
    println!("##");

    if let Some(ref_dir) = &args.reference_dir {
        // Per-assembly lookup: each group gets its own reference file.
        for (assembly_id, sequences) in sequences_by_assembly {
            let ref_file = find_reference_file(ref_dir, &assembly_id, &args.reference_default);
            print!("## o Loading reference: {} ... ", ref_file);
            let t = Instant::now();
            let genome_map = load_reference(&ref_file).unwrap_or_else(|e| {
                eprintln!("## Error: Could not open reference file '{}': {}", ref_file, e);
                std::process::exit(1);
            });
            println!("{} sequences loaded in {:.1}s", genome_map.len(), t.elapsed().as_secs_f32());
            let mut results = process_sequences(sequences, &genome_map, args.map_sequences, !args.boyer_moore, debug_mode, remapped_assembly.as_deref(), args.remove_duplicates);
            for record in results.iter_mut() {
                if record.validated.is_none() {
                    record.validated = Some("invalid".to_string());
                    validation_failed = true;
                }
            }
            processed_sequences.extend(results);
        }
    } else if let Some(ref_default) = &args.reference_default {
        // Single reference: flatten all groups into one batch so there is only
        // one validation+mapping pass (and one progress bar) regardless of how
        // many distinct assembly_id prefixes appear in the input identifiers.
        print!("## o Loading reference: {} ... ", ref_default);
        let t = Instant::now();
        let genome_map = load_reference(ref_default).unwrap_or_else(|e| {
            eprintln!("## Error: Could not open reference file '{}': {}", ref_default, e);
            std::process::exit(1);
        });
        println!("{} sequences loaded in {:.1}s", genome_map.len(), t.elapsed().as_secs_f32());
        let all_sequences: Vec<SequenceRecord> = sequences_by_assembly.into_values().flatten().collect();
        let mut results = process_sequences(all_sequences, &genome_map, args.map_sequences, !args.boyer_moore, debug_mode, remapped_assembly.as_deref(), args.remove_duplicates);
        for record in results.iter_mut() {
            if record.validated.is_none() {
                record.validated = Some("invalid".to_string());
                validation_failed = true;
            }
        }
        processed_sequences.extend(results);
    } else {
        panic!("No reference file provided for sequences without an assembly_id");
    }

    processed_sequences.sort_by_key(|record| (record.input_file.clone(), record.order));

    let elapsed = t_start.elapsed();
    let total_secs = elapsed.as_secs();
    println!("## Total runtime: {:02}:{:02}:{:02}", total_secs / 3600, (total_secs % 3600) / 60, total_secs % 60);

    for input_file in &args.input {
        let input_path = Path::new(input_file);
        let input_dir = input_path.parent().unwrap_or_else(|| Path::new("."));
        let input_dir_canonical = if input_dir.as_os_str().is_empty() {
            env::current_dir().unwrap_or_else(|_| PathBuf::from("."))
        } else {
            fs::canonicalize(input_dir).unwrap_or_else(|_| input_dir.to_path_buf())
        };

        let output_file = match args.output_dir {
            Some(ref dir) => {
                let output_dir_path = Path::new(&dir);
                let output_dir_canonical = fs::canonicalize(output_dir_path)
                    .unwrap_or_else(|_| output_dir_path.to_path_buf());

                if input_dir_canonical == output_dir_canonical {
                    input_path
                        .with_extension("discoord")
                        .to_string_lossy()
                        .to_string()
                } else {
                    output_dir_path
                        .join(input_path.file_name().unwrap_or_default())
                        .to_string_lossy()
                        .to_string()
                }
            }
            None => String::from(""),
        };
        let (is_gzip, format) = detect_format_and_compression(input_file)
            .expect("Failed to detect input format and compression");

        let records: Vec<SequenceRecord> = processed_sequences
            .iter()
            .filter(|record| record.input_file == *input_file)
            .cloned()
            .collect();

        output_results(&records, args.log_level.clone(), format!("Summary for {}", input_file));
        if args.output_dir.is_some() && !output_file.is_empty() {
            let output_records: Vec<SequenceRecord> = records
                .iter()
                .filter(|r| r.validated.as_deref() != Some("removed_remapped_duplicate"))
                .cloned()
                .collect();
            match format {
                "Fasta" => write_fasta_output(&output_records, &all_metadata, &output_file, is_gzip, false)
                    .expect("Failed to write output"),
                "Stockholm" => write_stockholm_output(&output_records, &all_metadata, &output_file, is_gzip, false)
                    .expect("Failed to write output"),
                "TabDelimited" | "CommaDelimited" => {
                    write_delimited_output(&output_records, &output_file, is_gzip, false, format)
                        .expect("Failed to write output")
                }
                _ => {
                    eprintln!("Unsupported format for writing output: {}", format);
                    process::exit(1);
                }
            }
        }
    }

    if validation_failed {
        process::exit(1);
    } else {
        process::exit(0);
    }
}
