use clap::ValueEnum;
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use rayon::prelude::*;
use aho_corasick::AhoCorasick;
use indicatif::{ProgressBar, ProgressStyle};
use std::time::Duration;
use dfam_stk_io::{IDVersion, SeqRow};
use std::collections::HashMap;
use std::fs::{File, OpenOptions};
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};
use std::path::Path;
use bio::io::fasta;
use memmap2::Mmap;

#[derive(Clone, Debug, ValueEnum, PartialEq)]
pub enum LogLevel {
    Summary,
    PerRecord,
    Detailed,
}

/// Controls for the interval pass that runs after validation and mapping.
///
/// discoord always detects and reports; these two switches only decide whether
/// it acts.  They are independent: either, both or neither may be set.
#[derive(Clone, Debug)]
pub struct IntervalOptions {
    /// Drop records whose span sits inside another record's span.
    pub remove_contained: bool,
    /// Fold each run of mutually overlapping records into one record whose
    /// sequence is re-extracted from the reference over the union span.
    pub merge_overlapping: bool,
    /// Overlap needed before two records join the same merge cluster.  A value
    /// of 1 or more is a base-pair count; a value between 0 and 1 is a fraction
    /// of the shorter of the two spans.
    pub min_overlap: f64,
}

impl Default for IntervalOptions {
    fn default() -> Self {
        IntervalOptions { remove_contained: false, merge_overlapping: false, min_overlap: 1.0 }
    }
}

#[derive(Clone, Debug)]
pub struct SequenceRecord {
    pub input_file: String,
    pub metadata_idx: usize,
    pub order: usize,
    pub original_id: Option<String>,
    pub assembly_id: Option<String>,
    pub sequence_id: String,
    /// 1-based fully closed, as written in the identifier. This crate audits
    /// identifiers in their published form, so it keeps that convention and
    /// says so in the name.
    pub start_1b: Option<u64>,
    pub end_1b: Option<u64>,
    pub orient: Option<char>,
    pub inferred_version: Option<IDVersion>,
    pub sequence: Vec<u8>,
    pub aligned_seq: Option<Vec<u8>>,
    pub validated: Option<String>,
}

impl SequenceRecord {
    /// Build a `SequenceRecord` from an already-parsed `SeqRow`.
    ///
    /// Gap characters (`.` and `-`) are stripped from `aligned_seq` to produce
    /// the ungapped `sequence` used for coordinate validation.  When
    /// `aligned_seq` is empty (e.g. during FASTA parsing before sequence lines
    /// are read), both `sequence` and `aligned_seq` are left empty/None and
    /// must be filled in by the caller.
    pub fn from_seq_row(row: &SeqRow, file_path: &str, order: usize, metadata_idx: usize) -> Self {
        let (sequence, aligned_seq) = if row.aligned_seq.is_empty() {
            (Vec::new(), None)
        } else {
            let seq: Vec<u8> = row.aligned_seq.bytes()
                .filter(|&b| b != b'.' && b != b'-')
                .collect();
            (seq, Some(row.aligned_seq.as_bytes().to_vec()))
        };
        SequenceRecord {
            input_file: file_path.to_string(),
            metadata_idx,
            order,
            original_id: Some(row.original_id.clone()),
            assembly_id: row.assembly_id.clone(),
            sequence_id: row.sequence_id.clone().unwrap_or_else(|| row.original_id.clone()),
            start_1b: row.span.and_then(|sp| sp.as_1b_closed()).map(|(s, _)| s),
            end_1b: row.span.and_then(|sp| sp.as_1b_closed()).map(|(_, e)| e),
            orient: row.orient,
            inferred_version: row.inferred_version.clone(),
            sequence,
            aligned_seq,
            validated: None,
        }
    }

    /// Format the sequence identifier for output, prepending the assembly ID
    /// when present: `assembly:sequence_id` or just `sequence_id`.
    pub fn format_id(&self) -> String {
        match &self.assembly_id {
            Some(assembly) => format!("{}:{}", assembly, self.sequence_id),
            None => self.sequence_id.clone(),
        }
    }

    pub fn print_record(&self) {
        println!(
            "Smitten::Identifier: original_id: {}, assembly_id: {}, sequence_id: {}, start: {}, end: {}, orient: {}, inferred_version: {:?}, validated: {}",
            self.original_id.as_deref().unwrap_or("Unknown"),
            self.assembly_id.as_deref().unwrap_or("None"),
            self.sequence_id,
            self.start_1b.unwrap_or(0),
            self.end_1b.unwrap_or(0),
            self.orient.unwrap_or('?'),
            self.inferred_version,
            self.validated.as_deref().unwrap_or(""),
        );
    }
}

pub fn find_reference_file(ref_dir: &str, assembly_id: &Option<String>, default_reference: &Option<String>) -> String {
    if let Some(id) = assembly_id {
        let two_bit_path = Path::new(ref_dir).join(format!("{}.2bit", id));
        if two_bit_path.exists() {
            return two_bit_path.to_string_lossy().to_string();
        }
        let fa_path = Path::new(ref_dir).join(format!("{}.fa", id));
        if fa_path.exists() {
            return fa_path.to_string_lossy().to_string();
        }
        panic!("No reference file found for assembly_id: {:?}", assembly_id);
    } else {
        if let Some(default_ref) = default_reference {
            return default_ref.clone();
        } else {
            panic!("No assembly_id provided and no default reference file specified.");
        }
    }
}

pub fn load_reference(path: &str) -> io::Result<HashMap<String, Vec<u8>>> {
    if path.ends_with(".2bit") {
        load_genome_from_2bit_parallel(path)
    } else {
        load_genome_from_fasta_parallel(path)
    }
}

/// Derive a canonical assembly name from a reference file path by stripping
/// the directory component and any standard genomic file suffixes
/// (case-insensitive, handles stacked extensions like `.fa.gz`).
pub fn derive_assembly_name(path: &str) -> String {
    let mut name = Path::new(path)
        .file_name()
        .and_then(|n| n.to_str())
        .unwrap_or(path)
        .to_string();

    let suffixes = [".gz", ".2bit", ".fasta", ".fna", ".fa", ".fas", ".fsa"];
    loop {
        let lower = name.to_lowercase();
        let mut stripped = false;
        for suffix in &suffixes {
            if lower.ends_with(suffix) {
                name.truncate(name.len() - suffix.len());
                stripped = true;
                break;
            }
        }
        if !stripped {
            break;
        }
    }
    name
}

pub fn process_sequences(
    sequences: Vec<SequenceRecord>,
    genome_map: &HashMap<String, Vec<u8>>,
    map_sequences: bool,
    use_aho_corasick: bool,
    debug_mode: bool,
    assembly_name: Option<&str>,
    remove_duplicates: bool,
    intervals: &IntervalOptions,
    log_level: &LogLevel,
) -> Vec<SequenceRecord> {
    let mut results = sequences;

    let t_validate = std::time::Instant::now();
    validate_sequences(&mut results, genome_map, debug_mode);
    let invalid_count = results.iter().filter(|r| r.validated.is_none()).count();
    let valid_count = results.len() - invalid_count;
    println!(
        "## o Validation: {}/{} validated in {:.2}s — {} need mapping",
        valid_count,
        results.len(),
        t_validate.elapsed().as_secs_f32(),
        invalid_count,
    );

    if map_sequences {
        if use_aho_corasick {
            aho_corasick_search_with_validation(&mut results, genome_map, debug_mode, assembly_name, remove_duplicates);
        } else {
            boyer_moore_search_with_validation(&mut results, genome_map, debug_mode, assembly_name, remove_duplicates);
        }
    }

    // Stamp the reference's assembly onto every record that resolved to
    // coordinates (validated, coordinate-fixed, or remapped) but carries no
    // assembly_id of its own, so output identifiers are emitted in full V2 form
    // (assembly_id:sequence_id:start-end_orient) wherever coordinates are known.
    //
    // Gated on `validated.is_some()`: a record that failed validation still holds
    // its originally-parsed (and now known-wrong) coordinates, so it must NOT be
    // dressed up as a full, resolved V2 identifier — it stays bare.  Records that
    // already name an assembly are left untouched.
    if let Some(name) = assembly_name {
        for record in results.iter_mut() {
            if record.assembly_id.is_none()
                && record.validated.is_some()
                && record.start_1b.is_some()
                && record.end_1b.is_some()
                && record.orient.is_some()
            {
                record.assembly_id = Some(name.to_string());
            }
        }
    }

    analyze_intervals(&mut results, genome_map, intervals, log_level);

    results
}

pub fn parse_stockholm(file_path: &str, is_gzip: bool) -> Result<(Vec<SequenceRecord>, Vec<String>), String> {
    // For non-gzip files we know the compressed == on-disk size and can show a
    // determinate bytes bar.  For gzip the on-disk size is the compressed size
    // so byte counts would overshoot; use a spinner instead.
    let file_size = if !is_gzip {
        std::fs::metadata(file_path).map(|m| m.len()).unwrap_or(0)
    } else {
        0
    };

    let progress: ProgressBar = if file_size > 0 {
        let pb = ProgressBar::new(file_size);
        pb.set_style(
            ProgressStyle::default_bar()
                .template("## o Parsing STK: [{bar:40.cyan/blue}] {bytes}/{total_bytes} ({elapsed})")
                .unwrap()
                .progress_chars("#>-"),
        );
        pb
    } else {
        let pb = ProgressBar::new_spinner();
        pb.set_style(
            ProgressStyle::default_spinner()
                .template("## o Parsing STK: {spinner:.cyan} {msg}")
                .unwrap(),
        );
        pb.set_message("0 records, 0 sequences");
        pb
    };
    progress.enable_steady_tick(Duration::from_millis(100));

    let file = File::open(file_path).map_err(|e| format!("Could not open file {}: {}", file_path, e))?;
    let reader: Box<dyn BufRead> = if is_gzip {
            let decoder = GzDecoder::new(file);
            Box::new(BufReader::new(decoder))
        } else {
            Box::new(BufReader::new(file))
        };
    let mut sequences = Vec::new();
    let mut metadata = Vec::new();
    let mut current_metadata = String::new();
    let mut current_record: Vec<(String,String)> = Vec::new();
    let mut metadata_idx = 0;
    let mut order = 0;
    let mut bytes_read: u64 = 0;

    for line in reader.lines() {
        let line = line.map_err(|e| format!("Error reading line: {}", e))?;
        bytes_read += line.len() as u64 + 1; // +1 approximates the newline
        if file_size > 0 {
            progress.set_position(bytes_read.min(file_size));
        }

        if line.starts_with("#") {
            current_metadata.push_str(&line);
            current_metadata.push('\n');
        } else if line.trim().is_empty() {
            // Ignore blank lines
        } else if line.starts_with("//") {
            if current_record.is_empty() {
                return Err("Unexpected '//' without sequences in the record".to_string());
            }

            metadata.push(current_metadata.clone());
            current_metadata.clear();

            for (name, seq) in current_record.drain(..) {
                let row = SeqRow::from_name_seq(&name, &seq);
                if row.sequence_id.is_none() {
                    println!("Failed to parse identifier: {} ... leaving unchanged", name);
                }
                sequences.push(SequenceRecord::from_seq_row(&row, file_path, order, metadata_idx));
                order += 1;
                // Update spinner message every 500 sequences (gzip path).
                if file_size == 0 && order % 500 == 0 {
                    progress.set_message(format!(
                        "{} records, {} sequences",
                        metadata_idx + 1, order
                    ));
                }
            }

            if file_size == 0 {
                progress.set_message(format!(
                    "{} records, {} sequences",
                    metadata_idx + 1, sequences.len()
                ));
            }

            metadata_idx += 1;
        } else {
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() == 2 {
                current_record.push((parts[0].to_string(), parts[1].to_string()));
            } else {
                return Err(format!("Malformed alignment line: {}", line));
            }
        }
    }

    progress.finish_and_clear();

    if !current_record.is_empty() || !current_metadata.is_empty() {
        return Err("Missing trailing '//' at the end of the Stockholm file".to_string());
    }

    Ok((sequences, metadata))
}

pub fn write_stockholm_output(
    records: &[SequenceRecord],
    metadata: &[String],
    output_path: &str,
    is_gzip: bool,
    append: bool,
) -> io::Result<()> {
    let file = OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(!append)
        .append(append)
        .open(output_path)?;

    let mut writer: Box<dyn Write> = if is_gzip {
        Box::new(BufWriter::new(GzEncoder::new(file, Compression::default())))
    } else {
        Box::new(BufWriter::new(file))
    };

    let mut grouped_records: HashMap<usize, Vec<&SequenceRecord>> = HashMap::new();
    for record in records {
        grouped_records
            .entry(record.metadata_idx)
            .or_default()
            .push(record);
    }

    for (metadata_idx, group) in grouped_records {
        if let Some(metadata_entry) = metadata.get(metadata_idx) {
            write!(writer, "{}", metadata_entry)?;
        }
        for record in group {
            let aligned_seq = record
                .aligned_seq
                .as_ref()
                .map(|seq| String::from_utf8_lossy(seq).to_string())
                .unwrap_or_else(|| String::from_utf8_lossy(&record.sequence).to_string());

            let v2_id = record.format_id();

            if record.start_1b.is_some() && record.end_1b.is_some() && record.orient.is_some() {
                writeln!(writer, "{}:{}-{}_{} {}", v2_id, record.start_1b.unwrap(), record.end_1b.unwrap(),
                        record.orient.unwrap(), aligned_seq)?;
            } else {
                writeln!(writer, "{} {}", v2_id, aligned_seq)?;
            }
        }
        writeln!(writer, "//")?;
    }

    writer.flush()?;
    Ok(())
}

pub fn write_fasta_output(
    records: &[SequenceRecord],
    metadata: &[String],
    output_path: &str,
    is_gzip: bool,
    append: bool,
) -> io::Result<()> {
    let file = OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(!append)
        .append(append)
        .open(output_path)?;

    let mut writer: Box<dyn Write> = if is_gzip {
        Box::new(BufWriter::new(GzEncoder::new(file, Compression::default())))
    } else {
        Box::new(BufWriter::new(file))
    };

    for record in records {
        let v2_id = record.format_id();
        let empty_string = String::new();
        let metadata_entry = metadata.get(record.metadata_idx).unwrap_or(&empty_string);
        if record.start_1b.is_some() && record.end_1b.is_some() && record.orient.is_some() {
            writeln!(writer, ">{}:{}-{}_{} {}", v2_id, record.start_1b.unwrap(), record.end_1b.unwrap(),
                    record.orient.unwrap(), metadata_entry)?;
        } else {
            writeln!(writer, ">{} {}", v2_id, metadata_entry)?;
        }

        writeln!(writer, "{}", String::from_utf8_lossy(&record.sequence))?;
    }

    writer.flush()?;
    Ok(())
}

pub fn write_delimited_output(
    records: &[SequenceRecord],
    output_path: &str,
    is_gzip: bool,
    append: bool,
    format: &str,
) -> io::Result<()> {
    let delimiter = match format {
        "TabDelimited" => '\t',
        "CommaDelimited" => ',',
        _ => {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "Unsupported format. Use 'TabDelimited' or 'CommaDelimited'.",
            ))
        }
    };

    let file = OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(!append)
        .append(append)
        .open(output_path)?;

    let mut writer: Box<dyn Write> = if is_gzip {
        Box::new(BufWriter::new(GzEncoder::new(file, Compression::default())))
    } else {
        Box::new(BufWriter::new(file))
    };

    for record in records {
        let assembly_id = record.assembly_id.clone().unwrap_or_default();
        let sequence_id = record.sequence_id.clone();
        let start = record.start_1b.map(|v| v.to_string()).unwrap_or_default();
        let end = record.end_1b.map(|v| v.to_string()).unwrap_or_default();
        let orient = record.orient.clone().unwrap_or_default();
        let sequence = String::from_utf8_lossy(&record.sequence);

        writeln!(
            writer,
            "{}{}{}{}{}{}{}{}{}{}{}",
            assembly_id, delimiter,
            sequence_id, delimiter,
            start, delimiter,
            end, delimiter,
            orient, delimiter,
            sequence
        )?;
    }

    writer.flush()?;
    Ok(())
}

pub fn detect_format_and_compression(path: &str) -> io::Result<(bool, &'static str)> {
    let file = match File::open(path) {
        Ok(file) => file,
        Err(e) if e.kind() == io::ErrorKind::NotFound => {
            eprintln!("Error: The file '{}' was not found.", path);
            return Err(e);
        }
        Err(e) => {
            eprintln!("Error: Could not open file '{}': {}", path, e);
            return Err(e);
        }
    };

    let mut buf_reader = BufReader::new(file);

    let is_gzip = match buf_reader.fill_buf() {
        Ok(data) => data.starts_with(b"\x1F\x8B"),
        Err(_) => false,
    };

    let mut reader: Box<dyn BufRead> = if is_gzip {
        let file = File::open(path)?;
        Box::new(BufReader::new(GzDecoder::new(file)))
    } else {
        Box::new(buf_reader)
    };

    let mut magic_buffer = [0; 4];
    reader.read_exact(&mut magic_buffer)?;
    let magic_be = u32::from_be_bytes(magic_buffer);
    let magic_le = u32::from_le_bytes(magic_buffer);
    if magic_be == 0x1A412743 || magic_le == 0x1A412743 {
        return Ok((is_gzip, "TwoBit"));
    }

    let mut first_line = String::from_utf8_lossy(&magic_buffer).to_string();
    let mut format = None;
    for _ in 0..15 {
        let mut line = String::new();

        if !first_line.is_empty() {
            line.push_str(&first_line);
            first_line.clear();
        }
        if reader.read_line(&mut line)? == 0 {
            break;
        }

        if line.trim().is_empty() {
            continue;
        }

        if line.starts_with(">") {
            format = Some("Fasta");
            break;
        } else if line.starts_with("# STOCKHOLM 1.0") {
            format = Some("Stockholm");
            break;
        }

        if line.contains('\t') && line.split('\t').count() > 1 {
            format = Some("TabDelimited");
            break;
        } else if line.contains(',') && line.split(',').count() > 1 {
            format = Some("CommaDelimited");
            break;
        }
    }

    match format {
        Some(fmt) => Ok((is_gzip, fmt)),
        None => Err(io::Error::new(io::ErrorKind::InvalidData, "Unknown file format")),
    }
}

pub fn validate_sequences(
    records: &mut [SequenceRecord],
    genome_map: &HashMap<String, Vec<u8>>,
    debug_mode: bool,
) {
    let mut fix_counts: HashMap<String, usize> = HashMap::new();

    for record in records.iter_mut() {
        let genome_sequence = match genome_map.get(&record.sequence_id) {
            Some(seq) => seq,
            None => continue,
        };

        if record.start_1b.is_none() {
            if genome_sequence == &record.sequence {
                record.validated = Some("valid".to_string());
                *fix_counts.entry(record.validated.clone().unwrap()).or_insert(0) += 1;
            }
            continue;
        }

        let range_length = match (record.end_1b, record.start_1b) {
            (Some(end), Some(start)) => end - start,
             _ => 0,
        };
        let fasta_sequence_length = record.sequence.len() as u64;

        let mut validation_str = String::new();
        if range_length == fasta_sequence_length {
            if debug_mode {
                println!(
                    "Detected half-open coordinates for record {}. Converting to one-based fully closed.",
                    record.original_id.as_ref().unwrap()
                );
            }
            validation_str.push_str("_halfopen");
            record.start_1b = record.start_1b.map(|start| start + 1);
        }

        let start = record.start_1b.unwrap() as usize - 1;
        let end = record.end_1b.unwrap() as usize;
        let fasta_sequence = &record.sequence;
        let rev_complement = reverse_complement(fasta_sequence);
        let mut located = false;

        if start < genome_sequence.len() && end <= genome_sequence.len() {
            let mut direct_match_orient: Option<char> = None;
            if &genome_sequence[start..end] == fasta_sequence {
                direct_match_orient = Some('+');
            }
            if &genome_sequence[start..end] == &rev_complement {
                if direct_match_orient.is_none() {
                    direct_match_orient = Some('-');
                }else {
                    direct_match_orient = Some('B');
                }
            }
            if direct_match_orient.is_some() {
                located = true;
                if direct_match_orient == Some('B') || direct_match_orient == record.orient {
                    if debug_mode {
                        println!("Direct match validated for: {:?}", record);
                    }
                } else
                {
                    validation_str.push_str("_orient");
                    record.orient = direct_match_orient;
                }
            }
        }

        if !located {
            let shifts: [isize; 6] = [-3, -2, -1, 1, 2, 3];
            let orig_len = end.saturating_sub(start);
            for shift in shifts.iter() {
                let shifted_start = if *shift < 0 {
                    start.saturating_sub((-*shift) as usize)
                } else {
                    start.saturating_add(*shift as usize)
                };

                let shifted_end = if *shift < 0  {
                    end.saturating_sub((-*shift) as usize)
                } else {
                    end.saturating_add(*shift as usize)
                };
                let new_len = shifted_end.saturating_sub(shifted_start);

                if new_len == orig_len && shifted_end  <= genome_sequence.len() {
                    if &genome_sequence[shifted_start..shifted_end] == fasta_sequence {
                        validation_str.push_str(&format!("{}{}{}",
                            if record.orient == Some('-') { "_orient" } else { "" },
                            if *shift >= 0 { "_plus" } else { "_minus" },
                            shift.abs()));
                        record.start_1b = Some((shifted_start + 1) as u64);
                        record.end_1b = Some(shifted_end as u64);
                        record.orient = if record.orient == Some('-') { Some('+') } else { Some('+') };
                        located = true;
                        break;
                    }

                    if &genome_sequence[shifted_start..shifted_end] == &rev_complement {
                        validation_str.push_str(&format!("{}{}{}",
                            if record.orient == Some('+') { "_orient" } else { "" },
                            if *shift >= 0 { "_plus" } else { "_minus" },
                            shift.abs()));
                        record.start_1b = Some((shifted_start + 1) as u64);
                        record.end_1b = Some(shifted_end as u64);
                        record.orient = if record.orient == Some('+') { Some('-') } else { Some('-') };
                        located = true;
                        break;
                    }
                }
            }
        }

        if located {
            if validation_str.is_empty() {
                record.validated = Some("valid".to_string());
            } else {
                record.validated = Some(format!("fixed{}",validation_str));
            }
            *fix_counts.entry(record.validated.clone().unwrap()).or_insert(0) += 1;
        }
    }
}

pub fn reverse_complement(dna: &[u8]) -> Vec<u8> {
    dna.iter()
        .rev()
        .map(|&base| match base {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            _ => base,
        })
        .collect()
}

pub fn parse_fasta(file_path: &str, is_gzip: bool) -> (Vec<SequenceRecord>, Vec<String>) {
    let file = File::open(file_path).expect("Could not open Fasta file");
    let reader: Box<dyn BufRead> = if is_gzip {
        let decoder = GzDecoder::new(file);
        Box::new(BufReader::new(decoder))
    } else {
        Box::new(BufReader::new(file))
    };

    let mut sequences = Vec::new();
    let mut metadata = Vec::new();
    let mut order = 0;

    let mut current_sequence: Option<SequenceRecord> = None;

    for line in reader.lines() {
        let line = line.expect("Error reading Fasta file");
        if line.starts_with('>') {
            if let Some(record) = current_sequence.take() {
                sequences.push(record);
            }

            let header = line[1..].trim();
            let mut parts = header.splitn(2, char::is_whitespace);
            let orig_id = parts.next().unwrap().to_string();
            let metadata_entry = parts.next().unwrap_or("").to_string();

            metadata.push(metadata_entry.clone());

            // Parse the identifier via SeqRow; sequence bytes are appended below.
            let row = SeqRow::from_name_seq(&orig_id, "");
            current_sequence = Some(SequenceRecord::from_seq_row(&row, file_path, order, order));

            order += 1;
        } else if let Some(record) = current_sequence.as_mut() {
            record.sequence.extend(line.trim().bytes());
        }
    }
    if let Some(record) = current_sequence.take() {
        sequences.push(record);
    }

    (sequences, metadata)
}

pub fn parse_delimited_file(file_path: &str) -> (Vec<SequenceRecord>, Vec<String>) {
    let file = File::open(file_path).expect("Could not open delimited file");
    let reader = BufReader::new(file);
    let mut records = Vec::new();
    let mut order = 0;

    for line in reader.lines() {
        let line = line.expect("Could not read line");
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() >= 2 {
            let id = parts[0].to_string();
            let seq = parts[1].to_string().into_bytes();

            let row = SeqRow::from_name_seq(&id, "");
            let mut record = SequenceRecord::from_seq_row(&row, file_path, order, 0);
            record.sequence = seq;

            order += 1;
            records.push(record);
        }
    }

    (records, vec![])
}

pub fn boyer_moore_search_with_validation(
    records: &mut [SequenceRecord],
    genome_map: &HashMap<String, Vec<u8>>,
    debug_mode: bool,
    remapped_assembly: Option<&str>,
    remove_duplicates: bool,
) {
    let invalid_indices: Vec<usize> = records
        .iter()
        .enumerate()
        .filter(|(_, r)| r.validated.is_none())
        .map(|(i, _)| i)
        .collect();

    let invalid_count = invalid_indices.len();
    let t_map = std::time::Instant::now();
    let progress = ProgressBar::new(invalid_count as u64);
    progress.set_style(
        ProgressStyle::default_bar()
            .template("## o Mapping: {spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {pos}/{len} ({eta})")
            .unwrap()
            .progress_chars("#>-"),
    );
    progress.enable_steady_tick(Duration::from_millis(100));

    // Phase 1 (parallel): scan the genome for every unvalidated record.
    // Results are (record_index, sorted_hits); order within the Vec is arbitrary.
    let mut all_hits: Vec<(usize, Vec<(usize, char, String)>)> = invalid_indices
        .par_iter()
        .map(|&idx| {
            let record = &records[idx];
            let pattern = &record.sequence;
            let rev_complement_pattern = reverse_complement(pattern);
            let original_sequence_id = record.sequence_id.clone();
            let mut found_positions: Vec<(usize, char, String)> = Vec::new();

            if let Some(target_sequence) = genome_map.get(&original_sequence_id) {
                found_positions.extend(boyer_moore_search(target_sequence, pattern)
                    .into_iter().map(|pos| (pos, '+', original_sequence_id.clone())));
                found_positions.extend(boyer_moore_search(target_sequence, &rev_complement_pattern)
                    .into_iter().map(|pos| (pos, '-', original_sequence_id.clone())));
            }

            if found_positions.is_empty() {
                for (seq_name, genome_sequence) in genome_map {
                    if seq_name == &original_sequence_id { continue; }
                    found_positions.extend(boyer_moore_search(genome_sequence, pattern)
                        .into_iter().map(|pos| (pos, '+', seq_name.clone())));
                    found_positions.extend(boyer_moore_search(genome_sequence, &rev_complement_pattern)
                        .into_iter().map(|pos| (pos, '-', seq_name.clone())));
                }
            }

            found_positions.sort_by(|a, b| {
                let a_same = a.2 == original_sequence_id;
                let b_same = b.2 == original_sequence_id;
                b_same.cmp(&a_same)
                    .then_with(|| {
                        let dist_a = record.start_1b.map_or(usize::MAX, |start| (a.0 as isize - start as isize).unsigned_abs());
                        let dist_b = record.start_1b.map_or(usize::MAX, |start| (b.0 as isize - start as isize).unsigned_abs());
                        dist_a.cmp(&dist_b)
                    })
                    .then_with(|| a.2.cmp(&b.2))
                    .then_with(|| a.1.cmp(&b.1))
            });

            progress.inc(1);
            (idx, found_positions)
        })
        .collect();

    progress.finish_and_clear();

    // Phase 2 (sequential, file order): claim positions deterministically.
    // Earlier records in the file get first pick; later ones are steered away
    // from already-claimed positions when alternatives exist.
    all_hits.sort_by_key(|(idx, _)| *idx);

    let mut occupied: std::collections::HashSet<(String, u64, u64, char)> = records
        .iter()
        .filter(|r| r.validated.is_some())
        .filter_map(|r| match (r.start_1b, r.end_1b, r.orient) {
            (Some(s), Some(e), Some(o)) => Some((r.sequence_id.clone(), s, e, o)),
            _ => None,
        })
        .collect();

    let mut mapped_count = 0usize;
    for (idx, found_positions) in all_hits {
        let record = &mut records[idx];
        if found_positions.is_empty() {
            continue;
        }
        let pat_len = record.sequence.len();
        let total_hits = found_positions.len();

        let first_free = found_positions.iter().find(|hit| {
            let s = hit.0 as u64 + 1;
            let e = (hit.0 + pat_len) as u64;
            !occupied.contains(&(hit.2.clone(), s, e, hit.1))
        });

        let chosen = match first_free {
            Some(hit) => {
                let s = hit.0 as u64 + 1;
                let e = (hit.0 + pat_len) as u64;
                occupied.insert((hit.2.clone(), s, e, hit.1));
                Some(hit.clone())
            }
            None => {
                if remove_duplicates { None } else { Some(found_positions[0].clone()) }
            }
        };

        match chosen {
            None => {
                record.validated = Some("removed_remapped_duplicate".to_string());
            }
            Some(best) => {
                record.start_1b = Some(best.0 as u64 + 1);
                record.end_1b = Some((best.0 + pat_len) as u64);
                record.orient = Some(best.1);
                record.sequence_id = best.2.clone();
                record.assembly_id = remapped_assembly.map(|s| s.to_string());
                record.validated = Some(if total_hits == 1 {
                    "fixed_remapped_unique".to_string()
                } else {
                    "fixed_remapped_ambig".to_string()
                });
                mapped_count += 1;
            }
        }

        if debug_mode {
            if let Some(v) = &record.validated {
                println!("Boyer-Moore {} fix for record: {:?}", v, record);
            } else {
                println!("Boyer-Moore failed to fix for record: {:?}", record);
            }
        }
    }
    println!("## o Mapping: {}/{} mapped in {:.2}s", mapped_count, invalid_count, t_map.elapsed().as_secs_f32());
}

/// Sort a hit list for one sequence: hits on the original chromosome first,
/// then by distance from the original coordinates, then by chromosome name and
/// position for deterministic tie-breaking.
fn sort_hits_by_proximity(
    hits: &mut Vec<(usize, char, String)>,
    original_seq_id: &str,
    record_start: Option<u64>,
) {
    hits.sort_by(|a, b| {
        b.2.eq(original_seq_id).cmp(&a.2.eq(original_seq_id))
            .then_with(|| {
                let dist_a = record_start.map_or(usize::MAX, |s| (a.0 as isize - s as isize).unsigned_abs());
                let dist_b = record_start.map_or(usize::MAX, |s| (b.0 as isize - s as isize).unsigned_abs());
                dist_a.cmp(&dist_b)
            })
            .then_with(|| a.2.cmp(&b.2))
            .then_with(|| a.1.cmp(&b.1))
    });
}

/// Scan every chromosome in `genome_map` for all sequences in `batch_indices`
/// simultaneously, using a single Aho-Corasick automaton built from all their
/// patterns (forward + reverse-complement).
///
/// This is O(genome_size × num_chromosomes) regardless of how many sequences
/// are in the batch — orders of magnitude faster than per-sequence whole-genome
/// searches when many sequences need mapping.
///
/// Returns `(record_index, sorted_hits)` for every sequence in the batch.
fn batch_genome_scan(
    batch_indices: &[usize],
    records: &[SequenceRecord],
    genome_map: &HashMap<String, Vec<u8>>,
    progress: &ProgressBar,
) -> Vec<(usize, Vec<(usize, char, String)>)> {
    if batch_indices.is_empty() {
        return Vec::new();
    }

    // Pattern layout: 2k = forward of batch[k], 2k+1 = reverse-complement.
    let patterns: Vec<Vec<u8>> = batch_indices
        .iter()
        .flat_map(|&idx| {
            let seq = &records[idx].sequence;
            [seq.clone(), reverse_complement(seq)]
        })
        .collect();

    let ac = AhoCorasick::new(patterns.iter().map(|p| p.as_slice()))
        .expect("failed to build batch Aho-Corasick automaton");

    // Search every chromosome in parallel; collect raw (batch_pos, start, strand, chrom) tuples.
    let raw_hits: Vec<(usize, usize, char, String)> = genome_map
        .par_iter()
        .flat_map_iter(|(chrom_name, chrom_seq)| {
            let mut local = Vec::new();
            for m in ac.find_overlapping_iter(chrom_seq) {
                let pat_idx = m.pattern().as_usize();
                let batch_pos = pat_idx / 2;
                let strand = if pat_idx % 2 == 0 { '+' } else { '-' };
                local.push((batch_pos, m.start(), strand, chrom_name.clone()));
            }
            local
        })
        .collect();

    // Group hits by their position within the batch.
    let mut hits_by_pos: Vec<Vec<(usize, char, String)>> = vec![Vec::new(); batch_indices.len()];
    for (batch_pos, start, strand, chrom) in raw_hits {
        if batch_pos < hits_by_pos.len() {
            hits_by_pos[batch_pos].push((start, strand, chrom));
        }
    }

    // Sort each sequence's hits and advance the progress bar.
    let mut results = Vec::with_capacity(batch_indices.len());
    for (bp, &idx) in batch_indices.iter().enumerate() {
        let record = &records[idx];
        let mut hits = std::mem::take(&mut hits_by_pos[bp]);
        sort_hits_by_proximity(&mut hits, &record.sequence_id, record.start_1b);
        results.push((idx, hits));
        progress.inc(1);
    }
    results
}

/// Aho-Corasick variant of the mapping step.
///
/// Identical selection logic to `boyer_moore_search_with_validation` but searches
/// both the forward pattern and its reverse complement in a **single pass** over
/// each chromosome, halving the genome I/O compared to two separate Boyer-Moore
/// calls.
pub fn aho_corasick_search_with_validation(
    records: &mut [SequenceRecord],
    genome_map: &HashMap<String, Vec<u8>>,
    debug_mode: bool,
    remapped_assembly: Option<&str>,
    remove_duplicates: bool,
) {
    let invalid_indices: Vec<usize> = records
        .iter()
        .enumerate()
        .filter(|(_, r)| r.validated.is_none())
        .map(|(i, _)| i)
        .collect();

    let invalid_count = invalid_indices.len();
    let t_map = std::time::Instant::now();
    let progress = ProgressBar::new(invalid_count as u64);
    progress.set_style(
        ProgressStyle::default_bar()
            .template("## o Mapping: {spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {pos}/{len} ({eta})")
            .unwrap()
            .progress_chars("#>-"),
    );
    progress.enable_steady_tick(Duration::from_millis(100));

    // Partition: sequences whose chromosome IS in genome_map can use a fast
    // per-sequence search (they usually hit on their own chromosome right away).
    // Sequences whose chromosome is NOT in genome_map must scan the whole genome;
    // batching them lets one Aho-Corasick pass cover BATCH_SIZE sequences at once,
    // reducing total genome traversals from N×chromosomes to ⌈N/BATCH_SIZE⌉×chromosomes.
    let (chrom_found, chrom_missing): (Vec<usize>, Vec<usize>) = invalid_indices
        .iter()
        .partition(|&&idx| genome_map.contains_key(&records[idx].sequence_id));

    // Phase 1a (parallel by sequence): per-sequence search for records whose
    // chromosome exists in the genome.  Whole-genome fallback is kept for the
    // rare case where the target chromosome yields no hit.
    let mut all_hits: Vec<(usize, Vec<(usize, char, String)>)> = chrom_found
        .par_iter()
        .map(|&idx| {
            let record = &records[idx];
            let pattern = &record.sequence;
            let rev_complement_pattern = reverse_complement(pattern);
            let original_sequence_id = record.sequence_id.clone();

            let ac = AhoCorasick::new([pattern.as_slice(), rev_complement_pattern.as_slice()])
                .expect("failed to build Aho-Corasick automaton");

            let mut found_positions: Vec<(usize, char, String)> = Vec::new();

            if let Some(target_sequence) = genome_map.get(&original_sequence_id) {
                for m in ac.find_overlapping_iter(target_sequence) {
                    let strand = if m.pattern().as_usize() == 0 { '+' } else { '-' };
                    found_positions.push((m.start(), strand, original_sequence_id.clone()));
                }
            }

            if found_positions.is_empty() {
                for (seq_name, genome_sequence) in genome_map {
                    if seq_name == &original_sequence_id { continue; }
                    for m in ac.find_overlapping_iter(genome_sequence) {
                        let strand = if m.pattern().as_usize() == 0 { '+' } else { '-' };
                        found_positions.push((m.start(), strand, seq_name.clone()));
                    }
                }
            }

            sort_hits_by_proximity(&mut found_positions, &original_sequence_id, record.start_1b);
            progress.inc(1);
            (idx, found_positions)
        })
        .collect();

    // Phase 1b: batch whole-genome scan for sequences whose chromosome is absent
    // from genome_map.  Each batch builds one multi-pattern AC automaton and
    // searches every chromosome once, so the cost is O(batch_count × genome_size)
    // rather than O(sequence_count × genome_size).
    const BATCH_SIZE: usize = 1000;
    for batch in chrom_missing.chunks(BATCH_SIZE) {
        let batch_results = batch_genome_scan(batch, records, genome_map, &progress);
        all_hits.extend(batch_results);
    }

    progress.finish_and_clear();

    // Phase 2 (sequential, file order): claim positions deterministically.
    // Earlier records in the file get first pick; later ones are steered away
    // from already-claimed positions when alternatives exist.
    all_hits.sort_by_key(|(idx, _)| *idx);

    let mut occupied: std::collections::HashSet<(String, u64, u64, char)> = records
        .iter()
        .filter(|r| r.validated.is_some())
        .filter_map(|r| match (r.start_1b, r.end_1b, r.orient) {
            (Some(s), Some(e), Some(o)) => Some((r.sequence_id.clone(), s, e, o)),
            _ => None,
        })
        .collect();

    let mut mapped_count = 0usize;
    for (idx, found_positions) in all_hits {
        let record = &mut records[idx];
        if found_positions.is_empty() {
            continue;
        }
        let pat_len = record.sequence.len();
        let total_hits = found_positions.len();

        let first_free = found_positions.iter().find(|hit| {
            let s = hit.0 as u64 + 1;
            let e = (hit.0 + pat_len) as u64;
            !occupied.contains(&(hit.2.clone(), s, e, hit.1))
        });

        let chosen = match first_free {
            Some(hit) => {
                let s = hit.0 as u64 + 1;
                let e = (hit.0 + pat_len) as u64;
                occupied.insert((hit.2.clone(), s, e, hit.1));
                Some(hit.clone())
            }
            None => {
                if remove_duplicates { None } else { Some(found_positions[0].clone()) }
            }
        };

        match chosen {
            None => {
                record.validated = Some("removed_remapped_duplicate".to_string());
            }
            Some(best) => {
                record.start_1b = Some(best.0 as u64 + 1);
                record.end_1b = Some((best.0 + pat_len) as u64);
                record.orient = Some(best.1);
                record.sequence_id = best.2.clone();
                record.assembly_id = remapped_assembly.map(|s| s.to_string());
                record.validated = Some(if total_hits == 1 {
                    "fixed_remapped_unique".to_string()
                } else {
                    "fixed_remapped_ambig".to_string()
                });
                mapped_count += 1;
            }
        }

        if debug_mode {
            if let Some(v) = &record.validated {
                println!("AhoCorasick {} fix for record: {:?}", v, record);
            } else {
                println!("AhoCorasick failed to fix for record: {:?}", record);
            }
        }
    }
    println!("## o Mapping: {}/{} mapped in {:.2}s", mapped_count, invalid_count, t_map.elapsed().as_secs_f32());
}

pub fn boyer_moore_search(text: &[u8], pattern: &[u8]) -> Vec<usize> {
    let m = pattern.len();
    let n = text.len();
    if m == 0 || m > n {
        return vec![];
    }

    let mut bad_char = [-1; 256];
    bad_char_heuristic(pattern, &mut bad_char);

    let mut positions = Vec::new();
    let mut s = 0;

    while s <= n - m {
        let mut j = (m - 1) as isize;

        while j >= 0 && pattern[j as usize] == text[s + j as usize] {
            j -= 1;
        }

        if j < 0 {
            positions.push(s);
            s += if s + m < n { m.saturating_sub(bad_char[text[s + m] as usize].max(0) as usize) } else { 1 };
        } else {
            s += (j - bad_char[text[s + j as usize] as usize]).max(1) as usize;
        }
    }

    positions
}

fn bad_char_heuristic(pattern: &[u8], bad_char: &mut [isize; 256]) {
    for i in 0..256 {
        bad_char[i] = -1;
    }
    for (i, &ch) in pattern.iter().enumerate() {
        bad_char[ch as usize] = i as isize;
    }
}

/// Overlap, in base pairs, that two spans must share before they are allowed to
/// join the same merge cluster.  `min_overlap` below 1 is read as a fraction of
/// the shorter span; 1 or more is read as a literal base-pair count.
fn required_overlap(len_a: u64, len_b: u64, min_overlap: f64) -> u64 {
    if min_overlap < 1.0 {
        let shorter = len_a.min(len_b) as f64;
        (((min_overlap.max(0.0)) * shorter).ceil() as u64).max(1)
    } else {
        min_overlap as u64
    }
}

/// Interval-level redundancy handling, run after validation and mapping.
///
/// Records are grouped by (input file, sequence id), so separate families and
/// separate reference sequences never interact.  Strand is ignored: two records
/// covering the same bases are redundant whichever way round they were written,
/// and a merged cluster takes the orientation of its longest member.
///
/// This always counts and reports.  `opts` only decides what it acts on, so
/// discoord tells you what is there whether or not you asked it to do anything
/// about it.
pub fn analyze_intervals(
    records: &mut [SequenceRecord],
    genome_map: &HashMap<String, Vec<u8>>,
    opts: &IntervalOptions,
    log_level: &LogLevel,
) {
    // A record earns a span only if its coordinates were resolved.  Anything
    // still carrying its originally-parsed (and known-wrong) numbers would
    // poison the comparison, so it sits the pass out.
    let spans: Vec<Option<(u64, u64, char)>> = records
        .iter()
        .map(|r| match (r.start_1b, r.end_1b, r.orient, r.validated.as_deref()) {
            (Some(s), Some(e), Some(o), Some(v))
                if s <= e && v != "invalid" && v != "removed_remapped_duplicate" =>
            {
                Some((s, e, o))
            }
            _ => None,
        })
        .collect();

    let mut groups: HashMap<(String, String), Vec<usize>> = HashMap::new();
    for (i, span) in spans.iter().enumerate() {
        if span.is_none() {
            continue;
        }
        groups
            .entry((records[i].input_file.clone(), records[i].sequence_id.clone()))
            .or_default()
            .push(i);
    }
    let eligible: usize = groups.values().map(|v| v.len()).sum();
    if eligible == 0 {
        return;
    }
    let mut keys: Vec<(String, String)> = groups.keys().cloned().collect();
    keys.sort();

    let mut contained: Vec<(usize, usize)> = Vec::new();
    let mut overlap_pairs: u64 = 0;
    let mut clusters: Vec<Vec<usize>> = Vec::new();

    for key in &keys {
        let idxs = &groups[key];

        // Start ascending, longest first at a tie, so a container is always
        // seen before anything it swallows; equal spans fall back to input
        // order, which keeps the earliest copy.
        let mut by_start = idxs.clone();
        by_start.sort_by_key(|&i| {
            let (s, e, _) = spans[i].unwrap();
            (s, std::cmp::Reverse(e), records[i].order)
        });

        let mut contained_here: std::collections::HashSet<usize> = std::collections::HashSet::new();
        let mut reach: Option<(u64, usize)> = None;
        for &i in &by_start {
            let (_, e, _) = spans[i].unwrap();
            match reach {
                // The record holding `reach` cannot itself be contained, so the
                // container reported here always survives the pass.
                Some((max_end, container)) if e <= max_end => {
                    contained.push((i, container));
                    contained_here.insert(i);
                }
                _ => reach = Some((e, i)),
            }
        }

        // Exact pair count by sweep: retire everything that ends before this
        // record starts, and whatever is still open overlaps it.
        let mut open: std::collections::BinaryHeap<std::cmp::Reverse<u64>> =
            std::collections::BinaryHeap::new();
        for &i in &by_start {
            let (s, e, _) = spans[i].unwrap();
            while let Some(&std::cmp::Reverse(end)) = open.peek() {
                if end < s {
                    open.pop();
                } else {
                    break;
                }
            }
            overlap_pairs += open.len() as u64;
            open.push(std::cmp::Reverse(e));
        }

        // Build clusters whether or not merging is switched on, so the report
        // can say what merging would do.  When containment removal is also on,
        // the records it drops are already out of the running.
        let mut live: Vec<usize> = idxs
            .iter()
            .copied()
            .filter(|i| !(opts.remove_contained && contained_here.contains(i)))
            .collect();
        live.sort_by_key(|&i| {
            let (s, e, _) = spans[i].unwrap();
            (s, e, records[i].order)
        });

        let mut cluster: Vec<usize> = Vec::new();
        let (mut cl_start, mut cl_end) = (0u64, 0u64);
        for &i in &live {
            let (s, e, _) = spans[i].unwrap();
            if cluster.is_empty() {
                cluster.push(i);
                cl_start = s;
                cl_end = e;
                continue;
            }
            // Measured against the cluster's running span, not just the previous
            // record, so a chain of tiled windows stays one cluster.
            let shared = (cl_end.min(e) as i64) - (cl_start.max(s) as i64) + 1;
            let needed = required_overlap(cl_end - cl_start + 1, e - s + 1, opts.min_overlap) as i64;
            if shared >= needed {
                cluster.push(i);
                cl_end = cl_end.max(e);
            } else {
                if cluster.len() > 1 {
                    clusters.push(cluster.clone());
                }
                cluster.clear();
                cluster.push(i);
                cl_start = s;
                cl_end = e;
            }
        }
        if cluster.len() > 1 {
            clusters.push(cluster);
        }
    }

    let clustered_records: usize = clusters.iter().map(|c| c.len()).sum();

    // ---- act ----------------------------------------------------------------
    let mut merged_detail: Vec<(String, u64, u64, char, usize)> = Vec::new();
    let mut alignment_lost = false;
    if opts.merge_overlapping {
        for cluster in &clusters {
            // The longest member sets orientation and keeps its description;
            // ties go to whichever came first in the file.
            let rep = *cluster
                .iter()
                .max_by_key(|&&i| {
                    let (s, e, _) = spans[i].unwrap();
                    (e - s, std::cmp::Reverse(records[i].order))
                })
                .unwrap();
            let new_start = cluster.iter().map(|&i| spans[i].unwrap().0).min().unwrap();
            let new_end = cluster.iter().map(|&i| spans[i].unwrap().1).max().unwrap();
            // Only worth warning about when the input actually carried an
            // alignment for merging to invalidate.
            let had_alignment = cluster.iter().any(|&i| records[i].aligned_seq.is_some());
            let new_order = cluster.iter().map(|&i| records[i].order).min().unwrap();
            let orient = records[rep].orient.unwrap_or('+');

            let sequence = match genome_map.get(&records[rep].sequence_id) {
                Some(g) if new_start >= 1 && new_end as usize <= g.len() => {
                    let sub = &g[(new_start - 1) as usize..new_end as usize];
                    if orient == '-' {
                        reverse_complement(sub)
                    } else {
                        sub.to_vec()
                    }
                }
                _ => {
                    eprintln!(
                        "## Warning: merged span {}:{}-{} lies outside the reference; cluster left unmerged",
                        records[rep].sequence_id, new_start, new_end
                    );
                    continue;
                }
            };

            for &i in cluster {
                if i != rep {
                    records[i].validated = Some("removed_merged".to_string());
                }
            }
            let sequence_id = records[rep].sequence_id.clone();
            let r = &mut records[rep];
            r.start_1b = Some(new_start);
            r.end_1b = Some(new_end);
            r.orient = Some(orient);
            r.order = new_order;
            r.sequence = sequence;
            // A merged span has no alignment columns; Stockholm output falls
            // back to the ungapped sequence, which no longer fits the MSA.
            r.aligned_seq = None;
            r.validated = Some("merged_overlapping".to_string());
            if had_alignment {
                alignment_lost = true;
            }
            merged_detail.push((sequence_id, new_start, new_end, orient, cluster.len()));
        }
    }

    if opts.remove_contained {
        for &(i, _) in &contained {
            records[i].validated = Some("removed_contained".to_string());
        }
    }

    // ---- report -------------------------------------------------------------
    println!(
        "## o Intervals: {} record(s) with coordinates in {} group(s)",
        eligible,
        keys.len()
    );
    println!("##     Overlapping pairs: {}", overlap_pairs);
    println!(
        "##     Contained records: {} ({})",
        contained.len(),
        if opts.remove_contained { "removed" } else { "reported only" }
    );
    println!(
        "##     Merge clusters: {} covering {} record(s) ({})",
        clusters.len(),
        clustered_records,
        if opts.merge_overlapping {
            format!("{} merged record(s) written", merged_detail.len())
        } else {
            "reported only".to_string()
        }
    );
    if alignment_lost {
        println!("##     Note: merged records dropped their alignment columns; the output is no longer a valid MSA");
    }

    if *log_level != LogLevel::Summary {
        if !contained.is_empty() {
            println!("##     Contained detail:");
            for &(i, container) in &contained {
                let (s, e, o) = spans[i].unwrap();
                let (cs, ce, co) = spans[container].unwrap();
                println!(
                    "##       {}:{}-{}_{} inside {}:{}-{}_{}",
                    records[i].sequence_id, s, e, o,
                    records[container].sequence_id, cs, ce, co
                );
            }
        }
        if !merged_detail.is_empty() {
            println!("##     Merge detail:");
            for (seq_id, s, e, orient, count) in &merged_detail {
                println!("##       {}:{}-{}_{} <- {} records", seq_id, s, e, orient, count);
            }
        }
    }
}

pub fn output_results(records: &[SequenceRecord], format: LogLevel, label: String) {
    match format {
        LogLevel::Summary | LogLevel::PerRecord => {
            let total_records = records.len();
            let mut fix_counts = HashMap::new();
            let mut fixed_count = 0;
            for record in records.iter() {
                let v = record.validated.as_deref();
                if v.is_some()
                    && !matches!(
                        v,
                        Some("valid")
                            | Some("invalid")
                            | Some("removed_remapped_duplicate")
                            | Some("removed_contained")
                            | Some("removed_merged")
                            | Some("merged_overlapping")
                    )
                {
                    fixed_count += 1;
                    *fix_counts.entry(record.validated.clone().unwrap()).or_insert(0) += 1;
                }
            }
            let valid_count = records.iter().filter(|r| r.validated.as_deref() == Some("valid")).count();
            let invalid_count = records.iter().filter(|r| r.validated.as_deref() == Some("invalid")).count();
            let count_of = |status: &str| {
                records.iter().filter(|r| r.validated.as_deref() == Some(status)).count()
            };
            let removed_dup_count = count_of("removed_remapped_duplicate");
            let removed_contained_count = count_of("removed_contained");
            let absorbed_count = count_of("removed_merged");
            let merged_count = count_of("merged_overlapping");

            println!("{}:", label);
            println!("  Total Sequences: {}", total_records);
            println!("     Accurate Coordinates: {}", valid_count);
            println!("     Repaired Coordinates: {}", fixed_count);
            for (fix_type, count) in fix_counts {
                println!("        {}: {}", fix_type, count);
            }
            if removed_dup_count > 0 {
                println!("     Removed Duplicate Sequences: {}", removed_dup_count);
            }
            if removed_contained_count > 0 {
                println!("     Removed Contained Records: {}", removed_contained_count);
            }
            if merged_count > 0 {
                println!(
                    "     Merged Records: {} (absorbing {} others)",
                    merged_count, absorbed_count
                );
            }
            println!("     Invalid Coordinates: {}", invalid_count);
        }
        LogLevel::Detailed => {
            println!("Detailed Report:");
            for record in records {
                record.print_record();
            }
        }
    }
}

/// Convert `(name, aligned_seq)` rows from a parsed Stockholm record into
/// `SequenceRecord`s suitable for `validate_sequences`.
///
/// Rows whose identifiers cannot be parsed by Smitten (e.g. bare consensus
/// labels) are silently skipped — they have no genomic coordinates to check.
/// Gap characters (`.` and `-`) are stripped from the sequence before
/// validation.
pub fn records_from_rows(rows: &[SeqRow], file_label: &str) -> Vec<SequenceRecord> {
    rows.iter()
        .enumerate()
        .filter(|(_, row)| row.sequence_id.is_some()) // skip unparseable identifiers
        .map(|(order, row)| SequenceRecord::from_seq_row(row, file_label, order, 0))
        .collect()
}

pub fn load_genome_from_fasta_parallel(path: &str) -> io::Result<HashMap<String, Vec<u8>>> {
    let reader = fasta::Reader::from_file(path).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    let sequences: Vec<(String, Vec<u8>)> = reader
        .records()
        .par_bridge()
        .map(|result| {
            let record = result.map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
            let name = record.id().to_string();
            let sequence = record.seq().to_ascii_uppercase().to_vec();
            Ok((name, sequence))
        })
        .collect::<Result<_, io::Error>>()?;

    Ok(sequences.into_iter().collect())
}

pub fn load_genome_from_2bit_parallel(path: &str) -> io::Result<HashMap<String, Vec<u8>>> {
    let file = File::open(path)?;
    let mmap = unsafe { Mmap::map(&file)? };
    let buffer = &mmap[..];

    let is_little_endian = match u32::from_be_bytes(buffer[0..4].try_into().unwrap()) {
        0x1A412743 => false,
        0x4327411A => true,
        _ => return Err(io::Error::new(io::ErrorKind::InvalidData, "Invalid 2bit signature")),
    };

    let read_u32 = |offset: usize| {
        let bytes: [u8; 4] = buffer[offset..offset + 4].try_into().unwrap();
        if is_little_endian { u32::from_le_bytes(bytes) } else { u32::from_be_bytes(bytes) }
    };

    let read_u64 = |offset: usize| {
        let bytes: [u8; 8] = buffer[offset..offset + 8].try_into().unwrap();
        if is_little_endian { u64::from_le_bytes(bytes) } else { u64::from_be_bytes(bytes) }
    };

    let version = read_u32(4);
    if version > 1 {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "Unsupported 2bit version"));
    }

    let seq_count = read_u32(8) as usize;

    let mut sequences = Vec::new();
    let mut offset = 16;

    for _ in 0..seq_count {
        let name_len = buffer[offset] as usize;
        offset += 1;

        let name = String::from_utf8(buffer[offset..offset + name_len].to_vec()).unwrap();
        offset += name_len;

        let seq_offset = if version == 0 {
            read_u32(offset) as u64
        } else {
            read_u64(offset)
        };

        offset += if version == 0 { 4 } else { 8 };

        sequences.push((name, seq_offset));
    }

    let genome_map: HashMap<String, Vec<u8>> = sequences
        .into_par_iter()
        .map(|(name, seq_offset)| {
            let dna_size = read_u32(seq_offset as usize) as usize;

            let n_block_count = read_u32((seq_offset + 4) as usize) as usize;
            let mut n_block_starts = Vec::with_capacity(n_block_count);
            let mut n_block_sizes = Vec::with_capacity(n_block_count);

            let mut current_offset = (seq_offset + 8) as usize;

            for _ in 0..n_block_count {
                let start = read_u32(current_offset) as usize;
                n_block_starts.push(start);
                current_offset += 4;
            }

            for _ in 0..n_block_count {
                let size = read_u32(current_offset) as usize;
                n_block_sizes.push(size);
                current_offset += 4;
            }

            let mask_block_count = read_u32((current_offset) as usize) as usize;
            current_offset = current_offset + (mask_block_count * 8) + 4;

            current_offset += 4;

            let mut genome = vec![b'N'; dna_size];
            for i in 0..((dna_size + 3) / 4) {
                let byte = buffer[current_offset + i];
                for j in 0..4 {
                    let pos = i * 4 + j;
                    if pos >= dna_size {
                        break;
                    }
                    genome[pos] = match (byte >> ((3 - j) * 2)) & 0b11 {
                        0 => b'T',
                        1 => b'C',
                        2 => b'A',
                        3 => b'G',
                        _ => b'N',
                    };
                }
            }

            for (&start, &size) in n_block_starts.iter().zip(n_block_sizes.iter()) {
                for pos in start..(start + size) {
                    if pos < genome.len() {
                        genome[pos] = b'N';
                    }
                }
            }

            (name, genome)
        })
        .collect();

    Ok(genome_map)
}

/// Initialise the global Rayon thread pool to use exactly `n` threads.
///
/// Call this once at program start before any parallel work is submitted.
/// Panics if the global pool has already been initialised.
pub fn init_thread_pool(n: usize) {
    rayon::ThreadPoolBuilder::new()
        .num_threads(n)
        .build_global()
        .expect("Failed to build global thread pool");
}

#[cfg(test)]
mod coordinate_tests {
    use super::*;

    /// `SeqRow` stores a half-open span; this crate audits identifiers in
    /// their written, 1-based closed form. The round trip must be exact.
    #[test]
    fn record_coordinates_are_as_written_in_the_identifier() {
        let row = SeqRow::from_name_seq("chr1:101-200_+", "ACGT");
        let rec = SequenceRecord::from_seq_row(&row, "x.stk", 0, 0);
        assert_eq!((rec.start_1b, rec.end_1b, rec.orient), (Some(101), Some(200), Some('+')));

        let bare = SequenceRecord::from_seq_row(&SeqRow::from_name_seq("consensus", "ACGT"), "x.stk", 0, 0);
        assert_eq!((bare.start_1b, bare.end_1b), (None, None));
    }
}

#[cfg(test)]
mod interval_tests {
    use super::*;

    fn rec(file: &str, seq_id: &str, s: u64, e: u64, o: char, order: usize) -> SequenceRecord {
        SequenceRecord {
            input_file: file.to_string(),
            metadata_idx: 0,
            order,
            original_id: Some(format!("{}:{}-{}_{}", seq_id, s, e, o)),
            assembly_id: None,
            sequence_id: seq_id.to_string(),
            start_1b: Some(s),
            end_1b: Some(e),
            orient: Some(o),
            inferred_version: None,
            sequence: Vec::new(),
            aligned_seq: None,
            validated: Some("valid".to_string()),
        }
    }

    fn genome(len: usize) -> HashMap<String, Vec<u8>> {
        let bases = b"ACGT";
        let mut x: u64 = 12345;
        let seq: Vec<u8> = (0..len)
            .map(|_| {
                x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                bases[(x >> 33) as usize % 4]
            })
            .collect();
        let mut m = HashMap::new();
        m.insert("chr1".to_string(), seq.clone());
        m.insert("chr2".to_string(), seq);
        m
    }

    fn opts(remove_contained: bool, merge_overlapping: bool, min_overlap: f64) -> IntervalOptions {
        IntervalOptions { remove_contained, merge_overlapping, min_overlap }
    }

    fn status(r: &SequenceRecord) -> &str {
        r.validated.as_deref().unwrap_or("")
    }

    fn run(records: &mut Vec<SequenceRecord>, o: &IntervalOptions) {
        analyze_intervals(records, &genome(2000), o, &LogLevel::Summary);
    }

    #[test]
    fn required_overlap_reads_below_one_as_a_fraction_of_the_shorter_span() {
        assert_eq!(required_overlap(100, 200, 0.5), 50);
        assert_eq!(required_overlap(200, 100, 0.5), 50);
        // Rounds up, so a fraction never silently admits a shorter overlap.
        assert_eq!(required_overlap(101, 999, 0.5), 51);
        // A fraction can never fall to zero and let touching spans merge.
        assert_eq!(required_overlap(100, 100, 0.0), 1);
    }

    #[test]
    fn required_overlap_reads_one_and_above_as_base_pairs() {
        assert_eq!(required_overlap(100, 200, 1.0), 1);
        assert_eq!(required_overlap(100, 200, 40.0), 40);
    }

    /// Exact duplicates are containment's degenerate case: the earliest record
    /// in the file is the keeper and every later copy is contained in it.
    #[test]
    fn identical_spans_keep_the_earliest_record() {
        let mut r = vec![
            rec("f.fa", "chr1", 100, 200, '+', 0),
            rec("f.fa", "chr1", 100, 200, '+', 1),
            rec("f.fa", "chr1", 100, 200, '+', 2),
        ];
        run(&mut r, &opts(true, false, 1.0));
        assert_eq!(status(&r[0]), "valid");
        assert_eq!(status(&r[1]), "removed_contained");
        assert_eq!(status(&r[2]), "removed_contained");
    }

    #[test]
    fn containment_ignores_strand() {
        let mut r = vec![
            rec("f.fa", "chr1", 100, 200, '+', 0),
            rec("f.fa", "chr1", 120, 180, '-', 1),
        ];
        run(&mut r, &opts(true, false, 1.0));
        assert_eq!(status(&r[0]), "valid");
        assert_eq!(status(&r[1]), "removed_contained");
    }

    /// A record that swallows a shorter one can itself be swallowed by a longer
    /// one; the survivor must be the outermost span, not the middle link.
    #[test]
    fn a_nested_chain_collapses_to_its_outermost_span() {
        let mut r = vec![
            rec("f.fa", "chr1", 150, 160, '+', 0),
            rec("f.fa", "chr1", 120, 180, '+', 1),
            rec("f.fa", "chr1", 100, 200, '+', 2),
        ];
        run(&mut r, &opts(true, false, 1.0));
        assert_eq!(status(&r[0]), "removed_contained");
        assert_eq!(status(&r[1]), "removed_contained");
        assert_eq!(status(&r[2]), "valid");
    }

    /// Overlap is measured against the cluster's running span, so a staircase of
    /// tiled windows stays one cluster even though the ends never touch.
    /// The surviving record is the longest member of the cluster, not the first
    /// one in the file.
    #[test]
    fn the_longest_member_is_the_survivor() {
        let mut r = vec![
            rec("f.fa", "chr1", 1, 50, '+', 0),
            rec("f.fa", "chr1", 40, 200, '+', 1),
            rec("f.fa", "chr1", 190, 240, '+', 2),
        ];
        run(&mut r, &opts(false, true, 1.0));
        assert_eq!(status(&r[0]), "removed_merged");
        assert_eq!(status(&r[1]), "merged_overlapping");
        assert_eq!((r[1].start_1b, r[1].end_1b), (Some(1), Some(240)));
        assert_eq!(status(&r[2]), "removed_merged");
    }

    #[test]
    fn tiled_windows_chain_into_a_single_cluster() {
        let mut r = vec![
            rec("f.fa", "chr1", 1, 100, '+', 0),
            rec("f.fa", "chr1", 51, 150, '+', 1),
            rec("f.fa", "chr1", 101, 200, '+', 2),
        ];
        run(&mut r, &opts(false, true, 1.0));
        // Equal lengths, so the tie falls to the earliest record in the file.
        assert_eq!(status(&r[0]), "merged_overlapping");
        assert_eq!((r[0].start_1b, r[0].end_1b), (Some(1), Some(200)));
        assert_eq!(status(&r[1]), "removed_merged");
        assert_eq!(status(&r[2]), "removed_merged");
    }

    #[test]
    fn a_fractional_threshold_refuses_a_weak_join() {
        let spans = || {
            vec![
                rec("f.fa", "chr1", 1, 100, '+', 0),
                rec("f.fa", "chr1", 99, 200, '+', 1),
            ]
        };
        // 2 bp shared, which clears the 1 bp default...
        let mut loose = spans();
        run(&mut loose, &opts(false, true, 1.0));
        let merged = loose.iter().find(|x| status(x) == "merged_overlapping").unwrap();
        assert_eq!((merged.start_1b, merged.end_1b), (Some(1), Some(200)));

        // ...but not half of the shorter span.
        let mut strict = spans();
        run(&mut strict, &opts(false, true, 0.5));
        assert!(strict.iter().all(|x| status(x) == "valid"));
    }

    #[test]
    fn a_base_pair_threshold_cuts_exactly_where_it_says() {
        let spans = || {
            vec![
                rec("f.fa", "chr1", 1, 100, '+', 0),
                rec("f.fa", "chr1", 99, 200, '+', 1),
            ]
        };
        let mut at = spans();
        run(&mut at, &opts(false, true, 2.0));
        assert!(at.iter().any(|x| status(x) == "merged_overlapping"));

        let mut over = spans();
        run(&mut over, &opts(false, true, 3.0));
        assert!(over.iter().all(|x| status(x) == "valid"));
    }

    /// Grouping keys on the input file as well as the sequence id, so two
    /// families passed in one invocation never bleed into each other.
    #[test]
    fn separate_input_files_never_interact() {
        let mut r = vec![
            rec("a.fa", "chr1", 100, 200, '+', 0),
            rec("b.fa", "chr1", 100, 200, '+', 0),
        ];
        run(&mut r, &opts(true, true, 1.0));
        assert_eq!(status(&r[0]), "valid");
        assert_eq!(status(&r[1]), "valid");
    }

    #[test]
    fn separate_reference_sequences_never_interact() {
        let mut r = vec![
            rec("f.fa", "chr1", 100, 200, '+', 0),
            rec("f.fa", "chr2", 100, 200, '+', 1),
        ];
        run(&mut r, &opts(true, true, 1.0));
        assert_eq!(status(&r[0]), "valid");
        assert_eq!(status(&r[1]), "valid");
    }

    #[test]
    fn a_merged_record_carries_the_reference_sequence_for_its_union_span() {
        let g = genome(2000);
        let mut r = vec![
            rec("f.fa", "chr1", 11, 20, '+', 0),
            rec("f.fa", "chr1", 16, 25, '+', 1),
        ];
        analyze_intervals(&mut r, &g, &opts(false, true, 1.0), &LogLevel::Summary);
        assert_eq!((r[0].start_1b, r[0].end_1b), (Some(11), Some(25)));
        assert_eq!(r[0].sequence, g["chr1"][10..25].to_vec());
    }

    /// The longest member sets the orientation, and a minus-strand result is
    /// reverse complemented out of the reference.
    #[test]
    fn the_longest_member_sets_the_strand_and_the_sequence_follows_it() {
        let g = genome(2000);
        let mut r = vec![
            rec("f.fa", "chr1", 25, 34, '+', 0),
            rec("f.fa", "chr1", 11, 30, '-', 1),
        ];
        analyze_intervals(&mut r, &g, &opts(false, true, 1.0), &LogLevel::Summary);
        let merged = r.iter().find(|x| status(x) == "merged_overlapping").unwrap();
        assert_eq!((merged.start_1b, merged.end_1b, merged.orient), (Some(11), Some(34), Some('-')));
        assert_eq!(merged.sequence, reverse_complement(&g["chr1"][10..34]));
    }

    #[test]
    fn a_merged_record_drops_its_alignment_columns() {
        let mut r = vec![
            rec("f.stk", "chr1", 11, 20, '+', 0),
            rec("f.stk", "chr1", 16, 25, '+', 1),
        ];
        r[0].aligned_seq = Some(b"ACGT----ACGT".to_vec());
        r[1].aligned_seq = Some(b"ACGT----ACGT".to_vec());
        run(&mut r, &opts(false, true, 1.0));
        assert!(r[0].aligned_seq.is_none());
    }

    /// Dropping contained records first must not move the union span, since a
    /// contained record contributes no bases the cluster does not already hold.
    #[test]
    fn removing_contained_records_leaves_the_merged_span_unchanged() {
        let spans = || {
            vec![
                rec("f.fa", "chr1", 1, 100, '+', 0),
                rec("f.fa", "chr1", 20, 60, '+', 1),
                rec("f.fa", "chr1", 51, 150, '+', 2),
            ]
        };
        let mut merge_only = spans();
        run(&mut merge_only, &opts(false, true, 1.0));

        let mut both = spans();
        run(&mut both, &opts(true, true, 1.0));

        let span_of = |v: &Vec<SequenceRecord>| {
            let m = v.iter().find(|x| status(x) == "merged_overlapping").unwrap();
            (m.start_1b, m.end_1b)
        };
        assert_eq!(span_of(&merge_only), (Some(1), Some(150)));
        assert_eq!(span_of(&both), span_of(&merge_only));
    }

    /// A record still holding its originally-parsed, known-wrong coordinates
    /// must not be compared against records whose coordinates were resolved.
    #[test]
    fn records_without_resolved_coordinates_sit_the_pass_out() {
        let mut r = vec![
            rec("f.fa", "chr1", 100, 200, '+', 0),
            rec("f.fa", "chr1", 120, 180, '+', 1),
        ];
        r[1].validated = Some("invalid".to_string());
        run(&mut r, &opts(true, true, 1.0));
        assert_eq!(status(&r[0]), "valid");
        assert_eq!(status(&r[1]), "invalid");
    }

    #[test]
    fn detection_alone_changes_nothing() {
        let mut r = vec![
            rec("f.fa", "chr1", 1, 100, '+', 0),
            rec("f.fa", "chr1", 20, 60, '+', 1),
            rec("f.fa", "chr1", 51, 150, '+', 2),
        ];
        run(&mut r, &opts(false, false, 1.0));
        assert!(r.iter().all(|x| status(x) == "valid"));
        assert_eq!((r[0].start_1b, r[0].end_1b), (Some(1), Some(100)));
    }
}
