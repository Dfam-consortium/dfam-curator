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

#[derive(Clone, Debug)]
pub struct SequenceRecord {
    pub input_file: String,
    pub metadata_idx: usize,
    pub order: usize,
    pub original_id: Option<String>,
    pub assembly_id: Option<String>,
    pub sequence_id: String,
    pub start: Option<u64>,
    pub end: Option<u64>,
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
            start: row.seq_start,
            end: row.seq_end,
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
            self.start.unwrap_or(0),
            self.end.unwrap_or(0),
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
    remapped_assembly: Option<&str>,
    remove_duplicates: bool,
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
            aho_corasick_search_with_validation(&mut results, genome_map, debug_mode, remapped_assembly, remove_duplicates);
        } else {
            boyer_moore_search_with_validation(&mut results, genome_map, debug_mode, remapped_assembly, remove_duplicates);
        }
    }

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

            if record.start.is_some() && record.end.is_some() && record.orient.is_some() {
                writeln!(writer, "{}:{}-{}_{} {}", v2_id, record.start.unwrap(), record.end.unwrap(),
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
        if record.start.is_some() && record.end.is_some() && record.orient.is_some() {
            writeln!(writer, ">{}:{}-{}_{} {}", v2_id, record.start.unwrap(), record.end.unwrap(),
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
        let start = record.start.map(|v| v.to_string()).unwrap_or_default();
        let end = record.end.map(|v| v.to_string()).unwrap_or_default();
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

        if record.start.is_none() {
            if genome_sequence == &record.sequence {
                record.validated = Some("valid".to_string());
                *fix_counts.entry(record.validated.clone().unwrap()).or_insert(0) += 1;
            }
            continue;
        }

        let range_length = match (record.end, record.start) {
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
            record.start = record.start.map(|start| start + 1);
        }

        let start = record.start.unwrap() as usize - 1;
        let end = record.end.unwrap() as usize;
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
                        record.start = Some((shifted_start + 1) as u64);
                        record.end = Some(shifted_end as u64);
                        record.orient = if record.orient == Some('-') { Some('+') } else { Some('+') };
                        located = true;
                        break;
                    }

                    if &genome_sequence[shifted_start..shifted_end] == &rev_complement {
                        validation_str.push_str(&format!("{}{}{}",
                            if record.orient == Some('+') { "_orient" } else { "" },
                            if *shift >= 0 { "_plus" } else { "_minus" },
                            shift.abs()));
                        record.start = Some((shifted_start + 1) as u64);
                        record.end = Some(shifted_end as u64);
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
                        let dist_a = record.start.map_or(usize::MAX, |start| (a.0 as isize - start as isize).unsigned_abs());
                        let dist_b = record.start.map_or(usize::MAX, |start| (b.0 as isize - start as isize).unsigned_abs());
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
        .filter_map(|r| match (r.start, r.end, r.orient) {
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
                record.start = Some(best.0 as u64 + 1);
                record.end = Some((best.0 + pat_len) as u64);
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
        sort_hits_by_proximity(&mut hits, &record.sequence_id, record.start);
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

            sort_hits_by_proximity(&mut found_positions, &original_sequence_id, record.start);
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
        .filter_map(|r| match (r.start, r.end, r.orient) {
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
                record.start = Some(best.0 as u64 + 1);
                record.end = Some((best.0 + pat_len) as u64);
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

pub fn output_results(records: &[SequenceRecord], format: LogLevel, label: String) {
    match format {
        LogLevel::Summary | LogLevel::PerRecord => {
            let total_records = records.len();
            let mut fix_counts = HashMap::new();
            let mut fixed_count = 0;
            for record in records.iter() {
                let v = record.validated.as_deref();
                if v.is_some() && v != Some("valid") && v != Some("invalid")
                    && v != Some("removed_remapped_duplicate")
                {
                    fixed_count += 1;
                    *fix_counts.entry(record.validated.clone().unwrap()).or_insert(0) += 1;
                }
            }
            let valid_count = records.iter().filter(|r| r.validated.as_deref() == Some("valid")).count();
            let invalid_count = records.iter().filter(|r| r.validated.as_deref() == Some("invalid")).count();
            let removed_dup_count = records.iter().filter(|r| r.validated.as_deref() == Some("removed_remapped_duplicate")).count();

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
