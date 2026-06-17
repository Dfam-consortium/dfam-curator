# dfam-curator

A toolkit for validating and curating transposable element family seed alignments
prior to submission to [Dfam](https://dfam.org) database. Provides command-line tools for
viewing, validating, editing, and repairing Stockholm-format alignment files,
along with a library of shared alignment and consensus-calling routines.

---
 

## Tools

- **[linup](#linup)** — MSA viewer and format converter
- **[stk](#stk)** — Stockholm file metadata tool (lint · edit · extract)
- **[discoord](#discoord)** — Sequence coordinate validator and repairer
- **[update-cache](#update-cache)** — Populate the stk-lint validation cache

---

## Installation

### Pre-built binaries

Download the latest release for your platform from the
[GitHub releases page](https://github.com/Dfam-consortium/dfam-curator/releases).
Extract the archive and place the binaries (`linup`, `stk`, `discoord`, `update-cache`)
somewhere on your `PATH`.

### Build from source

The `main` branch always reflects the current release, so building from source
is equivalent to installing a release tarball. Requires a recent stable
[Rust toolchain](https://rustup.rs).

```sh
git clone https://github.com/Dfam-consortium/dfam-curator.git
cd dfam-curator
cargo build --release
```

Compiled binaries are placed in `target/release/`.

---

## Quick Start

Initialize the dfam-curator data cache.  This stores reference data from the
NCBI Taxonomy, and the Dfam database used to validate new submissions.

```sh
./target/release/update-cache
```

Validate metadata in a Stockholm file of TE families that are ready for submission to Dfam:

```sh
./target/release/stk lint myfile.stk
```

A more rigorous validation involves checking that the sequence identifiers and sequences
in the Stockholm MSA can be mapped back to an assembly file.  To perform this validation:

```sh
./target/release/stk lint --genome myassembly.fa myfile.stk
```

Basic STK editing of metadata is provided using the `stk edit` command.  For more information
use '-h' with any of these subcommands.

---

## linup

MSA viewer and format converter. Reads a multiple sequence alignment in
Stockholm, FASTA/A2M, or crossmatch `.align` format (or rmblastn tabular
output with `--blast-tab`) and writes it in a requested format with optional
trimming, slicing, or reverse-complementing. A Rust port of
`RepeatModeler/util/Linup`.

```
USAGE:
    linup [OPTIONS] <INPUT>

ARGUMENTS:
    <INPUT>    Input alignment file

OPTIONS:
    --blast-tab              Treat input as rmblastn tabular output (requires --ref-seq)
    --ref-seq <FILE>         Reference FASTA used as the BLAST subject
    --format <FORMAT>        Output format [default: linup]
                               linup       Perl Linup-compatible pretty-print blocks
                               stockholm   Stockholm 1.0 with #=GC RF consensus line
                               msa         Aligned FASTA / A2M (gaps preserved)
                               fasta       Unaligned FASTA (gaps stripped)
                               consensus   Consensus sequence only
                               stats       Per-sequence Kimura divergence statistics
    -i, --include-ref        Include the reference row when calling the consensus
    --trim-left <N>          Trim N ungapped reference bp from the left
    --trim-right <N>         Trim N ungapped reference bp from the right
    --trim-ambig             Trim ambiguous (non-ACGT) bases from both ends of the consensus
    --sub-align <START-END>  Slice to a 1-based, fully-closed consensus coordinate range
    --min-len <N>            Minimum bases retained after --sub-align (drops shorter sequences)
    --revcomp                Reverse-complement the entire alignment before output
    --name <STRING>          Override the family name / ID in output
    --include-gaps           With --format consensus: retain gap characters
    --select <SELECT>        Select one record from a multi-record Stockholm file
                               Numeric value → 1-based record number
                               Non-numeric   → exact #=GF ID match
```

`--trim-left/--trim-right`, `--trim-ambig`, `--sub-align`, and `--revcomp`
are mutually exclusive. `--min-len` is only valid with `--sub-align`.

---

## stk

Dfam Stockholm file toolkit. Operates on one or more `.stk` files and provides
three subcommands.

```
USAGE:
    stk <SUBCOMMAND>
```

### stk lint

Validate Stockholm files and report structural and semantic diagnostics.

```
USAGE:
    stk lint [OPTIONS] <INPUT>...

ARGUMENTS:
    <INPUT>...    One or more Stockholm files to check

OPTIONS:
    --cache-dir <PATH>          Override the cache directory
                                  (default: $STK_CACHE_DIR, then ~/.cache/stk)
    --min-severity <LEVEL>      Minimum severity to report: error | warn | info  [default: info]
    --no-cache-warn             Suppress notices about missing tier-2 cache files
    --genome <FILE>             Reference genome (FASTA or .2bit) for coordinate validation
    --no-network                Skip live validation of PubMed IDs (RM) and DOIs (RD)
```

Exit status: 0 = clean, 1 = at least one ERROR, 2 = I/O failure.

Checks include: Stockholm structure, required annotation fields, record terminator,
duplicate IDs (with AC-aware severity), RF/consensus agreement, taxonomy names,
Dfam classification strings, live PubMed/DOI resolution (disable with `--no-network`),
and (with `--genome`) sequence coordinate validity.
Validation is split into tier-1 (always available) and tier-2 (requires a populated
cache; see [update-cache](#update-cache)).

### stk extract

Extract a single record from a multi-record Stockholm file.

```
USAGE:
    stk extract --select <SELECT> [OPTIONS] <INPUT>...

ARGUMENTS:
    <INPUT>...    One or more Stockholm files to search

OPTIONS:
    --select <SELECT>    Record to extract (required)
                           Numeric value → 1-based record number
                           Non-numeric   → exact #=GF ID match
    -o, --output <FILE>  Write output to FILE instead of stdout
```

```sh
stk extract --select MyFam families.stk
stk extract --select 3 families.stk
stk extract --select MyFam -o MyFam.stk families.stk
```

### stk edit

Edit `#=GF` annotation fields across records. Operations are applied in a
fixed order: `--delete`, then `--set`, then `--append`, then `--sub`.

```
USAGE:
    stk edit [OPTIONS] <INPUT>...

ARGUMENTS:
    <INPUT>...    One or more Stockholm files to edit

OPTIONS:
    --set <TAG> <VALUE>      Set (or replace) a GF field; repeatable
    --delete <TAG>           Remove all occurrences of a GF tag; repeatable
    --append <TAG> <VALUE>   Append a new GF field (for multi-valued tags); repeatable
    --sub <TAG> <EXPR>       Apply a regex substitution to all values of a GF tag
                               Format: /PATTERN/REPLACEMENT/[g]
                               First char is the delimiter (any char works)
                               Append /g to replace all matches; omit for first-match only
                               Capture groups use $1, $2, … in the replacement
                               For multi-valued fields (OC, CC) applied per line
    --update-consensus       Recompute the #=GC RF consensus from the aligned sequences
                               Replaces any existing RF line
    --select <SELECT>        Only edit matching records; others pass through unchanged
                               Numeric value → 1-based record number
                               Non-numeric   → exact #=GF ID match
    -o, --output <FILE>      Write output to FILE instead of stdout
```

Operations are applied in a fixed sequence: `--delete`, then `--set`, then `--append`, then `--sub`.

```sh
stk edit --set AU "Barbara McClintock" families.stk
stk edit --delete SE --set DE "Updated description" families.stk
stk edit --select MyFam --append OC "Mus musculus" families.stk
stk edit --set AU "Barbara McClintock" -o fixed.stk families.stk
stk edit --sub ID "/^(.*)-$/$1/" families.stk
stk edit --sub DE "/foo/bar/g" families.stk
stk edit --update-consensus families.stk
stk edit --select MyFam --update-consensus families.stk
```

---

## discoord

Validates and repairs sequence coordinates in FASTA or Stockholm files by
comparing extracted subsequences against a reference genome. Detects
half-open interval errors, small coordinate shifts, strand inversions, and
sequences that have migrated to a different genomic location.

```
USAGE:
    discoord [OPTIONS] <INPUT>...

ARGUMENTS:
    <INPUT>...    Input files (FASTA, Stockholm, or tab/comma-delimited)

OPTIONS:
    -d, --reference-dir <DIR>    Directory of per-assembly reference files
    -r, --reference-default <FILE>
                                 Default reference genome (FASTA or .2bit)
    -m, --map-sequences          Enable sequence mapping for unresolvable coordinates
                                   Uses Aho-Corasick by default (single-pass dual-strand)
    --boyer-moore                Use Boyer-Moore instead of Aho-Corasick with -m
    -l, --log-level <LEVEL>      Output verbosity: summary | per-record | detailed
                                   [default: summary]
    -x, --threads <N>            Number of threads for parallel processing
    -o, --output-dir <DIR>       Write corrected output files to DIR
                                   Use "." to write alongside each input file
    -u, --remove-duplicates      When mapping, remove a sequence if every genomic hit is
                                   already occupied by an earlier sequence in the file.
                                   Kept sequences are assigned to the first unoccupied hit.
                                   Removed sequences are reported as "removed_remapped_duplicate"
```

### Coordinate validation

Without `-m`, discoord checks each sequence's stored coordinates against the
reference genome and attempts common repairs:

| Result label | Meaning |
|---|---|
| `valid` | Coordinates are correct as given |
| `fixed_half_open` | Converted from half-open to fully-closed interval |
| `fixed_shift_±N` | Corrected by a small positional shift (±1–3 bp) |
| `fixed_revcomp` | Strand was inverted |
| `invalid` | Could not be resolved |

### Sequence mapping (`-m`)

When `-m` is set, sequences that remain unresolved after coordinate checking
are searched against the reference genome. The best match is chosen by:

1. Non-redundancy — avoids a location already occupied by another sequence in the set
2. Proximity — prefers the same chromosome and start coordinate closest to the original
3. Uniqueness — reports `fixed_remapped_unique` or `fixed_remapped_ambig`

With `-u` (`--remove-duplicates`), sequences whose every hit is already occupied are
dropped and reported as `removed_remapped_duplicate` rather than retained at an occupied
position.

Aho-Corasick (default) searches both strands in a single genome pass;
Boyer-Moore (`--boyer-moore`) searches each strand separately.

### Delimited input format

Tab- or comma-delimited files must have these columns (with a header row):

| Column | Description |
|---|---|
| `assembly_id` | Assembly or reference identifier (may be blank) |
| `sequence_id` | Sequence/chromosome name in the reference |
| `start` | 1-based start coordinate |
| `end` | 1-based fully-closed end coordinate |
| `sequence` | Sequence data |

---

## update-cache

Downloads and populates the tier-2 validation data used by `stk lint`.

```
USAGE:
    update-cache [OPTIONS] [SUBCOMMAND]

SUBCOMMANDS:
    all             Download all data, then show cache status  [default]
    taxonomy        Download NCBI taxonomy only
    classifications Download Dfam TP classifications only
    names           Download Dfam family names only
    info            Show cache status

OPTIONS:
    --cache-dir <PATH>    Override the cache directory
                            (default: $STK_CACHE_DIR, then ~/.cache/stk)
    --force               Re-download files even if they already exist
```

All four cache files are fetched automatically — no manual downloads required:

| File | Source |
|---|---|
| `taxonomy.tsv` | NCBI taxonomy scientific names (from taxdmp.zip) |
| `taxonomy-common.tsv` | Common name / synonym → scientific name mapping |
| `classification.tsv` | Dfam TP classification strings (from dfam.org) |
| `dfam-names.txt` | Dfam family names (from dfam.org) |

Requires `curl` and `unzip` to be available on `PATH`.

---

## Documentation

- **[Dfam_Seeds.md](Dfam_Seeds.md)** — Detailed guide to the Dfam seed alignment format:
  what a seed record is, the Stockholm 1.0 format, required and optional metadata fields,
  sequence identifier conventions (Smitten format), and annotated examples ranging from
  a minimal valid record to a full submission.

---

## Library crates

| Crate | Description |
|---|---|
| `dfam-stk-io` | Stockholm 1.0 streaming parser, `StkRecord` and `SeqRow` types, Smitten identifier integration |
| `dfam-coord` | Coordinate validation engine used by `discoord`: genome loading (FASTA / .2bit), range checking, offset correction, sequence mapping, and Stockholm/FASTA/delimited I/O |

---

## License

This project is released under the
[CC0 1.0 Universal Public Domain Dedication](LICENSE).
To the extent possible under law, the authors have waived all copyright and
related rights to this software.
