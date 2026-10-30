# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

## 0.6.9 - 2026-10-30
### Added
- stk convert --to comsa: write the (cleaned) Stockholm output as a CoMSA
  archive (.smsac) using comsa-rs 0.1.2. The archive decompresses to exactly
  what --to stk writes and is readable by CoMSA 1.1
- stk reads CoMSA archives as input: every subcommand that takes Stockholm
  (lint, extract, edit, convert, map-consensus) detects an archive by content
  and decompresses it one family at a time. `stk convert archive.smsac`
  decompresses to Stockholm by default
- Dfam Stockholm JSON, a JSON form of Dfam seed alignment files for programs
  that ingest them (spec: docs/dfam-stockholm-json.md; JSON Schema alongside).
  `stk convert --to json` writes it, and `stk convert` reads it back,
  recognised by content. Stockholm stays the interchange format: the other
  subcommands refuse JSON input and point to `stk convert`
- stk lint --format json: the findings as one JSON document, for programs
  that run lint (docs/stk-lint-json.md)
### Fixed
- stk lint reports a file with no Stockholm record at all as the ERROR
  `no_records`.  A FASTA, EMBL, plain-text or empty file used to lint clean
- stk exits with status 2 on any failure (unreadable or unsupported input, I/O
  error), as documented.  It used to exit 1, which for `stk lint` was
  indistinguishable from "found errors"
### Changed
- stk binaries now link comsa-rs (GPL-3.0-only) and are distributed under
  the GPL-3.0; the dfam-curator sources remain CC0

## 0.6.8 - 2026-10-28
### Fixed
- --fixed-cons <FASTA>: align the instances to a final consensus once and
  write it, unrefined, as #=GC RF (no bootstrap, refinement, extension or
  block repair); a single instance is enough
- stk lint: an RF line made only of 'X'/'x' and gaps is a legacy occupancy
  mask, not a called consensus, and is no longer compared against the
  consensus. This removes most of the warnings on unedited DR families
  exported from Dfam 4.0
- stk lint: AU names with middle initials ("Robert M. Hubley",
  "Arian F.A. Smit") are accepted; a '.' is flagged only in the first or
  last word, or in a middle word that is not initials

## 0.6.7 - 2026-10-20
### Changed
- Bumped to rmblast-lib 3.0.10, dfam-lib 0.2.4 and RepeatAfterMe
  (ram-core) 1.0.0

## 0.6.6 - 2026-10-19
### Fixed
- cons-core/te-composer: minus-strand instances were laid into the MSA
  backwards, so inputs in genome orientation (e.g. LTR_retriever records)
  produced N-riddled consensi. They are now placed in reference orientation
### Changed
- Bumped dfam-lib

## 0.6.5 - 2026-09-02
### Added
- te-composer: --min option for the Dfam consensus caller

## 0.6.4 - 2026-08-28
### Added
- te-composer documentation

## 0.6.3 - 2026-08-27
### Fixed
- Builds on non-x86 targets (aarch64-apple-darwin): parasail is x86-only,
  so other targets fall back to the reference aligner

## 0.6.2 - 2026-08-27
### Added
- te-composer: Refiner/autocons replacement for building consensus
  sequences
### Changed
- Moved to dfam-lib 0.2.0 spans and RepeatAfterMe 0.2.0; linup output is
  unchanged

## 0.6.1 - 2026-07-16
### Added
- Initial Repbase-to-Dfam classification translation table
- Conversion-time checks for species and name problems

## 0.6.0 - 2026-07-15
### Added
- stk lint: legacy Smitten identifier format check (seqid_legacy_format),
  RepeatModeler automated name checks (id_automated_name, automated_name),
  and cross-record summaries for AC update records and missing ORCIDs
- stk writing commands (edit, extract, convert --to stk) normalize output to
  Dfam conventions by default (gap characters to '.', 7-digit accessions
  widened to 9); --no-clean emits records verbatim
- stk repbase-import (internal): combines a Repbase IG MSA and IG family
  record into one Dfam Stockholm record
- Readers for the Repbase IG MSA and IG family-record formats, with format
  auto-detection
- examples/ directory
### Changed
- discoord emits full V2 identifiers (assembly stamped in) for every
  resolved record, not only remapped ones
### Fixed
- linup: the reference row now numbers from 1 instead of 0

## 0.5.0 - 2026-06-17
### Added
- stk lint validates the AU convention
  `[ORCID:####-####-####-####] firstname lastname[; ...]`
- stk lint validates the RD (DOI) and RM fields against DOI.org and PubMed
- tsd: target site duplication finder
### Changed
- Format conversion moved from linup to `stk convert`, with support for
  multi-record formats such as Stockholm
- Release archives unpack to `dfam-curator-<version>/bin/<tools>`
### Removed
- Windows build target

## 0.4.3 - 2026-06-10
- First tagged release: linup, stk and discoord; native release builds for
  Linux x86_64 and macOS Apple Silicon
