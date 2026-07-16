# Dfam Seed Alignments

## What Is a Dfam Seed Alignment?

A Dfam seed alignment is a curated multiple sequence alignment (MSA) of genomic instances
(copies) of a repetitive element family, together with the family's metadata.  It is the
primary submission artefact for a Dfam entry: everything from the consensus sequence to the
taxonomic scope of the family is derived from it.

Each family in the database is represented by one **seed alignment record**.  The record
contains two things:

1. **Metadata** — who described the family, what it is classified as, which organisms it
   occurs in, literature references, and so on.
2. **Multiple sequence alignment** — the actual genomic copies that were used to build or
   verify the consensus, aligned column-by-column so that homologous positions are in register.

The format used to store both of these together is **Stockholm 1.0**, a plain-text alignment
format that was developed for the Pfam and Rfam databases and later adopted by Dfam.

---

## Stockholm 1.0 Format Overview

A Stockholm file is a sequence of one or more records, each delimited by a header line and a
terminator line:

```
# STOCKHOLM 1.0
...annotation lines...
...sequence rows...
//
```

There are three kinds of annotation lines:

| Prefix  | Meaning                                     |
|---------|---------------------------------------------|
| `#=GF`  | Per-**f**ile (record-level) metadata field  |
| `#=GC`  | Per-**c**olumn annotation (one value per alignment column); Dfam uses `RF` and `MM` |
| `#=GS`  | Per-**s**equence metadata (not used in Dfam) |
| `#=GR`  | Per-**r**esidue annotation (not used in Dfam) |

Stockholm has **no free-text comment syntax**.  Every line beginning with `#` must use one of
the structured prefixes above — plain `# comment` lines are not valid inside a record.

### Sequence rows

Each aligned sequence is written as a name, whitespace, and the aligned sequence string:

```
GCA_000001405.15:chr1:100-200_+    ACGT..ACGT
```

- All sequence rows in a record must have the same aligned length.
- Dfam has adopted the **Smitten** identifier format (see
  [MSA Sequence Identifiers](#msa-sequence-identifiers) below).

#### Gap characters

Stockholm permits four gap characters — `-` (dash), `.` (period), `_` (underscore) and
`~` (tilde) — and `stk` reads all four.  **Dfam has standardized on `.` (period).**  A file
using any of the others is accepted but `stk lint` reports a warning (`seq_nonstandard_gap`),
one per sequence row, naming the character it found.

Whitespace is never a gap.

```
GCA_000001405.15:chr1:100-200_+    ACGT..ACGT      # Dfam convention
GCA_000001405.15:chr1:100-200_+    ACGT--ACGT      # read, but warned
```

#### One line per sequence — block format is not accepted

Stockholm allows a long alignment to be split into **blocks**: the alignment is broken into
column ranges, each sequence appears once per block, and the blocks are separated by blank
lines.  At this time Dfam does not permit the use of block format and requires that each
seed sequence appears on a single line.

### The `#=GC RF` line

The RF (reference/consensus) line is a per-column annotation that encodes the consensus
sequence of the family (derived from the MSA).  It uses IUPAC/IUB ambiguity codes for 
columns with mixed bases and `.` for gap-only columns.  In older versions of Dfam 
Stockholm files the characters `X` or `x` were used to indicate consensus occupancy 
columns for use with the HMMER hmmbuild tool.

By default the RF line is expected to be the consensus called from the alignment.  A curator
who has hand-built the consensus instead should say so with the [`CT` field](#ct--consensus-type).

### The `#=GC MM` line

The MM (model mask) line is an optional per-column annotation, defined by HMMER, that marks
alignment columns lying within a **masked range**:

| Character | Meaning                                              |
|-----------|------------------------------------------------------|
| `m`       | The column lies within a masked range                |
| `.`       | The column is not masked                             |

For a match state built from a masked column, `hmmbuild` emits background frequencies rather
than the frequencies observed in the alignment.  In other words, the column still contributes
to the length and architecture of the model, but carries no information content — it neither
adds to nor subtracts from a match score.  This is the mechanism to use when a region of the
family should be kept in the model for spacing reasons but should not itself drive matches —
for example a simple-repeat or otherwise low-complexity stretch that would generate spurious
hits.

```
#=GC RF    CAGACTTGATGCACGACGTAAACGTGACT
#=GC MM    ..........mmmmmmmm...........
```

Rules:

- Only `m` and `.` are valid characters (lower-case `m`; `stk lint` reports
  `mm_invalid_chars` otherwise).
- The line must be exactly as wide as the `#=GC RF` line and the sequence rows
  (`mm_length_mismatch`).
- The field is optional.  A record with no `MM` line is treated as fully unmasked.

Note that `stk` writes `#=GC RF` immediately before the sequence rows and any other `#=GC`
annotation, including `MM`, immediately after them.  Both placements are valid Stockholm and
are read correctly by HMMER.

Dfam imports only `RF` and `MM`.  Other Stockholm per-column annotations (`SS_cons`,
`PP_cons`, and so on) are passed through the file but ignored, and `stk lint` notes them with
an INFO (`unknown_gc_tag`).

---

## Minimal Example

The smallest valid Dfam seed record contains the five required metadata fields, an RF line,
and at least one sequence row:

```
# STOCKHOLM 1.0
#=GF DE    Short interspersed nuclear element, type 1
#=GF AU    Smith J
#=GF TP    Interspersed_Repeat;SINE/ID
#=GF OC    Mammalia
#=GF SQ    1
#=GC RF    CAGACTTGATGCACGACGTAAACGTGACT
seq1       CAGACTTGATGCACGACGTAAACGTGACT
//
```

This record has no accession (it has not yet been accepted into Dfam), no literature
references, and only two sequences.  It is the minimum that dfam-curator validation
tool `stk lint` will accept without errors.

---

## More Complete Example

A more complex submission will normally include an optional name (ID field), longer 
description (CC fields), publications (RN/RM fields), cross-database references (DR fields), 
and multiple well-annotated sequences.  Here's an example of a more complete record, 
with all of these:

```
# STOCKHOLM 1.0
#=GF DE    MER3822 internal sequence
#=GF AU    Robert Smith
#=GF TP    Interspersed_Repeat;LTR/ERV1
#=GF OC    Bos taurus
#=GF OC    Homo sapiens
#=GF SQ    3
#=GF RN    [1]
#=GF RM    12490638
#=GF DR    Repbase; MER4-int
#=GF CC    Internal sequence of the MER4 family of LTR retrotransposons.  This
#=GF CC    family can be described fully here using any number of CC lines.
#=GC RF                                    CAGACTTGATGCACGACGTAAACGTGACT.GGATCGCCCTAGATCGCGATCGATCGATCG
GCA_000001405.15:chr1:1234567-1234627_+    CAGACTTGATGCACGACGTAAACGTGACTCGGATCGCCCTAGATCGCGATCGATCGATCG
GCA_000001405.15:chr5:9876543-9876603_-    CAGACTTGATGCACGACGTAAACGTGACT.GGATCGCCCTAGATCGCGATCGATCGATCG
GCA_000001405.15:chr7:5555000-5555060_+    CAGACTTGATGCCCGACGTAAACGTGACT.GGATCGCCCTAGATCGCGATCGATCGATCG
//
```

Notes on this example:

- Upon acceptance to Dfam this record will be updated with a `AC` field providing the 
  accession and version suffix.
- Two `OC` lines together say: "this family is found in Bos taurus, and in Homo sapiens". 
  This encoding supports both clade level assignment as well as evidence for 
  horizontal transfer.
- The publication information currently supports PMID identifiers only.  To specify a
  publication provide `RN` (the reference-order tag) **first**, followed immediately by
  `RM` (the PubMed ID).  `RN` must precede all other fields in its block — `stk lint`
  will report an error if that ordering is violated, and `stk edit --set` will reorder
  fields automatically.  Additional publications may be encoded with `RN [2]` and so on.
  Dfam will soon be replacing the PubMed specification with a DOI one.  Once a family
  has been accepted to Dfam additional fields will be automatically populated with
  publication title, authors, and journal (`RT`, `RA`, `RL`).
- The RF line contains the consensus for the three sequences, and is called using specific
  consensus caller used by Dfam ( found in the dfam-curator toolkit ).
- Sequence names follow the Smitten format (assembly accession, sequence identifier, coordinates,
  strand).

---

## Required Fields

The following `#=GF` fields and the `#=GC RF` line are **required**.  `stk lint` will report
an `ERROR` for any record that is missing one.

| Field      | Meaning                                                  |
|------------|----------------------------------------------------------|
| `DE`       | Short description of the family (max 80 characters, one line) |
| `AU`       | Author(s) who created or curated this record             |
| `TP`       | Classification path in the RepeatMasker hierarchy        |
| `OC`       | Taxonomic scope — one or more scientific taxon names     |
| `SQ`       | Number of sequence rows (must equal the actual count)    |
| `#=GC RF`  | Consensus / reference annotation (one character per alignment column) |

---

## Field Reference

### `AC` — Accession

The Dfam accession number.  Format: `DF` or `DR` followed by 9 digits, with an optional
version suffix (e.g. `DF000000001`, `DR000123456.2`).

> **`AC` is assigned by Dfam, not by the submitter.**  It appears only on families that have
> already been released in the database.  **Do not put an `AC` field on a new submission** —
> and never invent one.  It is not a required field, and a record without one is perfectly
> valid.

- Accessions beginning with `DF` are curated families.
- Accessions beginning with `DR` or "Dfam Raw" families are uncurated families 
  that have been generated by de novo TE identification tools.
- Records submitted for the first time will not have an `AC` field yet.  Once assigned, `AC`
  becomes the stable, authoritative identifier for the family.
- The presence of an `AC` is what tells Dfam "this record **replaces** the already-released
  family with this accession".  An invented accession that happens to be well-formed (e.g.
  `DR0000001`) is therefore read as an update to a real, unrelated family rather than as a
  new submission.
- Because a plausible-looking invented accession cannot be told apart from a real one by format
  alone, `stk lint` reports an INFO line (`ac_update_records`) counting how many records in the
  file carry an `AC`, and naming them.  If you did not intend those records to update existing
  Dfam families, remove the `AC` field.

### `ID` — Identifier (name)

A short human-readable name for the family (max 45 characters).

- Providing a name is not required, but it is highly recommended for curated families or
  families that appear specifically in publications.  For example, a 2,000 family library
  produced by a de novo tool does not need family-level names, but a deeply curated library
  of 200 families that are described in a publication should have names.
- Must be **unique** within file and must not already exist in Dfam (unless the record
  carries an `AC` field marking it as an update to the existing Dfam family).
- Names generated by automated pipelines such as RepeatModeler (e.g.
  `rnd-3_family-101`) are almost always non-unique and will typically fail the uniqueness
  check.  Rename families to something descriptive before submission.
- Names should be limited to alphanumeric and common delmiters 
  such as `-`, `_`, `.`, and `:`.  In addition, names must not be purely numeric.

### `DE` — Description

A concise one-line description of the family (max 80 characters, one line only).

Examples:

```
#=GF DE    Short interspersed nuclear element, type 1
#=GF DE    Long terminal repeat of the HERV-K endogenous retrovirus

or for de novo generated families a one-line description might be:

#=GF DE    Repetitive element family identified in Gadus morhua
```

### `AU` — Author

Who created or curated this seed alignment.  Use **`First Last`** for each author,
with multiple authors separated by **semicolons**.  Commas are not allowed as
separators because they appear in some names.

```
#=GF AU    Barbara McClintock
#=GF AU    Pita Enriquez-Lopez; Weidong Bao
```

The older `Last Initial` convention (`Smith J`) is **no longer accepted** — spell out the
given name.  Multiple `AU` lines are allowed, and each line may itself carry several
semicolon-separated authors.

Curators who have an ORCID may prefix their name with it:

```
#=GF AU    ORCID:0000-0001-2345-6789 Barbara McClintock
#=GF AU    ORCID:0000-0002-1825-0097 Josiah Carberry; Weidong Bao
```

The ORCID is written as `ORCID:` immediately followed by the bare
`xxxx-xxxx-xxxx-xxxx` identifier (the final group may end in `X`, the ORCID check digit),
then a space, then the author's name.  The prefix applies only to the author it precedes,
so a mix of credited and uncredited authors on one line is fine.  An ORCID must not be
repeated within a single `AU` field.

ORCIDs are optional but encouraged: they are what lets Dfam credit a curator unambiguously
rather than guessing among everyone who shares their name.  `stk lint` reports an INFO line
(`orcid_missing`) counting how many records name at least one author without one.

Formats that `stk lint` flags as **errors**:
- Single-letter first name: `B McClintock`
- Old `Last Initial` style: `McClintock B`, `Smith J`
- The same ORCID used twice in one `AU` field

Formats that `stk lint` flags as **warnings**:
- Any period in a name word: `B. McClintock`, `Barbara J. McClintock`, `A.F.A. Smit`
- Comma-separated lists: `Barbara McClintock, Roy Britten`
- Missing space (single token): `McClintock`
- A colon outside an `ORCID:` prefix
- An `ORCID:` prefix with no name after it, or a malformed identifier

A single-letter middle initial *without* a period is accepted: `Barbara B McClintock`.

### `SE` — Source

Where the consensus or the seed sequences were originally obtained from (max 80 characters).

```
#=GF SE    RepeatModeler 2.0
#=GF SE    Repbase:AluSx
```

### `TP` — Type / Classification

The classification string describing where this family sits in the **RepeatMasker hierarchy**.
It is a semicolon-delimited path from the root class down to the family type.

```
#=GF TP    Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;Retrotransposon;Long_Terminal_Repeat_Element;Gypsy-ERV;Retroviridae;Orthoretrovirinae;ERV1
#=GF TP    Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;7SL-RNA_Promoter;No-core;L1-dependent;Alu
#=GF TP    Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Tc1-Mariner;Tc1
#=GF TP    Tandem_Repeat;Satellite;Centromeric
```
or shorthand aliases for each may be used:

```
#=GF TP    LTR/ERV1
#=GF TP    SINE/Alu
#=GF TP    DNA/TcMar-Tc1
#=GF TP    Satellite/centromeric
```

When using `stk lint` with the `--classification` option, the value is validated against the
official Dfam classification list found [at the Dfam website](https://www.dfam.org/classification/tree).
See the TSV file for the full list of valid classifications including the shorthand "type/subtype"
format for each.


### `OC` — Organism / Taxonomic Scope

One or more lines identifying the taxonomic scope of the family.  Each `OC` value must be a
**scientific name** (species or higher-rank taxon name) recognised by NCBI Taxonomy.

```
#=GF OC    Homo sapiens
#=GF OC    Mammalia
#=GF OC    Eukaryota
```

Multiple `OC` lines are valid in a single family record and may be used to identify
specify the scope of horizontally transfer.

When using `stk lint` with the `--taxonomy` option, each value is validated against a list of
NCBI scientific names, common names, and aliases.  

### `SQ` — Sequence count

The number of sequence rows in the record.  Must be a non-negative integer and must exactly
match the number of sequence rows that follow.

```
#=GF SQ    47
```

### `TD` — Target site duplication

The target site duplication (TSD) consensus sequence, written using IUB ambiguity codes.

```
#=GF TD    TTAAAA
#=GF TD    TA
```

When only the length of the typical TSD is known the length may be encoded as a string
of Ns.  For example a family with 5bp TSDs would be encoded as:

```
#=GF TD  NNNNN
```

### `CT` — Consensus type

Describes how the `#=GC RF` consensus was produced.  The field is optional; when it is
absent the consensus is assumed to have been called from the alignment by the Dfam
consensus caller.

The value must be one of a reserved set of words.  At present the only reserved word is:

| Word        | Meaning                                                            |
|-------------|--------------------------------------------------------------------|
| `handbuilt` | The consensus was authored or adjusted by a curator, and is not reproducible by calling the consensus from the alignment |

```
#=GF CT    handbuilt
```

A curator hand-edits a consensus when the called consensus is not the biologically correct
one — for example when a substitution in the called consensus introduces a spurious in-frame
stop codon, or a deletion introduces a frameshift, that the ancestral element would not have
had.  In such a record the `#=GC RF` line deliberately differs from what the consensus caller
would produce from the aligned copies.

Effect on `stk lint`:

- `stk lint` normally re-calls the consensus from the alignment and warns
  (`rf_consensus_mismatch`) when `#=GC RF` disagrees with it.  When `CT` is `handbuilt`
  that comparison is **skipped**, since a disagreement is the expected state.
- A `CT` value that is not a reserved word is an error (`ct_unknown`).  The consensus check
  still runs in that case.

Note that `stk edit --update-consensus` overwrites `#=GC RF` with a freshly called consensus,
which discards the curator's hand-built version.  Since the stored consensus is then a called
one, the `CT` field no longer describes it and is removed from the record (a note is written
to stderr).

### `KD` — Kimura Divergence

A numeric measure of sequence divergence within the seed alignment.

### `BM` — Build method

The command or pipeline used to build the alignment.

```
#=GF BM    RepeatModeler 2.0
```

### `RN` / `RM` / `RD` — Literature reference block

A structured literature citation.  Each reference block begins with `RN [N]` (a sequential
reference number in square brackets), followed by at least one of:

| Tag | Meaning                                      |
|-----|----------------------------------------------|
| `RN` | Reference number — **must come first**, format `[N]` |
| `RM` | PubMed ID                                   |
| `RD` | DOI (bare, e.g. `10.1093/nar/gkl1049`, or as a `https://doi.org/` URL) |

Every reference block must have `RN` followed by at least one of `RM` or `RD`.  `stk lint`
will report an error if `RM`/`RD` is present without `RN`.  `stk edit --set` will
automatically reorder fields into the correct sequence if they are supplied out of order.

By default `stk lint` makes live HTTP requests to verify that each PubMed ID and DOI
actually resolves.  Pass `--no-network` to skip this.  bioRxiv/medRxiv version suffixes
(e.g. `10.1101/2024.01.27.577580v2`) are automatically stripped before resolving since
`doi.org` only registers the base DOI.

> **Note:** The standard Stockholm tags `RT` (title), `RA` (authors), and `RL` (journal)
> are **not imported by Dfam** and will be silently ignored.  `stk lint` warns if they are
> present.  Use `RM` or `RD` to link publications.

### `DR` — Database cross-reference

A cross-reference to an external database.

```
#=GF DR    Repbase; AluSx
#=GF DR    WikiPedia; Alu_element
```

### `CC` — Comment / curatorial note

Free-form plain-text notes intended for human readers.  May appear multiple times.

```
#=GF CC    Internal portion of a full-length MLT1 element.
#=GF CC    Consensus matches RepBase entry MLT1A0 at >95% identity.
```

### `**` — Private / curator annotation

Curator working notes not intended for public display.

---

## Special Topics

### Family Identifiers: `AC` vs `ID`

Dfam uses two identifiers with very different semantics:

**`AC` (accession)** is the stable, permanent handle assigned by Dfam.  Once a record is
accepted, its accession never changes.  Use it in scripts, pipelines, and cross-references;
it is the only identifier guaranteed to remain valid across database updates.

**`ID` (name)** is a human-friendly label.  Names can be renamed when a family is
reclassified or when a naming conflict is resolved.  They are typically preserved by
demotion to a family alias. 

For new submissions, the `AC` field will be absent until the record is accepted.  Duplicate
`ID` values in the same file are errors unless one of the records carries an `AC` (which
marks it as an update to an existing entry). 

Automated pipeline names (e.g. `rnd-3_family-101` from RepeatModeler) must be replaced with
meaningful names before submission.  They are almost always non-unique and will fail the
`stk lint` uniqueness check.

### Taxonomic Labels (`OC` and NCBI Taxonomy)

Dfam uses scientific taxon names from **NCBI Taxonomy** to label the organisms in which a
family has been found.  Every `OC` value must be a recognised scientific name:

- Use **binomial names** for species: `Homo sapiens`, `Mus musculus`
- Use **higher-rank names** for broader scope: `Mammalia`, `Eukaryota`, `Metazoa`
- Do **not** use NCBI taxon IDs — use the scientific name.

When uncertain whether a name is valid, run `stk lint --taxonomy` which will suggest
correctly spelled names via fuzzy matching or common-name lookup.

Multiple `OC` lines are acceptable and encouraged when the family is documented in several
related species.

### Classification (`TP`)

The `TP` field places the family in the **RepeatMasker classification hierarchy**, a
semicolon-delimited path:

```
Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;Retrotransposon;Long_Terminal_Repeat_Element;Gypsy-ERV;Retroviridae;Orthoretrovirinae;ERV1
Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;7SL-RNA_Promoter;No-core;L1-dependent;Alu
Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Tc1-Mariner;Tc1
```

Use `stk lint --classification <file>` to validate the `TP` value against the official Dfam
classification list.  When in doubt about the appropriate classification, consult the
[Dfam classification browser](https://dfam.org/classification) or use a broad type such as
`Interspersed_Repeat;Unknown` as a temporary placeholder.

### MSA Sequence Identifiers

Every sequence row in a Dfam seed file uses a **Smitten-format identifier** that encodes the
precise genomic location of the copy:

```
GCA_000001405.15:chr1:1234567-1234627_+
```

The components are:

| Component             | Example                | Description                      |
|-----------------------|------------------------|----------------------------------|
| Assembly accession    | `GCA_000001405.15`     | NCBI/ENA assembly accession      |
| Chromosome / scaffold | `chr1`                 | Sequence name within the assembly |
| Coordinates           | `1234567-1234627`      | 1-based, fully-closed interval |
| Strand                | `+` or `-`             | Orientation of the copy          |

Using public database accessions (rather than local file paths or custom labels) is important
because Dfam's downstream algorithms rely on being able to retrieve the flanking sequence for
target-site duplication analysis and family extension.  Instances labelled with private or
opaque names cannot be used in these analyses.

For more details on the Smitten identifier specification, including alternative formats and
version history, see the [Smitten GitHub repository](https://github.com/Dfam-consortium/Smitten.git).

---

## Validating Your File

Use the `stk lint` command to check a seed file before submission:

```sh
# Tier-1 checks only (format and field validation):
stk lint my_families.stk

# Tier-2 checks (validates TP against known classifications,
#               OC against NCBI taxonomy, ID against existing Dfam names):
stk lint --classification dfam_classification.txt \
         --taxonomy ncbi_taxonomy.txt \
         --names dfam_names.txt \
         my_families.stk
```

Tier-1 errors are always reported.  Tier-2 checks require the corresponding reference files
and validate field values against external databases.

Common issues flagged by `stk lint`:

| Check                  | Severity | Meaning                                        |
|------------------------|----------|------------------------------------------------|
| `missing_required_field` | ERROR  | `DE`, `AU`, `TP`, `OC`, or `SQ` is absent     |
| `rf_missing`           | ERROR    | `#=GC RF` line is absent                       |
| `sq_mismatch`          | ERROR    | `SQ` count does not match actual sequence rows |
| `seq_invalid_chars`    | ERROR    | Sequence row contains a character that is neither an IUB code nor a gap |
| `seq_nonstandard_gap`  | WARN     | Sequence row uses `-`, `_` or `~` as a gap instead of `.` |
| `block_format`         | ERROR    | Blank line splits the alignment into interleaved blocks |
| `id_too_long`          | ERROR    | `ID` exceeds 45 characters                     |
| `de_too_long`          | ERROR    | `DE` exceeds 80 characters                     |
| `ac_format`            | ERROR    | `AC` does not match `DF/DR` + 7 or 9 digits    |
| `duplicate_id`         | ERROR    | Same `ID` used in two records without `AC`     |
| `ac_update_records`    | INFO     | Count of records carrying an `AC` (these update existing Dfam families) |
| `orcid_missing`        | INFO     | Count of records naming an author without an ORCID |
| `au_format`            | ERROR    | `AU` token uses abbreviated/single-letter first name, `Last Initial` style, or duplicate ORCID |
| `au_format`            | WARN     | `AU` token uses abbreviated initials with periods, commas, colons, or missing space |
| `ref_block_incomplete` | ERROR    | `RM`/`RD` present without `RN`, or `RN` present without `RM`/`RD` |
| `ref_block_order`      | ERROR    | `RM`/`RD` appears before its `RN` line         |
| `citation_fields_unused` | WARN   | `RT`/`RA`/`RL` present but not imported by Dfam |
| `pmid_unknown`         | ERROR    | PubMed ID not found at pubmed.ncbi.nlm.nih.gov (network check) |
| `doi_unknown`          | ERROR    | DOI could not be resolved via doi.org (network check) |
| `rf_consensus_mismatch`| WARN     | `#=GC RF` does not match the consensus called from the alignment (not checked when `CT` is `handbuilt`) |
| `ct_unknown`           | ERROR    | `CT` value is not a reserved consensus-type word |
| `mm_invalid_chars`     | ERROR    | `#=GC MM` contains a character other than `m` or `.` |
| `mm_length_mismatch`   | ERROR    | `#=GC MM` width differs from `#=GC RF` / the alignment |
| `unknown_gc_tag`       | INFO     | `#=GC` annotation other than `RF` or `MM` (not imported) |
| `oc_unknown`           | ERROR    | `OC` value not found in NCBI taxonomy          |
| `tp_unknown`           | ERROR    | `TP` value not found in Dfam classification    |
| `id_in_dfam`           | ERROR    | `ID` already exists in Dfam (add `AC` to update)|
| `unknown_gf_tag`       | INFO     | Unrecognised `#=GF` tag (possible typo)        |

To rebuild the `#=GC RF` line from the alignment before linting:

```sh
stk edit --update-consensus my_families.stk > my_families_rf.stk
stk lint my_families_rf.stk
```
