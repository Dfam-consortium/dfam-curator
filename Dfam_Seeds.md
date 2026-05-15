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
| `#=GC`  | Per-**c**olumn annotation (one value per alignment column) |
| `#=GS`  | Per-**s**equence metadata (not used in Dfam) |
| `#=GR`  | Per-**r**esidue annotation (not used in Dfam) |

Stockholm has **no free-text comment syntax**.  Every line beginning with `#` must use one of
the structured prefixes above — plain `# comment` lines are not valid inside a record.

### Sequence rows

Each aligned sequence is written as a name, whitespace, and the aligned sequence string:

```
GCA_000001405.15:chr1:100-200:+    ACGT..ACGT
```

- Gap characters are written as `.` (period).  The `-` character is **not valid** in Dfam
  seed files; `stk lint` will flag it as an error.
- All sequence rows in a record must have the same aligned length.
- Dfam has adopted the **Smitten** identifier format (see
  [MSA Sequence Identifiers](#msa-sequence-identifiers) below).

### The `#=GC RF` line

The RF (reference/consensus) line is a per-column annotation that encodes the consensus
sequence of the family (derived from the MSA).  It uses IUPAC/IUB ambiguity codes for 
columns with mixed bases and `.` for gap-only columns.  In older versions of Dfam 
Stockholm files the characters `X` or `x` were used to indicate consensus occupancy 
columns for use with the HMMER hmmbuild tool.

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
GCA_000001405.15:chr1:1234567-1234627:+    CAGACTTGATGCACGACGTAAACGTGACTCGGATCGCCCTAGATCGCGATCGATCGATCG
GCA_000001405.15:chr5:9876543-9876603:-    CAGACTTGATGCACGACGTAAACGTGACT.GGATCGCCCTAGATCGCGATCGATCGATCG
GCA_000001405.15:chr7:5555000-5555060:+    CAGACTTGATGCCCGACGTAAACGTGACT.GGATCGCCCTAGATCGCGATCGATCGATCG
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

- Accessions beginning with `DF` are curated families.
- Accessions beginning with `DR` or "Dfam Raw" families are uncurated families 
  that have been generated by de novo TE identification tools.
- Records submitted for the first time will not have an `AC` field yet.  Once assigned, `AC`
  becomes the stable, authoritative identifier for the family.
- Do not invent an accession; leave the field absent until one is officially assigned.

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

Who created or curated this seed alignment.  Free-form text; conventional format is
`Last F` (last name then initial) or multiple authors separated by semicolons.

```
#=GF AU    Smith R
#=GF AU    Jayapal P; Ruef W
```

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

### `KD` — Kimura Divergence

A numeric measure of sequence divergence within the seed alignment.

### `BM` — Build method

The command or pipeline used to build the alignment.

```
#=GF BM    RepeatModeler 2.0
```

### `RN` / `RM` / `RT` / `RA` / `RL` — Literature reference block

A structured literature citation.  Each reference block begins with `RN [N]` (a sequential
reference number in brackets), followed by:

| Tag | Meaning             |
|-----|---------------------|
| `RN` | Reference number (e.g. `[1]`) — **must come first** |
| `RM` | PubMed ID (required)           |
| `RT` | Title               |
| `RA` | Author(s)           |
| `RL` | Journal and volume  |

`RN` must be the first line in each reference block, followed by `RM` (required), then
optionally `RT`, `RA`, and `RL` in that order.  `stk lint` will report an error if any
of `RM`/`RT`/`RA`/`RL` appears before the `RN` line, and a warning if `RN` or `RM` is
present without the other.  `stk edit --set` will automatically reorder fields into the
correct sequence if they are supplied out of order.

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
GCA_000001405.15:chr1:1234567-1234627:+
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
| `seq_invalid_chars`    | ERROR    | Sequence row contains `-` or `~` instead of `.`|
| `id_too_long`          | ERROR    | `ID` exceeds 45 characters                     |
| `de_too_long`          | ERROR    | `DE` exceeds 80 characters                     |
| `ac_format`            | ERROR    | `AC` does not match `DF/DR` + 7 or 9 digits    |
| `duplicate_id`         | ERROR    | Same `ID` used in two records without `AC`     |
| `rf_consensus_mismatch`| WARN     | `#=GC RF` does not match the consensus called from the alignment |
| `oc_unknown`           | ERROR    | `OC` value not found in NCBI taxonomy          |
| `tp_unknown`           | ERROR    | `TP` value not found in Dfam classification    |
| `id_in_dfam`           | ERROR    | `ID` already exists in Dfam (add `AC` to update)|
| `unknown_gf_tag`       | INFO     | Unrecognised `#=GF` tag (possible typo)        |

To rebuild the `#=GC RF` line from the alignment before linting:

```sh
stk edit --update-consensus my_families.stk > my_families_rf.stk
stk lint my_families_rf.stk
```
