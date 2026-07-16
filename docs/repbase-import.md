# stk repbase-import (internal)

Converts a Repbase family into a Dfam Stockholm record by combining two
IG-derived input files:

- an **IG MSA** (the aligned sequences, `--msa`), and
- an **IG family record** (the metadata + consensus, `--record`).

This subcommand is intended for internal Repbase→Dfam conversion work and is
deliberately left out of the main `README.md`.

```
USAGE:
    stk repbase-import --msa <FILE> --record <FILE> [OPTIONS]

OPTIONS:
    --msa <FILE>         Repbase IG MSA (aligned FASTA with `; FRAGMENT` comments)
    --record <FILE>      Repbase IG family record (EMBL-style metadata + consensus)
    -o, --output <FILE>  Write output to FILE instead of stdout
    --no-clean           Do not normalize the record on write (see below)
```

```sh
stk repbase-import --msa Mariner-N5_CyaStr.aln --record Mariner-N5_CyaStr -o Mariner-N5_CyaStr.stk
```

## Input formats

Both formats are loosely based on IntelliGenetics/Stanford (IG) format: comment
lines begin with `;`, and there is no `1`/`2` sequence terminator.

### IG MSA (`--msa`)

Aligned FASTA with a `; FRAGMENT <start> -> <end>` comment (and a blank `;`
comment) preceding each row. The first record is the **consensus** (it becomes
the MSA reference / `#=GC RF`); the remaining records are the aligned instances.

```
; FRAGMENT 1 -> 383
;
Mariner-N5_CyaStr
CTGGATAATTTCGACC....----TAAGCC....ATTATCCAG
; FRAGMENT 1 -> 383
;
JAOVFP010000019.1_1
CTGGATAATTTCGACC....----TAAGCC....ATTATTCAG
```

`FRAGMENT a -> b` becomes the row's `seq_start`/`seq_end` (orientation reverse
when `a > b`). These coordinates are taken as-is and are **not** assumed
accurate — they are expected to be verified/repaired later by `discoord`.

The IG MSA can also be opened directly by `linup` and `stk convert` (it is a
recognized input format), which is how the "view an IG MSA" use case is served.

### IG family record (`--record`)

EMBL-style content where every metadata line is prefixed with `;` and carries a
two-letter tag (`ID`, `DE`, `KW`, `OS`, `OC`, `RN`/`RA`/`RT`/`RL`, `CC`, `SQ`, …;
`XX` lines are separators). After the `;SQ` line, a bare identifier line
introduces the wrapped (ungapped) consensus sequence.

```
;ID   Mariner-N5_CyaStr DNA   ; PLN   ; 6225 BP
;DE   DNA transposon from the Cyathus striatus genome, consensus.
;KW   Mariner/Tc1; DNA transposon; Transposable Element; nonautonomous;
;KW   Mariner-N5_CyaStr.
;OS   Cyathus striatus
;RN   [1]  ()
;RA   Bao,W.
;RT   DNA transposons from the Cyathus striatus genome.
;RL   Direct Submission to RR (8-Jul-2026)
;CC   ~96% identical to consensus.
;SQ   Sequence 6225 BP; ...
Mariner-N5_CyaStr
CTGGATAATTTCGACC...
```

## Field mapping

| Repbase (IG) | Dfam (Stockholm) | Notes |
|---|---|---|
| `ID` identifier token | `#=GF ID` | first whitespace token only |
| *(boilerplate)* | `#=GF DE` | always `Repbase TE family` (Repbase `DE` → `CC`) |
| `KW` (up to and incl. `Transposable Element`) | `#=GF TP` | via the classification table; unmapped → `TP` omitted + warning |
| `OS` species name | `#=GF OC` | taxon name only; the `;OC` lineage is dropped |
| `RN`/`RA`/`RT`/`RL` | reference block *or* `AU` + `**` | see **Reference handling** |
| `DE`, then `CC` | `#=GF CC` | |
| MSA consensus (first row) | `#=GC RF` | |
| MSA instances | sequence rows | id `name:start-end_orient` from `FRAGMENT` |

### Classification (`KW` → `#=GF TP`)

The lookup key is the `KW` tokens up to **and including** `Transposable Element`
(e.g. `Mariner/Tc1; DNA transposon; Transposable Element`). It is matched
(case-insensitively) against a seed table in `src/dfam/repbase.rs`
(`KW_CLASS_TABLE`). Extend that table with `(key, TP-lineage)` pairs as new
superfamilies are encountered. When a key is not found, `TP` is left unset and a
warning is emitted.

### Reference handling

A **Repbase Reports direct submission** (an `RL` containing "Repbase Reports" or
a standalone `RR` token) is not treated as a citation:

- its author (`RA`) becomes the Dfam curator → `#=GF AU`
- its date (from the `RL` parentheses, else the record `DT`) becomes a curator
  note → `#=GF ** Repbase Reports Submission Date: <date>`

All other references are emitted as normal Stockholm reference blocks
(`RN` renumbered sequentially, then `RT`/`RA`/`RL`).

## Consensus validation (advisory)

The importer runs two consensus checks and prints warnings to stderr. They do
not change the output — they exist to surface divergent examples for future
handling.

1. **MSA vs record** — the `--msa` consensus (first row, ungapped) and the
   `--record` file's own consensus should be the same sequence. Only the import
   can check this, since neither `lint` nor the output STK sees both files.
2. **RF vs called consensus** — the consensus (`#=GC RF`) should reproduce the
   consensus called from the aligned instances. This is the same test `stk lint`
   reports as `rf_consensus_mismatch`, run here (via the shared
   `lint::rf_consensus_status`) so the import alerts without a separate lint pass.

## Clean-on-write

By default the written record is normalized to Dfam conventions (gap characters
`-`, `_`, `~` → `.`; 7-digit accessions widened to 9 digits), consistent with
`stk edit`/`extract`/`convert`. Pass `--no-clean` to suppress this.

## Known gaps

- **`#=GF AU` format** — Repbase authors use `Last,Initial` (e.g. `Bao,W.`),
  which `stk lint` flags with `au_format`; this is left for the curator to clean
  up (and to add an ORCID). The importer does not rewrite the author.
- **Divergent consensus** — when the consensus genuinely differs (either between
  the two files or from the called consensus), the importer only *warns*; the
  correct handling is still to be designed.
