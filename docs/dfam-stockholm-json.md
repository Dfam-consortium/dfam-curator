# Dfam Stockholm JSON (version 1.0)

Dfam Stockholm JSON is a JSON form of a Dfam seed alignment file. It carries the
same content as the Stockholm 1.0 format described in [Dfam_Seeds.md](../Dfam_Seeds.md),
but names the fields instead of using two-letter tags. It is meant for programs that
ingest seed alignments: validate the Stockholm with `stk lint`, then convert it with
`stk convert --to json` and read the JSON.

```sh
stk lint families.stk
stk convert --to json -o families.json families.stk
```

Stockholm remains Dfam's interchange format, so `stk convert` is the only
subcommand that reads or writes JSON. It recognises a JSON input by its content (a
file whose first non-blank character is `{`), and `stk convert families.json`
writes it back out as Stockholm. The other subcommands (`lint`, `extract`, `edit`,
`map-consensus`) refuse a JSON input with an error that points to `convert`.

Lint in particular should run on the Stockholm that is submitted. Converting to JSON
fixes the record's layout (field order, `SQ`, reference blocks), so linting the
JSON would miss errors that the Stockholm source has.

A JSON Schema (draft 2020-12) for the format is in
[dfam-stockholm-json.schema.json](dfam-stockholm-json.schema.json). The schema
checks structure only. The Dfam content rules (required fields, the length of `ID`
and `DE`, the classification and taxonomy lists, and so on) are checked by `stk lint`
on the Stockholm.

---

## Document

```json
{
  "format": "dfam-stockholm-json",
  "version": "1.0",
  "records": [ { ... }, { ... } ]
}
```

| Member    | Type   | Meaning |
|-----------|--------|---------|
| `format`  | string | Always `"dfam-stockholm-json"` |
| `version` | string | `"MAJOR.MINOR"`. Readers accept any version with the major number they know; `stk` reads `1.x` |
| `records` | array  | One object per Stockholm record (`# STOCKHOLM 1.0` … `//`), in file order |

Unknown members are an error at every level, so a misspelt field name is reported
instead of being dropped. A later minor version may add optional members. A reader
written for 1.0 that wants to accept such documents should ignore members it does not
know rather than reject them.

---

## Record

Every member is optional except `sequences`. An absent member means the field was
absent from the Stockholm record. A field with no value is left out of the JSON
entirely; it is never written as `null` or `[]`.

| Member                    | Stockholm  | Type | Notes |
|---------------------------|------------|------|-------|
| `name`                    | `#=GF ID`  | string | |
| `accession`               | `#=GF AC`  | string | e.g. `"DF000000001.2"` |
| `description`             | `#=GF DE`  | string | |
| `authors`                 | `#=GF AU`  | array of [author](#author) | Every author on every `AU` line, in order |
| `source`                  | `#=GF SE`  | string | |
| `classification`          | `#=GF TP`  | string | Full path or shorthand alias, as written |
| `clades`                  | `#=GF OC`  | array of string | One entry per `OC` line |
| `target_site_duplication` | `#=GF TD`  | string | |
| `consensus_type`          | `#=GF CT`  | string | e.g. `"handbuilt"` |
| `kimura_divergence`       | `#=GF KD`  | number | |
| `build_method`            | `#=GF BM`  | string | |
| `references`              | `#=GF RN` blocks | array of [reference](#reference) | |
| `database_references`     | `#=GF DR`  | array of string | One entry per `DR` line, e.g. `"Repbase; AluSx"` |
| `comments`                | `#=GF CC`  | array of string | One entry per `CC` line |
| `features`                | `#=GF FT`  | array of [feature](#feature) | |
| `curator_notes`           | `#=GF **`  | array of string | One entry per `**` line |
| `other_gf`                | other `#=GF` | array of `{"tag", "value"}` | Tags with no named member, in document order |
| `reference_annotation`    | `#=GC RF`  | string | Gapped: one character per alignment column. See below |
| `model_mask`              | `#=GC MM`  | string | Gapped: one character per alignment column |
| `other_gc`                | other `#=GC` | object, tag → string | e.g. `{"SS_cons": "..."}` |
| `sequences`               | sequence rows | array of [sequence](#sequence) | In document order |

`reference_annotation` holds the `#=GC RF` line as written: one character per
alignment column, gaps included. It takes its name from HMMER, which defines `RF` as
the reference annotation. The member is deliberately not called `consensus`: in Dfam,
the family consensus is this line with the gaps removed. Older files may hold `x`/`X`
occupancy marks here instead of bases (see [The `#=GC RF` line](../Dfam_Seeds.md#the-gc-rf-line)).

`#=GF SQ` has no member. It is the length of `sequences`, and it is written from that
length on the way back to Stockholm.

String values are the Stockholm values with surrounding whitespace removed, which is
how `stk` reads them.

### Author

```json
{ "name": "Barbara McClintock", "orcid": "0000-0001-2345-6789" }
```

One `AU` line can hold several authors separated by `;`, and each may be prefixed
with `ORCID:<id>`. In the JSON each author is an object. `orcid` is the bare identifier
without the `ORCID:` prefix, and it is left out when the author has none. `name` is
the rest of the author's text, as written.

### Reference

```json
{ "number": 1, "pubmed_id": "12490638", "doi": "10.1093/nar/gkl1049" }
```

A reference is an `RN [n]` line and the `RM`, `RD`, `RT`, `RA` and `RL` lines that follow
it. `number` is required and the other members are optional.

| Member      | Stockholm |
|-------------|-----------|
| `number`    | `RN [n]` (the integer `n`) |
| `pubmed_id` | `RM` |
| `doi`       | `RD` |
| `title`     | `RT` (not imported by Dfam) |
| `authors`   | `RA` (not imported by Dfam) |
| `location`  | `RL` (not imported by Dfam) |

### Feature

```json
{
  "key": "CDS",
  "location": "join(312..602,1346..2053)",
  "qualifiers": [
    { "name": "product", "value": "Gypsy-13_AnMou-I_1p" },
    { "name": "pseudo",  "value": null }
  ]
}
```

A feature is an `FT` key + location line and the `/qualifier` lines below it.
`location` is the location string as written. Its grammar and coordinate system are
described under [FT — Feature table](../Dfam_Seeds.md#ft--feature-table). In short, it
is 1-based on the ungapped consensus, and descending coordinates mean the minus strand.
A qualifier's `value` has the surrounding quotes removed. A qualifier with no value,
such as `/pseudo`, has `"value": null`.

### Sequence

```json
{
  "id": "GCA_000001405.15:chr1:1234567-1234627_+",
  "aligned": "CAGACTTGATGCACGACG.TAAACGTG",
  "coordinates": {
    "assembly": "GCA_000001405.15",
    "sequence": "chr1",
    "start": 1234567,
    "end": 1234627,
    "strand": "+"
  }
}
```

| Member        | Meaning |
|---------------|---------|
| `id`          | The row name exactly as written. This is the authoritative identifier |
| `aligned`     | The gapped sequence, as wide as the alignment |
| `coordinates` | The [Smitten](https://github.com/Dfam-consortium/Smitten) identifier parsed out of `id`. It is present only when `id` has a coordinate range. It is written for the reader's convenience and is **ignored on input** |

In `coordinates`, `start` and `end` are 1-based and fully closed, with `start <= end` on
either strand. `assembly` and `strand` are left out when `id` does not give them. When
`id` has nested ranges, `coordinates` describes the outermost one; use `id` for the
full detail.

---

## Converting between Stockholm and JSON

**Stockholm → JSON.** `stk convert --to json` cleans each record exactly as `--to stk`
does (gaps to `.`, a 7-digit `AC` widened to 9 digits), unless `--no-clean` is given.
It refuses, with an error that names the record, any record whose structure the JSON
cannot hold:

- a second `ID`, `AC`, `DE`, `SE`, `TP`, `TD`, `CT`, `KD`, `BM` or `SQ` line
- an `RM`/`RD`/`RT`/`RA`/`RL` line before any `RN`, or a repeated one within a block
- an `RN` value that is not `[n]`
- an `FT` qualifier with no feature above it, or an `FT` feature with no location
- a `KD` that is not a number
- an `SQ` that does not match the number of sequence rows

`stk lint` reports each of these as an error, so a file that lints without errors
always converts. When a record is refused, no output file is left behind.

**JSON → Stockholm.** Stockholm read from JSON has the same content but a canonical
layout:

- `#=GF` lines come in this order: `ID AC DE AU SE TP OC SQ TD CT KD BM`, then each
  reference block (`RN RM RD RT RA RL`), then `DR CC FT **`, then `other_gf` in its
  own order.
- All authors are written on one `AU` line, joined with `; `.
- `SQ` is the number of sequences.
- Qualifier values are written quoted (`/name="value"`).
- `KD` is written in its shortest exact decimal form (`12.50` becomes `12.5`).

A JSON value that cannot be written on one Stockholm line is refused: a line break in
any string, or whitespace in a sequence `id` or in an `other_gf`/`other_gc` tag.

So Stockholm → JSON → Stockholm keeps every field and every row, and a JSON document
converts to Stockholm and back to exactly the same JSON.

Lines that `stk` does not keep when reading Stockholm (`#=GS`, `#=GR`, and `#` lines
that are not annotations) are not part of the JSON either.

---

## Example

The [more complete example](../Dfam_Seeds.md#more-complete-example) in Dfam_Seeds.md
converts to:

```json
{
  "format": "dfam-stockholm-json",
  "version": "1.0",
  "records": [
    {
      "description": "MER3822 internal sequence",
      "authors": [ { "name": "Robert Smith" } ],
      "classification": "Interspersed_Repeat;LTR/ERV1",
      "clades": [ "Bos taurus", "Homo sapiens" ],
      "references": [ { "number": 1, "pubmed_id": "12490638" } ],
      "database_references": [ "Repbase; MER4-int" ],
      "comments": [
        "Internal sequence of the MER4 family of LTR retrotransposons.  This",
        "family can be described fully here using any number of CC lines."
      ],
      "reference_annotation": "CAGACTTGATGCACGACGTAAACGTGACT.GGATCGCCCTAGATCGCGATCGATCGATCG",
      "sequences": [
        {
          "id": "GCA_000001405.15:chr1:1234567-1234627_+",
          "aligned": "CAGACTTGATGCACGACGTAAACGTGACTCGGATCGCCCTAGATCGCGATCGATCGATCG",
          "coordinates": { "assembly": "GCA_000001405.15", "sequence": "chr1", "start": 1234567, "end": 1234627, "strand": "+" }
        },
        {
          "id": "GCA_000001405.15:chr5:9876543-9876603_-",
          "aligned": "CAGACTTGATGCACGACGTAAACGTGACT.GGATCGCCCTAGATCGCGATCGATCGATCG",
          "coordinates": { "assembly": "GCA_000001405.15", "sequence": "chr5", "start": 9876543, "end": 9876603, "strand": "-" }
        },
        {
          "id": "GCA_000001405.15:chr7:5555000-5555060_+",
          "aligned": "CAGACTTGATGCCCGACGTAAACGTGACT.GGATCGCCCTAGATCGCGATCGATCGATCG",
          "coordinates": { "assembly": "GCA_000001405.15", "sequence": "chr7", "start": 5555000, "end": 5555060, "strand": "+" }
        }
      ]
    }
  ]
}
```

(`stk` writes one member per line; the layout above is compacted for reading.)

---

## Notes for readers

- A JSON document is read whole, so reading one takes memory proportional to its
  size. Writing streams one record at a time.
