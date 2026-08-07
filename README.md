# gtf_pyparser

A lightweight, dependency-free Python library for parsing GTF (Gene Transfer Format) files into a structured object model.

## Why gtf_pyparser?

Most GTF parsing libraries are either too heavy (requiring pandas, SQLite, or large C extensions) or return flat data structures that don't reflect the natural hierarchy of genomic annotation. `gtf_pyparser` gives you a clean `Gene → Transcript → features` object model with no dependencies beyond the Python standard library.

### new V0.3+
added intervall tree support. you can now query genes over a position.

### Roadmap

Currently working on direct access / indexing by name and position.

## Installation

```bash
pip install gtf_pyparser
```

Position lookups (`Gtf.build_intervals` / `Gtf.get_genes_at_position`) require the optional `intervaltree` dependency:

```bash
pip install gtf_pyparser[intervals]
```

## Quick start

```python
import gtf_pyparser

genes = gtf_pyparser.parse_gtf("Homo_sapiens.GRCh38.gtf")

gene = genes["ENSG00000139618"]
print(gene.symbol, gene.start, gene.end)

for transcript_id, transcript in gene.transcripts.items():
    exons = transcript.features.get("exon", [])
    introns = gtf_pyparser.get_intron(exons)
    print(f"  {transcript_id}: {len(exons)} exons, {len(introns)} introns")
```

## Coordinate convention

All coordinates are **0-based half-open** (start inclusive, end exclusive), consistent with BED and BAM format. GTF files are 1-based inclusive; the conversion is applied automatically during parsing.

## Data model

```
Gtf
├── dict: dict[primary_key, Gene]
│   └── Gene
│       ├── interval: Interval  // start, end, chr, strand, phase, attributes.
│       ├── attribute: dict
│       └── transcripts: dict[transcript_id, Transcript]
│           └── Transcript
│               ├── interval: Interval
│               ├── attribute: dict
│               └── features: dict[feature_type, list[Interval]]   // exon / CDS / UTR / ...
└── it: dict[chr, IntervalTree]   // lazily built, keyed by gene_id
```

## API

### Parsing

#### `gtf_pyparser.parse_gtf(gtf_file, primary_key="gene_id")`

Parse a GTF file into a dict-like container of `Gene` objects, keyed by `primary_key`.

```python
genes = gtf_pyparser.parse_gtf("annotation.gtf")
genes = gtf_pyparser.parse_gtf("annotation.gtf", primary_key="gene_name")
```

- **`gtf_file`** — path to the GTF file
- **`primary_key`** — attribute used to key the returned container (default: `"gene_id"`)
- **Returns** — a `Gtf` object supporting indexing, `in`, `len`, iteration, and `.keys()`/`.values()`/`.items()`, exactly like `dict[str, Gene]`

#### `gtf_pyparser.get_attr(string)`

Parse a raw GTF attribute string into a key-value dictionary. Useful when processing GTF lines outside of the full parser.

```python
attrs = gtf_pyparser.get_attr('gene_id "ENSG00000139618"; gene_name "BRCA2";')
# {"gene_id": "ENSG00000139618", "gene_name": "BRCA2"}
```

Repeated attribute keys (e.g. multiple `tag "..."` entries) are collected into a list instead of overwriting each other.

### Position lookups

*Requires the `intervaltree` package — see [Installation](#installation).*

`Gtf.get_genes_at_position(position, flanking=0)` returns the `gene_id` of every gene overlapping a genomic position. The underlying `IntervalTree` index is built lazily on first use (or explicitly via `Gtf.build_intervals()`).

```python
genes = gtf_pyparser.parse_gtf("annotation.gtf")

genes.get_genes_at_position({"chr": "chr1", "start": 100_000, "end": 100_050})
# ["ENSG00000139618"]

genes.get_genes_at_position({"chr": "chr1", "start": 100_000, "end": 100_050}, flanking=500)
```

- **`position`** — dict-like with `chr`, `start`, `end` keys, `start < end`
- **`flanking`** — bases to extend the query on both sides (default `0`)
- **Returns** — `list[str]` of overlapping `gene_id`s, or `None` if the chromosome isn't in the GTF or nothing overlaps

### Derived features

#### `gtf_pyparser.get_intron(exons)`

Derive intron intervals from a list of exon `Interval` objects belonging to a single transcript. Introns are numbered biologically: from 1 upward on the `+` strand, and from n down to 1 on the `-` strand.

```python
exons = transcript.features.get("exon", [])
introns = gtf_pyparser.get_intron(exons)
```

- **`exons`** — list of `Interval`, need not be pre-sorted
- **Returns** — `list[Interval]`, empty if fewer than two exons are provided

### Position classification

`Transcript.classify_position(position, strand, strand_aware=True)` and `Gene.classify_position(position, strand, strand_aware=True)` classify a 0-based genomic position relative to a transcript's (or every transcript of a gene's) exon/intron structure.

```python
gene.classify_position(150, "+")
# {"ENST00000001": "exon", "ENST00000002": "intron"}
```

Returns one of `"geneStart"`, `"geneEnd"`, `"exon"`, `"intron"`, `"junctionDonnor"`, `"junctionAcceptor"`, or `None` if the position falls outside the feature's span. `Gene.classify_position` returns a `dict[transcript_id, result]` covering every transcript on the gene.

### Data classes

#### `Interval`

Immutable (frozen) genomic coordinate record — replaced rather than mutated when coordinates need updating.

| Attribute   | Type          | Description                                              |
|-------------|---------------|----------------------------------------------------------|
| `chr`       | `str`         | Chromosome or sequence name                              |
| `start`     | `int`         | 0-based start (inclusive)                                |
| `end`       | `int`         | 0-based end (exclusive)                                  |
| `strand`    | `str`         | `"+"`, `"-"`, or `"."` for unstranded                   |
| `phase`     | `int` or `str`| Reading frame (0, 1, 2), or `"."` if not applicable     |
| `attribute` | `dict`        | Key-value pairs from the GTF attribute column            |

Other properties: `length`, `position` (dict with just chr/start/end/strand).

Methods:

- `overlaps(other, strand_aware=True, closed=False)` / `contains(position, strand, strand_aware=True, closed=False)` — half-open by default (touching intervals don't overlap; `position == end` isn't contained); pass `closed=True` to treat both endpoints as inclusive.
- `eq_pos(other)` — compare genomic position and strand only, ignoring phase/attribute.
- `clone()` — return a copy with its own independent `attribute` dict.
- `to_dict()` / `Interval.from_dict(dict_)` — round-trip through a plain dict.
- `Interval.from_position_str(string, is_one_based=False)` — parse a `"chr:start-end(strand)"` string; pass `is_one_based=True` to convert from 1-based inclusive (GTF/samtools region notation) to this class's 0-based half-open convention.

#### `Transcript`

Groups all GTF records sharing a `transcript_id`. The genomic span always reflects the union of all features seen so far.

| Attribute           | Type                          | Description                              |
|---------------------|-------------------------------|------------------------------------------|
| `transcript_id`     | `str`                         | Ensembl transcript ID or equivalent      |
| `transcript_symbol` | `str` or `None`               | Human-readable transcript symbol         |
| `interval`          | `Interval`                    | Current genomic span of the transcript   |
| `features`          | `dict[str, list[Interval]]`   | Feature type → list of intervals         |

Convenience properties: `start`, `end`, `chr`, `phase`, `length`, `attribute` (of the transcript's own interval), `exons`, `exons_positions`, `intron` (derived via `get_intron`), `cds`, `mrna`, `utr_5p`, `utr_3p`.

#### `Gene`

A gene locus containing one or more transcript isoforms.

| Attribute     | Type                      | Description                                                        |
|---------------|---------------------------|--------------------------------------------------------------------|
| `gene_id`     | `str`                     | Primary identifier (e.g. Ensembl gene ID)                          |
| `symbol`      | `str`                     | Gene symbol, resolved from `gene_symbol`, `gene_name`, or `gene_id`|
| `interval`    | `Interval`                | Genomic span from the GTF `gene` record                            |
| `transcripts` | `dict[str, Transcript]`   | Transcript ID → Transcript                                         |

Convenience properties: `start`, `end`, `chr`, `phase`, `length`, `attribute`, `biotype` (checks `biotype`, `gene_biotype`, then `gene_type`), `transcript_names`, `transcript_length`, `has_transcript`, `has_exon`, `exon` / `intron` (list of `(gene_id, transcript_id, ...)` tuples across all transcripts).

## Interoperating with easyfasta

`gtf_pyparser` has no dependency on [`easyfasta`](https://github.com/rLannes/easyfasta) (or any other sequence library) — the two interoperate purely by shape. `Interval` supports dict-style access (`interval["start"]`, `interval.get("strand")`), which is exactly the protocol `easyfasta.fai_common.query_position`/`query_iter`/`query_splice` expect from a "position": anything with `["chr"]`, `["start"]`, `["end"]`, and an optional `["strand"]`. That means `Interval` objects — including a transcript's `.exons` — can be passed straight through, no glue code or extra dependency required:

```python
import gtf_pyparser
from easyfasta import fai_common

genes = gtf_pyparser.parse_gtf("annotation.gtf")
transcript = genes["ENSG00000139618"].transcripts["ENST00000001"]

# a single feature
cds_seq = fai_common.query_position("genome.fa", transcript.cds[0])

# concatenate all exons into the spliced transcript sequence
mrna_seq = fai_common.query_splice("genome.fa", transcript.exons)
```

`Interval.to_dict()` / `Interval.from_dict()` round-trip through the same `chr`/`start`/`end`/`strand` shape, so plain dicts work interchangeably with `Interval` objects anywhere this protocol is expected.

## Logging

Progress and error messages are emitted under the `"gtf_pyparser"` logger. To suppress informational output:

```python
import logging
logging.getLogger("gtf_pyparser").setLevel(logging.WARNING)
```

## Versioning

`gtf_pyparser.__version__` reflects the installed package's version, derived from git tags via `setuptools_scm`. In an unbuilt source checkout (no install step run yet) it falls back to `"unknown"`.

## Development

```bash
pip install -e ".[dev]"
pytest
```

## License
MIT

Citation / acknowledgement: Romain Lannes 2026