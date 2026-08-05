"""
gtf_pyparser — lightweight GTF parsing library
==============================================

A small, dependency-free library for parsing GTF (Gene Transfer Format) files
into an in-memory object model. The primary use case is bioinformatics pipelines
that need to query gene structure — transcripts, exons, CDS boundaries, and
derived features such as introns — without pulling in a heavy framework.

Coordinate convention
---------------------
All coordinates in this library use **0-based half-open** intervals
(start inclusive, end exclusive), matching BED/BAM convention.
GTF files are 1-based inclusive; the conversion is applied automatically
during parsing.

Typical usage
-------------
Load a GTF file and iterate over genes and their transcripts::

    import gtf_pyparser

    genes = gtf_pyparser.parse_gtf("Homo_sapiens.GRCh38.gtf")

    gene = genes["ENSG00000139618"]
    print(gene.symbol, gene.start, gene.end)

    for transcript_id, transcript in gene.transcripts.items():
        exons = transcript.features.get("exon", [])
        introns = gtf_pyparser.get_intron(exons)
        print(f"  {transcript_id}: {len(exons)} exons, {len(introns)} introns")

Silencing log output
--------------------
The parser emits progress and error messages through the standard
:mod:`logging` module under the logger named ``"gtf_pyparser"``. To suppress
informational output::

    import logging
    logging.getLogger("gtf_pyparser").setLevel(logging.WARNING)

Public API
----------
Data classes:

- :class:`Interval`   — immutable genomic coordinate record
- :class:`Transcript` — transcript isoform with grouped feature intervals
- :class:`Gene`       — gene locus containing one or more transcripts

Functions:

- :func:`parse_gtf`   — parse a GTF file into a dict of Gene objects
- :func:`get_intron`  — derive intron intervals from a list of exon intervals
- :func:`get_attr`    — parse a raw GTF attribute string into a key-value dict

The following are considered internal and not part of the stable API:

- :meth:`Gene._parse_transcript_line` — incremental transcript builder
"""

from .gtf_pyparser import (
    Interval,
    Transcript,
    Gene,
    get_intron,
    get_attr,
    gtf_to_dict as parse_gtf,
)
try:
    from ._version import version as __version__
except ImportError:
    __version__ = "unknown"

__all__ = [
    "Interval",
    "Transcript",
    "Gene",
    "get_intron",
    "get_attr",
    "parse_gtf",
    "__version__",
]