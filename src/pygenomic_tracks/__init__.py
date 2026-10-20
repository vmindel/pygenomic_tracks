"""Genome-browser-style tracks from pandas DataFrames, drawn on matplotlib axes.

Every track function has the shape ``fn(df, ax, region, **column_args)``.
DataFrames need ``chrom``, ``start``, ``end`` (0-based, half-open; override
with ``cols=``), which is what pymisha's ``gextract``/``gscreen`` return.
"""
from ._coords import Region, as_region, select
from ._pack import pack_rows
from .annotations import annotation, cache_dir, cached, download, gene_region
from .datasets import (
    example_genes,
    example_gtf_path,
    example_motif_hits,
    example_peaks,
    example_region,
    example_signal,
)
from .genes import genes
from .intervals import intervals
from .io import from_bigwig, read_bed, read_gtf
from .layout import figure, format_genomic_axis
from .signal import signal

__version__ = "0.1.0"

__all__ = [
    "Region",
    "annotation",
    "as_region",
    "cache_dir",
    "cached",
    "download",
    "example_genes",
    "example_gtf_path",
    "example_motif_hits",
    "example_peaks",
    "example_region",
    "example_signal",
    "figure",
    "format_genomic_axis",
    "from_bigwig",
    "gene_region",
    "genes",
    "intervals",
    "pack_rows",
    "read_bed",
    "read_gtf",
    "select",
    "signal",
]
