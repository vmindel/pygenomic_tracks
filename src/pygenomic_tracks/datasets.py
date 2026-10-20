"""Small synthetic datasets for the tutorials and tests.

All tables describe one made-up locus, ``chr2:1,200,000-1,350,000``, and use
the same shape pymisha returns: ``chrom``, ``start``, ``end`` plus value
columns. They are generated deterministically, so the figures in the docs
are reproducible.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from ._coords import Region
from .io import read_gtf

_PEAKS = [(0.18, 40.0), (0.42, 25.0), (0.47, 18.0), (0.81, 55.0)]


def example_region() -> Region:
    """The example locus, ``chr2:1,200,000-1,350,000``."""
    return Region("chr2", 1_200_000, 1_350_000)


def example_signal(kind: str = "chec", bin_size: int = 50) -> pd.DataFrame:
    """Binned signal over the example locus.

    Parameters
    ----------
    kind : {"chec", "control"}
        ``"chec"`` has four peaks on a gamma background and a stretch of
        missing bins (NaN); ``"control"`` is background only.
    bin_size : int
        Bin width in bp.

    Returns
    -------
    pandas.DataFrame
        Columns ``chrom``, ``start``, ``end``, ``signal``, ``intervalID``,
        as from ``pymisha.gextract(..., iterator=bin_size)``.
    """
    if kind not in ("chec", "control"):
        raise ValueError(f"kind must be 'chec' or 'control', got {kind!r}")
    region = example_region()
    rng = np.random.default_rng(0 if kind == "chec" else 1)
    starts = np.arange(region.start, region.end, bin_size)
    x = (starts - region.start) / region.length
    sig = rng.gamma(1.2, 1.0, len(starts))
    if kind == "chec":
        for centre, height in _PEAKS:
            sig += height * np.exp(-((x - centre) ** 2) / (2 * 0.004 ** 2))
        sig[(x > 0.6) & (x < 0.63)] = np.nan
    return pd.DataFrame({"chrom": "chr2", "start": starts, "end": starts + bin_size,
                         "signal": sig.round(3), "intervalID": 1})


def example_peaks() -> pd.DataFrame:
    """Annotated regulatory elements: ``chrom``, ``start``, ``end``, ``name``, ``cre``, ``score``."""
    return pd.DataFrame({
        "chrom": "chr2",
        "start": [1_222_000, 1_261_000, 1_268_000, 1_270_500, 1_318_000, 1_325_000],
        "end": [1_228_000, 1_266_000, 1_273_000, 1_274_000, 1_324_000, 1_327_500],
        "name": ["P1", "E1", "E2", "C1", "P2", "E3"],
        "cre": ["promoter", "enhancer", "enhancer", "CTCF", "promoter", "enhancer"],
        "score": [880, 410, 350, 600, 950, 220],
    })


def example_motif_hits(n: int = 120) -> pd.DataFrame:
    """Motif hits: ``chrom``, ``start``, ``end``, ``motif``, ``score``, ``strand``.

    Hits are 12 bp, as from a PWM scan (e.g. ``pymisha.gscreen`` on a PWM
    virtual track, or moodspy positions).
    """
    region = example_region()
    rng = np.random.default_rng(2)
    pos = np.sort(rng.integers(region.start, region.end - 12, n))
    return pd.DataFrame({
        "chrom": "chr2",
        "start": pos,
        "end": pos + 12,
        "motif": rng.choice(["Msn2", "Gata2", "Ctcf"], n, p=[0.5, 0.3, 0.2]),
        "score": rng.uniform(8, 14, n).round(2),
        "strand": rng.choice(["+", "-"], n),
    })


def example_gtf_path() -> Path:
    """Path of the bundled synthetic GTF (UCSC style, 1-based) for the example locus."""
    return Path(__file__).resolve().parent / "data" / "example.gtf"


def example_genes() -> pd.DataFrame:
    """The bundled GTF read with :func:`~pygenomic_tracks.read_gtf` (exon and CDS rows)."""
    return read_gtf(example_gtf_path())
