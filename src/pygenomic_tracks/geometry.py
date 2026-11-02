"""What to draw, independent of how it is drawn.

These functions turn tables and a region into plain arrays and DataFrames
of steps, boxes and rows. They take no matplotlib ``Axes``: decisions that
depend on the screen (bin count, how much room a label needs) take the
resolution as an argument (``n_bins``, ``bp_per_px``).

The matplotlib tracks in this package are built on them, and so is
plotnine-genomic-tracks. This module is public API: changes to it are
breaking changes for other packages.
"""
from __future__ import annotations

from typing import Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from ._coords import (
    DEFAULT_COLS,
    Region,
    RegionLike,
    as_region,
    check_columns,
    chrom_mask,
    clipped,
    select,
)
from ._pack import pack_rows

__all__ = [
    "AGGREGATIONS",
    "Region",
    "as_region",
    "bin_signal",
    "clipped",
    "feature_rows",
    "gene_models",
    "genomic_unit",
    "pack_rows",
    "select",
    "signal_steps",
    "text_width_bp",
    "to_stairs",
    "widen",
]

AGGREGATIONS = ("max", "min", "mean")


# Signal


def to_stairs(starts, ends, values):
    """Edges and values for ``matplotlib.axes.Axes.stairs``.

    A NaN segment is inserted for every gap between consecutive intervals,
    so gaps stay empty instead of being bridged.

    Parameters
    ----------
    starts, ends, values : numpy.ndarray
        Non-overlapping intervals sorted by start.

    Returns
    -------
    edges : numpy.ndarray
        ``len(values) + 1`` edges.
    values : numpy.ndarray
    """
    n = len(starts)
    gap = starts[1:] > ends[:-1]
    n_gaps_before = np.concatenate([[0], np.cumsum(gap)])
    pos = np.arange(n) + n_gaps_before
    total = n + int(n_gaps_before[-1])
    edges = np.empty(total + 1)
    out = np.full(total, np.nan)
    out[pos] = values
    edges[pos] = starts
    after_gap = np.nonzero(gap)[0] + 1
    edges[pos[after_gap] - 1] = ends[after_gap - 1]
    edges[-1] = ends[-1]
    return edges, out


def bin_signal(starts, ends, values, region: Region, n_bins: int, agg: str):
    """Aggregate intervals onto ``n_bins`` equal bins spanning the region.

    Parameters
    ----------
    starts, ends, values : numpy.ndarray
        Intervals, clipped to ``region``; they may overlap.
    region : Region
    n_bins : int
    agg : {"max", "min", "mean"}
        ``"mean"`` is weighted by overlap length.

    Returns
    -------
    edges : numpy.ndarray
        ``n_bins + 1`` bin edges.
    values : numpy.ndarray
        One value per bin; NaN where no interval overlaps.
    """
    edges = np.linspace(region.start, region.end, n_bins + 1)
    width = region.length / n_bins
    b0 = np.clip(np.floor((starts - region.start) / width).astype(int), 0, n_bins - 1)
    b1 = np.clip(np.ceil((ends - region.start) / width).astype(int) - 1, 0, n_bins - 1)
    b1 = np.maximum(b0, b1)
    counts = b1 - b0 + 1
    rep = np.repeat(np.arange(len(starts)), counts)
    offset = np.arange(counts.sum()) - np.repeat(np.cumsum(counts) - counts, counts)
    bins = b0[rep] + offset
    v = values[rep]
    keep = ~np.isnan(v)
    bins, v, rep = bins[keep], v[keep], rep[keep]

    if agg == "max":
        out = np.full(n_bins, -np.inf)
        np.maximum.at(out, bins, v)
        out[np.isneginf(out)] = np.nan
    elif agg == "min":
        out = np.full(n_bins, np.inf)
        np.minimum.at(out, bins, v)
        out[np.isposinf(out)] = np.nan
    else:
        overlap = np.minimum(ends[rep], edges[bins + 1]) - np.maximum(starts[rep], edges[bins])
        wsum = np.bincount(bins, weights=overlap, minlength=n_bins)
        vsum = np.bincount(bins, weights=v * overlap, minlength=n_bins)
        with np.errstate(invalid="ignore", divide="ignore"):
            out = np.where(wsum > 0, vsum / wsum, np.nan)
    return edges, out


def signal_steps(
    df: pd.DataFrame,
    region: RegionLike,
    value: str = "signal",
    cols: Sequence[str] = DEFAULT_COLS,
    n_bins: Optional[int] = None,
    agg: str = "max",
) -> pd.DataFrame:
    """Steps of a signal track within a region.

    Rows are clipped to the region. They are aggregated onto ``n_bins``
    equal bins when there are more rows than bins, or when rows overlap;
    otherwise every row is one step at its real extent.

    Parameters
    ----------
    df : pandas.DataFrame
        Intervals with the ``cols`` columns and a numeric ``value`` column.
    region : Region, str or tuple
    value : str
        Column holding the signal.
    cols : tuple of str
        Names of the chrom, start and end columns.
    n_bins : int, optional
        Resolution, typically the panel width in pixels. ``None`` never
        aggregates and raises on overlapping rows.
    agg : {"max", "min", "mean"}
        Aggregation used when binning.

    Returns
    -------
    pandas.DataFrame
        ``start``, ``end``, ``value``, sorted by start, one row per step.
        Gaps between steps, empty bins and NaN values are left out.

    Raises
    ------
    ValueError
        If ``agg`` is unknown, or rows overlap with ``n_bins=None``.
    """
    if agg not in AGGREGATIONS:
        raise ValueError(f"agg must be one of {AGGREGATIONS}, got {agg!r}")
    region = as_region(region)
    check_columns(df, [value])
    sub = select(df, region, cols)
    starts, ends = clipped(sub, region, cols)
    values = sub[value].to_numpy(dtype=float)
    overlapping = len(starts) > 1 and bool(np.any(starts[1:] < ends[:-1]))

    if n_bins is not None and len(starts) and (len(starts) > n_bins or overlapping):
        edges, values = bin_signal(starts, ends, values, region, int(n_bins), agg)
        starts, ends = edges[:-1], edges[1:]
    elif overlapping:
        raise ValueError(
            "signal intervals overlap; pass n_bins to aggregate them, or merge them first"
        )

    out = pd.DataFrame({"start": starts, "end": ends, "value": values})
    return out[out["value"].notna()].reset_index(drop=True)


# Features


def widen(starts, ends, min_bp: float):
    """Expand features narrower than ``min_bp`` symmetrically around their centre."""
    narrow = (ends - starts) < min_bp
    mid = (starts + ends) / 2
    return np.where(narrow, mid - min_bp / 2, starts), np.where(narrow, mid + min_bp / 2, ends)


def text_width_bp(
    labels: Sequence[str], fontsize: float, bp_per_px: float, dpi: float
) -> np.ndarray:
    """Approximate width of each label in bp.

    Glyphs are taken as 0.6 em wide, which avoids measuring text with a
    renderer.
    """
    char_px = 0.6 * fontsize * dpi / 72
    n_chars = np.array([len(t) for t in labels], dtype=float)
    return n_chars * char_px * bp_per_px


def feature_rows(
    starts,
    ends,
    bp_per_px: float,
    labels: Optional[Sequence[str]] = None,
    fontsize: float = 7,
    dpi: float = 100,
    pack: bool = True,
    gap_px: float = 4.0,
) -> np.ndarray:
    """Row index per feature so that features, and their labels, do not overlap.

    Parameters
    ----------
    starts, ends : array-like
    bp_per_px : float
        Resolution of the drawing; sets how much room a label and the gap
        between features need.
    labels : sequence of str, optional
        Labels centred under each feature; packing leaves room for them.
    fontsize : float
        Label size in points.
    dpi : float
        Pixels per inch of the drawing.
    pack : bool
        ``False`` puts everything on row 0.
    gap_px : float
        Minimum horizontal space between neighbours in one row.

    Returns
    -------
    numpy.ndarray of int
    """
    starts = np.asarray(starts, dtype=float)
    ends = np.asarray(ends, dtype=float)
    if not pack or len(starts) == 0:
        return np.zeros(len(starts), dtype=int)
    lo, hi = starts, ends
    if labels is not None:
        half = text_width_bp(labels, fontsize, bp_per_px, dpi) / 2
        mid = (starts + ends) / 2
        lo, hi = np.minimum(starts, mid - half), np.maximum(ends, mid + half)
    return pack_rows(lo, hi, gap=gap_px * bp_per_px)


def gene_models(
    df: pd.DataFrame,
    region: RegionLike,
    cols: Sequence[str] = DEFAULT_COLS,
    feature: str = "feature",
    transcript: str = "transcript_id",
    name: Optional[str] = "gene_name",
    strand: Optional[str] = "strand",
    transcripts: str = "all",
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Transcripts and their exon/CDS parts within a region.

    Parameters
    ----------
    df : pandas.DataFrame
        One row per exon/CDS, as :func:`pygenomic_tracks.read_gtf` or
        :func:`pygenomic_tracks.annotation` return. Rows with another
        ``feature`` are ignored.
    region : Region, str or tuple
    cols, feature, transcript, name, strand :
        Column names, as in :func:`pygenomic_tracks.genes`.
    transcripts : {"all", "longest"}
        ``"longest"`` keeps one transcript per gene name (longest span).

    Returns
    -------
    spans : pandas.DataFrame
        One row per transcript overlapping the region, sorted by start,
        indexed by transcript id: ``start`` and ``end`` (unclipped),
        ``name``, and ``strand`` when ``strand`` is given.
    parts : pandas.DataFrame
        Exon and CDS rows of those transcripts: ``transcript``, ``start``,
        ``end``, ``feature`` and ``full``. ``full`` is True for CDS and
        for every part of a transcript without CDS; other parts are UTRs.
    """
    region = as_region(region)
    c, s, e = cols
    needed = list(cols) + [feature, transcript]
    needed += [x for x in (name, strand) if x is not None]
    check_columns(df, needed)
    if transcripts not in ("all", "longest"):
        raise ValueError(f"transcripts must be 'all' or 'longest', got {transcripts!r}")

    parts = df.loc[chrom_mask(df[c], region.chrom) & df[feature].isin(["exon", "CDS"])]
    span = parts.groupby(transcript, sort=False).agg(start=(s, "min"), end=(e, "max"))
    span = span[(span["end"] > region.start) & (span["start"] < region.end)]
    if name is not None:
        span["name"] = parts.groupby(transcript, sort=False)[name].first().loc[span.index]
    else:
        span["name"] = span.index.astype(str)
    if transcripts == "longest" and len(span):
        span = span.assign(_len=span["end"] - span["start"])
        span = span.sort_values("_len", ascending=False).groupby("name", sort=False).head(1)
        span = span.drop(columns="_len")
    span = span.sort_values("start", kind="stable")
    if strand is not None:
        span["strand"] = parts.groupby(transcript, sort=False)[strand].first().loc[span.index]

    parts = parts[parts[transcript].isin(span.index)]
    coding = set(parts.loc[parts[feature] == "CDS", transcript])
    out = pd.DataFrame(
        {
            "transcript": parts[transcript].to_numpy(),
            "start": parts[s].to_numpy(),
            "end": parts[e].to_numpy(),
            "feature": parts[feature].to_numpy(),
        }
    )
    out["full"] = (out["feature"] == "CDS") | ~out["transcript"].isin(coding)
    return span, out


# Axis


def genomic_unit(length: int) -> Tuple[str, float]:
    """Unit for labelling a window of ``length`` bp: ``("Mb", 1e6)``, ``("kb", 1e3)`` or ``("bp", 1.0)``."""
    if length >= 2_000_000:
        return "Mb", 1e6
    if length >= 2_000:
        return "kb", 1e3
    return "bp", 1.0
