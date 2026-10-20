"""Quantitative signal tracks (coverage, ChEC/ChIP-seq, bigwig-derived)."""
from __future__ import annotations

from typing import Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd
from matplotlib.patches import StepPatch
from matplotlib.ticker import MaxNLocator

from ._coords import DEFAULT_COLS, Region, RegionLike, as_region, check_columns, clipped, select
from .layout import finish_track, hide_y

_AGGS = ("max", "min", "mean")


def _to_stairs(starts, ends, values):
    """Edges/values for ``ax.stairs``, with NaN segments filling gaps between intervals."""
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


def _bin_signal(starts, ends, values, region: Region, n_bins: int, agg: str):
    """Aggregate intervals onto ``n_bins`` equal bins spanning the region.

    ``mean`` is weighted by overlap length; bins with no data are NaN.
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


def signal(
    df: pd.DataFrame,
    ax,
    region: RegionLike,
    value: str = "signal",
    cols: Sequence[str] = DEFAULT_COLS,
    n_bins: Union[int, str, None] = "auto",
    agg: str = "max",
    color: str = "C0",
    fill: bool = True,
    alpha: float = 1.0,
    linewidth: float = 0.8,
    ylim: Union[float, Tuple[float, float], None] = None,
    range_label: bool = True,
    name: Optional[str] = None,
    label: Optional[str] = None,
):
    """Draw a quantitative signal track.

    Each row is drawn as a step over its real extent, so variable-width bins
    are correct and gaps between rows, or NaN values, stay empty instead of
    being interpolated across.

    Parameters
    ----------
    df : pandas.DataFrame
        Intervals with ``chrom``, ``start``, ``end`` and a numeric ``value``
        column, e.g. ``pymisha.gextract`` output. Rows outside ``region`` are
        ignored.
    ax : matplotlib.axes.Axes
        Axes to draw on.
    region : Region, str or tuple
        ``"chr2:1,200,000-1,350,000"``, ``("chr2", 1200000, 1350000)`` or a
        :class:`Region`. 0-based, half-open.
    value : str
        Column holding the signal.
    cols : tuple of str
        Names of the chrom, start and end columns.
    n_bins : "auto", int or None
        ``"auto"`` aggregates to about one bin per pixel of ``ax`` when there
        are more rows than pixels, or when rows overlap. An int fixes the bin
        count. ``None`` draws every row as is and raises on overlapping rows.
    agg : {"max", "min", "mean"}
        Aggregation used when binning. ``"max"`` keeps narrow peaks visible;
        ``"mean"`` is weighted by overlap length.
    color, alpha : matplotlib colour and float
        Fill (or line) colour and opacity.
    fill : bool
        Filled area (default) or outline only.
    linewidth : float
        Line width when ``fill=False``.
    ylim : float or (float, float), optional
        A single ``ymax`` or ``(ymin, ymax)``. Pass the same value to several
        tracks to put them on one scale. By default the range is
        ``[min(0, data), max(0, data)]`` over every signal already drawn on
        ``ax``, so overlays share one scale and are not clipped.
    range_label : bool
        Hide the y-axis and print an IGV-style ``[min - max]`` label instead.
    name : str, optional
        Legend entry for this series; adds a legend to ``ax``. Useful when
        overlaying several signals on one axes.
    label : str, optional
        Track label printed to the left of the axes.

    Returns
    -------
    matplotlib.axes.Axes
        ``ax``, for chaining.

    Raises
    ------
    ValueError
        If ``agg`` is unknown, or rows overlap with ``n_bins=None``.
    KeyError
        If a required column is missing.
    """
    if agg not in _AGGS:
        raise ValueError(f"agg must be one of {_AGGS}, got {agg!r}")
    region = as_region(region)
    check_columns(df, [value])
    sub = select(df, region, cols)
    starts, ends = clipped(sub, region, cols)
    values = sub[value].to_numpy(dtype=float)

    if n_bins == "auto":
        n_bins = max(1, int(ax.get_window_extent().width))
    overlapping = len(starts) > 1 and bool(np.any(starts[1:] < ends[:-1]))

    if len(starts) == 0:
        edges, heights = np.array([]), np.array([])
    elif n_bins is not None and (len(starts) > n_bins or overlapping):
        edges, heights = _bin_signal(starts, ends, values, region, int(n_bins), agg)
    elif overlapping:
        raise ValueError(
            "signal intervals overlap; pass n_bins to aggregate them, or merge them first"
        )
    else:
        edges, heights = _to_stairs(starts, ends, values)

    if len(heights) and np.isfinite(heights).any():
        ax.stairs(
            heights, edges, baseline=0, fill=fill, color=color, alpha=alpha,
            linewidth=0 if fill else linewidth, label=name,
        )
        if ylim is None:
            # Range over every step series on the axes, so overlays share one scale.
            drawn = np.concatenate([p.get_data().values for p in ax.patches
                                    if isinstance(p, StepPatch)])
            lo, hi = min(0.0, np.nanmin(drawn)), max(0.0, np.nanmax(drawn))
        elif np.isscalar(ylim):
            lo, hi = 0.0, float(ylim)
        else:
            lo, hi = ylim
        if hi == lo:
            hi = lo + 1.0
        ax.set_ylim(lo, hi)

    finish_track(ax, region, label)
    if range_label:
        lo, hi = ax.get_ylim()
        hide_y(ax)
        text = f"[{lo:.3g} - {hi:.3g}]"
        existing = [t for t in ax.texts if t.get_gid() == "pygt-range"]
        if existing:
            existing[0].set_text(text)
        else:
            ax.text(0.005, 0.98, text, transform=ax.transAxes, ha="left", va="top",
                    fontsize=7, gid="pygt-range")
    else:
        ax.yaxis.set_major_locator(MaxNLocator(nbins=2))
    if name is not None:
        ax.legend(loc="upper right", frameon=False, fontsize=7)
    return ax
