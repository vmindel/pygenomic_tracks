"""Quantitative signal tracks (coverage, ChEC/ChIP-seq, bigwig-derived)."""
from __future__ import annotations

from typing import Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd
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
    label: Optional[str] = None,
):
    """Draw a signal track from interval rows with a numeric ``value`` column.

    Intervals are drawn as steps over their true extent; gaps between
    intervals are left empty rather than interpolated across.

    n_bins:
        "auto" aggregates to roughly one bin per pixel of the axes when there
        are more intervals than pixels. An int sets the bin count, None draws
        every interval as is (overlapping intervals then raise).
    agg:
        "max" (default, keeps narrow peaks visible), "min" or "mean"
        (overlap-weighted).
    ylim:
        (ymin, ymax), or a single ymax, to put several tracks on one scale.
    range_label:
        Replace the y-axis by an IGV-style "[min - max]" label.
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
            linewidth=0 if fill else linewidth,
        )
        if ylim is None:
            lo, hi = min(0.0, np.nanmin(heights)), max(0.0, np.nanmax(heights))
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
        ax.text(0.005, 0.98, f"[{lo:.3g} - {hi:.3g}]", transform=ax.transAxes,
                ha="left", va="top", fontsize=7)
    else:
        ax.yaxis.set_major_locator(MaxNLocator(nbins=2))
    return ax
