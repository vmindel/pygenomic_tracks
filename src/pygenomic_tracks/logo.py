"""Sequence logo tracks: letters scaled by attribution, ISM or PWM values."""
from __future__ import annotations

from functools import lru_cache
from typing import Mapping, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd
from matplotlib.font_manager import FontProperties
from matplotlib.patches import PathPatch
from matplotlib.path import Path
from matplotlib.textpath import TextPath
from matplotlib.ticker import MaxNLocator

from ._coords import DEFAULT_COLS, RegionLike, as_region
from .geometry import bin_signal, logo_stacks, position_totals, to_stairs
from .layout import finish_track, hide_y

DNA = ("A", "C", "G", "T")

BASE_COLORS = {"A": "#109648", "C": "#255C99", "G": "#F7B32B", "T": "#D62839", "U": "#D62839"}

_DEFAULT_FONT = FontProperties(family="DejaVu Sans", weight="bold")

# Horizontal space left between neighbouring letters, as a fraction of their width.
_GAP = 0.04


@lru_cache(maxsize=256)
def _unit_glyph(letter: str, font: FontProperties) -> Tuple[np.ndarray, np.ndarray]:
    """Outline of ``letter`` scaled so its ink fills the unit square."""
    path = TextPath((0, 0), letter, size=1, prop=font)
    verts = path.vertices.copy()
    lo, hi = verts.min(axis=0), verts.max(axis=0)
    verts = (verts - lo) / np.where(hi > lo, hi - lo, 1.0)
    codes = path.codes if path.codes is not None else np.full(len(verts), Path.LINETO)
    return verts, np.asarray(codes)


def _glyphs_path(letter, font, x0, x1, y0, y1, flip) -> Path:
    """One compound path holding ``letter`` stretched over each box."""
    verts, codes = _unit_glyph(letter, font)
    n = len(verts)
    ux = np.tile(verts[:, 0], len(x0))
    uy = np.tile(verts[:, 1], len(x0))
    uy = np.where(np.repeat(flip, n), 1.0 - uy, uy)
    x = np.repeat(x0, n) + ux * np.repeat(x1 - x0, n)
    y = np.repeat(y0, n) + uy * np.repeat(y1 - y0, n)
    return Path(np.column_stack([x, y]), np.tile(codes, len(x0)))


def logo_table(
    values,
    chrom: str,
    start: int,
    sequence: Optional[str] = None,
    bases: Sequence[str] = DNA,
) -> pd.DataFrame:
    """Turn a model's per-base output into a table for :func:`logo`.

    Parameters
    ----------
    values : array-like
        ``(L, len(bases))`` scores per position and base, e.g. hypothetical
        contributions or ISM effects. A ``(len(bases), L)`` array, the usual
        channels-first layout of a PyTorch model, is transposed. A
        one-dimensional ``(L,)`` array is one score per position and needs
        ``sequence``.
    chrom : str
    start : int
        0-based coordinate of the first position.
    sequence : str, optional
        Reference sequence of length ``L``. With a 2-D ``values`` the scores
        are multiplied by its one-hot encoding, which keeps only the
        reference base at each position (hypothetical to actual
        contributions). With a 1-D ``values`` each score is put on the
        reference base. Positions with other characters (``N``) get zeros.
        Case is ignored.
    bases : sequence of str
        Letters, in the column order of ``values``.

    Returns
    -------
    pandas.DataFrame
        ``chrom``, ``start``, ``end`` (1 bp rows) and one column per base.

    Raises
    ------
    ValueError
        If the shapes of ``values`` and ``sequence`` do not match.
    """
    bases = list(bases)
    arr = np.asarray(values, dtype=float)
    if arr.ndim == 2 and arr.shape[1] != len(bases) and arr.shape[0] == len(bases):
        arr = arr.T
    onehot = None
    if sequence is not None:
        seq = np.array(list(sequence.upper()))
        onehot = (seq[:, None] == np.array(bases)[None, :]).astype(float)
    if arr.ndim == 1:
        if onehot is None:
            raise ValueError("one score per position needs sequence= to place it on a base")
        if len(arr) != len(onehot):
            raise ValueError(f"values has {len(arr)} positions, sequence {len(onehot)}")
        arr = arr[:, None] * onehot
    elif arr.ndim == 2 and arr.shape[1] == len(bases):
        if onehot is not None:
            if len(onehot) != len(arr):
                raise ValueError(f"values has {len(arr)} positions, sequence {len(onehot)}")
            arr = arr * onehot
    else:
        raise ValueError(
            f"values must be (L,), (L, {len(bases)}) or ({len(bases)}, L); got {arr.shape}"
        )
    pos = int(start) + np.arange(len(arr))
    out = pd.DataFrame({"chrom": chrom, "start": pos, "end": pos + 1})
    for j, b in enumerate(bases):
        out[b] = arr[:, j]
    return out


def logo(
    df: pd.DataFrame,
    ax,
    region: RegionLike,
    bases: Sequence[str] = DNA,
    cols: Sequence[str] = DEFAULT_COLS,
    colors: Optional[Mapping[str, str]] = None,
    color: str = "0.35",
    letters: Union[bool, str] = "auto",
    min_letter_px: float = 2.0,
    flip_below: bool = False,
    font: Optional[FontProperties] = None,
    baseline: bool = True,
    ylim: Union[float, Tuple[float, float], None] = None,
    range_label: bool = True,
    label: Optional[str] = None,
):
    """Draw a sequence logo: letters scaled by attribution, ISM or PWM values.

    At each position the positive values are stacked upward from zero and
    the negative values downward, with the largest magnitude furthest from
    zero. Actual contributions (one non-zero base per position) give one
    letter per position; hypothetical contributions, ISM or PWMs give a
    stack.

    Parameters
    ----------
    df : pandas.DataFrame
        Rows with ``chrom``, ``start``, ``end`` (normally 1 bp each) and one
        numeric column per letter in ``bases``. :func:`logo_table` builds
        this from a model's array. Rows outside ``region`` are ignored and
        NaN counts as zero.
    ax : matplotlib.axes.Axes
    region : Region, str or tuple
    bases : sequence of str
        Columns to draw; each column name is the letter drawn. Any letters
        work, e.g. the 20 amino acids for a protein attribution.
    cols : tuple of str
        Names of the chrom, start and end columns.
    colors : dict, optional
        ``{letter: colour}``, merged over the defaults. For nucleotide
        ``bases`` the defaults are A green, C blue, G orange, T and U red;
        other alphabets have none, so every letter not named is ``color``.
    color : matplotlib colour
        Letters missing from ``colors``, and the bars of a zoomed-out track.
    letters : "auto", True or False
        ``"auto"`` draws letters when one position spans at least
        ``min_letter_px`` pixels of ``ax`` and otherwise the summed positive
        and negative values per position as bars, aggregated to screen
        resolution (``max`` of the positive sums, ``min`` of the negative
        ones). ``True`` and ``False`` force one or the other.
    min_letter_px : float
        Narrowest letter drawn under ``letters="auto"``.
    flip_below : bool
        Draw letters below zero upside down.
    font : matplotlib.font_manager.FontProperties, optional
        Letter font; DejaVu Sans Bold by default. Each glyph is stretched to
        fill its box, so widths do not depend on the font's metrics.
    baseline : bool
        Draw a thin line at zero when the range spans it.
    ylim : float or (float, float), optional
        A single ``ymax`` or ``(ymin, ymax)``; pass the same value to several
        logos to compare them. By default ``[min(0, data), max(0, data)]``.
    range_label : bool
        Hide the y-axis and print an IGV-style ``[min - max]`` label instead.
    label : str, optional
        Track label printed to the left of the axes.

    Returns
    -------
    matplotlib.axes.Axes
        ``ax``, for chaining.

    Raises
    ------
    KeyError
        If a column in ``cols`` or ``bases`` is missing.
    ValueError
        If ``letters`` is not ``"auto"``, ``True`` or ``False``.
    """
    if letters not in ("auto", True, False):
        raise ValueError(f"letters must be 'auto', True or False, got {letters!r}")
    region = as_region(region)
    nucleotide = set(bases) <= set("ACGTUN")
    colors = {**(BASE_COLORS if nucleotide else {}), **(colors or {})}
    font = font or _DEFAULT_FONT
    px = max(1.0, ax.get_window_extent().width)
    if letters == "auto":
        letters = px / region.length >= min_letter_px

    lo = hi = 0.0
    if letters:
        stacks = logo_stacks(df, region, bases=bases, cols=cols)
        if len(stacks):
            width = (stacks["end"] - stacks["start"]).to_numpy()
            stacks = stacks.assign(x0=stacks["start"] + width * _GAP / 2,
                                   x1=stacks["end"] - width * _GAP / 2)
            for base, g in stacks.groupby("base", sort=False):
                y0 = np.minimum(g["bottom"], g["top"]).to_numpy()
                y1 = np.maximum(g["bottom"], g["top"]).to_numpy()
                flip = (g["top"] < 0).to_numpy() & flip_below
                path = _glyphs_path(base, font, g["x0"].to_numpy(), g["x1"].to_numpy(),
                                    y0, y1, flip)
                ax.add_patch(PathPatch(path, facecolor=colors.get(base, color),
                                       edgecolor="none", gid=f"pygt-logo-{base}"))
            lo = min(0.0, stacks[["bottom", "top"]].min().min())
            hi = max(0.0, stacks[["bottom", "top"]].max().max())
    else:
        totals = position_totals(df, region, bases=bases, cols=cols)
        if len(totals):
            s, e = totals["start"].to_numpy(float), totals["end"].to_numpy(float)
            n_bins = int(px)
            for part, agg in (("pos", "max"), ("neg", "min")):
                v = totals[part].to_numpy(float)
                if len(s) > n_bins or (len(s) > 1 and np.any(s[1:] < e[:-1])):
                    edges, v = bin_signal(s, e, v, region, n_bins, agg)
                    ss, ee = edges[:-1], edges[1:]
                else:
                    ss, ee = s, e
                keep = np.isfinite(v) & (v != 0)
                if keep.any():
                    edg, h = to_stairs(ss[keep], ee[keep], v[keep])
                    ax.stairs(h, edg, baseline=0, fill=True, color=color, linewidth=0)
            lo = min(0.0, totals["neg"].min())
            hi = max(0.0, totals["pos"].max())

    if ylim is not None:
        lo, hi = (0.0, float(ylim)) if np.isscalar(ylim) else ylim
    if hi == lo:
        hi = lo + 1.0
    ax.set_ylim(lo, hi)
    if baseline and lo < 0 < hi:
        ax.axhline(0, color="0.6", linewidth=0.5, zorder=0)

    finish_track(ax, region, label)
    if range_label:
        hide_y(ax)
        ax.text(0.005, 0.98, f"[{lo:.3g} - {hi:.3g}]", transform=ax.transAxes,
                ha="left", va="top", fontsize=7, gid="pygt-range")
    else:
        ax.yaxis.set_major_locator(MaxNLocator(nbins=2))
    return ax
