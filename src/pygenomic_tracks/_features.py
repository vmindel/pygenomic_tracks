"""Shared machinery for feature tracks (intervals, genes): packing, labels, colours, strand."""
from __future__ import annotations

from typing import Dict, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize

from ._coords import Region
from ._pack import pack_rows

LABEL_PITCH = 1.8  # row pitch (in box-height units) when labels are drawn under features
PLAIN_PITCH = 1.0


def bp_per_px(ax, region: Region) -> float:
    return region.length / max(1.0, ax.get_window_extent().width)


def text_width_bp(ax, region: Region, labels: Sequence[str], fontsize: float) -> np.ndarray:
    # Approximate glyph width of 0.6 em; avoids a renderer round-trip per label.
    char_px = 0.6 * fontsize * ax.figure.dpi / 72
    n_chars = np.array([len(t) for t in labels], dtype=float)
    return n_chars * char_px * bp_per_px(ax, region)


def widen(starts, ends, min_bp: float):
    """Expand features narrower than ``min_bp`` symmetrically around their centre."""
    narrow = (ends - starts) < min_bp
    mid = (starts + ends) / 2
    return np.where(narrow, mid - min_bp / 2, starts), np.where(narrow, mid + min_bp / 2, ends)


def assign_rows(ax, region, starts, ends, labels, fontsize, pack, gap_px=4.0):
    """Row index per feature, accounting for label width when labels are drawn."""
    if not pack or len(starts) == 0:
        return np.zeros(len(starts), dtype=int)
    lo, hi = starts, ends
    if labels is not None:
        half = text_width_bp(ax, region, labels, fontsize) / 2
        mid = (starts + ends) / 2
        lo, hi = np.minimum(starts, mid - half), np.maximum(ends, mid + half)
    return pack_rows(lo, hi, gap=gap_px * bp_per_px(ax, region))


def set_row_limits(ax, n_rows: int, pitch: float, height: float, labelled: bool) -> None:
    bottom = (max(n_rows, 1) - 1) * pitch + height / 2 + (pitch - height if labelled else 0.2)
    ax.set_ylim(bottom, -height / 2 - 0.2)


def draw_labels(ax, xs, ys, labels, height, fontsize, color="0.15") -> None:
    for x, y, t in zip(xs, ys, labels):
        ax.text(x, y + height / 2 + 0.08, t, ha="center", va="top", fontsize=fontsize,
                color=color, clip_on=True)


def resolve_colors(
    values: Optional[pd.Series],
    n: int,
    color,
    palette,
    vmin: Optional[float],
    vmax: Optional[float],
) -> Tuple[list, Dict[str, object]]:
    """Per-feature colours plus a {category: colour} map for a legend (empty if continuous).

    - ``values`` None: every feature gets ``color``.
    - numeric ``values``: mapped through a colormap (``palette`` name, default viridis).
    - otherwise categorical: ``palette`` may be a dict, a list of colours or a colormap name.
    """
    if values is None:
        return [color] * n, {}

    if pd.api.types.is_numeric_dtype(values) and not pd.api.types.is_bool_dtype(values):
        v = values.to_numpy(dtype=float)
        norm = Normalize(
            vmin=np.nanmin(v) if vmin is None and n else vmin,
            vmax=np.nanmax(v) if vmax is None and n else vmax,
        )
        cmap = plt.get_cmap(palette or "viridis")
        return [cmap(norm(x)) for x in v], {}

    if isinstance(values.dtype, pd.CategoricalDtype):
        cats = [str(c) for c in values.cat.categories]
    else:
        cats = sorted(values.astype(str).unique())
    keys = values.astype(str).tolist()

    if isinstance(palette, dict):
        mapping = {c: palette.get(c, color) for c in cats}
    elif isinstance(palette, (list, tuple)):
        mapping = {c: palette[i % len(palette)] for i, c in enumerate(cats)}
    else:
        cmap = plt.get_cmap(palette or ("tab10" if len(cats) <= 10 else "tab20"))
        mapping = {c: cmap(i % cmap.N) for i, c in enumerate(cats)}
    return [mapping[k] for k in keys], mapping


def draw_chevrons(ax, region, starts, ends, ys, strands, color, spacing_px=40.0, size=3.5,
                  min_px=6.0, zorder=3) -> None:
    """Strand arrowheads along each feature; features narrower than ``min_px`` get none."""
    step = spacing_px * bp_per_px(ax, region)
    min_len = min_px * bp_per_px(ax, region)
    pts = {"+": ([], []), "-": ([], [])}
    for s, e, y, st in zip(starts, ends, ys, strands):
        if st not in pts or (e - s) < min_len:
            continue
        n = max(1, int((e - s) // step))
        xs = s + (np.arange(n) + 0.5) * (e - s) / n
        pts[st][0].extend(xs)
        pts[st][1].extend([y] * n)
    for st, marker in (("+", ">"), ("-", "<")):
        xs, ys_ = pts[st]
        if xs:
            ax.scatter(xs, ys_, marker=marker, s=size ** 2, c=[color], linewidths=0,
                       zorder=zorder, clip_on=True)
