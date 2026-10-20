"""Generic interval tracks: peaks, CREs, motif hits, any BED-like table."""
from __future__ import annotations

from typing import Optional, Sequence

import numpy as np
import pandas as pd
from matplotlib.collections import PatchCollection
from matplotlib.patches import Patch, Rectangle

from . import _features as F
from ._coords import DEFAULT_COLS, RegionLike, as_region, check_columns, clipped, select
from .layout import finish_track, hide_y


def intervals(
    df: pd.DataFrame,
    ax,
    region: RegionLike,
    cols: Sequence[str] = DEFAULT_COLS,
    color="0.35",
    color_by: Optional[str] = None,
    palette=None,
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    label_by: Optional[str] = None,
    strand: Optional[str] = None,
    pack: bool = True,
    height: float = 0.6,
    min_width_px: float = 1.0,
    edgecolor="none",
    fontsize: float = 7,
    legend: bool = False,
    label: Optional[str] = None,
):
    """Draw intervals as boxes, one colour by default.

    color_by:
        Optional column to colour by. Numeric columns (e.g. motif score) go
        through a colormap (``palette`` = colormap name, ``vmin``/``vmax``);
        anything else is treated as categories (``palette`` = dict, list of
        colours or colormap name).
    label_by:
        Optional column whose values are printed under each feature.
    strand:
        Optional strand column ("+"/"-"); draws arrowheads on features wide
        enough to hold them.
    pack:
        Stack overlapping features into rows; False draws everything on one row.
    min_width_px:
        Features narrower than this many pixels are widened so that short
        features (motif hits in a 1 Mb window) stay visible.
    """
    region = as_region(region)
    check_columns(df, [c for c in (color_by, label_by, strand) if c is not None])
    sub = select(df, region, cols)
    starts, ends = clipped(sub, region, cols)
    starts, ends = F.widen(starts, ends, min_width_px * F.bp_per_px(ax, region))
    n = len(sub)

    labels = sub[label_by].astype(str).tolist() if label_by is not None else None
    rows = F.assign_rows(ax, region, starts, ends, labels, fontsize, pack)
    pitch = F.LABEL_PITCH if labels is not None else F.PLAIN_PITCH
    ys = rows * pitch

    colors, mapping = F.resolve_colors(
        sub[color_by] if color_by is not None else None, n, color, palette, vmin, vmax
    )
    if n:
        boxes = [Rectangle((s, y - height / 2), e - s, height)
                 for s, e, y in zip(starts, ends, ys)]
        ax.add_collection(PatchCollection(boxes, facecolors=colors, edgecolors=edgecolor,
                                          linewidths=0.5, zorder=2))
    if strand is not None:
        F.draw_chevrons(ax, region, starts, ends, ys, sub[strand].astype(str), "white")
    if labels is not None:
        F.draw_labels(ax, (starts + ends) / 2, ys, labels, height, fontsize)

    F.set_row_limits(ax, int(rows.max()) + 1 if n else 1, pitch, height, labels is not None)
    finish_track(ax, region, label)
    hide_y(ax)
    if legend and mapping:
        ax.legend(handles=[Patch(color=c, label=k) for k, c in mapping.items()],
                  loc="upper left", bbox_to_anchor=(1.0, 1.0), frameon=False,
                  fontsize=fontsize)
    return ax
