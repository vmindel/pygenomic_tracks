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
    """Draw a feature track of boxes: peaks, CREs, motif hits, any BED-like table.

    Everything is drawn in one colour by default; colouring, labels and
    strand marks are opt-in.

    Parameters
    ----------
    df : pandas.DataFrame
        Intervals with ``chrom``, ``start``, ``end``; other columns are used
        only when named below. Rows outside ``region`` are ignored.
    ax : matplotlib.axes.Axes
        Axes to draw on.
    region : Region, str or tuple
        The window to draw, 0-based half-open.
    cols : tuple of str
        Names of the chrom, start and end columns.
    color : matplotlib colour
        Colour of every feature when ``color_by`` is not set, and of
        categories missing from a ``palette`` dict.
    color_by : str, optional
        Column to colour by. A numeric column (e.g. motif score) is mapped
        through a colormap; any other column is treated as categories.
    palette : str, list or dict, optional
        Numeric ``color_by``: a colormap name (default ``"viridis"``).
        Categorical: a ``{category: colour}`` dict (keeps colours stable
        across figures), a list of colours, or a colormap name (default
        ``tab10``/``tab20``).
    vmin, vmax : float, optional
        Colormap limits for a numeric ``color_by``; default to the data range
        in the window.
    label_by : str, optional
        Column printed under each feature. Packing leaves room for it.
    strand : str, optional
        Strand column (``"+"``/``"-"``). Draws arrowheads on features at
        least 6 px wide; narrower features get none.
    pack : bool
        Stack overlapping features (and their labels) into rows. ``False``
        draws everything on one row.
    height : float
        Box height as a fraction of a row.
    min_width_px : float
        Features narrower than this many pixels are widened around their
        centre, so 10 bp motifs stay visible in a 1 Mb window.
    edgecolor : matplotlib colour
        Box outline; ``"none"`` by default.
    fontsize : float
        Size of feature labels and legend text.
    legend : bool
        Add a legend of the categories (categorical ``color_by`` only).
    label : str, optional
        Track label printed to the left of the axes.

    Returns
    -------
    matplotlib.axes.Axes
        ``ax``, for chaining.
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
