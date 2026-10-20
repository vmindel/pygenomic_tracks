"""Gene model tracks from long-format (GTF-like) tables."""
from __future__ import annotations

import os
from typing import Optional, Sequence, Union

import numpy as np
import pandas as pd
from matplotlib.collections import LineCollection, PatchCollection
from matplotlib.patches import Rectangle

from . import _features as F
from .annotations import annotation
from ._coords import DEFAULT_COLS, RegionLike, as_region, check_columns, chrom_mask
from .layout import finish_track, hide_y


def genes(
    df: Union[pd.DataFrame, str, os.PathLike],
    ax,
    region: RegionLike,
    cols: Sequence[str] = DEFAULT_COLS,
    feature: str = "feature",
    transcript: str = "transcript_id",
    name: Optional[str] = "gene_name",
    strand: Optional[str] = "strand",
    transcripts: str = "all",
    color="0.2",
    height: float = 0.6,
    utr_ratio: float = 0.5,
    show_names: bool = True,
    pack: bool = True,
    fontsize: float = 7,
    label: Optional[str] = None,
):
    """Draw transcript models: intron line, exons, thicker CDS, strand arrowheads, names.

    Parameters
    ----------
    df : pandas.DataFrame, str or path
        A UCSC assembly name (``"hg38"``; downloaded and cached on first use,
        see :func:`annotation`), a GTF path, or a DataFrame with one row per
        exon/CDS such as :func:`read_gtf` returns.
    ax : matplotlib.axes.Axes
        Axes to draw on.
    region : Region, str or tuple
        The window to draw. Transcripts crossing its edges are drawn clipped.
    cols : tuple of str
        Names of the chrom, start and end columns.
    feature : str
        Column with the feature type. Rows other than ``"exon"`` and
        ``"CDS"`` are ignored.
    transcript : str
        Column grouping rows into transcripts.
    name : str or None
        Column with the gene name shown under each transcript; ``None``
        shows transcript ids.
    strand : str or None
        Strand column for arrowheads; ``None`` draws none.
    transcripts : {"all", "longest"}
        ``"longest"`` keeps one transcript per gene name (longest span).
    color : matplotlib colour
        Colour of the whole model.
    height : float
        Height of coding exons as a fraction of a row.
    utr_ratio : float
        Height of non-coding parts of coding transcripts relative to
        ``height``. Transcripts without CDS rows are drawn at full height.
    show_names : bool
        Print names under transcripts.
    pack : bool
        Stack overlapping transcripts into rows.
    fontsize : float
        Size of the names.
    label : str, optional
        Track label printed to the left of the axes.

    Returns
    -------
    matplotlib.axes.Axes
        ``ax``, for chaining.
    """
    df = annotation(df)
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
    span = span.sort_values("start", kind="stable")
    if strand is not None:
        span["strand"] = parts.groupby(transcript, sort=False)[strand].first().loc[span.index]

    t_start = np.clip(span["start"].to_numpy(float), region.start, region.end)
    t_end = np.clip(span["end"].to_numpy(float), region.start, region.end)
    labels = span["name"].astype(str).tolist() if show_names else None
    rows = F.assign_rows(ax, region, t_start, t_end, labels, fontsize, pack)
    pitch = F.LABEL_PITCH if show_names else F.PLAIN_PITCH
    y_of = dict(zip(span.index, rows * pitch))

    if len(span):
        ax.add_collection(LineCollection(
            [[(a, y), (b, y)] for a, b, y in zip(t_start, t_end, rows * pitch)],
            colors=color, linewidths=0.8, zorder=1,
        ))
        if strand is not None:
            F.draw_chevrons(ax, region, t_start, t_end, rows * pitch,
                            span["strand"].astype(str), color, size=3.0, zorder=1)

        parts = parts[parts[transcript].isin(span.index)]
        coding = set(parts.loc[parts[feature] == "CDS", transcript])
        boxes = []
        for tid, fs, fe, ft in zip(parts[transcript], parts[s], parts[e], parts[feature]):
            full = ft == "CDS" or tid not in coding
            h = height if full else height * utr_ratio
            boxes.append(Rectangle((fs, y_of[tid] - h / 2), fe - fs, h))
        ax.add_collection(PatchCollection(boxes, facecolors=color, edgecolors="none", zorder=2))

        if labels is not None:
            F.draw_labels(ax, (t_start + t_end) / 2, rows * pitch, labels, height, fontsize)

    F.set_row_limits(ax, int(rows.max()) + 1 if len(span) else 1, pitch, height, show_names)
    finish_track(ax, region, label)
    hide_y(ax)
    return ax
