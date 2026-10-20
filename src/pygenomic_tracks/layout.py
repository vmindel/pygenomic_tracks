"""Figure layout and genomic axis formatting."""
from __future__ import annotations

from typing import List, Optional, Sequence, Tuple, Union

import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, MaxNLocator

from ._coords import Region, RegionLike, as_region


def _unit(length: int) -> Tuple[str, float]:
    if length >= 2_000_000:
        return "Mb", 1e6
    if length >= 2_000:
        return "kb", 1e3
    return "bp", 1.0


def format_genomic_axis(ax, region: RegionLike, xlabel: bool = True) -> None:
    """Set the x-limits to ``region`` and label ticks in bp, kb or Mb.

    Track functions call this themselves; use it directly for axes you draw
    on with plain matplotlib.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
    region : Region, str or tuple
    xlabel : bool
        Also set the axis label, e.g. ``"chr2 (kb)"``.
    """
    region = as_region(region)
    unit, scale = _unit(region.length)
    ax.set_xlim(region.start, region.end)
    ax.xaxis.set_major_locator(MaxNLocator(nbins=5))
    ax.xaxis.set_major_formatter(
        FuncFormatter(lambda x, _: f"{x / scale:,.3f}".rstrip("0").rstrip("."))
    )
    if xlabel:
        ax.set_xlabel(f"{region.chrom} ({unit})")


def _is_shared_x(ax) -> bool:
    return len(ax.get_shared_x_axes().get_siblings(ax)) > 1


def finish_track(ax, region: Region, label: Optional[str]) -> None:
    """Common end-of-track styling: x-axis, spines, track label."""
    # Shared axes get their x-label from figure(); a standalone ax labels itself.
    format_genomic_axis(ax, region, xlabel=not _is_shared_x(ax))
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if label is not None:
        ax.set_ylabel(label, rotation=0, ha="right", va="center")


def hide_y(ax) -> None:
    ax.spines["left"].set_visible(False)
    ax.set_yticks([])


def figure(
    region: RegionLike,
    heights: Union[int, Sequence[float]] = 3,
    width: float = 8.0,
    track_height: float = 0.9,
    hspace: float = 0.15,
) -> Tuple[plt.Figure, List]:
    """Create stacked track axes sharing one genomic x-axis.

    Parameters
    ----------
    region : Region, str or tuple
        The window; sets the x-limits and the axis label of the bottom track.
    heights : int or sequence of float
        Number of tracks, or their relative heights.
    width : float
        Figure width in inches.
    track_height : float
        Inches per unit of ``heights``.
    hspace : float
        Vertical space between tracks, as a fraction of the mean track height.

    Returns
    -------
    fig : matplotlib.figure.Figure
    axes : list of matplotlib.axes.Axes
        One per track, top to bottom.
    """
    region = as_region(region)
    if isinstance(heights, int):
        heights = [1.0] * heights
    heights = list(heights)
    fig, axes = plt.subplots(
        len(heights),
        1,
        figsize=(width, sum(heights) * track_height + 0.6),
        sharex=True,
        squeeze=False,
        gridspec_kw={"height_ratios": heights, "hspace": hspace},
    )
    axes = list(axes[:, 0])
    for ax in axes:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    for ax in axes[:-1]:
        ax.spines["bottom"].set_visible(False)
        ax.tick_params(axis="x", length=0)
    format_genomic_axis(axes[-1], region)
    return fig, axes
