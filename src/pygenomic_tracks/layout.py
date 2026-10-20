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
    """Set x-limits to the region and label ticks in bp/kb/Mb."""
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
    """Stacked track axes sharing a genomic x-axis.

    ``heights`` is either the number of tracks or their relative heights.
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
