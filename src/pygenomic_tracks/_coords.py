"""Region parsing and DataFrame selection.

All coordinates are 0-based, half-open (BED convention).
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Sequence, Tuple, Union

import numpy as np
import pandas as pd

DEFAULT_COLS: Tuple[str, str, str] = ("chrom", "start", "end")

_REGION_RE = re.compile(r"^(?P<chrom>[^:\s]+):(?P<start>\d+)-(?P<end>\d+)$")


@dataclass(frozen=True)
class Region:
    chrom: str
    start: int
    end: int

    def __post_init__(self):
        if self.start < 0:
            raise ValueError(f"region start must be >= 0, got {self.start}")
        if self.end <= self.start:
            raise ValueError(f"region end ({self.end}) must be > start ({self.start})")

    @property
    def length(self) -> int:
        return self.end - self.start

    def __str__(self) -> str:
        return f"{self.chrom}:{self.start:,}-{self.end:,}"


RegionLike = Union[Region, str, Sequence]


def as_region(region: RegionLike) -> Region:
    """Coerce ``"chr2:1,200,000-1,350,000"``, ``("chr2", 1200000, 1350000)`` or a Region."""
    if isinstance(region, Region):
        return region
    if isinstance(region, str):
        m = _REGION_RE.match(region.replace(",", "").replace("_", "").strip())
        if m is None:
            raise ValueError(f"cannot parse region {region!r}; expected 'chrom:start-end'")
        return Region(m["chrom"], int(m["start"]), int(m["end"]))
    if isinstance(region, (tuple, list)) and len(region) == 3:
        chrom, start, end = region
        return Region(str(chrom), int(start), int(end))
    raise TypeError(f"cannot interpret {region!r} as a region")


def check_columns(df: pd.DataFrame, columns: Sequence[str]) -> None:
    missing = [c for c in columns if c not in df.columns]
    if missing:
        raise KeyError(f"DataFrame is missing column(s) {missing}; has {list(df.columns)}")


def chrom_mask(chroms: pd.Series, chrom: str) -> pd.Series:
    """Rows on ``chrom``; falls back to the name with the "chr" prefix toggled.

    misha databases differ in naming ("chr1" vs "1"), so a region typed in
    one style still finds data stored in the other.
    """
    if pd.api.types.is_numeric_dtype(chroms):
        chroms = chroms.astype(str)
    mask = chroms == chrom
    if not mask.any():
        alt = chrom[3:] if chrom.startswith("chr") else "chr" + chrom
        mask = chroms == alt
    return mask


def select(df: pd.DataFrame, region: Region, cols: Sequence[str] = DEFAULT_COLS) -> pd.DataFrame:
    """Rows of ``df`` overlapping ``region``, sorted by start."""
    c, s, e = cols
    check_columns(df, cols)
    mask = chrom_mask(df[c], region.chrom) & (df[e] > region.start) & (df[s] < region.end)
    return df.loc[mask].sort_values(s, kind="stable")


def clipped(sub: pd.DataFrame, region: Region, cols: Sequence[str] = DEFAULT_COLS):
    """Start/end arrays clipped to the region, as float."""
    _, s, e = cols
    starts = np.clip(sub[s].to_numpy(dtype=float), region.start, region.end)
    ends = np.clip(sub[e].to_numpy(dtype=float), region.start, region.end)
    return starts, ends
