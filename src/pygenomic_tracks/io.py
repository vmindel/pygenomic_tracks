"""Readers that return plain pandas DataFrames in 0-based, half-open coordinates."""
from __future__ import annotations

from typing import Optional, Sequence

import pandas as pd

from ._coords import RegionLike, as_region

BED_COLUMNS = [
    "chrom", "start", "end", "name", "score", "strand",
    "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts",
]

GTF_COLUMNS = ["chrom", "source", "feature", "start", "end", "score", "strand", "frame",
               "attributes"]


def read_bed(path, names: Optional[Sequence[str]] = None, **kwargs) -> pd.DataFrame:
    """Read a BED/bedGraph-like file. Columns beyond ``names`` get BED names, then ``colN``."""
    df = pd.read_csv(path, sep="\t", header=None, comment="#", dtype={0: str}, **kwargs)
    if names is None:
        names = BED_COLUMNS[: df.shape[1]] + [f"col{i}" for i in range(12, df.shape[1])]
    df.columns = list(names)
    return df


def read_gtf(
    path,
    features: Optional[Sequence[str]] = ("exon", "CDS"),
    attributes: Sequence[str] = ("gene_id", "gene_name", "transcript_id"),
) -> pd.DataFrame:
    """Read a GTF into one row per feature, ready for :func:`genes`.

    Start is shifted to 0-based. Each name in ``attributes`` becomes a column
    (NaN when absent). Cache the result (e.g. parquet) if you read a full
    GENCODE annotation repeatedly.
    """
    df = pd.read_csv(path, sep="\t", header=None, names=GTF_COLUMNS, comment="#",
                     dtype={"chrom": str})
    if features is not None:
        df = df.loc[df["feature"].isin(features)].copy()
    df["start"] = df["start"] - 1
    for a in attributes:
        df[a] = df["attributes"].str.extract(rf'(?:^|;)\s*{a} "([^"]*)"', expand=False)
    return df.drop(columns="attributes").reset_index(drop=True)


def from_bigwig(path, region: RegionLike, value: str = "signal") -> pd.DataFrame:
    """Intervals of a bigWig within ``region`` as a DataFrame (needs pyBigWig)."""
    try:
        import pyBigWig
    except ImportError as err:
        raise ImportError("from_bigwig needs pyBigWig: pip install 'pygenomic-tracks[bigwig]'") from err
    region = as_region(region)
    bw = pyBigWig.open(str(path))
    try:
        rows = bw.intervals(region.chrom, region.start, region.end) or []
    finally:
        bw.close()
    df = pd.DataFrame(rows, columns=["start", "end", value])
    df.insert(0, "chrom", region.chrom)
    return df
