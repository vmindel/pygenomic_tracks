"""Check a genomic table against a region before plotting it with pygenomic_tracks.

Reports what would make a track empty or wrong: missing columns,
chromosome-name mismatches, rows in the window, invalid coordinates,
1-based-looking starts and overlapping rows (which signal() bins).
Exits with 1 if the table cannot be plotted in the region.

    python check_table.py peaks.bed --region chr2:1,200,000-1,350,000
    python check_table.py extract.tsv --region chr2:1200000-1350000 --value chec.msn2
    python check_table.py genes.gtf --region chr2:1200000-1350000 --gtf

Reads .bed/.bedGraph (no header, BED column names), .tsv/.csv/.txt (header),
.parquet, .pkl and GTF (--gtf, or a .gtf/.gtf.gz path).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import pygenomic_tracks as pygt


def load(path: str, gtf: bool) -> pd.DataFrame:
    p = Path(path)
    name = p.name.lower()
    if gtf or name.endswith((".gtf", ".gtf.gz")):
        return pygt.read_gtf(p)
    if name.endswith((".bed", ".bed.gz", ".bedgraph", ".bedgraph.gz", ".narrowpeak", ".broadpeak")):
        return pygt.read_bed(p)
    if name.endswith(".parquet"):
        return pd.read_parquet(p)
    if name.endswith(".pkl"):
        return pd.read_pickle(p)
    sep = "," if name.endswith((".csv", ".csv.gz")) else "\t"
    return pd.read_csv(p, sep=sep)


def check(df: pd.DataFrame, region: pygt.Region, cols, value=None) -> list:
    """Return a list of (level, message); level is 'ERROR', 'WARN' or 'OK'."""
    out = []
    c, s, e = cols
    missing = [x for x in (*cols, *([value] if value else [])) if x not in df.columns]
    if missing:
        return [("ERROR", f"missing column(s) {missing}; columns are {list(df.columns)}")]

    chroms = df[c].astype(str)
    names = pd.unique(chroms)
    alt = region.chrom[3:] if region.chrom.startswith("chr") else "chr" + region.chrom
    if region.chrom in set(names):
        out.append(("OK", f"chromosome {region.chrom!r} present"))
    elif alt in set(names):
        out.append(("WARN", f"table uses {alt!r}, region uses {region.chrom!r}; "
                            "pygt matches them automatically"))
    else:
        out.append(("ERROR", f"no rows on {region.chrom!r} or {alt!r}; table chromosomes "
                             f"include {list(names[:10])}"))
        return out

    starts, ends = df[s].to_numpy(), df[e].to_numpy()
    bad = int(np.sum(ends <= starts))
    if bad:
        out.append(("WARN", f"{bad} rows have end <= start; they draw nothing"))
    if np.issubdtype(starts.dtype, np.number) and (starts < 0).any():
        out.append(("ERROR", f"{int((starts < 0).sum())} rows have negative start"))
    widths = ends - starts
    if len(widths) and (widths % 2 == 1).mean() > 0.9 and np.median(widths) > 1:
        out.append(("WARN", "over 90% of widths are odd: possibly 1-based closed coordinates "
                            "(GTF/R/IGV style); pygt expects 0-based half-open (start - 1)"))

    sub = pygt.select(df, region, cols)
    if sub.empty:
        lo = df.loc[chroms.isin([region.chrom, alt]), s].min()
        hi = df.loc[chroms.isin([region.chrom, alt]), e].max()
        out.append(("ERROR", f"0 rows overlap {region}; rows on this chromosome span {lo:,}-{hi:,}"))
        return out
    out.append(("OK", f"{len(sub):,} rows overlap {region}"))

    st, en = sub[s].to_numpy(), sub[e].to_numpy()
    n_overlap = int(np.sum(st[1:] < en[:-1]))
    if n_overlap:
        out.append(("WARN", f"{n_overlap} rows overlap their predecessor: signal() bins them "
                            "(agg=max by default); intervals() packs them into rows"))
    if value:
        v = pd.to_numeric(sub[value], errors="coerce")
        if v.isna().all():
            out.append(("ERROR", f"column {value!r} has no numeric values in the window"))
        else:
            nan = int(v.isna().sum())
            out.append(("OK", f"{value!r}: min {v.min():.3g}, max {v.max():.3g}, "
                              f"{nan} NaN rows (drawn as gaps)"))
            if len(sub) > 5000:
                out.append(("OK", f"{len(sub):,} rows: signal() aggregates to about one bin per pixel"))
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("table", help="BED/bedGraph/TSV/CSV/parquet/pickle/GTF file")
    ap.add_argument("--region", required=True, help="e.g. chr2:1,200,000-1,350,000")
    ap.add_argument("--value", help="numeric column for a signal track")
    ap.add_argument("--cols", default="chrom,start,end", help="chrom,start,end column names")
    ap.add_argument("--gtf", action="store_true", help="parse the table as GTF")
    args = ap.parse_args(argv)

    df = load(args.table, args.gtf)
    region = pygt.as_region(args.region)
    results = check(df, region, tuple(args.cols.split(",")), args.value)
    print(f"{args.table}: {len(df):,} rows, columns {list(df.columns)}")
    for level, msg in results:
        print(f"[{level}] {msg}")
    return 1 if any(level == "ERROR" for level, _ in results) else 0


if __name__ == "__main__":
    sys.exit(main())
