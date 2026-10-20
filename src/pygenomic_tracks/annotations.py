"""Download, prepare and cache gene annotations for UCSC assemblies.

Any UCSC assembly name works (hg38, hg19, mm10, mm39, sacCer3, dm6, ce11,
danRer11, ...): GTFs are fetched from the assembly's ``bigZips/genes``
directory, trying the tracks in ``DEFAULT_TRACKS`` order until one exists.
UCSC chromosome names ("chr1", "chrI") match misha databases built from UCSC
sequences.

The raw GTF and the prepared exon/CDS table are kept in :func:`cache_dir`.
Compute nodes without internet can use a cache filled beforehand on a login
node with ``python -m pygenomic_tracks download hg38``.
"""
from __future__ import annotations

import os
import pickle
import shutil
import sys
import tempfile
import urllib.error
import urllib.request
import warnings
from pathlib import Path
from typing import Dict, Optional, Sequence, Tuple, Union

import pandas as pd

from ._coords import Region
from .io import read_gtf

UCSC_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/{genome}/bigZips/genes/{genome}.{track}.gtf.gz"
DEFAULT_TRACKS: Tuple[str, ...] = ("ncbiRefSeq", "refGene", "ensGene", "knownGene")
COLUMNS = ["chrom", "start", "end", "strand", "feature", "transcript_id", "gene_id", "gene_name"]

_memory: Dict[Tuple[str, float], pd.DataFrame] = {}


def cache_dir() -> Path:
    """``$PYGT_CACHE_DIR``, else ``$XDG_CACHE_HOME/pygenomic_tracks``, else ``~/.cache/pygenomic_tracks``."""
    if os.environ.get("PYGT_CACHE_DIR"):
        return Path(os.environ["PYGT_CACHE_DIR"]).expanduser()
    base = os.environ.get("XDG_CACHE_HOME") or Path.home() / ".cache"
    return Path(base).expanduser() / "pygenomic_tracks"


def _atomic_write(dest: Path, write) -> None:
    """Write via a temp file and rename, so concurrent jobs never see a partial file."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=dest.parent, prefix=dest.name + ".", suffix=".part")
    try:
        with os.fdopen(fd, "wb") as fh:
            write(fh)
        os.replace(tmp, dest)
    except BaseException:
        if os.path.exists(tmp):
            os.unlink(tmp)
        raise


def _is_missing(err: urllib.error.URLError) -> bool:
    if isinstance(err, urllib.error.HTTPError):
        return err.code == 404
    return isinstance(err.reason, FileNotFoundError)


def _fetch(url: str, dest: Path) -> None:
    with urllib.request.urlopen(url, timeout=60) as resp:
        _atomic_write(dest, lambda fh: shutil.copyfileobj(resp, fh, 1 << 20))


def _prepare(gtf_path: Path) -> pd.DataFrame:
    df = read_gtf(gtf_path, attributes=("gene_id", "gene_name", "transcript_id"))
    df["gene_name"] = df["gene_name"].fillna(df["gene_id"])
    df = df[COLUMNS].sort_values(["chrom", "start"], kind="stable").reset_index(drop=True)
    for col in ("chrom", "strand", "feature"):
        df[col] = df[col].astype("category")
    return df


def _gtf_path(genome: str, track: str) -> Path:
    return cache_dir() / genome / f"{genome}.{track}.gtf.gz"


def _table_path(genome: str, track: str) -> Path:
    return cache_dir() / genome / f"{genome}.{track}.pkl"


def _cached_track(genome: str, tracks: Sequence[str]) -> Optional[str]:
    for track in tracks:
        if _table_path(genome, track).exists() or _gtf_path(genome, track).exists():
            return track
    return None


def download(genome: str, track: Optional[str] = None, force: bool = False) -> Path:
    """Fetch and prepare the annotation for ``genome``; returns the prepared table's path.

    ``track`` picks a UCSC gene track (e.g. "knownGene" for GENCODE on hg38);
    by default the first of ``DEFAULT_TRACKS`` that UCSC provides is used.
    """
    tracks = (track,) if track else DEFAULT_TRACKS
    found = None if force else _cached_track(genome, tracks)
    if found is None:
        for t in tracks:
            url = UCSC_URL.format(genome=genome, track=t)
            try:
                print(f"pygenomic_tracks: trying {url}", file=sys.stderr)
                _fetch(url, _gtf_path(genome, t))
                found = t
                break
            except urllib.error.URLError as err:
                if not _is_missing(err):
                    raise
        if found is None:
            raise FileNotFoundError(
                f"no GTF for {genome!r} among tracks {list(tracks)}; see "
                f"https://hgdownload.soe.ucsc.edu/goldenPath/{genome}/bigZips/genes/ "
                "and pass track=, or give read_gtf a local file"
            )

    table = _table_path(genome, found)
    if force or not table.exists():
        df = _prepare(_gtf_path(genome, found))
        _atomic_write(table, lambda fh: pickle.dump(df, fh, protocol=pickle.HIGHEST_PROTOCOL))
    return table


def _load_table(path: Path, gtf: Path) -> pd.DataFrame:
    try:
        with open(path, "rb") as fh:
            return pickle.load(fh)
    except Exception:
        # A pickle from another pandas version may not load; rebuild from the raw GTF.
        df = _prepare(gtf)
        _atomic_write(path, lambda fh: pickle.dump(df, fh, protocol=pickle.HIGHEST_PROTOCOL))
        return df


def annotation(
    source: Union[str, os.PathLike, pd.DataFrame],
    track: Optional[str] = None,
    download_missing: bool = True,
) -> pd.DataFrame:
    """Exon/CDS table for :func:`genes`, from a genome name, a GTF path or a DataFrame.

    Results are kept in memory, so repeated calls within a session are free.
    Do not modify the returned DataFrame in place; it is shared.
    """
    if isinstance(source, pd.DataFrame):
        return source
    path = Path(source).expanduser()
    if path.suffix in (".gtf", ".gz") or path.is_file():
        key = (str(path.resolve()), path.stat().st_mtime)
        if key not in _memory:
            _memory[key] = _prepare(path)
        return _memory[key]

    genome = str(source)
    tracks = (track,) if track else DEFAULT_TRACKS
    found = _cached_track(genome, tracks)
    if found is None or not _table_path(genome, found).exists():
        if not download_missing:
            raise FileNotFoundError(
                f"{genome!r} is not in {cache_dir()}; run "
                f"'python -m pygenomic_tracks download {genome}' on a node with internet"
            )
        download(genome, track)
        found = _cached_track(genome, tracks)
    table = _table_path(genome, found)
    key = (str(table), table.stat().st_mtime)
    if key not in _memory:
        _memory[key] = _load_table(table, _gtf_path(genome, found))
    return _memory[key]


def gene_region(
    source: Union[str, os.PathLike, pd.DataFrame],
    gene: str,
    flank: int = 0,
    track: Optional[str] = None,
) -> Region:
    """Region spanning all transcripts of ``gene`` (gene_name or gene_id), plus ``flank`` bp."""
    df = annotation(source, track=track)
    hits = df[(df["gene_name"] == gene) | (df["gene_id"] == gene)]
    if hits.empty:
        hits = df[df["gene_name"].str.lower() == gene.lower()]
    if hits.empty:
        raise KeyError(f"gene {gene!r} not found")
    chroms = hits["chrom"].astype(str).value_counts()
    if len(chroms) > 1:
        warnings.warn(f"{gene} is on {list(chroms.index)}; using {chroms.index[0]}")
    hits = hits[hits["chrom"].astype(str) == chroms.index[0]]
    return Region(chroms.index[0], max(0, int(hits["start"].min()) - flank),
                  int(hits["end"].max()) + flank)


def cached() -> pd.DataFrame:
    """Genomes and tracks currently in the cache."""
    rows = []
    root = cache_dir()
    for p in sorted(root.glob("*/*.gtf.gz")):
        genome, track = p.name[: -len(".gtf.gz")].split(".", 1)
        table = _table_path(genome, track)
        rows.append({"genome": genome, "track": track, "prepared": table.exists(),
                     "size_mb": round(p.stat().st_size / 1e6, 1)})
    return pd.DataFrame(rows, columns=["genome", "track", "prepared", "size_mb"])
