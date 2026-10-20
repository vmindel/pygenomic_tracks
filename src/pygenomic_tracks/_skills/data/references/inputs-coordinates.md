# Input tables and coordinates

## Shape

Every track takes a pandas DataFrame with three coordinate columns plus
whatever value columns the track uses:

| track | needs | optional columns, named by argument |
|---|---|---|
| `pygt.signal` | `chrom, start, end` + numeric `value=` column | — |
| `pygt.intervals` | `chrom, start, end` | `color_by=`, `label_by=`, `strand=` |
| `pygt.genes` | `chrom, start, end, feature, transcript_id` | `gene_name`, `strand` |

- Other columns (pymisha's `intervalID`, BED `score`, ...) are ignored.
- Column names differ? Pass `cols=("chr", "from", "to")` rather than
  renaming. The order is chrom, start, end.
- A missing column raises `KeyError` that lists the columns the DataFrame
  has.

## Coordinates

- **0-based, half-open** everywhere: BED, bedGraph, misha/pymisha, Python
  slicing. A 12 bp motif at the first base is `start=0, end=12`.
- GTF/GFF are 1-based and closed. Only `pygt.read_gtf` and
  `pygt.annotation` convert (start - 1). If you parse a GTF yourself,
  subtract 1 from `start`.
- Region strings are read as 0-based: `"chr2:1,200,000-1,350,000"` is
  `Region("chr2", 1200000, 1350000)`. Commas and underscores are ignored.
  A coordinate copied from the UCSC/IGV location bar is 1-based, so it is
  1 bp off. That doesn't matter at track scale but does when matching single
  bases.
- `region` is the drawing window. Rows overlapping it are drawn, clipped
  at its edges; the rest are ignored, so a whole-genome table can be passed
  directly. Filtering 10^7 rows costs a few hundred ms per call. For
  many figures from one large table, subset it once with `pygt.select`.

## Chromosome names

- Matching is exact string equality on `chrom`. Numeric chrom columns are
  compared as strings.
- **Silent failure:** a region on a chromosome the table does not have
  draws an empty track with no error. The one exception: when no row matches,
  the `chr`-toggled name is tried (`chr2` ↔ `2`, `chrX` ↔ `X`). Other
  mismatches still fail silently: `chrM` vs `MT`, Ensembl vs UCSC scaffold
  names, yeast `chrI` vs `chr1`.
- Check with `df.chrom.unique()` or `scripts/check_table.py`.

## Readers

```python
peaks = pygt.read_bed("peaks.narrowPeak")       # chrom start end name score strand ... (BED names)
bg = pygt.read_bed("cov.bedGraph", names=["chrom", "start", "end", "signal"])
bw = pygt.from_bigwig("cov.bw", region)          # needs pyBigWig; chrom start end signal
gtf = pygt.read_gtf("genes.gtf.gz")              # exon + CDS rows, 0-based
```

- `pygt.read_bed` names columns by BED position. A bedGraph read without
  `names=` gets its value column called `name`, so pass `names=`.
- `pygt.read_bed` does not skip `track`/`browser` header lines; strip them
  first or pass `skiprows=`.
- `pygt.from_bigwig` returns only the stored intervals in the window, at
  their native resolution, so a 1 bp bigWig over 1 Mb is 10^6 rows.
  `pygt.signal` bins it to about one bin per pixel, which is fine;
  reading it is the slow part.
- Any DataFrame works: pymisha output, `pd.read_parquet`, moodspy
  positions, pyranges after `.df` and renaming
  (`Chromosome/Start/End` → pass `cols=("Chromosome", "Start", "End")`).
