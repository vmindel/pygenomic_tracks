# With pymisha

[pymisha](https://pypi.org/project/pymisha/) returns genomic data as DataFrames that start
with `chrom, start, end`, which is the pygt default. Its output goes into the track
functions unchanged. The examples below need a misha database, so they are shown but
not executed by the docs build. The column layout shown is from pymisha 0.14 and its
example database.

## What gextract returns

With pymisha's example database:

```py
import pymisha as pm
import pygenomic_tracks as pygt

pm.gdb_init_examples()
d = pm.gextract("dense_track", intervals=pm.gintervals("1", 0, 20000), iterator=50)
print(d.head(3))
```

```text
  chrom  start  end  dense_track  intervalID
0     1      0   50     0.177778           1
1     1     50  100     0.160000           1
2     1    100  150     0.180000           1
```

This table goes straight into `pygt.signal(d, ax, "1:0-20000", value="dense_track")`.

## Track expressions → signal

With your own database:

```py
pm.gdb_init("/path/to/trackdb")

region = pygt.Region("chr2", 1_200_000, 1_350_000)
scope = pm.gintervals(region.chrom, region.start, region.end)
d = pm.gextract("chec.msn2", intervals=scope, iterator=50)

fig, axes = pygt.figure(region, heights=[1.2, 1])
pygt.signal(d, axes[0], region, value="chec.msn2", label="Msn2")
pygt.genes("mm10", axes[1], region, transcripts="longest")
```

- The value column is named after the expression. Pass `colnames=["msn2"]` to
  `gextract` for a cleaner name.
- `iterator=` is the bin size. 20–50 bp is plenty for windows of 100 kb or more, since
  `pygt.signal` aggregates to about one bin per pixel anyway. A 1 bp iterator over a
  large window makes `gextract` slow.
- `intervalID` is ignored.
- Bins without data in sparse tracks are NaN and show as gaps.

## Several tracks at once

```py
d = pm.gextract(["chec.msn2", "chec.free"], intervals=scope, iterator=50,
                colnames=["msn2", "free"])

fig, axes = pygt.figure(region, heights=[1, 1])
pygt.signal(d, axes[0], region, value="msn2", ylim=60, label="Msn2")
pygt.signal(d, axes[1], region, value="free", ylim=60, color="0.4", label="free MNase")
```

Virtual tracks work the same way: define them with `pm.gvtrack_create(...)` and
extract them by name.

## Logical expressions → intervals

`gscreen` returns the intervals where an expression is true, with adjacent bins merged:

```py
peaks = pm.gscreen("chec.msn2 > 20", intervals=scope)      # chrom, start, end
pygt.intervals(peaks, ax, region, color="tab:red", label="Msn2 > 20")
```

## Interval sets

Saved interval sets load with their extra columns, which can drive colour and labels:

```py
cres = pm.gintervals_load("cres", chrom="chr2")
pygt.intervals(cres, ax, region, color_by="type", label_by="name", legend=True)
```

## Chromosome names and coordinates

- misha and pygt are both **0-based, half-open**. Never add or subtract 1 when moving
  tables between them.
- Chromosome names come from the database. A database built from UCSC FASTA uses
  `chr1`, which matches `pygt.annotation("mm10")`. Some databases (including the
  pymisha example database) use `1`. When a table has no rows on the requested name,
  pygt tries it with the `chr` prefix added or removed, so `"chr1:0-20000"` finds rows
  stored as `1`.
- 2D tracks (Hi-C contacts) are not supported by any pygt track yet.
