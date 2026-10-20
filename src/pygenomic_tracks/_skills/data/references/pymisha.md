# pymisha workflows

pymisha returns DataFrames with `chrom, start, end` first, which is the
pygt default. No conversion is needed.

## gextract → signal

```python
import pymisha as pm
pm.gdb_init("/path/to/trackdb")

region = pygt.Region("chr2", 1_200_000, 1_350_000)
scope = pm.gintervals(region.chrom, region.start, region.end)
d = pm.gextract("chec.msn2", intervals=scope, iterator=50)
# columns: chrom, start, end, chec.msn2, intervalID
pygt.signal(d, ax, region, value="chec.msn2")
```

- The value column is named after the expression. Use
  `colnames=["msn2"]` in gextract for a clean name, or pass the expression
  string as `value=`.
- `iterator=` sets the bin size. For a 150 kb window, 20–50 bp is plenty;
  `pygt.signal` aggregates further to pixels anyway. A 1 bp iterator over
  Mb windows is slow in gextract, not in pygt.
- Sparse tracks give NaN bins; they show as gaps (see `references/signal.md`).
- Several expressions at once: `pm.gextract(["a", "b"], ...)`, then one
  `pygt.signal` call per column.
- Virtual tracks work the same: define with `pm.gvtrack_create(...)`, then
  extract by name.
- 2D tracks (Hi-C) are not supported by any pygt track.

## gscreen → intervals

```python
peaks = pm.gscreen("chec.msn2 > 20", intervals=scope)   # chrom, start, end
pygt.intervals(peaks, ax, region, color="tab:red")
```

`gscreen` returns merged runs of iterator bins where the expression holds.
Add value columns yourself (e.g. `pm.gextract` over `peaks` with an
`intervals_join="intervals"` join, or a pandas merge) for `color_by=`.

## Interval sets and annotations

- Saved intervals sets load as DataFrames: `pm.gintervals_load("cres")` →
  `pygt.intervals(...)`. Extra columns stored with the set are available
  for `color_by=` and `label_by=`.
- Chromosome names come from the misha database. A DB built from UCSC FASTA
  uses `chr1`, matching `pygt.annotation("hg38")`. Some DBs (and the
  pymisha example DB) use `1`; pygt tries the `chr`-toggled name
  automatically when a table has no rows on the requested name.

## Coordinates

misha and pygt are both 0-based, half-open. Do not add or subtract 1 when
moving tables between them.
