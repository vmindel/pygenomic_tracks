---
name: pygenomic-tracks
description: Use for any task involving the pygenomic_tracks library (import pygenomic_tracks as pygt) — genome-browser-style figures in Python with matplotlib from pandas DataFrames. Triggers on plotting a genomic region or locus; signal / coverage / ChEC-seq / ChIP-seq / ATAC tracks from pymisha gextract, bigWig or bedGraph; BED-like interval tracks such as peaks, CREs or motif hits; gene models from a GTF or a genome name (hg38, mm10, sacCer3); stacking tracks into one figure; Gviz or plotgardener-style plots in Python. This is a router skill — read the relevant file under references/ before writing pygt code.
---

# pygenomic_tracks

Draws genome-browser tracks on matplotlib axes from plain pandas
DataFrames. Every track function has the shape
`pygt.<track>(df, ax, region, **column_args)`: the DataFrame holds
`chrom`, `start`, `end` plus value columns, and arguments name which column
is used and how. pymisha's `gextract`/`gscreen` output goes in unchanged.
This file only routes; the reference files hold the signatures, defaults
and silent failures. Read the file for the task before writing code.

## Task to reference file

| task | read |
|---|---|
| **input tables**: required columns, 0-based coordinates, `cols=`, chromosome naming, BED/bedGraph/bigWig readers | `references/inputs-coordinates.md` |
| **signal tracks**: `pygt.signal`, binning/aggregation, NaN gaps, shared y-scales, overlays, negative values | `references/signal.md` |
| **interval tracks**: `pygt.intervals` for peaks, CREs, motif hits; colours, labels, strand, packing | `references/intervals.md` |
| **gene models and annotations**: `pygt.genes`, `pygt.annotation`, `pygt.gene_region`, cache, offline cluster nodes | `references/genes-annotations.md` |
| **figure layout**: `pygt.figure`, heights, own axes, highlights, several loci, saving PDF/PNG | `references/layout-figures.md` |
| **pymisha workflows**: `gextract` → signal, `gscreen` → intervals, virtual tracks, intervals sets | `references/pymisha.md` |
| "**why is my track empty / wrong?**": check a table against a region | `scripts/check_table.py` |

`scripts/check_table.py` (path relative to this skill's root) reads a
TSV/CSV/BED/parquet table and reports, for a region, missing columns,
chromosome-name mismatches, rows in the window, invalid or 1-based-looking
coordinates and overlaps. Run `python scripts/check_table.py --help`.

## Conventions

- `import pygenomic_tracks as pygt`. Track functions:
  `pygt.signal(df, ax, region, value="signal", ...)`,
  `pygt.intervals(df, ax, region, ...)`,
  `pygt.genes(df_or_genome, ax, region, ...)`. Each returns `ax`.
- `region` is `"chr2:1,200,000-1,350,000"`, `("chr2", 1200000, 1350000)` or
  `pygt.Region`. It must be passed to every track; it is never inferred.
- Coordinates are **0-based, half-open** (BED, misha, Python). GTF is
  1-based; `pygt.read_gtf` and `pygt.annotation` convert it.
- Default column names are `chrom`, `start`, `end`; rename with
  `cols=("chr", "from", "to")`. No bioframe or pyranges dependency.
- `pygt.figure(region, heights=[...])` returns `(fig, axes)` with a shared
  x-axis. Tracks also work on any single `ax`.
- `pygt.intervals` draws one colour unless `color_by=` is given. A numeric
  `color_by` uses a colormap, anything else is categorical.
- If a table has no rows on the region's chromosome, the name with the
  `chr` prefix toggled is tried (`chr2` and `2`). Anything else that does
  not match gives an empty track **without an error**.
