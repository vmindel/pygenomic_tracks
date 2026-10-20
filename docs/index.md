# pygenomic_tracks

**Genome-browser-style figures in Python**: signal, interval and gene tracks drawn on
matplotlib axes from plain pandas DataFrames. It covers what Gviz and plotgardener do
in R, built around tables instead of track objects.

![Example locus](images/demo.png)

```py
import pygenomic_tracks as pygt

region = "chr2:1,200,000-1,350,000"
fig, axes = pygt.figure(region, heights=[1.2, 0.6, 0.9, 1.0])
pygt.signal(chec, axes[0], region, value="signal", label="ChEC-seq")
pygt.intervals(cres, axes[1], region, color_by="cre", label_by="name", legend=True)
pygt.intervals(motif_hits, axes[2], region, color_by="score", palette="Greys")
pygt.genes("mm10", axes[3], region, transcripts="longest")
```

## Design

- **DataFrames in.** Every track is `pygt.<track>(df, ax, region, **column_args)`. A
  table needs `chrom`, `start`, `end` (0-based, half-open). Arguments name which other
  column is used and how: `value=`, `color_by=`, `label_by=`, `strand=`. No track
  classes, no config files.
- **pymisha-ready.** `gextract` and `gscreen` return this shape, so their output is
  plotted without conversion. BED files, bigWigs and moodspy motif positions work the
  same way.
- **Plain matplotlib.** Tracks draw on any `ax`. Highlights, annotations, insets and
  multi-panel layouts are ordinary matplotlib calls. `pygt.figure` sets up the common
  stacked layout.
- **Correct by default.** Signal is drawn as steps over the real interval extent, with
  gaps kept and aggregation to screen resolution (`max` by default, so narrow peaks
  survive). Short features such as motif hits are kept at least 1 px wide.
- **Annotations by name.** `pygt.genes("hg38", ...)` downloads and caches the UCSC
  annotation on first use. A one-line command fills the cache for offline cluster
  nodes.
- **Light dependencies.** numpy, pandas, matplotlib. pyBigWig is optional.
- **A [Claude Code skill](claude-skill.md)** ships with the package.

## Where next

| | |
|---|---|
| [Quick start](tutorials/quickstart.md) | the data model and a first four-track figure |
| [Signal tracks](tutorials/signal.md) | steps vs lines, binning, shared scales, overlays, ratios |
| [Interval tracks](tutorials/intervals.md) | peaks, CREs, motif hits; colour, labels, strand, packing |
| [Genes and annotations](tutorials/genes.md) | GTFs, genome names, the download cache, clusters |
| [Figure layout](tutorials/layout.md) | stacked tracks, highlights, panels, several loci, saving |
| [With pymisha](tutorials/pymisha.md) | from track expressions to figures |
| [API reference](api/tracks.md) | every public function |

Every tutorial code block with output is executed by the test suite. Its printed
output is checked and its figures are regenerated from the code shown.
