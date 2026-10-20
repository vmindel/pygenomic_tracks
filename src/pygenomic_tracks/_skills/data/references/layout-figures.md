# Figure layout

## pygt.figure

```python
fig, axes = pygt.figure(region, heights=[1.2, 0.5, 0.8, 1.0], width=8.0,
                        track_height=0.9, hspace=0.15)
```

- `heights`: relative heights, one per track, top to bottom (an int means
  that many equal tracks). Figure height = `sum(heights) * track_height + 0.6`
  inches.
- All axes share x. The bottom axis gets ticks in bp/kb/Mb (chosen from
  the window length) and the label `"chr2 (kb)"`. Setting xlim on any axes
  moves all of them.
- Call tracks **after** creating the figure at its final size. Auto-binning
  in `pygt.signal`, `min_width_px` and label packing measure the axes in
  pixels at call time. Resizing the figure afterwards (`fig.set_size_inches`)
  makes them stale; draw again instead.
- Track labels (`label=` on each track) are drawn left of the axes and
  extend outside the figure. Save with `bbox_inches="tight"`.

## Own axes

Any track works on any `ax`, e.g. inside a larger multi-panel figure:

```python
fig = plt.figure(figsize=(10, 5))
gs = fig.add_gridspec(3, 2, width_ratios=[3, 1])
ax_sig = fig.add_subplot(gs[0, 0])
ax_genes = fig.add_subplot(gs[1, 0], sharex=ax_sig)
pygt.signal(chec, ax_sig, region)
pygt.genes("mm10", ax_genes, region)
```

A standalone (non-shared) axes labels its own x-axis. Shared axes leave the
label to you: call `pygt.format_genomic_axis(bottom_ax, region)`.

## Annotating

Axes are in genomic coordinates, so plain matplotlib works:

```python
for ax in axes:
    ax.axvspan(1_261_000, 1_274_000, color="gold", alpha=0.2, lw=0)   # highlight
axes[0].axvline(1_320_500, color="k", lw=0.5, ls="--")
```

## Several loci

One figure per locus is simplest. For a side-by-side panel, build columns
without shared x across columns:

```python
fig, axes = plt.subplots(2, 2, figsize=(10, 3), sharex="col")
for col, reg in enumerate([r1, r2]):
    pygt.signal(chec, axes[0, col], reg, ylim=60)     # same ylim across loci
    pygt.genes("hg38", axes[1, col], reg, transcripts="longest")
```

## Saving

- PNG for slides: `fig.savefig("locus.png", dpi=300, bbox_inches="tight")`.
- PDF/SVG for papers. Set `matplotlib.rcParams["pdf.fonttype"] = 42` (and
  `"svg.fonttype": "none"`) **before** drawing, so text stays editable in
  Illustrator.
- Signal tracks are already reduced to about one step per pixel, so vector
  files stay small. `n_bins=None` on a big window produces huge PDFs.
