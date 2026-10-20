# Interval tracks (pygt.intervals)

```python
pygt.intervals(df, ax, region, cols=("chrom", "start", "end"), color="0.35",
               color_by=None, palette=None, vmin=None, vmax=None,
               label_by=None, strand=None, pack=True, height=0.6,
               min_width_px=1.0, edgecolor="none", fontsize=7,
               legend=False, label=None)
```

For any BED-like table: peaks, CREs, motif hits, pymisha `gscreen` output,
domains. **Default is one colour, no labels, no strand.** Everything else
is opt-in by naming a column.

## Colour

| want | call |
|---|---|
| one colour | `color="tab:red"` |
| by category (CRE class, motif name) | `color_by="cre"`, optionally `palette={"promoter": "#c0392b", ...}`, `legend=True` |
| by a number (motif score, -log10 p) | `color_by="score", palette="Greys", vmin=6, vmax=14` |

- Numeric vs categorical is decided by dtype. A numeric column you mean as
  categories (cluster ids 0..9) must be cast:
  `df.assign(cluster=df.cluster.astype(str))`.
- Categorical colours from a colormap are assigned in sorted category order
  **within the window**, so the same category can get a different colour in
  another region. For colours stable across figures, pass a dict
  `palette`. Categories missing from the dict get `color`.
- Numeric `vmin`/`vmax` default to the data range in the window. Fix them
  when comparing loci.
- There is no colourbar. Add one yourself if needed with
  `fig.colorbar(matplotlib.cm.ScalarMappable(norm, cmap), ax=ax)`.

## Short features (motif hits)

- `min_width_px=1` widens anything narrower than 1 px around its centre,
  so 12 bp motifs stay visible in a 150 kb window. Raise it to 2–3 for
  print. Widened boxes overstate width; for exact extents zoom in.
- Dense hits stack into several rows with `pack=True`. For hundreds of
  hits per window, `pack=False` (one row, overlaps drawn on top of each other)
  usually reads better. Combine with `color_by="motif"` to keep motifs
  apart, or draw one `pygt.intervals` call per motif on separate axes.
- `strand="strand"` draws arrowheads only on features at least 6 px wide,
  so motif hits in a large window get none.

## Labels and packing

- `label_by="name"` prints text under each feature. Packing reserves the
  label's estimated width, so labelled tracks need more rows. Give the
  track more height in `pygt.figure(heights=...)`.
- Packing is greedy first fit with a 4 px gap. Row order is not
  meaningful.
- Labels on hundreds of features are unreadable; filter the table first
  (e.g. top N by score) and label only those.
