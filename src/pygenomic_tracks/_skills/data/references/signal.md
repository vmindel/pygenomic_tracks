# Signal tracks (pygt.signal)

```python
pygt.signal(df, ax, region, value="signal", cols=("chrom", "start", "end"),
            n_bins="auto", agg="max", color="C0", fill=True, alpha=1.0,
            linewidth=0.8, ylim=None, range_label=True, name=None, label=None)
```

## Drawing model

- Each row is a step over `[start, end)` (`ax.stairs`), so variable-width
  bins are drawn to scale. **Do not** replace this with
  `ax.plot(df.start, df.signal)`, which shifts every bin by half its width
  and draws lines across gaps.
- Gaps between rows and NaN values are left empty. pymisha `gextract`
  returns NaN where a sparse track has no data, and these show as holes.
  Fill them with `df.fillna({"signal": 0})` if zeros are what you mean.

## Binning (n_bins, agg)

- `n_bins="auto"`: when there are more rows than pixels across `ax`, the
  window is cut into about one bin per pixel and aggregated. Fewer rows are
  drawn as they are.
- `agg="max"` (default) keeps narrow peaks: a 50 bp ChEC peak survives in a
  1 Mb window. `"mean"` is weighted by overlap length and suits
  coverage. `"min"` exists for depletion tracks.
- **Overlapping rows** (e.g. gextract over overlapping scope intervals,
  where the same coordinates appear twice) are always binned under
  `"auto"`. With `n_bins=None` they raise `ValueError`.
- `n_bins=None` draws every row. Only use it for small windows; it is slow
  and the file large for 10^5+ rows.
- Auto-binning reads the axes' pixel width at call time, so create the
  figure at its final size first. `pygt.figure` does this.

## Y-scale

- Default range: `[min(0, data), max(0, data)]` over **every signal already
  drawn on this ax**, so overlays share a scale. Negative values (log2
  ratios) get a range that spans 0.
- To compare tracks on separate axes, pass the same `ylim`: either a single
  ymax (`ylim=60`) or `(ymin, ymax)`. Without it, each track scales to its
  own maximum and visual heights are not comparable. This is the commonest
  misleading figure.
- `range_label=True` (default) hides the y-axis and prints `[0 - 59.2]`
  top-left, as in IGV. `range_label=False` keeps a two-tick y-axis.

## Overlays

```python
pygt.signal(chec, ax, region, color="C0", alpha=0.7, name="Msn2")
pygt.signal(ctrl, ax, region, color="0.4", fill=False, name="free MNase")
```

`name=` adds a legend entry. Draw the larger signal first, or use
`fill=False` for the second, so one does not hide the other.

## Several value columns

pymisha `gextract(["a", "b"], ...)` returns one column per expression. Pick
one with `value="a"`; expressions with characters such as `.` or spaces are
the column names unless `colnames=` was given to gextract.
