# Sequence logos (pygt.logo, pygt.logo_table)

```python
pygt.logo(df, ax, region, bases=("A", "C", "G", "T"), cols=("chrom", "start", "end"),
          colors=None, color="0.35", letters="auto", min_letter_px=2.0,
          flip_below=False, font=None, baseline=True, ylim=None,
          range_label=True, label=None)
pygt.logo_table(values, chrom, start, sequence=None, bases=("A", "C", "G", "T"))
```

## Input

- One row per position (normally 1 bp) with `chrom`, `start`, `end` and one
  numeric column per letter. Column names are the letters drawn; `bases=`
  picks them, so amino acids work (`bases=list("ACDEFGHIKLMNPQRSTVWY")`).
- NaN counts as zero. Overlapping rows raise `ValueError` (the same
  position twice would overdraw); deduplicate first.
- Model output is an array. Build the table with `pygt.logo_table`:
  - `(L, 4)` or channels-first `(4, L)` (PyTorch); a `(4, L)` array is
    transposed automatically. **Ambiguous when L == 4**: then it is read as
    `(L, 4)`.
  - `sequence=` with a 2-D array multiplies by the one-hot encoding:
    hypothetical contributions become actual ones. Do this explicitly; the
    track does not know which kind it was given.
  - a 1-D `(L,)` array (one score per position, e.g. DeepLIFT summed over
    bases) needs `sequence=` and is placed on the reference base.
  - `start` is the 0-based coordinate of the first position. A model input
    of length L centred on a summit at `s` starts at `s - L // 2`.
- From a batch `(N, 4, L)`, pass one example: `attr[i]`.

## Drawing

- Positives stack upward from 0, negatives downward, largest magnitude
  furthest from 0. `flip_below=True` draws negative letters upside down.
- `letters="auto"`: letters only when one position spans at least
  `min_letter_px` pixels of `ax`; otherwise bars of the summed positive
  and negative values per position, binned to screen resolution with
  `max`/`min`. On an 8 in figure that switch is at roughly 350 bp. Force
  with `letters=True` / `letters=False`. Like `pygt.signal`, it reads the
  axes' pixel width at call time, so create the figure at its final size.
- Default colours apply only to nucleotide alphabets (A green, C blue,
  G orange, T/U red). For other alphabets every letter not in `colors=` is
  `color`, so name the groups you care about.

## Scale

- Default range `[min(0, data), max(0, data)]`, per track. To compare
  models, alleles or conditions, pass the same `ylim` to each logo. Do not
  compare hypothetical and actual tracks by height without one.
- `range_label=True` prints `[min - max]` top-left and hides the y-axis.

## Geometry

`pygt.geometry.logo_stacks(df, region, bases)` returns the letter boxes
(`start`, `end`, `base`, `bottom`, `top`) and
`pygt.geometry.position_totals(df, region, bases)` the zoomed-out bars
(`pos`, `neg`), without drawing.
