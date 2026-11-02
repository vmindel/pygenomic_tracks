# Geometry

`pygenomic_tracks.geometry` computes what a track contains (signal steps,
gene parts, row assignments) without drawing it. The matplotlib tracks are
built on it, and other front-ends, such as plotnine-genomic-tracks, use it
so that binning, packing and annotation handling behave identically
everywhere.

Nothing here takes a matplotlib `Axes`. Where the screen matters, the
resolution is an argument: `n_bins` for signal, `bp_per_px` and `dpi` for
label-aware packing.

This module is public API; changes to it are breaking changes.

::: pygenomic_tracks.geometry
    options:
      members:
        - signal_steps
        - bin_signal
        - to_stairs
        - gene_models
        - feature_rows
        - text_width_bp
        - widen
        - pack_rows
        - genomic_unit
