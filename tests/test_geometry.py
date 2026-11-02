"""The public geometry API, without any matplotlib Axes."""
import numpy as np
import pandas as pd
import pytest

import pygenomic_tracks as pygt
from pygenomic_tracks import geometry as G


def intervals(starts, ends, values, chrom="chr1"):
    return pd.DataFrame({"chrom": chrom, "start": starts, "end": ends, "signal": values})


def test_signal_steps_raw_rows_are_clipped_and_keep_gaps():
    df = intervals([0, 10, 30, 90], [10, 20, 40, 120], [1.0, np.nan, 3.0, 4.0])
    steps = G.signal_steps(df, ("chr1", 5, 100))
    assert steps.to_dict("list") == {
        "start": [5.0, 30.0, 90.0],
        "end": [10.0, 40.0, 100.0],
        "value": [1.0, 3.0, 4.0],
    }


def test_signal_steps_bins_when_rows_outnumber_bins():
    df = intervals(np.arange(0, 100, 10), np.arange(10, 110, 10), np.arange(10.0))
    steps = G.signal_steps(df, ("chr1", 0, 100), n_bins=2, agg="max")
    assert steps["value"].tolist() == [4.0, 9.0]
    assert steps[["start", "end"]].to_numpy().tolist() == [[0, 50], [50, 100]]


def test_signal_steps_bins_overlapping_rows_and_rejects_them_unbinned():
    df = intervals([0, 5], [10, 15], [1.0, 2.0])
    assert G.signal_steps(df, ("chr1", 0, 20), n_bins=4)["value"].max() == 2.0
    with pytest.raises(ValueError, match="overlap"):
        G.signal_steps(df, ("chr1", 0, 20))


def test_signal_steps_matches_the_matplotlib_track():
    df, region = pygt.example_signal(), pygt.example_region()
    steps = G.signal_steps(df, region, n_bins=300, agg="mean")
    edges, heights = G.bin_signal(
        *G.clipped(G.select(df, region), region),
        df.pipe(G.select, region)["signal"].to_numpy(float),
        region, 300, "mean",
    )
    np.testing.assert_allclose(steps["value"], heights[~np.isnan(heights)])


def test_signal_steps_checks_arguments():
    df = intervals([0], [10], [1.0])
    with pytest.raises(ValueError, match="agg"):
        G.signal_steps(df, ("chr1", 0, 10), agg="median")
    with pytest.raises(KeyError):
        G.signal_steps(df, ("chr1", 0, 10), value="missing")


def test_signal_steps_toggles_chr_prefix():
    df = intervals([0], [10], [1.0], chrom="1")
    assert len(G.signal_steps(df, ("chr1", 0, 10))) == 1


def test_feature_rows_leave_room_for_labels():
    starts, ends = np.array([0.0, 120.0]), np.array([100.0, 220.0])
    # 1 bp per px: the features are 20 px apart
    assert G.feature_rows(starts, ends, bp_per_px=1).tolist() == [0, 0]
    # a 30-character label at 7 pt is ~175 px wide and needs its own row
    labels = ["x" * 30, "y"]
    assert G.feature_rows(starts, ends, bp_per_px=1, labels=labels).tolist() == [0, 1]
    assert G.feature_rows(starts, ends, bp_per_px=1, labels=labels, pack=False).tolist() == [0, 0]


def test_widen_keeps_centre():
    s, e = G.widen(np.array([10.0, 0.0]), np.array([12.0, 100.0]), 10)
    assert s.tolist() == [6.0, 0.0] and e.tolist() == [16.0, 100.0]


def test_gene_models():
    genes, region = pygt.example_genes(), pygt.example_region()
    spans, parts = G.gene_models(genes, region)
    assert spans["start"].is_monotonic_increasing
    assert {"start", "end", "name", "strand"} <= set(spans.columns)
    assert set(parts["transcript"]) == set(spans.index)
    assert set(parts["feature"]) <= {"exon", "CDS"}
    coding = set(parts.loc[parts["feature"] == "CDS", "transcript"])
    utr = parts[~parts["full"]]
    assert (utr["feature"] == "exon").all() and utr["transcript"].isin(coding).all()


def test_gene_models_longest_keeps_one_per_gene():
    genes, region = pygt.example_genes(), pygt.example_region()
    spans, _ = G.gene_models(genes, region, transcripts="longest")
    assert spans["name"].is_unique
    with pytest.raises(ValueError):
        G.gene_models(genes, region, transcripts="shortest")


@pytest.mark.parametrize(
    "length, unit", [(500, "bp"), (2_000, "kb"), (1_999_999, "kb"), (2_000_000, "Mb")]
)
def test_genomic_unit(length, unit):
    assert G.genomic_unit(length)[0] == unit


def test_old_private_names_still_import():
    from pygenomic_tracks.signal import _bin_signal, _to_stairs

    assert _bin_signal is G.bin_signal and _to_stairs is G.to_stairs
