import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from matplotlib.patches import StepPatch

import pygenomic_tracks as pygt
from pygenomic_tracks import geometry as G


@pytest.fixture
def ax():
    fig, ax = plt.subplots(figsize=(6, 1.5))
    yield ax
    plt.close(fig)


def letters(ax):
    return [p for p in ax.patches if str(p.get_gid()).startswith("pygt-logo-")]


def table(rows, start=100):
    return pygt.logo_table(np.array(rows, dtype=float), "chr1", start)


def test_logo_table_shapes_and_one_hot():
    hyp = np.arange(12, dtype=float).reshape(3, 4)
    t = pygt.logo_table(hyp, "chr1", 10)
    assert t[["start", "end"]].to_numpy().tolist() == [[10, 11], [11, 12], [12, 13]]
    assert t[list("ACGT")].to_numpy().tolist() == hyp.tolist()

    channels_first = pygt.logo_table(hyp.T, "chr1", 10)
    pd.testing.assert_frame_equal(channels_first, t)

    actual = pygt.logo_table(hyp, "chr1", 10, sequence="gNt")
    assert actual[list("ACGT")].to_numpy().tolist() == [[0, 0, 2, 0], [0, 0, 0, 0], [0, 0, 0, 11]]

    per_pos = pygt.logo_table([0.5, -1.0], "chr1", 0, sequence="AC")
    assert per_pos[list("ACGT")].to_numpy().tolist() == [[0.5, 0, 0, 0], [0, -1.0, 0, 0]]


def test_logo_table_rejects_bad_shapes():
    with pytest.raises(ValueError, match="sequence"):
        pygt.logo_table([1.0, 2.0], "chr1", 0)
    with pytest.raises(ValueError, match="positions"):
        pygt.logo_table(np.ones((3, 4)), "chr1", 0, sequence="AC")
    with pytest.raises(ValueError, match="must be"):
        pygt.logo_table(np.ones((3, 5)), "chr1", 0)


def test_logo_stacks_order_from_zero_outward():
    t = table([[0.2, 0.5, -0.1, -0.3]])
    s = G.logo_stacks(t, ("chr1", 0, 200))
    assert s["base"].tolist() == ["A", "C", "G", "T"]
    np.testing.assert_allclose(s["bottom"], [0, 0.2, 0, -0.1])
    np.testing.assert_allclose(s["top"], [0.2, 0.7, -0.1, -0.4])


def test_logo_stacks_drop_zeros_nan_and_clip():
    t = table([[1, 0, 0, 0], [0, np.nan, 0, 0], [0, 0, 0, -2]])
    s = G.logo_stacks(t, ("chr1", 100, 102))
    assert s[["start", "base", "top"]].values.tolist() == [[100.0, "A", 1.0]]
    tot = G.position_totals(t, ("chr1", 0, 200))
    assert tot["pos"].tolist() == [1, 0, 0] and tot["neg"].tolist() == [0, 0, -2]


def test_logo_stacks_reject_overlaps_and_missing_columns():
    t = table([[1, 0, 0, 0], [1, 0, 0, 0]])
    with pytest.raises(ValueError, match="overlap"):
        G.logo_stacks(pd.concat([t, t]), ("chr1", 0, 200))
    with pytest.raises(KeyError):
        G.logo_stacks(t.drop(columns="T"), ("chr1", 0, 200))


def test_logo_draws_one_patch_per_letter_and_sets_range(ax):
    t = table([[0.2, 0.5, -0.1, -0.3], [0, 1.0, 0, 0]])
    pygt.logo(t, ax, ("chr1", 100, 102))
    gids = sorted(p.get_gid() for p in letters(ax))
    assert gids == [f"pygt-logo-{b}" for b in "ACGT"]
    assert ax.get_ylim() == pytest.approx((-0.4, 1.0))
    c = next(p for p in ax.patches if p.get_gid() == "pygt-logo-C")
    ext = c.get_path().get_extents()
    assert ext.y1 == pytest.approx(1.0) and ext.x0 > 100 and ext.x1 < 102


def test_logo_flip_below_mirrors_negative_letters(ax):
    t = table([[0, 0, 0, -1.0]])
    pygt.logo(t, ax, ("chr1", 100, 101))
    up = letters(ax)[0].get_path().vertices[:, 1]
    ax.clear()
    pygt.logo(t, ax, ("chr1", 100, 101), flip_below=True)
    down = letters(ax)[0].get_path().vertices[:, 1]
    np.testing.assert_allclose(down, -1.0 - up)


def test_logo_falls_back_to_bars_when_zoomed_out(ax):
    hyp = pygt.example_attributions()
    region = ("chr2", 1_321_000, 1_322_000)
    pygt.logo(hyp, ax, region)
    assert not letters(ax)
    assert sum(isinstance(p, StepPatch) for p in ax.patches) == 2
    lo, hi = ax.get_ylim()
    tot = G.position_totals(hyp, region)
    assert lo == pytest.approx(tot["neg"].min()) and hi == pytest.approx(tot["pos"].max())

    ax.clear()
    pygt.logo(hyp, ax, ("chr2", 1_321_450, 1_321_520))
    assert letters(ax) and not any(isinstance(p, StepPatch) for p in ax.patches)


def test_logo_options(ax):
    t = table([[1.0, 0, 0, 0]])
    pygt.logo(t, ax, ("chr1", 100, 101), colors={"A": "black"}, ylim=2, range_label=False)
    patch = letters(ax)[0]
    assert patch.get_facecolor()[:3] == (0, 0, 0)
    assert ax.get_ylim() == (0, 2)
    with pytest.raises(ValueError, match="letters"):
        pygt.logo(t, ax, ("chr1", 100, 101), letters="yes")


def test_logo_protein_letters(ax):
    df = pd.DataFrame({"chrom": "p", "start": [0, 1, 2], "end": [1, 2, 3],
                       "S": [0.5, 0, 0], "D": [0, -0.4, 0], "G": [0, 0, 0.3]})
    pygt.logo(df, ax, ("p", 0, 3), bases=["S", "D", "G"], color="0.2", colors={"D": "red"})
    face = {p.get_gid(): p.get_facecolor()[:3] for p in ax.patches}
    # Nucleotide colours must not leak onto amino acids that share a letter (G).
    assert face == {"pygt-logo-S": (0.2,) * 3, "pygt-logo-D": (1, 0, 0),
                    "pygt-logo-G": (0.2,) * 3}


def test_logo_empty_region(ax):
    pygt.logo(table([[1.0, 0, 0, 0]]), ax, ("chr1", 0, 50))
    assert not ax.patches and ax.get_ylim() == (0, 1)
