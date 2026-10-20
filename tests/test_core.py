import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

import pygenomic_tracks as pygt
from pygenomic_tracks._coords import Region
from pygenomic_tracks.signal import _bin_signal, _to_stairs


@pytest.fixture
def ax():
    fig, ax = plt.subplots(figsize=(6, 1.5))
    yield ax
    plt.close(fig)


# --- regions and selection -------------------------------------------------

def test_region_parsing():
    assert pygt.as_region("chr2:1,200,000-1,350,000") == Region("chr2", 1_200_000, 1_350_000)
    assert pygt.as_region(("chrX", 10, 20)) == Region("chrX", 10, 20)
    with pytest.raises(ValueError):
        pygt.as_region("chr2:100-50")
    with pytest.raises(ValueError):
        pygt.as_region("chr2-100-200")


def test_select_overlap_and_chr_prefix():
    df = pd.DataFrame({"chrom": ["1", "1", "1", "2"], "start": [0, 90, 200, 100],
                       "end": [10, 110, 300, 150]})
    sub = pygt.select(df, Region("chr1", 100, 250))
    assert sub["start"].tolist() == [90, 200]


def test_select_categorical_chrom():
    df = pd.DataFrame({"chrom": pd.Categorical(["chr1", "chr2"]), "start": [0, 0], "end": [5, 5]})
    assert len(pygt.select(df, Region("chr2", 0, 10))) == 1


def test_missing_column_message():
    with pytest.raises(KeyError, match="missing"):
        pygt.select(pd.DataFrame({"chr": ["1"], "start": [0], "end": [1]}), Region("1", 0, 1))


# --- packing ----------------------------------------------------------------

def test_pack_rows():
    assert pygt.pack_rows([0, 20, 40], [10, 30, 50]).tolist() == [0, 0, 0]
    assert pygt.pack_rows([0, 5, 12], [10, 15, 20]).tolist() == [0, 1, 0]
    assert pygt.pack_rows([0, 11], [10, 20], gap=5).tolist() == [0, 1]


# --- signal ----------------------------------------------------------------

def test_to_stairs_inserts_nan_gaps():
    edges, vals = _to_stairs(np.array([0., 10, 30]), np.array([10., 20, 40]), np.array([1., 2, 3]))
    assert edges.tolist() == [0, 10, 20, 30, 40]
    assert vals[[0, 1, 3]].tolist() == [1, 2, 3]
    assert np.isnan(vals[2])


def test_bin_signal_max_and_weighted_mean():
    region = Region("c", 0, 100)
    starts = np.array([0., 25, 50])
    ends = np.array([25., 50, 100])
    vals = np.array([1., 3, 5])
    _, mx = _bin_signal(starts, ends, vals, region, 2, "max")
    assert mx.tolist() == [3, 5]
    _, mean = _bin_signal(starts, ends, vals, region, 2, "mean")
    assert mean.tolist() == [2, 5]


def test_bin_signal_empty_bins_are_nan():
    _, out = _bin_signal(np.array([0.]), np.array([10.]), np.array([1.]), Region("c", 0, 100), 10, "max")
    assert out[0] == 1 and np.isnan(out[1:]).all()


def test_signal_overlap_requires_binning(ax):
    df = pd.DataFrame({"chrom": "c", "start": [0, 5], "end": [10, 15], "signal": [1, 2]})
    with pytest.raises(ValueError, match="overlap"):
        pygt.signal(df, ax, ("c", 0, 20), n_bins=None)
    pygt.signal(df, ax, ("c", 0, 20))  # auto-binning handles it


def test_signal_draws_and_sets_limits(ax):
    df = pd.DataFrame({"chrom": "chr1", "start": np.arange(0, 1000, 10),
                       "end": np.arange(10, 1010, 10), "v": np.linspace(0, 9, 100)})
    pygt.signal(df, ax, "chr1:0-1000", value="v", n_bins=None, range_label=False)
    assert ax.get_xlim() == (0, 1000)
    assert ax.get_ylim()[1] == pytest.approx(9)


def test_signal_empty_region(ax):
    df = pd.DataFrame({"chrom": "chr1", "start": [0], "end": [10], "signal": [1.0]})
    pygt.signal(df, ax, "chr2:0-100")


# --- intervals ---------------------------------------------------------------

def test_intervals_single_colour_default(ax):
    df = pd.DataFrame({"chrom": "c", "start": [0, 50], "end": [20, 70], "cre": ["a", "b"]})
    pygt.intervals(df, ax, ("c", 0, 100))
    faces = ax.collections[0].get_facecolors()
    assert len(faces) == 2 and np.allclose(faces[0], faces[1])
    assert ax.get_legend() is None


def test_intervals_categorical_and_numeric_colour(ax):
    df = pd.DataFrame({"chrom": "c", "start": [0, 50], "end": [20, 70],
                       "cre": ["a", "b"], "score": [1.0, 9.0]})
    pygt.intervals(df, ax, ("c", 0, 100), color_by="cre", legend=True)
    faces = ax.collections[0].get_facecolors()
    assert not np.allclose(faces[0], faces[1])
    assert [t.get_text() for t in ax.get_legend().get_texts()] == ["a", "b"]

    ax.clear()
    pygt.intervals(df, ax, ("c", 0, 100), color_by="score", palette="Greys")
    faces = ax.collections[0].get_facecolors()
    assert faces[1][:3].sum() < faces[0][:3].sum()  # higher score is darker in Greys


def test_intervals_min_width_keeps_short_features_visible(ax):
    df = pd.DataFrame({"chrom": "c", "start": [500_000], "end": [500_010]})
    pygt.intervals(df, ax, ("c", 0, 1_000_000), min_width_px=2)
    width = ax.collections[0].get_paths()[0].get_extents().width
    assert width > 10


def test_intervals_labels_and_packing(ax):
    df = pd.DataFrame({"chrom": "c", "start": [0, 5], "end": [10, 15], "name": ["x", "y"]})
    pygt.intervals(df, ax, ("c", 0, 100), label_by="name")
    assert {t.get_text() for t in ax.texts} == {"x", "y"}
    ys = [p.get_extents().y0 for p in ax.collections[0].get_paths()]
    assert ys[0] != ys[1]


# --- genes ---------------------------------------------------------------------

def _gtf():
    return pd.DataFrame({
        "chrom": "chr1",
        "feature": ["exon", "exon", "CDS", "CDS", "exon", "exon"],
        "start": [100, 400, 150, 400, 100, 400],
        "end": [200, 500, 200, 450, 200, 700],
        "strand": "+",
        "transcript_id": ["t1", "t1", "t1", "t1", "t2", "t2"],
        "gene_name": "G",
    })


def test_genes_draws_transcripts(ax):
    pygt.genes(_gtf(), ax, "chr1:0-1000")
    assert [t.get_text() for t in ax.texts] == ["G", "G"]


def test_genes_longest(ax):
    pygt.genes(_gtf(), ax, "chr1:0-1000", transcripts="longest")
    assert len(ax.texts) == 1


def test_genes_outside_region(ax):
    pygt.genes(_gtf(), ax, "chr1:5000-6000")
    assert len(ax.texts) == 0


# --- io and layout ---------------------------------------------------------------

def test_read_gtf(tmp_path):
    p = tmp_path / "a.gtf"
    p.write_text(
        "#!genome-build test\n"
        'chr1\tsrc\tgene\t1\t100\t.\t+\t.\tgene_id "g1"; gene_name "G";\n'
        'chr1\tsrc\texon\t1\t100\t.\t+\t.\tgene_id "g1"; transcript_id "t1"; gene_name "G";\n'
        'chr1\tsrc\tCDS\t11\t50\t.\t+\t0\tgene_id "g1"; transcript_id "t1"; havana_gene_name "H";\n'
    )
    df = pygt.read_gtf(p)
    assert df["feature"].tolist() == ["exon", "CDS"]
    assert df["start"].tolist() == [0, 10] and df["end"].tolist() == [100, 50]
    assert df["transcript_id"].tolist() == ["t1", "t1"]
    assert df["gene_name"].iloc[0] == "G" and pd.isna(df["gene_name"].iloc[1])


def test_read_bed(tmp_path):
    p = tmp_path / "a.bed"
    p.write_text("chr1\t0\t10\tp1\t5\t+\nchr1\t20\t30\tp2\t7\t-\n")
    df = pygt.read_bed(p)
    assert list(df.columns) == ["chrom", "start", "end", "name", "score", "strand"]


def test_figure_shares_x():
    fig, axes = pygt.figure("chr1:1,000,000-3,500,000", heights=[1, 0.5])
    assert len(axes) == 2
    assert axes[0].get_xlim() == (1_000_000, 3_500_000)
    assert axes[-1].get_xlabel() == "chr1 (Mb)"
    plt.close(fig)


def test_signal_overlay_shares_range(ax):
    region = pygt.example_region()
    pygt.signal(pygt.example_signal("control"), ax, region, name="ctrl")
    low = ax.get_ylim()[1]
    pygt.signal(pygt.example_signal("chec"), ax, region, name="chec")
    assert ax.get_ylim()[1] > low
    pygt.signal(pygt.example_signal("control"), ax, region)
    assert ax.get_ylim()[1] > low  # a smaller third series does not shrink the range
    assert len([t for t in ax.texts if t.get_gid() == "pygt-range"]) == 1
    assert [t.get_text() for t in ax.get_legend().get_texts()] == ["ctrl", "chec"]


def test_signal_after_clear_resets_range(ax):
    region = pygt.example_region()
    pygt.signal(pygt.example_signal("chec"), ax, region)
    ax.clear()
    pygt.signal(pygt.example_signal("control"), ax, region)
    assert ax.get_ylim()[1] < 20
