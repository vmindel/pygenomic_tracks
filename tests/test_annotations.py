import gzip

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pytest

import pygenomic_tracks as pygt
from pygenomic_tracks import annotations
from pygenomic_tracks.__main__ import main

GTF = """\
chr1\ttest\ttranscript\t1001\t5000\t.\t+\t.\tgene_id "Myc"; transcript_id "NM_1"; gene_name "Myc";
chr1\ttest\texon\t1001\t1500\t.\t+\t.\tgene_id "Myc"; transcript_id "NM_1"; exon_number "1"; gene_name "Myc";
chr1\ttest\texon\t4001\t5000\t.\t+\t.\tgene_id "Myc"; transcript_id "NM_1"; exon_number "2"; gene_name "Myc";
chr1\ttest\tCDS\t1201\t1500\t.\t+\t0\tgene_id "Myc"; transcript_id "NM_1"; gene_name "Myc";
chr1\ttest\tCDS\t4001\t4600\t.\t+\t0\tgene_id "Myc"; transcript_id "NM_1"; gene_name "Myc";
chr2\ttest\texon\t101\t300\t.\t-\t.\tgene_id "NR_9"; transcript_id "NR_9";
"""


@pytest.fixture
def ucsc(tmp_path, monkeypatch):
    """Fake UCSC served from file:// with only the refGene track for genome 'tg1'."""
    server = tmp_path / "server"
    (server / "tg1").mkdir(parents=True)
    with gzip.open(server / "tg1" / "tg1.refGene.gtf.gz", "wt") as fh:
        fh.write(GTF)
    monkeypatch.setattr(annotations, "UCSC_URL",
                        server.as_uri() + "/{genome}/{genome}.{track}.gtf.gz")
    monkeypatch.setenv("PYGT_CACHE_DIR", str(tmp_path / "cache"))
    annotations._memory.clear()
    return tmp_path


def test_download_falls_back_to_available_track(ucsc):
    table = pygt.download("tg1")
    assert table.name == "tg1.refGene.pkl"
    assert (ucsc / "cache" / "tg1" / "tg1.refGene.gtf.gz").exists()


def test_unknown_genome_lists_tracks_tried(ucsc):
    with pytest.raises(FileNotFoundError, match="ncbiRefSeq"):
        pygt.download("nope1")


def test_annotation_table(ucsc):
    df = pygt.annotation("tg1")
    assert set(df["feature"]) == {"exon", "CDS"}
    assert df["start"].min() == 100  # 0-based
    # transcripts without gene_name fall back to gene_id
    assert df.loc[df["chrom"] == "chr2", "gene_name"].iloc[0] == "NR_9"
    assert pygt.annotation("tg1") is df  # in-memory cache


def test_offline_without_cache(ucsc):
    with pytest.raises(FileNotFoundError, match="python -m pygenomic_tracks download"):
        pygt.annotation("tg1", download_missing=False)


def test_gene_region(ucsc):
    assert pygt.gene_region("tg1", "Myc", flank=500) == pygt.Region("chr1", 500, 5500)
    assert pygt.gene_region("tg1", "MYC").start == 1000  # case-insensitive fallback
    with pytest.raises(KeyError):
        pygt.gene_region("tg1", "Sox2")


def test_genes_accepts_genome_name(ucsc):
    fig, ax = plt.subplots()
    pygt.genes("tg1", ax, pygt.gene_region("tg1", "Myc", flank=1000))
    assert [t.get_text() for t in ax.texts] == ["Myc"]
    plt.close(fig)


def test_corrupt_table_is_rebuilt(ucsc):
    table = pygt.download("tg1")
    table.write_bytes(b"not a pickle")
    annotations._memory.clear()
    assert len(pygt.annotation("tg1")) == 5


def test_annotation_from_gtf_path(tmp_path):
    p = tmp_path / "x.gtf"
    p.write_text(GTF)
    assert len(pygt.annotation(p)) == 5


def test_cli(ucsc, capsys):
    assert main(["download", "tg1"]) == 0
    assert main(["list"]) == 0
    out = capsys.readouterr().out
    assert "tg1" in out and "refGene" in out
