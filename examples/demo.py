"""Synthetic demo: signal, CREs, motif hits and gene models in one figure.

The tables mimic pymisha output (chrom/start/end + value columns), e.g.
    chec = pm.gextract("chec.msn2", intervals=pm.gintervals("chr2", 1_200_000, 1_350_000), iterator=50)
    hits = pm.gscreen("pwm.msn2 > 8", ...)
Run:  python examples/demo.py  ->  docs/demo.png
"""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import numpy as np
import pandas as pd

import pygenomic_tracks as pygt

rng = np.random.default_rng(0)
region = pygt.Region("chr2", 1_200_000, 1_350_000)

# 50 bp binned signal with a few peaks and one stretch of missing bins (NaN)
starts = np.arange(region.start, region.end, 50)
x = (starts - region.start) / region.length
sig = rng.gamma(1.2, 1.0, len(starts))
for centre, h in [(0.18, 40), (0.42, 25), (0.47, 18), (0.81, 55)]:
    sig += h * np.exp(-((x - centre) ** 2) / (2 * 0.004 ** 2))
sig[(x > 0.6) & (x < 0.63)] = np.nan
chec = pd.DataFrame({"chrom": "chr2", "start": starts, "end": starts + 50,
                     "chec.msn2": sig, "intervalID": 1})

cres = pd.DataFrame({
    "chrom": "chr2",
    "start": [1_222_000, 1_261_000, 1_268_000, 1_270_500, 1_318_000, 1_325_000],
    "end":   [1_228_000, 1_266_000, 1_273_000, 1_274_000, 1_324_000, 1_327_500],
    "cre":   ["promoter", "enhancer", "enhancer", "CTCF", "promoter", "enhancer"],
    "name":  ["P1", "E1", "E2", "C1", "P2", "E3"],
})

hit_pos = np.sort(rng.integers(region.start, region.end, 120))
motifs = pd.DataFrame({"chrom": "chr2", "start": hit_pos, "end": hit_pos + 12,
                       "score": rng.uniform(8, 14, len(hit_pos)),
                       "strand": rng.choice(["+", "-"], len(hit_pos))})


def transcript(tid, gene, strand, exons, cds):
    rows = [("chr2", "exon", s, e, strand, tid, gene) for s, e in exons]
    rows += [("chr2", "CDS", s, e, strand, tid, gene) for s, e in cds]
    return rows


gtf = pd.DataFrame(
    transcript("T1", "GENEA", "+",
               [(1_215_000, 1_216_500), (1_224_000, 1_225_000), (1_240_000, 1_243_000)],
               [(1_216_000, 1_216_500), (1_224_000, 1_225_000), (1_240_000, 1_241_500)])
    + transcript("T2", "GENEA", "+",
                 [(1_215_000, 1_216_500), (1_240_000, 1_243_000)],
                 [(1_216_000, 1_216_500), (1_240_000, 1_241_500)])
    + transcript("T3", "GENEB", "-",
                 [(1_290_000, 1_292_000), (1_300_000, 1_301_000), (1_330_000, 1_331_200)],
                 [(1_291_000, 1_292_000), (1_300_000, 1_301_000), (1_330_000, 1_330_600)])
    + transcript("T4", "LNC1", "+", [(1_255_000, 1_256_000), (1_262_000, 1_263_500)], []),
    columns=["chrom", "feature", "start", "end", "strand", "transcript_id", "gene_name"],
)

fig, axes = pygt.figure(region, heights=[1.2, 0.6, 0.9, 1.1])
pygt.signal(chec, axes[0], region, value="chec.msn2", label="Msn2 ChEC")
pygt.intervals(cres, axes[1], region, color_by="cre", label_by="name", legend=True,
               palette={"promoter": "#c0392b", "enhancer": "#e67e22", "CTCF": "#2c7fb8"},
               label="CREs")
pygt.intervals(motifs, axes[2], region, color_by="score", palette="Greys", vmin=6,
               strand="strand", label="Msn2 motifs")
pygt.genes(gtf, axes[3], region, label="Genes")

out = Path(__file__).resolve().parents[1] / "docs" / "demo.png"
fig.savefig(out, dpi=150, bbox_inches="tight")
print(out)
