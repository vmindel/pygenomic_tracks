"""The README figure: signal, overlay, CREs, motif hits and gene models.

Uses the bundled synthetic example locus. With pymisha the tables would come from
    chec = pm.gextract("chec.msn2", intervals=scope, iterator=50)
    hits = pm.gscreen("pwm.msn2 > 8", intervals=scope)
Run:  python examples/demo.py  ->  docs/images/demo.png
"""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import pygenomic_tracks as pygt

region = pygt.example_region()
cre_colours = {"promoter": "#c0392b", "enhancer": "#e67e22", "CTCF": "#2c7fb8"}

fig, axes = pygt.figure(region, heights=[1.2, 0.6, 0.9, 1.1])
pygt.signal(pygt.example_signal("chec"), axes[0], region, alpha=0.8, name="Msn2 ChEC",
            label="ChEC-seq")
pygt.signal(pygt.example_signal("control"), axes[0], region, color="0.3", fill=False,
            name="free MNase")
pygt.intervals(pygt.example_peaks(), axes[1], region, color_by="cre", palette=cre_colours,
               label_by="name", legend=True, label="CREs")
pygt.intervals(pygt.example_motif_hits(), axes[2], region, color_by="score",
               palette="Greys", vmin=6, label="Motif hits")
pygt.genes(pygt.example_gtf_path(), axes[3], region, label="Genes")

out = Path(__file__).resolve().parents[1] / "docs" / "images" / "demo.png"
fig.savefig(out, dpi=150, bbox_inches="tight")
print(out)
