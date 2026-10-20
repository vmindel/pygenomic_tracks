# Installation

pygenomic_tracks needs Python 3.10 or newer and depends only on numpy, pandas and
matplotlib.

```bash
pip install git+https://github.com/vmindel/pygenomic_tracks
```

From a clone:

```bash
git clone https://github.com/vmindel/pygenomic_tracks
cd pygenomic_tracks
pip install -e ".[test]"
pytest
```

Optional extras:

| extra | adds | for |
|---|---|---|
| `bigwig` | pyBigWig | `pygt.from_bigwig` |
| `test` | pytest | running the test suite |
| `docs` | mkdocs, mkdocs-material, mkdocstrings | building these docs |

```bash
pip install -e ".[bigwig]"
```

pymisha is not a dependency. Install it separately when you work with misha
databases.

## Commands

Installing the package adds two commands:

| command | does |
|---|---|
| `pygt download hg38 mm10` | download and prepare gene annotations into the cache |
| `pygt list` | show cached annotations |
| `pygt-install-skills` | install the [Claude Code skill](claude-skill.md) |

`python -m pygenomic_tracks ...` is equivalent to `pygt ...`.

## On a cluster

- **Annotations:** compute nodes often have no internet, so run `pygt download <genome>`
  on a login node once. See [Genes and annotations](tutorials/genes.md#clusters-without-internet).
- **Cache location:** the cache defaults to `~/.cache/pygenomic_tracks`. If your home
  quota is small, or to share one cache across the lab, set `PYGT_CACHE_DIR`.
- **No display:** batch jobs have no display. Matplotlib picks a non-interactive
  backend automatically when no display is available. Set `MPLBACKEND=Agg` if a job
  tries to open a window.

## Building these docs

```bash
pip install -e ".[docs]"
mkdocs serve                       # live preview at http://127.0.0.1:8000
python docs/run_docs.py            # check tutorial outputs and figures (also run by pytest)
python docs/run_docs.py --update   # refresh outputs and figures after an intended change
```
