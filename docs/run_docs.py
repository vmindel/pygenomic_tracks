"""Execute the tutorial code, check its printed output and (re)draw its figures.

Convention in docs/tutorials/*.md:

* ```` ```python ```` blocks are executed in order, one shared namespace per
  page, inside a temporary working directory, with the Agg backend.
* A ```` ```text ```` block directly after a python block is its expected
  stdout.
* Image links ``![...](../images/<name>.png)`` between a python block and the
  next python block are the figures that block draws, in order of creation.
  Figures are closed after every block.
* ```` ```py ````, ```` ```bash ```` and other blocks are shown but never
  executed (for things needing a network, a misha database, a cluster...).

    python docs/run_docs.py            # check all tutorials, exit 1 on problems
    python docs/run_docs.py --update   # rewrite text blocks and images from fresh output
"""

from __future__ import annotations

import contextlib
import io
import os
import re
import sys
import tempfile
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

TUTORIALS = Path(__file__).parent / "tutorials"
FENCE = re.compile(r"^```(\w*)[^\n]*\n(.*?)^```[ \t]*\n", re.S | re.M)
IMAGE = re.compile(r"!\[[^\]]*\]\(([^)\s]+\.png)\)")


def _blocks(text: str):
    """Yield (lang, body, start, end) for every fenced block."""
    for m in FENCE.finditer(text):
        yield m.group(1), m.group(2), m.start(), m.end()


def run_page(path: Path, save_images: bool = False):
    """Run a page.

    Returns one dict per python block: code, expected stdout (or None), actual
    stdout, the span of the expected block, and the image paths it should draw
    with the number of figures it did draw.
    """
    import numpy as np
    import pandas as pd

    pd.set_option("display.width", 100)
    pd.set_option("display.max_columns", 20)
    np.set_printoptions(precision=4, suppress=True, linewidth=100)
    plt.rcParams.update({"figure.dpi": 100, "savefig.dpi": 110})

    text = path.read_text()
    blocks = list(_blocks(text))
    py_idx = [i for i, b in enumerate(blocks) if b[0] == "python"]
    results = []
    namespace: dict = {"__name__": "__docs__"}
    cwd = os.getcwd()
    with tempfile.TemporaryDirectory() as tmp:
        os.chdir(tmp)
        try:
            for k, i in enumerate(py_idx):
                _, body, _, end = blocks[i]
                expected, span = None, None
                if i + 1 < len(blocks):
                    nlang, nbody, nstart, nend = blocks[i + 1]
                    if nlang == "text" and text[end:nstart].strip() == "":
                        expected, span = nbody, (nstart, nend)
                stop = blocks[py_idx[k + 1]][2] if k + 1 < len(py_idx) else len(text)
                images = [(path.parent / m).resolve() for m in IMAGE.findall(text[end:stop])]

                before = set(plt.get_fignums())
                buf = io.StringIO()
                with contextlib.redirect_stdout(buf), warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    exec(compile(body, f"{path.name}:block{k}", "exec"), namespace)
                new = [n for n in plt.get_fignums() if n not in before]
                if save_images:
                    for img, num in zip(images, new):
                        img.parent.mkdir(parents=True, exist_ok=True)
                        plt.figure(num).savefig(img, bbox_inches="tight")
                plt.close("all")
                results.append({"code": body, "expected": expected, "actual": buf.getvalue(),
                                "span": span, "images": images, "n_figures": len(new)})
        finally:
            os.chdir(cwd)
    return results


def _norm(s: str | None) -> str:
    return "\n".join(line.rstrip() for line in (s or "").strip("\n").splitlines())


def check_page(path: Path) -> list[str]:
    errors = []
    for r in run_page(path):
        first = r["code"].strip().splitlines()[0]
        where = f"{path.name}: block starting {first!r}"
        if r["expected"] is not None and _norm(r["expected"]) != _norm(r["actual"]):
            errors.append(f"{where}: output changed:\n--- expected\n{_norm(r['expected'])}\n"
                          f"--- actual\n{_norm(r['actual'])}")
        elif r["expected"] is None and r["actual"].strip():
            errors.append(f"{where}: prints output but has no text block")
        if len(r["images"]) != r["n_figures"]:
            errors.append(f"{where}: draws {r['n_figures']} figure(s) but links {len(r['images'])} image(s)")
        for img in r["images"]:
            if not img.is_file():
                errors.append(f"{where}: image {img.name} missing; run with --update")
    return errors


def update_page(path: Path) -> None:
    text = path.read_text()
    out = []
    pos = 0
    for r in run_page(path, save_images=True):
        if r["span"] is None:
            continue
        out.append(text[pos : r["span"][0]])
        out.append("```text\n" + _norm(r["actual"]) + "\n```\n")
        pos = r["span"][1]
    out.append(text[pos:])
    path.write_text("".join(out))


def main(argv: list[str]) -> int:
    pages = sorted(TUTORIALS.glob("*.md"))
    if "--update" in argv:
        for p in pages:
            update_page(p)
            print(f"updated {p.name}")
        return 0
    errors = [e for p in pages for e in check_page(p)]
    for e in errors:
        print(e, end="\n\n")
    print(f"{len(pages)} pages, {len(errors)} problems")
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
