"""Tutorial code must run, print what the docs say, and draw the figures it links."""

import importlib.util
from pathlib import Path

import pytest

DOCS = Path(__file__).parent.parent / "docs"
spec = importlib.util.spec_from_file_location("run_docs", DOCS / "run_docs.py")
run_docs = importlib.util.module_from_spec(spec)
spec.loader.exec_module(run_docs)


@pytest.mark.parametrize("page", sorted((DOCS / "tutorials").glob("*.md")), ids=lambda p: p.stem)
def test_tutorial(page):
    errors = run_docs.check_page(page)
    assert not errors, "\n\n".join(errors) + "\n\nRun `python docs/run_docs.py --update` if intended."


def test_docs_mention_only_real_api():
    import re

    import pygenomic_tracks as pygt

    names = set()
    for page in DOCS.rglob("*.md"):
        names |= set(re.findall(r"\bpygt\.([A-Za-z_]\w*)", page.read_text()))
    missing = sorted(n for n in names if not hasattr(pygt, n))
    assert not missing, f"docs mention pygt.{missing}, which is not exported"
