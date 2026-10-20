import re
import subprocess
import sys

import pandas as pd
import pytest

import pygenomic_tracks as pygt
from pygenomic_tracks._skills.install import bundled_skill_dir, install_skill, main

SKILL = bundled_skill_dir()
CHECK = SKILL / "scripts" / "check_table.py"


def test_skill_frontmatter_and_routes():
    text = (SKILL / "SKILL.md").read_text()
    front = re.match(r"^---\n(.*?)\n---\n", text, re.S)
    assert front, "SKILL.md needs YAML frontmatter"
    fields = dict(line.split(": ", 1) for line in front.group(1).splitlines())
    assert fields["name"] == "pygenomic-tracks"
    assert 100 < len(fields["description"]) <= 1024
    routed = set(re.findall(r"`((?:references|scripts)/[\w.-]+)`", text))
    assert routed, "router table is empty"
    for rel in routed:
        assert (SKILL / rel).is_file(), f"SKILL.md routes to missing {rel}"
    on_disk = {f"references/{p.name}" for p in (SKILL / "references").glob("*.md")}
    assert on_disk <= routed, f"reference files not routed from SKILL.md: {on_disk - routed}"


def test_skill_mentions_only_real_api():
    names = set()
    for path in [SKILL / "SKILL.md", *(SKILL / "references").glob("*.md")]:
        text = path.read_text()
        names |= set(re.findall(r"\bpygt\.([A-Za-z_]\w*)", text))
        for ref in re.findall(r"`(references/[\w.-]+)`", text):
            assert (SKILL / ref).is_file(), f"{path.name} points to missing {ref}"
    missing = sorted(n for n in names if not hasattr(pygt, n))
    assert not missing, f"skill mentions pygt.{missing}, which is not exported"


def test_install(tmp_path, capsys):
    dest = install_skill(tmp_path / "skills" / "pygenomic-tracks")
    assert (dest / "SKILL.md").read_text() == (SKILL / "SKILL.md").read_text()
    assert (dest / "references" / "signal.md").is_file()
    assert (dest / "scripts" / "check_table.py").is_file()
    with pytest.raises(FileExistsError):
        install_skill(dest)
    install_skill(dest, force=True)
    assert main(["--dest", str(dest)]) == 1
    assert main(["--dest", str(dest), "--force"]) == 0
    assert main(["--print-path"]) == 0
    assert capsys.readouterr().out.strip().endswith(str(SKILL))


def _run(*args):
    return subprocess.run([sys.executable, str(CHECK), *map(str, args)],
                          capture_output=True, text=True)


def test_check_table_ok(tmp_path):
    p = tmp_path / "signal.tsv"
    pygt.example_signal().to_csv(p, sep="\t", index=False)
    out = _run(p, "--region", "chr2:1,200,000-1,350,000", "--value", "signal")
    assert out.returncode == 0, out.stdout + out.stderr
    assert "rows overlap chr2" in out.stdout and "NaN rows" in out.stdout


def test_check_table_problems(tmp_path):
    p = tmp_path / "peaks.bed"
    pygt.example_peaks()[["chrom", "start", "end", "name"]].to_csv(p, sep="\t", header=False,
                                                                   index=False)
    out = _run(p, "--region", "chrX:0-1000")
    assert out.returncode == 1 and "no rows on 'chrX'" in out.stdout

    out = _run(p, "--region", "chr2:0-1000")
    assert out.returncode == 1 and "0 rows overlap" in out.stdout

    p2 = tmp_path / "ensembl.tsv"
    pygt.example_peaks().assign(chrom="2").to_csv(p2, sep="\t", index=False)
    out = _run(p2, "--region", "chr2:1,200,000-1,350,000")
    assert out.returncode == 0 and "matches them automatically" in out.stdout

    out = _run(p2, "--region", "chr2:1,200,000-1,350,000", "--value", "nope")
    assert out.returncode == 1 and "missing column" in out.stdout


def test_check_table_one_based_warning(tmp_path):
    df = pd.DataFrame({"chrom": "chr1", "start": [1, 101, 201], "end": [100, 200, 300]})
    p = tmp_path / "onebased.tsv"
    df.to_csv(p, sep="\t", index=False)
    out = _run(p, "--region", "chr1:0-1000")
    assert "1-based" in out.stdout
