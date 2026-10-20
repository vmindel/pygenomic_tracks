# Claude Code skill

pygenomic_tracks ships an [Agent Skill](https://docs.claude.com/en/docs/claude-code/skills)
for [Claude Code](https://claude.com/claude-code) that teaches the agent to use the
library correctly. It covers the data model, the API, the defaults that silently change
a figure, and recipes for pymisha workflows.

## Install

Claude Code doesn't scan installed Python packages for skills, so copy the skill once
into your personal skills directory:

```bash
pygt-install-skills
```

This installs it to `~/.claude/skills/pygenomic-tracks/`, where Claude Code can use it
in every project from the next session on.

| command | effect |
|---|---|
| `pygt-install-skills` | install to `~/.claude/skills/pygenomic-tracks` |
| `pygt-install-skills --dest .claude/skills/pygenomic-tracks` | install into the current project only |
| `pygt-install-skills --force` | overwrite, e.g. after upgrading pygenomic_tracks |
| `pygt-install-skills --print-path` | print the bundled copy's location |

Once it's installed, ask for plots in plain language. For example: "plot Msn2 ChEC
signal with motif hits and genes around GAL1 in sacCer3", "make a figure of this
gextract output with the CREs coloured by class", or "why is my gene track empty?".
Claude Code reads the skill automatically.

## What's inside

The skill is a short router, `SKILL.md`, that sends the agent to one topic file under
`references/`. The agent reads that file only when the task needs it, so the skill costs
almost no context until it's used.

| file | covers |
|---|---|
| `references/inputs-coordinates.md` | required columns, 0-based coordinates, `cols=`, chromosome names, readers |
| `references/signal.md` | step drawing, binning and `agg`, NaN gaps, shared `ylim`, overlays |
| `references/intervals.md` | default colour, `color_by` (categorical vs numeric), motif hits, labels, strand |
| `references/genes-annotations.md` | `genes`, genome names, cache, offline clusters, own GTFs |
| `references/layout-figures.md` | `figure`, own axes, highlights, several loci, saving PDFs |
| `references/pymisha.md` | `gextract` → signal, `gscreen` → intervals, interval sets |
| `scripts/check_table.py` | checks a table against a region: columns, chromosome names, rows in window, coordinate sanity |

Each file emphasises footguns that give a wrong figure without an error:

- comparing signal tracks that each scale to their own maximum;
- a numeric column meant as categories;
- category colours that change between loci unless a dict palette is used;
- chromosome names that don't match (`chrM` vs `MT`);
- 1-based coordinates from GTF, R or IGV;
- resizing a figure after drawing.

The test suite checks the skill:

- every routed file exists;
- every reference file is routed;
- every `pygt.<name>` it mentions is a real export;
- the install command works;
- `check_table.py` reports correctly on good and bad tables.
