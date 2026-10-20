# Command line

## pygt

```text
pygt download GENOME [GENOME ...] [--track TRACK] [--force]
pygt list
```

`download` fetches and prepares the gene annotation of each UCSC assembly into the cache
(see [Genes and annotations](../tutorials/genes.md)); `list` shows what is cached and
where. `python -m pygenomic_tracks` is the same command.

| option | effect |
|---|---|
| `--track` | UCSC gene track, e.g. `knownGene`, `ncbiRefSeq`, `ensGene` |
| `--force` | re-download and re-prepare |

The cache directory is `$PYGT_CACHE_DIR`, else `$XDG_CACHE_HOME/pygenomic_tracks`, else
`~/.cache/pygenomic_tracks`.

## pygt-install-skills

```text
pygt-install-skills [--dest DEST] [--force] [--print-path]
```

Copies the bundled [Claude Code skill](../claude-skill.md) to
`~/.claude/skills/pygenomic-tracks` (or `--dest`).
