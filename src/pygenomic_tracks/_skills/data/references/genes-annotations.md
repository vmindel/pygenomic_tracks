# Gene models and annotations

## Drawing (pygt.genes)

```python
pygt.genes(df_or_genome, ax, region, feature="feature", transcript="transcript_id",
           name="gene_name", strand="strand", transcripts="all", color="0.2",
           height=0.6, utr_ratio=0.5, show_names=True, pack=True, fontsize=7,
           label=None)
```

- First argument: a UCSC assembly name (`"hg38"`), a GTF path, or a
  DataFrame with **one row per exon/CDS** (long format). A table with one
  row per gene cannot draw exons; use `pygt.intervals` for it instead.
- Only `feature == "exon"` and `"CDS"` rows are drawn. A transcript with
  CDS rows draws its non-coding exon parts at `utr_ratio` of full height.
  One without CDS (lncRNA) is drawn at full height throughout.
- `transcripts="all"` draws every isoform. RefSeq/GENCODE loci then become
  tall stacks, so use `transcripts="longest"` (one per `gene_name`, longest
  span) for most figures.
- The transcript span comes from its exons on the region's chromosome.
  Transcripts crossing the window edge are clipped.

## Annotations by genome name

```python
ann = pygt.annotation("mm10")                         # DataFrame, cached
region = pygt.gene_region("mm10", "Myc", flank=50_000)
pygt.genes("mm10", ax, region, transcripts="longest")
```

- Source: UCSC `goldenPath/<genome>/bigZips/genes/<genome>.<track>.gtf.gz`.
  Any UCSC assembly name works. Tracks are tried in order
  `ncbiRefSeq, refGene, ensGene, knownGene`; force one with
  `track="knownGene"` (GENCODE on hg38/mm39).
- UCSC chromosome names (`chr1`, yeast `chrI`) match misha databases built
  from UCSC FASTA.
- Gene labels: `gene_name` from the GTF, falling back to `gene_id`. With an
  `ensGene` track, names may be Ensembl/systematic ids (yeast `YAL001C`),
  not symbols.
- `pygt.gene_region` matches gene name or id exactly, then
  case-insensitively by name. It warns if the gene is on several
  chromosomes (PAR genes) and uses the one with most rows.
- `pygt.annotation` returns the same in-memory object on every call.
  **Copy before modifying** (`ann = pygt.annotation("hg38").copy()`).

## Cache and clusters

- Cache directory: `pygt.cache_dir()`, which is `$PYGT_CACHE_DIR`, else
  `~/.cache/pygenomic_tracks`. It holds the raw GTF and a prepared pickle.
  hg38 is tens of MB.
- **Compute nodes often have no internet** (LSF, SLURM). Fill the cache on a
  login node:

  ```bash
  pygt download hg38 mm10 sacCer3
  pygt list
  ```

  In jobs, call `pygt.annotation("hg38", download_missing=False)`, so a
  missing genome raises at once with the command to run instead of
  hanging on a download.
- Share one cache per lab with `export PYGT_CACHE_DIR=/shared/path`. Writes
  are atomic (temp file + rename), so parallel jobs don't corrupt it.
- A pickle written by another pandas version that fails to load is rebuilt
  from the raw GTF automatically. `pygt download hg38 --force` refreshes
  everything.
- A network error (proxy 403, DNS) raises `URLError`. "No GTF for genome"
  raises `FileNotFoundError` naming the UCSC directory to check.

## Own GTF

```python
gtf = pygt.read_gtf("gencode.vM25.annotation.gtf.gz")   # exon + CDS, 0-based
pygt.genes(gtf, ax, region, transcripts="longest")
```

Reading a full GENCODE GTF takes tens of seconds; save the result with
`gtf.to_parquet(...)` once.
