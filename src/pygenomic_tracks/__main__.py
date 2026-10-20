"""Command line: ``python -m pygenomic_tracks download hg38 mm10`` or ``pygt list``."""
from __future__ import annotations

import argparse
import sys

from .annotations import cache_dir, cached, download


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(prog="pygt", description="Manage cached gene annotations.")
    sub = parser.add_subparsers(dest="cmd", required=True)

    dl = sub.add_parser("download", help="download and prepare annotations for UCSC assemblies")
    dl.add_argument("genomes", nargs="+", help="UCSC assembly names, e.g. hg38 mm10 sacCer3")
    dl.add_argument("--track", help="UCSC gene track, e.g. knownGene, ncbiRefSeq, ensGene")
    dl.add_argument("--force", action="store_true", help="re-download and re-prepare")

    sub.add_parser("list", help="show cached annotations")

    args = parser.parse_args(argv)
    if args.cmd == "download":
        for genome in args.genomes:
            print(download(genome, track=args.track, force=args.force))
    else:
        table = cached()
        print(f"cache: {cache_dir()}")
        print(table.to_string(index=False) if len(table) else "(empty)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
