from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

# Bound native libraries before importing numerical dependencies; avoid oversubscription.
for _key in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"]:
    os.environ[_key] = "1"
os.environ["MPLBACKEND"] = "Agg"

from . import analysis
from .common import load_config, outpaths, write_json, validate_config
from .data import ingest, download_demo
from .extensions import enrichment, pseudobulk_de
from .report import report

STAGES = {"ingest": ingest, "qc": analysis.qc, "normalize": analysis.normalize,
          "embed": analysis.embed, "markers": analysis.markers, "annotate": analysis.annotate,
          "enrichment": enrichment, "report": report}


def run_stage(stage, cfg):
    root = outpaths(cfg)
    started = time.perf_counter()
    print(f"[{stage}] start", flush=True)
    STAGES[stage](cfg)
    elapsed = round(time.perf_counter() - started, 3)
    # Separate files avoid a concurrent update race when Snakemake runs branches.
    write_json(root / f"provenance/timing_{stage}.json", {"stage": stage, "seconds": elapsed})
    print(f"[{stage}] completed in {elapsed:.1f}s", flush=True)


def main():
    parser = argparse.ArgumentParser(description="中文教学 scRNA-seq pipeline")
    sub = parser.add_subparsers(dest="command", required=True)
    download_parser = sub.add_parser("download")
    download_parser.add_argument("--out", default="data/raw/pbmc3k_raw.h5ad")
    for name in ["run", "stage"]:
        p = sub.add_parser(name)
        if name == "stage":
            p.add_argument("stage", choices=STAGES)
        options = p.add_mutually_exclusive_group()
        options.add_argument("--config", default="config/config.yaml")
        options.add_argument("--config-json", help=argparse.SUPPRESS)
    p = sub.add_parser("pseudobulk")
    p.add_argument("--input", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--cell-type", required=True)
    p.add_argument("--case", required=True)
    p.add_argument("--control", required=True)
    p.add_argument("--min-cells", type=int, default=20)
    p.add_argument("--paired", action="store_true")
    args = parser.parse_args()
    if args.command == "download":
        download_demo(args.out)
        return
    if args.command == "pseudobulk":
        pseudobulk_de(args.input, args.out, args.cell_type, args.case, args.control, args.min_cells, args.paired)
        return
    cfg = validate_config(json.loads(args.config_json)) if args.config_json else load_config(args.config)
    for stage in (STAGES if args.command == "run" else [args.stage]):
        run_stage(stage, cfg)


if __name__ == "__main__":
    main()
