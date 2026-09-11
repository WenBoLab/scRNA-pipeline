#!/usr/bin/env python3
"""Build or execute a modern 10x GEX Cell Ranger count job and hand off its counts."""
import argparse
import csv
import re
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

import yaml


def build_command(args):
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]*", args.sample_id):
        raise ValueError("sample-id may contain letters, digits, underscore and dash")
    if args.cores < 1 or args.memory_gb < 1:
        raise ValueError("cores and memory-gb must be positive")
    return ["cellranger", "count", f"--id={args.sample_id}",
            f"--fastqs={Path(args.fastqs).resolve()}", f"--sample={args.fastq_sample}",
            f"--transcriptome={Path(args.transcriptome).resolve()}",
            f"--create-bam={'true' if args.create_bam else 'false'}",
            f"--localcores={args.cores}", f"--localmem={args.memory_gb}"]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ["sample-id", "fastqs", "fastq-sample", "transcriptome", "donor", "batch"]:
        p.add_argument("--" + key, required=True)
    p.add_argument("--condition", default="unspecified")
    p.add_argument("--cores", type=int, default=8)
    p.add_argument("--memory-gb", type=int, default=64)
    p.add_argument("--outdir", default="results/cellranger")
    p.add_argument("--create-bam", action="store_true")
    p.add_argument("--execute", action="store_true", help="Run count; without this flag only display the command")
    p.add_argument("--run-analysis", action="store_true", help="After count, execute the teaching count-to-report pipeline")
    p.add_argument("--config", default="config/config.yaml")
    args = p.parse_args()
    cmd = build_command(args)
    print("Working directory:", Path(args.outdir).resolve())
    print(shlex.join(cmd))
    if not args.execute:
        print("Command preview only. Add --execute after checking sample chemistry and reference.")
        return
    if not shutil.which("cellranger"):
        p.error("cellranger is not on PATH; install it from the official 10x website")
    reference, fastqs = Path(args.transcriptome), Path(args.fastqs)
    if not reference.is_dir() or not (reference / "reference.json").is_file():
        p.error("transcriptome must be a Cell Ranger reference directory containing reference.json")
    if not fastqs.is_dir() or not any(fastqs.rglob("*.fastq.gz")):
        p.error("No .fastq.gz files found")
    root = Path(args.outdir).resolve(); root.mkdir(parents=True, exist_ok=True)
    subprocess.run(cmd, cwd=root, check=True)
    matrix = root / args.sample_id / "outs/filtered_feature_bc_matrix.h5"
    if not matrix.is_file():
        raise RuntimeError("Cell Ranger finished but expected filtered_feature_bc_matrix.h5 is missing")
    sheet = root / f"{args.sample_id}.samples.tsv"
    with sheet.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample_id", "path", "format", "donor", "condition", "batch"])
        writer.writerow([args.sample_id, str(matrix), "10x_h5", args.donor, args.condition, args.batch])
    cfg = yaml.safe_load(Path(args.config).read_text())
    cfg.update(project=args.sample_id, output_dir=str(root / f"{args.sample_id}_analysis"))
    cfg["input"] = {"mode": "samples", "samples": str(sheet)}
    cfg_path = root / f"{args.sample_id}.analysis.yaml"
    cfg_path.write_text(yaml.safe_dump(cfg, allow_unicode=True, sort_keys=False))
    analysis_cmd = [sys.executable, "-m", "scrna_learn.cli", "run", "--config", str(cfg_path)]
    print("Downstream command:", shlex.join(analysis_cmd))
    if args.run_analysis:
        subprocess.run(analysis_cmd, check=True)


if __name__ == "__main__":
    main()

