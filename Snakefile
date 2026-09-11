import csv
import json
from pathlib import Path

configfile: "config/config.yaml"

ROOT = config["output_dir"]
CFG = json.dumps(config, ensure_ascii=False)
CODE = sorted(str(p) for p in Path("src/scrna_learn").glob("*.py")) + ["uv.lock", "profiles/default/config.yaml"]

# Quote configuration with Snakemake's :q formatter; never concatenate user paths into shell code.
def source_inputs():
    if config["input"]["mode"] == "demo":
        return [config["input"]["path"]]
    sample_sheet = config["input"]["samples"]
    files = [sample_sheet]
    with open(sample_sheet) as handle:
        for row in csv.DictReader((line for line in handle if not line.startswith("#")), delimiter="\t"):
            path = Path(row["path"])
            if path.is_dir():
                files.extend(str(p) for p in sorted(path.iterdir()) if p.is_file())
            else:
                files.append(str(path))
    return files

rule all:
    input:
        f"{ROOT}/report.html",
        f"{ROOT}/provenance/artifact_checksums.json",
        f"{ROOT}/objects/06_annotated.h5ad"

if config["input"]["mode"] == "demo":
    rule download:
        input: "src/scrna_learn/pbmc3k.sha256", "resources/pbmc3k_counts.h5ad", "resources/pbmc3k_counts.source.json"
        output: config["input"]["path"]
        log: f"{ROOT}/logs/download.log"
        shell: "scrna-learn download --out {output:q} > {log:q} 2>&1"

rule ingest:
    input: raw=source_inputs(), code=CODE
    output:
        f"{ROOT}/objects/01_counts.h5ad",
        f"{ROOT}/provenance/inputs.json"
    params: cfg=CFG
    resources: mem_mb=2000
    log: f"{ROOT}/logs/01_ingest.log"
    shell: "scrna-learn stage ingest --config-json {params.cfg:q} > {log:q} 2>&1"

rule qc:
    input: f"{ROOT}/objects/01_counts.h5ad", CODE
    output:
        f"{ROOT}/objects/02_qc.h5ad",
        f"{ROOT}/tables/cell_qc.csv",
        f"{ROOT}/tables/qc_summary.json",
        f"{ROOT}/figures/qc.png"
    params: cfg=CFG
    resources: mem_mb=4000
    log: f"{ROOT}/logs/02_qc.log"
    shell: "scrna-learn stage qc --config-json {params.cfg:q} > {log:q} 2>&1"

rule normalize:
    input: f"{ROOT}/objects/02_qc.h5ad", CODE
    output:
        f"{ROOT}/objects/03_normalized.h5ad",
        f"{ROOT}/tables/gene_metrics.csv",
        f"{ROOT}/figures/highly_variable_genes.png"
    params: cfg=CFG
    resources: mem_mb=3000
    log: f"{ROOT}/logs/03_normalize.log"
    shell: "scrna-learn stage normalize --config-json {params.cfg:q} > {log:q} 2>&1"

rule embed:
    input: f"{ROOT}/objects/03_normalized.h5ad", CODE
    output:
        f"{ROOT}/objects/04_clustered.h5ad",
        f"{ROOT}/tables/cluster_counts.csv",
        f"{ROOT}/figures/pca_variance.png",
        f"{ROOT}/figures/umap_clusters.png"
    params: cfg=CFG
    resources: mem_mb=4000
    log: f"{ROOT}/logs/04_embed.log"
    shell: "scrna-learn stage embed --config-json {params.cfg:q} > {log:q} 2>&1"

rule markers:
    input: f"{ROOT}/objects/04_clustered.h5ad", CODE
    output:
        f"{ROOT}/objects/05_markers.h5ad",
        f"{ROOT}/tables/markers_all.csv",
        f"{ROOT}/tables/markers_positive.csv",
        f"{ROOT}/tables/markers_top10.csv"
    params: cfg=CFG
    resources: mem_mb=4000
    log: f"{ROOT}/logs/05_markers.log"
    shell: "scrna-learn stage markers --config-json {params.cfg:q} > {log:q} 2>&1"

rule annotate:
    input:
        f"{ROOT}/objects/05_markers.h5ad", CODE,
        config["annotation"]["markers"],
        ([config["annotation"]["manual_labels"]] if config["annotation"]["manual_labels"] else [])
    output:
        f"{ROOT}/objects/06_annotated.h5ad",
        f"{ROOT}/tables/annotation_evidence.csv",
        f"{ROOT}/tables/annotation_expression.csv",
        f"{ROOT}/tables/annotation_scores.csv",
        f"{ROOT}/tables/cell_metadata.csv",
        f"{ROOT}/tables/celltype_counts.csv",
        f"{ROOT}/figures/umap_celltypes.png",
        f"{ROOT}/figures/marker_dotplot.png"
    params: cfg=CFG
    resources: mem_mb=3000
    log: f"{ROOT}/logs/06_annotate.log"
    shell: "scrna-learn stage annotate --config-json {params.cfg:q} > {log:q} 2>&1"

rule enrichment:
    input:
        f"{ROOT}/objects/05_markers.h5ad", f"{ROOT}/tables/markers_positive.csv", CODE,
        ([config["enrichment"]["gmt"]] if config["enrichment"]["gmt"] else [])
    output:
        f"{ROOT}/tables/enrichment.csv",
        f"{ROOT}/tables/enrichment_top.csv",
        f"{ROOT}/tables/enrichment_status.json"
    params: cfg=CFG
    log: f"{ROOT}/logs/07_enrichment.log"
    shell: "scrna-learn stage enrichment --config-json {params.cfg:q} > {log:q} 2>&1"

rule report:
    input:
        rules.ingest.output, rules.qc.output, rules.normalize.output, rules.embed.output,
        rules.markers.output, rules.annotate.output, rules.enrichment.output, CODE
    output:
        f"{ROOT}/report.html",
        f"{ROOT}/run_summary.json",
        f"{ROOT}/provenance/run_manifest.json",
        f"{ROOT}/provenance/config_used.yaml",
        f"{ROOT}/provenance/artifact_checksums.json"
    params: cfg=CFG
    resources: mem_mb=2000
    log: f"{ROOT}/logs/08_report.log"
    shell: "scrna-learn stage report --config-json {params.cfg:q} > {log:q} 2>&1"
