from __future__ import annotations

import base64
import html
import importlib.metadata
import json
import platform
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import scanpy as sc
import yaml

from .common import outpaths, sha256, write_json


def report(cfg):
    root = outpaths(cfg)
    a = sc.read_h5ad(root / "objects/06_annotated.h5ad")
    qc = json.loads((root / "tables/qc_summary.json").read_text())
    enrichment = json.loads((root / "tables/enrichment_status.json").read_text())
    versions = {name: importlib.metadata.version(name) for name in [
        "scrna-seq-learning", "scanpy", "anndata", "numpy", "pandas", "scipy", "numba", "igraph", "snakemake", "pydeseq2"]}
    try:
        revision = subprocess.check_output(["git", "rev-parse", "HEAD"], stderr=subprocess.DEVNULL, text=True).strip()
        dirty = bool(subprocess.check_output(["git", "status", "--porcelain"], text=True).strip())
    except (subprocess.SubprocessError, FileNotFoundError):
        revision, dirty = "uncommitted", True
    manifest = {"created_at_utc": datetime.now(timezone.utc).isoformat(), "python": platform.python_version(),
                "platform": platform.platform(), "versions": versions, "config": cfg,
                "git_commit": revision, "git_dirty": dirty,
                "inputs": json.loads((root / "provenance/inputs.json").read_text()),
                "lock_sha256": sha256("uv.lock") if Path("uv.lock").is_file() else None,
                "source_sha256": {str(p): sha256(p) for p in sorted(Path(__file__).parent.glob("*.py"))},
                "workflow_sha256": {str(p): sha256(p) for p in [Path("Snakefile"), Path("profiles/default/config.yaml")]
                                    if p.is_file()}}
    write_json(root / "provenance/run_manifest.json", manifest)
    (root / "provenance/config_used.yaml").write_text(yaml.safe_dump(cfg, allow_unicode=True, sort_keys=False))
    summary = {**qc, "clusters": int(a.obs.cluster.nunique()), "hvg": int(a.var.highly_variable.sum()),
               "samples": int(a.obs.sample_id.nunique()), "donors": int(a.obs.donor.nunique()),
               "cluster_sizes": {str(k): int(v) for k, v in a.obs.cluster.value_counts().sort_index().items()},
               "annotation_status": sorted(a.obs.annotation_status.unique().tolist()),
               "cell_type_counts": {str(k): int(v) for k, v in a.obs.cell_type.value_counts().items()},
               "qc_per_sample": qc["samples"]}
    write_json(root / "run_summary.json", summary)
    def table(path, n=None):
        frame = pd.read_csv(root / path)
        if n:
            frame = frame.head(n)
        return '<div class="table">' + frame.to_html(index=False, border=0, float_format=lambda v: f"{v:.3g}", escape=True) + '</div>'
    def figure(name, caption):
        p = root / "figures" / name
        if not p.exists():
            return ""
        image = base64.b64encode(p.read_bytes()).decode()
        return f'<figure><img src="data:image/png;base64,{image}" alt="{html.escape(caption)}"><figcaption>{html.escape(caption)}</figcaption></figure>'
    reviewed = a.obs.annotation_status.eq("user_reviewed").all()
    dataset_label = ("10x Genomics PBMC3k（CC BY 4.0；通过 Scanpy raw-count 镜像获取）"
                     if cfg["input"]["mode"] == "demo" else "自定义输入；样本和文件校验信息见 provenance/inputs.json")
    status = "已使用人工复核标签" if reviewed else "当前细胞类型为候选注释，需人工复核"
    stats = [("输入细胞", a0 := qc["input_cells"]), ("保留细胞", a.n_obs), ("保留基因", a.n_vars), ("Leiden clusters", a.obs.cluster.nunique())]
    cards = "".join(f'<div class="card"><span>{label}</span><strong>{value:,}</strong></div>' for label, value in stats)
    body = f'''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>{html.escape(cfg['project'])} · scRNA-seq analysis report</title>
<style>body{{margin:0;background:#f3f6f6;color:#183335;font-family:system-ui,-apple-system,"Microsoft YaHei",sans-serif;line-height:1.75}}main{{max-width:1150px;margin:auto;padding:44px 26px}}header{{border-top:5px solid #19857d;padding-top:20px}}.eyebrow{{letter-spacing:.14em;font-size:12px;color:#547775}}h1{{font-size:38px;margin:8px 0}}h2{{font-size:23px;margin:36px 0 14px}}.cards{{display:grid;grid-template-columns:repeat(4,1fr);gap:16px;margin:26px 0}}.card,section{{background:white;padding:22px;border-radius:12px;box-shadow:0 1px 3px #18333510}}.card span{{display:block;color:#57716f;font-size:13px}}.card strong{{font-size:34px}}.note{{background:#e4f1ee;border-left:4px solid #19857d;padding:15px 20px}}.table{{overflow:auto}}table{{border-collapse:collapse;font-size:13px;width:100%}}td,th{{padding:9px 12px;border-bottom:1px solid #dce7e5;text-align:left;white-space:nowrap}}th{{color:#246f69;background:#f0f6f5}}figure{{margin:20px 0}}img{{max-width:100%;height:auto}}figcaption{{font-size:13px;color:#536c6b}}pre{{white-space:pre-wrap;font-size:12px;background:#f5f8f8;padding:18px}}footer{{font-size:12px;color:#647574;padding:24px 0}}@media(max-width:680px){{.cards{{grid-template-columns:repeat(2,1fr)}}h1{{font-size:27px}}}}</style>
<main><header><div class="eyebrow">SCRNA-SEQ LEARNING / REPRODUCIBLE ANALYSIS</div><h1>{html.escape(cfg['project'])} 分析报告</h1><p>从原始 UMI 计数到质控、聚类、marker 和细胞注释。</p></header>
<div class="cards">{cards}</div><p class="note">{status}。{a.obs.donor.nunique()} 位供者；cluster marker 属于探索性比较，不能据此推断疾病组间差异。</p>
<h2>01 · 质控与双细胞检测</h2><section><p>保留率 {a.n_obs / a0:.1%}。阈值与逐细胞去留原因见配置和 cell_qc.csv。双细胞标记是算法预测。</p>
<div class="table">{pd.DataFrame(qc['samples']).to_html(index=False, border=0, escape=True)}</div>{figure('qc.png','基因数、UMI 数与线粒体计数比例：过滤前和保留细胞。')}
{''.join(figure(str(p.relative_to(root / 'figures')), 'Scrublet 实测和模拟得分分布；竖线为当前阈值。') for p in sorted((root / 'figures').glob('doublets_*.png')))}</section>
<h2>02 · 特征选择与聚类</h2><section><p>使用 {int(a.var.highly_variable.sum())} 个高变基因构建 PCA 与邻居图。UMAP 距离不等同于发育时间。</p>
{figure('pca_variance.png','主成分解释方差：用于检查选择的维度。')}{figure('umap_clusters.png','Leiden 聚类与线粒体比例。')}{figure('tsne_clusters.png','同一 PCA 表征的 t-SNE 可视化。')}{figure('umap_samples.png','检查样本、条件与批次的分布。')}</section>
<h2>03 · 注释证据</h2><section>{figure('umap_celltypes.png','细胞类型标签：未经人工复核时为候选标签。')}{table('tables/annotation_evidence.csv')}
{figure('marker_dotplot.png','点大小代表表达细胞比例，颜色为按基因缩放的平均 log 表达。')}{figure('feature_markers.png','经典 marker 的 UMAP 表达分布。')}</section>
<h2>04 · Marker 与细胞组成</h2><section><p>下表展示前 30 行；完整结果在 markers_all.csv / markers_top10.csv。近似 log2FC 来自 log-normalized 表达。校正 p 值按每个 cluster 的基因检验计算。</p>{table('tables/markers_top10.csv',30)}{table('tables/celltype_counts.csv')}</section>
<h2>05 · 功能富集与分析边界</h2><section>{figure('enrichment.png','以全部受检基因为背景的 cluster marker ORA；每个 cluster 展示至多 5 个条目。')}<pre>{html.escape(json.dumps(enrichment,ensure_ascii=False,indent=2))}</pre><p>默认 PBMC3k 为单供者健康样本，不运行疾病差异表达、批次整合、拟时序或细胞通讯推断。多样本统计入口和适用条件见中文教程。</p></section>
<h2>06 · 可复现信息</h2><section><p>所有中间对象保存在 objects/。X 为 log-normalized 表达，layers["counts"] 为未归一化的 UMI 计数。PCA 缩放只发生在高变基因的独立对象中。</p><pre>{html.escape(json.dumps({'versions':versions,'config':cfg},ensure_ascii=False,indent=2))}</pre></section>
<footer>数据：{html.escape(dataset_label)}。来源、输入 SHA256 和代码版本见 provenance/run_manifest.json。生成时间 {manifest['created_at_utc']}。</footer></main></html>'''
    (root / "report.html").write_text(body)
    artifacts = {str(p.relative_to(root)): sha256(p) for p in sorted(root.rglob("*"))
                 if p.is_file() and p.name != "artifact_checksums.json" and "logs" not in p.relative_to(root).parts
                 and not p.name.startswith("timing")}
    write_json(root / "provenance/artifact_checksums.json", artifacts)
