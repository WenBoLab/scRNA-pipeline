"""Build the readable teaching notebook; its execution is verified before delivery."""
from pathlib import Path
import nbformat as nbf

cells = []
def md(text): cells.append(nbf.v4.new_markdown_cell(text))
def code(text): cells.append(nbf.v4.new_code_cell(text))

md("""# 从 PBMC3k 学习 scRNA-seq 全流程

这本 Notebook 使用真实的 10x PBMC3k 未归一化计数，逐步运行与 Snakemake 相同的实现。
从上到下执行，不要先运行后面的格子。每一步请先读它回答的问题，再看输出。

数据：1 位健康供者，2,700 个 PBMC；[10x 来源](https://www.10xgenomics.com/datasets/3-k-pbm-cs-from-a-healthy-donor-1-standard-1-1-0)。
完整代码在 `src/scrna_learn/`，方法解释在 `docs/02_walkthrough.md`。
这份数据用于学习细胞分群，不包含可用于病例对照比较的独立生物学重复。""")
code("""from pathlib import Path
import os
for key in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMBA_NUM_THREADS']:
    os.environ[key] = '1'
root = Path.cwd()
if not (root / 'pyproject.toml').exists() and (root.parent / 'pyproject.toml').exists():
    root = root.parent
assert (root / 'pyproject.toml').exists(), '请在本仓库中打开 Notebook'
os.chdir(root)
import json
import numpy as np
import pandas as pd
import scanpy as sc
from IPython.display import display, Image
from scrna_learn.common import load_config
from scrna_learn.cli import run_stage
cfg = load_config('config/config.yaml')
cfg['output_dir'] = 'results/notebook_pbmc3k'
out = Path(cfg['output_dir'])
print('Scanpy', sc.__version__)
display(pd.DataFrame([cfg['qc']]).drop(columns='sample_overrides'))""")
md("""## 1. 输入：细胞 × 基因的矩阵

`ingest` 会校验演示文件的 SHA256，确认整数 UMI，加入样本和供者信息。
`obs` 是逐细胞 metadata；`var` 是逐基因 metadata。先检查细胞数、基因数和原始计数是否合理。""")
code("""run_stage('ingest', cfg)
counts = sc.read_h5ad(out / 'objects/01_counts.h5ad')
print('cells × genes:', counts.shape)
print('total UMI:', int(counts.X.sum()))
display(counts.obs.head())
display(counts.var.head())""")
md("""## 2. 质控：每一个被移除的细胞都有原因

基础指标由 `sc.pp.calculate_qc_metrics` 计算。当前 PBMC 教学阈值是基因数 ≥200 且 <2500、线粒体比例 <5%。
Scrublet 在每个 capture 的未归一化计数上预测双细胞。请检查自动阈值，而不是默认它一定正确。
先看到全部细胞的去留表，再看总体数量。""")
code("""run_stage('qc', cfg)
qc = pd.read_csv(out / 'tables/cell_qc.csv')
display(qc.groupby('decision', dropna=False).size().rename('cells').to_frame())
display(qc.loc[~qc.retained, ['cell_id', 'n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'decision']].head(10))
display(Image(filename=str(out / 'figures/qc.png')))
for fig in sorted((out / 'figures').glob('doublets_*.png')):
    display(Image(filename=str(fig)))""")
md("""## 3. 标准化与高变基因

核心运算是 `normalize_total(target_sum=1e4)` 和 `log1p`，原始计数保存在 `layers['counts']`。
这里 `highly_variable_genes(flavor='seurat')` 使用 log 表达；不能与要求 counts 的 `seurat_v3` 混淆。
全对象保留全部通过质控的基因，只有构建 PCA 时使用 HVG 子对象。""")
code("""run_stage('normalize', cfg)
a = sc.read_h5ad(out / 'objects/03_normalized.h5ad')
print('保留细胞/基因:', a.shape, 'HVG:', int(a.var.highly_variable.sum()))
print('counts 层是否仍为整数:', np.allclose(a.layers['counts'].data, np.rint(a.layers['counts'].data)))
display(a.var.loc[a.var.highly_variable, ['means', 'dispersions_norm']].head(10))
display(Image(filename=str(out / 'figures/highly_variable_genes.png')))""")
md("""## 4. PCA、邻居图、Leiden 和 UMAP

PCA 提取主要表达差异；邻居图把相似的细胞联系起来；Leiden 在图上分群；UMAP 用于观察。
当前配置使用 40 个 PC、15 个邻居和 1.0 分辨率。
检查线粒体比例是否与某些群重合。UMAP 岛之间的距离不等同于发育时间。""")
code("""run_stage('embed', cfg)
a = sc.read_h5ad(out / 'objects/04_clustered.h5ad')
display(a.obs.cluster.value_counts().sort_index().rename('cells').to_frame())
display(Image(filename=str(out / 'figures/pca_variance.png')))
display(Image(filename=str(out / 'figures/umap_clusters.png')))""")
md("""## 5. 找 marker：不是进行疾病差异推断

对每个 cluster 与其余细胞做 Wilcoxon 探索性比较，BH 校正，查看效应量及表达比例。
输入是全基因的 log-normalized 表达，既不是缩放矩阵，也不是只有 HVG 的矩阵。
这些细胞先被用于聚类，再用于 marker 筛选，所以 p 值不能当作独立验证。""")
code("""run_stage('markers', cfg)
markers = pd.read_csv(out / 'tables/markers_top10.csv')
display(markers.groupby('group').head(3))""")
md("""## 6. 细胞注释：同时看谱系与竞争证据

程序输出候选标签，不将其包装为验证过的分类结果。
例如，NKG7/PRF1 也可以出现在细胞毒 T 细胞中，需要一起看 CD3D/CD3E。
`T/NK unresolved` 提醒你继续检查，不能为了让图好看而强行命名。
点图的大小代表表达细胞比例；颜色是按基因缩放后的平均表达。""")
code("""run_stage('annotate', cfg)
evidence = pd.read_csv(out / 'tables/annotation_evidence.csv')
display(evidence)
display(Image(filename=str(out / 'figures/marker_dotplot.png')))
display(Image(filename=str(out / 'figures/umap_celltypes.png')))""")
md("""## 7. 功能富集

使用固定的 Reactome_2022 人类基因集快照，在本地进行 hypergeometric ORA。
背景是所有受检基因，查询是阳性 marker 前 200 个，校正包含所有符合大小范围的条目。
富集是功能线索，不等同于通路活性或因果机制。""")
code("""run_stage('enrichment', cfg)
display(pd.read_csv(out / 'tables/enrichment_top.csv').head(20))
fig = out / 'figures/enrichment.png'
if fig.exists():
    display(Image(filename=str(fig)))""")
md("""## 8. 保存报告并检查可复现性

报告自带图片，可离线打开。配置、软件版本、输入摘要与阶段输出均被记录。
下面核对最终对象中的 counts 层确实来自原始输入。""")
code("""run_stage('report', cfg)
final = sc.read_h5ad(out / 'objects/06_annotated.h5ad')
expected = counts[final.obs_names, final.var_names].X
assert (final.layers['counts'] != expected).nnz == 0
assert np.isfinite(final.obsm['X_umap']).all()
print('原始 UMI 保存验证通过；UMAP 坐标有效。')
print('报告:', str(out / 'report.html'))
display(pd.Series(json.loads((out / 'run_summary.json').read_text())))""")
md("""## 完成后做三个练习

1. 选择一个 cluster，写出支持和反对其候选标签的 marker 证据。
2. 将分辨率改为 0.5，并使用另一个输出目录，比较合并后的细胞群，而不是比较编号。
3. 解释为什么一个供者的 2,607 个保留细胞不能成为病例对照研究中的 2,607 个独立样本。

进阶入口见 `docs/04_multisample.md`；FASTQ 上游见 `docs/05_fastq.md`；R/Seurat 对应关系见 `docs/07_seurat_mapping.md`。""")

nb = nbf.v4.new_notebook(cells=cells)
nb.metadata.kernelspec = {"display_name": "Python 3 (ipykernel)", "language": "python", "name": "python3"}
nb.metadata.language_info = {"name": "python", "version": "3.12"}
path = Path("notebooks/01_pbmc3k_walkthrough.ipynb")
path.parent.mkdir(parents=True, exist_ok=True)
nbf.write(nb, path)
print(path)
