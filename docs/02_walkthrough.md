# 按问题学习每个分析步骤

## 1. 输入矩阵是什么

AnnData 的行是细胞，列是基因。`obs` 保存逐细胞信息，`var` 保存逐基因信息，`obsm` 保存 PCA、UMAP 等细胞坐标。
10x 的 MTX 在文件层面通常是基因 × 细胞，`read_10x_mtx` 会读成 Scanpy 所需方向。

```python
import scanpy as sc
a = sc.read_h5ad("results/pbmc3k/objects/01_counts.h5ad")
print(a.shape)
print(a.obs.head())
print(a.var.head())
```

先确认矩阵是非负整数计数，而不是已经取对数或缩放过的表达值。
每个 capture 使用唯一 `sample_id`，细胞名加样本后缀，避免不同样本的相同 barcode 被混淆。
多样本默认取共同测量的基因交集，基因数量损失会记录在输入审计中；不同参考注释应先统一。

## 2. 哪些细胞需要审查

主要指标是检测到的基因数、UMI 总数和线粒体 UMI 比例。PBMC 教学配置为：

| 参数 | 默认值 | 具体边界 |
|---|---:|---|
| `min_genes` | 200 | 基因数 ≥ 200 |
| `max_genes` | 2500 | 基因数 < 2500 |
| `max_pct_mt` | 5 | 线粒体计数比例 < 5% |
| `min_cells_per_gene` | 3 | 至少在 3 个保留细胞中检测到 |

这些是 PBMC3k 的教学起点。不同组织、损伤状态、细胞核数据不能照抄阈值。
线粒体比例偏高可能提示应激或低质量，但它不是单独判定“死细胞”的实验金标准。

`cell_qc.csv` 保存原始细胞的所有判断。一个细胞可能同时触发多个失败原因，因此不能把各失败类别简单相加。
多样本逐 capture 过滤，支持 `qc.sample_overrides`。

## 3. 为什么检测双细胞

一个液滴可能含两个细胞，混合表达会干扰聚类与注释。
本实现对基础质控后、仍未归一化的每个 capture 分别调用 `sc.pp.scrublet`。
输出实测分数、模拟双细胞分数和实际阈值。

默认 `expected_rate: 0.04` 是演示设定，不是这份数据的已知真实双细胞率。
自动阈值也可能不理想，应结合上样信息、得分分布和 marker 混合模式复核。
本算法更容易发现不同类型构成的异型双细胞，不能保证移除所有同型双细胞。
关闭检测时报告会明确显示没有运行，不能将 False 标记解释为已证明是单细胞。

## 4. 标准化与原始计数怎样共存

对每个细胞按总 UMI 缩放到 10,000，再计算 `log(1+x)`，减少测序深度造成的尺度差别。
这不意味着所有技术因素都已消除。

```python
a.layers["counts"] = a.X.copy()
sc.pp.normalize_total(a, target_sum=1e4)
sc.pp.log1p(a)
```

| 对象位置 | 内容 | 用途 |
|---|---|---|
| `a.layers["counts"]` | 原始整数 UMI，经过细胞/基因筛选 | pseudobulk、重新建模 |
| `a.X` | 全部保留基因的 log-normalized 表达 | marker、表达图、注释 |
| `a.var["highly_variable"]` | 哪些基因是 HVG | 构建低维空间 |
| `a.obsm["X_pca"]` | PCA 坐标 | 邻居图；保留未整合表示 |
| `a.obsm["X_pca_harmony"]` | 可选的整合 PCA 坐标 | 多样本图结构 |
| `a.obsm["X_umap"]` | UMAP 坐标 | 可视化 |

## 5. 为什么选择高变基因

并非每个基因都对细胞间差异同等有用。选取 2,000 个高变基因帮助构建较紧凑的表示。
本实现使用 `flavor="seurat"`，其输入是 log-normalized 表达。
`seurat_v3` 是不同方法，要求原始计数；不能只更改参数名称而继续沿用同一输入。
多样本在 `sample_id` 层面选择 HVG，减少被单一样本特有高变基因主导的机会。

全对象不会裁剪为只有 HVG。我们复制 HVG 子对象做缩放和 PCA，然后把 PCA 坐标写回全基因对象。
这样能避免在寻找 marker 时遗漏非 HVG 的重要基因。

## 6. PCA、邻居图与聚类各做什么

PCA 将表达差异压缩到较少的维度。本教学配置最多使用 40 个 PC，具体数目也受细胞数和 HVG 数限制。
检查解释方差图；不要把固定“40”当成所有数据都适合的答案。

邻居图基于 PCA 中的相近细胞构建。Leiden 在这个图上寻找群落。
分辨率控制聚类粒度，但更高分辨率不等于更准确。
本版选择 `resolution: 1.0` 作为 PBMC3k 教学起点，随后需要用 marker 判断分出来的差别。

UMAP 和 t-SNE 将细胞放到二维平面。聚类来自 PCA 邻居图，不是在 UMAP 图片上手动画圈。
图上的整体方向、不同岛之间的远近不能直接转换为分化顺序或谱系关系。

## 7. Marker 的统计含义

`rank_genes_groups` 对每个 cluster 与其余细胞进行 Wilcoxon 比较，并在该 cluster 的全部受检基因内做 BH 校正。
输入为未缩放的 log-normalized 全基因矩阵。`logfoldchanges` 是根据平均 log 表达换算的近似 log2FC。
默认候选阳性 marker 要求近似 log2FC ≥ 0.25 且校正 p < 0.05，随后按得分输出前十个。

这些细胞先参与了聚类，又被用于检验 cluster 差异，因此结果用于探索和描述。
细胞很多时很小的表达差别也会有很小的 p 值，应同时看效应量、表达比例、多个 marker 和已知生物学。
疾病或处理效应改用有真实独立重复的样本级设计，见进阶章节。

## 8. 建议完成的练习

1. 查看 20 个被移除的细胞，解释各自去留原因。
2. 将分辨率改为 0.5，使用新的 `output_dir` 重跑，比较哪些群合并了；不要只比较 cluster 编号。
3. 比较 CD3D、NKG7、GNLY 的表达，解释为什么 NKG7 阳性不能单独证明是 NK 细胞。
4. 从全部 marker 中找一个非 HVG，说明为什么差异分析没有只使用 HVG。
5. 更改线粒体阈值，在新输出目录重跑，检查保留数量及细胞组成如何变化。
6. 在全基因原始 counts 层核对任意一个细胞，证明标准化未覆盖原始计数。

相关 API：[QC 与聚类教程](https://scanpy.readthedocs.io/en/stable/tutorials/basics/clustering.html)、[HVG 输入约定](https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.highly_variable_genes.html)、[Scrublet](https://scanpy.readthedocs.io/en/1.11.x/api/generated/scanpy.pp.scrublet.html)、[marker 检验](https://scanpy.readthedocs.io/en/1.11.x/generated/scanpy.tl.rank_genes_groups.html)。

