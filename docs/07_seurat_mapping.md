# 给熟悉 R / Seurat 的学习者

本仓库的可执行主线采用 Scanpy。下面是概念对应表，不表示两套实现的默认参数或结果完全等价。

| 学习目标 | Seurat 常见操作 | 本仓库 Scanpy 操作 |
|---|---|---|
| 读取 10x | `Read10X` | `sc.read_10x_mtx` / `sc.read_10x_h5` |
| 容器 | Seurat object | AnnData |
| 原始 UMI | RNA assay counts layer | `layers["counts"]` |
| 逐细胞 metadata | `object[[]]` | `adata.obs` |
| QC 指标 | `nFeature_RNA`, `nCount_RNA`, `PercentageFeatureSet` | `calculate_qc_metrics` |
| 归一化 | `NormalizeData` | `normalize_total` + `log1p` |
| 高变基因 | `FindVariableFeatures` | `highly_variable_genes` |
| 缩放 | `ScaleData` | 在独立 HVG 对象上 `scale` |
| PCA | `RunPCA` | `tl.pca` |
| 图与聚类 | `FindNeighbors`, `FindClusters` | `pp.neighbors`, `tl.leiden` |
| UMAP | `RunUMAP` | `tl.umap` |
| Marker | `FindAllMarkers` | `tl.rank_genes_groups` |
| 查看表达 | `FeaturePlot`, `DotPlot` | `pl.umap(color=...)`, `pl.dotplot` |

可以对照 [Seurat 官方 PBMC3k 教程](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html) 学习。
本项目没有把未经运行的 R 脚本列为已验证主线；避免让学习者误以为安装环境和统计结果已在两种语言中逐项验证。

