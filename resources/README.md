# 教学数据与来源

## PBMC3k 原始计数

`pbmc3k_counts.h5ad` 为 10x PBMC3k 的 2,700 × 32,738 原始 UMI 计数（6,390,631 UMI），从本项目已核验的分析前输入无损重新封装，以便离线学习。没有 QC、标准化或插补。
`pbmc3k_counts.source.json` 记录原始镜像 URL、前后 SHA256、重新封装方式、矩阵形状和数据许可。
完整来源和引用见 [数据说明](../docs/01_dataset.md)。数据由 10x Genomics 发布，采用 [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)，重新封装由本项目完成。

## Reactome 教学基因集

`reactome_2022.gmt` 是 Maayan Lab Enrichr 提供的 `Reactome_2022` 人类基因集，包含 1,818 个条目。
本项目于 2026-09-08 获取并保留原始文件，下载 URL、SHA256 与来源记录在 `reactome_source.json`。

这是一份按库名和文件摘要固定的教学快照，**不宣称是最新 Reactome 版本**。
不将在线服务器的 enrichment 结果复制到本流程；统计由本地实现按实测基因背景重新计算。
原文件第二列为空，本地解析保留该格式，来源以旁边的 JSON 为准。

- [Enrichr](https://maayanlab.cloud/Enrichr/)
- [基因集下载](https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=Reactome_2022)
- [Reactome 数据许可：CC0](https://reactome.org/license)

使用这些数据时请认可 Reactome 和 Enrichr 的工作。
Enrichr 的方法来源之一：Kuleshov et al. *Enrichr: a comprehensive gene set enrichment analysis web server 2016 update*. Nucleic Acids Research. DOI: [10.1093/nar/gkw377](https://doi.org/10.1093/nar/gkw377)。
