# 换成自己的数据与样本级统计

## 1. 准备真实样本表

复制 `config/config.yaml` 为 `config/my_project.yaml`，修改：

```yaml
project: my_project
output_dir: results/my_project
input:
  mode: samples
  samples: config/my_samples.tsv
```

样本表必须包含下面六列，每个 capture/library 一行。支持 Cell Ranger 目录、H5 和原始计数 H5AD：

| sample_id | path | format | donor | condition | batch |
|---|---|---|---|---|---|
| S01 | data/raw/S01/filtered_feature_bc_matrix | 10x_mtx | D01 | control | B1 |
| S02 | data/raw/S02/filtered_feature_bc_matrix.h5 | 10x_h5 | D02 | case | B1 |

表中两行只是格式示例，不足以做本仓库的病例对照统计。
`sample_id` 表示一次独立捕获；`donor` 表示独立生物学来源；同一人的技术重复必须用相同 `donor`。
`condition` 是研究分组；`batch` 来自实验记录，不可通过 UMAP 形状反推后随意填入。

```bash
uv run snakemake --configfile config/my_project.yaml --cores 2
```

如果 H5AD 的 `.X` 不是原始计数，需要先在自己的副本中把正确的 counts 层放到 `.X`，并验证其来源。
单靠数值像整数不能证明它是原始 UMI；程序只能做结构与数值检查，不能替你核实数据产生过程。
10x `raw_feature_bc_matrix` 包含空液滴，不能直接当作已确认细胞输入这条 filtered-count 主线。

## 2. 何时做 Harmony

至少两个真实 batch，且需要先审查实验设计。启用方式：

```yaml
embedding:
  integration: harmony
  batch_key: batch
```

保留同一节的其他默认字段。程序把校正结果放在 `X_pca_harmony`，使用该表示建图；`X_pca`、log 表达和原始 counts 保留。
只有一个 batch 时拒绝执行；如果每个 batch 都只对应一种 condition，则拒绝用 Harmony 解决完全混杂的设计。
校正后的 UMAP 混合均匀不等于生物学正确，仍应按 batch、donor、condition 和 cell type 分别观察。
本教学入口不宣称自动选择了最优整合方法，也不提供大规模整合性能基准。

## 3. 有独立供者后做 pseudobulk

程序按 **donor × condition × cell_type** 求和原始 UMI，把技术重复 capture 合并。
每个组合至少 20 个细胞才保留；这也是教学默认值，不是所有研究的统计门槛。
默认每组至少 3 个独立 donor，防止把 3,000 个同一人的细胞当作 3,000 个独立重复；真正的样本量仍需根据研究问题和功效确定。

先人工复核注释并写入配置，再运行：

```bash
uv run scrna-learn pseudobulk \
  --input results/my_project/objects/06_annotated.h5ad \
  --out results/my_project/pseudobulk_T \
  --cell-type "T cells" \
  --case case --control control
```

同一位 donor 在两个条件下都有观测时增加 `--paired`，设计改为 `~ donor + condition`。
未配对入口使用 `~ condition`；禁止在有重复 donor 的情况下误用独立组设计。
每个 pseudobulk 样本的细胞数及未达到最低细胞数的组合都会记录。

PyDESeq2 输入整数样本计数，至少 3 个 pseudobulk 样本中计数 ≥10 的基因进入模型。
输出的正 log2FC 表示 `case` 高于 `control`。程序不会把 Harmony 校正坐标或 log-normalized 值当作 DESeq2 输入。

该入口面向简单平衡教学设计。真实研究还可能需要年龄、性别、配对、批次等协变量以及缺失设计处理。
仅有聚合并不能自动消除混杂；应根据实验设计显式修改模型，而不是机械套用本入口。

依据：[Squair et al., Nature Communications, 2021](https://doi.org/10.1038/s41467-021-25960-2) 强调生物学重复的重要性。
[PyDESeq2 0.5.2 API](https://pydeseq2.readthedocs.io/en/v0.5.2/api/docstrings/pydeseq2.dds.DeseqDataSet.html) 提供本实现使用的公式设计接口。

## 4. GMT 功能富集

`enrichment.gmt` 接受 GMT 文件：每行是通路名、来源字段和若干基因 symbol，以 tab 分隔。
使用真实数据库时保留来源、版本/快照日期、物种和许可，不能把自行列的几个 marker 假装成 GO 或 KEGG。
本仓库提供 hypergeometric ORA；它不是 GSEA，不使用排序全基因统计量。

默认已经提供并启用 `resources/reactome_2022.gmt`，来自 Enrichr 的 Reactome_2022 人类库，包含 1,818 个条目。
这是有日期和 SHA256 的固定教学快照，不是最新数据库声明，详见 `resources/README.md`。
换成非人数据时应替换同物种基因集，或将 `gmt` 设为 `null` 明确跳过。

背景是本次 marker 检验的全部保留基因，查询是每个 cluster 的候选阳性 marker 前 200 个。
基因集先与背景取交集，默认保留 5–500 个实测基因的条目。
对所有合格条目做 BH 校正，包括没有交集、p=1 的条目，避免只校正命中条目造成偏差。
校正范围是每个 cluster 内的条目；跨很多 cluster 挑选显著结果时还应考虑更广泛的多重比较。

富集提示某组 marker 的功能特征，不证明通路活性，也不证明该通路导致某种疾病。
同一批细胞已经用于聚类与 marker 筛选，因此通路结果同样是探索性结果。

本次 PBMC3k 示例中，一些群的前列结果是翻译相关条目。应检查驱动这些结果的核糖体基因和 marker 列表，不能仅凭最小 p 值命名细胞状态。
例如 `Viral mRNA Translation` 这样的数据库条目也含有人类翻译机制的基因；出现该名称不说明供者存在病毒感染。

## 5. 不自动接入的模块

| 分析 | 需要的额外条件 |
|---|---|
| 空液滴与 ambient RNA 校正 | 包含空液滴的 raw matrix、实验背景与适当模型 |
| 细胞周期校正 | 合适的周期基因集，明确周期是否是感兴趣的真实信号 |
| 拟时序 | 连续状态、生物学起点、独立证据；成熟 PBMC 的岛状 UMAP 不能替代 |
| RNA velocity | 正确计算的 spliced/unspliced 等输入及模型适用性 |
| 细胞通讯 | 合适配体受体库、组织语境、样本级重复与验证 |
| 肿瘤 CNV 推断 | 合适参考细胞、基因组位置、组织语境与独立验证 |

这些分支不是每份 scRNA-seq 数据都必须做的“全流程步骤”。本仓库完整运行的是通用计数矩阵分析主线，并对需要额外数据的分支说明适用条件。
