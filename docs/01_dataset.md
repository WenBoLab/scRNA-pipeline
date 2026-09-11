# 为什么用 PBMC3k

## 来源和选择理由

| 项目 | 信息 |
|---|---|
| 名称 | 3k PBMCs from a Healthy Donor |
| 发布者 | 10x Genomics |
| 物种与材料 | 人，外周血单个核细胞 |
| 供者数 | 1 位健康供者 |
| 官方检出细胞数 | 2,700 |
| 本镜像矩阵形状 | 2,700 cells × 32,738 genes |
| 原始分析 | Cell Ranger 1.1.0；Scanpy 镜像由 hg19 计数矩阵转换 |
| 发布时间 | 2016-05-26 |
| 数据许可 | CC BY 4.0 |
| 仓库附带计数大小 | 7,742,219 bytes，约 7.4 MiB；无损重新封装 |
| 最初下载的镜像大小 | 5,855,727 bytes，约 5.6 MiB |
| 原始数据页面 | [10x 官方页面](https://www.10xgenomics.com/datasets/3-k-pbm-cs-from-a-healthy-donor-1-standard-1-1-0) |
| 教学镜像说明 | [scanpy.datasets.pbmc3k](https://scanpy.readthedocs.io/en/stable/generated/scanpy.datasets.pbmc3k.html) |

这份数据体积小，免疫细胞类型便于通过经典 marker 理解，并且可以与 Scanpy、Seurat 的官方入门教程交叉学习。这里的选择依据是教学可理解性和运行成本，数据较旧这一点已明确保留；它不代表最新测序化学版本的性能。

默认输入是随仓库附带的 `resources/pbmc3k_counts.h5ad`，包含未归一化的原始 UMI 计数。
首次数据于 2026-09-08 从 Scanpy 使用的 `pbmc3k_raw.h5ad` 镜像获得；已核对 Scanpy 1.11.5 的 `datasets.pbmc3k` 实现，其镜像由官方 MTX 计数读取后写成 H5AD。
2026-09-09 镜像服务返回 HTTP 502，因此本项目将已校验、未经 QC 的输入对象无损重新封装，恢复原始 barcode，仅保留原始计数和基因 ID，随仓库分发。
逐项比较确认全部 2,700 × 32,738 个计数及基因顺序相同，总 UMI 为 6,390,631，没有做过滤、标准化或插补。

重新封装的 H5AD 字节摘要与镜像文件不同。程序依据 `resources/pbmc3k_counts.source.json` 校验随附文件，并把来源记录写入运行结果。
文件摘要是本项目实测值，不冒充 10x 发布的官方摘要。

随附计数 SHA256：

```text
f8b708c02cb9fd81bfd30bc14c88b5d62549394cedadd487323f8cc7bae6f47e
```

最初镜像文件 SHA256（另存于 `src/scrna_learn/pbmc3k.sha256`）：

```text
89a96f1beaa2dd83a687666d3f19a4513ac27a2a2d12581fcd77afed7ea653a1
```

文件名中的 `raw` 指未归一化的 UMI 计数；这仍然是经过 Cell Ranger cell calling 的 filtered matrix，并不是包括全部空液滴的 raw droplet matrix。因此，这条教学主线无法重新完成空液滴识别或完整估计 ambient RNA。

## 能学什么，不能得出什么

可以学习数据结构、细胞质控、双细胞检测、标准化、高变基因、PCA、聚类、UMAP、marker 和细胞注释。
PBMC 主要用于免疫细胞教学，不能代表其他组织的完整细胞谱系。

不能把一个供者的细胞随机分为“患者”和“对照”，然后将细胞数当作样本数。
也不能为了展示批次校正而把这些细胞人为拆成两个批次，并把它当作真实整合案例。
PBMC3k 没有疾病对照、时间序列或剪接计数，本仓库不会用它给出疾病机制、RNA velocity 或发育轨迹结论。

## 数据授权和引用

代码许可为 MIT；数据和本仓库中基于其生成的图表遵循原始数据的 CC BY 4.0。
使用时标注：10x Genomics, *3k PBMCs from a Healthy Donor*, 2016，附上数据页面，说明图表经过本流程处理。
不要把数据许可证写成代码的 MIT。

参考学习：[Seurat guided clustering tutorial](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html)。本仓库实现与参数独立维护，聚类数量和 UMAP 方向无需与官方截图完全相同。
