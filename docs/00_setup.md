# 安装、运行与排错

## 1. 先准备什么

建议使用 Linux，或 Windows 上的 WSL2 Ubuntu。macOS 可运行矩阵分析；Cell Ranger 上游按 10x 的 Linux 系统要求准备。
本仓库固定 Python 3.12，所有依赖记录在 `uv.lock`。建议为 PBMC3k 矩阵教学预留 8 GB 内存和 2 GB 可用磁盘；这是教学环境建议，不是最低硬件保证。首次安装要联网。

安装 [uv](https://docs.astral.sh/uv/getting-started/installation/) 后，打开终端并进入包含 `pyproject.toml` 的目录：

```bash
cd scrna-seq-learning
uv sync --locked --group tutorial
uv run scrna-learn --help
```

`uv sync` 会在项目内创建独立 `.venv`。不需要把依赖装进系统 Python，也不需要先安装 R。
`--group tutorial` 增加 Jupyter 和 Notebook 执行工具；只运行命令行流程时可用 `uv sync --locked`。

## 2. 第一次运行

```bash
uv run snakemake --cores 2 --dry-run
uv run snakemake --cores 2
```

第一条只列任务；第二条会校验并复制仓库附带的原始计数、分析并生成报告。`--cores 2` 限制并发任务数，单个分析进程内部使用一个数值计算线程，避免重复争抢 CPU。资源声明 `mem_mb` 用于调度，不能代替操作系统的硬内存限制。

运行完成后，用浏览器打开 `results/pbmc3k/report.html`。报告图片内嵌，不需要联网。
再次执行同一命令，Snakemake 应提示没有工作需要执行。修改配置或代码后会重新计算。
目前配置作为整体追踪，调整一个参数可能保守地重跑多个阶段，这是为了保证结果不会意外沿用旧配置。

`profiles/default/config.yaml` 会被 Snakemake 自动加载，使用 SQLite 保存任务状态与锁。
这是工作流的一部分，请保留该目录；本版锁定的 Snakemake 支持此后端。
如果在没有任务运行时因异常中断留下锁，可先运行 `uv run snakemake --unlock`，再重新执行。

## 3. Notebook 学习方式

```bash
uv run --group tutorial jupyter lab notebooks/01_pbmc3k_walkthrough.ipynb
```

选择项目 Python 内核，从第一格向下执行。Notebook 使用 `results/notebook_pbmc3k/`，每个阶段都调用同一份分析实现，避免教程和自动化流程逐渐变成两个不同版本。每一步都会展示关键对象和检查结果。

如果已有 Jupyter，无法找到项目内核，可在本机执行：

```bash
uv run --group tutorial python -m ipykernel install --user --name scrna-learning --display-name "scRNA Learning"
```

## 4. 文件去哪里找

| 文件 | 你应检查的内容 |
|---|---|
| `report.html` | 总览、图表、分析边界 |
| `tables/cell_qc.csv` | 每个细胞的指标、过滤原因、双细胞状态 |
| `tables/qc_summary.json` | 每个样本的过滤数量和阈值 |
| `tables/gene_metrics.csv` | 高变基因标记、表达均值等 |
| `tables/markers_all.csv` | 所有通过基因过滤的基因的 cluster 检验 |
| `tables/markers_top10.csv` | 每个 cluster 的候选阳性 marker |
| `tables/annotation_evidence.csv` | 候选标签、竞争标签、marker 得分和缺失情况 |
| `tables/annotation_expression.csv` | marker 在 cluster 内的平均表达与表达比例 |
| `tables/celltype_counts.csv` | 每个样本中的细胞类型数量和比例 |
| `objects/06_annotated.h5ad` | 最终对象，可继续分析 |
| `logs/` | Snakemake 每一步的独立日志 |
| `provenance/` | 配置快照、版本、输入与结果 SHA256、耗时 |

## 5. 常见错误

| 情况 | 处理方法 |
|---|---|
| 缺少演示数据 | 确认完整解压或 clone 了仓库，`resources/pbmc3k_counts.h5ad` 和旁边的来源 JSON 均存在；默认流程无需访问镜像下载站 |
| Checksum mismatch | 当前文件与已验证的演示文件不同。删除本地这一个损坏的下载文件后重试；不要关掉校验 |
| No mitochondrial genes | 确認基因名是 symbol，人与小鼠的前缀分别常见为 `MT-`、`mt-`；Ensembl ID 需可靠映射 |
| Input is not integer UMI counts | 传入了归一化表达或缩放对象。选择 Cell Ranger 计数或含真实原始 UMI 的文件 |
| 太多细胞被过滤 | 查看 QC 分布、组织类型和失败原因，再调整阈值；不要为得到预期细胞数盲目放宽 |
| Scrublet 无阈值或报错 | 检查是否单个 capture 的原始计数、细胞是否太少，查看得分分布后再显式设置阈值 |
| 内存不足 | 先降低 `n_hvg`，避免把所有基因转为 dense；多样本项目改用更大机器 |
| Snakemake 显示 incomplete | 查看失败阶段日志，修复原因后用 `--rerun-incomplete` 重新执行 |
| 手动标签不匹配 cluster | 重新检查当前分辨率的 cluster 列表，完整填写一次每个编号；编号不具有跨运行生物学含义 |

不要把所有 warnings 都当成失败。HVG 独立对象在 `scale` 时转为 dense 是本实现的有意选择；原始计数与全基因 log 表达保留为 sparse。程序异常退出、缺文件、非有限值则应停止解读结果。
