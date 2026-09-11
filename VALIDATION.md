# 实际验证记录

验证日期：2026-09-09。以下数字来自本仓库真实执行产生的输出，参考文件在 `examples/pbmc3k/`。

## 真实数据主线

使用随仓库分发、经过 SHA256 校验的 PBMC3k 原始 UMI 计数，完成 Snakemake 的 10 个任务（含数据准备和最终目标）。
实际流程返回码为 0；随后 dry-run 返回码为 0，并显示 `Nothing to be done`。
在另一个新进程中重复 dry-run，仍确认全部目标存在且最新。工作流状态使用随仓库配置的 SQLite 后端。
本次实测全流程约 124 秒；共享运行环境的耗时不代表其他电脑的性能保证。

| 指标 | 实测值 |
|---|---:|
| 输入细胞 | 2,700 |
| 输入基因 | 32,738 |
| 输入 UMI 总数 | 6,390,631 |
| 基础 QC 后细胞 | 2,638 |
| 移除的预测双细胞 | 31 |
| 最终细胞 | 2,607 |
| 保留基因 | 13,611 |
| 高变基因 | 2,000 |
| Leiden clusters | 9 |
| Reactome cluster × term 检验 | 14,022 |
| 每个 cluster 内 BH 校正后 p < 0.05 的条目总数 | 304 |

候选标签包括 T、B、NK、两类单核细胞、树突细胞、血小板和 `T/NK unresolved`。
9 个 cluster 不等于 9 种已确认的细胞类型；标签尚未经过实验或外部参考验证。
Scrublet 的预测标记也不能视为已知真值。富集属于探索性 marker ORA。

## Notebook 与命令行一致性

教学 Notebook 的 9 个代码单元全部按顺序实际执行，输出已保存在 `.ipynb` 中。
命令行与 Notebook 独立运行得到相同细胞顺序、基因顺序、cluster、候选标签和原始 counts；UMAP 通过 `numpy.allclose` 比较。
最终对象的 `layers['counts']` 与分析前输入按细胞和基因索引后的 UMI 逐项相同；所有 UMAP 坐标有限。
结果清单内的 36 个文件 SHA256 全部匹配。

执行采用进程内 IPython 并捕获真实输出。本运行环境限制 Jupyter 内核的本地 socket，因此没有验证浏览器中 Jupyter 服务与内核通信。
这项限制不影响上述 9 个单元中的分析代码实际执行；在个人电脑打开 Notebook 的命令见 [安装教程](docs/00_setup.md)。

## 软件测试

`pytest -q`：**14 passed**。测试包括：

- 拒绝负数、非有限值和非整数计数。
- 10x MTX、10x H5 的细胞/基因方向和计数读入。
- 完全禁用网络请求时校验、复制和读取随附真实 PBMC3k。
- 在明确标注的合成数据上验证从输入到报告的流程以及 UMI 保存。
- 同一供者技术重复的 pseudobulk 聚合、计数守恒与设计检查。
- 在明确标注的合成重复样本上实际运行 PyDESeq2，检出预设方向的表达变化。
- 在明确标注的合成批次上实际运行 Harmony，并保留原始 counts。
- GMT 背景集合及包括零命中条目在内的 BH 校正。
- Cell Ranger 命令参数构造与字面路径传递。

测试中的合成结果用于验证软件，不作为真实 PBMC 生物学结果。运行中的依赖弃用提醒与 HVG 缩放的稀疏矩阵转 dense 提醒已保留在日志。

## 环境与复现

本次实际环境为 Linux、Python 3.12.14、Scanpy 1.11.5、Snakemake 9.26.1。
完整依赖由 `uv.lock` 固定；实际包版本和源代码 SHA256 见 `examples/pbmc3k/provenance/run_manifest.json`。
随机种子为 0；不同硬件或数值计算实现可能影响浮点结果，不能只凭 UMAP 的方向判定运行错误。

```bash
uv sync --locked --group tutorial
uv run pytest -q
uv run snakemake --cores 2
uv run snakemake --cores 2 --dry-run
```

## 尚未实测的范围

- 没有执行真实 FASTQ → Cell Ranger 比对：需要额外安装软件并准备匹配化学版本的 FASTQ 和参考。
- PBMC3k 只有一位健康供者；未用它声称完成真实多供者整合或疾病差异分析。
- 此处记录的是最初的本地验证；后续 GitHub 云端检查状态以 [Actions](https://github.com/flee70973-coder/F/actions/workflows/ci.yml) 的实际运行记录为准。
- macOS、Windows/WSL2 和大规模数据未做本次性能与兼容性实测。

发布步骤见 [GitHub 指南](docs/06_github.md)。机器可读验证摘要与原始验证日志见 `examples/pbmc3k/validation/`。
