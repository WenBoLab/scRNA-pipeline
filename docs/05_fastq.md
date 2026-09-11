# FASTQ 上游入口

## 两条入口的区别

默认演示从官方 filtered UMI count matrix 起步，已经完成比对、barcode/UMI 处理和 cell calling。
若拿到自己的现代 10x Chromium 单独 Gene Expression 文库 FASTQ，可使用本节入口连接上游。
仅有本仓库文件时不会声称重新运行过 Cell Ranger。

PBMC3k 2016 数据使用旧 GemCode/Cell Ranger 1.1.0 化学版本，read 布局与现代 10x 3′ 文库不同。
下面的现代 Cell Ranger 命令**不是**重处理这份历史数据 FASTQ 的兼容性承诺。
如需完整重现该历史上游，必须准备与化学版本、参考和软件匹配的环境。

## 准备工作

在满足 [10x 系统要求](https://www.10xgenomics.com/support/software/cell-ranger/latest/resources/cr-system-requirements) 的 Linux 机器上安装 Cell Ranger，
并准备匹配物种的 Cell Ranger transcriptome reference。
软件和参考不包含在 Python lock 中，也不打包入本仓库。

确认 FASTQ 命名、sample prefix、物种与参考版本。BCL 到 FASTQ 的 demultiplexing 使用测序平台支持的工具完成，本入口假定已有 FASTQ。
Flex、Multiome 和某些 multiplex assays 需要 `multi` 或 `arc` 等不同入口，不应套用这里的单 GEX `count`。

## 预览命令

```bash
uv run python scripts/cellranger_count.py \
  --sample-id sample_A \
  --fastqs /absolute/path/to/fastqs \
  --fastq-sample sample_A \
  --transcriptome /absolute/path/to/refdata-gex-GRCh38 \
  --donor donor_A --condition control --batch batch1
```

默认只打印即将运行的命令。实际运行时增加 `--execute`；增加 `--run-analysis` 可在计数结束后继续到报告：

```bash
uv run python scripts/cellranger_count.py \
  --sample-id sample_A \
  --fastqs /absolute/path/to/fastqs \
  --fastq-sample sample_A \
  --transcriptome /absolute/path/to/refdata-gex-GRCh38 \
  --donor donor_A --condition control --batch batch1 \
  --cores 8 --memory-gb 64 --execute --run-analysis
```

示例中的 8 cores/64 GB 来自常见官方命令形式，需按当前 Cell Ranger 要求、参考和规模调整。
默认显式设置 `--create-bam=false`；需要 BAM 时使用 `--create-bam`。
Cell Ranger 命令使用参数数组执行，不通过 shell 拼接样本路径。

成功后创建计数输出、单样本 TSV 和下游 YAML，路径打印在终端。
先阅读 Cell Ranger 的 `web_summary.html`：确认比对和细胞识别质量，不能只看最终 UMAP。
随后按组织类型调整 YAML 的 QC 与 marker，再分析自己的数据。

命令依据：[Cell Ranger count 官方说明](https://www.10xgenomics.com/support/software/cell-ranger/10.0/analysis/running-pipelines/cr-gex-count)。
本仓库验证命令构造和下游接口；由于没有安装 Cell Ranger、现代 FASTQ 与参考，本次没有执行实际 FASTQ 比对。

