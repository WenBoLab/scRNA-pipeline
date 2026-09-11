# PBMC3k 实际运行示例

这些图表和统计来自本项目对真实 10x PBMC3k 的运行，参数见 `provenance/config_used.yaml`，验证范围见仓库根目录的 `VALIDATION.md`。

先在本机浏览器打开 `report.html`；图片已内嵌。GitHub 文件预览通常显示 HTML 源码，下载后再打开即可。
`figures/` 适合直接查看，`tables/` 保留核心证据表，`validation/` 记录测试与运行结果。
这里是精简示例；完整中间 H5AD 和全部 marker 检验表由运行流程重新生成，另附完整结果压缩包。

原始数据由 10x Genomics 发布：*3k PBMCs from a Healthy Donor*（2016），[数据页面](https://www.10xgenomics.com/datasets/3-k-pbm-cs-from-a-healthy-donor-1-standard-1-1-0)，采用 [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)。
本项目对数据进行了质控、标准化、聚类、候选注释与可视化；图表和衍生结果保留数据来源授权要求。代码许可另见根目录 MIT LICENSE。
Reactome_2022 来源与许可见 `resources/README.md`。
