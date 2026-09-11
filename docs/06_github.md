# 发布到 GitHub

本项目使用已有公开仓库 [flee70973-coder/F](https://github.com/flee70973-coder/F)，项目名称为 `scrna-seq-learning`。

## 第一次使用 GitHub

1. 打开仓库首页，向下滚动阅读 README；点击“学习顺序”中的链接进入中文教程。
2. 要下载全部代码，点击绿色 **Code → Download ZIP**，解压后按安装教程运行。
3. 在 **Actions → Offline tests** 查看自动检查；绿色勾表示该次检查通过，红色叉可展开查看失败步骤。
4. 想保留一份到自己的 GitHub 账号，使用右上角 **Fork**；遇到可复现的问题可通过 **Issues** 提交。

后面的命令用于把项目发布到你自己新建的仓库，无需在当前仓库重复执行。

## 使用 GitHub CLI

安装 Git 和 [GitHub CLI](https://cli.github.com/)，在项目根目录运行：

```bash
gh auth login
git init -b main
git add .
git commit -m "Add reproducible scRNA-seq teaching pipeline and PBMC3k example"
gh repo create scrna-seq-learning --public --source=. --remote=origin --push
```

Git 首次提交若提示缺少姓名/邮箱，请按你自己的身份配置。不要使用其他人的 Git 身份。
已经是 Git 仓库时跳过初始化；已有提交时不必重复提交。已有同名远程仓库时不要覆盖，使用明确的仓库地址和分支。
希望先内部检查可把 `--public` 改为 `--private`，准备好后再调整公开状态。

依据：[GitHub 官方：推送本地代码](https://docs.github.com/en/migrations/importing-source-code/using-the-command-line-to-import-source-code/adding-locally-hosted-code-to-github)。

## 仓库上传哪些内容

代码、中文教程、Notebook、配置、锁文件和精简示例结果进入 Git。
`data/raw/`、完整 `results/`、`.venv/`、FASTQ、BAM 和密钥文件被排除。
已许可分发的 PBMC3k 教学计数约 7.4 MiB，随 `resources/` 一起进入 Git；初学者 clone 后无需依赖原镜像站。个人项目的原始数据不上传到 Git。

可将 `examples/pbmc3k/figures/umap_celltypes.png` 放在仓库首页展示，并在 README 介绍学习顺序。
如果分发完整结果，使用 GitHub Release 附件或合适的数据存储，并附上数据许可和校验信息。

## GitHub Actions

| 工作流 | 触发方式 | 内容 |
|---|---|---|
| Offline tests | push / pull request | 软件断言、明确标注的合成数据、Snakemake 端到端运行 |
| Real PBMC3k demo | Actions 页面手动 Run workflow | 校验随附真实 PBMC3k，运行所有默认步骤，上传结果 |

这里“离线”指测试数据不需要联网；首次安装 Python 依赖仍然需要网络。
合成测试结果仅验证软件行为，不能作为真实 PBMC 生物学发现。
工作流只读仓库，不自动推送分析结果或改写主分支。
远程 Actions 的执行状态必须在发布后查看；本地运行通过不能冒充 GitHub 云端通过。

## 适合教学协作的方式

为新分析或参数变更开分支，在 PR 中写清为什么调整、输入是否改变、哪些输出改变，以及如何验证。
学员提问可使用仓库的 issue 模板，附操作系统、命令和错误日志，帮助复现。
不要只上传最终 PNG 而省略参数和数据来源。
