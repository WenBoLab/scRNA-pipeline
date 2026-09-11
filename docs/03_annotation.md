# 用证据给细胞命名

## 先定大类，再讨论亚型

这份人 PBMC marker 面板是教学起点，不能直接移植到小鼠、肝组织、肿瘤上皮细胞或成纤维细胞。
面板结合 [Seurat PBMC3k 教程](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html) 中的经典标记及常用免疫细胞判读思路；它不是从训练数据估计准确率的分类模型。

| 候选大类 | 面板 marker | 判读重点 |
|---|---|---|
| T cells | CD3D、CD3E、TRAC、IL7R | CD3 谱系证据；IL7R 并不覆盖所有 T 亚型 |
| B cells | MS4A1、CD79A、CD79B、CD37 | 多个 B 谱系 marker 共同支持 |
| NK cells | GNLY、NKG7、PRF1、KLRD1 | 需同时检查 CD3，细胞毒 T 也可表达 NKG7、PRF1 |
| CD14 monocytes | CD14、LYZ、S100A8、S100A9 | LYZ 广泛见于髓系，需联合判断 |
| FCGR3A monocytes | FCGR3A、MS4A7、LST1、IFITM3 | FCGR3A 也可见于 NK，结合髓系 marker |
| Dendritic cells | FCER1A、CD1C、CST3、CLEC10A | CST3 单独不足以证明 DC |
| Platelets | PPBP、PF4、GNG11、SDPR | 小 cluster 更应检查低 RNA、混合表达与过滤影响 |

旧版 hg19 矩阵可能没有面板中的部分基因，例如本次 TRAC 未在可用基因名中检出。缺失会记录在证据表，不将“没有这个列”误解为“生物学上不表达”。

## 程序如何生成候选标签

每个面板至少有两个可用 marker 才参加打分。`score_genes` 将该面板的平均 log 表达与表达量相近的对照基因比较，随后计算 cluster 平均分。
最低得分和领先第二名的差值阈值分别在配置中设置。
分数不达标时保留 `Uncertain`；这些阈值是可解释的教学启发规则，不是概率或置信区间。

对 NK 候选另检查 CD3D/CD3E 的平均表达细胞比例；若大于 0.5，则保留 `T/NK unresolved`，提醒检查细胞毒 T 与混合群。
这条辅助规则同样不能替代参考图谱或人工核对，也不是将每个细胞独立分类。

默认将整个 cluster 赋予候选标签。真实细胞类型可能在同一个 cluster 内混合，因此 dotplot、逐细胞表达图和双细胞指标必须一起看。

## 怎样进行复核

打开 `annotation_evidence.csv`、`annotation_expression.csv`、`markers_top10.csv`，并查看 marker 点图与表达 UMAP。
至少回答：主 marker 是否覆盖多数细胞？是否存在互相冲突的谱系？是否由少量异常细胞拉高均值？有没有足够证据进一步区分亚型？

如果证据只支持 T 大类，就写 T cells；不要为了标签丰富而直接写 exhausted、Treg 或特定记忆亚型。
当前主线不进行 T 亚型鉴定。如果需要，可将 T 大类子集重新做 HVG、PCA、图与 marker，并加入相应谱系面板。对子集不应把主线的 UMAP 坐标当作重新聚类。

完成复核后，在自己的 TSV 中填写当前运行的每个 cluster：

```text
cluster	cell_type	evidence
0	T cells	CD3D/CD3E 表达并结合该 cluster 的其他 marker
```

以上只演示格式，实际文件必须完整覆盖当前所有 cluster，不能重复编号或空着证据列。
将路径填写到 `annotation.manual_labels` 后重跑。程序保留原始候选标签，同时增加复核标签。
更改输入、筛选、随机种子或聚类参数后，应重新检查编号对应关系，避免套用旧标签。

