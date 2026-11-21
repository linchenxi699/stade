# StaDE 使用手册

本手册面向使用者，说明数据准备、命令行与 Python API 的最小可用路径。统计定义见 [method.md](method.md)。

## 1. 输入数据

### 1.1 表达矩阵

制表符分隔（`.tsv` / `.tsv.gz`）或逗号分隔（`.csv`）：

```
gene_id            S1      S2      S3
ENSG00000000003    115     117     168
ENSG00000000005    0       0       2
```

- 第一行为表头，第一列为基因标识；
- 计数矩阵用于过滤与归一化，对数表达矩阵也能直接用（加 `--kind lognorm` 或在 Python 中自行跳过 `prep` 的变换）；
- 不允许重复的基因或样本标识（函数会直接报错，避免静默合并）。

### 1.2 样本元数据

```
file_id      group
S1           Tumor
S2           Normal
```

用 `--meta-index` 指定与表达矩阵列名对齐的那一列，用 `--group-col` 指定分组列，
用 `--case-col Tumor` 指定哪一组是"病例组"。

## 2. 命令行

### 2.1 ``stade prep``

低表达过滤 → TMM → log2 CPM，落盘 `bundle.npz`、`logcpm.tsv.gz`、`filtered_counts.tsv.gz`、`groups.tsv`。

```bash
stade prep --counts counts.tsv.gz --meta meta.tsv --meta-index file_id \
           --group-col group --case-col Tumor --tmm --outdir work
```

| 参数 | 默认 | 说明 |
|---|---|---|
| `--min-count` | 10 | 过滤阈值（计数） |
| `--min-prop` | 0.25 | 至少多少比例的样本达到阈值 |
| `--prior-count` | 0.5 | log2(CPM+prior) 的先验计数 |
| `--tmm` | 关 | 是否做 TMM 归一化 |

### 2.2 ``stade de``

经验贝叶斯差异表达 + bootstrap 稳定性 + 三法一致性 + 经验零分布。

```bash
stade de --bundle work/bundle.npz --case-col Tumor --n-boot 200 \
         --alpha 0.05 --min-abs-logfc 1.0 --min-stability 0.7 --outdir results
```

输出：

| 文件 | 内容 |
|---|---|
| `de_ebayes.tsv` | 原始经验贝叶斯结果 |
| `de_stability.tsv` | 加稳定性分量、`robust` 标记、local FDR |
| `de_consensus.tsv` | 三法一致性 |
| `bootstrap_t.tsv.gz` | 每次 bootstrap 的 t 统计量（供富集稳定性使用） |
| `de_summary.json` | 关键数字汇总 |
| `figs/*.png` | PCA / MA / 火山 / 稳定性 / 一致性 |
| `provenance_de.json` | 参数、输入 SHA-256、环境、摘要 |

### 2.3 ``stade enrich``

```bash
stade enrich --de results/de_stability.tsv \
             --genesets data/raw/genesets/KEGG_2021_Human.txt \
             --library KEGG --robust-only --n-perm 1000 \
             --n-boot 200 --boot-t results/bootstrap_t.tsv.gz \
             --symbol-map data/processed/gene_meta.tsv --symbol-col gene_name --use-symbols \
             --outdir results/enrichment
```

- `--genesets` 支持 `.gmt` 或 Enrichr 文本格式（自动识别）；
- `--robust-only` 只用"显著且稳定"的基因做 ORA；
- `--n-boot` + `--boot-t` 计算基因集选择频率；
- `--symbol-map` 把 Ensembl ID 映射为基因符号后再做富集（推荐）。

### 2.4 ``stade surv``

```bash
stade surv --bundle work/bundle.npz --clinical clinical.tsv \
           --genes candidates.txt --cut median --sig 5 --outdir results/survival
```

`--cut` 可取 `median`、`quartile`、`optimal`（最优切点扫描，仅用于发现阶段）。

## 3. Python API

```python
import pandas as pd
from stade.normalization import log_cpm, filter_low_expression
from stade.stats import fit_ebayes, significance
from stade.stability import bootstrap_de, gene_stability
from stade.enrichment import load_enrichr_tsv, ora_weighted, gsea

counts = pd.read_csv("counts.tsv.gz", sep="\t", index_col=0)
groups = pd.read_csv("meta.tsv", sep="\t", index_col=0)["group"]
expr = log_cpm(filter_low_expression(counts, groups=groups), use_tmm=True)

design = pd.DataFrame({"Intercept": 1.0,
                       "group": (groups == "Tumor").astype(float)}, index=expr.columns)
res = fit_ebayes(expr, design, coef="group")
de = significance(res.table, alpha=0.05, min_abs_logfc=1.0)

boot = bootstrap_de(expr, groups, "Tumor", n_boot=200, seed=20261215)
de = de.join(gene_stability(boot, de["logFC"]))

kegg = load_enrichr_tsv("KEGG_2021_Human.txt")
print(ora_weighted(de["stability"], kegg, universe=set(de.index)).head())
print(gsea(de["t"], kegg, n_perm=1000).head())
```

## 4. 零下载试跑

仓库自带一个**模拟**表达矩阵（基因符号与 KEGG 库为真实，表达值为模拟），用于快速验证流程：

```bash
python scripts/99_make_demo_dataset.py
stade prep --counts data/demo/demo_counts.tsv.gz --meta data/demo/demo_sample_meta.tsv \
           --meta-index sample_id --group-col group --case-col Tumor --tmm --outdir demo_work
stade de --bundle demo_work/bundle.npz --case-col Tumor --n-boot 50 --outdir demo_out
stade enrich --de demo_out/de_stability.tsv --genesets data/raw/genesets/KEGG_2021_Human.txt \
             --library KEGG_demo --robust-only --n-perm 300 --outdir demo_out/enrichment
```

## 5. 常见问题

**Q：没有基因符号，只有 Ensembl ID，能做富集吗？**
A：可以。用 `--symbol-map` 指向含 `gene_id`/`gene_name` 列的注释表；框架会在富集前把 ID 换成符号，并自动去重。

**Q：稳定性分数怎么选阈值？**
A：默认 0.7。先看 `de_stability.tsv` 里稳定性的分布（通常双峰），再结合下游用途；生存/预后筛选建议更严格（≥0.8）。

**Q：bootstrap 要跑多少次？**
A：过夜任务可跑 500–1000 次；笔记本上 200 次已能给出稳定的选择频率估计（本文案例 200 次耗时约 1 分钟）。

**Q：结果能复现吗？**
A：所有随机过程都有 `--seed`；每个子命令输出 `provenance_*.json`，记录输入文件 SHA-256 与全部参数。
