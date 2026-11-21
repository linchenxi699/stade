# StaDE · 稳定性感知的转录组差异表达与功能富集分析框架

[![CI](https://github.com/your-account/stade/actions/workflows/ci.yml/badge.svg)](https://github.com/your-account/stade/actions/workflows/ci.yml)
[![version](https://img.shields.io/badge/version-0.1.0-blue.svg)](#)
[![Python](https://img.shields.io/badge/python-%E2%89%A53.10-blue.svg)](#)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![Tests](https://img.shields.io/badge/tests-18%20passed-brightgreen.svg)](tests)

> **StaDE** = **Sta**bility-aware **D**ifferential **E**xpression.
> 一套可复现、可解释、以"稳定性"为一等公民的公共转录组分析框架。

公共转录组数据的二次分析普遍存在一个被忽视的问题：**显著性不等于可重复性**。同一个数据集，换一种检验方法、换一次随机重抽样、换一个预处理参数，差异基因列表与富集通路就可能明显漂移。StaDE 把"稳定性"从一句讨论空话变成可以计算、可以作图、可以直接进入富集统计量的量。

---

## 1. 它解决什么问题

传统流程给出的是"哪些基因显著"。StaDE 额外回答三个问题：

1. **这个基因的显著性有多稳？** 分层 bootstrap 重抽样 200 次，统计每个基因被反复选中的频率、方向一致性与效应量一致性，合成 0–1 的稳定性分数。
2. **这条通路的富集有多稳？** 把基因级稳定性传入超几何检验与 GSEA，并统计通路在 bootstrap 排序中的入选频率，融合成"显著且可重复"的综合证据。
3. **换一个队列还能复现吗？** 内置 GEO 队列解析与跨队列验证，输出复现率随稳定性分层的单调曲线。

## 2. 核心特性

| 模块 | 能力 | 关键实现 |
|---|---|---|
| 差异表达 | limma 风格经验贝叶斯 moderated t，均值—方差趋势先验，稳健方差收缩 | `stade.stats.fit_ebayes` |
| 方法一致性 | 同时给出 Welch t、Wilcoxon 秩和检验与三法一致性评分 | `stade.stats` / `stade.stability.consensus_table` |
| **基因稳定性** | 分层 bootstrap，输出选择频率、符号一致性、效应量一致性与综合稳定性分数 | `stade.stability.bootstrap_de` / `gene_stability` |
| 零分布校准 | 置换零分布 + Efron 经验零分布 local FDR | `stade.stability.empirical_null` |
| 经典富集 | 超几何 ORA、加权 Kolmogorov-Smirnov GSEA（基因集置换） | `stade.enrichment.ora` / `gsea` |
| **稳定性加权富集** | 稳定性加权 ORA（泊松—二项）、稳定性加权 GSEA、bootstrap 基因集选择频率、综合证据分数 | `stade.enrichment.ora_weighted` / `bootstrap_enrichment` / `integrated_enrichment` |
| 生存分析 | 单/多因素 Cox、Kaplan-Meier + log-rank、DE×预后联合排序 | `stade.survival` |
| 跨队列验证 | GEO series matrix 解析、logFC 一致性、预后复现 | `stade.geo` + `scripts/04*` |
| 可复现性 | 每步自动落盘 provenance JSON（输入 SHA-256、参数、环境、时间戳） | `stade.io.Provenance` |
| 出图 | PCA、MA、火山、稳定性、热图、富集气泡、KM、森林图、验证散点 | `stade.plotting` |

## 3. 安装

```bash
git clone https://github.com/your-account/stade.git
cd stade
python -m venv .venv
# Windows: .venv\Scripts\activate    macOS/Linux: source .venv/bin/activate
pip install -e ".[dev]"
python -m pytest -q        # 18 passed
```

仅依赖 `numpy / pandas / scipy / matplotlib / statsmodels / lifelines`，不需要 R、不需要 Bioconductor。

## 4. 零下载试跑（1 分钟）

仓库自带一个**模拟**表达矩阵（基因符号与 KEGG 基因集为真实，表达值为负二项模拟），用于快速验证安装是否成功：

```bash
python scripts/99_make_demo_dataset.py

stade prep   --counts data/demo/demo_counts.tsv.gz --meta data/demo/demo_sample_meta.tsv \
             --meta-index sample_id --group-col group --case-col Tumor --tmm --outdir demo_work
stade de     --bundle demo_work/bundle.npz --case-col Tumor --n-boot 50 --outdir demo_out
stade enrich --de demo_out/de_stability.tsv --genesets data/raw/genesets/KEGG_2021_Human.txt \
             --library KEGG_demo --robust-only --n-perm 300 --outdir demo_out/enrichment
```

也可以在 Python 里跑一个完整小例子：`python examples/demo_simulated.py`。

## 5. 真实数据案例：TCGA-LIHC

```bash
python scripts/01_download_data.py all         # GDC STAR counts + Enrichr 基因集 + Xena 生存表 + GSE14520
python scripts/02_build_matrix.py all          # 拼接矩阵 + 临床随访

stade prep --counts data/processed/counts.tsv.gz --meta data/processed/sample_meta.tsv \
           --meta-index file_id --group-col group --case-col Tumor --tmm --outdir work

stade de --bundle work/bundle.npz --case-col Tumor --n-boot 200 --outdir results

for LIB in KEGG_2021_Human GO_Biological_Process_2023 MSigDB_Hallmark_2020; do
  stade enrich --de results/de_stability.tsv --genesets data/raw/genesets/$LIB.txt \
               --library $LIB --robust-only --n-perm 1000 \
               --symbol-map data/processed/gene_meta.tsv --symbol-col gene_name --use-symbols \
               --outdir results/enrichment
done

python scripts/03b_finish_survival.py
python scripts/04_validate_geo.py
python scripts/04b_validate_prognosis.py
```

### 5.1 结果摘要（用本项目代码实际跑出）

| 步骤 | 结果 |
|---|---|
| 有效基因 / 样本 | 20,253 / 421（371 肿瘤 + 50 癌旁） |
| 差异基因（FDR<0.05, \|log2FC\|≥1） | 3,937（上调 3,152 / 下调 785） |
| 稳健差异基因（稳定性 ≥0.7） | 3,510（89.2%） |
| 经验零分布 | π0 ≈ 0.922，零分布 SD = 4.12 |
| KEGG 代表通路 | 下调：视黄醇代谢（12.4×, padj=3.8e-19）、药物代谢、外源物代谢、补体凝血、胆汁分泌；上调：细胞周期 |
| Hallmark 代表 | 上调 E2F 靶基因（4.5×, padj=4.7e-55）、G2/M、Myc；下调外源物代谢、脂肪酸代谢、凝血、胆汁酸代谢 |
| 生存（365 样本 / 130 事件） | 3,824 个基因 FDR<0.05；最强 UTP11 HR=1.75、YBX1 1.72、KPNA2 1.72、G6PD 1.64（每 SD） |
| 外部验证（GSE14520） | 10,749 个基因 log2FC ρ=0.691；稳健基因方向一致率 87.4%、验证显著率 23.2%（背景 7.2%） |
| **稳定性分层的复现率** | 2.2% → 10.1% → 8.9% → **24.1%**（稳定性 0–0.2 / 0.2–0.5 / 0.5–0.8 / 0.8–1.0） |

最后一行是全文最硬的证据：复现率随稳定性显著上升（最高区间是最低区间的十倍以上），说明稳定性分数捕捉的是跨平台可移植的真实信号，而非对 bootstrap 噪声的过拟合。

## 6. 作为 Python 库

```python
import pandas as pd
from stade.normalization import log_cpm, filter_low_expression
from stade.stats import fit_ebayes, significance
from stade.stability import bootstrap_de, gene_stability
from stade.enrichment import load_enrichr_tsv, ora_weighted

counts = pd.read_csv("counts.tsv.gz", sep="\t", index_col=0)
groups = pd.read_csv("meta.tsv", sep="\t", index_col=0)["group"]
expr = log_cpm(filter_low_expression(counts, groups=groups), use_tmm=True)

design = pd.DataFrame({"Intercept": 1.0,
                       "group": (groups == "Tumor").astype(float)}, index=expr.columns)
res = fit_ebayes(expr, design, coef="group")
de  = significance(res.table, alpha=0.05, min_abs_logfc=1.0)

boot = bootstrap_de(expr, groups, "Tumor", n_boot=200, seed=20261215)
de   = de.join(gene_stability(boot, de["logFC"]))

kegg = load_enrichr_tsv("KEGG_2021_Human.txt")
print(ora_weighted(de["stability"], kegg, universe=set(de.index)).head())
```

## 7. 稳定性是怎么定义的

对基因 *g*，分层 bootstrap 重抽样 *B* 次（保持两组样本量不变）：

- 选择频率 `sel_freq = mean_b[ padj_b < α 且同向且 |logFC_b| ≥ τ ]`
- 符号一致性 `sign_consistency = mean_b[ sign(logFC_b) == sign(logFC) ]`
- 效应量一致性 `cv_score = clip(1 - sd_b(logFC_b) / |mean_b(logFC_b)|, 0, 1)`
- 综合分数 `stability = sel_freq × (0.7·sign_consistency + 0.3·cv_score)`

乘积形式以"是否被重复选中"作为门控：从未被选中的基因稳定性为 0，避免加权求和把分数压缩到 0.3–0.6 的窄区间。

## 8. 目录结构

```
stade/                核心库
  io.py               数据读取、SHA-256 溯源
  normalization.py    过滤、TMM、log-CPM、分位数归一化
  stats.py            经验贝叶斯 moderated t、Welch、Wilcoxon、BH
  stability.py        bootstrap、稳定性分量、经验零分布、三法一致性
  enrichment.py       ORA、稳定性加权 ORA、GSEA、bootstrap 富集、综合证据
  survival.py         Cox、KM、候选基因联合排序
  geo.py              GEO series matrix 解析
  plotting.py         全部出图函数
  cli.py              命令行入口
scripts/              数据下载、矩阵构建、生存、验证、论文生成、demo 数据
tests/                18 个单元/性质测试
docs/                 method.md（统计定义）、usage.md（使用手册）
examples/             demo_simulated.py（合成数据完整流程）
data/demo/            自带模拟数据集（可直接跑）
data/raw/genesets/    随仓库提供的 KEGG / GO / Hallmark 基因集
```

## 9. 可复现性

- 所有随机过程都有显式 `--seed`，默认固定；
- 每个子命令写出 `provenance_*.json`：输入 SHA-256、全部参数、Python 与依赖版本、开始/结束时间、关键结果摘要；
- 下载脚本断点续传并按文件大小校验，失败自动重试；
- 论文里的每一个数字都可以在 `results/` 对应文件中追溯。

## 10. 已知局限

1. bootstrap 只反映样本抽样变异，不含批次与平台变异；
2. 稳定性加权 ORA 使用正态近似，基因集很小或权重极端时精度下降；
3. GSEA 采用基因集置换，强共调控基因集可能偏乐观；
4. 生存模型默认不校正临床混杂，HR 应视为发现阶段效应量；
5. 稳定性阈值 0.7 为经验值，应结合具体数据的分布选择。

详见 [docs/method.md](docs/method.md) 第 8 节。

## 11. 引用

```bibtex
@software{stade2026,
  title  = {StaDE: stability-aware differential expression and functional enrichment for public transcriptomes},
  author = {Lin, Chenxi},
  year   = {2026},
  note   = {Course project, Fuzhou University},
  url    = {https://github.com/your-account/stade}
}
```

## 12. 许可证

MIT License，见 [LICENSE](LICENSE)。欢迎提 issue 与 PR，见 [CONTRIBUTING.md](CONTRIBUTING.md)。
