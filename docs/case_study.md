# 案例研究：TCGA-LIHC 肝细胞癌转录组

本页记录论文所用案例的设计、数据来源与可复现命令。所有数字均由本仓库代码在本机实际运行得到。

## 1. 研究设计

| 项目 | 内容 |
|---|---|
| 发现队列 | TCGA-LIHC（GDC STAR-Counts），371 例原发肿瘤 + 50 例癌旁 |
| 验证队列 | GSE14520（GPL3921），225 例肿瘤 + 220 例癌旁 |
| 分组 | Tumor vs Normal |
| 主要终点 | 差异表达（FDR<0.05 且 \|log2FC\|≥1）、bootstrap 稳定性、通路富集 |
| 次要终点 | 总生存（OS）、跨队列差异方向与预后价值复现 |

## 2. 数据获取

```bash
python scripts/01_download_data.py all    # 约 1.8 GB，断点续传
python scripts/02_build_matrix.py all     # counts.tsv.gz + sample_meta.tsv + clinical.tsv
```

- GDC API：`/files` 过滤 `TCGA-LIHC` + `Gene Expression Quantification` + `STAR - Counts`，逐文件下载后按 `unstranded` 列拼接；
- 临床：GDC `/cases` 取 `days_to_death` / `days_to_last_follow_up`，与 UCSC Xena `LIHC_survival.txt` 合并；
- 基因集：Enrichr `geneSetLibrary?mode=text`（KEGG 2021、GO BP 2023、Hallmark 2020、Reactome 2022）；
- 验证队列：GEO FTP 的 GPL3921 series matrix、`GSE14520_Extra_Supplement.txt.gz` 与 `GPL3921.annot.gz`。

## 3. 分析流程

```bash
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
python scripts/04c_validation_summary.py
```

## 4. 关键结果

### 4.1 差异表达与稳定性

| 指标 | 数值 |
|---|---|
| 过滤后基因 | 20,253 |
| FDR<0.05 且 \|log2FC\|≥1 | 3,937（上调 3,152 / 下调 785） |
| 稳定性 ≥0.7 | 3,510（89.2%） |
| 显著基因稳定性中位数 | 0.976 |
| 先验自由度 d0 / 后验总自由度 | 3.42 / 422.4 |
| 经验零分布 SD / π0 | 4.12 / 0.922 |

### 4.2 富集（综合证据排序）

- **KEGG**：下调视黄醇代谢（12.4×, padj=3.8e-19）、药物代谢、外源物代谢、补体与凝血级联、胆汁分泌；上调细胞周期。
- **GO BP**：上调 DNA 代谢、DNA 损伤应答、DNA 修复、细胞周期正调控、姐妹染色单体分离、有丝分裂纺锤体；下调类固醇代谢、环氧化酶 P450 通路。
- **Hallmark**：上调 E2F 靶基因（4.5×, padj=4.7e-55）、G2/M 检查点、Myc 靶基因；下调外源物代谢、脂肪酸代谢、凝血、胆汁酸代谢。

### 4.3 生存

365 例肿瘤样本有可用 OS（130 事件），中位随访约 19.6 个月；3,824 个基因 FDR<0.05。
最强单基因（每 SD）：UTP11 HR=1.75、YBX1 1.72、KPNA2 1.72、G6PD 1.64。

### 4.4 外部验证（GSE14520）

| 指标 | 数值 |
|---|---|
| 可比基因 | 10,749 |
| log2FC Spearman ρ | 0.691 |
| 稳健基因方向一致率 | 87.4% |
| 稳健基因在验证队列显著率 | 23.2%（全体基因 7.2%） |
| 预后基因可检验 / 复现 | 2,405 / 96（FDR<0.05），方向一致率 73.0% |

**稳定性 ↔ 复现率**（TCGA 稳定性分层 → GSE14520 复现显著率）：

| 稳定性区间 | 可比基因 | 方向一致率 | 复现显著率 |
|---|---|---|---|
| 0–0.2 | 7,727 | 64.7% | 2.2% |
| 0.2–0.5 | 435 | 79.8% | 10.1% |
| 0.5–0.8 | 404 | 84.7% | 8.9% |
| 0.8–1.0 | 2,183 | 87.8% | **24.1%** |

## 5. 复现说明

- 全流程固定随机种子（默认 20261215）；
- 每步输出 `provenance_*.json`，含输入 SHA-256、参数、依赖版本与关键摘要；
- 论文中的每个数字都可在 `results/` 对应表格中检索到；
- 结果表体积较大，仓库不提交，按上述命令可完整重建。
