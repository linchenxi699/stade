# StaDE 方法说明

本文档说明 StaDE 中每个统计量的定义、实现选择与已知局限，便于审稿人与使用者复核。

## 1. 数据约定

- 表达矩阵：行 = 基因（Ensembl ID 或基因符号），列 = 样本。
- 分组向量：与矩阵列同序的 `pandas.Series`。
- 计数矩阵用于过滤与归一化；对数表达矩阵用于线性模型。

## 2. 预处理

| 步骤 | 规则 | 可调参数 |
|---|---|---|
| 低表达过滤 | 在至少一个分组中，≥ `min_count` 计数的样本比例 ≥ `min_prop` | `min_count=10`, `min_prop=0.25` |
| 归一化 | TMM（trimmed mean of M-values），log-ratio 截尾 30%，绝对表达截尾 5%，因子几何均值为 1 | `logratio_trim=0.3`, `sum_trim=0.05` |
| 变换 | `log2(CPM + prior)` | `prior=0.5` |

## 3. 经验贝叶斯 moderated t

对基因 *g*：

```
y_g = X β_g + ε_g,        ε_g ~ N(0, σ_g² I)
β̂_g = (XᵀX)⁻¹ Xᵀ y_g
s_g² = ‖y_g − Xβ̂_g‖² / (n − rank(X))
```

先验方差 `s0²(·)` 是平均表达 `ā_g` 的平滑函数（分箱 + 局部线性平滑），
先验自由度 `d0` 由矩估计得到：

```
d0 = 2 · trigamma⁻¹( Var(log s_g²) − trigamma(d/2) )
```

后验方差与 moderated t：

```
σ̃_g² = (d0 · s0²(ā_g) + d · s_g²) / (d0 + d)
t_g   = β̂_g / sqrt(σ̃_g² · (XᵀX)⁻¹_jj)  ~  t_{d+d0}
```

实现上：
- `trend=True` 迭代 3 次（先验依赖均值），`trend=False` 使用常数先验；
- `robust=True` 对远低于趋势的方差设下界，避免过度收缩；
- 分箱数 `n_bins=40`，每箱最少 5 个基因，二次局部线性平滑。

## 4. 稳定性

分层 bootstrap，第 *b* 次重抽样在每组内独立有放回抽取同样本量：

```
sel_freq_g   = (1/B) Σ_b 1[ padj_{g,b} < α ∧ sign(logFC_{g,b}) = sign(logFC_g) ∧ |logFC_{g,b}| ≥ τ ]
sign_cons_g  = (1/B) Σ_b 1[ sign(logFC_{g,b}) = sign(logFC_g) ]
cv_score_g   = clip(1 − SD_b(logFC_{g,b}) / |mean_b(logFC_{g,b})|, 0, 1)
stability_g  = sel_freq_g × (0.7 × sign_cons_g + 0.3 × cv_score_g)
```

- 乘积形式使“从未被选中”的基因稳定性为 0，避免加权求和造成的区间压缩；
- 默认 `α=0.05`、`τ=1`、`B=200`；
- 同时输出 bootstrap logFC 的 5%/95% 分位区间。

## 5. 多重检验

- BH：标准 step-up 过程，NaN 安全；
- 经验零分布：以 |z| 的 90% 分位内数据估计零分布中心（中位数），
  以 `q90 / Φ⁻¹(0.95)` 估计尺度，`π0` 由中央 50% 区域的 `f_emp/f0` 中位数估计；
- local FDR：`lfdr = π0 · f0(z) / f_emp(z)`，`f_emp` 为 Silverman 带宽的高斯核密度。

## 6. 富集

### 6.1 ORA

背景 `U` 为实际检验基因集合；查询集 `Q`；基因集 `S ⊆ U`：

```
P = Hypergeom(N=|U|, K=|S|, n=|Q|).sf(|Q∩S| − 1)
```

### 6.2 稳定性加权 ORA

以背景平均稳定性 `p0` 为零假设入选概率，权重和 `W = Σ_{g∈S} w_g`：

```
E[W] = |S| · p0,   Var[W] = |S| · p0 (1 − p0)
z = (W − E[W]) / sqrt(Var[W]),   p = 1 − Φ(z)
```

`baseline="zero"` 时 `p0=0`，即检验“集合内是否存在稳定性证据”。

### 6.3 GSEA

加权 Kolmogorov–Smirnov：

```
P_hit(i) = |r_i|^p · w_i / Σ_{j∈S} |r_j|^p · w_j
P_miss(i) = −1 / (N − |S|)
ES = max_k |Σ_{i≤k} (P_hit − P_miss)|
```

`w_i` 为可选稳定性权重（默认全 1）。p 值由基因集置换（同大小随机集合）得到，
`NES = ES / mean|ES_null|`。**注意**：基因集置换未完全保留基因间相关结构，
在强连锁区域可能偏乐观，正式发表时应辅以表型置换或报告效应量。

### 6.4 综合证据

```
evidence = 0.35 · rank(−log10 padj_ORA) + 0.35 · rank(−log10 padj_GSEA) + 0.30 · selection_freq
```

`selection_freq` 为 bootstrap 排序重复中该基因集进入显著前 15 名的比例。

## 7. 生存分析

- Cox：`h(t|x) = h0(t) exp(βᵀx)`，表达按基因标准化（每 SD）；
- KM：按中位数（或分位数/最优切点）二分组，log-rank 检验；
- 联合候选：`score = |logFC| · |log HR| · stability + 0.1 · (−log10 p_DE − log10 p_Cox)`，
  仅在 `padj_DE < 0.05`、`p_Cox < 0.05` 且 `stability ≥ 0.7` 的基因中排序。

## 8. 已知局限

1. bootstrap 只反映样本抽样变异，不包含批次/平台变异；
2. 稳定性加权 ORA 使用正态近似，极小基因集或极端权重时精度下降；
3. GSEA 采用基因集置换，未完全建模基因间相关；
4. 生存分析未系统校正临床混杂因素，HR 应视为发现阶段效应量；
5. 稳定性阈值（0.7）为经验值，应结合具体数据的稳定性分布选择。
