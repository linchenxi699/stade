# Changelog

All notable changes to StaDE are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project uses
[Semantic Versioning](https://semver.org/).

## [0.1.0] - 2026-12-15

### Added

- Empirical-Bayes moderated t-test (`fit_ebayes`) with a binned local-linear
  mean-variance trend prior, robust variance flooring and limma-style B statistic.
- TMM normalisation, log-CPM transformation, low-expression filtering and
  library-size outlier detection.
- Stratified bootstrap differential expression (`bootstrap_de`) and the composite
  gene stability score (`gene_stability`): selection frequency x (0.7 x sign
  consistency + 0.3 x effect-size consistency).
- BH correction, group-label permutation null and Efron-style empirical null with
  local FDR.
- Method concordance table (limma-EB / Welch / Wilcoxon).
- Over-representation analysis, stability-weighted ORA with a Poisson-binomial
  approximation, weighted Kolmogorov-Smirnov GSEA, bootstrap gene-set selection
  frequency and an integrated evidence score.
- Univariate/multivariate Cox models, Kaplan-Meier + log-rank, joint DE x prognosis
  candidate ranking.
- GEO series-matrix parser and cross-cohort validation helpers.
- `stade` command line interface: `prep`, `de`, `enrich`, `surv`.
- Provenance JSON with input SHA-256, parameters, environment and timestamps.
- 18 unit/property tests, a simulated end-to-end demo and a bundled demo dataset.
- Publication-grade figures: PCA, MA, volcano, stability, heatmap, enrichment dot
  plot, KM curves, forest plot, validation scatter.

### Known limitations

- Bootstrap captures sampling variability only (not batch or platform effects).
- Stability-weighted ORA relies on a normal approximation for the Poisson-binomial.
- GSEA uses gene-set permutation and may be optimistic for co-regulated genes.
- Survival models do not adjust for clinical confounders by default.
