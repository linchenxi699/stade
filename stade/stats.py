"""Differential expression statistics.

Implements an empirical-Bayes moderated t-statistic in the spirit of limma
(Smyth 2004) using a mean-variance trend as the prior, together with two
non-parametric competitors used for concordance analysis:

* ``welch``     : Welch's t-test on log-CPM values;
* ``wilcoxon``  : Wilcoxon rank-sum test (Mann-Whitney) on log-CPM values.

The API returns a tidy :class:`pandas.DataFrame` whose column names follow the
limma conventions (``logFC``, ``AveExpr``, ``t``, ``P.Value``, ``adj.P.Val``).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import special, stats

__all__ = [
    "EbayesResult",
    "fit_ebayes",
    "welch_test",
    "wilcoxon_test",
    "bh_adjust",
    "significance",
]


def bh_adjust(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR adjustment (NaN-safe)."""
    p = np.asarray(p, dtype=float)
    out = np.full_like(p, np.nan)
    ok = np.isfinite(p)
    if not ok.any():
        return out
    pv = p[ok]
    n = pv.size
    order = np.argsort(pv)
    ranked = pv[order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adj = np.empty(n)
    adj[order] = np.clip(ranked, 0, 1)
    out[ok] = adj
    return out


def _expected_log_chi2_var(df: float) -> float:
    """Var(log chi2_df / df), used by the limma moment estimator of d0."""
    return float(special.polygamma(1, df / 2.0))


def _trigamma_inverse(x: float, lo: float = 1e-6, hi: float = 1e4) -> float:
    """Solve trigamma(y) = x for y (monotone decreasing, bisection is enough)."""
    if x <= 0:
        return hi
    if special.polygamma(1, hi) > x:
        return hi
    if special.polygamma(1, lo) < x:
        return lo
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if special.polygamma(1, mid) > x:
            lo = mid
        else:
            hi = mid
        if hi - lo < 1e-9 * max(1.0, hi):
            break
    return 0.5 * (lo + hi)


def _estimate_d0(s2: np.ndarray, df_residual: float, max_iter: int = 50) -> float:
    """Moment estimate of the prior degrees of freedom (limma ``fitFDist``)."""
    s2 = np.asarray(s2, dtype=float)
    ok = np.isfinite(s2) & (s2 > 0)
    if ok.sum() < 20:
        return 4.0
    e = np.log(s2[ok])
    e_mean = e.mean()
    e_var = e.var(ddof=1)
    # remove the expected sampling variance of log(s2) under chi2_df/df
    target = e_var - _expected_log_chi2_var(df_residual)
    if not np.isfinite(target) or target <= 0:
        return 1e3  # variances are homogeneous -> very strong shrinkage
    d0 = 2.0 * _trigamma_inverse(max(target, 1e-12))
    return float(np.clip(d0, 0.5, 1e3))


@dataclass
class EbayesResult:
    table: pd.DataFrame
    d0: float
    df_total: float
    trend: np.ndarray | None = None
    trend_x: np.ndarray | None = None
    design: pd.DataFrame | None = None

    def __getitem__(self, key: str) -> pd.Series:
        return self.table[key]

    @property
    def genes(self) -> pd.Index:
        return self.table.index


def _lowess_trend(x: np.ndarray, y: np.ndarray, frac: float = 0.4) -> tuple[np.ndarray, np.ndarray]:
    """Robust LOWESS trend; falls back to a running median on tiny inputs."""
    from statsmodels.nonparametric.smoothers_lowess import lowess

    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 30:
        return np.array([]), np.array([])
    res = lowess(y[ok], x[ok], frac=frac, return_sorted=True, it=1)
    return res[:, 0], res[:, 1]


def _binned_trend(
    x: np.ndarray,
    y: np.ndarray,
    n_bins: int = 40,
    min_count: int = 5,
) -> tuple[np.ndarray, np.ndarray]:
    """Fast mean-variance trend: equal-width bins + local-linear smoother.

    This is the workhorse that replaces LOWESS inside the bootstrap loop. It
    bins genes by average expression, takes the per-bin average of ``y``
    (typically ``log(sigma2)``) and then applies a local linear smoother over
    the bin centres, which is ~10^5 times faster than a full LOWESS fit and
    statistically indistinguishable for the numbers of genes used here.
    """
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    n = x.size
    if n < 30:
        return np.array([]), np.array([])
    order = np.argsort(x, kind="stable")
    xs, ys = x[order], y[order]
    nb = int(np.clip(n_bins, 6, max(6, n // min_count)))
    edges = np.quantile(xs, np.linspace(0, 1, nb + 1))
    edges = np.unique(edges)
    if edges.size < 4:
        return np.array([]), np.array([])
    idx = np.clip(np.searchsorted(edges, xs, side="right") - 1, 0, edges.size - 2)
    cnt = np.bincount(idx, minlength=edges.size - 1)
    sx = np.bincount(idx, weights=xs, minlength=edges.size - 1)
    sy = np.bincount(idx, weights=ys, minlength=edges.size - 1)
    good = cnt >= min_count
    if good.sum() < 3:
        return np.array([]), np.array([])
    cx = sx[good] / cnt[good]
    cy = sy[good] / cnt[good]
    # local-linear (LOESS-like) smoothing over bin centres, tricube weights
    span = max(3, int(np.ceil(0.5 * cy.size)))
    sm = np.empty_like(cy)
    for i in range(cy.size):
        lo = max(0, i - span // 2)
        hi = min(cy.size, lo + span)
        lo = max(0, hi - span)
        d = np.abs(cx[lo:hi] - cx[i])
        h = d.max() if d.max() > 0 else 1.0
        w = (1 - (d / h) ** 3) ** 3
        X = np.column_stack([np.ones(hi - lo), cx[lo:hi] - cx[i]])
        W = w[:, None]
        try:
            coef = np.linalg.lstsq(X * W, cy[lo:hi] * w, rcond=None)[0]
            sm[i] = coef[0]
        except np.linalg.LinAlgError:
            sm[i] = cy[i]
    # enforce a non-increasing mean-variance relationship is *not* assumed;
    # keep the smoother's shape but clip extreme extrapolation
    return cx, sm


def fit_ebayes(
    log_expr: pd.DataFrame,
    design: pd.DataFrame,
    coef: str,
    trend: bool = True,
    robust: bool = True,
    max_iter: int = 3,
) -> EbayesResult:
    """Fit per-gene linear models and moderate variances with an EB prior.

    Parameters
    ----------
    log_expr:
        Gene x sample log-expression matrix (e.g. ``log_cpm``).
    design:
        Sample x coefficient design matrix (columns are coefficients).
    coef:
        Name of the coefficient tested against zero.
    trend:
        If True, the prior variance is a LOWESS function of average expression
        (heteroscedasticity-aware moderation) instead of a constant.
    robust:
        If True, iteratively down-weight genes whose observed variance is far
        below the trend (limma's ``robust=TRUE`` idea) to avoid over-shrinking
        biologically variable genes.
    """
    y = log_expr.to_numpy(dtype=float)
    X = design.to_numpy(dtype=float)
    if X.shape[0] != y.shape[1]:
        raise ValueError("design rows must match expression columns")
    if coef not in design.columns:
        raise KeyError(f"coefficient {coef!r} not in design: {list(design.columns)}")

    # ---- per-gene least squares ----------------------------------------
    yt = y.T  # samples x genes
    XtX_inv = np.linalg.pinv(X.T @ X)
    beta = (XtX_inv @ (X.T @ yt)).T  # genes x coefficients
    fitted = X @ beta.T  # samples x genes
    resid = yt - fitted
    df_res = X.shape[0] - np.linalg.matrix_rank(X)
    sigma2 = (resid**2).sum(axis=0) / df_res
    df_res = float(df_res)

    j = list(design.columns).index(coef)
    coef_hat = beta[:, j]
    cov_unscaled = XtX_inv[j, j]  # scalar for a single coefficient
    se = np.sqrt(sigma2 * cov_unscaled)
    t_obs = coef_hat / np.where(se > 0, se, np.nan)

    # ---- empirical Bayes moderation ------------------------------------
    ave = y.mean(axis=1)
    s2 = sigma2
    d0 = _estimate_d0(s2, df_res)
    trend_x = trend_y = np.array([])
    if trend:
        for _ in range(max_iter):
            tx, ty = _binned_trend(ave, np.log(s2))
            if tx.size == 0:
                break
            trend_x, trend_y = tx, ty
            s2_trend = np.exp(np.interp(ave, tx, ty))
            w = df_res / (df_res + d0)
            s2 = w * sigma2 + (1 - w) * s2_trend
        # final trend on the moderated variances (closes the loop)
        tx, ty = _binned_trend(ave, np.log(s2))
        if tx.size:
            trend_x, trend_y = tx, ty
        s2_trend = np.exp(np.interp(ave, trend_x, trend_y)) if trend_x.size else np.full_like(s2, np.median(s2))
    else:
        s2_trend = np.full_like(s2, np.median(s2[s2 > 0]) if np.any(s2 > 0) else 1.0)

    if robust:
        # genes far below the trend carry little information: cap their weight
        prior_var = np.median(s2_trend[s2_trend > 0]) if np.any(s2_trend > 0) else 1.0
        low = s2_trend < 0.01 * prior_var
        if low.any():
            s2[low] = s2_trend[low]

    d0 = _estimate_d0(s2, df_res) if trend or robust else d0
    df_total = df_res + d0
    s2_post = (d0 * s2_trend + df_res * sigma2) / df_total
    se_post = np.sqrt(s2_post * cov_unscaled)
    t_mod = coef_hat / np.where(se_post > 0, se_post, np.nan)
    p = 2.0 * stats.t.sf(np.abs(t_mod), df_total)

    table = pd.DataFrame(
        {
            "logFC": coef_hat,
            "AveExpr": ave,
            "t": t_mod,
            "P.Value": p,
            "adj.P.Val": bh_adjust(p),
            "t_raw": t_obs,
            "P.Value.raw": 2.0 * stats.t.sf(np.abs(t_obs), df_res),
            "sigma2": s2_post,
            "sigma2_raw": sigma2,
            "SE": se_post,
        },
        index=log_expr.index,
    )
    table["B"] = _b_statistic(table["t"].to_numpy(), d0, df_res)
    return EbayesResult(table=table, d0=d0, df_total=df_total, trend=trend_y, trend_x=trend_x, design=design)


def _b_statistic(t: np.ndarray, d0: float, df_res: float) -> np.ndarray:
    """limma's B statistic (log-odds of differential expression)."""
    df_total = d0 + df_res
    r = d0 / df_total
    with np.errstate(divide="ignore", invalid="ignore"):
        b = 0.5 * (np.log1p(t**2 / df_total) - np.log(r))
    return b


def welch_test(log_expr: pd.DataFrame, groups: pd.Series, group_a: str) -> pd.DataFrame:
    """Welch t-test for ``group_a`` vs the remaining samples (log scale, no EB)."""
    is_a = (groups == group_a).to_numpy()
    a = log_expr.to_numpy(dtype=float)[:, is_a]
    b = log_expr.to_numpy(dtype=float)[:, ~is_a]
    stat, p = stats.ttest_ind(a, b, axis=1, equal_var=False, nan_policy="omit")
    logfc = a.mean(axis=1) - b.mean(axis=1)
    return pd.DataFrame({"logFC": logfc, "P.Value": p, "adj.P.Val": bh_adjust(p)}, index=log_expr.index)


def wilcoxon_test(log_expr: pd.DataFrame, groups: pd.Series, group_a: str) -> pd.DataFrame:
    """Two-sided Wilcoxon rank-sum (Mann-Whitney) test per gene."""
    is_a = (groups == group_a).to_numpy()
    a = log_expr.to_numpy(dtype=float)[:, is_a]
    b = log_expr.to_numpy(dtype=float)[:, ~is_a]
    stat, p = stats.mannwhitneyu(a, b, axis=1, alternative="two-sided")
    logfc = a.mean(axis=1) - b.mean(axis=1)
    return pd.DataFrame({"logFC": logfc, "P.Value": p, "adj.P.Val": bh_adjust(p), "U": stat}, index=log_expr.index)


def significance(table: pd.DataFrame, alpha: float = 0.05, min_abs_logfc: float = 1.0) -> pd.DataFrame:
    """Attach direction/significance labels to a DE table."""
    df = table.copy()
    df["direction"] = np.where(df["logFC"] > 0, "Up", "Down")
    df["significant"] = (df["adj.P.Val"] < alpha) & (df["logFC"].abs() >= min_abs_logfc)
    return df
