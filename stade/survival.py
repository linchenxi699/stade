"""Univariate/multivariate survival analysis for candidate genes.

Thin, explicit wrappers around ``lifelines`` plus the data plumbing needed for
TCGA-style clinical tables (``OS.time`` in days, ``OS`` event indicator).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test

__all__ = [
    "prepare_survival",
    "cut_by_expression",
    "km_analysis",
    "cox_univariate",
    "cox_multivariate",
    "candidate_signature",
]


def prepare_survival(
    clinical: pd.DataFrame,
    expression: pd.DataFrame,
    time_col: str = "OS.time",
    event_col: str = "OS",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Intersect clinical follow-up and expression samples; return both tables."""
    clin = clinical.copy()
    if time_col not in clin or event_col not in clin:
        raise KeyError(f"clinical table needs {time_col!r} and {event_col!r}")
    clin = clin[[time_col, event_col] + [c for c in clin.columns if c not in (time_col, event_col)]]
    clin = clin.dropna(subset=[time_col, event_col])
    clin[time_col] = pd.to_numeric(clin[time_col], errors="coerce")
    clin[event_col] = pd.to_numeric(clin[event_col], errors="coerce")
    clin = clin.dropna(subset=[time_col, event_col])
    clin = clin[clin[time_col] > 0]
    samples = [s for s in expression.columns if s in clin.index]
    clin = clin.loc[samples]
    expr = expression.loc[:, samples]
    return clin, expr


def cut_by_expression(expr: pd.Series, method: str = "median", q: float = 0.5) -> pd.Series:
    """Dichotomise a continuous gene score into High/Low groups."""
    if method == "median":
        thr = expr.median()
    elif method == "quartile":
        thr = expr.quantile(q)
    else:
        raise ValueError(f"unknown cut method {method!r}; use km_analysis(cut='optimal') for scanning")
    return pd.Series(np.where(expr > thr, "High", "Low"), index=expr.index, name="group")


@dataclass
class KMResult:
    group: pd.Series
    logrank_p: float
    median_high: float
    median_low: float
    n_high: int
    n_low: int
    hr: float
    hr_ci: tuple[float, float]


def km_analysis(
    time: pd.Series,
    event: pd.Series,
    score: pd.Series,
    cut: str = "median",
    minutes: bool = True,
) -> KMResult:
    """Kaplan-Meier grouping + log-rank test (and the Cox HR as a by-product)."""
    score = score.reindex(time.index)
    if cut == "optimal":
        opt = _scan_optimal(time, event, score)
        group = pd.Series(np.where(score > opt.threshold, "High", "Low"), index=time.index)
    else:
        group = cut_by_expression(score.dropna(), method=cut)
    t = time.reindex(group.index)
    e = event.reindex(group.index).astype(int)
    hi = group == "High"
    lr = logrank_test(t[hi], t[~hi], e[hi], e[~hi])
    kmf = KaplanMeierFitter()
    kmf.fit(t[hi], e[hi])
    mh = float(kmf.median_survival_time_)
    kmf.fit(t[~hi], e[~hi])
    ml = float(kmf.median_survival_time_)
    cph = CoxPHFitter()
    df = pd.DataFrame({"T": t, "E": e, "High": hi.astype(float)})
    try:
        cph.fit(df, duration_col="T", event_col="E")
        hr = float(np.exp(cph.params_["High"]))
        ci = tuple(float(x) for x in np.exp(cph.confidence_intervals_.loc["High"]).to_numpy())
    except Exception:  # noqa: BLE001 - no events or singular
        hr, ci = np.nan, (np.nan, np.nan)
    scale = 1 / 30.425 if minutes else 1.0  # days -> months
    return KMResult(group, float(lr.p_value), mh * scale, ml * scale, int(hi.sum()), int((~hi).sum()), hr, ci)


@dataclass
class _OptCut:
    threshold: float
    p: float


def _scan_optimal(time: pd.Series, event: pd.Series, score: pd.Series, min_frac: float = 0.25) -> _OptCut:
    vals = np.sort(score.dropna().unique())
    lo = int(np.floor(len(vals) * min_frac))
    hi = int(np.ceil(len(vals) * (1 - min_frac)))
    if hi <= lo:
        return _OptCut(float(np.median(score.dropna())), 1.0)
    cand = vals[lo:hi: max(1, (hi - lo) // 80)]
    best = _OptCut(float(np.median(vals)), 1.0)
    for thr in cand:
        hi_mask = score > thr
        if hi_mask.sum() < 5 or (~hi_mask).sum() < 5:
            continue
        p = logrank_test(
            time[hi_mask], time[~hi_mask], event[hi_mask], event[~hi_mask]
        ).p_value
        if p < best.p:
            best = _OptCut(float(thr), float(p))
    return best


def cox_univariate(
    clinical: pd.DataFrame,
    expr: pd.DataFrame,
    genes: list[str] | None = None,
    time_col: str = "OS.time",
    event_col: str = "OS",
    zscore: bool = True,
) -> pd.DataFrame:
    """Per-gene univariate Cox model on z-scored expression."""
    genes = list(expr.index) if genes is None else [g for g in genes if g in expr.index]
    rows = []
    for g in genes:
        x = expr.loc[g].astype(float)
        if zscore:
            sd = x.std()
            x = (x - x.mean()) / (sd if sd > 0 else 1.0)
        df = pd.DataFrame({"T": clinical[time_col].values, "E": clinical[event_col].astype(int).values, "x": x.values}).dropna()
        if df["E"].sum() < 5 or df["x"].std() == 0:
            continue
        cph = CoxPHFitter()
        try:
            cph.fit(df, duration_col="T", event_col="E")
            s = cph.summary.loc["x"]
            rows.append(
                {
                    "gene": g,
                    "HR": float(np.exp(s["coef"])),
                    "HR_lower": float(np.exp(s["coef lower 95%"])),
                    "HR_upper": float(np.exp(s["coef upper 95%"])),
                    "z": float(s["z"]),
                    "pvalue": float(s["p"]),
                    "n": int(len(df)),
                    "events": int(df["E"].sum()),
                }
            )
        except Exception:  # noqa: BLE001
            continue
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    from .stats import bh_adjust

    out["padj"] = bh_adjust(out["pvalue"].to_numpy())
    return out.sort_values("pvalue").reset_index(drop=True)


def cox_multivariate(
    clinical: pd.DataFrame,
    expr: pd.DataFrame,
    genes: list[str],
    time_col: str = "OS.time",
    event_col: str = "OS",
    covariates: list[str] | None = None,
) -> tuple[pd.DataFrame, object]:
    """Multivariate Cox model for a small gene panel (+ optional covariates)."""
    genes = [g for g in genes if g in expr.index]
    data = expr.loc[genes].T.astype(float)
    data = (data - data.mean()) / data.std().replace(0, 1.0)
    data.columns = [f"expr_{g}" for g in genes]
    data["T"] = clinical[time_col].reindex(data.index)
    data["E"] = clinical[event_col].reindex(data.index).astype(float)
    if covariates:
        for c in covariates:
            if c in clinical.columns:
                data[c] = clinical[c].reindex(data.index)
    data = data.dropna()
    cph = CoxPHFitter(penalizer=0.05)
    cph.fit(data, duration_col="T", event_col="E")
    return cph.summary, cph


def candidate_signature(
    de: pd.DataFrame,
    surv: pd.DataFrame,
    top_n: int = 12,
    require_stability: float = 0.6,
    de_alpha: float = 0.05,
    cox_alpha: float = 0.05,
) -> pd.DataFrame:
    """Joint ranking of DE and prognostic evidence (discovery-level only).

    Only genes that are both differentially expressed (``adj.P.Val`` <
    ``de_alpha``) and prognostic (``cox_p`` < ``cox_alpha``) are eligible; the
    ranking score combines the two -log10(p) values.
    """
    a = de[["logFC", "adj.P.Val", "stability"]].copy() if "stability" in de else de[["logFC", "adj.P.Val"]].copy()
    a["direction"] = np.where(a["logFC"] > 0, "Up", "Down")
    b = surv.set_index("gene")[["HR", "pvalue", "padj"]].rename(columns={"pvalue": "cox_p", "padj": "cox_padj"})
    m = a.join(b, how="inner")
    if "stability" in m:
        m = m[m["stability"] >= require_stability]
    m = m[(m["adj.P.Val"] < de_alpha) & (m["cox_p"] < cox_alpha)]
    # effect-size aware ranking: differential magnitude x prognostic magnitude,
    # with a stability multiplier; the two -log10(p) terms keep significance in
    # the ranking but no longer dominate it.
    stab = m["stability"] if "stability" in m else pd.Series(1.0, index=m.index)
    m["score"] = (m["logFC"].abs() * np.log(m["HR"].clip(lower=1e-6)).abs() * stab) + 0.1 * (
        -np.log10(np.clip(m["adj.P.Val"], 1e-300, 1)) + -np.log10(np.clip(m["cox_p"], 1e-300, 1))
    )
    return m.sort_values("score", ascending=False).head(top_n)
