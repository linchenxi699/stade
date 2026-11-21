"""Unit and property tests for StaDE (pytest).

Run with ``python -m pytest -q`` from the repository root.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

import stade
from stade.enrichment import GeneSet, GeneSetCollection, gsea, load_enrichr_tsv, ora, ora_weighted
from stade.normalization import cpm, filter_low_expression, log_cpm, tmm_factors
from stade.stability import bootstrap_de, empirical_null, gene_stability
from stade.stats import bh_adjust, fit_ebayes, significance, welch_test, wilcoxon_test


# ---------------------------------------------------------------------------
# fixtures
# ---------------------------------------------------------------------------
def _simulate(n_genes: int = 400, n_per_group: int = 12, n_de: int = 40, seed: int = 0) -> tuple[pd.DataFrame, pd.Series]:
    rng = np.random.default_rng(seed)
    lib = rng.uniform(4e6, 8e6, size=2 * n_per_group)
    base = rng.gamma(0.6, 120, size=n_genes)
    counts = np.outer(base, lib / 1e6)
    counts *= rng.gamma(5.0, 0.2, size=counts.shape)
    effect = np.zeros(n_genes)
    effect[:n_de] = rng.choice([-1, 1], size=n_de) * rng.uniform(1.2, 2.5, size=n_de)
    is_case = np.r_[np.ones(n_per_group, bool), np.zeros(n_per_group, bool)]
    counts[:, is_case] *= 2 ** effect[:, None]
    genes = [f"G{i:04d}" for i in range(n_genes)]
    samples = [f"S{i:02d}" for i in range(2 * n_per_group)]
    counts_df = pd.DataFrame(np.round(counts).astype(int), index=genes, columns=samples)
    groups = pd.Series(np.where(is_case, "Tumor", "Normal"), index=samples, name="group")
    return counts_df, groups


@pytest.fixture(scope="module")
def sim() -> tuple[pd.DataFrame, pd.Series]:
    return _simulate()


@pytest.fixture(scope="module")
def logexpr(sim) -> pd.DataFrame:
    counts, _ = sim
    return log_cpm(counts, use_tmm=True)


# ---------------------------------------------------------------------------
# normalisation
# ---------------------------------------------------------------------------
def test_cpm_sums_to_1e6(sim):
    counts, _ = sim
    m = cpm(counts)
    assert np.allclose(m.sum(axis=0), 1e6, rtol=1e-9)


def test_tmm_geometric_mean_is_one(sim):
    counts, _ = sim
    f = tmm_factors(counts)
    assert np.exp(np.log(f).mean()) == pytest.approx(1.0, rel=1e-9)
    assert np.isfinite(f).all()


def test_filter_keeps_expressed_genes(sim):
    counts, groups = sim
    out = filter_low_expression(counts, min_count=10, min_prop=0.5, groups=groups)
    assert 0 < out.shape[0] <= counts.shape[0]


# ---------------------------------------------------------------------------
# multiple testing
# ---------------------------------------------------------------------------
def test_bh_matches_statsmodels():
    from statsmodels.stats.multitest import multipletests

    rng = np.random.default_rng(1)
    p = rng.uniform(0, 1, 200)
    p[:20] = rng.uniform(0, 1e-4, 20)
    ours = bh_adjust(p)
    ref = multipletests(p, method="fdr_bh")[1]
    assert np.allclose(ours, ref, atol=1e-12)


def test_bh_handles_nan():
    p = np.array([0.01, np.nan, 0.5, 0.02])
    out = bh_adjust(p)
    assert np.isnan(out[1])
    assert np.all((out[[0, 2, 3]] >= 0) & (out[[0, 2, 3]] <= 1))


# ---------------------------------------------------------------------------
# DE model
# ---------------------------------------------------------------------------
def test_ebayes_recovers_planted_effects(logexpr, sim):
    _, groups = sim
    design = pd.DataFrame({"Intercept": 1.0, "group": (groups == "Tumor").astype(float)}, index=logexpr.columns)
    res = fit_ebayes(logexpr, design, "group")
    de = significance(res.table, alpha=0.05, min_abs_logfc=1.0)
    assert (de.index[:40][de.loc[de.index[:40], "significant"]]).size >= 25
    assert de.loc[de.index[100:150], "significant"].sum() <= 3
    assert res.d0 > 0


def test_ebayes_shrinks_variance(logexpr, sim):
    _, groups = sim
    design = pd.DataFrame({"Intercept": 1.0, "group": (groups == "Tumor").astype(float)}, index=logexpr.columns)
    res = fit_ebayes(logexpr, design, "group")
    assert (res.table["sigma2"] < res.table["sigma2_raw"]).mean() > 0.5


def test_ebayes_without_trend(logexpr, sim):
    _, groups = sim
    design = pd.DataFrame({"Intercept": 1.0, "group": (groups == "Tumor").astype(float)}, index=logexpr.columns)
    res = fit_ebayes(logexpr, design, "group", trend=False, robust=False)
    assert np.isfinite(res.table["t"]).all()


def test_welch_and_wilcoxon_agree_on_direction(logexpr, sim):
    _, groups = sim
    w = welch_test(logexpr, groups, "Tumor")
    m = wilcoxon_test(logexpr, groups, "Tumor")
    assert np.corrcoef(w["logFC"], m["logFC"])[0, 1] > 0.99


# ---------------------------------------------------------------------------
# stability
# ---------------------------------------------------------------------------
def test_bootstrap_reproducible_and_shaped(logexpr, sim):
    _, groups = sim
    b1 = bootstrap_de(logexpr, groups, "Tumor", n_boot=4, seed=1)
    b2 = bootstrap_de(logexpr, groups, "Tumor", n_boot=4, seed=1)
    assert b1.logfc.shape == (logexpr.shape[0], 4)
    assert np.allclose(b1.logfc.to_numpy(), b2.logfc.to_numpy())
    assert b1.tstat is not None


def test_gene_stability_ranges(logexpr, sim):
    _, groups = sim
    b = bootstrap_de(logexpr, groups, "Tumor", n_boot=6, seed=3)
    design = pd.DataFrame({"Intercept": 1.0, "group": (groups == "Tumor").astype(float)}, index=logexpr.columns)
    res = fit_ebayes(logexpr, design, "group")
    st = gene_stability(b, res.table["logFC"])
    for col in ("sel_freq", "sign_consistency", "stability"):
        assert st[col].between(0, 1).all()
    assert st.loc[st.index[:40], "stability"].median() > st.loc[st.index[100:], "stability"].median()


def test_empirical_null_calibrated(logexpr):
    rng = np.random.default_rng(0)
    t = rng.standard_t(20, size=5000)
    out = empirical_null(t)
    sd = out["null_sd"].dropna()
    assert sd.size > 0
    assert sd.between(0.5, 2.0).all()
    assert out["lfdr"].dropna().between(0, 1).all()


# ---------------------------------------------------------------------------
# enrichment
# ---------------------------------------------------------------------------
def test_hypergeometric_known_value():
    # universe of 100 genes, set of 10, query of 10 with 5 overlapping
    from scipy import stats

    expected = stats.hypergeom.sf(4, 100, 10, 10)
    universe = {f"g{i}" for i in range(100)}
    members = sorted(universe)[:10]
    gs = GeneSetCollection({"set10": GeneSet("set10", frozenset(members))})
    query = set(members[:5]) | {f"g{i}" for i in range(90, 95)}
    out = ora(query, gs, universe, min_size=5, max_size=50)
    assert out.loc[0, "overlap"] == 5
    assert out.loc[0, "pvalue"] == pytest.approx(expected, rel=1e-6)


def test_ora_detects_enriched_set():
    universe = {f"g{i}" for i in range(1000)}
    members = list(universe)[:50]
    gs = GeneSetCollection({"sig": GeneSet("sig", frozenset(members)), "noise": GeneSet("noise", frozenset(list(universe)[500:550]))})
    out = ora(set(members), gs, universe, min_size=5, max_size=500)
    assert out.iloc[0]["term"] == "sig"
    assert out.iloc[0]["padj"] < 1e-10


def test_weighted_ora_prefers_stable_genes():
    universe = {f"g{i}" for i in range(500)}
    members = list(universe)[:30]
    gs = GeneSetCollection({"sig": GeneSet("sig", frozenset(members))})
    w = pd.Series(0.05, index=sorted(universe))
    w.loc[members] = 0.98
    out = ora_weighted(w, gs, universe, min_size=5, max_size=500)
    assert out.iloc[0]["weight_ratio"] > 1.5
    assert out.iloc[0]["pvalue"] < 1e-6
    out0 = ora_weighted(w, gs, universe, min_size=5, max_size=500, baseline="zero")
    assert out0.iloc[0]["pvalue"] < 1e-6


def test_gsea_detects_leading_block():
    genes = [f"g{i}" for i in range(300)]
    ranking = pd.Series(np.r_[np.linspace(5, 1, 20), np.linspace(0.9, -0.9, 260), np.linspace(-1, -5, 20)], index=genes)
    gs = GeneSetCollection({"top": GeneSet("top", frozenset(genes[:15])), "bottom": GeneSet("bottom", frozenset(genes[-15:]))})
    out = gsea(ranking, gs, n_perm=200, min_size=5, max_size=100, seed=3)
    assert out.iloc[0]["term"] in {"top", "bottom"}
    assert out["ES"].abs().max() > 0.8


def test_enrichr_reader(tmp_path):
    p = tmp_path / "lib.txt"
    p.write_text("Pathway A\t\tTP53\tEGFR\nPathway B (GO:0001234)\t\tMYC\tEGFR\n", encoding="utf-8")
    gs = load_enrichr_tsv(p)
    assert len(gs) == 2
    assert gs["Pathway A"].genes == frozenset({"TP53", "EGFR"})
    assert gs["Pathway B"].description == "GO:0001234"


# ---------------------------------------------------------------------------
# provenance
# ---------------------------------------------------------------------------
def test_provenance_roundtrip(tmp_path):
    import json

    (tmp_path / "x.txt").write_text("hello", encoding="utf-8")
    prov = stade.Provenance(command="test", parameters={"alpha": 0.05})
    prov.add_input("fixture", tmp_path / "x.txt", with_hash=False)
    prov.finish(ok=True)
    prov.write(tmp_path / "p.json")
    obj = json.loads((tmp_path / "p.json").read_text(encoding="utf-8"))
    assert obj["command"] == "test"
    assert obj["parameters"]["alpha"] == 0.05
    assert obj["environment"]["python"]
