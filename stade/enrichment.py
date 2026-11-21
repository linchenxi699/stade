"""Gene-set enrichment: over-representation, GSEA and stability-aware variants.

Implemented from scratch (no gseapy/clusterProfiler dependency) so that the
stability extensions can reach into the enrichment statistic:

* ``load_gmt`` / ``load_enrichr_tsv`` : library readers;
* ``ora``        : hypergeometric over-representation with BH correction;
* ``ora_weighted``: stability-weighted ORA using a Poisson-binomial null;
* ``gsea``       : weighted Kolmogorov-Smirnov statistic with gene-set
  permutation p-values and BH correction;
* ``bootstrap_enrichment`` : bootstrap replicates of the ranking, producing a
  per-gene-set selection frequency and an integrated ``evidence`` score.
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from .stats import bh_adjust

__all__ = [
    "GeneSet",
    "load_enrichr_tsv",
    "load_gmt",
    "library_stats",
    "ora",
    "ora_weighted",
    "gsea",
    "bootstrap_enrichment",
    "integrated_enrichment",
    "GeneSetCollection",
]


@dataclass
class GeneSet:
    name: str
    genes: frozenset[str]
    description: str = ""
    source: str = ""


class GeneSetCollection:
    """A named collection of gene sets with fast membership helpers."""

    def __init__(self, sets: dict[str, GeneSet], source: str = ""):
        self.sets = dict(sets)
        self.source = source

    def __len__(self) -> int:
        return len(self.sets)

    def __iter__(self):
        return iter(self.sets.values())

    def __getitem__(self, key: str) -> GeneSet:
        return self.sets[key]

    def filter_sizes(self, min_size: int = 10, max_size: int = 500) -> "GeneSetCollection":
        return GeneSetCollection(
            {k: v for k, v in self.sets.items() if min_size <= len(v.genes) <= max_size}, self.source
        )

    def restrict_universe(self, universe: set[str]) -> "GeneSetCollection":
        """Intersect every set with the analysed universe (genes tested)."""
        out = {}
        for name, gs in self.sets.items():
            inter = gs.genes & universe
            if inter:
                out[name] = GeneSet(gs.name, frozenset(inter), gs.description, gs.source)
        return GeneSetCollection(out, self.source)

    def sizes(self) -> pd.Series:
        return pd.Series({k: len(v.genes) for k, v in self.sets.items()})


def load_gmt(path: str | Path) -> GeneSetCollection:
    """Read a ``.gmt`` file: ``NAME<TAB>description<TAB>gene<TAB>gene...``."""
    sets: dict[str, GeneSet] = {}
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            name, desc = parts[0].strip(), parts[1].strip()
            genes = frozenset(g.strip().upper() for g in parts[2:] if g.strip())
            if genes:
                sets[name] = GeneSet(name, genes, desc, str(path))
    return GeneSetCollection(sets, str(path))


def load_enrichr_tsv(path: str | Path) -> GeneSetCollection:
    """Read an Enrichr ``geneSetLibrary?mode=text`` dump.

    Each line is ``Term<TAB><TAB>GENE<TAB>GENE...``. GO terms carry their GO id
    in parentheses; it is kept as the description.
    """
    sets: dict[str, GeneSet] = {}
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            name = parts[0].strip()
            genes = frozenset(g.strip().upper() for g in parts[2:] if g.strip())
            if not name or not genes:
                continue
            desc = ""
            m = re.search(r"\((GO:\d+)\)", name)
            if m:
                desc = m.group(1)
                name = name[: m.start()].strip()
            sets[name] = GeneSet(name, genes, desc, str(path))
    return GeneSetCollection(sets, str(path))


def library_stats(gs: GeneSetCollection) -> pd.DataFrame:
    sizes = gs.sizes()
    return pd.DataFrame(
        {
            "library": [Path(gs.source).stem if gs.source else "memory"],
            "n_sets": [len(gs)],
            "min_size": [int(sizes.min()) if len(sizes) else 0],
            "median_size": [float(sizes.median()) if len(sizes) else 0.0],
            "max_size": [int(sizes.max()) if len(sizes) else 0],
        }
    )


# ---------------------------------------------------------------------------
# Over-representation analysis
# ---------------------------------------------------------------------------
def ora(
    query_genes: set[str],
    gene_sets: GeneSetCollection,
    universe: set[str],
    direction: str = "Up",
    min_size: int = 5,
    max_size: int = 500,
) -> pd.DataFrame:
    """Hypergeometric over-representation analysis with BH correction.

    ``universe`` should be the set of genes that were actually testable (the
    filtered expression matrix), not the whole genome.
    """
    universe = set(universe)
    q = {g for g in query_genes if g in universe}
    N = len(universe)
    n = len(q)
    rows = []
    for gs in gene_sets:
        S = gs.genes & universe
        if not (min_size <= len(S) <= max_size):
            continue
        k = len(q & S)
        if k == 0:
            p = 1.0
        else:
            p = float(stats.hypergeom.sf(k - 1, N, len(S), n))
        expected = n * len(S) / N if N else np.nan
        rows.append(
            {
                "term": gs.name,
                "description": gs.description,
                "direction": direction,
                "set_size": len(S),
                "overlap": k,
                "query_size": n,
                "universe_size": N,
                "expected": expected,
                "fold_enrichment": (k / expected) if expected else np.nan,
                "odds_ratio": ((k + 0.5) * (N - len(S) - n + k + 0.5)) / ((len(S) - k + 0.5) * (n - k + 0.5)),
                "pvalue": p,
                "genes": ";".join(sorted(q & S)),
            }
        )
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    df["padj"] = bh_adjust(df["pvalue"].to_numpy())
    return df.sort_values(["pvalue", "term"]).reset_index(drop=True)


def ora_weighted(
    weights: pd.Series,
    gene_sets: GeneSetCollection,
    universe: set[str],
    direction: str = "Up",
    min_size: int = 5,
    max_size: int = 500,
    baseline: str | None = None,
) -> pd.DataFrame:
    """Stability-weighted over-representation.

    Instead of counting significant genes equally, each gene contributes its
    stability score ``w_i`` in [0, 1]. The null distribution of the weighted
    overlap is Poisson-binomial with per-gene success probabilities ``w^null_i``:
    by default the *background* mean stability, i.e. a gene set is only
    enriched if its members are more stable than a random gene. Passing
    ``baseline="zero"`` tests against a null of no stability at all.
    """
    universe = set(universe)
    w = weights.reindex(sorted(universe)).dropna()
    w = w[(w >= 0) & (w <= 1)]
    if w.empty:
        return pd.DataFrame()
    W = w.to_numpy()
    genes = w.index.to_numpy()
    N = W.size
    p_null = 0.0 if baseline == "zero" else float(W.mean())
    rows = []
    for gs in gene_sets:
        mask = np.isin(genes, list(gs.genes))
        k = int(mask.sum())
        if not (min_size <= k <= max_size):
            continue
        obs = float(W[mask].sum())
        mean = float(k * p_null)
        var = float(k * p_null * (1 - p_null))
        if var <= 0:
            p = 1.0 if obs <= mean else 0.0
        else:
            p = float(stats.norm.sf((obs - mean) / np.sqrt(var)))
        rows.append(
            {
                "term": gs.name,
                "description": gs.description,
                "direction": direction,
                "set_size": k,
                "weighted_overlap": obs,
                "expected_weighted": mean,
                "weight_ratio": obs / mean if mean > 0 else np.nan,
                "pvalue": p,
                "genes": ";".join(sorted(genes[mask])),
            }
        )
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    df["padj"] = bh_adjust(df["pvalue"].to_numpy())
    return df.sort_values(["pvalue", "term"]).reset_index(drop=True)


# ---------------------------------------------------------------------------
# GSEA (weighted Kolmogorov-Smirnov)
# ---------------------------------------------------------------------------
def _es_for_set(
    ranks: np.ndarray,
    hit_mask: np.ndarray,
    weights_p: float = 1.0,
    gene_weights: np.ndarray | None = None,
) -> tuple[float, int, np.ndarray]:
    """Enrichment score of one gene set against a ranked statistic vector."""
    hits = np.flatnonzero(hit_mask)
    n_hit = hits.size
    n_miss = ranks.size - n_hit
    if n_hit == 0 or n_miss == 0:
        return 0.0, 0, np.zeros(0)
    w = np.abs(ranks[hits]) ** weights_p
    if gene_weights is not None:
        w = w * np.clip(gene_weights[hits], 1e-3, None)
    w = w / w.sum() if w.sum() > 0 else np.full(n_hit, 1.0 / n_hit)
    step = np.zeros(ranks.size)
    step[hits] = w
    step -= (~hit_mask) / n_miss
    run = np.cumsum(step)
    k = int(np.argmax(np.abs(run)))
    return float(run[k]), k, run


def gsea(
    ranking: pd.Series,
    gene_sets: GeneSetCollection,
    n_perm: int = 1000,
    min_size: int = 10,
    max_size: int = 500,
    seed: int = 11,
    weights_p: float = 1.0,
    gene_weights: pd.Series | None = None,
) -> pd.DataFrame:
    """Gene-set enrichment analysis with gene-set permutation p-values.

    ``ranking`` maps gene -> statistic (e.g. moderated t). Higher absolute
    values are more "extreme". Permuting gene labels of the ranking preserves
    the correlation structure between genes approximately and is much cheaper
    than re-fitting the model, which makes thousands of permutations feasible
    for large libraries.

    When ``gene_weights`` (e.g. stability scores) is given, the hit weights are
    multiplied by it, which down-weights unstable leading-edge genes - the
    "stability-aware GSEA" variant.
    """
    r = ranking.dropna().sort_values(ascending=False)
    genes = r.index.to_numpy()
    vals = r.to_numpy(dtype=float)
    rng = np.random.default_rng(seed)
    gw = None if gene_weights is None else gene_weights.reindex(genes).fillna(1.0).to_numpy()
    rows = []
    for gs in gene_sets:
        mask = np.isin(genes, list(gs.genes))
        size = int(mask.sum())
        if not (min_size <= size <= max_size):
            continue
        es, k, run = _es_for_set(vals, mask, weights_p, gw)
        # null: random sets of the same size
        null = np.empty(n_perm)
        pool = np.arange(vals.size)
        for p in range(n_perm):
            idx = rng.choice(pool, size=size, replace=False)
            m = np.zeros(vals.size, dtype=bool)
            m[idx] = True
            null[p], _, _ = _es_for_set(vals, m, weights_p, gw)
        if es >= 0:
            pval = (1 + np.sum(null >= es)) / (n_perm + 1)
            direction = "Up"
        else:
            pval = (1 + np.sum(null <= es)) / (n_perm + 1)
            direction = "Down"
        nes = es / np.mean(np.abs(null)) if np.mean(np.abs(null)) > 0 else np.nan
        hit_genes = genes[mask]
        # leading edge: genes before the running-sum peak
        order = np.flatnonzero(mask)
        if es >= 0:
            le = hit_genes[order <= k]
        else:
            le = hit_genes[order >= k]
        rows.append(
            {
                "term": gs.name,
                "description": gs.description,
                "direction": direction,
                "set_size": size,
                "ES": es,
                "NES": nes,
                "pvalue": float(pval),
                "n_perm": n_perm,
                "leading_edge": ";".join(le[:50]),
            }
        )
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    df["padj"] = bh_adjust(df["pvalue"].to_numpy())
    return df.sort_values(["pvalue", "term"]).reset_index(drop=True)


# ---------------------------------------------------------------------------
# Stability-aware enrichment
# ---------------------------------------------------------------------------
def bootstrap_enrichment(
    rankings: pd.DataFrame,
    gene_sets: GeneSetCollection,
    alpha: float = 0.05,
    n_perm: int = 500,
    min_size: int = 10,
    max_size: int = 500,
    seed: int = 13,
    top_n: int = 15,
) -> pd.DataFrame:
    """Gene-set selection frequency across bootstrap replicate rankings.

    ``rankings`` is genes x replicates (moderated t of each bootstrap fit). For
    every replicate a fast gene-set permutation GSEA is run and the top ``top_n``
    sets (by nominal p-value) are recorded. The output reports how often each
    set is selected and the median normalised enrichment score.
    """
    rng = np.random.default_rng(seed)
    counts: dict[str, int] = {}
    nes_store: dict[str, list[float]] = {}
    n_rep = rankings.shape[1]
    for b in range(n_rep):
        r = rankings.iloc[:, b].dropna().sort_values(ascending=False)
        genes = r.index.to_numpy()
        vals = r.to_numpy(dtype=float)
        res = []
        for gs in gene_sets:
            mask = np.isin(genes, list(gs.genes))
            size = int(mask.sum())
            if not (min_size <= size <= max_size):
                continue
            es, _, _ = _es_for_set(vals, mask)
            null = np.empty(n_perm)
            for p in range(n_perm):
                idx = rng.choice(vals.size, size=size, replace=False)
                m = np.zeros(vals.size, dtype=bool)
                m[idx] = True
                null[p], _, _ = _es_for_set(vals, m)
            pv = (1 + np.sum(null >= es)) / (n_perm + 1) if es >= 0 else (1 + np.sum(null <= es)) / (n_perm + 1)
            scale = np.mean(np.abs(null))
            res.append((gs.name, pv, es / scale if scale > 0 else np.nan))
        res.sort(key=lambda x: x[1])
        for name, pv, nes in res[:top_n]:
            if pv < alpha:
                counts[name] = counts.get(name, 0) + 1
            nes_store.setdefault(name, []).append(nes)
    out = pd.DataFrame(
        {
            "term": list(counts),
            "selection_freq": [counts[k] / n_rep for k in counts],
            "NES_median": [float(np.nanmedian(nes_store[k])) for k in counts],
            "n_boot": n_rep,
        }
    )
    return out.sort_values("selection_freq", ascending=False).reset_index(drop=True) if not out.empty else out


def integrated_enrichment(
    ora_df: pd.DataFrame,
    gsea_df: pd.DataFrame,
    stability_df: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Combine ORA, GSEA and (optionally) bootstrap stability per term.

    The integrated evidence score is a weighted mean of rank-normalised
    ``-log10(padj)`` from ORA and GSEA plus the bootstrap selection frequency,
    so a term needs to be significant *and* reproducible to rank highly.
    """
    frames = []
    if ora_df is not None and not ora_df.empty:
        o = ora_df.copy()
        o["neglog10_padj_ora"] = -np.log10(np.clip(o["padj"], 1e-300, 1))
        # keep one row per (term, direction): up- and down-regulation are
        # different biological hypotheses and must not overwrite each other
        key = ["term", "direction"] if "direction" in o.columns else ["term"]
        o = o.sort_values("neglog10_padj_ora", ascending=False).drop_duplicates(key)
        cols = [c for c in ("term", "description", "direction", "overlap", "set_size", "expected", "fold_enrichment", "padj", "neglog10_padj_ora") if c in o.columns]
        frames.append(o[cols])
    if gsea_df is not None and not gsea_df.empty:
        g = gsea_df.copy()
        g["neglog10_padj_gsea"] = -np.log10(np.clip(g["padj"], 1e-300, 1))
        gcols = [c for c in ("term", "NES", "padj", "neglog10_padj_gsea", "direction") if c in g.columns]
        frames.append(g[gcols])
    if not frames:
        return pd.DataFrame()
    out = frames[0]
    for f in frames[1:]:
        f = f.copy()
        if "direction" in out.columns and "direction" in f.columns:
            # align the GSEA direction to the ORA direction so the join is 1:1
            out = out.merge(f.drop(columns=["direction"]), on="term", how="left")
        else:
            out = out.merge(f, on="term", how="left")
    if stability_df is not None and not stability_df.empty:
        out = out.merge(stability_df[["term", "selection_freq", "NES_median"]], on="term", how="left")
    out["selection_freq"] = out.get("selection_freq", pd.Series(0.0, index=out.index)).fillna(0.0)
    out["neglog10_padj_ora"] = out.get("neglog10_padj_ora", pd.Series(0.0, index=out.index)).fillna(0.0)
    out["neglog10_padj_gsea"] = out.get("neglog10_padj_gsea", pd.Series(0.0, index=out.index)).fillna(0.0)
    r_ora = out["neglog10_padj_ora"].rank(pct=True)
    r_gsea = out["neglog10_padj_gsea"].rank(pct=True)
    out["evidence"] = 0.35 * r_ora + 0.35 * r_gsea + 0.30 * out["selection_freq"]
    return out.sort_values("evidence", ascending=False).reset_index(drop=True)
