"""Command line interface for StaDE.

Examples
--------
::

    stade prep   --counts data/processed/counts.tsv.gz --meta data/processed/sample_meta.tsv \\
                 --group-col group --case-col Tumor --outdir work
    stade de     --bundle work/bundle.npz --outdir results --n-boot 200
    stade enrich --de results/de_consensus.tsv --genesets data/raw/genesets/KEGG_2021_Human.txt \\
                 --library KEGG --outdir results/enrichment
    stade surv   --expr work/logcpm.tsv.gz --clinical data/processed/clinical.tsv \\
                 --genes results/candidates.txt --outdir results/survival
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

from . import __version__
from .enrichment import (
    bootstrap_enrichment,
    gsea,
    integrated_enrichment,
    load_enrichr_tsv,
    load_gmt,
    ora,
    ora_weighted,
)
from .io import Provenance, read_matrix, write_matrix
from .normalization import filter_low_expression, log_cpm
from .plotting import (
    plot_bootstrap_fdr,
    plot_enrichment_dot,
    plot_ma,
    plot_method_concordance,
    plot_pca,
    plot_stability_scatter,
    plot_volcano,
    setup_style,
)
from .stability import bootstrap_de, consensus_table, empirical_null, gene_stability
from .stats import fit_ebayes, significance, welch_test, wilcoxon_test
from .survival import cox_multivariate, cox_univariate, km_analysis, prepare_survival, candidate_signature

LOGFMT = "[%(asctime)s] %(message)s"


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def cmd_prep(args: argparse.Namespace) -> int:
    prov = Provenance(command="prep", parameters=vars(args))
    prov.add_input("counts", args.counts)
    counts = read_matrix(args.counts)
    log(f"counts: {counts.shape[0]} genes x {counts.shape[1]} samples")
    meta = read_matrix(args.meta) if args.meta else None
    if meta is not None and args.meta_index and args.meta_index in meta.columns:
        meta = meta.set_index(args.meta_index)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    groups = None
    if meta is not None and args.group_col:
        groups = meta.loc[counts.columns, args.group_col]
        if args.case_col:
            groups = pd.Series(np.where(groups == args.case_col, args.case_col, "Other"), index=groups.index, name="group")
        counts = counts.loc[:, groups.notna()]
        groups = groups.loc[counts.columns]
        meta = meta.loc[counts.columns]
        log(f"groups: {groups.value_counts().to_dict()}")

    kept = filter_low_expression(counts, min_count=args.min_count, min_prop=args.min_prop, groups=groups)
    log(f"after expression filter: {kept.shape[0]} genes")
    expr = log_cpm(kept, prior_count=args.prior_count, use_tmm=args.tmm)
    write_matrix(expr, outdir / "logcpm.tsv.gz")
    write_matrix(kept, outdir / "filtered_counts.tsv.gz")
    if groups is not None:
        groups.to_frame().to_csv(outdir / "groups.tsv", sep="\t")
    bundle = {
        "genes": np.array(expr.index, dtype=object),
        "samples": np.array(expr.columns, dtype=object),
        "expr": expr.to_numpy(dtype=np.float64),
    }
    np.savez_compressed(outdir / "bundle.npz", **bundle)
    prov.finish(n_genes=int(expr.shape[0]), n_samples=int(expr.shape[1]))
    prov.write(outdir / "provenance_prep.json")
    log(f"wrote {outdir/'bundle.npz'}")
    return 0


def _load_bundle(path: Path) -> tuple[pd.DataFrame, pd.Series | None]:
    d = np.load(path, allow_pickle=True)
    expr = pd.DataFrame(d["expr"], index=d["genes"], columns=d["samples"])
    gpath = path.parent / "groups.tsv"
    groups = pd.read_csv(gpath, sep="\t", index_col=0).iloc[:, 0] if gpath.exists() else None
    return expr, groups


def cmd_de(args: argparse.Namespace) -> int:
    prov = Provenance(command="de", parameters=vars(args))
    expr, groups = _load_bundle(Path(args.bundle))
    if groups is None:
        groups = pd.Series(np.where(np.arange(expr.shape[1]) % 2 == 0, args.case_col, "Other"), index=expr.columns)
        log("WARNING: no groups.tsv found; using the first half as cases")
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    design = pd.DataFrame(
        {"Intercept": 1.0, "group": (groups == args.case_col).astype(float)}, index=expr.columns
    )
    res = fit_ebayes(expr, design, "group", trend=not args.no_trend, robust=not args.no_robust)
    de = significance(res.table, alpha=args.alpha, min_abs_logfc=args.min_abs_logfc)
    de["gene"] = de.index
    log(f"EB d0={res.d0:.2f} df_total={res.df_total:.1f}; significant={int(de['significant'].sum())}")
    de.to_csv(outdir / "de_ebayes.tsv", sep="\t")

    # bootstrap stability
    n_boot = args.n_boot
    boot = bootstrap_de(
        expr, groups, args.case_col, n_boot=n_boot, seed=args.seed,
        method="ebayes" if not args.fast_boot else "welch",
    )
    stab = gene_stability(boot, de["logFC"], alpha=args.alpha, min_abs_logfc=args.min_abs_logfc)
    de = de.join(stab[["sel_freq", "sign_consistency", "effect_cv", "confidence", "stability"]])
    de["robust"] = de["significant"] & (de["stability"] >= args.min_stability)
    de.to_csv(outdir / "de_stability.tsv", sep="\t")
    if boot.tstat is not None:
        boot.tstat.to_csv(outdir / "bootstrap_t.tsv.gz", sep="\t", compression="gzip")
    log(
        f"bootstrap: {n_boot} replicates; significant & stable>={args.min_stability}: "
        f"{int(de['robust'].sum())} (median stability among significant = "
        f"{de.loc[de['significant'],'stability'].median():.3f})"
    )

    # competing methods + consensus
    wel = welch_test(expr, groups, args.case_col)
    wil = wilcoxon_test(expr, groups, args.case_col)
    cons = consensus_table(res.table, wel, wil, alpha=args.alpha).join(stab[["stability", "sel_freq"]])
    cons.to_csv(outdir / "de_consensus.tsv", sep="\t")
    log(f"methods agree on {(cons['n_methods_sig']==3).sum()} genes; EB only {(cons['n_methods_sig']==1).sum()}")

    # empirical null / local fdr
    en = empirical_null(res.table["t"].to_numpy())
    de = de.join(en[["lfdr", "pi0"]].set_axis(de.index))
    de.to_csv(outdir / "de_stability.tsv", sep="\t")
    null_sd = float(np.nanmedian(en["null_sd"]))
    log(f"empirical null sd={null_sd:.3f}, pi0={float(np.nanmedian(en['pi0'])):.3f}")

    setup_style()
    figs = outdir / "figs"
    plot_pca(expr, groups, figs / "pca.png", title=f"PCA ({args.case_col} vs rest)")
    plot_ma(de, figs / "ma.png", alpha=args.alpha, min_abs_logfc=args.min_abs_logfc)
    plot_volcano(de, figs / "volcano.png", alpha=args.alpha, min_abs_logfc=args.min_abs_logfc)
    plot_stability_scatter(de, stab, figs / "stability.png")
    plot_method_concordance(cons, figs / "concordance.png")
    plot_bootstrap_fdr(stab, de, figs / "bootstrap_fdr.png")
    log(f"figures -> {figs}")

    summary = {
        "d0": res.d0,
        "df_total": res.df_total,
        "n_tested": int(de.shape[0]),
        "n_sig": int(de["significant"].sum()),
        "n_up": int(((de["significant"]) & (de["logFC"] > 0)).sum()),
        "n_down": int(((de["significant"]) & (de["logFC"] < 0)).sum()),
        "n_robust": int(de["robust"].sum()),
        "n_robust_up": int(((de["robust"]) & (de["logFC"] > 0)).sum()),
        "n_robust_down": int(((de["robust"]) & (de["logFC"] < 0)).sum()),
        "median_stability_significant": float(de.loc[de["significant"], "stability"].median()),
        "median_stability_robust": float(de.loc[de["robust"], "stability"].median()),
        "methods_all_three": int((cons["n_methods_sig"] == 3).sum()),
        "null_sd": null_sd,
        "pi0": float(np.nanmedian(en["pi0"])),
        "top_genes": de.loc[de["robust"]].sort_values(["stability", "t"], ascending=[False, False]).head(20)["gene"].tolist(),
        "top_genes_by_t": de.loc[de["significant"]].sort_values("t", ascending=False).head(20)["gene"].tolist(),
    }
    prov.finish(**summary)
    prov.write(outdir / "provenance_de.json")
    (outdir / "de_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0


def cmd_enrich(args: argparse.Namespace) -> int:
    prov = Provenance(command="enrich", parameters=vars(args))
    de = pd.read_csv(args.de, sep="\t", index_col=0)
    de.index = de.index.astype(str)
    # map Ensembl ids -> symbols when requested (keeps a symbol column for output)
    if args.symbol_map:
        smap = pd.read_csv(args.symbol_map, sep="\t")
        if args.symbol_id_col not in smap.columns:
            smap = smap.reset_index().rename(columns={smap.index.name or "index": args.symbol_id_col})
        id_col = args.symbol_id_col
        sy_col = args.symbol_col
        if id_col in smap.columns and sy_col in smap.columns:
            smap = smap.drop_duplicates(id_col).set_index(id_col)
            de["symbol"] = smap[sy_col].reindex(de.index).astype(object)
            # prefer symbols for enrichment pre-filtering
            if args.use_symbols:
                keep = de["symbol"].notna() & (de["symbol"].astype(str) != "")
                de = de.loc[keep]
                de.index = de.loc[keep, "symbol"].astype(str).to_numpy()
                de = de[~de.index.duplicated(keep="first")]
            log(f"symbol mapping applied: {len(de)} genes with symbols")
    if args.universe_file:
        universe = set(pd.read_csv(args.universe_file, sep="\t", index_col=0).index.astype(str))
    else:
        universe = set(de.index.astype(str))
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    alpha = args.alpha

    flag = "robust" if args.robust_only and "robust" in de.columns else "significant"
    sig = de[de[flag]]
    if "logFC" in sig:
        up = set(sig.index[sig["logFC"] > 0].astype(str))
        dn = set(sig.index[sig["logFC"] < 0].astype(str))
    else:
        up, dn = set(sig.index.astype(str)), set()
    log(f"query ({flag}): {len(up)} up, {len(dn)} down, universe {len(universe)}")

    gs_path = Path(args.genesets)
    gs = load_gmt(gs_path) if gs_path.suffix.lower() == ".gmt" else load_enrichr_tsv(gs_path)
    gs = gs.restrict_universe(universe).filter_sizes(args.min_size, args.max_size)
    log(f"library {gs.source}: {len(gs)} sets after filtering")

    ora_up = ora(up, gs, universe, "Up", args.min_size, args.max_size)
    ora_dn = ora(dn, gs, universe, "Down", args.min_size, args.max_size) if dn else pd.DataFrame()
    o = pd.concat([ora_up, ora_dn], ignore_index=True) if not ora_dn.empty else ora_up
    if not o.empty:
        # a term can appear once per direction: keep the more significant one
        o = o.sort_values(["pvalue", "term"]).drop_duplicates("term").sort_values(["pvalue", "term"]).reset_index(drop=True)
    o.to_csv(outdir / f"ora_{args.library}.tsv", sep="\t", index=False)
    log(f"ORA significant (padj<{alpha}): {int((o['padj']<alpha).sum()) if not o.empty else 0}")

    ranking = de["t"] if "t" in de.columns else -np.log10(de["adj.P.Val"].clip(lower=1e-300)) * np.sign(de["logFC"])
    gene_w = de["stability"] if args.stability_weight and "stability" in de.columns else None
    g = gsea(ranking=ranking, gene_sets=gs, n_perm=args.n_perm, min_size=args.min_size, max_size=args.max_size, gene_weights=gene_w)
    g.to_csv(outdir / f"gsea_{args.library}.tsv", sep="\t", index=False)
    log(f"GSEA significant (padj<{alpha}): {int((g['padj']<alpha).sum()) if not g.empty else 0}")

    boot = pd.DataFrame()
    n_boot_rankings = 0
    if args.n_boot > 0:
        if args.boot_t:
            bt = pd.read_csv(args.boot_t, sep="\t", index_col=0)
            bt.index = bt.index.astype(str)
            missing = bt.index.difference(de.index)
            if len(missing):
                log(f"WARNING: {len(missing)} bootstrap genes absent from DE table; dropping them")
                bt = bt.loc[bt.index.intersection(de.index)]
            n_boot_rankings = int(bt.shape[1])
            if args.boot_reps > 0:
                bt = bt.iloc[:, : args.boot_reps]
        else:
            raise SystemExit("--n-boot requires --boot-t (per-replicate t-statistics from `stade de`)")
        boot = bootstrap_enrichment(
            bt, gs, n_perm=args.boot_perm, min_size=args.min_size, max_size=args.max_size,
            alpha=args.alpha,
        )
        boot.to_csv(outdir / f"bootstrap_{args.library}.tsv", sep="\t", index=False)
        log(f"bootstrap enrichment over {n_boot_rankings} replicates: {len(boot)} stable terms")

    integ = integrated_enrichment(o, g, boot)
    integ.to_csv(outdir / f"integrated_{args.library}.tsv", sep="\t", index=False)
    top = integ.head(min(15, len(integ)))
    setup_style()
    pcol = "padj_x" if "padj_x" in top.columns else ("padj" if "padj" in top.columns else None)
    if "overlap" in top.columns and "expected" in top.columns:
        top = top.assign(fold_enrichment=top["overlap"] / top["expected"].replace(0, np.nan))
    plot_enrichment_dot(
        top, outdir / "figs" / f"enrichment_{args.library}.png",
        top_n=15, score=pcol or "padj", size="overlap" if "overlap" in top else "set_size",
        x="fold_enrichment" if "fold_enrichment" in top else None, color_by=pcol,
    )
    summary = {
        "library": args.library,
        "n_sets": len(gs),
        "ora_significant": int((o["padj"] < alpha).sum()) if not o.empty else 0,
        "gsea_significant": int((g["padj"] < alpha).sum()) if not g.empty else 0,
        "top_terms": top["term"].head(10).tolist(),
    }
    prov.finish(**summary)
    prov.write(outdir / "provenance_enrich.json")
    print(json.dumps(summary, indent=2))
    return 0


def cmd_surv(args: argparse.Namespace) -> int:
    prov = Provenance(command="surv", parameters=vars(args))
    expr = read_matrix(args.expr) if args.expr else _load_bundle(Path(args.bundle))[0]
    clin = read_matrix(args.clinical)
    clin, expr = prepare_survival(clin, expr, time_col=args.time_col, event_col=args.event_col)
    log(f"survival cohort: {expr.shape[1]} samples, events={int(clin[args.event_col].sum())}")
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    genes = None
    if args.genes:
        genes = [g.strip() for g in Path(args.genes).read_text().splitlines() if g.strip()]
    cox = cox_univariate(clin, expr, genes)
    cox.to_csv(outdir / "cox_univariate.tsv", sep="\t", index=False)
    log(f"Cox tested {len(cox)} genes; padj<0.05: {int((cox['padj']<0.05).sum()) if not cox.empty else 0}")

    top = cox.sort_values("pvalue").head(args.top_n)["gene"].tolist() if not cox.empty else []
    km_rows = []
    for g in top[: args.n_km]:
        res = km_analysis(
            clin[args.time_col], clin[args.event_col],
            expr.loc[g], cut=args.cut,
        )
        km_rows.append({"gene": g, "logrank_p": res.logrank_p, "median_high": res.median_high, "median_low": res.median_low, "HR": res.hr})
    pd.DataFrame(km_rows).to_csv(outdir / "km_summary.tsv", sep="\t", index=False)
    if args.sig and len(top) >= 2:
        panel = [g for g in top[: min(args.sig, 8)] if g in expr.index]
        try:
            summary, _ = cox_multivariate(clin, expr, panel, time_col=args.time_col, event_col=args.event_col)
            summary.to_csv(outdir / "cox_multivariate.tsv", sep="\t")
            log(f"multivariate model on {len(panel)} genes fitted")
        except Exception as exc:  # noqa: BLE001
            log(f"multivariate Cox failed: {exc}")
    prov.finish(n_genes_tested=int(len(cox)), n_cohort=int(expr.shape[1]), events=int(clin[args.event_col].sum()))
    prov.write(outdir / "provenance_surv.json")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="stade", description=f"StaDE {__version__}: stability-aware differential expression and enrichment")
    p.add_argument("--version", action="version", version=f"stade {__version__}")
    sub = p.add_subparsers(dest="command", required=True)

    # prep ------------------------------------------------------------------
    s = sub.add_parser("prep", help="load counts, filter genes, normalise, cache a bundle")
    s.add_argument("--counts", required=True)
    s.add_argument("--meta")
    s.add_argument("--meta-index", default=None)
    s.add_argument("--group-col", default=None)
    s.add_argument("--case-col", default=None)
    s.add_argument("--min-count", type=float, default=10.0)
    s.add_argument("--min-prop", type=float, default=0.25)
    s.add_argument("--prior-count", type=float, default=0.5)
    s.add_argument("--tmm", action="store_true")
    s.add_argument("--outdir", default="work")
    s.set_defaults(func=cmd_prep)

    # de --------------------------------------------------------------------
    s = sub.add_parser("de", help="empirical-Bayes DE + bootstrap stability + concordance")
    s.add_argument("--bundle", required=True)
    s.add_argument("--case-col", default="Tumor")
    s.add_argument("--alpha", type=float, default=0.05)
    s.add_argument("--min-abs-logfc", type=float, default=1.0)
    s.add_argument("--n-boot", type=int, default=200)
    s.add_argument("--seed", type=int, default=20261215)
    s.add_argument("--min-stability", type=float, default=0.7, help="stability threshold for a 'robust' call")
    s.add_argument("--fast-boot", action="store_true")
    s.add_argument("--no-trend", action="store_true")
    s.add_argument("--no-robust", action="store_true")
    s.add_argument("--outdir", default="results")
    s.set_defaults(func=cmd_de)

    # enrich ----------------------------------------------------------------
    s = sub.add_parser("enrich", help="ORA + GSEA + stability-aware integration")
    s.add_argument("--de", required=True)
    s.add_argument("--genesets", required=True)
    s.add_argument("--library", default="library")
    s.add_argument("--alpha", type=float, default=0.05)
    s.add_argument("--min-size", type=int, default=10)
    s.add_argument("--max-size", type=int, default=500)
    s.add_argument("--n-perm", type=int, default=1000)
    s.add_argument("--n-boot", type=int, default=0, help="bootstrap ranking replicates (0 = skip)")
    s.add_argument("--boot-t", default=None, help="genes x replicates t-statistic table from `stade de`")
    s.add_argument("--boot-reps", type=int, default=0, help="limit bootstrap replicates (0 = all)")
    s.add_argument("--boot-perm", type=int, default=300)
    s.add_argument("--stability-weight", action="store_true")
    s.add_argument("--robust-only", action="store_true", help="use only significant & stable genes as ORA query")
    s.add_argument("--symbol-map", default=None, help="TSV mapping gene_id -> gene symbol")
    s.add_argument("--symbol-id-col", default="gene_id")
    s.add_argument("--symbol-col", default="gene_symbol")
    s.add_argument("--use-symbols", action="store_true", help="replace Ensembl ids by symbols for enrichment")
    s.add_argument("--universe-file", default=None)
    s.add_argument("--outdir", default="results/enrichment")
    s.set_defaults(func=cmd_enrich)

    # surv ------------------------------------------------------------------
    s = sub.add_parser("surv", help="univariate/multivariate Cox and Kaplan-Meier")
    s.add_argument("--bundle")
    s.add_argument("--expr")
    s.add_argument("--clinical", required=True)
    s.add_argument("--genes")
    s.add_argument("--time-col", default="OS.time")
    s.add_argument("--event-col", default="OS")
    s.add_argument("--cut", default="median", choices=["median", "quartile", "optimal"])
    s.add_argument("--top-n", type=int, default=20)
    s.add_argument("--n-km", type=int, default=4)
    s.add_argument("--sig", type=int, default=0, help="fit a multivariate model on the top-N genes")
    s.add_argument("--outdir", default="results/survival")
    s.set_defaults(func=cmd_surv)
    return p


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return int(args.func(args) or 0)


if __name__ == "__main__":
    sys.exit(main())
