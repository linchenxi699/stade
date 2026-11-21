"""Create a small, self-contained demo dataset so the pipeline runs without downloads.

The expression values are *simulated* (negative-binomial counts on a realistic
gene x sample layout) while gene symbols and the gene-set library are real.
Use it to check that the pipeline works end to end; for the published case study
run ``scripts/01_download_data.py`` instead.

Usage::

    python scripts/99_make_demo_dataset.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "demo"
OUT.mkdir(parents=True, exist_ok=True)
GSETS = ROOT / "data" / "raw" / "genesets"


def pick_genes(n_genes: int, seed: int = 0) -> list[str]:
    """Take a reproducible subset of real human gene symbols from the KEGG library."""
    pool: list[str] = []
    src = GSETS / "KEGG_2021_Human.txt"
    if src.exists():
        for line in src.read_text(encoding="utf-8", errors="replace").splitlines():
            parts = line.split("\t")
            pool.extend(g.strip().upper() for g in parts[2:] if g.strip())
    pool = sorted(set(pool))
    if len(pool) < n_genes:  # fall back to placeholder symbols
        pool += [f"GENE{i:05d}" for i in range(n_genes - len(pool))]
    rng = np.random.default_rng(seed)
    keep = rng.choice(len(pool), size=min(n_genes, len(pool)), replace=False)
    return [pool[i] for i in sorted(keep)]


def main() -> int:
    n_genes, n_case, n_ctrl = 1200, 30, 28
    genes = pick_genes(n_genes)
    rng = np.random.default_rng(20261215)

    lib = rng.uniform(6e6, 1.1e7, size=n_case + n_ctrl)
    base = rng.gamma(0.7, 90, size=len(genes))
    counts = np.outer(base, lib / 1e6)
    counts *= rng.gamma(6.0, 0.18, size=counts.shape)

    # planted biology: a proliferative programme up in tumours and a
    # hepatocyte-metabolism programme down, plus 40 sparse "noise" genes
    up_set, down_set = genes[:60], genes[60:120]
    up_fc = rng.uniform(1.6, 2.8, size=len(up_set))
    dn_fc = -rng.uniform(1.6, 2.8, size=len(down_set))
    noise = rng.choice(np.arange(200, len(genes)), size=40, replace=False)
    idx = {g: i for i, g in enumerate(genes)}
    counts[[idx[g] for g in up_set], :n_case] *= 2 ** up_fc[:, None]
    counts[[idx[g] for g in down_set], :n_case] *= 2 ** dn_fc[:, None]
    spk = rng.uniform(2, 4, size=40)
    for j, g in enumerate(noise):
        hot = rng.choice(n_case, size=3, replace=False)
        counts[idx[genes[g]], hot] *= 2 ** spk[j]

    samples = [f"DEMO-T{i:02d}" for i in range(n_case)] + [f"DEMO-N{i:02d}" for i in range(n_ctrl)]
    counts_df = pd.DataFrame(np.round(counts).astype(int), index=genes, columns=samples)
    counts_df.to_csv(OUT / "demo_counts.tsv.gz", sep="\t", compression="gzip")

    meta = pd.DataFrame(
        {
            "sample_id": samples,
            "group": ["Tumor"] * n_case + ["Normal"] * n_ctrl,
        }
    )
    meta.to_csv(OUT / "demo_sample_meta.tsv", sep="\t", index=False)
    print(f"demo dataset -> {OUT}")
    print(f"  {counts_df.shape[0]} genes x {counts_df.shape[1]} samples "
          f"({n_case} tumour / {n_ctrl} normal)")
    print(f"  planted: {len(up_set)} up, {len(down_set)} down, {len(noise)} sparse noise genes")
    print("  NOTE: expression values are simulated; gene symbols are real.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
