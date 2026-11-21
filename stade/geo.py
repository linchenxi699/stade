"""GEO series-matrix reader (plain text, no GEOquery).

Supports the ``GSE14520`` layout: ``!Sample_*`` header lines followed by a
``!series_matrix_table_begin`` block with probe x sample values.
"""
from __future__ import annotations

import gzip
import re
from pathlib import Path

import numpy as np
import pandas as pd

__all__ = ["read_series_matrix", "load_gse14520"]


def _split(line: str) -> list[str]:
    return [v.strip().strip('"') for v in line.rstrip("\n").split("\t")]


def read_series_matrix(path: str | Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (expression dataframe probe x sample, sample metadata dataframe)."""
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    meta: dict[str, list[str]] = {}
    rows: list[list[str]] = []
    in_table = False
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            if in_table:
                rows.append(_split(line))
                continue
            if line.startswith("!Sample_"):
                parts = _split(line)
                key = parts[0].lstrip("!")
                meta.setdefault(key, []).append(parts[1:])
    expr = pd.DataFrame(rows[1:], columns=rows[0]).set_index("ID_REF")
    expr = expr.apply(pd.to_numeric, errors="coerce")
    # expand multi-row sample characteristics into columns
    samples = meta.get("Sample_geo_accession", [[]])[0]
    sm = pd.DataFrame(index=pd.Index(samples, name="gsm"))
    for key, values in meta.items():
        if key in ("Sample_geo_accession",):
            continue
        if len(values) == 1 and len(values[0]) == len(samples):
            sm[key.replace("Sample_", "")] = values[0]
        elif len(values) > 1:
            # e.g. characteristics split across several lines: merge per column
            if all(len(v) == len(samples) for v in values):
                merged = []
                for j in range(len(samples)):
                    parts = [v[j] for v in values if v[j]]
                    merged.append("; ".join(parts))
                sm[key.replace("Sample_", "")] = merged
    sm["title"] = meta.get("Sample_title", [[""] * len(samples)])[0]
    return expr, sm


def load_gse14520(
    matrix_path: str | Path,
    supplement_path: str | Path | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Load GSE14520: returns (expression genes x samples, sample meta, survival)."""
    expr, sm = read_series_matrix(matrix_path)
    # GPL3921 probes are mapped to symbols in the platform annotation; the
    # GSE14520 matrix stores probe ids, so collapse duplicates by mean.
    expr.index.name = "probe"
    title = sm["title"].astype(str)
    sm["group"] = np.where(title.str.contains("Tumor", case=False, na=False), "Tumor", "Normal")
    sm.loc[title.str.contains("Non-Tumor", case=False, na=False), "group"] = "Normal"
    c = sm.get("characteristics_ch1", pd.Series("", index=sm.index)).astype(str)
    if "Tissue" in c.str.cat():
        sm["group"] = np.where(c.str.contains("Non-Tumor"), "Normal", np.where(c.str.contains("Tumor"), "Tumor", sm["group"]))
    sm["sample_id"] = sm.index
    expr = expr.loc[:, [s for s in expr.columns if s in sm.index]]

    surv = pd.DataFrame(index=sm.index)
    if supplement_path is not None and Path(supplement_path).exists():
        sup = pd.read_csv(supplement_path, sep="\t", dtype=str)
        sup = sup.dropna(subset=["Affy_GSM"]).drop_duplicates("Affy_GSM").set_index("Affy_GSM")
        surv["OS"] = pd.to_numeric(sup["Survival status"].reindex(surv.index), errors="coerce")
        surv["OS.time"] = pd.to_numeric(sup["Survival months"].reindex(surv.index), errors="coerce") * 30.44
        surv["recurrence"] = pd.to_numeric(sup["Recurr status"].reindex(surv.index), errors="coerce")
        surv["stage"] = sup["TNM staging"].reindex(surv.index)
        surv["gender"] = sup["Gender"].reindex(surv.index)
        surv["age"] = pd.to_numeric(sup["Age"].reindex(surv.index), errors="coerce")
        surv["lcs_id"] = sup["LCS ID"].reindex(surv.index)
    return expr, sm, surv
