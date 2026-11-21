"""Assemble the GDC STAR-count files into an expression matrix + metadata.

The GDC augmented STAR files contain a ``gene_id`` / ``gene_name`` /
``gene_type`` header block followed by 9 numeric columns; ``unstranded`` is the
raw read count used by ``HTSeq``-style DE analysis.

Outputs (in ``data/processed``)
-------------------------------
counts.tsv.gz        gene x sample raw counts (Ensembl id, no version suffix)
gene_meta.tsv        gene_id, gene_name, gene_type, avg_count, detected_in
sample_meta.tsv      file_id, case_id, sample_id, sample_type, group
clinical.tsv         OS.time (days), OS (0/1), plus GDC clinical fields
"""
from __future__ import annotations

import concurrent.futures as cf
import io
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import requests

# repository root: derived from this file so the scripts work from any checkout
ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
OUT = ROOT / "data" / "processed"
OUT.mkdir(parents=True, exist_ok=True)
GDC_API = "https://api.gdc.cancer.gov"
UA = {"User-Agent": "Mozilla/5.0 (academic course project; python-requests)"}


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def parse_one(path: Path) -> tuple[str, np.ndarray, list[str], list[str]] | None:
    """Return (gene_id series, counts vector, gene_name list, gene_type list)."""
    try:
        df = pd.read_csv(path, sep="\t", comment="#", dtype=str, low_memory=False)
    except Exception as exc:  # noqa: BLE001
        log(f"parse FAILED {path.name}: {exc}")
        return None
    if "unstranded" not in df.columns or "gene_id" not in df.columns:
        return None
    n_top = df["gene_id"].str.startswith("N_")
    df = df.loc[~n_top].copy()
    gid = df["gene_id"].str.split(".").str[0]
    counts = pd.to_numeric(df["unstranded"], errors="coerce").fillna(0).to_numpy(dtype=np.int64)
    return gid.to_numpy(), counts, df["gene_name"].to_numpy(), df["gene_type"].to_numpy()


def build_matrix() -> None:
    index = pd.read_csv(RAW / "gdc_index.tsv", sep="\t")
    index = index[index["sample_type"].isin(["Primary Tumor", "Solid Tissue Normal"])].copy()
    index["group"] = np.where(index["sample_type"] == "Primary Tumor", "Tumor", "Normal")
    files = [(row.file_id, RAW / "gdc" / f"{row.file_id}.tsv") for row in index.itertuples()]
    files = [(fid, p) for fid, p in files if p.exists()]
    log(f"assembling {len(files)} samples")
    results: dict[str, tuple[np.ndarray, list[str], list[str]]] = {}
    gene_order: list[str] | None = None
    with cf.ThreadPoolExecutor(max_workers=8) as ex:
        futs = {ex.submit(parse_one, p): fid for fid, p in files}
        for fut in cf.as_completed(futs):
            fid = futs[fut]
            r = fut.result()
            if r is None:
                continue
            gid, counts, names, types = r
            if gene_order is None:
                gene_order = list(gid)
                gene_names, gene_types = names, types
            elif list(gid) != gene_order:
                # align by id, filling missing genes with 0
                pos = {g: i for i, g in enumerate(gid)}
                counts = np.array([counts[pos[g]] if g in pos else 0 for g in gene_order], dtype=np.int64)
            results[fid] = (counts, names, types)
    if not results or gene_order is None:
        raise RuntimeError("no GDC files parsed")
    sample_ids = list(results)
    mat = np.column_stack([results[s][0] for s in sample_ids])
    counts_df = pd.DataFrame(mat, index=pd.Index(gene_order, name="gene_id"), columns=sample_ids)
    log(f"raw matrix: {counts_df.shape}")
    counts_df.to_csv(OUT / "counts.tsv.gz", sep="\t", compression="gzip")

    gene_meta = pd.DataFrame({"gene_id": gene_order, "gene_name": gene_names, "gene_type": gene_types})
    gene_meta["avg_count"] = counts_df.mean(axis=1).to_numpy()
    gene_meta["detected_in"] = (counts_df > 0).sum(axis=1).to_numpy()
    # keep a unique symbol column: prefer the canonical symbol when duplicated
    gene_meta["gene_symbol"] = gene_meta["gene_name"].where(gene_meta["gene_name"].astype(bool), gene_meta["gene_id"])
    gene_meta.to_csv(OUT / "gene_meta.tsv", sep="\t", index=False)

    sm = index[index["file_id"].isin(sample_ids)][["file_id", "case_id", "sample_id", "sample_type", "group"]]
    sm.to_csv(OUT / "sample_meta.tsv", sep="\t", index=False)
    log(f"sample meta: {len(sm)}  tumor={int((sm.group=='Tumor').sum())} normal={int((sm.group=='Normal').sum())}")


def fetch_clinical() -> None:
    """Clinical + follow-up from the GDC API (one row per case)."""
    fields = [
        "submitter_id",
        "project.project_id",
        "demographic.gender",
        "demographic.age_at_index",
        "demographic.race",
        "diagnoses.age_at_diagnosis",
        "diagnoses.ajcc_pathologic_stage",
        "diagnoses.ajcc_pathologic_t",
        "diagnoses.ajcc_pathologic_n",
        "diagnoses.ajcc_pathologic_m",
        "diagnoses.primary_diagnosis",
        "diagnoses.tumor_grade",
        "diagnoses.vascular_invasion",
        "diagnoses.child_pugh_classification",
        "diagnoses.afp_level",
        "diagnoses.days_to_death",
        "diagnoses.days_to_last_follow_up",
        "diagnoses.vital_status",
    ]
    params = {
        "filters": json.dumps({"op": "in", "content": {"field": "project.project_id", "value": ["TCGA-LIHC"]}}),
        "fields": ",".join(fields),
        "size": 1000,
        "format": "JSON",
        "expand": "diagnoses,demographic",
    }
    r = requests.get(f"{GDC_API}/cases", params=params, headers=UA, timeout=180)
    r.raise_for_status()
    hits = r.json()["data"]["hits"]
    rows = []
    for h in hits:
        dem = h.get("demographic") or {}
        diags = h.get("diagnoses") or [{}]
        d = diags[0] if diags else {}
        rows.append(
            {
                "case_id": h.get("submitter_id"),
                "gender": dem.get("gender"),
                "age_at_index": dem.get("age_at_index"),
                "race": dem.get("race"),
                "age_at_diagnosis_days": d.get("age_at_diagnosis"),
                "stage": d.get("ajcc_pathologic_stage"),
                "t_stage": d.get("ajcc_pathologic_t"),
                "n_stage": d.get("ajcc_pathologic_n"),
                "m_stage": d.get("ajcc_pathologic_m"),
                "grade": d.get("tumor_grade"),
                "vascular_invasion": d.get("vascular_invasion"),
                "child_pugh": d.get("child_pugh_classification"),
                "afp": d.get("afp_level"),
                "days_to_death": d.get("days_to_death"),
                "days_to_last_follow_up": d.get("days_to_last_follow_up"),
                "vital_status": d.get("vital_status"),
            }
        )
    clin = pd.DataFrame(rows).drop_duplicates("case_id").set_index("case_id")
    log(f"clinical cases: {len(clin)}")
    clin.to_csv(OUT / "clinical_gdc.tsv", sep="\t")

    # ---- merge the curated Xena survival table (OS / DSS / DFI / PFI) ----
    xena_p = RAW / "LIHC_survival.txt"
    if xena_p.exists():
        x = pd.read_csv(xena_p, sep="\t")
        x = x.drop_duplicates("_PATIENT").set_index("_PATIENT")
        for col in ("OS", "OS.time", "DSS", "DSS.time", "PFI", "PFI.time"):
            if col in x.columns:
                clin[col] = x[col].reindex(clin.index)
        log(f"xena survival merged: {clin['OS.time'].notna().sum()} cases with OS time")

    surv = pd.DataFrame(index=clin.index)
    if "OS" in clin.columns and clin["OS"].notna().any():
        surv["OS"] = pd.to_numeric(clin["OS"], errors="coerce")
        surv["OS.time"] = pd.to_numeric(clin["OS.time"], errors="coerce")
    else:
        surv["OS"] = (clin["vital_status"].str.lower() == "dead").astype(float)
        surv["OS.time"] = np.where(
            clin["vital_status"].str.lower() == "dead",
            pd.to_numeric(clin["days_to_death"], errors="coerce"),
            pd.to_numeric(clin["days_to_last_follow_up"], errors="coerce"),
        )
    surv["vital_status"] = clin["vital_status"]
    surv["stage"] = clin["stage"]
    surv["age"] = pd.to_numeric(clin["age_at_index"], errors="coerce")
    surv["gender"] = clin["gender"]
    surv = surv.dropna(subset=["OS.time", "OS"])
    surv = surv[surv["OS.time"] > 0]
    surv.to_csv(OUT / "clinical.tsv", sep="\t")
    log(f"clinical with OS: {len(surv)}  events={int(surv['OS'].sum())}")


def check_survival_file() -> None:
    p = RAW / "LIHC_survival.txt"
    if p.exists():
        head = p.read_text(encoding="utf-8", errors="replace").splitlines()[:3]
        log(f"xena survival file: {len(head)} preview lines -> {head[0][:120]}")


def main() -> int:
    step = sys.argv[1] if len(sys.argv) > 1 else "all"
    if step in ("all", "matrix"):
        t0 = time.time()
        build_matrix()
        log(f"matrix done in {(time.time()-t0)/60:.1f} min")
    if step in ("all", "clinical"):
        fetch_clinical()
        check_survival_file()
    return 0


if __name__ == "__main__":
    sys.exit(main())
