"""Download TCGA-LIHC STAR gene counts via the GDC API + Enrichr gene-set libraries.

Outputs
-------
data/raw/gdc/            per-sample *.rna_seq.augmented_star_gene_counts.tsv
data/raw/gdc_index.tsv   file_id -> case / sample_type mapping
data/raw/genesets/*.txt  Enrichr text libraries (KEGG, GO BP, Hallmark)
data/raw/geo/GSE14520*   validation cohort (GPL3921 series matrix)

The script is resumable: files already on disk with a matching size are skipped.
"""
from __future__ import annotations

import concurrent.futures as cf
import json
import os
import sys
import time
from pathlib import Path

import requests

# repository root: derived from this file so the scripts work from any checkout
ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
GDC_DIR = RAW / "gdc"
GS_DIR = RAW / "genesets"
GEO_DIR = RAW / "geo"
for d in (GDC_DIR, GS_DIR, GEO_DIR):
    d.mkdir(parents=True, exist_ok=True)

GDC_API = "https://api.gdc.cancer.gov"
UA = {"User-Agent": "Mozilla/5.0 (academic course project; python-requests)"}

ENRICHR_LIBS = {
    "KEGG_2021_Human": GS_DIR / "KEGG_2021_Human.txt",
    "GO_Biological_Process_2023": GS_DIR / "GO_Biological_Process_2023.txt",
    "MSigDB_Hallmark_2020": GS_DIR / "MSigDB_Hallmark_2020.txt",
    "Reactome_2022": GS_DIR / "Reactome_2022.txt",
}


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def query_files() -> list[dict]:
    flt = {
        "op": "and",
        "content": [
            {"op": "in", "content": {"field": "cases.project.project_id", "value": ["TCGA-LIHC"]}},
            {"op": "in", "content": {"field": "data_type", "value": ["Gene Expression Quantification"]}},
            {"op": "in", "content": {"field": "analysis.workflow_type", "value": ["STAR - Counts"]}},
        ],
    }
    params = {
        "filters": json.dumps(flt),
        "fields": "file_id,file_size,cases.submitter_id,cases.samples.sample_type,cases.samples.submitter_id",
        "size": 2000,
        "format": "JSON",
    }
    r = requests.get(f"{GDC_API}/files", params=params, headers=UA, timeout=120)
    r.raise_for_status()
    hits = r.json()["data"]["hits"]
    rows = []
    for h in hits:
        case = (h.get("cases") or [{}])[0]
        samples = case.get("samples") or [{}]
        rows.append(
            {
                "file_id": h["file_id"],
                "file_size": h["file_size"],
                "case_id": case.get("submitter_id", ""),
                "sample_type": samples[0].get("sample_type", ""),
                "sample_id": samples[0].get("submitter_id", ""),
            }
        )
    return rows


def write_index(rows: list[dict]) -> None:
    idx = RAW / "gdc_index.tsv"
    with open(idx, "w", encoding="utf-8", newline="") as fh:
        fh.write("file_id\tfile_size\tcase_id\tsample_type\tsample_id\n")
        for r in rows:
            fh.write(f"{r['file_id']}\t{r['file_size']}\t{r['case_id']}\t{r['sample_type']}\t{r['sample_id']}\n")
    log(f"index rows: {len(rows)} -> {idx}")


def download_one(row: dict, retries: int = 4) -> bool:
    dest = GDC_DIR / f"{row['file_id']}.tsv"
    if dest.exists() and dest.stat().st_size == row["file_size"]:
        return True
    for attempt in range(retries):
        try:
            with requests.get(f"{GDC_API}/data/{row['file_id']}", headers=UA, stream=True, timeout=180) as r:
                if r.status_code != 200:
                    raise RuntimeError(f"HTTP {r.status_code}")
                tmp = dest.with_suffix(".part")
                with open(tmp, "wb") as fh:
                    for chunk in r.iter_content(1 << 20):
                        fh.write(chunk)
                tmp.replace(dest)
            ok = dest.stat().st_size == row["file_size"]
            if ok:
                return True
        except Exception as exc:  # noqa: BLE001
            time.sleep(1.5 * (attempt + 1))
            if attempt == retries - 1:
                log(f"FAILED {row['file_id']}: {type(exc).__name__} {exc}")
    return False


def download_all(rows: list[dict], workers: int = 8) -> None:
    done = 0
    t0 = time.time()
    with cf.ThreadPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(download_one, r): r for r in rows}
        for fut in cf.as_completed(futs):
            done += 1
            if done % 20 == 0 or done == len(rows):
                el = time.time() - t0
                log(f"gdc files {done}/{len(rows)}  elapsed {el/60:.1f} min")


def download_genesets() -> None:
    for name, dest in ENRICHR_LIBS.items():
        if dest.exists() and dest.stat().st_size > 1000:
            log(f"gene set cached: {dest.name}")
            continue
        try:
            r = requests.get(
                "https://maayanlab.cloud/Enrichr/geneSetLibrary",
                params={"mode": "text", "libraryName": name},
                headers=UA,
                timeout=180,
            )
            r.raise_for_status()
            dest.write_text(r.text, encoding="utf-8")
            n = sum(1 for line in r.text.splitlines() if line.strip())
            log(f"gene set {name}: {n} sets, {len(r.text)/1024:.0f} KB")
        except Exception as exc:  # noqa: BLE001
            log(f"gene set FAILED {name}: {exc}")


def download_geo() -> None:
    base = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE14nnn/GSE14520"
    for url, name in (
        (f"{base}/matrix/GSE14520-GPL3921_series_matrix.txt.gz", "GSE14520-GPL3921_series_matrix.txt.gz"),
        (f"{base}/suppl/GSE14520_Extra_Supplement.txt.gz", "GSE14520_Extra_Supplement.txt.gz"),
        (f"{base}/suppl/GSE14520_README.txt", "GSE14520_README.txt"),
    ):
        dest = GEO_DIR / name
        if dest.exists() and dest.stat().st_size > 1000:
            log(f"geo cached: {name}")
            continue
        try:
            with requests.get(url, headers=UA, stream=True, timeout=300) as r:
                r.raise_for_status()
                with open(dest, "wb") as fh:
                    for chunk in r.iter_content(1 << 20):
                        fh.write(chunk)
            log(f"geo {name}: {dest.stat().st_size/1e6:.2f} MB")
        except Exception as exc:  # noqa: BLE001
            log(f"geo FAILED {name}: {exc}")


def download_survival() -> None:
    url = "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/survival%2FLIHC_survival.txt"
    dest = RAW / "LIHC_survival.txt"
    if dest.exists() and dest.stat().st_size > 1000:
        log("survival cached")
        return
    try:
        r = requests.get(url, headers=UA, timeout=120)
        r.raise_for_status()
        dest.write_text(r.text, encoding="utf-8")
        log(f"survival: {len(r.text.splitlines())} lines")
    except Exception as exc:  # noqa: BLE001
        log(f"survival FAILED: {exc}")


def main() -> int:
    step = sys.argv[1] if len(sys.argv) > 1 else "all"
    if step in ("all", "files"):
        rows = query_files()
        write_index(rows)
        log(f"files to fetch: {len(rows)}  total {sum(r['file_size'] for r in rows)/1e9:.2f} GB")
        download_all(rows)
    if step in ("all", "meta"):
        download_survival()
        download_genesets()
        download_geo()
    log("done")
    return 0


if __name__ == "__main__":
    sys.exit(main())
