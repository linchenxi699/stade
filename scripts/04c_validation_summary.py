"""TCGA-LIHC <-> GSE14520 cross-cohort validation summary (numbers for the paper)."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

# repository root: derived from this file so the scripts work from any checkout
ROOT = Path(__file__).resolve().parents[1]
VAL = ROOT / "results" / "validation"
sys.path.insert(0, str(ROOT))


def main() -> int:
    cmp = pd.read_csv(VAL / "concordance_pairs.tsv", sep="\t", index_col=0)
    out: dict = {"n_common": int(len(cmp))}
    out["spearman_all"] = float(stats.spearmanr(cmp["tcga_logFC"], cmp["geo_logFC"]).statistic)
    sig = cmp[cmp["tcga_sig"]]
    out["n_tcga_sig"] = int(len(sig))
    out["spearman_sig"] = float(stats.spearmanr(sig["tcga_logFC"], sig["geo_logFC"]).statistic)
    rob = cmp[cmp["tcga_robust"]]
    out["n_tcga_robust"] = int(len(rob))
    out["sign_agreement_robust"] = float((np.sign(rob["tcga_logFC"]) == np.sign(rob["geo_logFC"])).mean())
    out["geo_sig_rate_robust"] = float(rob["geo_sig"].mean())
    out["geo_sig_rate_all"] = float(cmp["geo_sig"].mean())
    out["pearson_robust"] = float(rob["tcga_logFC"].corr(rob["geo_logFC"]))
    # replication rates by TCGA stability decile
    if "tcga_stability" in cmp:
        pass
    pv = VAL / "prognostic_validation.json"
    if pv.exists():
        out["prognostic"] = json.loads(pv.read_text(encoding="utf-8"))
    (VAL / "validation_summary.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
