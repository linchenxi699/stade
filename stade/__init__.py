"""StaDE - stability-aware differential expression & enrichment for public transcriptomes."""

from __future__ import annotations

__version__ = "0.1.0"

from .enrichment import (  # noqa: F401
    GeneSetCollection,
    bootstrap_enrichment,
    gsea,
    integrated_enrichment,
    load_enrichr_tsv,
    load_gmt,
    ora,
    ora_weighted,
)
from .io import ExpressionDataset, Provenance  # noqa: F401
from .normalization import filter_low_expression, log_cpm, tmm_factors  # noqa: F401
from .stability import bootstrap_de, consensus_table, empirical_null, gene_stability  # noqa: F401
from .stats import fit_ebayes, significance, welch_test, wilcoxon_test  # noqa: F401

__all__ = [
    "__version__",
    "ExpressionDataset",
    "Provenance",
    "filter_low_expression",
    "log_cpm",
    "tmm_factors",
    "fit_ebayes",
    "significance",
    "welch_test",
    "wilcoxon_test",
    "bootstrap_de",
    "gene_stability",
    "consensus_table",
    "empirical_null",
    "GeneSetCollection",
    "load_enrichr_tsv",
    "load_gmt",
    "ora",
    "ora_weighted",
    "gsea",
    "bootstrap_enrichment",
    "integrated_enrichment",
]
