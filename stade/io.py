"""Data loading, metadata handling and provenance recording.

Design goals
------------
* one canonical long/wide representation (:class:`ExpressionDataset`);
* explicit, hashable provenance for every artefact (:class:`Provenance`);
* zero silent coercion: missing values and duplicate identifiers are reported.
"""
from __future__ import annotations

import gzip
import hashlib
import json
import platform
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

__all__ = [
    "ExpressionDataset",
    "Provenance",
    "read_matrix",
    "write_matrix",
    "file_sha256",
]


def file_sha256(path: str | Path, chunk: int = 1 << 20) -> str:
    """Streaming SHA-256 of a file (used for data provenance fingerprints)."""
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while block := fh.read(chunk):
            h.update(block)
    return h.hexdigest()


def _open_text(path: str | Path):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")


def read_matrix(path: str | Path, index_col: int | None = 0) -> pd.DataFrame:
    """Read a delimited matrix (``.tsv``/``.csv``, optionally gzip-compressed).

    ``index_col=None`` keeps a default RangeIndex, which is handy when the file
    has no meaningful row identifier yet (e.g. GDC sample metadata).
    """
    path = Path(path)
    sep = "\t" if path.name.endswith((".tsv", ".tsv.gz", ".txt", ".txt.gz")) else ","
    return pd.read_csv(path, sep=sep, index_col=index_col, low_memory=False)


def write_matrix(df: pd.DataFrame, path: str | Path, sep: str = "\t", float_format: str = "%.6g") -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep=sep, float_format=float_format)


@dataclass
class ExpressionDataset:
    """A gene x sample expression matrix plus optional sample metadata.

    Parameters
    ----------
    counts:
        Rows are gene identifiers, columns are sample identifiers. Values are
        expected to be raw (or pseudo-raw) counts for the DEB/limma-style
        pipeline; a ``kind`` attribute records what they actually are.
    meta:
        Sample-level table indexed by the same sample identifiers.
    gene_meta:
        Optional gene-level annotations (e.g. gene symbols, biotype).
    kind:
        ``"counts"`` for integer-like counts, ``"lognorm"`` for already
        normalised/log expression values.
    """

    counts: pd.DataFrame
    meta: pd.DataFrame | None = None
    gene_meta: pd.DataFrame | None = None
    kind: str = "counts"
    name: str = "dataset"
    source: dict = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.counts.columns.has_duplicates:
            dup = self.counts.columns[self.counts.columns.duplicated()].unique().tolist()
            raise ValueError(f"{self.name}: duplicated sample identifiers: {dup[:5]}")
        if self.counts.index.has_duplicates:
            dup = self.counts.index[self.counts.index.duplicated()].unique().tolist()
            raise ValueError(f"{self.name}: duplicated gene identifiers: {dup[:5]}")
        if self.meta is not None:
            missing = set(self.counts.columns) - set(self.meta.index)
            if missing:
                raise ValueError(f"{self.name}: {len(missing)} samples missing from metadata")

    # -- basic views -------------------------------------------------------
    @property
    def n_genes(self) -> int:
        return int(self.counts.shape[0])

    @property
    def n_samples(self) -> int:
        return int(self.counts.shape[1])

    def subset(self, genes: Iterable[str] | None = None, samples: Iterable[str] | None = None) -> "ExpressionDataset":
        c = self.counts
        if genes is not None:
            c = c.loc[[g for g in genes if g in c.index]]
        if samples is not None:
            c = c.loc[:, [s for s in samples if s in c.columns]]
        meta = None if self.meta is None else self.meta.loc[c.columns]
        gm = None if self.gene_meta is None else self.gene_meta.loc[[g for g in c.index if g in self.gene_meta.index]]
        return ExpressionDataset(c, meta, gm, self.kind, self.name, dict(self.source))

    def group_vector(self, column: str) -> pd.Series:
        if self.meta is None or column not in self.meta.columns:
            raise KeyError(f"metadata column {column!r} not found")
        return self.meta.loc[self.counts.columns, column]

    def describe(self) -> str:
        g, s = self.n_genes, self.n_samples
        parts = [f"{self.name}: {g} genes x {s} samples", f"kind={self.kind}"]
        if self.meta is not None:
            for col in self.meta.columns:
                if self.meta[col].dtype == object and self.meta[col].nunique() <= 12:
                    vc = self.meta[col].value_counts().to_dict()
                    parts.append(f"  {col}: {vc}")
        nulls = int(self.counts.isna().sum().sum())
        if nulls:
            parts.append(f"  WARNING: {nulls} missing values")
        return "\n".join(parts)


@dataclass
class Provenance:
    """JSON-serialisable record of how an analysis artefact was produced."""

    tool: str = "stade"
    version: str = "0.1.0"
    command: str = ""
    inputs: dict = field(default_factory=dict)
    parameters: dict = field(default_factory=dict)
    started: str = ""
    finished: str = ""
    notes: dict = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not self.started:
            self.started = time.strftime("%Y-%m-%dT%H:%M:%S%z")

    def finish(self, **notes) -> "Provenance":
        self.finished = time.strftime("%Y-%m-%dT%H:%M:%S%z")
        self.notes.update(notes)
        return self

    def add_input(self, label: str, path: str | Path, with_hash: bool = True) -> None:
        p = Path(path)
        entry = {"path": str(p), "bytes": p.stat().st_size if p.exists() else None}
        if with_hash and p.exists() and p.stat().st_size < 2 << 30:
            entry["sha256"] = file_sha256(p)
        self.inputs[label] = entry

    def to_dict(self) -> dict:
        return {
            "tool": self.tool,
            "version": self.version,
            "command": self.command,
            "started": self.started,
            "finished": self.finished,
            "inputs": self.inputs,
            "parameters": self.parameters,
            "notes": self.notes,
            "environment": {
                "python": sys.version.split()[0],
                "platform": platform.platform(),
                "numpy": np.__version__,
                "pandas": pd.__version__,
            },
        }

    def write(self, path: str | Path) -> None:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(self.to_dict(), indent=2, ensure_ascii=False, default=str), encoding="utf-8")
