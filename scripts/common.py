"""Shared helpers for the KDM1A/ITGA5 analysis scripts."""

import gzip
import io
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
DATA = Path(os.environ.get("KDM1A_DATA", ROOT / "data"))
RESULTS = Path(os.environ.get("KDM1A_RESULTS", ROOT / "results"))

# Fixed gene memberships (defined in the discovery step and used unchanged
# downstream). ECM29 and Adaptor9 drop ITGA5, the only shared member.
ECM30 = [
    "PDGFRB", "COL7A1", "COL17A1", "ITGAV", "LRP1", "IGFBP4", "TGM2", "TNC",
    "SPARC", "FBLN1", "COL5A1", "SNAI2", "LAMC2", "COL1A1", "FBN1", "ECM1",
    "ITGA3", "NOTCH2", "L1CAM", "MYL9", "COL6A2", "TGFBI", "LAMC1", "PLAUR",
    "LAMB3", "FSTL1", "GLIPR1", "RHOB", "ITGA2", "ITGA5",
]
ADAPTOR10 = ["ITGA5", "ITGB1", "PXN", "TLN1", "VCL", "FERMT2", "ILK", "PTK2",
             "BCAR1", "ACTN1"]
ADHESION13 = ["L1CAM", "ITGB4", "ITGA6", "ITGA3", "COL17A1", "LAMA3", "LAMB3",
              "CDH2", "FSCN1", "LAMC2", "TNC", "TGFBI", "SERPINE1"]
INTEGRIN_A = ["ITGA1", "ITGA2", "ITGA2B", "ITGA3", "ITGA4", "ITGA5", "ITGA6",
              "ITGA7", "ITGA8", "ITGA9", "ITGA10", "ITGA11", "ITGAD", "ITGAE",
              "ITGAL", "ITGAM", "ITGAV", "ITGAX"]
INTEGRIN_B = ["ITGB1", "ITGB2", "ITGB3", "ITGB4", "ITGB5", "ITGB6", "ITGB7",
              "ITGB8"]

FOLD_THRESHOLD = 1.2
LOG2_THRESHOLD = float(np.log2(FOLD_THRESHOLD))


def ecm29():
    return [g for g in ECM30 if g not in ADAPTOR10]


def adaptor9():
    return [g for g in ADAPTOR10 if g != "ITGA5"]


def bh(pvalues):
    """Benjamini-Hochberg adjusted p values, NaN preserved."""
    p = np.asarray(pvalues, dtype=float)
    q = np.full(p.shape, np.nan)
    idx = np.flatnonzero(np.isfinite(p))
    if idx.size == 0:
        return q
    order = idx[np.argsort(p[idx], kind="mergesort")]
    m = order.size
    scaled = p[order] * m / np.arange(1, m + 1)
    scaled = np.minimum.accumulate(scaled[::-1])[::-1]
    q[order] = np.clip(scaled, 0.0, 1.0)
    return q


def welch_contrast(matrix, control, case):
    """Gene-wise Welch test between two sample groups, genome-wide BH."""
    a = matrix[control].to_numpy(dtype=float)
    b = matrix[case].to_numpy(dtype=float)
    t, p = stats.ttest_ind(b, a, axis=1, equal_var=False, nan_policy="omit")
    out = pd.DataFrame({
        "log2FC": np.nanmean(b, axis=1) - np.nanmean(a, axis=1),
        "t": t,
        "p": p,
    }, index=matrix.index)
    out["q"] = bh(out["p"].to_numpy())
    out["n_control"] = len(control)
    out["n_case"] = len(case)
    return out


def read_series_matrix(path):
    """Read the value table of a GEO series matrix file."""
    with gzip.open(path, "rt", errors="replace") as fh:
        text = fh.read()
    block = text.split("!series_matrix_table_begin\n")[1]
    block = block.split("!series_matrix_table_end")[0]
    table = pd.read_csv(io.StringIO(block), sep="\t", index_col=0)
    table.index.name = "probe"
    return table.apply(pd.to_numeric, errors="coerce")


def read_platform_annotation(path):
    """Probe to gene symbol map from a GEO platform annotation file."""
    with gzip.open(path, "rt", errors="replace") as fh:
        text = fh.read()
    block = text.split("!platform_table_begin\n")[1]
    block = block.split("!platform_table_end")[0]
    table = pd.read_csv(io.StringIO(block), sep="\t", dtype=str, low_memory=False)
    keep = table[["ID", "Gene symbol"]].dropna()
    keep["Gene symbol"] = (keep["Gene symbol"].str.split("///").str[0]
                           .str.strip().str.upper())
    keep = keep[keep["Gene symbol"] != ""]
    return keep.drop_duplicates("ID").set_index("ID")["Gene symbol"]


def collapse_to_symbols(matrix, mapping):
    """Map probes to symbols, keeping the probe with the highest mean signal."""
    frame = matrix.copy()
    frame["gene"] = frame.index.map(mapping)
    frame = frame.dropna(subset=["gene"])
    frame["mean"] = frame.drop(columns="gene").mean(axis=1)
    frame = frame.sort_values("mean", ascending=False).drop_duplicates("gene")
    return frame.drop(columns="mean").set_index("gene").sort_index()


def read_gmt(path):
    sets = {}
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) > 2:
                sets[fields[0]] = [g for g in fields[2:] if g]
    return sets


def percentile_scores(matrix, gene_sets):
    """Mean within-profile percentile rank of each gene set.

    Ranks are computed across all measured genes in a profile, so ties
    (including zeros) receive their average rank. Values are scaled to 0-1.
    """
    values = matrix.to_numpy(dtype=float)
    ranks = stats.rankdata(values, axis=0, method="average")
    position = {g: i for i, g in enumerate(matrix.index)}
    ng = matrix.shape[0]
    out = {}
    for name, genes in gene_sets.items():
        rows = [position[g] for g in genes if g in position]
        if rows:
            out[name] = ranks[rows].mean(axis=0) / ng
    return pd.DataFrame(out, index=matrix.columns)


def partial_spearman(a, b, controls):
    """Spearman correlation after residualising ranks on ranked controls."""
    ra = stats.rankdata(a)
    rb = stats.rankdata(b)
    design = np.column_stack([np.ones(len(ra))] + [stats.rankdata(c) for c in controls])
    ra = ra - design @ np.linalg.lstsq(design, ra, rcond=None)[0]
    rb = rb - design @ np.linalg.lstsq(design, rb, rcond=None)[0]
    if ra.std() < 1e-12 or rb.std() < 1e-12:
        return np.nan
    return float(np.corrcoef(ra, rb)[0, 1])


def fisher_z_test(coefficients):
    """One-sample test of Fisher-transformed coefficients against zero."""
    r = np.asarray(coefficients, dtype=float)
    r = r[np.isfinite(r)]
    r = np.clip(r, -0.999999, 0.999999)
    z = np.arctanh(r)
    if z.size < 2:
        return np.nan, np.nan
    t, p = stats.ttest_1samp(z, 0.0)
    return float(t), float(p)


def bootstrap_mean_ci(values, n_boot=10000, seed=37004, alpha=0.05):
    rng = np.random.default_rng(seed)
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size < 2:
        return np.nan, np.nan
    draws = rng.choice(values, size=(n_boot, values.size), replace=True).mean(axis=1)
    return float(np.quantile(draws, alpha / 2)), float(np.quantile(draws, 1 - alpha / 2))


def write_json(path, payload):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
