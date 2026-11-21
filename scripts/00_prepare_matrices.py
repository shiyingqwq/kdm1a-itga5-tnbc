"""Step 0: build processed expression matrices from the public GEO files.

Writes one genes x samples table per dataset plus a dataset manifest. Only the
KDM1A arms are kept where a series contains additional perturbations.

Usage:
    python scripts/00_prepare_matrices.py [--data DIR]
"""

import argparse
import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import common as C


def read_sample_titles(path):
    titles, accessions = [], []
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!Sample_title"):
                titles = [v.strip('"') for v in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_geo_accession"):
                accessions = [v.strip('"') for v in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!series_matrix_table_begin"):
                break
    return pd.DataFrame({"gsm": accessions, "title": titles})


def median_ratio_normalise(counts):
    positive = (counts > 0).all(axis=1)
    geometric = np.exp(np.log(counts[positive]).mean(axis=1))
    size_factors = counts[positive].div(geometric, axis=0).median(axis=0)
    return counts.div(size_factors, axis=1), size_factors


def write(matrix, name, out_dir, control, case, role):
    matrix.to_csv(out_dir / f"{name}.csv.gz", compression="gzip")
    return {"dataset": name, "control": ";".join(control), "case": ";".join(case),
            "role": role, "n_genes": matrix.shape[0], "n_samples": matrix.shape[1]}


def prepare_gse58492(geo, out_dir, annotation_unused=None):
    raw = pd.read_csv(geo / "GSE58492" / "GSE58492_siLSD1_results.txt.gz",
                      sep="\t", skiprows=2, header=None)
    raw.columns = ["gene", "siC-1", "siC-2", "siL1-1", "siL1-2", "fc_siL1",
                   "padj_siL1", "siL2-1", "siL2-2", "fc_siL2", "padj_siL2"]
    raw["gene"] = raw["gene"].str.upper()
    raw = raw.drop_duplicates("gene").set_index("gene")
    samples = ["siC-1", "siC-2", "siL1-1", "siL1-2", "siL2-1", "siL2-2"]
    counts = raw[samples].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    normalised, _ = median_ratio_normalise(counts)
    matrix = np.log2(normalised + 1.0)

    author = pd.DataFrame({
        "siL1_log2FC": np.log2(pd.to_numeric(raw["fc_siL1"], errors="coerce").clip(lower=1e-12)),
        "siL1_q": pd.to_numeric(raw["padj_siL1"], errors="coerce"),
        "siL2_log2FC": np.log2(pd.to_numeric(raw["fc_siL2"], errors="coerce").clip(lower=1e-12)),
        "siL2_q": pd.to_numeric(raw["padj_siL2"], errors="coerce"),
        "mean_count": counts.mean(axis=1),
    })
    author.to_csv(out_dir / "GSE58492_author_deseq.csv")
    matrix.to_csv(out_dir / "GSE58492.csv.gz", compression="gzip")
    return {"dataset": "GSE58492", "control": "siC-1;siC-2",
            "case": "siL1-1;siL1-2;siL2-1;siL2-2", "role": "discovery",
            "n_genes": matrix.shape[0], "n_samples": matrix.shape[1]}


def prepare_array(geo, name, platform_file, out_dir, select, control_regex, case_regex,
                  role):
    matrix_path = geo / name / f"{name}_series_matrix.txt.gz"
    meta = read_sample_titles(matrix_path)
    if select:
        meta = meta[meta["gsm"].isin(select)].reset_index(drop=True)
    chosen = meta[meta["title"].str.contains(control_regex, case=False, regex=True)
                  | meta["title"].str.contains(case_regex, case=False, regex=True)]
    control = chosen.loc[chosen["title"].str.contains(control_regex, case=False), "gsm"].tolist()
    case = chosen.loc[chosen["title"].str.contains(case_regex, case=False), "gsm"].tolist()
    case = [g for g in case if g not in control]

    signal = C.read_series_matrix(matrix_path)[control + case]
    annotation = C.read_platform_annotation(geo / platform_file)
    matrix = C.collapse_to_symbols(np.log2(signal.clip(lower=1.0)), annotation)
    matrix.to_csv(out_dir / f"{name}.csv.gz", compression="gzip")
    return {"dataset": name, "control": ";".join(control), "case": ";".join(case),
            "role": role, "n_genes": matrix.shape[0], "n_samples": matrix.shape[1]}


def prepare_gse154382(geo, out_dir):
    files = sorted((geo / "GSE154382").glob("GSM*_Genes_ReadCount.txt.gz"))
    columns = {}
    for path in files:
        label = path.name.split("_", 2)[2].replace("_Genes_ReadCount.txt.gz", "")
        table = pd.read_csv(path, sep="\t", header=None, index_col=0,
                            names=["gene", label])
        columns[label] = table[label]
    counts = pd.DataFrame(columns).fillna(0.0)
    counts.index = counts.index.str.upper()
    counts = counts.groupby(level=0).sum()
    normalised, _ = median_ratio_normalise(counts)
    matrix = np.log2(normalised + 1.0)
    control = [c for c in matrix.columns if "cont" in c]
    case = [c for c in matrix.columns if "KDM1A" in c]
    matrix.to_csv(out_dir / "GSE154382.csv.gz", compression="gzip")
    return {"dataset": "GSE154382", "control": ";".join(control), "case": ";".join(case),
            "role": "independent_validation", "n_genes": matrix.shape[0],
            "n_samples": matrix.shape[1]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", type=Path, default=C.DATA)
    args = ap.parse_args()
    geo = args.data / "geo"
    out_dir = args.data / "processed"
    out_dir.mkdir(parents=True, exist_ok=True)

    rows = [prepare_gse58492(geo, out_dir)]
    rows.append(prepare_array(geo, "GSE30775", "GPL570.annot.gz", out_dir, None,
                              "control siRNA", r"LSD1 siRNA", "independent_validation"))
    rows.append(prepare_array(geo, "GSE72687", "GPL570.annot.gz", out_dir,
                              [f"GSM29113{i}" for i in range(78, 84)],
                              "scramble", r"LSD1 shRNA", "independent_validation"))
    rows.append(prepare_array(geo, "GSE95165", "GSE95165/GPL10558.annot.gz", out_dir, None,
                              "control siRNA", r"LSD1 siRNA", "context"))
    rows.append(prepare_gse154382(geo, out_dir))

    manifest = pd.DataFrame(rows)
    manifest.to_csv(out_dir / "datasets.tsv", sep="\t", index=False)
    print(manifest.to_string(index=False))


if __name__ == "__main__":
    sys.exit(main())
