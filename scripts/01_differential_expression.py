"""Step 1: KDM1A knockdown contrasts and the ECM30 discovery set.

Combines the GSE58492 author DESeq contrasts (siL1, siL2) with a Welch
comparison of GSE30775, defines the shared response at BH q < 0.05 and a
1.2-fold change, and reconstructs ECM30, ECM29, Adaptor10, Adaptor9 and
Prolif20.

Usage:
    python scripts/01_differential_expression.py [--data DIR] [--out DIR]
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import common as C

CONTRASTS = ["GSE58492_siL1", "GSE58492_siL2", "GSE30775_siKDM1A"]


def read_manifest(data_dir):
    return pd.read_csv(data_dir / "processed" / "datasets.tsv", sep="\t")


def load_processed(data_dir, name):
    return pd.read_csv(data_dir / "processed" / f"{name}.csv.gz", index_col=0)


def authors_plus_welch(data_dir, manifest):
    author = pd.read_csv(data_dir / "processed" / "GSE58492_author_deseq.csv", index_col=0)
    row = manifest.set_index("dataset").loc["GSE30775"]
    matrix = load_processed(data_dir, "GSE30775")
    contrast = C.welch_contrast(matrix, row["control"].split(";"), row["case"].split(";"))

    combined = pd.DataFrame(index=author.index)
    combined["GSE58492_siL1"] = author["siL1_log2FC"]
    combined["GSE58492_siL2"] = author["siL2_log2FC"]
    combined["GSE58492_siL1_q"] = author["siL1_q"]
    combined["GSE58492_siL2_q"] = author["siL2_q"]
    combined["GSE30775_siKDM1A"] = contrast["log2FC"].reindex(combined.index)
    combined["GSE30775_siKDM1A_q"] = contrast["q"].reindex(combined.index)
    return combined, author, contrast


def shared_response(combined):
    effects = combined[CONTRASTS].to_numpy(dtype=float)
    qvalues = combined[[c + "_q" for c in CONTRASTS]].to_numpy(dtype=float)
    tested = np.isfinite(effects).all(axis=1) & np.isfinite(qvalues).all(axis=1)
    significant = ((qvalues < 0.05).all(axis=1)
                   & (np.abs(effects) >= C.LOG2_THRESHOLD).all(axis=1))
    induced = tested & significant & (effects > 0).all(axis=1)
    reduced = tested & significant & (effects < 0).all(axis=1)
    return tested, induced, reduced


def reconstruct_ecm30(author, hallmark):
    emt = set(hallmark["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"])
    adhesion = set(C.ADHESION13)
    minimum = np.minimum(author["siL1_log2FC"], author["siL2_log2FC"])
    eligible = ((author["siL1_q"] < 0.05) & (author["siL2_q"] < 0.05)
                & (minimum >= 0.25) & (author["mean_count"] >= 5.0)
                & author.index.isin(emt | adhesion))
    ranked = minimum[eligible].sort_values(ascending=False)
    members = ranked.head(30)
    detail = pd.DataFrame({
        "min_log2FC": ranked,
        "mean_count": author.loc[ranked.index, "mean_count"],
        "in_hallmark_emt": [g in emt for g in ranked.index],
        "in_adhesion13": [g in adhesion for g in ranked.index],
        "selected": [g in members.index for g in ranked.index],
    })
    return members, detail


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", type=Path, default=C.DATA)
    ap.add_argument("--out", type=Path, default=C.RESULTS)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    manifest = read_manifest(args.data)
    hallmark = C.read_gmt(args.data / "msigdb" / "h.all.v2024.1.Hs.symbols.gmt")
    combined, author, df30775 = authors_plus_welch(args.data, manifest)
    tested, induced, reduced = shared_response(combined)

    table = combined[CONTRASTS].copy()
    table["tested_all_three"] = tested
    table["shared_induced"] = induced
    table["shared_reduced"] = reduced
    table.to_csv(args.out / "contrasts.csv")

    shared = combined.loc[induced | reduced, CONTRASTS].copy()
    shared["direction"] = np.where(induced[induced | reduced], "up", "down")
    shared.sort_values("direction").to_csv(args.out / "shared_response.csv")

    members, detail = reconstruct_ecm30(author, hallmark)
    detail.to_csv(args.out / "ecm30_candidates.csv")

    prolif_genes = (set(hallmark["HALLMARK_E2F_TARGETS"])
                    | set(hallmark["HALLMARK_G2M_CHECKPOINT"]))
    prolif = df30775[(df30775["q"] < 0.05) & (df30775["log2FC"] <= -0.5)]
    prolif = prolif[prolif.index.isin(prolif_genes)].sort_values("log2FC").head(20)

    sets = {
        "ECM30": list(members.index),
        "ECM29": [g for g in members.index if g not in C.ADAPTOR10],
        "Adaptor10": C.ADAPTOR10,
        "Adaptor9": C.adaptor9(),
        "Prolif20": list(prolif.index),
        "Adhesion13": C.ADHESION13,
        "Integrin_A": C.INTEGRIN_A,
        "Integrin_B": C.INTEGRIN_B,
    }
    (args.out / "gene_sets.json").write_text(json.dumps(sets, indent=2), encoding="utf-8")

    summary = {
        "genes_with_estimates_in_all_three": int(tested.sum()),
        "shared_induced": int(induced.sum()),
        "shared_reduced": int(reduced.sum()),
        "ecm30_eligible": int(detail.shape[0]),
        "ecm30_selected": int(len(members)),
        "prolif20_selected": int(len(prolif)),
        "fold_change_threshold": C.FOLD_THRESHOLD,
        "log2_threshold": round(C.LOG2_THRESHOLD, 4),
    }
    C.write_json(args.out / "summary_step1.json", summary)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    sys.exit(main())
