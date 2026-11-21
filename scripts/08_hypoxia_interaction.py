"""Step 8: oxygen and HIF1A/EPAS1 double-knockout effects.

Wild-type and double-knockout cultures were profiled at 20% and 1% oxygen.
Genotype and oxygen contrasts are computed inside each level of the other
factor, so conditions are never pooled. A two-by-two ordinary least squares
model gives the genotype-by-oxygen interaction. Tests run on log2(FPKM + 1)
for genes with a maximum abundance of at least one, and are BH-adjusted
within each comparison.

Inputs:
    data/hypoxia/fpkm.csv    genes x samples, FPKM
    data/hypoxia/design.csv  columns: sample, genotype, oxygen
        genotype in {wild_type, double_knockout}; oxygen in {20, 1}

Usage:
    python scripts/08_hypoxia_interaction.py [--data DIR] [--out DIR]
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

import common as C

HIF_TARGETS = ["VEGFA", "SLC2A1", "PGK1", "LDHA", "BNIP3", "PDK1", "ANGPTL4",
               "CA9", "ADM", "NDRG1", "P4HA1", "EGLN3"]


def welch(matrix, first, second):
    a = matrix[first].to_numpy(dtype=float)
    b = matrix[second].to_numpy(dtype=float)
    t, p = stats.ttest_ind(b, a, axis=1, equal_var=False, nan_policy="omit")
    table = pd.DataFrame({"log2FC": np.nanmean(b, axis=1) - np.nanmean(a, axis=1),
                          "t": t, "p": p}, index=matrix.index)
    table["q"] = C.bh(table["p"].to_numpy())
    return table


def interaction(matrix, design):
    values = matrix[design["sample"]].to_numpy(dtype=float).T
    genotype = (design["genotype"] == "double_knockout").to_numpy(dtype=float)
    oxygen = (design["oxygen"].astype(float) == 1).to_numpy(dtype=float)
    design_matrix = np.column_stack([np.ones(len(design)), genotype, oxygen,
                                     genotype * oxygen])
    beta, *_ = np.linalg.lstsq(design_matrix, values, rcond=None)
    residual = values - design_matrix @ beta
    dof = values.shape[0] - design_matrix.shape[1]
    sigma2 = (residual ** 2).sum(axis=0) / dof
    covariance = np.linalg.inv(design_matrix.T @ design_matrix)
    se = np.sqrt(sigma2 * covariance[3, 3])
    with np.errstate(divide="ignore", invalid="ignore"):
        t = beta[3] / se
    p = 2 * stats.t.sf(np.abs(t), dof)
    table = pd.DataFrame({"interaction": beta[3], "t": t, "p": p}, index=matrix.index)
    table["q"] = C.bh(table["p"].to_numpy())
    return table


def signed_hif_score(matrix):
    present = [g for g in HIF_TARGETS if g in matrix.index]
    block = matrix.loc[present]
    centred = block.sub(block.mean(axis=1), axis=0)
    return centred.mean(axis=0)


def samples(design, genotype=None, oxygen=None):
    mask = pd.Series(True, index=design.index)
    if genotype is not None:
        mask &= design["genotype"] == genotype
    if oxygen is not None:
        mask &= design["oxygen"].astype(float) == oxygen
    return design.loc[mask, "sample"].tolist()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", type=Path, default=C.DATA)
    ap.add_argument("--out", type=Path, default=C.RESULTS)
    args = ap.parse_args()

    fpkm = pd.read_csv(args.data / "hypoxia" / "fpkm.csv", index_col=0)
    design = pd.read_csv(args.data / "hypoxia" / "design.csv")

    # Duplicate symbols are averaged before the log transform.
    fpkm = fpkm.groupby(level=0).mean()
    matrix = np.log2(fpkm.clip(lower=0) + 1.0)
    matrix = matrix.loc[fpkm.max(axis=1) >= 1.0]

    comparisons = {
        "knockout_vs_wild_type_20pct": welch(
            matrix, samples(design, "wild_type", 20), samples(design, "double_knockout", 20)),
        "knockout_vs_wild_type_1pct": welch(
            matrix, samples(design, "wild_type", 1), samples(design, "double_knockout", 1)),
        "hypoxia_vs_normoxia_wild_type": welch(
            matrix, samples(design, "wild_type", 20), samples(design, "wild_type", 1)),
        "hypoxia_vs_normoxia_double_knockout": welch(
            matrix, samples(design, "double_knockout", 20), samples(design, "double_knockout", 1)),
        "interaction": interaction(matrix, design),
    }
    for name, table in comparisons.items():
        table.insert(0, "comparison", name)
    combined = pd.concat(comparisons.values())
    combined.index.name = "gene"
    combined.to_csv(args.out / "hypoxia_gene_tests.csv")

    score = signed_hif_score(matrix)
    score_frame = pd.DataFrame({"sample": score.index,
                                "signed_hif_score": score.to_numpy()})
    score_frame = score_frame.merge(design, on="sample", how="left")
    score_frame.to_csv(args.out / "hypoxia_hif_score.csv", index=False)

    focal = [g for g in ["ITGA5", "ITGB1", "FN1", "VEGFA", "CA9"] if g in matrix.index]
    ko = comparisons["knockout_vs_wild_type_1pct"]
    hyp = comparisons["hypoxia_vs_normoxia_wild_type"]
    inter = comparisons["interaction"]
    summary = {
        "genes_tested": int(matrix.shape[0]),
        "samples": int(design.shape[0]),
        "hif_targets_measured": int(len([g for g in HIF_TARGETS if g in matrix.index])),
        "hif_score_normoxia": float(score.loc[samples(design, oxygen=20)].mean()),
        "hif_score_hypoxia": float(score.loc[samples(design, oxygen=1)].mean()),
        "focal": {
            gene: {
                "log2FC_hypoxia_wild_type": float(hyp.loc[gene, "log2FC"]),
                "q_hypoxia_wild_type": float(hyp.loc[gene, "q"]),
                "log2FC_dko_vs_wt_1pct": float(ko.loc[gene, "log2FC"]),
                "q_dko_vs_wt_1pct": float(ko.loc[gene, "q"]),
                "interaction": float(inter.loc[gene, "interaction"]),
                "q_interaction": float(inter.loc[gene, "q"]),
            } for gene in focal
        },
    }
    C.write_json(args.out / "summary_step8.json", summary)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    sys.exit(main())
