"""Step 5: baseline expression against CRISPR fitness in breast models.

Reads DepMap Public 24Q4 model annotations, expression and Chronos gene
effects, keeps breast lineage models and the nested TNBC subset, and correlates
each gene's expression with its own gene effect. More negative gene effects
indicate a larger fitness cost.

Inputs:
    data/depmap/Model.csv
    data/depmap/CRISPRGeneEffect.csv
    data/depmap/OmicsExpressionProteinCodingGenesTPMLogp1.csv

Usage:
    python scripts/05_depmap_dependency.py [--data DIR] [--out DIR]
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

import common as C

GENES = ["KDM1A", "ITGA5", "ITGB1", "ITGA2", "FN1"]


def read_gene_columns(path, genes, id_column="ModelID"):
    header = pd.read_csv(path, nrows=0).columns.tolist()
    wanted = [id_column]
    for gene in genes:
        for column in header:
            if column.split(" (")[0] == gene:
                wanted.append(column)
    return pd.read_csv(path, usecols=wanted).set_index(id_column)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", type=Path, default=C.DATA)
    ap.add_argument("--out", type=Path, default=C.RESULTS)
    ap.add_argument("--genes", nargs="*", default=GENES)
    args = ap.parse_args()

    depmap = args.data / "depmap"
    models = pd.read_csv(depmap / "Model.csv")
    breast = models[models["OncotreeLineage"].astype(str).str.contains("Breast", case=False)]
    tnbc_codes = breast["OncotreeCode"].astype(str).str.upper().str.contains("TNBC|TRIPLE")
    tnbc = breast[tnbc_codes]

    expression = read_gene_columns(
        depmap / "OmicsExpressionProteinCodingGenesTPMLogp1.csv", args.genes)
    fitness = read_gene_columns(depmap / "CRISPRGeneEffect.csv", args.genes)
    expression.columns = [c.split(" (")[0] for c in expression.columns]
    fitness.columns = [c.split(" (")[0] for c in fitness.columns]

    frames = []
    for cohort, members in (("breast", breast["ModelID"]), ("tnbc", tnbc["ModelID"])):
        ids = [m for m in members if m in expression.index and m in fitness.index]
        for gene in args.genes:
            if gene not in expression.columns or gene not in fitness.columns:
                continue
            pair = pd.DataFrame({"expression": expression.loc[ids, gene],
                                 "fitness": fitness.loc[ids, gene]}).dropna()
            if pair.shape[0] < 10:
                continue
            rho, p = stats.spearmanr(pair["expression"], pair["fitness"])
            frames.append({"cohort": cohort, "gene": gene, "n_models": pair.shape[0],
                           "rho": rho, "p": p})

    table = pd.DataFrame(frames)
    table["q"] = np.nan
    for cohort, block in table.groupby("cohort"):
        table.loc[block.index, "q"] = C.bh(block["p"].to_numpy())
    table.to_csv(args.out / "depmap_dependency.csv", index=False)

    summary = {
        "breast_models": int(breast.shape[0]),
        "tnbc_models": int(tnbc.shape[0]),
        "genes_tested": sorted(table["gene"].unique().tolist()) if not table.empty else [],
    }
    C.write_json(args.out / "summary_step5.json", summary)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    sys.exit(main())
