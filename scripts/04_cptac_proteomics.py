"""Step 4: candidate protein abundance against a leave-self-out ECM score.

For every measured protein the ECM30 score is rebuilt from the other members,
so the tested protein never contributes to its own score. Correlations are
Spearman coefficients with BH correction across all tested proteins; the TNBC
subgroup is analysed descriptively with the same score.

Inputs:
    data/cptac/protein_matrix.csv    samples x proteins, normalised log ratios
    data/cptac/sample_metadata.csv   columns: sample_id, tnbc_group

Usage:
    python scripts/04_cptac_proteomics.py [--data DIR] [--out DIR]
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

import common as C


def ecm_score(matrix, members, exclude=None):
    genes = [g for g in members if g in matrix.columns and g != exclude]
    block = matrix[genes]
    z = (block - block.mean()) / block.std(ddof=0)
    return z.mean(axis=1), len(genes)


def correlate(matrix, candidates, members):
    rows = []
    for protein in candidates:
        score, n_members = ecm_score(matrix, members, exclude=protein)
        pair = pd.DataFrame({"protein": matrix[protein], "score": score}).dropna()
        if pair.shape[0] < 10:
            continue
        rho, p = stats.spearmanr(pair["protein"], pair["score"])
        rows.append({"protein": protein, "n_samples": pair.shape[0],
                     "n_score_members": n_members, "rho": rho, "p": p})
    table = pd.DataFrame(rows)
    table["q"] = C.bh(table["p"].to_numpy())
    return table.sort_values("p").reset_index(drop=True)


def bootstrap_and_leave_one_out(matrix, protein, members, n_boot=10000, seed=37004):
    score, _ = ecm_score(matrix, members, exclude=protein)
    pair = pd.DataFrame({"protein": matrix[protein], "score": score}).dropna()
    rng = np.random.default_rng(seed)
    draws = []
    values = pair.to_numpy()
    for _ in range(n_boot):
        pick = rng.integers(0, values.shape[0], values.shape[0])
        draws.append(stats.spearmanr(values[pick, 0], values[pick, 1])[0])
    loo = []
    for i in range(values.shape[0]):
        keep = np.ones(values.shape[0], dtype=bool)
        keep[i] = False
        loo.append(stats.spearmanr(values[keep, 0], values[keep, 1])[0])
    return {
        "n_samples": int(values.shape[0]),
        "rho": float(stats.spearmanr(values[:, 0], values[:, 1])[0]),
        "ci_low": float(np.nanquantile(draws, 0.025)),
        "ci_high": float(np.nanquantile(draws, 0.975)),
        "loo_min": float(np.nanmin(loo)),
        "loo_max": float(np.nanmax(loo)),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", type=Path, default=C.DATA)
    ap.add_argument("--out", type=Path, default=C.RESULTS)
    ap.add_argument("--candidate", default="ITGA5")
    args = ap.parse_args()

    matrix = pd.read_csv(args.data / "cptac" / "protein_matrix.csv", index_col=0)
    meta = pd.read_csv(args.data / "cptac" / "sample_metadata.csv")
    matrix = matrix.reindex(meta["sample_id"])
    matrix.index = meta["sample_id"].to_numpy()

    candidates = [g for g in C.ECM30 if g in matrix.columns]
    full = correlate(matrix, candidates, C.ECM30)

    tnbc_samples = meta.loc[meta["tnbc_group"].astype(str).str.lower() == "tnbc", "sample_id"]
    subgroup = correlate(matrix.loc[matrix.index.isin(tnbc_samples)],
                         candidates, C.ECM30)
    subgroup = subgroup.rename(columns={"rho": "rho_tnbc", "p": "p_tnbc", "q": "q_tnbc"})
    merged = full.merge(subgroup[["protein", "n_samples", "rho_tnbc", "p_tnbc", "q_tnbc"]],
                        on="protein", how="left",
                        suffixes=("", "_subgroup"))
    merged.to_csv(args.out / "cptac_protein_correlations.csv", index=False)

    detail = {}
    if args.candidate in matrix.columns:
        detail = bootstrap_and_leave_one_out(matrix, args.candidate, C.ECM30)
        tnbc_matrix = matrix.loc[matrix.index.isin(tnbc_samples)]
        detail["tnbc"] = bootstrap_and_leave_one_out(tnbc_matrix, args.candidate, C.ECM30)

    summary = {
        "samples": int(matrix.shape[0]),
        "tnbc_samples": int(len(tnbc_samples)),
        "tests_in_family": int(full.shape[0]),
        args.candidate: detail,
    }
    C.write_json(args.out / "summary_step4.json", summary)
    print(json.dumps(summary, indent=2, default=float))


if __name__ == "__main__":
    sys.exit(main())
