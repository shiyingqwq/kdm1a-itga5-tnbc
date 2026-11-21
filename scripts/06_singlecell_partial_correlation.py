"""Step 6: within-patient partial correlation of epithelial programs.

Each program is scored per cell against the full measured background, then
ranked. Ranks are residualised on ranked detected-gene counts within patient,
and the patient-level coefficients are tested with a one-sample test of
Fisher-transformed values. Cell-level values are descriptive.

Input:
    data/singlecell/<cohort>.csv.gz
        one row per cell; columns patient, compartment and gene symbols

Usage:
    python scripts/06_singlecell_partial_correlation.py --cohort GSE176078
"""

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

import common as C

PAIRS = [("ECM29", "Adaptor10"), ("Adaptor10", "Prolif20"), ("ECM29", "Prolif20")]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", type=Path, default=C.DATA)
    ap.add_argument("--out", type=Path, default=C.RESULTS)
    ap.add_argument("--cohort", default="GSE176078")
    ap.add_argument("--compartment", default=None)
    args = ap.parse_args()

    frame = pd.read_csv(args.data / "singlecell" / f"{args.cohort}.csv.gz")
    if args.compartment and "compartment" in frame.columns:
        frame = frame[frame["compartment"] == args.compartment]

    patients = frame["patient"].astype(str)
    genes = frame.drop(columns=[c for c in ("patient", "compartment") if c in frame.columns])
    matrix = genes.T
    matrix.index.name = "gene"

    sets = json.loads((args.out / "gene_sets.json").read_text(encoding="utf-8"))
    programs = {name: sets[name] for name in
                ("ECM30", "ECM29", "Adaptor10", "Adaptor9", "Prolif20")}
    scores = C.percentile_scores(matrix, programs)
    scores["patient"] = patients.to_numpy()
    scores["detected"] = (genes > 0).sum(axis=1).to_numpy()

    per_patient = []
    for patient, block in scores.groupby("patient"):
        for first, second in PAIRS:
            rho = C.partial_spearman(block[first], block[second], [block["detected"]])
            per_patient.append({"patient": patient, "pair": f"{first}-{second}",
                                "n_cells": block.shape[0], "rho": rho})
    table = pd.DataFrame(per_patient)

    rows = []
    for pair, block in table.groupby("pair"):
        values = block["rho"].dropna()
        t, p = C.fisher_z_test(values)
        lower, upper = C.bootstrap_mean_ci(values, seed=37004)
        rows.append({"cohort": args.cohort, "pair": pair, "n_patients": int(values.size),
                     "mean_rho": float(values.mean()), "t": t, "p": p,
                     "ci_low": lower, "ci_high": upper})
    result = pd.DataFrame(rows)
    result["q"] = C.bh(result["p"].to_numpy())

    table.to_csv(args.out / f"singlecell_patient_{args.cohort}.csv", index=False)
    result.to_csv(args.out / f"singlecell_tests_{args.cohort}.csv", index=False)
    C.write_json(args.out / f"summary_step6_{args.cohort}.json", {
        "cells": int(frame.shape[0]),
        "patients": int(table["patient"].nunique()),
        "tests": result.to_dict(orient="records"),
    })
    print(result.to_string(index=False))


if __name__ == "__main__":
    sys.exit(main())
