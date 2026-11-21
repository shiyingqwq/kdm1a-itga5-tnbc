"""Step 2: program scores and the independent ECM30 comparison.

Scores each program as the mean within-profile percentile rank of its measured
members and compares knockdown with control samples by Welch t test. ECM30 was
derived in GSE58492 and Prolif20 in GSE30775, so those two comparisons are
discovery; the remaining nine form the validation family for BH correction.

Usage:
    python scripts/02_program_scores.py [--data DIR] [--out DIR]
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

import common as C

PROGRAMS = ["ECM30", "ECM29", "Adaptor10", "Adaptor9", "Prolif20"]
DISCOVERY = {("GSE58492", "ECM30"), ("GSE30775", "Prolif20")}
FAMILY_PROGRAMS = ["ECM30", "Prolif20"]


def comparisons(manifest):
    """Dataset, control and case samples for each analysed contrast."""
    out = []
    for _, entry in manifest.iterrows():
        control = entry["control"].split(";")
        case = entry["case"].split(";")
        if entry["dataset"] == "GSE58492":
            out.append(("GSE58492_siL1", "GSE58492", control, case[:2]))
            out.append(("GSE58492_siL2", "GSE58492", control, case[2:]))
        else:
            out.append((entry["dataset"], entry["dataset"], control, case))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", type=Path, default=C.DATA)
    ap.add_argument("--out", type=Path, default=C.RESULTS)
    args = ap.parse_args()

    manifest = pd.read_csv(args.data / "processed" / "datasets.tsv", sep="\t")
    sets = json.loads((args.out / "gene_sets.json").read_text(encoding="utf-8"))

    rows, sample_rows = [], []
    for label, dataset, control, case in comparisons(manifest):
        matrix = pd.read_csv(args.data / "processed" / f"{dataset}.csv.gz", index_col=0)
        present = {p: [g for g in sets[p] if g in matrix.index] for p in PROGRAMS}
        scores = C.percentile_scores(matrix, present)

        for program in PROGRAMS:
            a = scores.loc[control, program]
            b = scores.loc[case, program]
            t, p = stats.ttest_ind(b, a, equal_var=False)
            role = "discovery" if (dataset, program) in DISCOVERY else "validation"
            rows.append({
                "comparison": label,
                "dataset": dataset,
                "program": program,
                "role": role,
                "n_members_measured": len(present[program]),
                "n_members_defined": len(sets[program]),
                "mean_control": a.mean(),
                "mean_case": b.mean(),
                "difference": b.mean() - a.mean(),
                "t": t,
                "p": p,
            })
            for sample in control + case:
                sample_rows.append({"comparison": label, "program": program,
                                    "sample": sample,
                                    "group": "control" if sample in control else "knockdown",
                                    "score": scores.loc[sample, program]})

    tests = pd.DataFrame(rows)
    tests["q_validation_family"] = np.nan
    mask = (tests["role"] == "validation") & tests["program"].isin(FAMILY_PROGRAMS)
    tests.loc[mask, "q_validation_family"] = C.bh(tests.loc[mask, "p"].to_numpy())

    tests.to_csv(args.out / "program_tests.csv", index=False)
    pd.DataFrame(sample_rows).to_csv(args.out / "program_scores.csv", index=False)

    ecm = tests[tests["comparison"] == "GSE30775"].set_index("program").loc["ECM30"]
    contrasts = pd.read_csv(args.out / "contrasts.csv", index_col=0)
    effects = contrasts["GSE30775_siKDM1A"].reindex(sets["ECM30"]).dropna()

    summary = {
        "ecm30_in_GSE30775": {
            "members_measured": int(ecm["n_members_measured"]),
            "difference": round(float(ecm["difference"]), 4),
            "p": float(ecm["p"]),
            "q_validation_family": float(ecm["q_validation_family"]),
            "concordant_members": int((effects > 0).sum()),
        },
        "validation_family_size": int(mask.sum()),
    }
    C.write_json(args.out / "summary_step2.json", summary)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    sys.exit(main())
