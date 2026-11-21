"""Step 7: spatial associations in Visium sections.

Two quantities are computed per section and tested at the patient level:

1. ECM29 versus Adaptor10, scored as within-spot percentile ranks and
   correlated with total UMI depth as the only control. This is the spatial
   counterpart of the single-cell test in step 6. Ranks are taken over every
   feature in the counts table, which therefore keeps features with no counts.
2. A neighbourhood test: FN1 averaged over the six nearest spots is
   correlated with the focal adhesion score of the spot itself (geometric mean
   of ITGA5 and ITGB1), controlling depth and the ECM29 score.

Inputs per section:
    data/spatial/<section>/counts.csv.gz   spots x genes, raw UMI counts
    data/spatial/<section>/coords.csv      barcode, x, y
    data/spatial/<section>/patient.txt     one patient identifier

Usage:
    python scripts/07_spatial_neighborhood.py [--data DIR] [--out DIR]
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

import common as C

NEIGHBOURS = 6
FOCAL_GENES = ["ITGA5", "ITGB1"]


def neighbour_mean(values, coords, k=NEIGHBOURS):
    tree = cKDTree(coords)
    out = np.full(len(values), np.nan)
    for i, point in enumerate(coords):
        distances, index = tree.query(point, k=k + 1)
        others = [j for j in index if j != i][:k]
        if others:
            out[i] = np.nanmean(values[others])
    return out


def section_analysis(path):
    counts = pd.read_csv(path / "counts.csv.gz", index_col=0)
    coords = pd.read_csv(path / "coords.csv").set_index("barcode").reindex(counts.index)
    patient = (path / "patient.txt").read_text(encoding="utf-8").strip()

    depth = counts.sum(axis=1).to_numpy()
    library = counts.div(np.maximum(depth, 1.0), axis=0) * 1e4
    log_counts = np.log1p(library)

    scores = C.percentile_scores(log_counts.T, {"ECM29": C.ecm29(),
                                                "Adaptor10": C.ADAPTOR10})
    program = C.partial_spearman(scores["ECM29"], scores["Adaptor10"], [depth])

    result = {"patient": patient, "spots": int(counts.shape[0]),
              "rho_ecm_adaptor": program}

    available = [g for g in FOCAL_GENES if g in log_counts.columns]
    if len(available) < 2 or "FN1" not in log_counts.columns:
        result["rho_neighbourhood"] = np.nan
        return result

    focal = np.log1p(np.exp(np.log1p(library[available]).mean(axis=1)) - 1)
    ligand = neighbour_mean(log_counts["FN1"].to_numpy(), coords[["x", "y"]].to_numpy())
    frame = pd.DataFrame({"focal": focal, "ligand": ligand, "depth": depth,
                          "ecm": scores["ECM29"].to_numpy()}).dropna()
    result["rho_neighbourhood"] = (
        C.partial_spearman(frame["ligand"], frame["focal"], [frame["depth"], frame["ecm"]])
        if frame.shape[0] >= 50 else np.nan)
    return result


def summarise(values):
    values = pd.Series(values).dropna()
    if values.empty:
        return {}
    t, p = C.fisher_z_test(values)
    return {"patients": int(values.size), "median_rho": float(values.median()),
            "mean_rho": float(values.mean()), "t": t, "p": p,
            "positive_patients": int((values > 0).sum())}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", type=Path, default=C.DATA)
    ap.add_argument("--out", type=Path, default=C.RESULTS)
    args = ap.parse_args()

    sections = sorted(p for p in (args.data / "spatial").iterdir() if p.is_dir())
    rows = [r for r in (section_analysis(p) for p in sections) if r]
    table = pd.DataFrame(rows, columns=["patient", "spots", "rho_ecm_adaptor",
                                        "rho_neighbourhood"])
    table.to_csv(args.out / "spatial_patient_associations.csv", index=False)
    if table.empty:
        C.write_json(args.out / "summary_step7.json", {"sections": 0})
        print("no eligible section")
        return

    summary = {
        "sections": int(table.shape[0]),
        "spots": int(table["spots"].sum()),
        "ecm_adaptor": summarise(table["rho_ecm_adaptor"]),
        "neighbourhood": summarise(table["rho_neighbourhood"]),
    }
    C.write_json(args.out / "summary_step7.json", summary)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    sys.exit(main())
