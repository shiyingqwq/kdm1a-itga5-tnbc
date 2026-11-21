"""Step 3: over-representation of the dual-siRNA induced list.

The induced list requires BH q < 0.05 in both GSE58492 siRNAs and a 1.2-fold
change in the same direction. The background is every gene with an estimate in
both siRNAs. Terms are restricted to 10-500 genes after intersecting with the
background, and BH correction is applied within each library.

Usage:
    python scripts/03_pathway_overrepresentation.py [--data DIR] [--out DIR]
"""

import argparse
import json
import sys
from pathlib import Path

import pandas as pd
from scipy import stats

import common as C

LIBRARIES = {
    "Reactome_2022": "enrichment_enrichr_Reactome_2022.gmt",
    "KEGG_2021_Human": "enrichment_enrichr_KEGG_2021_Human.gmt",
    "Hallmark_2024_1": "h.all.v2024.1.Hs.symbols.gmt",
}
MIN_TERM = 10
MAX_TERM = 500


def induced_list(author):
    mask = ((author["siL1_q"] < 0.05) & (author["siL2_q"] < 0.05)
            & (author["siL1_log2FC"] >= C.LOG2_THRESHOLD)
            & (author["siL2_log2FC"] >= C.LOG2_THRESHOLD))
    return set(author.index[mask])


def ora(genes, background, sets):
    in_background = genes & background
    records = []
    for term, members in sets.items():
        term_background = set(members) & background
        size = len(term_background)
        if size < MIN_TERM or size > MAX_TERM:
            continue
        hits = len(term_background & genes)
        p = stats.hypergeom.sf(hits - 1, len(background), size, len(in_background))
        records.append({
            "term": term,
            "term_size_in_background": size,
            "hits": hits,
            "expected": len(in_background) * size / len(background),
            "p": p,
            "genes": ";".join(sorted(term_background & genes)),
        })
    table = pd.DataFrame(records)
    if table.empty:
        return table
    table["q"] = C.bh(table["p"].to_numpy())
    return table.sort_values("p").reset_index(drop=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", type=Path, default=C.DATA)
    ap.add_argument("--out", type=Path, default=C.RESULTS)
    args = ap.parse_args()

    author = pd.read_csv(args.data / "processed" / "GSE58492_author_deseq.csv", index_col=0)
    background = set(author.index[author["siL1_q"].notna() & author["siL2_q"].notna()])
    genes = induced_list(author)

    frames, summary = [], {}
    for name, filename in LIBRARIES.items():
        sets = C.read_gmt(args.data / "msigdb" / filename)
        table = ora(genes, background, sets)
        if table.empty:
            continue
        table.insert(0, "library", name)
        frames.append(table)
        summary[name] = {
            "list_size": len(genes),
            "list_in_background": len(genes & background),
            "background_size": len(background),
            "library_genes_in_background": len(background & set().union(*map(set, sets.values()))),
            "terms_tested": int(len(table)),
            "terms_q_below_0_05": int((table["q"] < 0.05).sum()),
            "top_terms": table.head(5)["term"].tolist(),
        }

    combined = pd.concat(frames, ignore_index=True)
    combined.to_csv(args.out / "pathway_ora.csv", index=False)

    ecm = combined[combined["term"].str.contains("Extracellular matrix organization",
                                                 case=False, na=False)]
    summary["extracellular_matrix_organization"] = (
        ecm[["library", "hits", "term_size_in_background", "p", "q"]]
        .to_dict(orient="records"))
    C.write_json(args.out / "summary_step3.json", summary)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    sys.exit(main())
