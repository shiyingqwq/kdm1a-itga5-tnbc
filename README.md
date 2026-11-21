# KDM1A depletion and the ECM-integrin response in TNBC

Analysis code for the manuscript "KDM1A depletion increases matrix and integrin
transcripts in MDA-MB-231 cells: links to epithelial adhesion states in TNBC".

The code covers the perturbation contrasts, the ECM30 program and its
independent test, pathway over-representation, the CPTAC protein analysis,
DepMap fitness, single-cell partial correlations, spatial associations and the
oxygen experiment. Analysis stops at the result tables; figures are rendered
from them.

## Scripts

| Script | Analysis | Input |
| --- | --- | --- |
| `00_prepare_matrices.py` | Builds gene-level matrices from GEO files | GSE58492, GSE30775, GSE72687, GSE95165, GSE154382 |
| `01_differential_expression.py` | Three knockdown contrasts, shared response, ECM30/ECM29/Adaptor10/Adaptor9/Prolif20 | processed matrices |
| `02_program_scores.py` | Percentile program scores and the independent ECM30 comparison | processed matrices |
| `03_pathway_overrepresentation.py` | Hypergeometric ORA for the dual-siRNA induced list | Reactome, KEGG, Hallmark |
| `04_cptac_proteomics.py` | Leave-self-out ECM score against candidate proteins | CPTAC breast tumours |
| `05_depmap_dependency.py` | Expression against Chronos gene effect | DepMap 24Q4 |
| `06_singlecell_partial_correlation.py` | Within-patient partial correlations, ECM29 against Adaptor10 and two control pairs | GSE118389, GSE176078 |
| `07_spatial_neighborhood.py` | ECM29 against Adaptor10 across in-tissue spots; neighbouring FN1 against focal ITGA5/ITGB1 | GSE210616 sections |
| `08_hypoxia_interaction.py` | Genotype and oxygen contrasts within each level, plus their interaction | GSE108833 |

## Data

The input files are public. `data/README.md` lists the accessions, the file
each script expects, the directory layout and the download commands. Raw GEO
files go under `data/geo/`; pathway libraries under `data/msigdb/`.

Two paths can be set by environment variable:

```
KDM1A_DATA      default ./data
KDM1A_RESULTS   default ./results
```

## Running

```
python -m pip install -r requirements.txt
python scripts/00_prepare_matrices.py
python scripts/01_differential_expression.py
python scripts/02_program_scores.py
python scripts/03_pathway_overrepresentation.py
python scripts/04_cptac_proteomics.py
python scripts/05_depmap_dependency.py
python scripts/06_singlecell_partial_correlation.py --cohort GSE118389
python scripts/06_singlecell_partial_correlation.py --cohort GSE176078
python scripts/07_spatial_neighborhood.py
python scripts/08_hypoxia_interaction.py
```

Steps 0 to 3 need the GEO series files for the five knockdown datasets, two
platform annotations and three gene-set libraries. Later steps need the
corresponding cohort files listed in `data/README.md`, which also describes how
each public object is reduced to the table the script reads.

Outputs are written to `results/`:

- `contrasts.csv`, `shared_response.csv`, `ecm30_candidates.csv`
- `gene_sets.json`, `program_tests.csv`, `program_scores.csv`
- `pathway_ora.csv`
- `cptac_protein_correlations.csv`, `depmap_dependency.csv`
- `singlecell_tests_*.csv`
- `spatial_patient_associations.csv`
- `hypoxia_gene_tests.csv`, `hypoxia_hif_score.csv`

## Output of a full run

The values below come from one run of steps 0 to 8 on the public data
described in `data/README.md`.

| Quantity | Value |
| --- | --- |
| Shared induced / reduced transcripts | 182 / 101 |
| ECM30 members from 36 candidates | 30 |
| GSE30775 ECM30 test | p = 0.0037, q = 0.0334, 25 of 30 concordant |
| GSE118389 ECM29-Adaptor10, 810 cells, six patients | mean partial rho = 0.300, p = 0.0043, q = 0.0128 |
| GSE210616 ECM29-Adaptor10, 29,298 spots, 22 sections | positive in 22 of 22, median rho = 0.293 |
| ITGA5, double knockout versus wild type at 1% oxygen | log2FC = -0.840, p = 0.0143 |
| ITGA5 genotype-by-oxygen interaction | -1.105, q = 0.0639 |

## Definitions used throughout

Effect estimates are log2 fold changes, knockdown minus control. Genome-wide
Benjamini-Hochberg q values are computed within each contrast. The shared
response requires q < 0.05, at least a 1.2-fold change (absolute log2FC >=
0.263) and the same direction in all three contrasts.

A program score is the mean within-profile percentile rank of the measured
members against all measured genes in that profile; ties, including zeros, take
average ranks. Scores are computed per sample, cell or spot and never pooled
across platforms.

Significance families follow the analysis in which the test was defined. The
nine validation tests for the ECM30 and Prolif20 programs are corrected
together; pathway terms are corrected within each library; the remaining
families are noted in the relevant script.

## Reproducibility

Run on Python 3.10 or later with the versions in `requirements.txt`. The
GSE176078 cohort of step 6 starts from the CELLxGENE object; `data/README.md`
gives the cell export. The comparison background in step 1 is the set of genes
with a finite estimate in all three contrasts, on the platform annotation
releases listed in `data/README.md`.

Gene-level q values in step 8 are adjusted over the 13,241 genes that pass the
abundance filter. Pathway q values in `pathway_ora.csv` are corrected across
every eligible term within each library.

## Code and archive

The repository is at https://github.com/shiyingqwq/kdm1a-itga5-tnbc. Tagged
releases are archived on Zenodo, and the DOI for the archived version is given
on the release page.

## License

MIT. See `LICENSE`.
