# Input files

Data are downloaded from the source repositories. The layout below is what the
scripts expect; `--data` points at this directory.

```
data/
  geo/
    GSE58492/GSE58492_siLSD1_results.txt.gz
    GSE30775/GSE30775_series_matrix.txt.gz
    GSE72687/GSE72687_series_matrix.txt.gz
    GSE95165/GSE95165_series_matrix.txt.gz
    GSE95165/GPL10558.annot.gz
    GSE154382/GSM*_Genes_ReadCount.txt.gz
    GPL570.annot.gz
  msigdb/
    h.all.v2024.1.Hs.symbols.gmt
    enrichment_enrichr_Reactome_2022.gmt
    enrichment_enrichr_KEGG_2021_Human.gmt
  cptac/
    protein_matrix.csv
    sample_metadata.csv
  depmap/
    Model.csv
    CRISPRGeneEffect.csv
    OmicsExpressionProteinCodingGenesTPMLogp1.csv
  singlecell/
    GSE118389.csv.gz
    GSE176078.csv.gz
  spatial/
    <section>/counts.csv.gz
    <section>/coords.csv
    <section>/patient.txt
  hypoxia/
    fpkm.csv
    design.csv
```

## GEO

```
base=https://ftp.ncbi.nlm.nih.gov/geo
curl -O $base/series/GSE58nnn/GSE58492/suppl/GSE58492_siLSD1_results.txt.gz
curl -O $base/series/GSE30nnn/GSE30775/matrix/GSE30775_series_matrix.txt.gz
curl -O $base/series/GSE72nnn/GSE72687/matrix/GSE72687_series_matrix.txt.gz
curl -O $base/series/GSE95nnn/GSE95165/matrix/GSE95165_series_matrix.txt.gz
curl -O $base/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz
curl -O $base/platforms/GPLnnn/GPL10558/annot/GPL10558.annot.gz
```

GSE154382 count files are the `GSM*_Genes_ReadCount.txt.gz` members of the
series supplementary tar archive.

## Gene sets

Hallmark from MSigDB:

```
curl -O https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs/h.all.v2024.1.Hs.symbols.gmt
```

Reactome 2022 and KEGG 2021 Human from Enrichr:

```
curl -o enrichment_enrichr_Reactome_2022.gmt \
  "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=Reactome_2022"
curl -o enrichment_enrichr_KEGG_2021_Human.gmt \
  "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=KEGG_2021_Human"
```

## CPTAC

Proteomic Data Commons, study PDC000120, release S060, version 5.4. Use the
normalised protein ratio matrix for the 122 tumours and the clinical table for
the subtype labels.

`protein_matrix.csv` is samples in rows with proteins as columns.
`sample_metadata.csv` has two columns, `sample_id` and `tnbc_group`, with the
value `TNBC` for the 28 triple-negative tumours.

## DepMap

Public 24Q4 release:

- `Model.csv`
- `CRISPRGeneEffect.csv`
- `OmicsExpressionProteinCodingGenesTPMLogp1.csv`

## Single cell

Two cohorts, one file per accession, selected with `--cohort`. Layout: one row
per cell with columns `patient`, `compartment` and gene symbols, author cell
labels kept. Values are log2(TPM + 1) for GSE118389 and log1p counts for
GSE176078; the script only ranks within a cell, so the two scales are never
compared.

- `GSE118389.csv.gz`: gene-level TPM from the series supplementary file. Keep
  the 810 author-annotated epithelial cells from six patients and set
  `compartment` to `epithelial`.
- `GSE176078.csv.gz`: cancer epithelial cells from seven IHC-defined TNBC
  donors in the CELLxGENE object `36ab5d3d-158e-4988-9700-e14fc7102fe4`. Take
  at most 500 barcodes per donor in lexicographic order; the cap is independent
  of expression and yields 3,124 cells. The CAF comparison in step 6 uses the
  same object with up to 100 barcodes per donor and compartment.

## Spatial

GSE210616, 10x Visium. Each section provides
`*_filtered_feature_bc_matrix.h5` and `*_tissue_positions_list.csv.gz`.

Build one directory per patient from the first available section:

- keep in-tissue spots with at least 200 detected genes and 500 UMIs;
- sum rows that map to the same gene symbol;
- write raw integer counts to `counts.csv.gz`, spots in rows with a `barcode`
  column, genes in columns;
- write `coords.csv` with `barcode`, `x` and `y` from `pixel_x` and `pixel_y`;
- write the patient identifier to `patient.txt`.

Keep the full feature set of the filtered matrix, 36,591 features, including
features that carry no counts in a given section. The per-spot percentile
background is every feature in `counts.csv.gz`, so the scores depend on it.

The 22 sections used here contain 29,298 spots after filtering. Counts are raw;
step 7 converts to counts per 10,000 UMIs before scoring.

## Hypoxia

`GSE108833_DifferentialExpression_FPKM.txt.gz` has `GENE_NAME` and 12 libraries.
Keep the wild-type and HIF1A/EPAS1 double-knockout columns at 20% and 1% oxygen
(eight libraries), average duplicate symbols on the FPKM scale and write
`fpkm.csv` with genes in rows. `design.csv` has `sample`, `genotype`
(`wild_type` or `double_knockout`) and `oxygen` (`20` or `1`). The ZMYND8
knockout columns are not used.
