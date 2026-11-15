# 1. Differential expression analysis (microarray)

Identifies genes that differ between RA and control samples.

**Pipeline** (`GSE10500/limma_GSE10500.R`, `GSE93272/limma_GSE93272.R`; identical logic, different Affymetrix chip):
read `.CEL` files → raw-data boxplot → RMA normalisation → boxplot → limma (`lmFit` → `contrasts.fit` → `eBayes`) → BH-adjusted `topTable` → fold changes (arithmetic and geometric) → gene annotation → threshold filtering → volcano plot and heatmaps.

Sample groups are defined in each `covdesc.txt` (columns `sample`, `treatment`).

## Run
```bash
cd GSE10500            # or GSE93272
mkdir data_files       # put the dataset's .CEL files here (from GEO)
Rscript limma_GSE10500.R
```

## DEG thresholds used for the reported results (thesis Table 4.4)

| Dataset | \|logFC\| | p-value | adj. p (BH) |
|---|---|---|---|
| GSE10500 | > 1 | < 0.05 | not applied (3 vs 5 samples) |
| GSE93272 | > 0.76 | < 0.05 | < 0.1 |

The scripts also write several alternative threshold sets (e.g. \|logFC\| > 1.5) used for exploration.

## Included results
- `results/GSE10500_DEGs_absLogFC1_p0.05.csv`, `..._absLogFC1.5_p0.05.csv`: filtered DEG tables (probe-set level, with symbols and expression values)
- `results/figures/`: raw vs normalised boxplots, volcano plot, heatmap, fold-change histogram

Note: `limma_GSE93272.R` was adapted from the GSE10500 script by changing the chip annotation packages and sample file.
