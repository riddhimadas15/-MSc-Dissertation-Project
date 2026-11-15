# 4. From formulations to bioactives to targets

Builds the "drug side" of the analysis.

1. **Formulations and herbs**: RA (*Amavata*) formulations and their herbs were taken from the Dravya database and the literature. Compiled in `data/formulations.csv` (formulation → herb → bioactive).
2. **Bioactives**: identified with IMPPAT and CTD; drug-likeness/ADME checked with SwissADME (Lipinski's rule of five). These are web-tool steps, documented in thesis §3.2.2–3.2.3.
3. **Targets**: known targets from PubChem and predicted targets from SwissTargetPrediction (SMILES from PubChem). Notes on each herb, bioactive and its targets are in `data/herb_bioactive_target_notes.txt`.
4. **Parsing**: `extract_targets_from_notes.py` converts those notes into a tidy bioactive–gene table (`data/drug_gene_targets.csv`, an early intermediate). The curated bioactive → target table used by the RAveda tool is `data/ayurvedic_targets.csv`.
5. **RA relevance**: `pubchem_swisstarget/*.R` intersect PubChem and SwissTargetPrediction targets with the RA gene list (`genes_RA.csv`).

```bash
python extract_targets_from_notes.py   # data/herb_bioactive_target_notes.txt -> data/drug_gene_targets.csv (overwrites the shipped file)
```
