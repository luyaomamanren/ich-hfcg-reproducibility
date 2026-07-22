# Reproducibility and claim-audit checklist

| Requirement | Implementation | Public output |
|---|---|---|
| Reprocess principal GEO cohorts from original deposits | GSE24265 CEL/RMA/blocked limma; GSE163256 raw counts/technical-replicate collapse/patient-blocked edgeR-voom; GSE266873 raw 10x/sample-aware Seurat and pseudobulk | `R/01_*`, `R/02_*`, `R/03_*`, complete tables and sample maps |
| Publish all analytical code | Numbered R workflow plus deterministic transfer/LD-query helpers | `R/`, `tools/`, `config/` |
| Record the computational environment | R version, package table, per-stage session information, `renv.lock` | `results/logs/`, `results/tables/R_package_versions.tsv`, `renv.lock` |
| Retain intermediate results | Design tables, replicate map, normalized-analysis inputs, harmonised variants, regional eQTL/GWAS extracts | `data/intermediate/` |
| Report MR instruments and strength | Harmonised IV table with alleles, betas, SEs, F statistic and strength flag | `results/tables/MR_instruments_harmonised.tsv` |
| Perform LD clumping | 1000 Genomes European LD, r2 < 0.001, 500 kb | `data/intermediate/mr_ld_clumping.tsv` |
| Heterogeneity and horizontal pleiotropy | Explicitly marked non-estimable because one independent IV remained per gene | `results/tables/MR_diagnostics.tsv` |
| Steiger directionality | Variance-explained comparison with documented frequency proxy | `results/tables/MR_diagnostics.tsv` |
| Multiple testing | Benjamini-Hochberg correction across four genes | `results/tables/MR_results.tsv` |
| Reduce LD-confounding risk | Approximate-Bayes-factor regional colocalisation and p12 prior sensitivity | `results/tables/colocalisation_*.tsv` |
| Correct Figure 4A | Log-odds effect-estimate axis with zero null line; odds ratios appear only after exponentiation | `results/figures/Figure_4_revised.*` and `source_data/` |
| Correct causal language | TLR4/STAT3 are genetically prioritized candidate regulators; HMOX1/NLRP3 are downstream expression/cell-state candidates | Revised manuscript and `REVISION_NOTES.md` |
| Audit file identity | Byte size and SHA-256 for input, intermediate, result and code files | `MANIFEST_SHA256.tsv` |

Blank or `NA` diagnostic cells are accompanied by an estimability note. They are not failed
computations and must not be replaced with zeroes.
