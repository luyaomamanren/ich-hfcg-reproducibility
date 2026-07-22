# Result-table dictionary

- `GSE24265_DE_from_CEL.tsv`: complete gene-level patient-blocked limma output after raw-CEL RMA.
- `GSE163256_DE_from_counts.tsv`: complete gene-level patient-blocked edgeR/voom output after technical-replicate collapse.
- `cross_cohort_gene_concordance.tsv`: gene-matched bulk estimates and direction comparison.
- `GSE266873_cell_qc_before.tsv`: cell-level metrics before prespecified filtering.
- `GSE266873_sample_qc_summary.tsv`: per-sample cell retention after filtering.
- `GSE266873_cluster_marker_module_scores.tsv`: auditable scores used to assign major cell classes.
- `GSE266873_cell_metadata.tsv`: cell barcode, sample, time group, cluster and assigned class.
- `GSE266873_cell_composition.tsv`: sample-level cell-class proportions.
- `GSE266873_candidate_expression.tsv`: sample-by-time-by-cell-class mean normalized expression.
- `GSE266873_pseudobulk_DE.tsv`: complete sample-level pseudobulk quasi-likelihood contrasts.
- `MR_instruments_harmonised.tsv`: retained independent instruments with alleles, effects and F statistics.
- `MR_results.tsv`: Wald-ratio estimates, confidence intervals, odds ratios and four-gene FDR.
- `MR_diagnostics.tsv`: instrument count, strength, Steiger results and estimability notes.
- `MR_data_provenance.tsv`: exposure, outcome, LD-reference and assumption record.
- `colocalisation_ABF.tsv`: posterior probabilities for all prespecified p12 values.
- `colocalisation_parameters.tsv`: prior parameters, matching threshold and genome builds.
- `R_package_versions.tsv`: installed package versions used in the recorded environment.

`NA` in the MR diagnostics is intentional when a statistic cannot be defined with one
independent instrument; the corresponding note column explains the reason.
