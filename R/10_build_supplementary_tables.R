source("R/_setup.R")
suppressPackageStartupMessages(library(writexl))

read_if_present <- function(path) {
  if (!file.exists(path)) return(data.frame(note = paste("Not generated:", basename(path))))
  read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
}

mr_tables <- list(
  `S5.1 MR estimates` = read_if_present(file.path(results_dir, "MR_results.tsv")),
  `S5.2 Harmonised IVs` = read_if_present(file.path(results_dir, "MR_instruments_harmonised.tsv")),
  `S5.3 LD clumping` = read_if_present(file.path(intermediate_dir, "mr_ld_clumping.tsv")),
  `S5.4 Diagnostics` = read_if_present(file.path(results_dir, "MR_diagnostics.tsv")),
  `S5.5 Provenance` = read_if_present(file.path(results_dir, "MR_data_provenance.tsv")),
  `S6.1 Colocalisation` = read_if_present(file.path(results_dir, "colocalisation_ABF.tsv")),
  `S6.2 Coloc parameters` = read_if_present(file.path(results_dir, "colocalisation_parameters.tsv"))
)
write_xlsx(mr_tables, file.path(results_dir, "TableS5_S6_MR_colocalisation_complete.xlsx"))

expression_tables <- list(
  `GSE24265 full DE` = read_if_present(file.path(results_dir, "GSE24265_DE_from_CEL.tsv")),
  `GSE24265 candidates` = subset(
    read_if_present(file.path(results_dir, "GSE24265_DE_from_CEL.tsv")),
    gene %in% c("TLR4", "STAT3", "HMOX1", "NLRP3")
  ),
  `GSE163256 full DE` = read_if_present(file.path(results_dir, "GSE163256_DE_from_counts.tsv")),
  `GSE163256 candidates` = subset(
    read_if_present(file.path(results_dir, "GSE163256_DE_from_counts.tsv")),
    gene %in% c("TLR4", "STAT3", "HMOX1", "NLRP3")
  ),
  `Cross-cohort candidates` = subset(
    read_if_present(file.path(results_dir, "cross_cohort_gene_concordance.tsv")),
    gene %in% c("TLR4", "STAT3", "HMOX1", "NLRP3")
  ),
  `GSE266873 candidates` = read_if_present(file.path(results_dir, "GSE266873_candidate_expression.tsv")),
  `GSE266873 composition` = read_if_present(file.path(results_dir, "GSE266873_cell_composition.tsv")),
  `GSE266873 sample QC` = read_if_present(file.path(results_dir, "GSE266873_sample_qc_summary.tsv")),
  `GSE266873 marker scores` = read_if_present(file.path(results_dir, "GSE266873_cluster_marker_module_scores.tsv")),
  `GSE266873 pseudobulk DE` = read_if_present(file.path(results_dir, "GSE266873_pseudobulk_DE.tsv"))
)
write_xlsx(expression_tables, file.path(results_dir, "TableS2_S4_raw_GEO_reanalysis.xlsx"))
save_session_info("10_supplementary_tables")
