scripts <- c(
  "R/00_download_inputs.R",
  "R/01_gse24265_microarray.R",
  "R/02_gse163256_bulk_rnaseq.R",
  "R/03_gse266873_scrna.R",
  "R/03b_gse266873_pseudobulk.R",
  "R/04_fetch_eqtl_regions.R",
  "R/04b_fetch_finngen_regions.R",
  "R/05_mr_diagnostics.R",
  "R/06_colocalisation.R",
  "R/07_build_figure4.R",
  "R/08_export_environment.R",
  "R/10_build_supplementary_tables.R",
  "R/09_build_manifests.R"
)
for (script in scripts) {
  message("\n===== ", script, " =====")
  status <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", script))
  if (!identical(status, 0L)) stop("Pipeline stage failed: ", script)
}
