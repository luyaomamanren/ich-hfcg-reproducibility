if (!requireNamespace("renv", quietly = TRUE)) stop("The renv package is required")
.libPaths(c(normalizePath("r-library", mustWork = FALSE), .libPaths()))
required <- c(
  "affy", "hgu133plus2cdf", "hgu133plus2.db", "AnnotationDbi", "org.Hs.eg.db",
  "data.table", "edgeR", "limma", "ggplot2", "Seurat", "Matrix", "patchwork",
  "GenomicRanges", "Rsamtools", "rtracklayer", "R.utils", "TwoSampleMR", "readxl",
  "writexl", "digest", "renv"
)
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Cannot snapshot missing packages: ", paste(missing, collapse = ", "))
renv::snapshot(packages = required, prompt = FALSE, force = TRUE)
