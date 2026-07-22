source("R/_setup.R")
suppressPackageStartupMessages(library(digest))

inventory <- function(root, label) {
  files <- list.files(root, recursive = TRUE, full.names = TRUE, all.files = TRUE)
  files <- files[file.info(files)$isdir %in% FALSE]
  normalized <- normalizePath(files, winslash = "/", mustWork = FALSE)
  files <- files[!grepl("\\.parts/|\\.part$|/test_(ftp_)?|hg38ToHg19\\.over\\.chain|GSE266873/10x/|/results/objects/|/tools/node_modules/|\\.inspect\\.ndjson$", normalized)]
  data.frame(
    collection = label,
    path = substring(normalizePath(files, winslash = "/"), nchar(project_root) + 2L),
    bytes = file.info(files)$size,
    sha256 = vapply(files, digest, character(1), algo = "sha256", file = TRUE),
    stringsAsFactors = FALSE
  )
}

manifest <- do.call(rbind, list(
  inventory(file.path(project_root, "data", "raw"), "raw_input"),
  inventory(file.path(project_root, "data", "intermediate"), "intermediate"),
  inventory(file.path(project_root, "results"), "result"),
  inventory(file.path(project_root, "R"), "R_code"),
  inventory(file.path(project_root, "tools"), "transfer_helper")
))
write_tsv(manifest, file.path(project_root, "MANIFEST_SHA256.tsv"))
save_session_info("09_manifest")
