options(stringsAsFactors = FALSE, scipen = 999)
set.seed(20260722)

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(project_root, "config", "genes.tsv"))) {
  stop("Run scripts from the repository root.")
}

local_library <- file.path(project_root, "r-library")
if (dir.exists(local_library)) .libPaths(c(local_library, .libPaths()))

dirs <- c(
  "data/raw", "data/extracted", "data/intermediate", "data/legacy_reference",
  "results/tables", "results/figures", "results/objects", "results/logs",
  "source_data"
)
invisible(lapply(file.path(project_root, dirs), dir.create, recursive = TRUE, showWarnings = FALSE))
raw_data_dir <- file.path(project_root, "data", "raw")
intermediate_dir <- file.path(project_root, "data", "intermediate")
results_dir <- file.path(project_root, "results", "tables")
figures_dir <- file.path(project_root, "results", "figures")
objects_dir <- file.path(project_root, "results", "objects")
logs_dir <- file.path(project_root, "results", "logs")

write_tsv <- function(x, path) {
  data.table::fwrite(x, path, sep = "\t", na = "NA", quote = FALSE)
}

save_session_info <- function(label) {
  path <- file.path(project_root, "results", "logs", paste0("sessionInfo_", label, ".txt"))
  capture.output(sessionInfo(), file = path)
  invisible(path)
}

assert_columns <- function(x, required, label = deparse(substitute(x))) {
  missing <- setdiff(required, names(x))
  if (length(missing)) stop(label, " lacks columns: ", paste(missing, collapse = ", "))
}

parse_geo_soft_titles <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con))
  accession <- character()
  title <- character()
  current <- NA_character_
  repeat {
    lines <- readLines(con, n = 10000, warn = FALSE)
    if (!length(lines)) break
    for (line in lines) {
      if (grepl("^\\^SAMPLE = ", line)) current <- sub("^\\^SAMPLE = ", "", line)
      if (grepl("^!Sample_title = ", line) && !is.na(current)) {
        accession <- c(accession, current)
        title <- c(title, sub("^!Sample_title = ", "", line))
      }
    }
  }
  unique(data.frame(accession = accession, title = title))
}

theme_manuscript <- function(base_size = 9) {
  ggplot2::theme_classic(base_size = base_size, base_family = "sans") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0),
      axis.title = ggplot2::element_text(face = "bold"),
      legend.title = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold")
    )
}
