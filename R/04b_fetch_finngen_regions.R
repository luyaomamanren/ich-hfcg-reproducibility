source("R/_setup.R")
suppressPackageStartupMessages({
  library(Rsamtools)
  library(GenomicRanges)
})

genes <- read.delim(file.path(project_root, "config", "genes.tsv"), stringsAsFactors = FALSE)
url <- paste0(
  "http://storage.googleapis.com/finngen-public-data-r10/summary_stats/",
  "finngen_R10_I9_ICH.gz"
)
columns <- c(
  "chromosome", "position", "ref", "alt", "rsids", "nearest_genes",
  "pval", "mlogp", "beta", "sebeta", "af_alt", "af_alt_cases", "af_alt_controls"
)
tabix <- TabixFile(url)
out_dir <- file.path(raw_data_dir, "genetics", "ich")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

for (i in seq_len(nrow(genes))) {
  row <- genes[i, ]
  region <- GRanges(row$chr_hg38,
                    IRanges(max(1, row$start_hg38 - 1000000), row$end_hg38 + 1000000))
  lines <- scanTabix(tabix, param = region)[[1]]
  con <- textConnection(lines)
  x <- read.delim(con, header = FALSE, col.names = columns, check.names = FALSE,
                  stringsAsFactors = FALSE)
  close(con)
  x <- x[!duplicated(paste(x$position, x$ref, x$alt, sep = ":")), ]
  path <- file.path(out_dir, paste0("finngen_R10_I9_ICH_", row$gene, "_region.tsv"))
  write_tsv(x, path)
  message(row$gene, ": ", nrow(x), " FinnGen ICH association rows")
  Sys.sleep(2)
}
save_session_info("04b_finngen_acquisition")
