source("R/_setup.R")
suppressPackageStartupMessages({
  library(Rsamtools)
  library(GenomicRanges)
})

genes <- read.delim(file.path(project_root, "config", "genes.tsv"), stringsAsFactors = FALSE)
url <- paste0(
  "ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/",
  "QTS000015/QTD000356/QTD000356.all.tsv.gz"
)
columns <- c(
  "molecular_trait_id", "chromosome", "position", "ref", "alt", "variant",
  "ma_samples", "maf", "pvalue", "beta", "se", "type", "ac", "an", "r2",
  "molecular_trait_object_id", "gene_id", "median_tpm", "rsid"
)
tabix <- TabixFile(url)
dir.create(file.path(raw_data_dir, "genetics", "eqtl"), recursive = TRUE, showWarnings = FALSE)

for (i in seq_len(nrow(genes))) {
  row <- genes[i, ]
  region <- GRanges(row$chr_hg38,
                    IRanges(max(1, row$start_hg38 - 1000000), row$end_hg38 + 1000000))
  lines <- scanTabix(tabix, param = region)[[1]]
  con <- textConnection(lines)
  x <- read.delim(con, header = FALSE, col.names = columns, check.names = FALSE,
                  stringsAsFactors = FALSE)
  close(con)
  x <- x[x$gene_id == row$ensembl_gene_id & !duplicated(x$rsid), ]
  path <- file.path(raw_data_dir, "genetics", "eqtl", paste0("QTD000356_", row$gene, ".tsv"))
  write_tsv(x, path)
  message(row$gene, ": ", nrow(x), " cis-eQTL association rows")
  Sys.sleep(2)
}
save_session_info("04_eqtl_acquisition")
