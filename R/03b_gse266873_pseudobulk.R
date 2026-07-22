source("R/_setup.R")
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(edgeR)
  library(ggplot2)
})

accession <- "GSE266873"
stage_path <- file.path(intermediate_dir, "GSE266873_annotation_stage.rds")
if (!file.exists(stage_path)) stop("Run R/03_gse266873_scrna.R first")
stage <- readRDS(stage_path)
meta <- stage$metadata
sample_map <- stage$sample_design
extract_dir <- stage$extract_dir

find_one <- function(gsm, pattern) {
  hits <- list.files(extract_dir, pattern = paste0(gsm, ".*", pattern),
                     recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  if (length(hits) != 1L) stop("Expected one ", pattern, " file for ", gsm,
                               "; found ", length(hits))
  hits
}

cell_types <- sort(unique(meta$cell_type))
pb_by_sample <- lapply(sample_map$sample, function(gsm) {
  counts_i <- ReadMtx(
    mtx = find_one(gsm, "matrix\\.mtx(\\.gz)?$"),
    cells = find_one(gsm, "barcodes\\.tsv(\\.gz)?$"),
    features = find_one(gsm, "(features|genes)\\.tsv(\\.gz)?$"),
    feature.column = 2
  )
  meta_i <- meta[meta$sample == gsm, , drop = FALSE]
  raw_barcodes <- sub(paste0("^", gsm, "_"), "", rownames(meta_i))
  keep_index <- match(raw_barcodes, colnames(counts_i))
  if (anyNA(keep_index)) stop("Unmatched filtered barcodes for ", gsm)
  counts_i <- counts_i[, keep_index, drop = FALSE]
  colnames(counts_i) <- rownames(meta_i)
  aggregated <- vapply(cell_types, function(celltype) {
    Matrix::rowSums(counts_i[, meta_i$cell_type == celltype, drop = FALSE])
  }, numeric(nrow(counts_i)))
  rownames(aggregated) <- rownames(counts_i)
  colnames(aggregated) <- paste(gsm, cell_types, sep = "__")
  aggregated
})
reference_genes <- rownames(pb_by_sample[[1]])
if (!all(vapply(pb_by_sample, function(x) identical(rownames(x), reference_genes), logical(1)))) {
  stop("Feature rows differ across the nine raw 10x matrices")
}
pb <- do.call(cbind, pb_by_sample)
rm(pb_by_sample)
invisible(gc())
pb_meta <- do.call(rbind, strsplit(colnames(pb), "__", fixed = TRUE))
pb_meta <- data.frame(sample = pb_meta[, 1], cell_type = pb_meta[, 2], stringsAsFactors = FALSE)
pb_meta <- merge(pb_meta, sample_map, by = "sample", all.x = TRUE, sort = FALSE)
rownames(pb_meta) <- paste(pb_meta$sample, pb_meta$cell_type, sep = "__")
pb <- pb[, rownames(pb_meta), drop = FALSE]

run_pb_de <- function(celltype) {
  keep_cols <- pb_meta$cell_type == celltype
  group_counts <- table(pb_meta$group[keep_cols])
  if (length(group_counts) < 3L || any(group_counts < 2L)) return(NULL)
  y <- DGEList(pb[, keep_cols, drop = FALSE])
  keep_genes <- filterByExpr(y, group = pb_meta$group[keep_cols])
  if (sum(keep_genes) < 100L) return(NULL)
  y <- calcNormFactors(y[keep_genes, , keep.lib.sizes = FALSE])
  group <- factor(pb_meta$group[keep_cols], levels = c("0-6h", "6-24h", "24-48h"))
  design <- model.matrix(~0 + group)
  colnames(design) <- levels(group)
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design, robust = TRUE)
  contrasts <- list(`6-24h_vs_0-6h` = c(-1, 1, 0), `24-48h_vs_0-6h` = c(-1, 0, 1))
  do.call(rbind, lapply(names(contrasts), function(label) {
    tab <- topTags(glmQLFTest(fit, contrast = contrasts[[label]]), n = Inf, sort.by = "none")$table
    tab$gene <- rownames(tab); tab$contrast <- label; tab$cell_type <- celltype
    tab
  }))
}
de <- do.call(rbind, Filter(Negate(is.null), lapply(cell_types, run_pb_de)))
write_tsv(de, file.path(results_dir, "GSE266873_pseudobulk_DE.tsv"))

plot_data <- cbind(as.data.frame(stage$umap), meta[rownames(stage$umap), , drop = FALSE])
colnames(plot_data)[1:2] <- c("UMAP_1", "UMAP_2")
p_umap <- ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2, colour = cell_type)) +
  geom_point(size = 0.15, alpha = 0.85) +
  facet_wrap(~group, nrow = 1) +
  labs(title = "GSE266873: annotation by prespecified marker modules",
       colour = NULL, x = "UMAP 1", y = "UMAP 2") +
  theme_classic(base_size = 11) +
  theme(legend.position = "right", strip.background = element_blank(),
        strip.text = element_text(face = "bold"))
ggsave(file.path(figures_dir, "GSE266873_UMAP_celltypes.pdf"), p_umap,
       width = 10.5, height = 4.0, device = cairo_pdf)
ggsave(file.path(figures_dir, "GSE266873_UMAP_celltypes.tiff"), p_umap,
       width = 10.5, height = 4.0, dpi = 600, compression = "lzw")
ggsave(file.path(figures_dir, "GSE266873_UMAP_celltypes_preview.png"), p_umap,
       width = 10.5, height = 4.0, dpi = 180, bg = "white")

saveRDS(stage,
        file.path(objects_dir, "GSE266873_reduced_reproducibility_object.rds"),
        compress = "gzip")
save_session_info("03b_gse266873_pseudobulk")
