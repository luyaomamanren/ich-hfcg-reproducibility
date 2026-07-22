source("R/_setup.R")
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(edgeR)
  library(ggplot2)
})

accession <- "GSE266873"
raw_dir <- file.path(project_root, "data", "raw", accession)
extract_dir <- file.path(project_root, "data", "extracted", accession)
dir.create(extract_dir, recursive = TRUE, showWarnings = FALSE)
tar_path <- file.path(raw_dir, paste0(accession, "_RAW.tar"))
direct_dir <- file.path(raw_dir, "10x")
direct_files <- list.files(direct_dir, recursive = TRUE, full.names = TRUE)
if (sum(grepl("matrix\\.mtx\\.gz$", direct_files)) == 9L) {
  extract_dir <- direct_dir
} else {
  if (!file.exists(tar_path)) stop("Missing ", tar_path, "; run R/00_download_inputs.R")
  utils::untar(tar_path, exdir = extract_dir)
}

sample_map <- data.frame(
  sample = paste0("GSM", 8255340:8255348),
  group = rep(c("0-6h", "6-24h", "24-48h"), each = 3),
  hour = c(2, 3, 2, 12, 12, 12, 26, 28, 32),
  replicate = rep(1:3, times = 3),
  stringsAsFactors = FALSE
)
write_tsv(sample_map, file.path(intermediate_dir, "GSE266873_sample_design.tsv"))

find_one <- function(gsm, pattern) {
  hits <- list.files(extract_dir, pattern = paste0(gsm, ".*", pattern),
                     recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  if (length(hits) != 1L) stop("Expected one ", pattern, " file for ", gsm,
                               "; found ", length(hits))
  hits
}

objects <- lapply(seq_len(nrow(sample_map)), function(i) {
  gsm <- sample_map$sample[i]
  counts <- ReadMtx(
    mtx = find_one(gsm, "matrix\\.mtx(\\.gz)?$"),
    cells = find_one(gsm, "barcodes\\.tsv(\\.gz)?$"),
    features = find_one(gsm, "(features|genes)\\.tsv(\\.gz)?$"),
    feature.column = 2
  )
  obj <- CreateSeuratObject(counts, project = gsm, min.cells = 3, min.features = 200)
  obj$sample <- gsm
  obj$group <- sample_map$group[i]
  obj$hour <- sample_map$hour[i]
  obj$replicate <- sample_map$replicate[i]
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj
})
names(objects) <- sample_map$sample
obj <- merge(objects[[1]], y = objects[-1], add.cell.ids = names(objects))
obj <- JoinLayers(obj)
rm(objects)
invisible(gc())
obj$group <- factor(obj$group, levels = c("0-6h", "6-24h", "24-48h"))

qc_before <- obj@meta.data
qc_before$barcode <- rownames(qc_before)
write_tsv(qc_before, file.path(results_dir, "GSE266873_cell_qc_before.tsv"))
n_before <- table(obj$sample)
obj <- subset(obj, subset = nFeature_RNA >= 200 & nFeature_RNA <= 6000 &
                nCount_RNA >= 500 & percent.mt < 20)
invisible(gc())
n_after <- table(obj$sample)
qc_summary <- data.frame(
  sample = names(n_before),
  cells_before = as.integer(n_before),
  cells_after = as.integer(n_after[names(n_before)]),
  retention_fraction = as.integer(n_after[names(n_before)]) / as.integer(n_before)
)
write_tsv(qc_summary, file.path(results_dir, "GSE266873_sample_qc_summary.tsv"))
obj <- NormalizeData(obj, verbose = FALSE)
obj <- FindVariableFeatures(obj, nfeatures = 3000, verbose = FALSE)
obj <- ScaleData(obj, features = VariableFeatures(obj), vars.to.regress = "percent.mt", verbose = FALSE)
obj <- RunPCA(obj, npcs = 30, verbose = FALSE)
obj <- FindNeighbors(obj, dims = 1:20, verbose = FALSE)
obj <- FindClusters(obj, resolution = 0.5, random.seed = 20260722, verbose = FALSE)
obj <- RunUMAP(obj, dims = 1:20, seed.use = 20260722, verbose = FALSE)

umap_table <- tibble::rownames_to_column(as.data.frame(Embeddings(obj, "umap")), "barcode")
umap_meta <- tibble::rownames_to_column(obj[[]], "barcode")
umap_table <- merge(
  umap_table,
  umap_meta[, c("barcode", "sample", "group", "hour", "seurat_clusters")],
  by = "barcode",
  all.x = TRUE,
  sort = FALSE
)
write_tsv(umap_table, file.path(intermediate_dir, "GSE266873_UMAP_coordinates.tsv"))

marker_sets <- list(
  `Microglia/macrophage` = c("C1QA", "C1QB", "AIF1", "APOE"),
  Monocyte = c("FCN1", "CTSS", "LILRB1", "CTSD"),
  Neutrophil = c("FCGR3B", "CSF3R", "S100A8", "S100A9"),
  `T/NK` = c("CD3D", "CD3E", "TRBC1", "NKG7"),
  B = c("MS4A1", "CD79A", "CD37", "CD74"),
  Platelet = c("PPBP", "PF4", "GP9"),
  Erythroid = c("HBB", "HBA1", "HBA2"),
  Endothelial = c("PECAM1", "VWF", "KDR")
)
avg <- AverageExpression(
  obj,
  assays = "RNA",
  features = unique(unlist(marker_sets)),
  layer = "data",
  group.by = "seurat_clusters",
  verbose = FALSE
)$RNA
cluster_ids <- colnames(avg)
scores <- sapply(marker_sets, function(genes) {
  genes <- intersect(genes, rownames(avg))
  if (!length(genes)) return(rep(NA_real_, ncol(avg)))
  colMeans(avg[genes, , drop = FALSE])
})
score_table <- data.frame(cluster = cluster_ids, scores, check.names = FALSE)
write_tsv(score_table, file.path(results_dir, "GSE266873_cluster_marker_module_scores.tsv"))
labels <- colnames(scores)[max.col(scores, ties.method = "first")]
names(labels) <- sub("^g(?=[0-9])", "", cluster_ids, perl = TRUE)
obj$cell_type <- unname(labels[as.character(obj$seurat_clusters)])
if (anyNA(obj$cell_type)) stop("Cluster-to-cell-type mapping produced missing labels")

meta <- obj@meta.data
meta$barcode <- rownames(meta)
write_tsv(meta, file.path(results_dir, "GSE266873_cell_metadata.tsv"))
composition <- as.data.frame(prop.table(table(meta$sample, meta$cell_type), margin = 1))
colnames(composition) <- c("sample", "cell_type", "proportion")
composition <- merge(composition, sample_map, by = "sample", all.x = TRUE)
write_tsv(composition, file.path(results_dir, "GSE266873_cell_composition.tsv"))

candidate_genes <- intersect(c("TLR4", "STAT3", "HMOX1", "NLRP3"), rownames(obj))
candidate_summary <- aggregate(
  FetchData(obj, vars = candidate_genes),
  by = list(sample = obj$sample, group = obj$group, hour = obj$hour, cell_type = obj$cell_type),
  FUN = mean
)
write_tsv(candidate_summary, file.path(results_dir, "GSE266873_candidate_expression.tsv"))

pca_embedding <- Embeddings(obj, "pca")
umap_embedding <- Embeddings(obj, "umap")
annotation_stage <- list(
  metadata = meta,
  umap = umap_embedding,
  pca = pca_embedding,
  sample_design = sample_map,
  cluster_annotation = labels,
  marker_sets = marker_sets,
  extract_dir = extract_dir
)
saveRDS(annotation_stage,
        file.path(intermediate_dir, "GSE266873_annotation_stage.rds"),
        compress = "gzip")
save_session_info("03_gse266873_annotation")
message("Annotation stage complete; run R/03b_gse266873_pseudobulk.R in a fresh R process.")
quit(save = "no", status = 0)
rm(obj)
invisible(gc())

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
de <- do.call(rbind, Filter(Negate(is.null), lapply(sort(unique(pb_meta$cell_type)), run_pb_de)))
write_tsv(de, file.path(results_dir, "GSE266873_pseudobulk_DE.tsv"))

plot_data <- cbind(as.data.frame(umap_embedding), meta[rownames(umap_embedding), , drop = FALSE])
p_umap <- ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2, colour = cell_type)) +
  geom_point(size = 0.15, alpha = 0.85) +
  facet_wrap(~group, nrow = 1) +
  labs(
    title = "GSE266873: annotation by prespecified marker modules",
    colour = NULL,
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  theme_classic(base_size = 11) +
  theme(legend.position = "right", strip.background = element_blank(), strip.text = element_text(face = "bold"))
ggsave(file.path(figures_dir, "GSE266873_UMAP_celltypes.pdf"), p_umap,
       width = 10.5, height = 4.0, device = cairo_pdf)
ggsave(file.path(figures_dir, "GSE266873_UMAP_celltypes.tiff"), p_umap,
       width = 10.5, height = 4.0, dpi = 600, compression = "lzw")
ggsave(file.path(figures_dir, "GSE266873_UMAP_celltypes_preview.png"), p_umap,
       width = 10.5, height = 4.0, dpi = 180, bg = "white")
reproducible_object <- list(
  metadata = meta,
  umap = umap_embedding,
  pca = pca_embedding,
  sample_design = sample_map,
  cluster_annotation = labels,
  marker_sets = marker_sets
)
saveRDS(reproducible_object,
        file.path(objects_dir, "GSE266873_reduced_reproducibility_object.rds"),
        compress = "gzip")
save_session_info("03_gse266873")
