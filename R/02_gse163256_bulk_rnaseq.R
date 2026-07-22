source("R/_setup.R")
suppressPackageStartupMessages({library(data.table); library(edgeR); library(limma); library(ggplot2); library(patchwork)})

raw_dir <- file.path(project_root, "data", "raw", "GSE163256")
soft_path <- file.path(raw_dir, "GSE163256_family.soft.gz")
files <- c(
  Monocyte = file.path(raw_dir, "GSE163256_monos_counts.csv.gz"),
  Neutrophil = file.path(raw_dir, "GSE163256_neuts_counts.csv.gz")
)
if (!all(file.exists(files)) || !file.exists(soft_path)) stop("Missing GSE163256 inputs")
titles <- parse_geo_soft_titles(soft_path)

parse_soft_internal_ids <- function(path) {
  lines <- readLines(gzfile(path), warn = FALSE)
  sample_starts <- grep("^\\^SAMPLE = ", lines)
  sample_ends <- c(sample_starts[-1] - 1L, length(lines))
  out <- lapply(seq_along(sample_starts), function(i) {
    block <- lines[sample_starts[i]:sample_ends[i]]
    accession <- sub("^\\^SAMPLE = ", "", block[1])
    title_line <- grep("^!Sample_title = ", block, value = TRUE)[1]
    descriptions <- sub("^!Sample_description = ", "",
                        grep("^!Sample_description = ", block, value = TRUE))
    internal <- descriptions[grepl("^(ICH|Control|Healthy|Batch)", descriptions,
                                   ignore.case = TRUE)][1]
    data.frame(accession = accession,
               title = sub("^!Sample_title = ", "", title_line),
               internal_id = internal, stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}
sample_lookup <- parse_soft_internal_ids(soft_path)

read_count_file <- function(path, cell_type) {
  x <- fread(path)
  gene_col <- names(x)[1]
  genes <- x[[gene_col]]
  x[[gene_col]] <- NULL
  mat <- as.matrix(x)
  storage.mode(mat) <- "integer"
  rownames(mat) <- genes
  internal_ids <- colnames(mat)
  meta <- sample_lookup[match(internal_ids, sample_lookup$internal_id), ]
  if (anyNA(meta$title)) stop("Could not match all count columns to GEO sample titles in ", path)
  meta$cell_type <- cell_type
  list(counts = mat, meta = meta)
}

parts <- Map(read_count_file, files, names(files))
common_genes <- Reduce(intersect, lapply(parts, function(z) rownames(z$counts)))
counts <- do.call(cbind, lapply(parts, function(z) z$counts[common_genes, , drop = FALSE]))
meta <- do.call(rbind, lapply(parts, `[[`, "meta"))

meta$source <- ifelse(grepl("from hematoma", meta$title, ignore.case = TRUE), "Hematoma",
                      ifelse(grepl("from HC blood", meta$title, ignore.case = TRUE), "HealthyBlood", "ICHBlood"))
meta$patient <- sub(".*patient number ([^ ]+) on day.*", "\\1", meta$title, ignore.case = TRUE)
meta$day <- as.integer(sub(".* on day ([0-9]+).*", "\\1", meta$title, ignore.case = TRUE))
meta$replicate <- as.integer(sub(".*replicate ([0-9]+).*", "\\1", meta$title, ignore.case = TRUE))
meta$bio_id <- paste(meta$cell_type, meta$source, meta$patient, meta$day, sep = "|")
if (anyNA(meta$day)) stop("Failed to parse day from GEO titles")
rownames(meta) <- colnames(counts)

collapsed <- rowsum(t(counts), group = meta$bio_id, reorder = FALSE)
collapsed <- t(collapsed)
bio_meta <- meta[match(colnames(collapsed), meta$bio_id), c("bio_id","cell_type","source","patient","day")]
rownames(bio_meta) <- bio_meta$bio_id
write_tsv(cbind(data.frame(sample = rownames(meta)), meta), file.path(project_root, "data", "intermediate", "GSE163256_technical_replicate_map.tsv"))
write_tsv(bio_meta, file.path(project_root, "data", "intermediate", "GSE163256_biological_sample_design.tsv"))

keep_samples <- bio_meta$source %in% c("Hematoma", "ICHBlood")
y <- DGEList(collapsed[, keep_samples, drop = FALSE])
m <- bio_meta[keep_samples, , drop = FALSE]
m$source <- factor(m$source, levels = c("ICHBlood", "Hematoma"))
m$cell_type <- factor(m$cell_type)
m$day_factor <- factor(m$day)
keep_genes <- filterByExpr(y, group = interaction(m$source, m$cell_type))
y <- calcNormFactors(y[keep_genes, , keep.lib.sizes = FALSE], method = "TMM")
design <- model.matrix(~ 0 + source + cell_type + day_factor, data = m)
v0 <- voom(y, design, plot = FALSE)
corfit <- duplicateCorrelation(v0, design, block = m$patient)
v <- voom(y, design, plot = FALSE, block = m$patient, correlation = corfit$consensus.correlation)
fit <- lmFit(v, design, block = m$patient, correlation = corfit$consensus.correlation)
contrast <- rep(0, ncol(design)); names(contrast) <- colnames(design)
contrast["sourceHematoma"] <- 1; contrast["sourceICHBlood"] <- -1
fit <- contrasts.fit(fit, contrasts = matrix(contrast, ncol = 1, dimnames = list(names(contrast), "Hematoma_vs_ICHBlood")))
fit <- eBayes(fit, robust = TRUE, trend = TRUE)
de <- topTable(fit, coef = 1, number = Inf, sort.by = "P")
de$gene <- rownames(de); de$direction <- ifelse(de$logFC > 0, "up_in_hematoma", "down_in_hematoma")
de <- de[, c("gene", setdiff(names(de), "gene"))]
write_tsv(de, file.path(project_root, "results", "tables", "GSE163256_DE_from_counts.tsv"))
saveRDS(v$E, file.path(project_root, "results", "objects", "GSE163256_voom_expression.rds"))

discovery_path <- file.path(project_root, "results", "tables", "GSE24265_DE_from_CEL.tsv")
if (file.exists(discovery_path)) {
  discovery <- fread(discovery_path)
  joined <- merge(discovery[, .(gene, logFC_discovery = logFC, FDR_discovery = adj.P.Val)],
                  de[, c("gene", "logFC", "adj.P.Val")], by = "gene")
  setnames(joined, c("logFC", "adj.P.Val"), c("logFC_validation", "FDR_validation"))
  joined$concordant <- sign(joined$logFC_discovery) == sign(joined$logFC_validation)
  joined$max_FDR_across_cohorts <- pmax(joined$FDR_discovery, joined$FDR_validation)
  write_tsv(joined, file.path(project_root, "results", "tables", "cross_cohort_gene_concordance.tsv"))

  core <- de[de$gene %in% c("TLR4", "STAT3", "HMOX1", "NLRP3"), ]
  core$gene <- factor(core$gene, levels = rev(c("TLR4", "STAT3", "HMOX1", "NLRP3")))
  core$se_from_moderated_t <- abs(core$logFC / core$t)
  core$low <- core$logFC - 1.96 * core$se_from_moderated_t
  core$high <- core$logFC + 1.96 * core$se_from_moderated_t
  core$label <- sprintf("FDR = %.3g", core$adj.P.Val)
  p_core <- ggplot(core, aes(logFC, gene)) +
    geom_vline(xintercept = 0, colour = "grey55", linetype = 2, linewidth = 0.45) +
    geom_errorbar(aes(xmin = low, xmax = high), orientation = "y", width = 0.16,
                  linewidth = 0.65, colour = "#244A73") +
    geom_point(size = 2.5, shape = 21, fill = "#D97A35", colour = "#244A73") +
    geom_text(aes(x = high, label = label), hjust = -0.08, size = 2.5) +
    scale_x_continuous(expand = expansion(mult = c(0.08, 0.70))) +
    labs(x = expression(log[2] * " fold change"), y = NULL,
         title = "A  Raw-count validation") +
    theme_manuscript() + theme(axis.text.y = element_text(face = "italic"))

  core_joined <- joined[joined$gene %in% c("TLR4", "STAT3", "HMOX1", "NLRP3"), ]
  rho <- suppressWarnings(cor(joined$logFC_discovery, joined$logFC_validation,
                              method = "spearman", use = "complete.obs"))
  p_concordance <- ggplot(joined, aes(logFC_discovery, logFC_validation)) +
    geom_hline(yintercept = 0, colour = "grey82", linewidth = 0.35) +
    geom_vline(xintercept = 0, colour = "grey82", linewidth = 0.35) +
    geom_point(alpha = 0.18, size = 0.65, colour = "#52789E") +
    geom_point(data = core_joined, size = 2.2, shape = 21,
               fill = "#D97A35", colour = "#244A73") +
    geom_text(data = core_joined, aes(label = gene), size = 2.5, vjust = -0.75,
              check_overlap = TRUE) +
    annotate("text", x = -Inf, y = Inf, label = sprintf("Spearman rho = %.2f", rho),
             hjust = -0.05, vjust = 1.3, size = 2.7) +
    labs(x = "GSE24265 log2 fold change", y = "GSE163256 log2 fold change",
         title = "B  Cross-cohort directionality") +
    theme_manuscript()
  figure3 <- p_core + p_concordance + plot_layout(widths = c(1.05, 0.95))
  write_tsv(core, file.path(project_root, "source_data", "Figure3A_GSE163256_candidates.tsv"))
  write_tsv(joined, file.path(project_root, "source_data", "Figure3B_cross_cohort_concordance.tsv"))
  ggsave(file.path(figures_dir, "Figure_3_raw_reanalysis.pdf"), figure3,
         width = 7.2, height = 3.8, device = cairo_pdf)
  ggsave(file.path(figures_dir, "Figure_3_raw_reanalysis.tiff"), figure3,
         width = 7.2, height = 3.8, dpi = 600, compression = "lzw")
  ggsave(file.path(figures_dir, "Figure_3_raw_reanalysis_preview.png"), figure3,
         width = 7.2, height = 3.8, dpi = 180, bg = "white")
}
save_session_info("02_gse163256")
