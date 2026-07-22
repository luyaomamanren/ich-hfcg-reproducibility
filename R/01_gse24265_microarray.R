source("R/_setup.R")
suppressPackageStartupMessages({
  library(affy)
  library(limma)
  library(AnnotationDbi)
  library(hgu133plus2.db)
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

raw_tar <- file.path(project_root, "data", "raw", "GSE24265", "GSE24265_RAW.tar")
if (!file.exists(raw_tar)) stop("Missing GSE24265_RAW.tar; run R/00_download_inputs.R")
cel_dir <- file.path(project_root, "data", "extracted", "GSE24265", "CEL")
dir.create(cel_dir, recursive = TRUE, showWarnings = FALSE)
untar(raw_tar, exdir = cel_dir)

gz_files <- list.files(cel_dir, pattern = "\\.CEL\\.gz$", full.names = TRUE, ignore.case = TRUE)
for (gz_path in gz_files) {
  cel_path <- sub("\\.gz$", "", gz_path, ignore.case = TRUE)
  if (!file.exists(cel_path)) {
    input <- gzfile(gz_path, "rb")
    output <- file(cel_path, "wb")
    repeat {
      block <- readBin(input, "raw", n = 1024 * 1024)
      if (!length(block)) break
      writeBin(block, output)
    }
    close(input); close(output)
  }
}

cel_files <- list.files(cel_dir, pattern = "\\.CEL$", full.names = TRUE, ignore.case = TRUE)
if (length(cel_files) != 11L) stop("Expected 11 CEL files, found ", length(cel_files))
raw_affy <- ReadAffy(filenames = cel_files)
eset <- rma(raw_affy)
expr_probe <- exprs(eset)

gsm <- sub(".*(GSM[0-9]+).*", "\\1", basename(sampleNames(eset)))
sample_map <- data.frame(
  accession = c("GSM596842","GSM596843","GSM596844","GSM596845","GSM596846","GSM596847",
                "GSM596848","GSM596849","GSM596850","GSM596851","GSM596852"),
  subject = c("27","27","27","28","28","28","30","30","31","31","31"),
  region = c("PH","CG","CW","PH","CG","CW","PH","CG","PH","CG","CW")
)
sample_map$tissue <- ifelse(sample_map$region == "PH", "PH", "Contralateral")
sample_map <- sample_map[match(gsm, sample_map$accession), ]
if (anyNA(sample_map$accession)) stop("CEL-to-GEO accession mapping failed")
rownames(sample_map) <- colnames(expr_probe)
write_tsv(sample_map, file.path(project_root, "data", "intermediate", "GSE24265_sample_design.tsv"))

annotation <- AnnotationDbi::select(
  hgu133plus2.db, keys = rownames(expr_probe), keytype = "PROBEID", columns = "SYMBOL"
)
annotation <- annotation[!is.na(annotation$SYMBOL) & annotation$SYMBOL != "", ]
probe_iqr <- apply(expr_probe, 1, IQR)
annotation$iqr <- probe_iqr[annotation$PROBEID]
annotation <- annotation[order(annotation$SYMBOL, -annotation$iqr, annotation$PROBEID), ]
annotation <- annotation[!duplicated(annotation$SYMBOL), ]
expr_gene <- expr_probe[annotation$PROBEID, , drop = FALSE]
rownames(expr_gene) <- annotation$SYMBOL

tissue <- factor(sample_map$tissue, levels = c("Contralateral", "PH"))
design <- model.matrix(~ 0 + tissue)
colnames(design) <- levels(tissue)
corfit <- duplicateCorrelation(expr_gene, design, block = sample_map$subject)
fit <- lmFit(expr_gene, design, block = sample_map$subject, correlation = corfit$consensus.correlation)
fit <- contrasts.fit(fit, makeContrasts(PH_vs_Contralateral = PH - Contralateral, levels = design))
fit <- eBayes(fit, robust = TRUE, trend = TRUE)
de <- topTable(fit, coef = "PH_vs_Contralateral", number = Inf, sort.by = "P")
de$gene <- rownames(de)
de$direction <- ifelse(de$logFC > 0, "up_in_PH", "down_in_PH")
de <- de[, c("gene", setdiff(names(de), "gene"))]
write_tsv(de, file.path(project_root, "results", "tables", "GSE24265_DE_from_CEL.tsv"))
saveRDS(expr_gene, file.path(project_root, "results", "objects", "GSE24265_RMA_gene_expression.rds"))

pca <- prcomp(t(expr_gene), scale. = FALSE)
pca_df <- cbind(sample_map, PC1 = pca$x[, 1], PC2 = pca$x[, 2])
write_tsv(pca_df, file.path(project_root, "data", "intermediate", "GSE24265_PCA_scores.tsv"))
p <- ggplot(pca_df, aes(PC1, PC2, color = tissue, shape = subject, label = accession)) +
  geom_point(size = 3) + geom_text(vjust = -0.7, size = 2.4, show.legend = FALSE) +
  scale_x_continuous(expand = expansion(mult = 0.18)) +
  scale_y_continuous(expand = expansion(mult = 0.18)) +
  scale_color_manual(values = c(Contralateral = "#4C78A8", PH = "#D1495B")) +
  labs(title = "GSE24265 raw-CEL RMA quality control", color = NULL) + theme_manuscript()
ggsave(file.path(project_root, "results", "figures", "GSE24265_PCA.pdf"), p, width = 5.5, height = 4.2)
ggsave(file.path(project_root, "results", "figures", "GSE24265_PCA.tiff"), p, width = 5.5, height = 4.2,
       dpi = 600, compression = "lzw")

de$significance <- ifelse(de$adj.P.Val < 0.05, de$direction, "not_FDR_significant")
candidates <- de[de$gene %in% c("TLR4", "STAT3", "HMOX1", "NLRP3"), ]
p_volcano <- ggplot(de, aes(logFC, -log10(P.Value), colour = significance)) +
  geom_point(alpha = 0.55, size = 0.8) +
  geom_point(data = candidates, size = 2.2, shape = 21, fill = "#D97A35", colour = "#244A73") +
  geom_text(data = candidates, aes(label = gene), colour = "#202020", size = 2.6,
            vjust = -0.75, check_overlap = TRUE) +
  scale_colour_manual(values = c(up_in_PH = "#C34A36", down_in_PH = "#3C78A8",
                                 not_FDR_significant = "grey72")) +
  labs(x = expression(log[2] * " fold change"), y = expression(-log[10] * " P"),
       title = "B  Differential expression") +
  theme_manuscript() + theme(legend.position = "none")
p_pca <- p + labs(title = "A  Raw-CEL RMA quality control")
figure2 <- p_pca + p_volcano + plot_layout(widths = c(1, 1))
write_tsv(de, file.path(project_root, "source_data", "Figure2_raw_CEL_differential_expression.tsv"))
ggsave(file.path(figures_dir, "Figure_2_raw_reanalysis.pdf"), figure2,
       width = 7.2, height = 3.8, device = cairo_pdf)
ggsave(file.path(figures_dir, "Figure_2_raw_reanalysis.tiff"), figure2,
       width = 7.2, height = 3.8, dpi = 600, compression = "lzw")
ggsave(file.path(figures_dir, "Figure_2_raw_reanalysis_preview.png"), figure2,
       width = 7.2, height = 3.8, dpi = 180, bg = "white")
save_session_info("01_gse24265")
