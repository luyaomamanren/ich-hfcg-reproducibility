source("R/_setup.R")
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

mr <- read.delim(file.path(results_dir, "MR_results.tsv"), check.names = FALSE)
mr <- mr[mr$method %in% c("Wald ratio", "Inverse variance weighted"), ]
mr$gene <- factor(mr$gene, levels = rev(c("TLR4", "STAT3", "HMOX1", "NLRP3")))
mr$low <- mr$beta - 1.96 * mr$se
mr$high <- mr$beta + 1.96 * mr$se
mr$label <- ifelse(is.na(mr$p_FDR), "FDR NA", sprintf("FDR = %.3g", mr$p_FDR))

p_a <- ggplot(mr, aes(beta, gene)) +
  geom_vline(xintercept = 0, linewidth = 0.45, colour = "grey55", linetype = 2) +
  geom_errorbar(aes(xmin = low, xmax = high), orientation = "y", width = 0.17,
                linewidth = 0.65, colour = "#244A73") +
  geom_point(size = 2.5, shape = 21, stroke = 0.6, fill = "#D97A35", colour = "#244A73") +
  geom_text(aes(x = high, label = label), hjust = -0.08, size = 2.7) +
  scale_x_continuous(expand = expansion(mult = c(0.08, 0.32))) +
  labs(x = "Effect estimate (log odds scale)", y = NULL,
       title = "A  Mendelian randomization") +
  theme_manuscript() +
  theme(plot.title = element_text(face = "bold", size = 10),
        axis.text.y = element_text(face = "italic"))

coloc <- read.delim(file.path(results_dir, "colocalisation_ABF.tsv"), check.names = FALSE)
coloc <- coloc[coloc$p12 == 1e-5, ]
coloc$gene <- factor(coloc$gene, levels = levels(mr$gene))
write_tsv(mr, file.path(project_root, "source_data", "Figure4A_MR_effect_estimates.tsv"))
write_tsv(coloc, file.path(project_root, "source_data", "Figure4B_colocalisation.tsv"))
p_b <- ggplot(coloc, aes(PP.H4, gene)) +
  geom_vline(xintercept = 0.8, linewidth = 0.45, colour = "grey55", linetype = 2) +
  geom_segment(data = coloc[!is.na(coloc$PP.H4), ],
               aes(x = 0, xend = PP.H4, yend = gene), linewidth = 1.0, colour = "#7CA6A1") +
  geom_point(data = coloc[!is.na(coloc$PP.H4), ], size = 2.5, shape = 21,
             fill = "#D97A35", colour = "#244A73", stroke = 0.6) +
  geom_text(data = coloc[is.na(coloc$PP.H4), ], aes(x = 0.5, label = "Not estimable"),
            colour = "grey40", size = 2.8, fontface = "italic") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(x = expression("Colocalization posterior " * PP[H4]), y = NULL,
       title = "B  Shared-signal assessment") +
  theme_manuscript() +
  theme(plot.title = element_text(face = "bold", size = 10),
        axis.text.y = element_text(face = "italic"))

figure <- p_a + p_b + plot_layout(widths = c(1.15, 0.85))
ggsave(file.path(figures_dir, "Figure_4_revised.pdf"), figure,
       width = 7.2, height = 3.6, device = cairo_pdf)
ggsave(file.path(figures_dir, "Figure_4_revised.tiff"), figure,
       width = 7.2, height = 3.6, dpi = 600, compression = "lzw")
ggsave(file.path(figures_dir, "Figure_4_revised_preview.png"), figure,
       width = 7.2, height = 3.6, dpi = 180, bg = "white")

diagnostics <- read.delim(file.path(results_dir, "MR_diagnostics.tsv"), check.names = FALSE)
diagnostics$gene <- factor(diagnostics$gene, levels = c("NLRP3", "HMOX1", "STAT3", "TLR4"))
write_tsv(diagnostics, file.path(project_root, "source_data", "Supplementary_Figure_S4_MR_diagnostics.tsv"))
p_s4a <- ggplot(diagnostics, aes(gene, min_F)) +
  geom_hline(yintercept = 10, colour = "grey55", linetype = 2, linewidth = 0.45) +
  geom_segment(aes(xend = gene, y = 1, yend = min_F), linewidth = 1.1, colour = "#7CA6A1") +
  geom_point(size = 2.7, shape = 21, fill = "#D97A35", colour = "#244A73", stroke = 0.6) +
  scale_y_log10(breaks = c(1, 10, 100, 1000, 10000)) +
  labs(x = NULL, y = "Minimum F statistic (log scale)", title = "A  Instrument strength") +
  theme_manuscript() +
  theme(plot.title = element_text(face = "bold", size = 10),
        axis.text.x = element_text(face = "italic"))

diag_text <- data.frame(
  y = 4:1,
  gene = c("HMOX1", "NLRP3", "STAT3", "TLR4"),
  text = paste0("1 IV; Q: not estimable; Egger: not estimable; Steiger: correct")
)
p_s4b <- ggplot(diag_text, aes(x = 0, y = y)) +
  geom_text(aes(label = gene), hjust = 0, fontface = "italic", size = 3.0) +
  geom_text(aes(x = 0.18, label = text), hjust = 0, size = 2.75, colour = "grey25") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0.5, 4.5), clip = "off") +
  labs(title = "B  Sensitivity-test estimability") +
  theme_void(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 10, hjust = 0),
        plot.margin = margin(8, 8, 8, 8))

supp4 <- p_s4a + p_s4b + plot_layout(widths = c(0.75, 1.45))
ggsave(file.path(figures_dir, "Supplementary_Figure_S4_corrected.pdf"), supp4,
       width = 7.2, height = 3.6, device = cairo_pdf)
ggsave(file.path(figures_dir, "Supplementary_Figure_S4_corrected.tiff"), supp4,
       width = 7.2, height = 3.6, dpi = 600, compression = "lzw")
ggsave(file.path(figures_dir, "Supplementary_Figure_S4_corrected_preview.png"), supp4,
       width = 7.2, height = 3.6, dpi = 180, bg = "white")
save_session_info("07_figure4")
