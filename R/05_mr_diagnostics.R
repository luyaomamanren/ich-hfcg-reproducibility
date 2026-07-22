source("R/_setup.R")
suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(readxl)
})

genes <- read.delim(file.path(project_root, "config", "genes.tsv"), check.names = FALSE)
genetic_dir <- file.path(project_root, "data", "raw", "genetics")

read_eqtl <- function(gene) {
  x <- read.delim(file.path(genetic_dir, "eqtl", paste0("QTD000356_", gene, ".tsv")),
                  check.names = FALSE, stringsAsFactors = FALSE)
  x <- x[x$type == "SNP" & grepl("^rs[0-9]+$", x$rsid), ]
  x <- x[!duplicated(x$rsid), ]
  x$gene <- rep(gene, nrow(x))
  x
}
eqtl <- do.call(rbind, lapply(genes$gene, read_eqtl))
write_tsv(eqtl, file.path(intermediate_dir, "QTD000356_eqtl_hg38_combined.tsv"))

read_ich <- function(gene) {
  x <- read.delim(
    file.path(genetic_dir, "ich", paste0("finngen_R10_I9_ICH_", gene, "_region.tsv")),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  x$gene <- gene
  x
}
ich <- do.call(rbind, lapply(genes$gene, read_ich))
write_tsv(ich, file.path(intermediate_dir, "finngen_R10_ICH_regional_combined.tsv"))

harmonise_gene <- function(gene) {
  e <- eqtl[eqtl$gene == gene, ]
  o <- ich[ich$gene == gene, ]
  if (!nrow(e) || !nrow(o)) return(NULL)
  m <- merge(e, o, by = c("chromosome", "position"), suffixes = c(".exposure", ".outcome"))
  m <- m[
    nchar(m$ref.exposure) == 1L & nchar(m$alt.exposure) == 1L &
      nchar(m$ref.outcome) == 1L & nchar(m$alt.outcome) == 1L,
  ]
  same <- toupper(m$alt.exposure) == toupper(m$alt.outcome) &
    toupper(m$ref.exposure) == toupper(m$ref.outcome)
  flip <- toupper(m$alt.exposure) == toupper(m$ref.outcome) &
    toupper(m$ref.exposure) == toupper(m$alt.outcome)
  m <- m[same | flip, ]
  flip <- flip[same | flip]
  m$beta_outcome_aligned <- ifelse(flip, -m$beta.outcome, m$beta.outcome)
  m$eaf_outcome_aligned <- ifelse(flip, 1 - m$af_alt, m$af_alt)
  m$effect_allele <- toupper(m$alt.exposure)
  m$other_allele <- toupper(m$ref.exposure)
  palindromic <- paste0(m$effect_allele, m$other_allele) %in% c("AT", "TA", "CG", "GC")
  m <- m[!(palindromic & m$maf > 0.42), ]
  m <- m[order(m$pvalue), ]
  m <- m[!duplicated(m$rsid), ]
  m$gene <- gene
  m
}
harmonised_all <- do.call(rbind, Filter(Negate(is.null), lapply(genes$gene, harmonise_gene)))
write_tsv(harmonised_all, file.path(intermediate_dir, "harmonised_eqtl_ich_all.tsv"))

legacy <- as.data.frame(read_excel(
  file.path(project_root, "data", "legacy_reference", "TableS5_MR_results_sensitivity.xlsx"),
  sheet = "harmonised_data"
))
legacy$gene <- sub("_expression$", "", legacy$exposure)

harmonise_legacy_gene <- function(gene) {
  e <- legacy[legacy$gene == gene & legacy$mr_keep %in% TRUE, ]
  o <- ich[ich$gene == gene, ]
  index <- vapply(e$SNP, function(snp) {
    hits <- which(vapply(strsplit(o$rsids, "[,;]"), function(z) snp %in% z, logical(1)))
    if (length(hits)) hits[1] else NA_integer_
  }, integer(1))
  keep <- !is.na(index)
  e <- e[keep, ]
  o <- o[index[keep], ]
  same <- toupper(e$effect_allele.exposure) == toupper(o$alt) &
    toupper(e$other_allele.exposure) == toupper(o$ref)
  flip <- toupper(e$effect_allele.exposure) == toupper(o$ref) &
    toupper(e$other_allele.exposure) == toupper(o$alt)
  e <- e[same | flip, ]
  o <- o[same | flip, ]
  flip <- flip[same | flip]
  if (!nrow(e)) return(NULL)
  maf_proxy <- ifelse(flip, 1 - o$af_alt, o$af_alt)
  maf_proxy <- pmin(maf_proxy, 1 - maf_proxy)
  palindromic <- paste0(
    toupper(e$effect_allele.exposure), toupper(e$other_allele.exposure)
  ) %in% c("AT", "TA", "CG", "GC")
  keep_pal <- !(palindromic & maf_proxy > 0.42)
  data.frame(
    gene = gene,
    rsid = e$SNP,
    pvalue = e$pval.exposure,
    position = o$position,
    beta.exposure = e$beta.exposure,
    se = e$se.exposure,
    samplesize.exposure = e$samplesize.exposure,
    maf = maf_proxy,
    beta_outcome_aligned = ifelse(flip, -o$beta, o$beta),
    sebeta = o$sebeta,
    eaf_outcome_aligned = ifelse(flip, 1 - o$af_alt, o$af_alt),
    effect_allele = toupper(e$effect_allele.exposure),
    other_allele = toupper(e$other_allele.exposure),
    outcome_ref = o$ref,
    outcome_alt = o$alt,
    outcome_pval = o$pval,
    frequency_note = "FinnGen EUR effect-allele frequency used as proxy for exposure MAF",
    stringsAsFactors = FALSE
  )[keep_pal, ]
}
legacy_harmonised <- do.call(
  rbind,
  Filter(Negate(is.null), lapply(genes$gene, harmonise_legacy_gene))
)
write_tsv(
  legacy_harmonised,
  file.path(intermediate_dir, "harmonised_eqtlgen_finngen_instruments.tsv")
)
candidates <- legacy_harmonised[legacy_harmonised$pvalue < 5e-8, ]
candidate_export <- candidates[, c("gene", "rsid", "pvalue", "position")]
write_tsv(candidate_export, file.path(intermediate_dir, "mr_candidates_preclump.tsv"))
clump_file <- file.path(intermediate_dir, "mr_ld_clumping.tsv")
python <- Sys.getenv(
  "PYTHON_FOR_DOWNLOADS",
  unset = "C:/Users/32574/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe"
)
status <- system2(python, c(
  file.path(project_root, "tools", "ensembl_ld_clump.py"),
  "--input", file.path(intermediate_dir, "mr_candidates_preclump.tsv"),
  "--output", clump_file, "--r2", "0.001"
))
if (status != 0) stop("LD clumping failed")
clump <- read.delim(clump_file, check.names = FALSE, stringsAsFactors = FALSE)
instruments <- merge(
  candidates,
  clump[clump$keep == "TRUE", c("gene", "rsid", "lead_rsid", "r2")],
  by = c("gene", "rsid")
)
instruments$F_statistic <- (instruments$beta.exposure / instruments$se)^2
instruments$strong_instrument <- instruments$F_statistic > 10
instruments <- instruments[instruments$strong_instrument, ]
write_tsv(instruments, file.path(results_dir, "MR_instruments_harmonised.tsv"))

run_mr <- function(d) {
  k <- nrow(d)
  if (!k) return(NULL)
  methods <- list()
  if (k == 1L) {
    b <- d$beta_outcome_aligned / d$beta.exposure
    se_ratio <- abs(d$sebeta / d$beta.exposure)
    methods[[1]] <- data.frame(
      method = "Wald ratio", nsnp = 1L, beta = b, se = se_ratio,
      pval = 2 * pnorm(-abs(b / se_ratio))
    )
  } else {
    ivw <- mr_ivw(d$beta.exposure, d$beta_outcome_aligned, d$se, d$sebeta)
    methods[[1]] <- data.frame(
      method = "Inverse variance weighted", nsnp = k,
      beta = ivw$b, se = ivw$se, pval = ivw$pval
    )
    wm <- mr_weighted_median(d$beta.exposure, d$beta_outcome_aligned, d$se, d$sebeta)
    methods[[2]] <- data.frame(
      method = "Weighted median", nsnp = k,
      beta = wm$b, se = wm$se, pval = wm$pval
    )
    if (k >= 3L) {
      egger <- mr_egger_regression(
        d$beta.exposure, d$beta_outcome_aligned, d$se, d$sebeta,
        default_parameters()
      )
      methods[[3]] <- data.frame(
        method = "MR Egger", nsnp = k,
        beta = egger$b, se = egger$se, pval = egger$pval
      )
    }
  }
  ans <- do.call(rbind, methods)
  ans$gene <- d$gene[1]
  ans$OR <- exp(ans$beta)
  ans$OR_low <- exp(ans$beta - 1.96 * ans$se)
  ans$OR_high <- exp(ans$beta + 1.96 * ans$se)
  ans
}
mr <- do.call(rbind, Filter(Negate(is.null), lapply(split(instruments, instruments$gene), run_mr)))
primary <- mr$method %in% c("Wald ratio", "Inverse variance weighted")
mr$p_FDR <- NA_real_
mr$p_FDR[primary] <- p.adjust(mr$pval[primary], method = "BH")
write_tsv(mr, file.path(results_dir, "MR_results.tsv"))

diagnostics_for_gene <- function(d) {
  k <- nrow(d)
  if (!k) return(NULL)
  ratio <- d$beta_outcome_aligned / d$beta.exposure
  ratio_se <- abs(d$sebeta / d$beta.exposure)
  w <- 1 / ratio_se^2
  ivw_beta <- sum(w * ratio) / sum(w)
  q <- if (k >= 2L) sum(w * (ratio - ivw_beta)^2) else NA_real_
  q_p <- if (k >= 2L) pchisq(q, df = k - 1, lower.tail = FALSE) else NA_real_
  egger_intercept <- egger_se <- egger_p <- NA_real_
  if (k >= 3L) {
    fit <- summary(lm(d$beta_outcome_aligned ~ d$beta.exposure, weights = 1 / d$sebeta^2))
    egger_intercept <- coef(fit)[1, 1]
    egger_se <- coef(fit)[1, 2]
    egger_p <- coef(fit)[1, 4]
  }
  r2_exp <- 2 * d$maf * (1 - d$maf) * d$beta.exposure^2
  r_out <- get_r_from_lor(
    d$beta_outcome_aligned, d$eaf_outcome_aligned,
    ncase = 4056, ncontrol = 371717, prevalence = 0.003
  )
  r2_out <- r_out^2
  data.frame(
    gene = d$gene[1], nsnp = k, mean_F = mean(d$F_statistic), min_F = min(d$F_statistic),
    cochran_Q = q, Q_df = if (k >= 2L) k - 1L else NA_integer_, Q_pval = q_p,
    egger_intercept = egger_intercept, egger_intercept_se = egger_se,
    egger_intercept_pval = egger_p,
    steiger_R2_exposure = sum(r2_exp), steiger_R2_outcome = sum(r2_out),
    steiger_correct_direction = sum(r2_exp) > sum(r2_out),
    heterogeneity_note = if (k < 2L) "Not estimable: one independent IV" else "Estimated",
    pleiotropy_note = if (k < 3L) {
      "MR-Egger intercept not estimable: fewer than three IVs"
    } else {
      "Estimated"
    },
    stringsAsFactors = FALSE
  )
}
diagnostics <- do.call(
  rbind,
  Filter(Negate(is.null), lapply(split(instruments, instruments$gene), diagnostics_for_gene))
)
write_tsv(diagnostics, file.path(results_dir, "MR_diagnostics.tsv"))

provenance <- data.frame(
  item = c(
    "MR exposure", "Colocalisation exposure", "Outcome", "Genome build", "LD reference", "Clumping threshold",
    "Exposure threshold", "Outcome cases", "Outcome controls",
    "Steiger assumed population prevalence"
  ),
  value = c(
    "eQTLGen whole-blood cis-eQTL instrument rows; per-SNP n=19,772-31,684",
    "GTEx whole blood eQTL Catalogue QTD000356; n=670",
    "FinnGen release 10 I9_ICH",
    "GRCh38 for exposure and outcome; Ensembl GRCh37 endpoint used only for 1000G LD by rsID",
    "1000 Genomes phase 3 EUR via Ensembl REST API",
    "r2 < 0.001, 500-kb maximum public-API window within each cis region",
    "P < 5e-8", "4056", "371717", "0.003"
  ),
  stringsAsFactors = FALSE
)
write_tsv(provenance, file.path(results_dir, "MR_data_provenance.tsv"))
save_session_info("05_mr")
