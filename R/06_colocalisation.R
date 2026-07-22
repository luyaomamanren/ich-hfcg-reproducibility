source("R/_setup.R")

dat <- read.delim(file.path(intermediate_dir, "harmonised_eqtl_ich_all.tsv"),
                  check.names = FALSE, stringsAsFactors = FALSE)
genes <- read.delim(file.path(project_root, "config", "genes.tsv"), stringsAsFactors = FALSE)

logsumexp <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}
logdiffexp <- function(a, b) {
  if (!is.finite(a) || b >= a) return(-Inf)
  a + log1p(-exp(b - a))
}
log_abf <- function(beta, variance, prior_variance) {
  0.5 * (log(variance / (variance + prior_variance)) +
           (beta^2 / variance) * prior_variance / (variance + prior_variance))
}

coloc_one <- function(d, p12 = 1e-5) {
  d <- d[is.finite(d$beta.exposure) & is.finite(d$se) &
           is.finite(d$beta_outcome_aligned) & is.finite(d$sebeta) &
           d$se > 0 & d$sebeta > 0, ]
  if (nrow(d) < 50L) {
    return(data.frame(nsnp = nrow(d), PP.H0 = NA, PP.H1 = NA, PP.H2 = NA,
                      PP.H3 = NA, PP.H4 = NA, top_shared_snp = NA,
                      top_shared_position_hg38 = NA, p12 = p12,
                      note = "Insufficient matched variants (<50)"))
  }
  l1 <- log_abf(d$beta.exposure, d$se^2, 0.15^2)
  l2 <- log_abf(d$beta_outcome_aligned, d$sebeta^2, 0.20^2)
  lsum1 <- logsumexp(l1); lsum2 <- logsumexp(l2); lshared <- logsumexp(l1 + l2)
  logs <- c(
    H0 = 0,
    H1 = log(1e-4) + lsum1,
    H2 = log(1e-4) + lsum2,
    H3 = log(1e-4) + log(1e-4) + logdiffexp(lsum1 + lsum2, lshared),
    H4 = log(p12) + lshared
  )
  posterior <- exp(logs - logsumexp(logs))
  top <- which.max(l1 + l2)
  data.frame(
    nsnp = nrow(d), PP.H0 = posterior["H0"], PP.H1 = posterior["H1"],
    PP.H2 = posterior["H2"], PP.H3 = posterior["H3"], PP.H4 = posterior["H4"],
    top_shared_snp = d$rsid[top], top_shared_position_hg38 = d$position[top],
    p12 = p12, note = "Approximate-Bayes-factor colocalisation", stringsAsFactors = FALSE
  )
}

priors <- c(1e-6, 1e-5, 1e-4)
results <- do.call(rbind, lapply(genes$gene, function(gene) {
  d <- dat[dat$gene == gene, ]
  ans <- do.call(rbind, lapply(priors, function(p12) coloc_one(d, p12)))
  ans$gene <- gene
  ans
}))
results$shared_given_associated <- with(results, PP.H4 / (PP.H3 + PP.H4))
write_tsv(results, file.path(results_dir, "colocalisation_ABF.tsv"))

method <- data.frame(
  parameter = c("Exposure prior SD", "Outcome log-OR prior SD", "p1", "p2", "p12 sensitivity",
                "Minimum matched variants", "Genome builds"),
  value = c("0.15", "0.20", "1e-4", "1e-4", "1e-6, 1e-5, 1e-4", "50",
            "eQTL Catalogue and FinnGen R10 both queried on GRCh38"),
  stringsAsFactors = FALSE
)
write_tsv(method, file.path(results_dir, "colocalisation_parameters.tsv"))
save_session_info("06_colocalisation")
