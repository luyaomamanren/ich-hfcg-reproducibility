# Reanalysis decisions and claim corrections

## Raw-data scope

- **GSE24265:** all 11 CEL files were RMA-normalized. Perihematomal tissue was compared with
  contralateral grey/white matter using limma with patient blocking.
- **GSE163256:** both raw count matrices were used. The 399 sequencing columns were mapped back
  to GEO records through `Sample_description`, technical replicates were summed, and 274
  biological sample/cell-type/day combinations were retained. Hematoma was compared with ICH
  blood using a patient-blocked edgeR/voom model controlling for cell type and day.
- **GSE266873:** the nine raw 10x matrices are processed as nine biological samples in three
  post-ICH time groups. Inference is performed at sample level by pseudobulk rather than treating
  cells as independent replicates.
- **GSE166638:** single-patient longitudinal data are supportive only and are not a principal
  replication cohort.

## Expression evidence

In GSE24265, HMOX1, STAT3 and NLRP3 were nominally higher in perihematomal tissue but none passed
genome-wide FDR correction; TLR4 was unchanged. In GSE163256, HMOX1, STAT3 and NLRP3 were higher
in hematoma-derived myeloid cells, whereas TLR4 was modestly lower. This is expression evidence,
not proof of upstream causality.

## MR audit and rerun

The original 16-20 SNP sets per gene were not independent after formal 1000 Genomes EUR
clumping. At r2 < 0.001, one instrument remained for each gene. Consequently:

- causal estimates are Wald ratios, not multi-instrument IVW estimates;
- Cochran heterogeneity and MR-Egger intercept tests are mathematically not estimable;
- leave-one-out analysis is not meaningful with one instrument;
- all four retained instruments had F > 10;
- Steiger comparisons supported the exposure-to-outcome direction under the stated assumptions;
- no gene survived Benjamini-Hochberg correction in FinnGen R10 I9_ICH.

Primary corrected estimates were: HMOX1 OR 1.188 (95% CI 0.858-1.644; FDR 0.598), NLRP3 OR
1.019 (0.790-1.313; FDR 0.887), STAT3 OR 0.644 (0.403-1.030; FDR 0.265), and TLR4 OR 0.985
(0.865-1.122; FDR 0.887).

## Colocalisation

At the default shared-signal prior p12 = 1e-5, posterior support for a shared eQTL/ICH signal was
low for NLRP3 (PP.H4 = 0.0057) and TLR4 (PP.H4 = 0.0115). GTEx whole blood QTD000356 contained
no regional STAT3 or HMOX1 eQTL rows, so colocalisation for these genes is reported as not
estimable rather than imputed or fabricated.

## Required terminology

- TLR4 and STAT3 are referred to as **prespecified genetically prioritized candidate
  regulators**. The corrected analysis did not confirm them as causal or therapeutic targets.
- HMOX1 and NLRP3 are **expression- and cell-state-associated downstream candidates**, not
  established causal nodes.
- Figure 4A reports **Effect estimate (log odds scale)**, so negative beta values are valid. Any
  odds-ratio presentation is exponentiated and therefore strictly positive.
