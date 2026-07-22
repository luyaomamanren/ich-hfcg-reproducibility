# Raw-data reanalysis of human intracerebral haemorrhage transcriptomic and genetic evidence

This repository is the reproducibility companion for the revised manuscript. It starts from
the public raw or minimally processed source files rather than manuscript-exported result
tables. The primary cohorts are GSE24265 (Affymetrix CEL files), GSE163256 (raw gene-count
matrices with technical replicates), and GSE266873 (10x single-cell matrices). GSE166638 is
treated as supportive single-patient context and is not a primary validation cohort.

## Reproduce the analysis

Use R 4.5.2 from the repository root:

```r
source("R/run_all.R")
```

`R/00_download_inputs.R` retrieves the public inputs and records checksums. R performs the
expression, MR, colocalisation and figure analyses. Small Python standard-library helpers are
included for resumable byte transfer and for deterministic Ensembl 1000 Genomes LD queries;
their inputs, thresholds and complete outputs are retained for audit.

The numbered scripts perform:

1. raw-data acquisition and checksum recording;
2. RMA normalization and paired/blocked limma analysis of GSE24265;
3. technical-replicate collapsing and patient-blocked edgeR/voom analysis of GSE163256;
4. sample-aware Seurat processing of GSE266873 (`R/03_*`) followed in a fresh,
   memory-isolated R process by raw-matrix pseudobulk aggregation (`R/03b_*`);
5. instrument selection, 1000 Genomes EUR LD clumping, F statistics, harmonisation, MR,
   heterogeneity, horizontal-pleiotropy and Steiger directionality analyses;
6. regional approximate-Bayes-factor colocalisation with prior sensitivity analysis;
7. manuscript-ready Figure 4 using an effect-estimate axis in panel A;
8. session, package, intermediate-result and SHA-256 manifests.

## Interpretation guardrails

- TLR4 and STAT3 are described only as *genetically prioritized candidate regulators* when
  supported by the fully rerun genetic analysis; they are not presented as validated therapies.
- HMOX1 and NLRP3 are expression- and cell-state-associated downstream candidates, not proven
  causal nodes.
- MR-Egger intercepts and Cochran heterogeneity cannot be estimated when LD clumping leaves too
  few independent instruments. Such cells are explicitly reported as not estimable rather than
  silently omitted.
- Colocalisation is used to distinguish a shared regional signal from two associations generated
  by distinct variants in linkage disequilibrium; it does not itself prove biological mediation.

## Repository map

- `R/`: complete analytical code.
- `tools/`: byte-transfer and public-LD query helpers.
- `config/`: prespecified genes and genomic regions.
- `data/intermediate/`: auditable transformations and design tables.
- `results/tables/`: complete numerical results and supplementary workbooks.
- `results/figures/`: vector, TIFF and preview outputs.
- `results/logs/`: R version and per-stage `sessionInfo()` records.
- `MANIFEST_SHA256.tsv`: size and checksum of every input, intermediate, result and script.

Raw public files and large serialized objects are intentionally excluded from Git. Their exact
URLs, sizes and SHA-256 hashes allow byte-identical reconstruction. The archival release bundle
contains code, intermediate tabular results, figures, logs and manifests.

## Public data sources

- GEO: GSE24265, GSE163256, GSE266873 and supportive GSE166638.
- eQTL Catalogue: GTEx whole-blood dataset QTD000356 (n = 670).
- ICH GWAS: FinnGen release 10 endpoint I9_ICH (4,056 cases and 371,717 controls).
- LD reference: 1000 Genomes phase 3 Europeans through the Ensembl GRCh37 REST service.

## Licence

Analysis code is released under the MIT Licence. Source datasets retain their original repository
terms; manuscript text and newly created tables/figures can be cited from the associated archive.
