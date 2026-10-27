# Multi-Algorithm Immune Deconvolution and LASSO Machine Learning in Oral Squamous Cell Carcinoma (OSCC)

[![R Version](https://img.shields.io/badge/R-4.5.3-blue.svg)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![OSF Preregistration](https://img.shields.io/badge/OSF-Preregistered-brightgreen.svg)](https://osf.io/7nr2z/)

## Overview
This repository contains the complete bioinformatics analytical pipeline, statistical survival models, and machine-learning code accompanying the manuscript submitted to ... :

> **"Multi-Algorithm Immune Deconvolution and LASSO Machine Learning Reveal Macrophage Exhaustion and a 12-Gene Prognostic Signature in Oral Squamous Cell Carcinoma"**

Challenging traditional immuno-oncology paradigms, this study integrates three orthogonal transcriptomic deconvolution architectures (**CIBERSORT**, **quanTIseq**, and **xCell**) across anatomically restricted primary Oral Squamous Cell Carcinoma (OSCC) tumors from The Cancer Genome Atlas (`TCGA-HNSC`, $N = 220$; $n = 219$ with complete survival follow-up). We demonstrate that physical M2 macrophage depletion—reflecting metabolic stromal exhaustion—acts as a potent, stage-independent driver of early mortality, and we distill a clinically translatable **12-gene LASSO prognostic risk signature**.

---

## Key Statistical Findings
* **Primary CIBERSORT Survival & Multivariate Cox Modeling:** Low M2 macrophage fractional abundance predicts rapid early mortality (log-rank $p = 0.001$) and acts as a major independent risk factor (HR = 3.70; 95% CI: 2.16–6.30; $p < 0.001$) after adjusting for patient age, gender, and AJCC pathologic stage.
* **Orthogonal Multi-Algorithm Validation:** 
  * **quanTIseq (Absolute Physical Cell Density):** Independently confirms the mortality hazard of physical M2 depletion (log-rank $p = 0.037$; Multivariate HR = 2.25; 95% CI: 1.22–4.15; $p = 0.00976$).
  * **xCell (ssGSEA Transcriptomic Signaling):** Reveals that preserved M1 pro-inflammatory signaling independently protects overall survival (log-rank $p = 0.017$; Multivariate HR = 1.84; 95% CI: 1.14–2.98; $p = 0.0129$).
* **Independence from TNM Clinical Staging:** Two-sided Wilcoxon rank-sum testing confirms no significant difference in M2 infiltration between Early-Stage (Stages I–II) and Late-Stage (Stages III–IV) OSCC tumors ($p = 0.21$).
* **12-Gene LASSO Prognostic Signature:** L1-penalized Cox regression across 60,660 transcripts isolated an optimal 12-gene panel anchored by canonical myeloid and metabolic regulators (*CD33*, *VSIG4*, *ALOX5*, *NR1H3*, *CD52*, and *HTRA4*).

---

## Repository Files

* `Survival_Analysis.R`: Master R script executing data filtering, multi-algorithm survival analysis (Kaplan-Meier & log-rank), multivariate Cox proportional-hazards models, Wilcoxon stage correlation, Pearson transcriptomic correlation funnel, and 10-fold cross-validated LASSO Cox regression.
* `OSCC_Clinical_Data.csv`: Anatomically filtered TCGA-OSCC clinical metadata ($N = 220$ primary oral cavity tumors).
* `TIMER_Results.csv`: Multi-algorithm immune deconvolution output matrix generated via TIMER3 (CIBERSORT, quanTIseq, and xCell estimates).
* `Table_1_LASSO_12_Gene_Signature.csv`: Final 12-gene LASSO Cox regression risk coefficients and HGNC symbols.
* `LASSO Final Risk Signature.txt`: Raw R console output of the non-zero LASSO coefficients at $\lambda_{\text{min}}$.
* `Multivariate Validarion QUANTISEQ M2 and xCell M1.txt`: Full R console statistical outputs for the orthogonal multivariate Cox regression models ($n = 185$, events = 72).
* `Top 10 Correlated Genes to M2.txt`: Top-ranking transcripts co-expressed with CIBERSORT M2 macrophage abundance.
* `Figure_1A` to `Figure_6` (`.png`): High-resolution (300 DPI) survival curves, forest plots, clinical stage boxplot, and LASSO coefficient barplot.

*> **Note on Large Genomic Files (>100 MB):** Due to GitHub file size limits, the serialized full transcriptomic matrix (`OSCC_Genomic_Data.rds`, 126.4 MB; 60,660 transcripts × 220 patients) is not stored directly in this git tree. It is automatically retrieved from the NCI Genomic Data Commons (GDC) via the `TCGAbiolinks` acquisition workflow or available upon request.*

---

## System Requirements & R Packages

All analyses were conducted in **R (version 4.5.3)**. Install the required CRAN and Bioconductor dependencies prior to running the pipeline:

```R
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# CRAN packages
install.packages(c("survival", "survminer", "ggpubr", "glmnet"))

# Bioconductor packages
BiocManager::install(c("TCGAbiolinks", "SummarizedExperiment", "org.Hs.eg.db"))
