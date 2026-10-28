# =========================================================================
# SCRIPT 01: TCGA-HNSC Data Acquisition, Anatomical Filtering & TIMER3 Prep
# TITLE: Multi-Algorithm Immune Deconvolution & LASSO Pipeline in OSCC
# DESCRIPTION: Queries GDC via TCGAbiolinks, isolates primary oral cavity 
#              squamous cell carcinoma (OSCC, N = 220), saves clinical/genomic 
#              objects, and formats the TPM matrix for TIMER3 deconvolution.
# =========================================================================

# -------------------------------------------------------------------------
# 0. PACKAGE INSTALLATION & SETUP
# -------------------------------------------------------------------------
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (bpkg in c("TCGAbiolinks", "SummarizedExperiment")) {
  if (!require(bpkg, character.only = TRUE, quietly = TRUE)) {
    BiocManager::install(bpkg, update = FALSE, ask = FALSE)
  }
}

library(TCGAbiolinks)
library(SummarizedExperiment)

dir.create("C:/TCGA", showWarnings = FALSE)
setwd("C:/TCGA")


# -------------------------------------------------------------------------
# PHASE 1: QUERY & DOWNLOAD TCGA-HNSC TRANSCRIPTOMIC DATA (STAR - COUNTS)
# -------------------------------------------------------------------------
cat("Querying NCI Genomic Data Commons (GDC) for TCGA-HNSC...\n")

query_hnsc <- GDCquery(
  project = "TCGA-HNSC",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = "Primary Tumor"
)

# Download and prepare SummarizedExperiment object (Initial N = 566)
GDCdownload(query_hnsc, directory = "C:/TCGA/TCGA-HNSC")
hnsc_se <- GDCprepare(query_hnsc, directory = "C:/TCGA/TCGA-HNSC")


# -------------------------------------------------------------------------
# PHASE 2: RIGOROUS ANATOMICAL FILTERING FOR PRIMARY ORAL CAVITY (OSCC)
# -------------------------------------------------------------------------
# Define strict anatomical oral cavity topographical sites
oral_cavity_sites <- c(
  "Cheek mucosa",
  "Tongue, NOS",
  "Mouth, NOS",
  "Ventral surface of tongue, NOS",
  "Floor of mouth, NOS",
  "Gum, NOS",
  "Hard palate, NOS"
)

clinical_all <- as.data.frame(colData(hnsc_se))

# Filter strictly for Oral Cavity primary sites and unique patients (N = 220)
oscc_mask <- clinical_all$tissue_or_organ_of_origin %in% oral_cavity_sites
oscc_se <- hnsc_se[, oscc_mask]
oscc_se <- oscc_se[, !duplicated(oscc_se$bcr_patient_barcode)]

cat("Total primary OSCC patients isolated:", ncol(oscc_se), "\n")


# -------------------------------------------------------------------------
# PHASE 3: EXPORT PURIFIED CLINICAL CSV & SERIALIZED GENOMIC OBJECT (.RDS)
# -------------------------------------------------------------------------
# 1. Flatten list columns in clinical metadata and export to CSV
oscc_clinical <- as.data.frame(colData(oscc_se))
list_cols <- sapply(oscc_clinical, is.list)
oscc_clinical[list_cols] <- lapply(oscc_clinical[list_cols], function(x) sapply(x, paste, collapse = "; "))
rownames(oscc_clinical) <- oscc_clinical$bcr_patient_barcode

write.csv(oscc_clinical, "C:/TCGA/OSCC_Clinical_Data.csv", row.names = TRUE)

# 2. Save purified OSCC SummarizedExperiment object (60,660 transcripts x 220 patients)
saveRDS(oscc_se, "C:/TCGA/OSCC_Genomic_Data.rds")


# -------------------------------------------------------------------------
# PHASE 4: TPM NORMALIZATION, GENE MAPPING & ROUNDING FOR TIMER3 DECONVOLUTION
# -------------------------------------------------------------------------
# 1. Extract Transcripts Per Million (TPM) assay and HGNC gene symbols
tpm_matrix <- assay(oscc_se, "tpm_unstrand")
colnames(tpm_matrix) <- oscc_se$bcr_patient_barcode
gene_symbols <- rowData(oscc_se)$gene_name

# 2. Filter valid gene symbols and aggregate duplicate gene identifiers by mean
valid_genes <- !is.na(gene_symbols) & gene_symbols != ""
tpm_valid <- tpm_matrix[valid_genes, ]
symbols_valid <- gene_symbols[valid_genes]

cat("Aggregating duplicate gene symbols by mean expression...\n")
tpm_aggregated <- rowsum(tpm_valid, group = symbols_valid) / as.vector(table(symbols_valid))

# 3. Round TPM values to 2 decimal places to meet web-server size constraints
tpm_rounded <- round(tpm_aggregated, 2)

# 4. Export formatted mixture matrix for TIMER3 web server upload
write.csv(tpm_rounded, "C:/TCGA/TIMER_Mixture_Rounded.csv", row.names = TRUE, quote = FALSE)

cat("\nPreprocessing complete! Upload 'TIMER_Mixture_Rounded.csv' to TIMER3 (HNSC) and save output as 'TIMER_Results.csv'.\n")