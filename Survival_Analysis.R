```R
# =========================================================================
# TITLE: Multi-Algorithm Immune Deconvolution & LASSO Pipeline in OSCC
# JOURNAL: Oral Surgery, Oral Medicine, Oral Pathology and Oral Radiology (OOOO)
# COHORT: TCGA-HNSC Anatomically Filtered Oral Cavity Cohort (N = 220; n = 219 OS)
# =========================================================================

# -------------------------------------------------------------------------
# 0. PACKAGE INSTALLATION & ENVIRONMENT SETUP
# -------------------------------------------------------------------------
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (pkg in c("survival", "survminer", "ggpubr", "glmnet")) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) install.packages(pkg)
}
for (bpkg in c("SummarizedExperiment", "org.Hs.eg.db")) {
  if (!require(bpkg, character.only = TRUE, quietly = TRUE)) BiocManager::install(bpkg, update = FALSE, ask = FALSE)
}

library(survival)
library(survminer)
library(ggpubr)
library(glmnet)
library(SummarizedExperiment)
library(org.Hs.eg.db)

# Set working directory (Adjust path if running locally in C:/TCGA)
if (dir.exists("C:/TCGA")) setwd("C:/TCGA")
dir.create("Results", showWarnings = FALSE)

# -------------------------------------------------------------------------
# PHASE 1: CLINICAL & TIMER3 IMMUNE DECONVOLUTION DATA PREPARATION
# -------------------------------------------------------------------------
clinical <- read.csv("OSCC_Clinical_Data.csv", row.names = 1)

clinical$OS_time <- ifelse(clinical$vital_status == "Dead", 
                           clinical$days_to_death, 
                           clinical$days_to_last_follow_up)
clinical$OS_time_months <- as.numeric(clinical$OS_time) / 30.41
clinical$OS_status <- ifelse(clinical$vital_status == "Dead", 1, 0)
clinical <- clinical[!is.na(clinical$OS_time_months) & clinical$OS_time_months > 0, ]

immune_data <- read.csv("TIMER_Results.csv", row.names = 1, check.names = FALSE)
immune_data <- as.data.frame(t(immune_data))

immune_data$M1_CIB   <- as.numeric(as.character(immune_data$`Macrophage M1_CIBERSORT`))
immune_data$M2_Macs  <- as.numeric(as.character(immune_data$`Macrophage M2_CIBERSORT`))
immune_data$M1_quant <- as.numeric(as.character(immune_data$`Macrophage M1_QUANTISEQ`))
immune_data$M2_quant <- as.numeric(as.character(immune_data$`Macrophage M2_QUANTISEQ`))
immune_data$M1_xcell <- as.numeric(as.character(immune_data$`Macrophage M1_XCELL`))
immune_data$M2_xcell <- as.numeric(as.character(immune_data$`Macrophage M2_XCELL`))

merged_data <- merge(clinical, immune_data, by = "row.names")
rownames(merged_data) <- merged_data$Row.names

# -------------------------------------------------------------------------
# PHASE 2: PRIMARY KAPLAN-MEIER SURVIVAL ANALYSIS (CIBERSORT M1 & M2)
# -------------------------------------------------------------------------
# Figure 1A: CIBERSORT M1
cut_cib_m1 <- surv_cutpoint(merged_data, time = "OS_time_months", event = "OS_status", variables = "M1_CIB")
merged_data$M1_CIB_Group <- ifelse(merged_data$M1_CIB >= cut_cib_m1$cutpoint$cutpoint, "High", "Low")
fit_cib_M1 <- survfit(Surv(OS_time_months, OS_status) ~ M1_CIB_Group, data = merged_data)

plot_cib_M1 <- ggsurvplot(fit_cib_M1, data = merged_data, pval = TRUE, conf.int = TRUE, risk.table = TRUE,
                          palette = c("red", "blue"), title = "OSCC Survival: M1 Macrophage Infiltration (CIBERSORT)",
                          xlab = "Time in Months", ylab = "Overall Survival Probability")

# Figure 1B: CIBERSORT M2
cut_cib_m2 <- surv_cutpoint(merged_data, time = "OS_time_months", event = "OS_status", variables = "M2_Macs")
merged_data$M2_Group <- ifelse(merged_data$M2_Macs >= cut_cib_m2$cutpoint$cutpoint, "High", "Low")
fit_cib_M2 <- survfit(Surv(OS_time_months, OS_status) ~ M2_Group, data = merged_data)

plot_cib_M2 <- ggsurvplot(fit_cib_M2, data = merged_data, pval = TRUE, conf.int = TRUE, risk.table = TRUE,
                          palette = c("red", "blue"), title = "OSCC Survival: M2 Macrophage Infiltration (CIBERSORT)",
                          xlab = "Time in Months", ylab = "Overall Survival Probability")

# -------------------------------------------------------------------------
# PHASE 3: PRIMARY MULTIVARIATE COX REGRESSION (CIBERSORT M2 - FIGURE 2)
# -------------------------------------------------------------------------
cox_model <- coxph(Surv(OS_time_months, OS_status) ~ M2_Group + age_at_index + gender + ajcc_pathologic_stage, data = merged_data)
summary(cox_model)
ggforest(cox_model, data = merged_data, main = "Multivariate Hazard Ratios: OSCC Cohort", fontsize = 1.0)

# -------------------------------------------------------------------------
# PHASE 4: ORTHOGONAL MULTI-ALGORITHM VALIDATION (quanTIseq & xCell)
# -------------------------------------------------------------------------
# Figure 3A: quanTIseq M1
cut_quant_m1 <- surv_cutpoint(merged_data, time = "OS_time_months", event = "OS_status", variables = "M1_quant")
merged_data$M1_quant_Group <- ifelse(merged_data$M1_quant >= cut_quant_m1$cutpoint$cutpoint, "High", "Low")
fit_quant_M1 <- survfit(Surv(OS_time_months, OS_status) ~ M1_quant_Group, data = merged_data)
plot_quant_M1 <- ggsurvplot(fit_quant_M1, data = merged_data, pval = TRUE, conf.int = TRUE, risk.table = TRUE,
                            palette = c("red", "blue"), title = "OSCC Survival: M1 Macrophage Infiltration (quanTIseq)",
                            xlab = "Time in Months", ylab = "Overall Survival Probability")

# Figure 3B: quanTIseq M2
cut_quant_m2 <- surv_cutpoint(merged_data, time = "OS_time_months", event = "OS_status", variables = "M2_quant")
merged_data$M2_quant_Group <- ifelse(merged_data$M2_quant >= cut_quant_m2$cutpoint$cutpoint, "High", "Low")
fit_quant_M2 <- survfit(Surv(OS_time_months, OS_status) ~ M2_quant_Group, data = merged_data)
plot_quant_M2 <- ggsurvplot(fit_quant_M2, data = merged_data, pval = TRUE, conf.int = TRUE, risk.table = TRUE,
                            palette = c("red", "blue"), title = "OSCC Survival: M2 Macrophage Infiltration (quanTIseq)",
                            xlab = "Time in Months", ylab = "Overall Survival Probability")

# Figure 3C: xCell M1
cut_xcell_m1 <- surv_cutpoint(merged_data, time = "OS_time_months", event = "OS_status", variables = "M1_xcell")
merged_data$M1_xcell_Group <- ifelse(merged_data$M1_xcell >= cut_xcell_m1$cutpoint$cutpoint, "High", "Low")
fit_xcell_M1 <- survfit(Surv(OS_time_months, OS_status) ~ M1_xcell_Group, data = merged_data)
plot_xcell_M1 <- ggsurvplot(fit_xcell_M1, data = merged_data, pval = TRUE, conf.int = TRUE, risk.table = TRUE,
                            palette = c("red", "blue"), title = "OSCC Survival: M1 Macrophage Infiltration (xCell)",
                            xlab = "Time in Months", ylab = "Overall Survival Probability")

# Figure 3D: xCell M2
cut_xcell_m2 <- surv_cutpoint(merged_data, time = "OS_time_months", event = "OS_status", variables = "M2_xcell")
merged_data$M2_xcell_Group <- ifelse(merged_data$M2_xcell >= cut_xcell_m2$cutpoint$cutpoint, "High", "Low")
fit_xcell_M2 <- survfit(Surv(OS_time_months, OS_status) ~ M2_xcell_Group, data = merged_data)
plot_xcell_M2 <- ggsurvplot(fit_xcell_M2, data = merged_data, pval = TRUE, conf.int = TRUE, risk.table = TRUE,
                            palette = c("red", "blue"), title = "OSCC Survival: M2 Macrophage Infiltration (xCell)",
                            xlab = "Time in Months", ylab = "Overall Survival Probability")

# Figure 4A & 4B: Orthogonal Multivariate Cox Models (quanTIseq M2 & xCell M1)
cox_quant <- coxph(Surv(OS_time_months, OS_status) ~ M2_quant_Group + age_at_index + gender + ajcc_pathologic_stage, data = merged_data)
summary(cox_quant)
ggforest(cox_quant, data = merged_data, main = "Multivariate Hazard Ratios: quanTIseq M2 Density", fontsize = 1.0)

cox_xcell <- coxph(Surv(OS_time_months, OS_status) ~ M1_xcell_Group + age_at_index + gender + ajcc_pathologic_stage, data = merged_data)
summary(cox_xcell)
ggforest(cox_xcell, data = merged_data, main = "Multivariate Hazard Ratios: xCell M1 Signaling", fontsize = 1.0)

# -------------------------------------------------------------------------
# PHASE 5: CLINICAL STAGE CORRELATION (WILCOXON RANK-SUM TEST - FIGURE 5)
# -------------------------------------------------------------------------
merged_data$Stage_Group <- ifelse(merged_data$ajcc_pathologic_stage %in% c("Stage I", "Stage II"), "Early (I-II)",
                                  ifelse(merged_data$ajcc_pathologic_stage %in% c("Stage III", "Stage IVA", "Stage IVB", "Stage IVC"), "Late (III-IV)", NA))
plot_data <- merged_data[!is.na(merged_data$Stage_Group), ]

stage_plot <- ggboxplot(plot_data, x = "Stage_Group", y = "M2_Macs", color = "Stage_Group", 
                        palette = c("#00AFBB", "#FC4E07"), add = "jitter", add.params = list(alpha = 0.5),
                        title = "M2 Macrophage Infiltration by OSCC Clinical Stage",
                        xlab = "Pathological Stage", ylab = "M2 Macrophage Fraction (CIBERSORT)") +
              stat_compare_means(method = "wilcox.test", label.y = max(plot_data$M2_Macs, na.rm = TRUE) * 1.05)
print(stage_plot)

# -------------------------------------------------------------------------
# PHASE 6 & 7: GENOMIC CORRELATION FUNNEL & LASSO COX REGRESSION
# -------------------------------------------------------------------------
if (file.exists("OSCC_Genomic_Data.rds")) {
  genomic_data <- readRDS("OSCC_Genomic_Data.rds")
  genomic_matrix <- assay(genomic_data)
  genomic_t <- as.data.frame(t(genomic_matrix))
  
  genomic_t$Patient_ID <- substr(rownames(genomic_t), 1, 12)
  merged_data$Patient_ID <- substr(merged_data$bcr_patient_barcode, 1, 12)
  
  ml_data <- merge(merged_data[, c("Patient_ID", "OS_time_months", "OS_status", "M2_Macs")], 
                   genomic_t, by = "Patient_ID")
  
  gene_columns <- ml_data[, 5:ncol(ml_data)]
  correlations <- cor(ml_data$M2_Macs, gene_columns, use = "pairwise.complete.obs")
  
  cor_df <- data.frame(Gene = colnames(correlations), Correlation = as.vector(correlations))
  cor_df <- cor_df[!is.na(cor_df$Correlation), ]
  cor_df <- cor_df[order(-cor_df$Correlation), ]
  
  top_100_genes <- cor_df$Gene[1:100]
  ml_clean <- ml_data[!is.na(ml_data$OS_time_months) & !is.na(ml_data$OS_status), ]
  
  x_matrix <- as.matrix(ml_clean[, top_100_genes])
  y_survival <- Surv(ml_clean$OS_time_months, ml_clean$OS_status)
  
  set.seed(123)
  lasso_cv <- cv.glmnet(x = x_matrix, y = y_survival, family = "cox", alpha = 1)
  best_lambda <- lasso_cv$lambda.min
  lasso_coefs <- as.matrix(coef(lasso_cv, s = best_lambda))
  final_signature <- lasso_coefs[lasso_coefs[, 1] != 0, , drop = FALSE]
  print(final_signature)
}

# -------------------------------------------------------------------------
# PHASE 8: FIGURE 6 - 12-GENE LASSO PROGNOSTIC SIGNATURE BARPLOT
# -------------------------------------------------------------------------
signature_df <- data.frame(
  Gene = c("ENSG00000225331", "CD33", "VSIG4", "ALOX5", 
           "GALM", "CD52", "NR1H3", "ENSG00000257838", 
           "C19orf38", "HTRA4", "MIR3648-1", "ENSG00000229671"),
  Coefficient = c(0.01691193, 0.002886999, 0.000249641, 0.00002169772,
                  -0.000004034941, -0.00001273053, -0.000297142, -0.0005140829,
                  -0.001921628, -0.002426121, -0.02332472, -0.03579118)
)
signature_df$Effect <- ifelse(signature_df$Coefficient > 0, "Increases Risk (Hazard)", "Protective")

lasso_bar_plot <- ggbarplot(signature_df, x = "Gene", y = "Coefficient",
                            fill = "Effect", color = "Effect",
                            palette = c("#FC4E07", "#00AFBB"),
                            sort.val = "asc", sort.by.groups = FALSE,
                            orientation = "horiz",
                            title = "LASSO-Derived 12-Gene Prognostic Signature in OSCC",
                            xlab = "Transcript / Gene Symbol",
                            ylab = "LASSO Cox Regression Coefficient (Weight)") +
                  geom_hline(yintercept = 0, linetype = "dashed", color = "black")
print(lasso_bar_plot)
