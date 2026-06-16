suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(survival)
  library(patchwork)
})

##############################################
####### Figure5 setup
##############################################

rm(list = ls())

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(survival)
  library(survminer)
  library(glmnet)
  library(ggplot2)
})

############################################################
## 1. Paths
############################################################

workdir <- "E:/code/400T/TCGA/TCGA-BRCA"

expr_file <- file.path(workdir, "TCGA-BRCA.star_fpkm.tsv.gz")
surv_file <- file.path(workdir, "TCGA-BRCA.survival.tsv.gz")
mut_file  <- file.path(workdir, "TCGA-BRCA.somaticmutation_wxs.tsv")

outdir <- file.path(workdir, "TCGA_BRCA_tumor_only_Cox_LASSO_mutation")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

############################################################
## 2. Candidate genes
############################################################

genes <- c(
  "ARID1B", "MUC17", "MUC5AC", "LRP1", "C3", "LAMA1",
  "MUC16", "LAMA5", "ADAMTSL3", "KMT2D", "SMARCA4",
  "HSPG2", "SETD2", "CCL2", "JAK1", "CXCR4", "FN1",
  "MUC5B", "TGFB1", "CREBBP", "COL5A1", "SMARCB1",
  "CFB", "MMP14", "COL1A1", "ITGAV", "CTLA4", "KRT18",
  "FCGBP", "EP300", "NFKBIA", "DNMT3A", "ARID1A", "FAP",
  "CCR2", "NFKB1", "LAG3", "STAT3", "MMP9", "IKBKB",
  "COL6A1", "TGFBR1", "MUC4", "PDCD1", "MMP2", "NLRP3",
  "TET2"
)

############################################################
## 3. Read the expression matrix and keep only tumor samples labeled 01
############################################################

expr_raw <- fread(expr_file, data.table = FALSE)

cat("Expression matrix dimension:\n")
print(dim(expr_raw))
print(colnames(expr_raw)[1:10])

colnames(expr_raw)[1] <- "Gene"

expr_raw <- expr_raw %>%
  mutate(
    Ensembl_ID = as.character(Gene),
    Ensembl_ID = str_replace(Ensembl_ID, "\\..*$", "")
  )

# Identify TCGA sample columns in the expression matrix
all_sample_cols <- colnames(expr_raw)[str_detect(colnames(expr_raw), "^TCGA-")]

# Keep only tumor samples: positions 14-15 of the TCGA barcode are 01
tumor_sample_cols <- all_sample_cols[
  substr(all_sample_cols, 14, 15) == "01"
]

cat("All expression samples:", length(all_sample_cols), "\n")
cat("Tumor expression samples:", length(tumor_sample_cols), "\n")

if (length(tumor_sample_cols) == 0) {
  stop("表达矩阵中没有识别到 TCGA 01 肿瘤样本，请检查样本名格式。")
}

############################################################
## 4. Convert Ensembl IDs to gene symbols
############################################################

if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install("org.Hs.eg.db")
}

suppressPackageStartupMessages({
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

expr_raw$Symbol <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = expr_raw$Ensembl_ID,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

print(head(expr_raw[, c("Gene", "Ensembl_ID", "Symbol")]))

############################################################
## 5. Extract candidate gene expression
############################################################

expr_raw2 <- expr_raw %>%
  filter(Symbol %in% genes)

matched_genes <- sort(unique(expr_raw2$Symbol))
missing_genes <- setdiff(genes, matched_genes)

cat("Matched genes:\n")
print(matched_genes)

cat("Missing genes:\n")
print(missing_genes)

if (nrow(expr_raw2) == 0) {
  stop("没有匹配到候选基因。")
}

# Keep only Symbol plus tumor sample columns
expr_raw2 <- expr_raw2 %>%
  dplyr::select(Symbol, all_of(tumor_sample_cols))

expr_raw2[, tumor_sample_cols] <- lapply(
  expr_raw2[, tumor_sample_cols],
  as.numeric
)

# Average values when multiple Ensembl IDs map to the same symbol
expr_mat_df <- expr_raw2 %>%
  group_by(Symbol) %>%
  summarise(
    across(
      all_of(tumor_sample_cols),
      ~ mean(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )

expr_mat <- expr_mat_df %>%
  column_to_rownames("Symbol") %>%
  as.data.frame()

cat("Candidate gene expression matrix:\n")
print(dim(expr_mat))
print(rownames(expr_mat)[1:10])

############################################################
## 6. Transpose the expression matrix: samples x genes
############################################################

expr_t <- as.data.frame(t(expr_mat))

expr_t <- expr_t %>%
  rownames_to_column("sample_id") %>%
  mutate(
    sample_type = substr(sample_id, 14, 15),
    patient_id = substr(sample_id, 1, 12)
  ) %>%
  filter(sample_type == "01")

genes_expr <- intersect(genes, colnames(expr_t))

cat("Genes after transpose:", length(genes_expr), "\n")
print(genes_expr)

if (length(genes_expr) < 2) {
  stop("转置后可用基因少于 2 个。")
}

# Average multiple tumor samples from the same patient
expr_t_avg <- expr_t %>%
  dplyr::select(patient_id, all_of(genes_expr)) %>%
  group_by(patient_id) %>%
  summarise(
    across(
      all_of(genes_expr),
      ~ mean(as.numeric(.x), na.rm = TRUE)
    ),
    .groups = "drop"
  )

# FPKM log2 transformation
expr_t_avg <- expr_t_avg %>%
  mutate(
    across(
      all_of(genes_expr),
      ~ log2(.x + 1)
    )
  )

cat("Tumor expression data after transpose:\n")
print(dim(expr_t_avg))
print(head(expr_t_avg[, 1:min(8, ncol(expr_t_avg))]))

############################################################
## 7. Read survival data and keep only tumor samples labeled 01
############################################################

surv_raw <- fread(surv_file, data.table = FALSE)

cat("Survival columns:\n")
print(colnames(surv_raw))
print(head(surv_raw))

if (!all(c("OS", "OS.time") %in% colnames(surv_raw))) {
  stop("survival 文件中没有 OS 或 OS.time。")
}

sample_col <- if ("sample" %in% colnames(surv_raw)) {
  "sample"
} else if ("sample_id" %in% colnames(surv_raw)) {
  "sample_id"
} else {
  colnames(surv_raw)[1]
}

surv_df <- surv_raw %>%
  mutate(
    sample_id = as.character(.data[[sample_col]]),
    sample_type = substr(sample_id, 14, 15),
    patient_id = substr(sample_id, 1, 12),
    OS_time = as.numeric(OS.time),
    OS_event = as.numeric(OS)
  ) %>%
  filter(
    sample_type == "01",
    !is.na(OS_time),
    !is.na(OS_event),
    OS_time > 0
  ) %>%
  dplyr::select(patient_id, OS_time, OS_event) %>%
  distinct(patient_id, .keep_all = TRUE)

cat("Tumor survival data:\n")
print(dim(surv_df))
print(head(surv_df))

############################################################
## 8. Merge expression and survival data
############################################################

cox_df <- inner_join(
  surv_df,
  expr_t_avg,
  by = "patient_id"
)

genes_use <- intersect(genes, colnames(cox_df))

cat("Final sample number:", nrow(cox_df), "\n")
cat("Final gene number:", length(genes_use), "\n")
print(genes_use)

if (length(genes_use) < 2) {
  stop("可用于 Cox/LASSO 的候选基因少于 2 个。")
}

write.csv(
  cox_df,
  file.path(outdir, "00_TCGA_BRCA_tumor_only_expression_survival_candidate_genes.csv"),
  row.names = FALSE
)

############################################################
## 9. Univariate Cox analysis
############################################################

univ_cox_res <- lapply(genes_use, function(gene) {
  
  tmp <- cox_df %>%
    dplyr::select(OS_time, OS_event, all_of(gene)) %>%
    mutate(
      OS_time = as.numeric(OS_time),
      OS_event = as.numeric(OS_event),
      gene_expr = as.numeric(.data[[gene]])
    ) %>%
    filter(
      !is.na(OS_time),
      !is.na(OS_event),
      !is.na(gene_expr),
      OS_time > 0
    )
  
  if (nrow(tmp) < 10 || length(unique(tmp$gene_expr)) < 2) {
    return(data.frame(
      gene = gene,
      HR = NA,
      coef = NA,
      lower95 = NA,
      upper95 = NA,
      pvalue = NA,
      n = nrow(tmp),
      event = sum(tmp$OS_event == 1, na.rm = TRUE),
      note = "insufficient variation"
    ))
  }
  
  fit <- tryCatch(
    coxph(Surv(OS_time, OS_event) ~ gene_expr, data = tmp),
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    return(data.frame(
      gene = gene,
      HR = NA,
      coef = NA,
      lower95 = NA,
      upper95 = NA,
      pvalue = NA,
      n = nrow(tmp),
      event = sum(tmp$OS_event == 1, na.rm = TRUE),
      note = "cox failed"
    ))
  }
  
  fit_sum <- summary(fit)
  
  data.frame(
    gene = gene,
    HR = fit_sum$coefficients[1, "exp(coef)"],
    coef = fit_sum$coefficients[1, "coef"],
    lower95 = fit_sum$conf.int[1, "lower .95"],
    upper95 = fit_sum$conf.int[1, "upper .95"],
    pvalue = fit_sum$coefficients[1, "Pr(>|z|)"],
    n = nrow(tmp),
    event = sum(tmp$OS_event == 1, na.rm = TRUE),
    note = "ok",
    stringsAsFactors = FALSE
  )
}) %>%
  bind_rows() %>%
  mutate(
    FDR = p.adjust(pvalue, method = "BH"),
    sig = case_when(
      pvalue < 0.1 & HR > 1 ~ "Risk",
      pvalue < 0.1 & HR < 1 ~ "Protective",
      TRUE ~ "NS"
    )
  ) %>%
  arrange(pvalue)

write.csv(
  univ_cox_res,
  file.path(outdir, "01_univariate_cox_candidate_genes.csv"),
  row.names = FALSE
)

print(univ_cox_res)




##############################################
#######Figure5A
##############################################

forest_df_noMUC17 <- univ_cox_res %>%
  dplyr::filter(
    !is.na(HR),
    !is.na(lower95),
    !is.na(upper95),
    gene != "MUC17"
  ) %>%
  dplyr::mutate(
    sig = dplyr::case_when(
      pvalue < 0.1 & HR > 1 ~ "Risk",
      pvalue < 0.1 & HR < 1 ~ "Protective",
      TRUE ~ "NS"
    )
  ) %>%
  dplyr::arrange(HR) %>%
  dplyr::mutate(
    gene = factor(gene, levels = gene)
  )


p_forest_noMUC17 <- ggplot(
  forest_df_noMUC17,
  aes(
    x = HR,
    y = gene,
    color = sig
  )
) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    color = "grey40",
    linewidth = 0.5
  ) +
  geom_errorbarh(
    aes(
      xmin = lower95,
      xmax = upper95
    ),
    height = 0.25,
    linewidth = 0.55
  ) +
  geom_point(size = 2.5) +
  scale_x_log10() +
  scale_color_manual(
    values = c(
      "Risk" = "#E64B35",
      "Protective" = "#4DBBD5",
      "NS" = "#7A7A7A"
    )
  ) +
  labs(
    x = "Hazard ratio per log2(FPKM + 1) unit",
    y = NULL,
    color = NULL,
    title = "Univariate Cox regression",
    subtitle = "TCGA-BRCA tumor samples only; MUC17 excluded"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 14
    ),
    plot.subtitle = element_text(
      hjust = 0.5,
      size = 9,
      color = "grey30"
    ),
    axis.text.y = element_text(
      size = 8,
      color = "black"
    ),
    axis.text.x = element_text(
      color = "black"
    ),
    legend.position = "right"
  )

p_forest_noMUC17







##############################################
#######Figure5C
##############################################


library(survival)
library(survminer)
library(ggplot2)

## 1. Ensure risk_group order: Low risk as the reference group
cox_lasso_df$risk_group <- factor(
  cox_lasso_df$risk_group,
  levels = c("Low risk", "High risk")
)

## Use this if the original groups are Low / High
# cox_lasso_df$risk_group <- factor(
#   cox_lasso_df$risk_group,
#   levels = c("Low", "High"),
#   labels = c("Low risk", "High risk")
# )

## 2. Fit KM curves
fit_km <- survfit(
  Surv(OS_time, OS_event) ~ risk_group,
  data = cox_lasso_df
)

## 3. Cox regression: High risk vs Low risk
cox_fit <- coxph(
  Surv(OS_time, OS_event) ~ risk_group,
  data = cox_lasso_df
)

cox_sum <- summary(cox_fit)

HR <- cox_sum$conf.int[1, "exp(coef)"]
HR_low <- cox_sum$conf.int[1, "lower .95"]
HR_high <- cox_sum$conf.int[1, "upper .95"]
cox_p <- cox_sum$coefficients[1, "Pr(>|z|)"]

## 4. Generate HR text
hr_label <- paste0(
  "High vs Low: HR=", sprintf("%.2f", HR),
  " (95% CI ",
  sprintf("%.2f", HR_low), "-",
  sprintf("%.2f", HR_high), ")",
  ", P=",
  ifelse(
    cox_p < 0.001,
    "<0.001",
    sprintf("%.3f", cox_p)
  )
)

## 5. Draw KM plot
p_km <- ggsurvplot(
  fit_km,
  data = cox_lasso_df,
  pval = TRUE,
  conf.int = FALSE,
  risk.table = TRUE,
  risk.table.height = 0.25,
  palette = c("#4DBBD5", "#E64B35"),
  xlab = "Time (days)",
  ylab = "Overall survival probability",
  legend.title = "",
  legend.labs = c("Low risk", "High risk"),
  legend = "right",
  title = "LASSO risk score"
)

## 6. Main survival plot: square panel with HR in the lower-left corner
p_km$plot <- p_km$plot +
  theme_bw() +
  theme(
    aspect.ratio = 1,
    plot.title = element_text(
      hjust = 0.5,
      size = 11
    ),
    axis.text = element_text(
      size = 9,
      color = "black"
    ),
    axis.title = element_text(
      size = 10
    ),
    legend.title = element_text(
      size = 9
    ),
    legend.text = element_text(
      size = 8
    )
  ) +
  annotate(
    "text",
    x = max(cox_lasso_df$OS_time, na.rm = TRUE) * 0.04,
    y = 0.08,
    label = hr_label,
    hjust = 0,
    vjust = 0,
    size = 2.2
  )

## 7. Format the risk table
p_km$table <- p_km$table +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      size = 8,
      color = "black"
    ),
    axis.text.y = element_text(
      size = 8,
      color = "black"
    ),
    axis.title.x = element_text(
      size = 9
    ),
    axis.title.y = element_text(
      size = 9
    )
  )

## 8. Output the plot
p_km

## 9. Save PDF
pdf(
  file = "E:/code/400T/LASOO_COX/TCGA_BRCA_tumor_only_Cox_LASSO_mutation/LASSO_risk_score_KM_square.pdf",
  width = 5.5,
  height = 6.5
)

print(p_km)

dev.off()










library(survival)
library(survminer)
library(ggplot2)
library(coxphf)

## If coxphf is not installed, run this once first
## install.packages("coxphf")

## 1. Ensure group order: Wildtype as the reference group
surv_mut_integrated_df$mutation_integrated_group <- factor(
  surv_mut_integrated_df$mutation_integrated_group,
  levels = c(
    "Wildtype",
    "HR<1 mutant only",
    "HR>1 mutant only",
    "Both mutant"
  )
)

## 2. Check sample and event counts in each group
aggregate(
  OS_event ~ mutation_integrated_group,
  data = surv_mut_integrated_df,
  FUN = function(x) c(
    n = length(x),
    event = sum(x == 1, na.rm = TRUE),
    censored = sum(x == 0, na.rm = TRUE)
  )
)

## 3. Fit KM curves
fit_integrated_km <- survfit(
  Surv(OS_time, OS_event) ~ mutation_integrated_group,
  data = surv_mut_integrated_df
)

## 4. Firth Cox regression: calculate HR for each group relative to Wildtype
cox_integrated_firth <- coxphf(
  Surv(OS_time, OS_event) ~ mutation_integrated_group,
  data = surv_mut_integrated_df
)

## 5. Extract Firth Cox HR results
hr_df <- data.frame(
  group = names(cox_integrated_firth$coefficients),
  HR = exp(cox_integrated_firth$coefficients),
  lower95 = exp(cox_integrated_firth$ci.lower),
  upper95 = exp(cox_integrated_firth$ci.upper),
  pvalue = cox_integrated_firth$prob
)

## Remove variable name prefixes and keep only group names
hr_df$group <- gsub(
  "mutation_integrated_group",
  "",
  hr_df$group
)

## Simplify display names
hr_df$group <- gsub(
  " mutant only",
  "",
  hr_df$group
)

hr_df$group <- gsub(
  " mutant",
  "",
  hr_df$group
)

## 6. Generate HR text for the plot
hr_text <- paste0(
  "Ref: Wildtype\n",
  paste0(
    hr_df$group,
    ": HR=", sprintf("%.2f", hr_df$HR),
    ", P=", ifelse(
      hr_df$pvalue < 0.001,
      "<0.001",
      sprintf("%.3f", hr_df$pvalue)
    ),
    collapse = "\n"
  )
)

## To show 95% CI, replace the hr_text above with this
# hr_text <- paste0(
#   "Ref: Wildtype\n",
#   paste0(
#     hr_df$group,
#     ": HR=", sprintf("%.2f", hr_df$HR),
#     " (95% CI ",
#     sprintf("%.2f", hr_df$lower95),
#     "-",
#     sprintf("%.2f", hr_df$upper95),
#     "), P=",
#     ifelse(
#       hr_df$pvalue < 0.001,
#       "<0.001",
#       sprintf("%.3f", hr_df$pvalue)
#     ),
#     collapse = "\n"
#   )
# )

## 7. Draw KM plot
p_integrated_km <- ggsurvplot(
  fit_integrated_km,
  data = surv_mut_integrated_df,
  pval = TRUE,
  conf.int = FALSE,
  risk.table = TRUE,
  risk.table.height = 0.28,
  palette = c(
    "#7A7A7A",
    "#4DBBD5",
    "#E64B35",
    "#8E44AD"
  ),
  xlab = "Time (days)",
  ylab = "Overall survival probability",
  legend.title = "Mutation group",
  legend.labs = c(
    "Wildtype",
    "HR<1 mutant only",
    "HR>1 mutant only",
    "Both mutant"
  ),
  title = "Survival analysis by integrated mutation groups"
)

## 8. Main survival plot: square panel with Firth HR in the lower-left corner
p_integrated_km$plot <- p_integrated_km$plot +
  theme_bw() +
  theme(
    aspect.ratio = 1,
    plot.title = element_text(
      hjust = 0.5,
      size = 11
    ),
    axis.text = element_text(
      size = 9,
      color = "black"
    ),
    axis.title = element_text(
      size = 10
    ),
    legend.title = element_text(
      size = 9
    ),
    legend.text = element_text(
      size = 8
    )
  ) +
  annotate(
    "text",
    x = max(surv_mut_integrated_df$OS_time, na.rm = TRUE) * 0.04,
    y = 0.08,
    label = hr_text,
    hjust = 0,
    vjust = 0,
    size = 2.0
  )

## 9. Format the risk table
p_integrated_km$table <- p_integrated_km$table +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      size = 8,
      color = "black"
    ),
    axis.text.y = element_text(
      size = 8,
      color = "black"
    ),
    axis.title.x = element_text(
      size = 9
    ),
    axis.title.y = element_text(
      size = 9
    )
  )

## 10. Output the plot
p_integrated_km

## 11. Save PDF
pdf(
  file = "E:/code/400T/LASOO_COX/TCGA_BRCA_tumor_only_Cox_LASSO_mutation/integrated_KM_square_firth.pdf",
  width = 5.5,
  height = 6.5
)

print(p_integrated_km)

dev.off()




##############################################
#######Figure5D-E
##############################################


library(survival)
library(survminer)
library(ggplot2)
library(coxphf)

## 1. Ensure group order
surv_mut_only_df$mutation_integrated_group <- factor(
  surv_mut_only_df$mutation_integrated_group,
  levels = c(
    "HR<1 mutant only",
    "HR>1 mutant only",
    "Both mutant"
  )
)

## 2. Check event counts in each group
aggregate(
  OS_event ~ mutation_integrated_group,
  data = surv_mut_only_df,
  FUN = function(x) c(
    n = length(x),
    event = sum(x == 1, na.rm = TRUE),
    censored = sum(x == 0, na.rm = TRUE)
  )
)

## 3. Fit KM curves
fit_mut_only_km <- survfit(
  Surv(OS_time, OS_event) ~ mutation_integrated_group,
  data = surv_mut_only_df
)

get_pairwise_hr_firth <- function(data, group_col, group1, group2) {
  
  tmp_df <- data[data[[group_col]] %in% c(group1, group2), ]
  
  tmp_df[[group_col]] <- factor(
    tmp_df[[group_col]],
    levels = c(group1, group2)
  )
  
  fit <- coxphf(
    as.formula(
      paste0("Surv(OS_time, OS_event) ~ ", group_col)
    ),
    data = tmp_df
  )
  
  comparison_name <- paste0(
    gsub(" mutant only| mutant", "", group2),
    " vs ",
    gsub(" mutant only| mutant", "", group1)
  )
  
  data.frame(
    comparison = comparison_name,
    HR = exp(fit$coefficients),
    lower95 = fit$ci.lower,
    upper95 = fit$ci.upper,
    pvalue = fit$prob
  )
}

hr_pairwise_df <- rbind(
  get_pairwise_hr_firth(
    surv_mut_only_df,
    "mutation_integrated_group",
    "HR<1 mutant only",
    "HR>1 mutant only"
  ),
  get_pairwise_hr_firth(
    surv_mut_only_df,
    "mutation_integrated_group",
    "Both mutant",
    "HR<1 mutant only"
  ),
  get_pairwise_hr_firth(
    surv_mut_only_df,
    "mutation_integrated_group",
    "Both mutant",
    "HR>1 mutant only"
  )
)

hr_pairwise_df

## View the HR result table
hr_pairwise_df

## 6. Generate HR text for the plot
hr_text <- paste0(
  hr_pairwise_df$comparison,
  ": HR=", sprintf("%.2f", hr_pairwise_df$HR),
  " (95% CI ",
  sprintf("%.2f", hr_pairwise_df$lower95),
  "-",
  sprintf("%.2f", hr_pairwise_df$upper95),
  "), P=",
  ifelse(
    hr_pairwise_df$pvalue < 0.001,
    "<0.001",
    sprintf("%.3f", hr_pairwise_df$pvalue)
  ),
  collapse = "\n"
)

## 7. Draw KM plot
p_mut_only_km <- ggsurvplot(
  fit_mut_only_km,
  data = surv_mut_only_df,
  pval = TRUE,
  conf.int = FALSE,
  risk.table = TRUE,
  risk.table.height = 0.25,
  palette = c("#4DBBD5", "#E64B35", "#8E44AD"),
  xlab = "Time (days)",
  ylab = "Overall survival probability",
  legend.title = "Mutation group",
  legend.labs = c("HR<1", "HR>1", "Both"),
  legend = "right",
  title = "Survival analysis among mutant patients"
)

## 8. Main survival plot: square panel with Firth HR in the lower-left corner
p_mut_only_km$plot <- p_mut_only_km$plot +
  theme_bw() +
  theme(
    aspect.ratio = 1,
    plot.title = element_text(hjust = 0.5, size = 11),
    axis.text = element_text(size = 9, color = "black"),
    axis.title = element_text(size = 10),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8)
  ) +
  annotate(
    "text",
    x = max(surv_mut_only_df$OS_time, na.rm = TRUE) * 0.04,
    y = 0.07,
    label = hr_text,
    hjust = 0,
    vjust = 0,
    size = 1.9
  )

## 9. Format the risk table
p_mut_only_km$table <- p_mut_only_km$table +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 8, color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),
    axis.title.x = element_text(size = 9),
    axis.title.y = element_text(size = 9)
  )

## 10. Output the plot
p_mut_only_km




##############################################
#######Figure5B
##############################################





suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(survival)
  library(ggplot2)
})

############################################################
## Apply fixed LASSO/GEO coefficients to TCGA-BRCA
##
## Risk score = sum(coef_i * z-scored log2(FPKM + 1)_i)
## Z-score is calculated within the TCGA-BRCA cohort.
############################################################

tcga_dir <- "E:/code/400T/TCGA/TCGA-BRCA"
expr_file <- file.path(tcga_dir, "TCGA-BRCA.star_fpkm.tsv.gz")
surv_file <- file.path(tcga_dir, "TCGA-BRCA.survival.tsv.gz")
clinical_file <- file.path(tcga_dir, "TCGA-BRCA.clinical.tsv.gz")

outdir <- "E:/code/400T/TCGA_fixed_risk_score_results"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

## Coefficients from the previous LASSO result.
## These are fixed external coefficients and are not refitted in TCGA.
risk_coef <- c(
  ARID1B =  1.16870213,
  JAK1   = -0.31956337,
  MUC4   =  0.31647725,
  CFB    = -0.24106792,
  PDCD1  = -0.22657677,
  NFKBIA = -0.22434813,
  FCGBP  = -0.07262415
)

genes <- names(risk_coef)

############################################################
## 1. Read TCGA-BRCA FPKM expression matrix
############################################################

expr_raw <- fread(expr_file, data.table = FALSE)
colnames(expr_raw)[1] <- "Gene"

expr_raw <- expr_raw %>%
  mutate(
    Ensembl_ID = as.character(Gene),
    Ensembl_ID = str_replace(Ensembl_ID, "\\..*$", "")
  )

all_sample_cols <- colnames(expr_raw)[str_detect(colnames(expr_raw), "^TCGA-")]
tumor_sample_cols <- all_sample_cols[substr(all_sample_cols, 14, 15) == "01"]

cat("Expression total samples: ", length(all_sample_cols), "\n", sep = "")
cat("Expression primary tumor samples: ", length(tumor_sample_cols), "\n", sep = "")

############################################################
## 2. Map Ensembl IDs to gene symbols without org.Hs.eg.db
############################################################

## Ensembl GRCh38 IDs for the 7 fixed model genes.
## Version suffixes in the expression matrix are removed above.
gene_map <- data.frame(
  Symbol = c("ARID1B", "JAK1", "MUC4", "CFB", "PDCD1", "NFKBIA", "FCGBP"),
  Ensembl_ID = c(
    "ENSG00000049618",
    "ENSG00000162434",
    "ENSG00000145113",
    "ENSG00000243649",
    "ENSG00000188389",
    "ENSG00000100906",
    "ENSG00000275395"
  ),
  stringsAsFactors = FALSE
)

expr_raw <- expr_raw %>%
  left_join(gene_map, by = "Ensembl_ID")

expr_gene <- expr_raw %>%
  filter(Symbol %in% genes) %>%
  dplyr::select(Symbol, all_of(tumor_sample_cols))

matched_genes <- sort(unique(expr_gene$Symbol))
missing_genes <- setdiff(genes, matched_genes)

cat("Matched genes:\n")
print(matched_genes)
cat("Missing genes:\n")
print(missing_genes)

if (length(matched_genes) != length(genes)) {
  stop("Not all risk-score genes were matched in TCGA expression matrix.")
}

expr_gene[, tumor_sample_cols] <- lapply(expr_gene[, tumor_sample_cols], as.numeric)

## If multiple Ensembl IDs map to the same symbol, use mean expression.
expr_mat <- expr_gene %>%
  group_by(Symbol) %>%
  summarise(
    across(all_of(tumor_sample_cols), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  column_to_rownames("Symbol") %>%
  as.data.frame()

############################################################
## 3. Transpose expression and average repeated patient samples
############################################################

expr_t <- as.data.frame(t(expr_mat)) %>%
  rownames_to_column("sample_id") %>%
  mutate(
    sample_type = substr(sample_id, 14, 15),
    patient_id = substr(sample_id, 1, 12)
  ) %>%
  filter(sample_type == "01")

expr_t_avg <- expr_t %>%
  dplyr::select(patient_id, all_of(genes)) %>%
  group_by(patient_id) %>%
  summarise(
    across(all_of(genes), ~ mean(as.numeric(.x), na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    across(all_of(genes), ~ log2(.x + 1))
  )

############################################################
## 4. Read TCGA-BRCA survival data
############################################################

surv_raw <- fread(surv_file, data.table = FALSE)

sample_col <- if ("sample" %in% colnames(surv_raw)) {
  "sample"
} else if ("sample_id" %in% colnames(surv_raw)) {
  "sample_id"
} else {
  colnames(surv_raw)[1]
}

surv_df <- surv_raw %>%
  mutate(
    sample_id = as.character(.data[[sample_col]]),
    sample_type = substr(sample_id, 14, 15),
    patient_id = substr(sample_id, 1, 12),
    OS_time = as.numeric(OS.time),
    OS_event = as.numeric(OS)
  ) %>%
  filter(
    sample_type == "01",
    !is.na(OS_time),
    !is.na(OS_event),
    OS_time > 0
  ) %>%
  dplyr::select(patient_id, OS_time, OS_event) %>%
  distinct(patient_id, .keep_all = TRUE)

clinical_raw <- fread(clinical_file, data.table = FALSE)

clean_unknown <- function(x) {
  x <- as.character(x)
  x[x %in% c("", "NA", "not reported", "Not Reported", "Not Available",
             "not available", "Unknown", "unknown", "--")] <- NA
  x
}

clinical_df <- clinical_raw %>%
  mutate(
    sample_id = as.character(sample),
    sample_type = substr(sample_id, 14, 15),
    patient_id = substr(sample_id, 1, 12),
    age = as.numeric(clean_unknown(age_at_index.demographic)),
    age_group = case_when(
      !is.na(age) & age < 60 ~ "<60",
      !is.na(age) & age >= 60 ~ ">=60",
      TRUE ~ NA_character_
    ),
    gender = str_to_title(clean_unknown(gender.demographic)),
    race_raw = str_to_lower(clean_unknown(race.demographic)),
    race = case_when(
      str_detect(race_raw, "white") ~ "White",
      str_detect(race_raw, "black") ~ "Black",
      str_detect(race_raw, "asian") ~ "Asian",
      !is.na(race_raw) ~ "Other",
      TRUE ~ NA_character_
    ),
    stage_raw = str_to_upper(clean_unknown(ajcc_pathologic_stage.diagnoses)),
    stage = case_when(
      str_detect(stage_raw, "STAGE IV") ~ "Stage IV",
      str_detect(stage_raw, "STAGE III") ~ "Stage III",
      str_detect(stage_raw, "STAGE II") ~ "Stage II",
      str_detect(stage_raw, "STAGE I") ~ "Stage I",
      TRUE ~ NA_character_
    ),
    T_raw = str_to_upper(clean_unknown(ajcc_pathologic_t.diagnoses)),
    N_raw = str_to_upper(clean_unknown(ajcc_pathologic_n.diagnoses)),
    M_raw = str_to_upper(clean_unknown(ajcc_pathologic_m.diagnoses)),
    T_stage = str_extract(T_raw, "T[0-4]"),
    N_stage = str_extract(N_raw, "N[0-3]"),
    M_stage = str_extract(M_raw, "M[0-1]")
  ) %>%
  filter(sample_type == "01") %>%
  dplyr::select(patient_id, age, age_group, gender, race, stage, T_stage, N_stage, M_stage) %>%
  distinct(patient_id, .keep_all = TRUE) %>%
  mutate(
    age_group = factor(age_group, levels = c("<60", ">=60")),
    gender = factor(gender),
    race = factor(race, levels = c("White", "Black", "Asian", "Other")),
    stage = factor(stage, levels = c("Stage I", "Stage II", "Stage III", "Stage IV")),
    T_stage = factor(T_stage, levels = c("T1", "T2", "T3", "T4")),
    N_stage = factor(N_stage, levels = c("N0", "N1", "N2", "N3")),
    M_stage = factor(M_stage, levels = c("M0", "M1"))
  )

############################################################
## 5. Merge expression and survival
############################################################

tcga_df <- inner_join(surv_df, expr_t_avg, by = "patient_id") %>%
  left_join(clinical_df, by = "patient_id")

cat("TCGA patients with expression and OS: ", nrow(tcga_df), "\n", sep = "")

############################################################
## 6. Calculate fixed-coefficient risk score
############################################################

z_mat <- scale(as.matrix(tcga_df[, genes, drop = FALSE]))
colnames(z_mat) <- genes

tcga_df$risk_score <- as.vector(z_mat[, genes, drop = FALSE] %*% risk_coef[genes])
tcga_df$risk_group <- ifelse(
  tcga_df$risk_score >= median(tcga_df$risk_score, na.rm = TRUE),
  "High risk",
  "Low risk"
)
tcga_df$risk_group <- factor(tcga_df$risk_group, levels = c("Low risk", "High risk"))

write.csv(
  tcga_df,
  file.path(outdir, "TCGA_BRCA_fixed_LASSO_risk_score_data.csv"),
  row.names = FALSE
)

############################################################
## 7. Cox model: fixed risk score and risk group
############################################################

cox_score <- coxph(
  Surv(OS_time, OS_event) ~ risk_score,
  data = tcga_df
)

cox_group <- coxph(
  Surv(OS_time, OS_event) ~ risk_group,
  data = tcga_df
)

cox_score_sum <- summary(cox_score)
cox_group_sum <- summary(cox_group)

score_res <- data.frame(
  model = "Continuous fixed risk score",
  variable = "risk_score",
  HR = cox_score_sum$conf.int[1, "exp(coef)"],
  lower95 = cox_score_sum$conf.int[1, "lower .95"],
  upper95 = cox_score_sum$conf.int[1, "upper .95"],
  pvalue = cox_score_sum$coefficients[1, "Pr(>|z|)"]
)

group_res <- data.frame(
  model = "Median risk group",
  variable = rownames(cox_group_sum$coefficients)[1],
  HR = cox_group_sum$conf.int[1, "exp(coef)"],
  lower95 = cox_group_sum$conf.int[1, "lower .95"],
  upper95 = cox_group_sum$conf.int[1, "upper .95"],
  pvalue = cox_group_sum$coefficients[1, "Pr(>|z|)"]
)

write.csv(
  rbind(score_res, group_res),
  file.path(outdir, "TCGA_BRCA_fixed_risk_score_cox_results.csv"),
  row.names = FALSE
)

############################################################
## 8. Multivariate Cox with the 7 model genes in TCGA
############################################################

tcga_gene_z_df <- cbind(
  tcga_df[, c("patient_id", "OS_time", "OS_event")],
  as.data.frame(z_mat[, genes, drop = FALSE])
)

multi_formula <- as.formula(
  paste0("Surv(OS_time, OS_event) ~ ", paste(genes, collapse = " + "))
)

multi_gene_fit <- coxph(
  multi_formula,
  data = tcga_gene_z_df
)

multi_sum <- summary(multi_gene_fit)
multi_gene_res <- data.frame(
  gene = rownames(multi_sum$coefficients),
  TCGA_multivariate_HR = multi_sum$conf.int[, "exp(coef)"],
  lower95 = multi_sum$conf.int[, "lower .95"],
  upper95 = multi_sum$conf.int[, "upper .95"],
  pvalue = multi_sum$coefficients[, "Pr(>|z|)"],
  fixed_LASSO_coef = risk_coef[rownames(multi_sum$coefficients)],
  row.names = NULL
)

write.csv(
  multi_gene_res,
  file.path(outdir, "TCGA_BRCA_7_gene_multivariate_cox_results.csv"),
  row.names = FALSE
)

sink(file.path(outdir, "TCGA_BRCA_fixed_risk_score_summary.txt"))
cat("Risk score formula:\n")
print(risk_coef)
cat("\nSample number:\n")
print(nrow(tcga_df))
cat("\nContinuous risk score Cox:\n")
print(summary(cox_score))
cat("\nMedian risk group Cox:\n")
print(summary(cox_group))
cat("\n7-gene multivariate Cox in TCGA:\n")
print(summary(multi_gene_fit))
sink()

############################################################
## 9. Survival Cox forest plots, not Kaplan-Meier curves
############################################################

format_p <- function(p) {
  ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
}

extract_cox_rows <- function(fit, model_name, variable_label = NULL) {
  fit_sum <- summary(fit)
  variable <- if (is.null(variable_label)) {
    rownames(fit_sum$coefficients)
  } else {
    rep(variable_label, nrow(fit_sum$coefficients))
  }
  data.frame(
    model = model_name,
    term = rownames(fit_sum$coefficients),
    variable = variable,
    HR = fit_sum$conf.int[, "exp(coef)"],
    lower95 = fit_sum$conf.int[, "lower .95"],
    upper95 = fit_sum$conf.int[, "upper .95"],
    pvalue = fit_sum$coefficients[, "Pr(>|z|)"],
    row.names = NULL
  )
}

term_label <- function(term) {
  term <- str_replace(term, "^risk_score$", "Risk score, continuous")
  term <- str_replace(term, "^risk_groupHigh risk$", "High vs low risk group")
  term <- str_replace(term, "^age_group>=60$", "Age: >=60 vs <60")
  term <- str_replace(term, "^gender", "Gender: ")
  term <- str_replace(term, "^race", "Race: ")
  term <- str_replace(term, "^stage", "Pathologic stage: ")
  term <- str_replace(term, "^T_stage", "Pathologic T: ")
  term <- str_replace(term, "^N_stage", "Pathologic N: ")
  term <- str_replace(term, "^M_stage", "Pathologic M: ")
  term
}

term_order <- function(term) {
  case_when(
    term == "risk_score" ~ 1,
    term == "risk_groupHigh risk" ~ 2,
    term == "age_group>=60" ~ 3,
    str_detect(term, "^gender") ~ 4,
    str_detect(term, "^race") ~ 5,
    str_detect(term, "^stage") ~ 6,
    str_detect(term, "^T_stage") ~ 7,
    str_detect(term, "^N_stage") ~ 8,
    str_detect(term, "^M_stage") ~ 9,
    TRUE ~ 99
  )
}

fit_univariate_cox <- function(data, var, label) {
  tmp <- data %>%
    dplyr::select(OS_time, OS_event, all_of(var)) %>%
    filter(!is.na(.data[[var]]))
  
  if (nrow(tmp) < 20 || sum(tmp$OS_event, na.rm = TRUE) < 5) {
    return(NULL)
  }
  if (is.factor(tmp[[var]]) || is.character(tmp[[var]])) {
    tmp[[var]] <- droplevels(factor(tmp[[var]]))
    if (nlevels(tmp[[var]]) < 2) {
      return(NULL)
    }
  }
  
  fit <- coxph(as.formula(paste0("Surv(OS_time, OS_event) ~ ", var)), data = tmp)
  extract_cox_rows(fit, model_name = "Univariate Cox", variable_label = label)
}

plot_forest <- function(plot_df, title, subtitle, outfile, width = 7.2, height = 4.8) {
  plot_df <- plot_df %>%
    filter(is.finite(HR), is.finite(lower95), is.finite(upper95)) %>%
    mutate(
      label = term_label(term),
      plot_order = term_order(term),
      sig = case_when(
        pvalue < 0.1 & HR > 1 ~ "Risk",
        pvalue < 0.1 & HR < 1 ~ "Protective",
        TRUE ~ "NS"
      )
    ) %>%
    arrange(plot_order, term) %>%
    mutate(
      label = factor(label, levels = rev(label)),
      sig = factor(sig, levels = c("Risk", "Protective", "NS"))
    )
  
  p <- ggplot(
    plot_df,
    aes(
      x = HR,
      y = label,
      color = sig
    )
  ) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey45", linewidth = 0.4) +
    geom_errorbar(
      aes(xmin = lower95, xmax = upper95),
      orientation = "y",
      width = 0.18,
      linewidth = 0.55
    ) +
    geom_point(size = 2.5) +
    scale_x_log10() +
    scale_color_manual(
      values = c(
        "Risk" = "#E64B35",
        "Protective" = "#4DBBD5",
        "NS" = "#7A7A7A"
      )
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Hazard ratio (log scale)",
      y = NULL,
      color = NULL
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey30"),
      axis.text.y = element_text(size = 8, color = "black"),
      axis.text.x = element_text(color = "black"),
      legend.position = "right"
    )
  
  pdf(file.path(outdir, outfile), width = width, height = height)
  print(p)
  dev.off()
}

univariate_vars <- c(
  risk_group = "Risk group",
  age_group = "Age group",
  gender = "Gender",
  stage = "Pathologic stage"
)

univariate_clinical_res <- bind_rows(
  lapply(
    names(univariate_vars),
    function(v) fit_univariate_cox(tcga_df, v, univariate_vars[[v]])
  )
)

write.csv(
  univariate_clinical_res,
  file.path(outdir, "TCGA_BRCA_survival_univariate_clinical_cox_results.csv"),
  row.names = FALSE
)

plot_forest(
  univariate_clinical_res,
  "Univariate Cox regression",
  "TCGA-BRCA tumor samples only; fixed risk score and clinical variables",
  "TCGA_BRCA_survival_univariate_clinical_Cox_forest.pdf",
  width = 8.8,
  height = max(4.5, 0.35 * nrow(univariate_clinical_res) + 1.6)
)

multivariable_vars <- c("risk_group", "age_group", "gender", "stage")
multivariable_vars <- multivariable_vars[sapply(multivariable_vars, function(v) {
  tmp <- tcga_df[[v]]
  tmp <- tmp[!is.na(tmp)]
  if (is.factor(tmp) || is.character(tmp)) {
    length(unique(tmp)) >= 2
  } else {
    length(tmp) >= 20
  }
})]

multivariable_df <- tcga_df %>%
  dplyr::select(OS_time, OS_event, all_of(multivariable_vars)) %>%
  filter(complete.cases(.))

multivariable_formula <- as.formula(
  paste0("Surv(OS_time, OS_event) ~ ", paste(multivariable_vars, collapse = " + "))
)

multivariable_fit <- coxph(multivariable_formula, data = multivariable_df)
multivariable_clinical_res <- extract_cox_rows(
  multivariable_fit,
  model_name = "Multivariable Cox"
)

write.csv(
  multivariable_clinical_res,
  file.path(outdir, "TCGA_BRCA_survival_multivariable_clinical_cox_results.csv"),
  row.names = FALSE
)

plot_forest(
  multivariable_clinical_res,
  "Multivariable Cox regression",
  "TCGA-BRCA tumor samples only; adjusted for fixed risk group and clinical variables",
  "TCGA_BRCA_survival_multivariable_clinical_Cox_forest.pdf",
  width = 8.8,
  height = max(4.2, 0.35 * nrow(multivariable_clinical_res) + 1.6)
)

cat("Done. Results saved to: ", outdir, "\n", sep = "")