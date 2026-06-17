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
