##############################################
#######FigureS5A
##############################################



lasso_genes <- univ_cox_res %>%
  filter(!is.na(pvalue), pvalue < 0.1) %>%
  pull(gene)

if (length(lasso_genes) < 2) {
  lasso_genes <- genes_use
}

cat("Genes used for LASSO:\n")
print(lasso_genes)

x <- as.matrix(cox_df[, lasso_genes, drop = FALSE])
y <- Surv(cox_df$OS_time, cox_df$OS_event)

keep <- complete.cases(x) &
  !is.na(cox_df$OS_time) &
  !is.na(cox_df$OS_event)

x <- x[keep, , drop = FALSE]
y <- y[keep]
cox_lasso_df <- cox_df[keep, ]

set.seed(123)

cvfit <- cv.glmnet(
  x = x,
  y = y,
  family = "cox",
  alpha = 1,
  nfolds = 10
)

pdf(file.path(outdir, "03_LASSO_Cox_cvfit.pdf"), width = 6, height = 5)
plot(cvfit)
dev.off()

png(file.path(outdir, "03_LASSO_Cox_cvfit.png"), width = 1800, height = 1500, res = 300)
plot(cvfit)
dev.off()

coef_min <- coef(cvfit, s = "lambda.min")

lasso_res <- data.frame(
  gene = rownames(coef_min),
  coef = as.numeric(coef_min)
) %>%
  filter(coef != 0) %>%
  arrange(desc(abs(coef)))

coef_1se <- coef(cvfit, s = "lambda.1se")

lasso_res_1se <- data.frame(
  gene = rownames(coef_1se),
  coef = as.numeric(coef_1se)
) %>%
  filter(coef != 0) %>%
  arrange(desc(abs(coef)))

write.csv(lasso_res, file.path(outdir, "04_LASSO_selected_genes_lambda_min.csv"), row.names = FALSE)
write.csv(lasso_res_1se, file.path(outdir, "04_LASSO_selected_genes_lambda_1se.csv"), row.names = FALSE)

print(lasso_res)
print(lasso_res_1se)

if (nrow(lasso_res) > 0) {
  final_lasso_res <- lasso_res
  final_lambda <- "lambda.min"
} else if (nrow(lasso_res_1se) > 0) {
  final_lasso_res <- lasso_res_1se
  final_lambda <- "lambda.1se"
} else {
  stop("LASSO 没有筛选到任何基因。")
}

final_genes <- final_lasso_res$gene
final_coef <- final_lasso_res$coef



##############################################
#######FigureS5B
##############################################
############################################################
#LASSO coefficient plot
############################################################

p_lasso_coef <- final_lasso_res %>%
  mutate(
    gene = reorder(gene, coef),
    direction = ifelse(coef > 0, "Risk", "Protective")
  ) %>%
  ggplot(aes(x = gene, y = coef, fill = direction)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.3) +
  coord_flip() +
  scale_fill_manual(
    values = c(
      "Risk" = "#E64B35",
      "Protective" = "#4DBBD5"
    )
  ) +
  labs(
    x = NULL,
    y = "LASSO coefficient",
    title = paste0("LASSO Cox selected genes (", final_lambda, ")")
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text = element_text(color = "black", face = "bold")
  )

ggsave(file.path(outdir, "05_LASSO_selected_gene_coefficients.pdf"), p_lasso_coef, width = 5.5, height = 4)
ggsave(file.path(outdir, "05_LASSO_selected_gene_coefficients.png"), p_lasso_coef, width = 5.5, height = 4, dpi = 300)







##############################################
#######FigureS5C
##############################################

############################################################
## 2. Mutant-only KM
############################################################

x_max_mut_only <- ceiling(
  max(surv_mut_only_df$OS_time, na.rm = TRUE) / x_break
) * x_break

p_mut_only_km <- ggsurvplot(
  fit_mut_only_km,
  data = surv_mut_only_df,
  pval = TRUE,
  conf.int = FALSE,
  risk.table = TRUE,
  risk.table.height = 0.25,
  risk.table.fontsize = 4,
  palette = c(
    "#4DBBD5",
    "#E64B35",
    "#8E44AD"
  ),
  xlim = c(0, x_max_mut_only),
  break.time.by = x_break,
  xlab = "Time (days)",
  ylab = "Overall survival probability",
  legend = "bottom",
  legend.title = "Mutation group",
  legend.labs = c(
    "HR<1",
    "HR>1",
    "Both"
  ),
  title = "Survival analysis among mutant patients",
  ggtheme = theme_bw(),
  tables.theme = theme_bw()
)

p_mut_only_km$plot <- p_mut_only_km$plot +
  theme_bw() %+replace%
  theme(
    aspect.ratio = 1,
    plot.title = element_text(
      hjust = 0.5,
      size = 11,
      face = "bold"
    ),
    axis.text = element_text(
      size = 9,
      color = "black"
    ),
    axis.title = element_text(
      size = 10,
      color = "black"
    ),
    legend.title = element_text(
      size = 9,
      color = "black"
    ),
    legend.text = element_text(
      size = 8,
      color = "black"
    ),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.6
    ),
    panel.grid.major = element_line(
      color = "grey90",
      linewidth = 0.3
    ),
    panel.grid.minor = element_line(
      color = "grey92",
      linewidth = 0.2
    )
  )

p_mut_only_km$table <- p_mut_only_km$table +
  theme_bw() %+replace%
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
      size = 9,
      color = "black"
    ),
    axis.title.y = element_text(
      size = 9,
      color = "black"
    ),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.6
    ),
    panel.grid.major = element_line(
      color = "grey90",
      linewidth = 0.3
    ),
    panel.grid.minor = element_blank()
  )

p_mut_only_km

pdf(
  file = "E:/code/400T/LASOO_COX/TCGA_BRCA_tumor_only_Cox_LASSO_mutation/mut_only_KM_square.pdf",
  width = 5.5,
  height = 6.8
)

print(p_mut_only_km)

dev.off()



##############################################
#######FigureS5D
##############################################

library(survival)
library(survminer)
library(ggplot2)
library(dplyr)
library(coxphf)

## Separate KM plots by mutation variant type, with Wildtype reference and HR labels:
## 1) Wildtype + HR<1 mutant only variant types
## 2) Wildtype + HR>1 mutant only variant types
##
## Required object in environment:
##   surv_mut_integrated_df with patient_id, OS_time, OS_event, mutation_integrated_group
##
## MAF-like mutation file with Variant_Classification:
##   E:/code/400T/TCGA/TCGA_maf_mutation.tsv

maf_file <- "E:/code/400T/TCGA/TCGA_maf_mutation.tsv"
mutation_type_col <- "Variant_Classification"
patient_id_col <- "patient_id"
sample_col <- "Tumor_Sample_Barcode"
gene_col <- "Hugo_Symbol"

standardize_mutation_group <- function(x) {
  x_chr <- trimws(as.character(x))
  out <- x_chr
  out[x_chr %in% c("Wildtype", "WT", "wildtype")] <- "Wildtype"
  out[x_chr %in% c("HR<1 mutant only", "HR<1", "Protective", "Protective mutant only")] <- "HR<1 mutant only"
  out[x_chr %in% c("HR>1 mutant only", "HR>1", "Risk", "Risk mutant only")] <- "HR>1 mutant only"
  out[x_chr %in% c("Both mutant", "Both", "Dual-effect", "Dual effect", "Dual-effect mutant")] <- "Both mutant"
  factor(out, levels = c("Wildtype", "HR<1 mutant only", "HR>1 mutant only", "Both mutant"))
}

format_pvalue <- function(pvalue) {
  ifelse(is.na(pvalue), "NA", ifelse(pvalue < 0.001, "<0.001", sprintf("%.3f", pvalue)))
}

read_mutation_type_df <- function(maf_file) {
  maf_raw <- read.delim(
    maf_file,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  required_cols <- c(sample_col, gene_col, mutation_type_col)
  missing_cols <- setdiff(required_cols, colnames(maf_raw))
  if (length(missing_cols) > 0) {
    stop(paste0(
      "Missing columns in mutation file: ", paste(missing_cols, collapse = ", "),
      "\nAvailable columns: ", paste(colnames(maf_raw), collapse = ", ")
    ))
  }
  
  maf_raw %>%
    mutate(
      Tumor_Sample_Barcode = as.character(.data[[sample_col]]),
      patient_id = substr(Tumor_Sample_Barcode, 1, 12),
      sample_type = substr(Tumor_Sample_Barcode, 14, 15),
      Gene = as.character(.data[[gene_col]]),
      mutation_type = as.character(.data[[mutation_type_col]])
    ) %>%
    filter(
      sample_type == "01",
      !is.na(patient_id), patient_id != "",
      !is.na(Gene), Gene != "",
      !is.na(mutation_type), mutation_type != ""
    ) %>%
    distinct(patient_id, Gene, mutation_type)
}

get_target_genes <- function(target_group) {
  if (target_group == "HR<1 mutant only" && exists("hr_lt1_genes")) return(hr_lt1_genes)
  if (target_group == "HR>1 mutant only" && exists("hr_gt1_genes")) return(hr_gt1_genes)
  NULL
}

collapse_patient_mutation_type <- function(x) {
  x <- sort(unique(na.omit(as.character(x))))
  if (length(x) == 0) return(NA_character_)
  if (length(x) == 1) return(x)
  "Multiple"
}

get_survival_source_df <- function() {
  if (exists("surv_mut_integrated_df")) return(surv_mut_integrated_df)
  stop(
    paste0(
      "surv_mut_integrated_df not found. Wildtype is not available in surv_mut_only_df. ",
      "Please run the upstream code that creates surv_mut_integrated_df first."
    )
  )
}

build_variant_type_surv_df <- function(data, target_group, mut_type_df) {
  data <- data %>%
    mutate(
      mutation_integrated_group_raw = as.character(mutation_integrated_group),
      mutation_integrated_group = standardize_mutation_group(mutation_integrated_group)
    )
  
  wildtype_df <- data %>%
    filter(mutation_integrated_group == "Wildtype") %>%
    filter(!is.na(OS_time), !is.na(OS_event)) %>%
    mutate(mutation_type = "Wildtype")
  
  mutant_base_df <- data %>%
    filter(mutation_integrated_group == target_group) %>%
    filter(!is.na(OS_time), !is.na(OS_event))
  
  target_genes <- get_target_genes(target_group)
  mut_df_use <- mut_type_df %>%
    filter(patient_id %in% mutant_base_df[[patient_id_col]])
  
  if (!is.null(target_genes)) {
    mut_df_use <- mut_df_use %>% filter(Gene %in% target_genes)
  }
  
  patient_type_df <- mut_df_use %>%
    group_by(patient_id) %>%
    summarise(
      mutation_type = collapse_patient_mutation_type(mutation_type),
      .groups = "drop"
    ) %>%
    filter(!is.na(mutation_type))
  
  mutant_df <- mutant_base_df %>%
    inner_join(patient_type_df, by = "patient_id")
  
  bind_rows(wildtype_df, mutant_df) %>%
    mutate(
      mutation_type = factor(
        mutation_type,
        levels = c("Wildtype", setdiff(sort(unique(mutation_type)), "Wildtype"))
      )
    )
}

get_hr_vs_wildtype <- function(df) {
  type_levels <- setdiff(levels(droplevels(df$mutation_type)), "Wildtype")
  
  if (length(type_levels) == 0) return(data.frame())
  
  do.call(rbind, lapply(type_levels, function(tp) {
    tmp_df <- df %>%
      filter(mutation_type %in% c("Wildtype", tp)) %>%
      mutate(mutation_type_for_model = factor(mutation_type, levels = c("Wildtype", tp)))
    
    n_by_group <- table(tmp_df$mutation_type_for_model)
    event_by_group <- tapply(tmp_df$OS_event, tmp_df$mutation_type_for_model, function(x) sum(x == 1, na.rm = TRUE))
    
    if (length(n_by_group) < 2 || any(n_by_group == 0) || sum(event_by_group, na.rm = TRUE) == 0) {
      return(data.frame(
        comparison = paste0(tp, " vs Wildtype"),
        HR = NA_real_, lower95 = NA_real_, upper95 = NA_real_, pvalue = NA_real_
      ))
    }
    
    fit <- tryCatch(
      coxphf(
        Surv(OS_time, OS_event) ~ mutation_type_for_model,
        data = tmp_df,
        maxit = 1000,
        maxstep = 0.25
      ),
      error = function(e) NULL
    )
    
    if (is.null(fit)) {
      return(data.frame(
        comparison = paste0(tp, " vs Wildtype"),
        HR = NA_real_, lower95 = NA_real_, upper95 = NA_real_, pvalue = NA_real_
      ))
    }
    
    data.frame(
      comparison = paste0(tp, " vs Wildtype"),
      HR = exp(fit$coefficients),
      lower95 = fit$ci.lower,
      upper95 = fit$ci.upper,
      pvalue = fit$prob,
      row.names = NULL
    )
  }))
}

build_hr_text <- function(hr_df) {
  if (nrow(hr_df) == 0) return("Ref: Wildtype")
  
  paste0(
    "Ref: Wildtype\n",
    paste0(
      hr_df$comparison,
      ": HR=", sprintf("%.2f", hr_df$HR),
      ", P=", format_pvalue(hr_df$pvalue),
      collapse = "\n"
    )
  )
}

plot_km_by_variant_type <- function(data, target_group, mut_type_df) {
  df <- build_variant_type_surv_df(data, target_group, mut_type_df)
  
  if (nrow(df) == 0) {
    stop(paste0("No available rows for ", target_group, ". Check patient_id matching and mutation type data."))
  }
  
  cat("\nWildtype + ", target_group, "\n", sep = "")
  print(
    aggregate(
      x = list(OS_event = df$OS_event),
      by = list(Mutation_type = df$mutation_type),
      FUN = function(x) c(
        n = length(x),
        event = sum(x == 1, na.rm = TRUE),
        censored = sum(x == 0, na.rm = TRUE)
      )
    )
  )
  
  hr_df <- get_hr_vs_wildtype(df)
  print(hr_df)
  hr_text <- build_hr_text(hr_df)
  
  fit <- survfit(Surv(OS_time, OS_event) ~ mutation_type, data = df)
  
  p <- ggsurvplot(
    fit,
    data = df,
    pval = TRUE,
    conf.int = FALSE,
    risk.table = TRUE,
    risk.table.height = 0.20,
    surv.plot.height = 0.72,
    legend = "right",
    legend.labs = levels(df$mutation_type),
    xlab = "Time (days)",
    ylab = "Overall survival probability",
    legend.title = "Mutation type",
    title = paste0("Wildtype + ", target_group),
    palette = "Dark2"
  )
  
  p$plot <- p$plot +
    theme_bw() +
    theme(
      aspect.ratio = 1,
      plot.title = element_text(hjust = 0.5, size = 11),
      axis.text = element_text(size = 10, color = "black"),
      axis.title = element_text(size = 12),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9)
    ) +
    annotate(
      "text",
      x = max(df$OS_time, na.rm = TRUE) * 0.04,
      y = 0.07,
      label = hr_text,
      hjust = 0,
      vjust = 0,
      size = 2.4
    )
  
  p$table <- p$table +
    scale_y_discrete(labels = function(x) gsub("^mutation_type=", "", x)) +
    theme_bw() +
    theme(
      axis.text.x = element_text(size = 9, color = "black"),
      axis.text.y = element_text(size = 9, color = "black"),
      axis.title.x = element_text(size = 10),
      axis.title.y = element_text(size = 10)
    )
  
  p
}

survival_source_df <- get_survival_source_df()
survival_source_df$mutation_integrated_group <- standardize_mutation_group(
  survival_source_df$mutation_integrated_group
)

cat("Available mutation_integrated_group values after standardization:\n")
print(table(survival_source_df$mutation_integrated_group, useNA = "ifany"))

mut_type_df <- read_mutation_type_df(maf_file)
cat("Mutation type counts in MAF file:\n")
print(sort(table(mut_type_df$mutation_type), decreasing = TRUE))

p_hr_lt1_variant_type <- plot_km_by_variant_type(
  survival_source_df,
  target_group = "HR<1 mutant only",
  mut_type_df = mut_type_df
)

p_hr_gt1_variant_type <- plot_km_by_variant_type(
  survival_source_df,
  target_group = "HR>1 mutant only",
  mut_type_df = mut_type_df
)

p_hr_lt1_variant_type
p_hr_gt1_variant_type

combined_variant_type_km <- arrange_ggsurvplots(
  list(p_hr_lt1_variant_type, p_hr_gt1_variant_type),
  ncol = 2,
  nrow = 1,
  risk.table.height = 0.25
)

combined_variant_type_km

## Save two separate PDF files
pdf(
  file = "E:/code/400T/TCGA/HR_lt1_mutant_only_variant_type_KM_with_wildtype_HR.pdf",
  width = 10,
  height = 8
)
print(p_hr_lt1_variant_type)
dev.off()

pdf(
  file = "E:/code/400T/TCGA/HR_gt1_mutant_only_variant_type_KM_with_wildtype_HR.pdf",
  width = 10,
  height = 8
)
print(p_hr_gt1_variant_type)
dev.off()

## Optional combined export
## pdf(
##   file = "E:/code/400T/TCGA/combined_variant_type_KM_by_HR_group_with_wildtype_HR.pdf",
##   width = 14,
##   height = 8
## )
## print(combined_variant_type_km)
## dev.off()