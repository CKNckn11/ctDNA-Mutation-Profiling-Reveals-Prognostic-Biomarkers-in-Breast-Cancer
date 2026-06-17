library(dplyr)
library(ggplot2)
library(cowplot)

load("E:/code/400T/all.RData")

df <- clinical_df_all2 %>%
  distinct(clinical_id, .keep_all = TRUE) %>%
  mutate(
    pT_table = ifelse(is.na(pT_raw) | pT_raw %in% c("", "/", "NA"), NA, paste0("T", pT_raw)),
    pT_table = factor(pT_table, levels = c("T1", "T2", "T3", "T4")),
    cN_table = factor(cN_group, levels = c("negative", "positive"), labels = c("Negative", "Positive")),
    Stage_table = factor(Stage_group, levels = c("I", "II", "III", "IV")),
    Subtype_table = yunzhong_subtype,
    Subtype_table = ifelse(Subtype_table == "Luminal", NA, Subtype_table),
    Subtype_table = factor(Subtype_table, levels = c("Luminal A", "Luminal B", "HER2-enriched", "Triple Negative")),
    log_mut_count = log2(MAF_Mutation_Count + 1),
    log2_pre_library_concentration = ifelse(
      is.na(log2_pre_library_concentration) & !is.na(pre_library_concentration),
      log2(pre_library_concentration + 1),
      log2_pre_library_concentration
    )
  )

out_dir <- "E:/code/400T/TCGA/Figure"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

clinical_cols <- c(
  "T1" = "#F3B56B",
  "T2" = "#E84A3C",
  "T3" = "#8B1E2D",
  "T4" = "#5B1020",
  "Negative" = "#67A9CF",
  "Positive" = "#D63B37",
  "cfDNA-low" = "#79C7BE",
  "cfDNA-high" = "#D92727"
)

subtype_cols <- c(
  "Healthy" = "#9E9E9E",
  "Luminal A" = "#5DAECC",
  "Luminal B" = "#D46A6A",
  "HER2-enriched" = "#B8B8B8",
  "Triple Negative" = "#D9A84E"
)

stage_cols <- c("I" = "#F8C991", "II" = "#F4B36A", "III" = "#E6423A", "IV" = "#8B1E2D")

fmt_p_num <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) "<0.001" else signif(p, 2)
}

bar_count_plot <- function(data, var, title, fill_values = clinical_cols) {
  x <- data %>%
    filter(!is.na(.data[[var]])) %>%
    count(.data[[var]], name = "n") %>%
    mutate(label = as.character(n))
  names(x)[1] <- "group"

  ggplot(x, aes(group, n, fill = group)) +
    geom_col(width = 0.68, color = "white", linewidth = 0.4) +
    geom_text(aes(label = label), vjust = -0.25, size = 3.0) +
    scale_fill_manual(values = fill_values, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
    labs(title = title, x = NULL, y = "Number") +
    theme_classic(base_size = 10) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 11),
      axis.text.x = element_text(angle = 35, hjust = 1),
      legend.position = "none",
      axis.line = element_line(linewidth = 0.6),
      axis.ticks = element_line(linewidth = 0.6)
    )
}

add_pairwise_brackets <- function(p, dd, xvar, yvar, ref_group, compare_groups) {
  ymax <- max(dd[[yvar]], na.rm = TRUE)
  ymin <- min(dd[[yvar]], na.rm = TRUE)
  step <- (ymax - ymin) * 0.12
  if (step == 0 || is.na(step)) step <- 0.2
  base_y <- ymax + step * 0.35
  tick <- step * 0.12
  levels_x <- levels(dd[[xvar]])

  for (i in seq_along(compare_groups)) {
    g <- compare_groups[i]
    ref_vals <- dd %>% filter(.data[[xvar]] == ref_group) %>% pull(.data[[yvar]])
    cmp_vals <- dd %>% filter(.data[[xvar]] == g) %>% pull(.data[[yvar]])
    pval <- tryCatch(wilcox.test(ref_vals, cmp_vals)$p.value, error = function(e) NA_real_)
    x1 <- which(levels_x == ref_group)
    x2 <- which(levels_x == g)
    y <- base_y + (i - 1) * step

    p <- p +
      annotate("segment", x = x1, xend = x2, y = y, yend = y, linewidth = 0.4) +
      annotate("segment", x = x1, xend = x1, y = y, yend = y - tick, linewidth = 0.4) +
      annotate("segment", x = x2, xend = x2, y = y, yend = y - tick, linewidth = 0.4) +
      annotate("text", x = (x1 + x2) / 2, y = y + tick, label = fmt_p_num(pval), size = 2.6)
  }

  p + coord_cartesian(ylim = c(ymin, base_y + length(compare_groups) * step + step * 0.4), clip = "off")
}

box_plot_ref_p <- function(data, xvar, yvar, title, ylab, fill_values, ref_group, compare_groups) {
  dd <- data %>% filter(!is.na(.data[[xvar]]), !is.na(.data[[yvar]]))

  p <- ggplot(dd, aes(.data[[xvar]], .data[[yvar]], fill = .data[[xvar]])) +
    geom_boxplot(width = 0.62, outlier.shape = NA, alpha = 0.72, color = "black", linewidth = 0.45) +
    geom_jitter(aes(color = .data[[xvar]]), width = 0.18, size = 1.15, alpha = 0.65, show.legend = FALSE) +
    scale_fill_manual(values = fill_values, drop = FALSE) +
    scale_color_manual(values = fill_values, drop = FALSE) +
    labs(title = title, x = NULL, y = ylab) +
    theme_classic(base_size = 10) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 11),
      axis.text.x = element_text(angle = 35, hjust = 1),
      legend.position = "none",
      axis.line = element_line(linewidth = 0.6),
      axis.ticks = element_line(linewidth = 0.6)
    )

  add_pairwise_brackets(p, dd, xvar, yvar, ref_group, compare_groups)
}

box_plot_left_ref_brackets <- function(data, xvar, yvar, title, ylab, fill_values, left_ref_group) {
  dd <- data %>% filter(!is.na(.data[[xvar]]), !is.na(.data[[yvar]]))
  dd[[xvar]] <- factor(dd[[xvar]], levels = names(fill_values))

  p <- ggplot(dd, aes(.data[[xvar]], .data[[yvar]], fill = .data[[xvar]])) +
    geom_boxplot(width = 0.62, outlier.shape = NA, alpha = 0.72, color = "black", linewidth = 0.45) +
    geom_jitter(aes(color = .data[[xvar]]), width = 0.18, size = 1.15, alpha = 0.65, show.legend = FALSE) +
    scale_fill_manual(values = fill_values, drop = FALSE) +
    scale_color_manual(values = fill_values, drop = FALSE) +
    labs(title = title, x = NULL, y = ylab) +
    theme_classic(base_size = 10) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 11),
      axis.text.x = element_text(angle = 35, hjust = 1),
      legend.position = "none",
      axis.line = element_line(linewidth = 0.6),
      axis.ticks = element_line(linewidth = 0.6)
    )

  levels_x <- levels(dd[[xvar]])
  ref_vals <- dd %>% filter(.data[[xvar]] == left_ref_group) %>% pull(.data[[yvar]])
  compare_groups <- rev(setdiff(levels_x, left_ref_group))
  y_positions <- c(7.80, 7.35, 6.90, 6.45)
  tick <- 0.08

  for (i in seq_along(compare_groups)) {
    g <- compare_groups[i]
    cmp_vals <- dd %>% filter(.data[[xvar]] == g) %>% pull(.data[[yvar]])
    pval <- tryCatch(wilcox.test(ref_vals, cmp_vals)$p.value, error = function(e) NA_real_)
    x1 <- which(levels_x == left_ref_group)
    x2 <- which(levels_x == g)
    y <- y_positions[i]

    p <- p +
      annotate("segment", x = x1, xend = x2, y = y, yend = y, linewidth = 0.45) +
      annotate("segment", x = x1, xend = x1, y = y, yend = y - tick, linewidth = 0.45) +
      annotate("segment", x = x2, xend = x2, y = y, yend = y - tick, linewidth = 0.45) +
      annotate("text", x = (x1 + x2) / 2, y = y + 0.05, label = fmt_p_num(pval), size = 2.5)
  }

  p + coord_cartesian(ylim = c(0, 8.3), clip = "off")
}

healthy_df <- read.csv(
  "C:/Users/cknckn11/Desktop/健康人样本检测结果健康人.csv",
  header = FALSE,
  stringsAsFactors = FALSE,
  fileEncoding = "GBK"
)

colnames(healthy_df) <- c("sample_id", "sample_no", "blank", "group", "pre_library_concentration", "qc")

healthy_conc_df <- healthy_df %>%
  filter(qc == "pass", !is.na(pre_library_concentration)) %>%
  transmute(
    Subtype_table = factor("Healthy", levels = names(subtype_cols)),
    log2_pre_library_concentration = log2(as.numeric(pre_library_concentration) + 1)
  )

tumor_conc_df <- df %>%
  filter(!is.na(Subtype_table), !is.na(log2_pre_library_concentration)) %>%
  transmute(
    Subtype_table = factor(as.character(Subtype_table), levels = names(subtype_cols)),
    log2_pre_library_concentration = log2_pre_library_concentration
  )

conc_plot_df <- bind_rows(healthy_conc_df, tumor_conc_df)
##############################################
#######Figure2A
##############################################
pA <- bar_count_plot(df, "pT_table", "Pathological tumor size", clinical_cols)
##############################################
#######Figure2B
##############################################
pB <- bar_count_plot(df, "cN_table", "Clinical nodal status", clinical_cols)
##############################################
#######Figure2C
##############################################
pC <- box_plot_ref_p(
  df,
  "Stage_table",
  "log_mut_count",
  "Mutation count by stage",
  "log2(mutation count + 1)",
  stage_cols,
  ref_group = "I",
  compare_groups = c("II", "III", "IV")
)
##############################################
#######Figure2D
##############################################
pD <- box_plot_ref_p(
  df,
  "Subtype_table",
  "Avg_VAF",
  "Mean VAF by subtype",
  "Mean VAF",
  subtype_cols[names(subtype_cols) != "Healthy"],
  ref_group = "Luminal A",
  compare_groups = c("Luminal B", "HER2-enriched", "Triple Negative")
)
##############################################
#######Figure2E
##############################################
pE <- box_plot_left_ref_brackets(
  conc_plot_df,
  "Subtype_table",
  "log2_pre_library_concentration",
  "Pre-library concentration by subtype",
  "log2(pre-library concentration + 1)",
  subtype_cols,
  left_ref_group = "Healthy"
)
fig <- plot_grid(
  pA, pB,
  pC, pD,
  pE,
  labels = LETTERS[1:6],
  label_fontface = "bold",
  label_size = 15,
  ncol = 2,
  align = "hv"
)

pdf(file.path(out_dir, "FigureS1_clinical_ctDNA_supplement.pdf"), width = 8.0, height = 11.0)
print(fig)
dev.off()

png(file.path(out_dir, "FigureS1_clinical_ctDNA_supplement.png"), width = 8.0, height = 11.0, units = "in", res = 600)
print(fig)
dev.off()

cat("Saved:\n")
cat(file.path(out_dir, "FigureS1_clinical_ctDNA_supplement.pdf"), "\n")
cat(file.path(out_dir, "FigureS1_clinical_ctDNA_supplement.png"), "\n")
