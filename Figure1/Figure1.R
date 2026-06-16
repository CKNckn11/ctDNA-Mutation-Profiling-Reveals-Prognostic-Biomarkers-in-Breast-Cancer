############################################################
## Clinical characteristics pie plots
## Stage / Subtype / ER / PR / HER2 / Ki67 / Age / cfDNA burden
## Labels in legend only
############################################################

library(dplyr)
library(ggplot2)
library(patchwork)

############################################################
## 1. Resolve paths and import data
############################################################

# Resolve the script directory so the script can be run from any working directory.
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  normalizePath(getwd())
}

script_dir <- get_script_dir()
data_dir <- file.path(script_dir, "data")
outdir <- file.path(script_dir, "output")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# Load the exported Figure 1 data object first. This file contains the
# clinical_df_all object used below and is the preferred local data source.
data_rdata <- file.path(data_dir, "Figure1_data.RData")
if (file.exists(data_rdata)) {
  load(data_rdata)
} else {
  # Fallback for sharing or rerunning the script when only CSV files are present.
  clinical_csv <- file.path(data_dir, "clinical_df_all.csv")
  if (!file.exists(clinical_csv)) {
    stop("Missing Figure 1 data. Expected either: ", data_rdata, " or ", clinical_csv)
  }
  clinical_df_all <- read.csv(clinical_csv, stringsAsFactors = FALSE, check.names = FALSE)
}

if (!exists("clinical_df_all")) {
  stop("The imported Figure 1 data must contain clinical_df_all.")
}

plot_df <- clinical_df_all

# Validate the columns required for the clinical pie plots before plotting.
required_cols <- c(
  "Stage_group",
  "ER_group",
  "PR_group",
  "HER2_group",
  "Ki67_status_group",
  "Age_group",
  "cfDNA_group"
)
missing_cols <- setdiff(required_cols, colnames(plot_df))
if (length(missing_cols) > 0) {
  stop("Missing required columns in clinical_df_all: ", paste(missing_cols, collapse = ", "))
}

############################################################
## 2. Prepare variables
############################################################

plot_df <- plot_df %>%
  mutate(
    Stage_group = as.character(Stage_group),
    ER_group = as.character(ER_group),
    PR_group = as.character(PR_group),
    HER2_group = as.character(HER2_group),
    Ki67_status_group = as.character(Ki67_status_group),
    Age_group = as.character(Age_group),
    cfDNA_group = as.character(cfDNA_group)
  )

## Prefer yunzhong_subtype for subtype
if ("yunzhong_subtype" %in% colnames(plot_df)) {
  plot_df <- plot_df %>%
    mutate(Subtype_plot = as.character(yunzhong_subtype))
} else if ("Breast_subtype" %in% colnames(plot_df)) {
  plot_df <- plot_df %>%
    mutate(Subtype_plot = as.character(Breast_subtype))
} else {
  stop("Could not find yunzhong_subtype or Breast_subtype.")
}

## Clean invalid values and remove stage 0
plot_df <- plot_df %>%
  mutate(
    Stage_group = ifelse(Stage_group %in% c("", "/", "NA", "N/A", "0"), NA, Stage_group),
    Subtype_plot = ifelse(Subtype_plot %in% c("", "/", "NA", "N/A"), NA, Subtype_plot),
    ER_group = ifelse(ER_group %in% c("", "/", "NA", "N/A"), NA, ER_group),
    PR_group = ifelse(PR_group %in% c("", "/", "NA", "N/A"), NA, PR_group),
    HER2_group = ifelse(HER2_group %in% c("", "/", "NA", "N/A"), NA, HER2_group),
    Ki67_status_group = ifelse(Ki67_status_group %in% c("", "/", "NA", "N/A"), NA, Ki67_status_group),
    Age_group = ifelse(Age_group %in% c("", "/", "NA", "N/A"), NA, Age_group),
    cfDNA_group = ifelse(cfDNA_group %in% c("", "/", "NA", "N/A"), NA, cfDNA_group)
  )

############################################################
## 3. Colors
############################################################

stage_cols <- c(
  "I" = "#FDD49E",
  "II" = "#FDBF6F",
  "III" = "#EF3B2C",
  "IV" = "#99000D"
)

subtype_cols <- c(
  "Luminal" = "#E7B64F",
  "Luminal A" = "#4DBBD5",
  "Luminal B" = "#E64B35",
  "HER2-enriched" = "#9E9E9E",
  "HER2-equivocal" = "#F0B75E",
  "Triple Negative" = "#3C5488",
  "TNBC" = "#3C5488"
)

binary_cols <- c(
  "negative" = "#6BAED6",
  "positive" = "#DE2D26",
  "Negative" = "#6BAED6",
  "Positive" = "#DE2D26"
)

her2_cols <- c(
  "negative" = "#6BAED6",
  "positive" = "#DE2D26",
  "equivocal" = "#F0B75E",
  "Negative" = "#6BAED6",
  "Positive" = "#DE2D26",
  "Equivocal" = "#F0B75E"
)

ki67_cols <- c(
  "low" = "#6BAED6",
  "high" = "#DE2D26",
  "Low" = "#6BAED6",
  "High" = "#DE2D26"
)

age_cols <- c(
  "<60" = "#91D1C2",
  ">=60" = "#DC0000",
  "< 60" = "#91D1C2",
  "≥60" = "#DC0000"
)

cfdna_cols <- c(
  "cfDNA-low" = "#6BAED6",
  "cfDNA-high" = "#DE2D26"
)

############################################################
## 4. Pie chart function: put counts in the legend
############################################################

plot_pie_legend_percent <- function(data, var, title, fill_cols = NULL) {
  
  pie_df <- data %>%
    dplyr::filter(
      !is.na(.data[[var]]),
      .data[[var]] != "",
      .data[[var]] != "/"
    ) %>%
    dplyr::count(.data[[var]], name = "n") %>%
    dplyr::mutate(
      group = as.character(.data[[var]]),
      percent = n / sum(n) * 100,
      legend_label = paste0(group, " (n=", n, ", ", sprintf("%.1f", percent), "%)")
    )
  
  pie_df$group <- factor(pie_df$group, levels = pie_df$group)
  
  p <- ggplot(
    pie_df,
    aes(x = "", y = n, fill = group)
  ) +
    geom_col(
      width = 1,
      color = "white",
      linewidth = 0.55
    ) +
    coord_polar(theta = "y") +
    labs(
      title = title,
      fill = NULL
    ) +
    theme_void(base_size = 11) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 13
      ),
      legend.position = "right",
      legend.text = element_text(
        size = 7.8,
        color = "black"
      ),
      legend.key.size = unit(0.38, "cm"),
      plot.margin = margin(4, 4, 4, 4)
    )
  
  if (!is.null(fill_cols)) {
    use_cols <- fill_cols[names(fill_cols) %in% as.character(pie_df$group)]
    p <- p +
      scale_fill_manual(
        values = use_cols,
        labels = setNames(pie_df$legend_label, pie_df$group)
      )
  } else {
    p <- p +
      scale_fill_discrete(
        labels = setNames(pie_df$legend_label, pie_df$group)
      )
  }
  
  return(p)
}

############################################################
## 5. Draw pie charts
############################################################

p_stage <- plot_pie_legend_percent(
  plot_df,
  "Stage_group",
  "Stage",
  stage_cols
)

p_subtype <- plot_pie_legend_percent(
  plot_df,
  "Subtype_plot",
  "Subtype",
  subtype_cols
)

p_er <- plot_pie_legend_percent(
  plot_df,
  "ER_group",
  "ER",
  binary_cols
)

p_pr <- plot_pie_legend_percent(
  plot_df,
  "PR_group",
  "PR",
  binary_cols
)

p_her2 <- plot_pie_legend_percent(
  plot_df,
  "HER2_group",
  "HER2",
  her2_cols
)

p_ki67 <- plot_pie_legend_percent(
  plot_df,
  "Ki67_status_group",
  "Ki67",
  ki67_cols
)

p_age <- plot_pie_legend_percent(
  plot_df,
  "Age_group",
  "Age",
  age_cols
)

p_cfdna <- plot_pie_legend_percent(
  plot_df,
  "cfDNA_group",
  "cfDNA burden",
  cfdna_cols
)

############################################################
## 6. Combine the 8-panel figure
############################################################

p_all_pie8 <- 
  (p_stage + p_subtype + p_er + p_pr) /
  (p_her2 + p_ki67 + p_age + p_cfdna) +
  patchwork::plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(
      face = "bold",
      size = 15
    )
  )

p_all_pie8


