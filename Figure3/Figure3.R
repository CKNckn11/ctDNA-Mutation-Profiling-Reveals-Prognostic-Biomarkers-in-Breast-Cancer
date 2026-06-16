##############################################
#######Figure3B
##############################################

head(maf_sub@clinical.data)

# Check
for (f in c("cohort", "Age", "pT", "pN", "Stage", "ER", "Ki67")) {
  cat("\n", f, "\n")
  print(table(maf_sub@clinical.data[[f]], useNA = "ifany"))
}

# Then use this to draw the oncoplot:

variant_colors <- c(
  Missense_Mutation = "#377EB8",
  Nonsense_Mutation = "#E64B35",
  Frame_Shift_Del = "#4DBBD5",
  Frame_Shift_Ins = "#00A087",
  In_Frame_Del = "#F39B7F",
  In_Frame_Ins = "#8491B4",
  Splice_Site = "#91D1C2",
  Translation_Start_Site = "#DC0000",
  Nonstop_Mutation = "#7E6148",
  Multi_Hit = "#8DD3C7"
)


ann_colors <- list(
  cohort = c(
    "Kunming" = "#4DB3C8",
    "Tsinghua" = "#B94747"
  ),
  Age = c(
    "<60" = "#E7A6B0",
    ">=60" = "#B94747"
  ),
  pT = c(
    "T0/is" = "#F2E2F2",
    "T1" = "#D8B7D8",
    "T2" = "#B987C0",
    "T3" = "#8B4B9C",
    "T4" = "#5A236E"
  ),
  pN = c(
    "N0" = "#D9D9D9",
    "N1" = "#B6A6D8",
    "N2" = "#8C6BB1",
    "N3" = "#5E3C99"
  ),
  Stage = c(
    "0" = "#F0F0F0",
    "I" = "#DDECC9",
    "II" = "#A8DDB5",
    "III" = "#43A2CA",
    "IV" = "#0868AC"
  ),
  ER = c(
    "Negative" = "#9ECAE1",
    "Positive" = "#DE2D26"
  ),
  Ki67 = c(
    "Negative" = "#9ECAE1",
    "Positive" = "#DE2D26"
  ),
  yunzhong_subtype =  c(
    "Luminal" = "#E2B34C",          # Yellow
    "Luminal A" = "#6DB7E3",        # Blue
    "Luminal B" = "#D96B63",        # Red
    "HER2-enriched" = "#9E9E9E",    # Gray
    "Triple Negative" = "#4F79B7"   # Dark blue
  )
)
clinicalFeatures <- c("cohort", "Age", "pT", "pN", "Stage", "ER", "Ki67",'yunzhong_subtype')

clinicalFeatures_use <- clinicalFeatures[
  sapply(clinicalFeatures, function(x) {
    vals <- maf_sub@clinical.data[[x]]
    vals <- vals[!is.na(vals) & vals != ""]
    length(unique(vals)) > 0
  })
]

ann_colors_use <- ann_colors[clinicalFeatures_use]
genes_oncoplot <- c(
  "MUC12",
  "MUC5AC",
  "MUC5B",
  "MUC16",
  "LRP1",
  "KMT2D",
  "MUC17",
  "LAMA5",
  "HSPG2",
  "ADAMTSL3",
  "FCGBP",
  "LAMA1",
  "CREBBP"
)


pdf(
  file.path(outdir, "Figure3B.pdf"),
  width = 12,
  height = 9,
  useDingbats = FALSE
)

oncoplot(
  maf = maf_sub,
  genes = genes_oncoplot,
  #top = 30,
  clinicalFeatures = clinicalFeatures_use,
  colors = variant_colors,
  annotationColor = ann_colors_use,
  sortByAnnotation = TRUE,
  showTumorSampleBarcodes = FALSE,
  drawRowBar = TRUE,
  drawColBar = TRUE,
  removeNonMutated = T,
  titleText = "Clinical annotated cfDNA mutation landscape"
)

dev.off()







##############################################
#######Figure3C
##############################################
############################################################
## 1. Sample-level estimated TMB
############################################################

exome_size_mb <- 38

sample_burden_df <- mut_data %>%
  dplyr::filter(
    !is.na(.data[[sample_col]]),
    !is.na(.data[[gene_col]]),
    .data[[gene_col]] != ""
  ) %>%
  dplyr::group_by(sample_id = .data[[sample_col]]) %>%
  dplyr::summarise(
    mutation_count = dplyr::n(),
    mutated_genes = dplyr::n_distinct(.data[[gene_col]]),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    estimated_TMB = mutation_count / exome_size_mb,
    mutation_burden = mutation_count
  )

head(sample_burden_df)

############################################################
## 2. Merge mutation, subtype, and estimated TMB data
############################################################

mut_subtype_df <- mut_data %>%
  dplyr::select(
    sample_id = all_of(sample_col),
    Gene = all_of(gene_col)
  ) %>%
  dplyr::mutate(
    sample_id = as.character(sample_id),
    Gene = toupper(as.character(Gene))
  ) %>%
  dplyr::filter(
    !is.na(sample_id),
    !is.na(Gene),
    Gene != ""
  ) %>%
  dplyr::left_join(clinical_subtype_df, by = "sample_id") %>%
  dplyr::left_join(sample_burden_df, by = "sample_id") %>%
  dplyr::filter(!is.na(subtype))


############################################################
## 3. Top 5 among 47 candidate genes within each subtype
############################################################

candidate_genes <- toupper(c(
  "ARID1B", "MUC17", "MUC5AC", "LRP1", "C3", "LAMA1",
  "MUC16", "LAMA5", "ADAMTSL3", "KMT2D", "SMARCA4",
  "HSPG2", "SETD2", "CCL2", "JAK1", "CXCR4", "FN1",
  "MUC5B", "TGFB1", "CREBBP", "COL5A1", "SMARCB1",
  "CFB", "MMP14", "COL1A1", "ITGAV", "CTLA4", "KRT18",
  "FCGBP", "EP300", "NFKBIA", "DNMT3A", "ARID1A", "FAP",
  "CCR2", "NFKB1", "LAG3", "STAT3", "MMP9", "IKBKB",
  "COL6A1", "TGFBR1", "MUC4", "PDCD1", "MMP2", "NLRP3",
  "TET2"
))

subtype_n_df <- clinical_subtype_df %>%
  dplyr::count(subtype, name = "total_samples")

top5_gene_subtype_df <- mut_subtype_df %>%
  dplyr::filter(Gene %in% candidate_genes) %>%
  dplyr::distinct(
    sample_id,
    subtype,
    Gene,
    estimated_TMB,
    mutation_burden
  ) %>%
  dplyr::group_by(subtype, Gene) %>%
  dplyr::summarise(
    mutated_samples = dplyr::n_distinct(sample_id),
    median_estimated_TMB = median(estimated_TMB, na.rm = TRUE),
    mean_estimated_TMB = mean(estimated_TMB, na.rm = TRUE),
    median_mutation_burden = median(mutation_burden, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::left_join(subtype_n_df, by = "subtype") %>%
  dplyr::mutate(
    mutation_frequency = mutated_samples / total_samples * 100
  ) %>%
  dplyr::filter(mutated_samples >= 1) %>%
  dplyr::group_by(subtype) %>%
  dplyr::arrange(
    desc(mutation_frequency),
    desc(mutated_samples),
    desc(median_estimated_TMB),
    .by_group = TRUE
  ) %>%
  dplyr::slice_head(n = 5) %>%
  dplyr::ungroup()

top5_gene_subtype_df


p_top5_candidate_dot <- ggplot(
  top5_gene_subtype_df,
  aes(
    x = subtype,
    y = Gene
  )
) +
  geom_point(
    aes(
      size = mutation_frequency,
      color = median_estimated_TMB
    ),
    alpha = 0.9
  ) +
  scale_size_continuous(
    range = c(3, 10),
    name = "Mutation frequency (%)"
  ) +
  scale_color_gradient(
    low = "lightblue",
    high = "#D73027",
    name = "Median estimated TMB"
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = "Top 5 mutated candidate genes across breast cancer subtypes",
    subtitle = "Dot color represents the median estimated TMB of samples carrying each gene mutation"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey30"),
    axis.text.x = element_text(angle = 35, hjust = 1, color = "black", face = "bold"),
    axis.text.y = element_text(color = "black", face = "italic"),
    legend.position = "right"
  )

p_top5_candidate_dot






##############################################
#######Figure3D
##############################################

library(dplyr)
library(ggplot2)
library(stringr)

# Convert enrichResult to data.frame.
# If ego_bp is not loaded in the current session, read the exported GO table.
if (exists("ego_bp")) {
  go_df <- as.data.frame(ego_bp)
} else {
  go_file <- file.path("data", "expanded_TME_inflammation_curated_GO_enrichment.tsv")
  go_df <- read.delim(go_file, header = TRUE, sep = "\t", check.names = FALSE)
}

# Convert GeneRatio like "8/44" to numeric value
ratio_to_numeric <- function(x) {
  if (is.numeric(x)) {
    return(x)
  }
  sapply(strsplit(as.character(x), "/"), function(z) {
    if (length(z) == 2) {
      as.numeric(z[1]) / as.numeric(z[2])
    } else {
      as.numeric(z[1])
    }
  })
}

plot_df <- go_df %>%
  mutate(
    GeneRatio_num = ratio_to_numeric(GeneRatio),
    Count = if ("Count" %in% colnames(.)) Count else Gene_Count,
    neg_log10_FDR = -log10(p.adjust)
  ) %>%
  arrange(p.adjust) %>%
  slice_head(n = 14) %>%   # number of terms shown
  arrange(GeneRatio_num)

# Draw GO BP dotplot
panel_d <- ggplot(
  plot_df,
  aes(
    x = GeneRatio_num,
    y = factor(Description, levels = Description)
  )
) +
  geom_point(
    aes(size = Count, color = neg_log10_FDR),
    alpha = 0.95
  ) +
  scale_color_gradient(
    low = "#2C5A99",
    high = "#9D2F2F",
    name = expression(-log[10](FDR))
  ) +
  scale_size_continuous(
    name = "Gene count",
    range = c(3, 9)
  ) +
  labs(
    x = "Gene ratio",
    y = NULL
  ) +
  theme_classic(base_size = 10) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    plot.margin = margin(5, 10, 5, 5)
  )

panel_d

ggsave(
  file.path("output", "Figure3D_GO_BP_dotplot.pdf"),
  panel_d,
  width = 6.5,
  height = 5,
  useDingbats = FALSE
)

ggsave(
  file.path("output", "Figure3D_GO_BP_dotplot.png"),
  panel_d,
  width = 6.5,
  height = 5,
  dpi = 300
)

# Optional short file name for quick testing in the current working directory
ggsave(
  "GO_BP_dotplot.pdf",
  panel_d,
  width = 6.5,
  height = 5,
  useDingbats = FALSE
)

