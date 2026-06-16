suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tidyr)
  library(patchwork)
})

##############################################
####### Figure4 setup
##############################################

###############################################################
##clin
clinical_df <- read.table(
  "E:/code/400T/TCGA/TCGA-BRCA/TCGA-BRCA.clinical.tsv.gz",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

dim(clinical_df)
colnames(clinical_df)
head(clinical_df)


clinical_surv_df <- clinical_df %>%
  dplyr::transmute(
    patient_id = substr(submitter_id, 1, 12),
    OS_event = ifelse(vital_status.demographic == "Dead", 1, 0),
    OS_time = ifelse(
      !is.na(days_to_death.demographic),
      days_to_death.demographic,
      days_to_last_follow_up.diagnoses
    ),
    vital_status = vital_status.demographic
  ) %>%
  dplyr::filter(
    !is.na(OS_time),
    OS_time > 0,
    !is.na(OS_event)
  )


clinical_anno <- clinical_df %>%
  dplyr::transmute(
    Tumor_Sample_Barcode = substr(sample, 1, 16),
    Gender = gender.demographic,
    Age = as.numeric(age_at_index.demographic),
    Stage = ajcc_pathologic_stage.diagnoses,
    T_stage = ajcc_pathologic_t.diagnoses,
    N_stage = ajcc_pathologic_n.diagnoses,
    M_stage = ajcc_pathologic_m.diagnoses,
    Grade = tumor_grade.diagnoses,
    Sample_Type = sample_type.samples
  ) %>%
  dplyr::distinct(Tumor_Sample_Barcode, .keep_all = TRUE)






library(dplyr)
library(maftools)

## 1. Remove Grade because all Grade values are NA here

clinical_anno_use <- clinical_anno %>%
  transmute(
    Tumor_Sample_Barcode,
    Age = ifelse(Age < 60, "<60", ">=60"),
    Gender = Gender,
    
    pT = case_when(
      grepl("Tis|T0", T_stage, ignore.case = TRUE) ~ "T0/is",
      grepl("T1", T_stage, ignore.case = TRUE) ~ "T1",
      grepl("T2", T_stage, ignore.case = TRUE) ~ "T2",
      grepl("T3", T_stage, ignore.case = TRUE) ~ "T3",
      grepl("T4", T_stage, ignore.case = TRUE) ~ "T4",
      TRUE ~ NA_character_
    ),
    
    pN = case_when(
      grepl("N0", N_stage, ignore.case = TRUE) ~ "N0",
      grepl("N1", N_stage, ignore.case = TRUE) ~ "N1",
      grepl("N2", N_stage, ignore.case = TRUE) ~ "N2",
      grepl("N3", N_stage, ignore.case = TRUE) ~ "N3",
      TRUE ~ NA_character_
    ),
    
    pM = case_when(
      grepl("M0", M_stage, ignore.case = TRUE) ~ "M0",
      grepl("M1", M_stage, ignore.case = TRUE) ~ "M1",
      TRUE ~ NA_character_
    ),
    
    Stage = case_when(
      grepl("Stage I[A-C]*$", Stage, ignore.case = TRUE) ~ "I",
      grepl("Stage II[A-C]*$", Stage, ignore.case = TRUE) ~ "II",
      grepl("Stage III[A-C]*$", Stage, ignore.case = TRUE) ~ "III",
      grepl("Stage IV", Stage, ignore.case = TRUE) ~ "IV",
      TRUE ~ NA_character_
    ),
    
    Sample_Type = Sample_Type
  )

## If the result above is > 0, continue with the code below
## If it is 0, barcode lengths are inconsistent and should be changed to 12-digit patient IDs

## 3. Regenerate the MAF object with clinical annotations
maf <- read.maf(
  maf = maf_df,
  clinicalData = clinical_anno_use,
  verbose = TRUE
)

## 4. Check whether clinical data were successfully added to the MAF object
colnames(maf@clinical.data)
head(maf@clinical.data)

## 5. Set mutation colors
mut_col <- c(
  Missense_Mutation = "#377EB8",
  Nonsense_Mutation = "#E64B35",
  Frame_Shift_Del = "#4DBBD5",
  Frame_Shift_Ins = "#00A087",
  In_Frame_Del = "#F39B7F",
  In_Frame_Ins = "#8491B4",
  Splice_Site = "#91D1C2",
  Translation_Start_Site = "#DC0000",
  Nonstop_Mutation = "#7E6148"
)

## 6. Set clinical annotation colors
ann_colors <- list(
  #Age = c(
  #  "<60" = "#E7A6B0",
  #  ">=60" = "#B94747"
  #),
  
  Gender = c(
    "female" = "#E64B35",
    "male" = "#4DBBD5",
    "not reported" = "#BDBDBD"
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
  
  pM = c(
    "M0" = "#BDBDBD",
    "M1" = "#E64B35"
  ),
  
  Stage = c(
    "I" = "#DDECC9",
    "II" = "#A8DDB5",
    "III" = "#43A2CA",
    "IV" = "#0868AC"
  ),
  
  Sample_Type = c(
    "Primary Tumor" = "#4DBBD5",
    "Solid Tissue Normal" = "#7A7A7A",
    "Metastatic" = "#E64B35"
  )
)

## 7. Select clinical annotations that actually exist and have values
clinicalFeatures <- c(
  #"Age",
  "Gender",
  "pT",
  "pN",
  "pM",
  "Stage",
  "Sample_Type"
)

clinicalFeatures_use <- clinicalFeatures[
  sapply(clinicalFeatures, function(x) {
    vals <- maf@clinical.data[[x]]
    vals <- vals[!is.na(vals) & vals != ""]
    length(unique(vals)) > 0
  })
]

clinicalFeatures_use

## 8. Keep only colors that are actually used
ann_colors_use <- ann_colors[names(ann_colors) %in% clinicalFeatures_use]

## 9. Draw oncoplot
oncoplot(
  maf = maf,
  genes = genes_oncoplot,
  colors = mut_col,
  removeNonMutated = TRUE,
  clinicalFeatures = clinicalFeatures_use,
  annotationColor = ann_colors_use,
  sortByAnnotation = TRUE,
  showTumorSampleBarcodes = FALSE,
  drawRowBar = TRUE,
  drawColBar = TRUE,
  titleText = "Clinical annotated TCGA-BRCA mutation landscape"
)



pdf(
  file = "E:/code/400T/TCGA/TCGA-BRCA/TCGA_BRCA_clinical_annotated_oncoplot.pdf",
  width = 12,
  height = 8
)

oncoplot(
  maf = maf,
  genes = genes_oncoplot,
  colors = mut_col,
  removeNonMutated = TRUE,
  clinicalFeatures = clinicalFeatures_use,
  annotationColor = ann_colors_use,
  sortByAnnotation = TRUE,
  showTumorSampleBarcodes = FALSE,
  drawRowBar = TRUE,
  drawColBar = TRUE,
  fontSize = 0.6,
  titleText = "Clinical annotated TCGA-BRCA mutation landscape"
)

dev.off()






##############################################
####### Shared dotplot function
##############################################

library(dplyr)
library(ggplot2)

############################################################
## 1. Sample-level estimated TMB
############################################################

exome_size_mb <- 38

sample_burden_df <- maf_df %>%
  dplyr::filter(
    !is.na(Tumor_Sample_Barcode),
    !is.na(Hugo_Symbol),
    Hugo_Symbol != ""
  ) %>%
  dplyr::group_by(sample_id = Tumor_Sample_Barcode) %>%
  dplyr::summarise(
    mutation_count = dplyr::n(),
    mutated_genes = dplyr::n_distinct(Hugo_Symbol),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    estimated_TMB = mutation_count / exome_size_mb,
    mutation_burden = mutation_count
  )

head(sample_burden_df)



stage_df <- clinical_anno %>%
  dplyr::select(
    sample_id = Tumor_Sample_Barcode,
    Stage
  ) %>%
  dplyr::filter(
    !is.na(Stage),
    Stage != ""
  ) %>%
  dplyr::distinct(sample_id, .keep_all = TRUE)

table(stage_df$Stage, useNA = "ifany")



stage_df <- clinical_df %>%
  dplyr::transmute(
    sample_id = substr(sample, 1, 16),
    Stage_raw = ajcc_pathologic_stage.diagnoses
  ) %>%
  dplyr::mutate(
    Stage = dplyr::case_when(
      grepl("Stage I[A-C]*$", Stage_raw, ignore.case = TRUE) ~ "Stage I",
      grepl("Stage II[A-C]*$", Stage_raw, ignore.case = TRUE) ~ "Stage II",
      grepl("Stage III[A-C]*$", Stage_raw, ignore.case = TRUE) ~ "Stage III",
      grepl("Stage IV", Stage_raw, ignore.case = TRUE) ~ "Stage IV",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::select(sample_id, Stage) %>%
  dplyr::filter(
    !is.na(Stage),
    Stage != ""
  ) %>%
  dplyr::distinct(sample_id, .keep_all = TRUE)

table(stage_df$Stage, useNA = "ifany")







mut_stage_df <- maf_df %>%
  dplyr::select(
    sample_id = Tumor_Sample_Barcode,
    Gene = Hugo_Symbol
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
  dplyr::left_join(stage_df, by = "sample_id") %>%
  dplyr::left_join(sample_burden_df, by = "sample_id") %>%
  dplyr::filter(!is.na(Stage))

head(mut_stage_df)



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


stage_n_df <- stage_df %>%
  dplyr::count(Stage, name = "total_samples")

top5_gene_stage_df <- mut_stage_df %>%
  dplyr::filter(Gene %in% candidate_genes) %>%
  dplyr::distinct(
    sample_id,
    Stage,
    Gene,
    estimated_TMB,
    mutation_burden
  ) %>%
  dplyr::group_by(Stage, Gene) %>%
  dplyr::summarise(
    mutated_samples = dplyr::n_distinct(sample_id),
    median_estimated_TMB = median(estimated_TMB, na.rm = TRUE),
    mean_estimated_TMB = mean(estimated_TMB, na.rm = TRUE),
    median_mutation_burden = median(mutation_burden, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::left_join(stage_n_df, by = "Stage") %>%
  dplyr::mutate(
    mutation_frequency = mutated_samples / total_samples * 100
  ) %>%
  dplyr::filter(mutated_samples >= 1) %>%
  dplyr::group_by(Stage) %>%
  dplyr::arrange(
    desc(mutation_frequency),
    desc(mutated_samples),
    desc(median_estimated_TMB),
    .by_group = TRUE
  ) %>%
  dplyr::slice_head(n = 5) %>%
  dplyr::ungroup()

top5_gene_stage_df


stage_order <- c(
  "Stage I",
  "Stage II",
  "Stage III",
  "Stage IV"
)

stage_order <- stage_order[
  stage_order %in% unique(as.character(top5_gene_stage_df$Stage))
]

stage_order <- c(
  stage_order,
  setdiff(unique(as.character(top5_gene_stage_df$Stage)), stage_order)
)

gene_order_by_stage <- top5_gene_stage_df %>%
  dplyr::mutate(
    Stage = as.character(Stage),
    Gene = as.character(Gene)
  ) %>%
  dplyr::mutate(
    Stage = factor(Stage, levels = stage_order)
  ) %>%
  dplyr::arrange(
    Stage,
    desc(mutation_frequency),
    desc(mutated_samples),
    desc(median_estimated_TMB)
  ) %>%
  dplyr::pull(Gene) %>%
  unique()

top5_gene_stage_df_plot <- top5_gene_stage_df %>%
  dplyr::mutate(
    Stage = factor(as.character(Stage), levels = stage_order),
    Gene = factor(as.character(Gene), levels = rev(gene_order_by_stage))
  )



p_top5_stage_tmb_dot <- ggplot(
  top5_gene_stage_df_plot,
  aes(
    x = Stage,
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
    title = "Top 5 mutated candidate genes across TCGA-BRCA stages",
    subtitle = "Dot color represents the median estimated TMB of samples carrying each gene mutation"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 15
    ),
    plot.subtitle = element_text(
      hjust = 0.5,
      size = 10,
      color = "grey30"
    ),
    axis.text.x = element_text(
      angle = 35,
      hjust = 1,
      color = "black",
      face = "bold"
    ),
    axis.text.y = element_text(
      color = "black",
      face = "italic"
    ),
    legend.position = "right"
  )

p_top5_stage_tmb_dot



ggsave(
  filename = "E:/code/400T/TCGA/TCGA-BRCA/TCGA_BRCA_stage_TMB_dotplot.pdf",
  plot = p_top5_stage_tmb_dot,
  width = 6,
  height = 7,
  device = cairo_pdf
)





library(dplyr)
library(ggplot2)

############################################################
## 1. Sample-level estimated TMB
############################################################

exome_size_mb <- 38

sample_burden_df <- maf_df %>%
  dplyr::filter(
    !is.na(Tumor_Sample_Barcode),
    !is.na(Hugo_Symbol),
    Hugo_Symbol != ""
  ) %>%
  dplyr::group_by(sample_id = Tumor_Sample_Barcode) %>%
  dplyr::summarise(
    mutation_count = dplyr::n(),
    mutated_genes = dplyr::n_distinct(Hugo_Symbol),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    estimated_TMB = mutation_count / exome_size_mb,
    mutation_burden = mutation_count
  )

############################################################
## 2. Prepare pT groups
############################################################

pt_df <- clinical_anno %>%
  dplyr::select(
    sample_id = Tumor_Sample_Barcode,
    pT
  ) %>%
  dplyr::filter(
    !is.na(pT),
    pT != ""
  ) %>%
  dplyr::distinct(sample_id, .keep_all = TRUE)

table(pt_df$pT, useNA = "ifany")

############################################################
## 3. Merge mutation, pT, and estimated TMB data
############################################################

mut_pt_df <- maf_df %>%
  dplyr::select(
    sample_id = Tumor_Sample_Barcode,
    Gene = Hugo_Symbol
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
  dplyr::left_join(pt_df, by = "sample_id") %>%
  dplyr::left_join(sample_burden_df, by = "sample_id") %>%
  dplyr::filter(!is.na(pT))

############################################################
## 4. Top 5 candidate genes within each pT group
############################################################

pt_n_df <- pt_df %>%
  dplyr::count(pT, name = "total_samples")

top5_gene_pt_df <- mut_pt_df %>%
  dplyr::filter(Gene %in% candidate_genes) %>%
  dplyr::distinct(
    sample_id,
    pT,
    Gene,
    estimated_TMB,
    mutation_burden
  ) %>%
  dplyr::group_by(pT, Gene) %>%
  dplyr::summarise(
    mutated_samples = dplyr::n_distinct(sample_id),
    median_estimated_TMB = median(estimated_TMB, na.rm = TRUE),
    mean_estimated_TMB = mean(estimated_TMB, na.rm = TRUE),
    median_mutation_burden = median(mutation_burden, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::left_join(pt_n_df, by = "pT") %>%
  dplyr::mutate(
    mutation_frequency = mutated_samples / total_samples * 100
  ) %>%
  dplyr::filter(mutated_samples >= 1) %>%
  dplyr::group_by(pT) %>%
  dplyr::arrange(
    desc(mutation_frequency),
    desc(mutated_samples),
    desc(median_estimated_TMB),
    .by_group = TRUE
  ) %>%
  dplyr::slice_head(n = 5) %>%
  dplyr::ungroup()

top5_gene_pt_df

############################################################
## 5. Order pT and Gene
############################################################

pt_order <- c("T0/is", "T1", "T2", "T3", "T4")

pt_order <- pt_order[
  pt_order %in% unique(as.character(top5_gene_pt_df$pT))
]

pt_order <- c(
  pt_order,
  setdiff(unique(as.character(top5_gene_pt_df$pT)), pt_order)
)

gene_order_by_pt <- top5_gene_pt_df %>%
  dplyr::mutate(
    pT = as.character(pT),
    Gene = as.character(Gene)
  ) %>%
  dplyr::mutate(
    pT = factor(pT, levels = pt_order)
  ) %>%
  dplyr::arrange(
    pT,
    desc(mutation_frequency),
    desc(mutated_samples),
    desc(median_estimated_TMB)
  ) %>%
  dplyr::pull(Gene) %>%
  unique()

top5_gene_pt_df_plot <- top5_gene_pt_df %>%
  dplyr::mutate(
    pT = factor(as.character(pT), levels = pt_order),
    Gene = factor(as.character(Gene), levels = rev(gene_order_by_pt))
  )

############################################################
## 6. Draw pT TMB dotplot
############################################################

p_top5_pt_tmb_dot <- ggplot(
  top5_gene_pt_df_plot,
  aes(
    x = pT,
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
    title = "Top 5 mutated candidate genes across TCGA-BRCA pT stages",
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

p_top5_pt_tmb_dot
ggsave(
  filename = "E:/code/400T/TCGA/TCGA-BRCA/TCGA_BRCA_PT_stage_TMB_dotplot.pdf",
  plot = p_top5_pt_tmb_dot,
  width = 6,
  height = 7,
  device = cairo_pdf
)



############################################################################

############################################################
## 1. Prepare pN groups
############################################################

pn_df <- clinical_anno %>%
  dplyr::select(
    sample_id = Tumor_Sample_Barcode,
    pN
  ) %>%
  dplyr::filter(
    !is.na(pN),
    pN != ""
  ) %>%
  dplyr::distinct(sample_id, .keep_all = TRUE)

table(pn_df$pN, useNA = "ifany")

############################################################
## 2. Merge mutation, pN, and estimated TMB data
############################################################

mut_pn_df <- maf_df %>%
  dplyr::select(
    sample_id = Tumor_Sample_Barcode,
    Gene = Hugo_Symbol
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
  dplyr::left_join(pn_df, by = "sample_id") %>%
  dplyr::left_join(sample_burden_df, by = "sample_id") %>%
  dplyr::filter(!is.na(pN))

############################################################
## 3. Top 5 candidate genes within each pN group
############################################################

pn_n_df <- pn_df %>%
  dplyr::count(pN, name = "total_samples")

top5_gene_pn_df <- mut_pn_df %>%
  dplyr::filter(Gene %in% candidate_genes) %>%
  dplyr::distinct(
    sample_id,
    pN,
    Gene,
    estimated_TMB,
    mutation_burden
  ) %>%
  dplyr::group_by(pN, Gene) %>%
  dplyr::summarise(
    mutated_samples = dplyr::n_distinct(sample_id),
    median_estimated_TMB = median(estimated_TMB, na.rm = TRUE),
    mean_estimated_TMB = mean(estimated_TMB, na.rm = TRUE),
    median_mutation_burden = median(mutation_burden, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::left_join(pn_n_df, by = "pN") %>%
  dplyr::mutate(
    mutation_frequency = mutated_samples / total_samples * 100
  ) %>%
  dplyr::filter(mutated_samples >= 1) %>%
  dplyr::group_by(pN) %>%
  dplyr::arrange(
    desc(mutation_frequency),
    desc(mutated_samples),
    desc(median_estimated_TMB),
    .by_group = TRUE
  ) %>%
  dplyr::slice_head(n = 5) %>%
  dplyr::ungroup()

top5_gene_pn_df

############################################################
## 4. Order pN and Gene
############################################################

pn_order <- c("N0", "N1", "N2", "N3")

pn_order <- pn_order[
  pn_order %in% unique(as.character(top5_gene_pn_df$pN))
]

pn_order <- c(
  pn_order,
  setdiff(unique(as.character(top5_gene_pn_df$pN)), pn_order)
)

gene_order_by_pn <- top5_gene_pn_df %>%
  dplyr::mutate(
    pN = as.character(pN),
    Gene = as.character(Gene)
  ) %>%
  dplyr::mutate(
    pN = factor(pN, levels = pn_order)
  ) %>%
  dplyr::arrange(
    pN,
    desc(mutation_frequency),
    desc(mutated_samples),
    desc(median_estimated_TMB)
  ) %>%
  dplyr::pull(Gene) %>%
  unique()

top5_gene_pn_df_plot <- top5_gene_pn_df %>%
  dplyr::mutate(
    pN = factor(as.character(pN), levels = pn_order),
    Gene = factor(as.character(Gene), levels = rev(gene_order_by_pn))
  )

############################################################
## 5. Draw pN TMB dotplot
############################################################

p_top5_pn_tmb_dot <- ggplot(
  top5_gene_pn_df_plot,
  aes(
    x = pN,
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
    title = "Top 5 mutated candidate genes across TCGA-BRCA pN stages",
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

p_top5_pn_tmb_dot

ggsave(
  filename = "E:/code/400T/TCGA/TCGA-BRCA/TCGA_BRCA_pN_TMB_dotplot.pdf",
  plot = p_top5_pn_tmb_dot,
  width = 6,
  height = 7,
  device = cairo_pdf
)


###################################################

############################################################
## TMB dotplot grouped by Sample_Type / type
############################################################

library(dplyr)
library(ggplot2)

############################################################
## 1. Sample-level estimated TMB
############################################################

exome_size_mb <- 38

sample_burden_df <- maf_df %>%
  dplyr::filter(
    !is.na(Tumor_Sample_Barcode),
    !is.na(Hugo_Symbol),
    Hugo_Symbol != ""
  ) %>%
  dplyr::group_by(sample_id = Tumor_Sample_Barcode) %>%
  dplyr::summarise(
    mutation_count = dplyr::n(),
    mutated_genes = dplyr::n_distinct(Hugo_Symbol),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    estimated_TMB = mutation_count / exome_size_mb,
    mutation_burden = mutation_count
  )

############################################################
## 2. Prepare Sample_Type groups
############################################################

type_df <- clinical_anno %>%
  dplyr::select(
    sample_id = Tumor_Sample_Barcode,
    Sample_Type
  ) %>%
  dplyr::filter(
    !is.na(Sample_Type),
    Sample_Type != ""
  ) %>%
  dplyr::distinct(sample_id, .keep_all = TRUE)

table(type_df$Sample_Type, useNA = "ifany")

############################################################
## 3. Merge mutation, Sample_Type, and estimated TMB data
############################################################

mut_type_df <- maf_df %>%
  dplyr::select(
    sample_id = Tumor_Sample_Barcode,
    Gene = Hugo_Symbol
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
  dplyr::left_join(type_df, by = "sample_id") %>%
  dplyr::left_join(sample_burden_df, by = "sample_id") %>%
  dplyr::filter(!is.na(Sample_Type))

############################################################
## 4. Top 5 candidate genes within each Sample_Type
############################################################

type_n_df <- type_df %>%
  dplyr::count(Sample_Type, name = "total_samples")

top5_gene_type_df <- mut_type_df %>%
  dplyr::filter(Gene %in% candidate_genes) %>%
  dplyr::distinct(
    sample_id,
    Sample_Type,
    Gene,
    estimated_TMB,
    mutation_burden
  ) %>%
  dplyr::group_by(Sample_Type, Gene) %>%
  dplyr::summarise(
    mutated_samples = dplyr::n_distinct(sample_id),
    median_estimated_TMB = median(estimated_TMB, na.rm = TRUE),
    mean_estimated_TMB = mean(estimated_TMB, na.rm = TRUE),
    median_mutation_burden = median(mutation_burden, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::left_join(type_n_df, by = "Sample_Type") %>%
  dplyr::mutate(
    mutation_frequency = mutated_samples / total_samples * 100
  ) %>%
  dplyr::filter(mutated_samples >= 1) %>%
  dplyr::group_by(Sample_Type) %>%
  dplyr::arrange(
    desc(mutation_frequency),
    desc(mutated_samples),
    desc(median_estimated_TMB),
    .by_group = TRUE
  ) %>%
  dplyr::slice_head(n = 5) %>%
  dplyr::ungroup()

top5_gene_type_df

############################################################
## 5. Order Sample_Type and Gene
############################################################

type_order <- c(
  "Primary Tumor",
  "Metastatic",
  "Solid Tissue Normal"
)

type_order <- type_order[
  type_order %in% unique(as.character(top5_gene_type_df$Sample_Type))
]

type_order <- c(
  type_order,
  setdiff(unique(as.character(top5_gene_type_df$Sample_Type)), type_order)
)

gene_order_by_type <- top5_gene_type_df %>%
  dplyr::mutate(
    Sample_Type = as.character(Sample_Type),
    Gene = as.character(Gene)
  ) %>%
  dplyr::mutate(
    Sample_Type = factor(Sample_Type, levels = type_order)
  ) %>%
  dplyr::arrange(
    Sample_Type,
    desc(mutation_frequency),
    desc(mutated_samples),
    desc(median_estimated_TMB)
  ) %>%
  dplyr::pull(Gene) %>%
  unique()

top5_gene_type_df_plot <- top5_gene_type_df %>%
  dplyr::mutate(
    Sample_Type = factor(as.character(Sample_Type), levels = type_order),
    Gene = factor(as.character(Gene), levels = rev(gene_order_by_type))
  )

############################################################
## 6. Draw Sample_Type TMB dotplot
############################################################

p_top5_type_tmb_dot <- ggplot(
  top5_gene_type_df_plot,
  aes(
    x = Sample_Type,
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
    title = "Sample type-specific mutation burden of candidate genes",
    subtitle = "Dot color represents the median estimated TMB of samples carrying each gene mutation"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 15
    ),
    plot.subtitle = element_text(
      hjust = 0.5,
      size = 10,
      color = "grey30"
    ),
    axis.text.x = element_text(
      angle = 35,
      hjust = 1,
      color = "black",
      face = "bold"
    ),
    axis.text.y = element_text(
      color = "black",
      face = "italic"
    ),
    legend.position = "right"
  )

p_top5_type_tmb_dot

ggsave(
  filename = "E:/code/400T/TCGA/TCGA-BRCA/TCGA_BRCA_SampleType_TMB_dotplot.pdf",
  plot = p_top5_type_tmb_dot,
  width = 6,
  height = 7,
  device = cairo_pdf
)