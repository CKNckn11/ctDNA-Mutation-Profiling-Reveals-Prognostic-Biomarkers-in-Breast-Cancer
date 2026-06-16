
mut <- read.table(
  "E:/code/400T/TCGA/TCGA-BRCA/TCGA-BRCA.somaticmutation_wxs.tsv",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

head(mut)
colnames(mut)
dim(mut)


library(dplyr)

maf_df <- mut %>%
  transmute(
    Hugo_Symbol = gene,
    Chromosome = gsub("chr","",chrom),
    Start_Position = start,
    End_Position = end,
    Reference_Allele = ref,
    Tumor_Seq_Allele2 = alt,
    Tumor_Sample_Barcode = substr(sample,1,16),  # 用sample列
    Variant_Classification = effect
  )

maf_df$Variant_Classification <- dplyr::recode(
  maf_df$Variant_Classification,
  "missense_variant" = "Missense_Mutation",
  "synonymous_variant" = "Silent",
  "stop_gained" = "Nonsense_Mutation",
  "frameshift_variant" = "Frame_Shift_Ins",
  "splice_region_variant" = "Splice_Site",
  .default = "Missense_Mutation"
)

# 添加必需列
maf_df$Variant_Type <- "SNP"

maf <- read.maf(maf = maf_df)
genes_oncoplot=candidate_genes
library(maftools)

oncoplot(
  maf = maf,
  #top = 30,
  colors = mut_col,
  genes =genes_oncoplot ,
  removeNonMutated = TRUE,
  #clinicalFeatures = clinicalFeatures_use,
  annotationColor = ann_colors_use,
  sortByAnnotation = TRUE,
  showTumorSampleBarcodes = FALSE,
  drawRowBar = TRUE,
  drawColBar = TRUE,
  titleText = "Clinical annotated cfDNA mutation landscape"
)





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


maf <- read.maf(
  maf = maf_df,
  clinicalData = clinical_anno_use,
  verbose = TRUE
)

## 
colnames(maf@clinical.data)
head(maf@clinical.data)

## 
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

##
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

## 
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

##
ann_colors_use <- ann_colors[names(ann_colors) %in% clinicalFeatures_use]

##
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

##############################################
#######FigureS3A
##############################################

pdf(
  file = "E:/code/400T/TCGA/TCGA-BRCA/TCGA_BRCA_clinical_annotated_oncoplot.pdf",
  width = 12,
  height = 8
)

oncoplot(
  maf = maf,
  #genes = genes_oncoplot,
  top = 40,
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
#######FigureSB
##############################################
library(maftools)
library(dplyr)
library(ggplot2)
library(patchwork)
## BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
library(BSgenome.Hsapiens.UCSC.hg38)

#outdir <- "Figure"
#dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

maf_df <- maf@data
## -----------------------------
## -----------------------------
snv_df <- maf_df %>%
  dplyr::filter(
    nchar(Reference_Allele) == 1,
    nchar(Tumor_Seq_Allele2) == 1,
    Reference_Allele %in% c("A", "C", "G", "T"),
    Tumor_Seq_Allele2 %in% c("A", "C", "G", "T"),
    Reference_Allele != Tumor_Seq_Allele2
  ) %>%
  dplyr::mutate(
    Chromosome = as.character(Chromosome),
    Start_Position = as.numeric(Start_Position)
  )

## -----------------------------
## -----------------------------
snv_df <- snv_df %>%
  dplyr::mutate(
    chr_use = ifelse(grepl("^chr", Chromosome), Chromosome, paste0("chr", Chromosome))
  )

## -----------------------------
## -----------------------------
snv_df$tri_context <- BSgenome::getSeq(
  BSgenome.Hsapiens.UCSC.hg38,
  names = snv_df$chr_use,
  start = snv_df$Start_Position - 1,
  end = snv_df$Start_Position + 1
) %>%
  as.character() %>%
  toupper()

## -----------------------------
## -----------------------------
comp_base <- function(x) {
  chartr("ACGT", "TGCA", x)
}

revcomp_tri <- function(x) {
  sapply(x, function(z) {
    paste0(rev(strsplit(comp_base(z), "")[[1]]), collapse = "")
  })
}

snv_df <- snv_df %>%
  dplyr::mutate(
    ref = Reference_Allele,
    alt = Tumor_Seq_Allele2,
    
    tri_use = ifelse(ref %in% c("C", "T"), tri_context, revcomp_tri(tri_context)),
    ref_use = ifelse(ref %in% c("C", "T"), ref, comp_base(ref)),
    alt_use = ifelse(ref %in% c("C", "T"), alt, comp_base(alt)),
    
    mut_type = paste0(ref_use, ">", alt_use),
    left_base = substr(tri_use, 1, 1),
    right_base = substr(tri_use, 3, 3),
    context = paste0(left_base, ref_use, right_base)
  ) %>%
  dplyr::filter(
    mut_type %in% c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"),
    left_base %in% c("A", "C", "G", "T"),
    right_base %in% c("A", "C", "G", "T")
  )

## -----------------------------
## -----------------------------
all_context <- expand.grid(
  left_base = c("A", "C", "G", "T"),
  right_base = c("A", "C", "G", "T"),
  mut_type = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"),
  stringsAsFactors = FALSE
) %>%
  dplyr::mutate(
    context = paste0(left_base, substr(mut_type, 1, 1), right_base)
  )

spectrum_df <- snv_df %>%
  dplyr::count(mut_type, context, name = "count") %>%
  dplyr::right_join(all_context, by = c("mut_type", "context")) %>%
  dplyr::mutate(
    count = ifelse(is.na(count), 0, count),
    relative_contribution = count / sum(count)
  )

spectrum_df$mut_type <- factor(
  spectrum_df$mut_type,
  levels = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
)

spectrum_df$context <- factor(
  spectrum_df$context,
  levels = unique(all_context$context)
)

## -----------------------------
## -----------------------------
pB <- ggplot(
  spectrum_df,
  aes(x = context, y = relative_contribution, fill = mut_type)
) +
  geom_col(width = 0.75, color = "black", linewidth = 0.2) +
  facet_grid(. ~ mut_type, scales = "free_x", space = "free_x") +
  scale_fill_manual(
    values = c(
      "C>A" = "#2CA9E1",
      "C>G" = "black",
      "C>T" = "#E64B35",
      "T>A" = "grey70",
      "T>C" = "#8BC34A",
      "T>G" = "#F4A3A3"
    )
  ) +
  labs(
    x = NULL,
    y = "Relative contribution"
  ) +
  theme_bw(base_size = 10) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "grey80", color = "black"),
    strip.text = element_text(color = "black", size = 9),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
    axis.text.y = element_text(size = 8),
    panel.grid = element_blank()
  )




##############################################
#######FigureSC
##############################################

variant_df <- maf_df %>%
  dplyr::count(Variant_Classification, name = "count") %>%
  dplyr::arrange(count)

pC <- ggplot(
  variant_df,
  aes(x = count, y = reorder(Variant_Classification, count), fill = Variant_Classification)
) +
  geom_col(width = 0.75) +
  scale_fill_manual(
    values = c(
      "Missense_Mutation" = "#49A148",
      "Nonsense_Mutation" = "#E41A1C",
      "Frame_Shift_Ins" = "#6A3D9A",
      "Frame_Shift_Del" = "#377EB8",
      "In_Frame_Del" = "#FF7F00",
      "In_Frame_Ins" = "#A65628"
    )
  ) +
  labs(
    title = "Variant type",
    x = NULL,
    y = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5),
    axis.text.y = element_text(size = 9)
  )


snv_type_df <- snv_df %>%
  dplyr::mutate(
    SNV_type = paste0(Reference_Allele, ">", Tumor_Seq_Allele2),
    SNV_type = case_when(
      SNV_type == "G>T" ~ "C>A",
      SNV_type == "G>C" ~ "C>G",
      SNV_type == "G>A" ~ "C>T",
      SNV_type == "A>T" ~ "T>A",
      SNV_type == "A>G" ~ "T>C",
      SNV_type == "A>C" ~ "T>G",
      TRUE ~ SNV_type
    )
  ) %>%
  dplyr::filter(SNV_type %in% c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")) %>%
  dplyr::count(SNV_type, name = "count")

snv_type_df$SNV_type <- factor(
  snv_type_df$SNV_type,
  levels = c("T>G", "T>A", "T>C", "C>T", "C>G", "C>A")
)


##############################################
#######FigureSD
##############################################


pD <- ggplot(
  snv_type_df,
  aes(x = count, y = SNV_type, fill = SNV_type)
) +
  geom_col(width = 0.75) +
  scale_fill_manual(
    values = c(
      "C>A" = "#3C91CF",
      "C>G" = "#405BA7",
      "C>T" = "#EF4B43",
      "T>C" = "#F5B82E",
      "T>A" = "#55A868",
      "T>G" = "#F28E2B"
    )
  ) +
  labs(
    title = "SNV type",
    x = NULL,
    y = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5),
    axis.text.y = element_text(size = 9)
  )
