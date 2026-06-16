
##############################################
#######Figure2A
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
  file.path(outdir, "Figure2A.pdf"),
  width = 12,
  height = 9,
  useDingbats = FALSE
)

oncoplot(
  maf = maf_sub,
  #genes = genes_oncoplot,
  top = 30,
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
#######Figure2B
##############################################
############################################################
## Dotplot:
## Top 5 mutated genes within each subtype
## Remove Luminal and keep Luminal A / Luminal B
## Dot size  = mutation frequency (%)
## Dot color = median estimated TMB of mutated samples
############################################################

library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)

############################################################
## 1. Extract mutation and clinical data from maf_sub
############################################################

mut_data <- maf_sub@data %>%
  as.data.frame()

clinical_data <- maf_sub@clinical.data %>%
  as.data.frame()

sample_col <- "Tumor_Sample_Barcode"
gene_col <- "Hugo_Symbol"
subtype_col <- "yunzhong_subtype"

############################################################
## 2. Sample-level estimated TMB
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
## 3. Prepare subtype information
## Remove standalone Luminal and keep Luminal A / Luminal B
############################################################

clinical_subtype_df <- clinical_data %>%
  dplyr::select(
    sample_id = all_of(sample_col),
    subtype = all_of(subtype_col)
  ) %>%
  dplyr::mutate(
    sample_id = as.character(sample_id),
    subtype = as.character(subtype),
    subtype = ifelse(subtype %in% c("", "/", "NA", "N/A"), NA, subtype)
  ) %>%
  dplyr::filter(!is.na(subtype)) %>%
  dplyr::filter(subtype != "Luminal") %>%
  dplyr::distinct(sample_id, subtype)

table(clinical_subtype_df$subtype)

############################################################
## 4. Merge mutation, subtype, and estimated TMB data
## Note: candidate genes are not restricted here; this is the genome-wide top 5
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
## 5. Total sample count for each subtype
############################################################

subtype_n_df <- clinical_subtype_df %>%
  dplyr::count(subtype, name = "total_samples")

print(subtype_n_df)

############################################################
## 6. Top 5 genes within each subtype across all genes
############################################################

top5_gene_subtype_df <- mut_subtype_df %>%
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

print(top5_gene_subtype_df)

#write.csv(
#  top5_gene_subtype_df,
#  file.path(outdir, "Subtype_top5_all_genes_no_Luminal_dotplot_data.csv"),
#  row.names = FALSE
#)

############################################################
## 7. Order subtype and Gene
## Subtype order:
## Luminal A -> Luminal B -> Triple Negative -> HER2-enriched
############################################################

subtype_order <- c(
  "Luminal A",
  "Luminal B",
  "Triple Negative",
  "HER2-enriched"
)

subtype_order <- subtype_order[
  subtype_order %in% unique(as.character(top5_gene_subtype_df$subtype))
]

subtype_order <- c(
  subtype_order,
  setdiff(unique(as.character(top5_gene_subtype_df$subtype)), subtype_order)
)

gene_order_by_subtype <- top5_gene_subtype_df %>%
  dplyr::mutate(
    subtype = as.character(subtype),
    Gene = as.character(Gene)
  ) %>%
  dplyr::mutate(
    subtype = factor(subtype, levels = subtype_order)
  ) %>%
  dplyr::arrange(
    subtype,
    desc(mutation_frequency),
    desc(mutated_samples),
    desc(median_estimated_TMB)
  ) %>%
  dplyr::pull(Gene) %>%
  unique()

top5_gene_subtype_df_plot <- top5_gene_subtype_df %>%
  dplyr::mutate(
    subtype = factor(as.character(subtype), levels = subtype_order),
    Gene = factor(as.character(Gene), levels = rev(gene_order_by_subtype))
  )

############################################################
## 8. Draw dotplot
############################################################

p_top5_allgenes_dot <- ggplot(
  top5_gene_subtype_df_plot,
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
    title = "Top 5 mutated genes across breast cancer subtypes",
    subtitle = "Dot size indicates mutation frequency; dot color indicates median estimated TMB"
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
      face = "bold",
      size = 11
    ),
    axis.text.y = element_text(
      color = "black",
      face = "italic",
      size = 10
    ),
    axis.line = element_line(linewidth = 0.5),
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    plot.margin = margin(10, 15, 10, 10)
  )

p_top5_allgenes_dot






##############################################
#######Figure2C
##############################################

library(maftools)

pdf(
  file.path(outdir, "Figure3C_maf_summary.pdf"),
  width = 8,
  height = 5,
  useDingbats = FALSE
)

plotmafSummary(
  maf = maf_sub,
  rmOutlier = TRUE,
  addStat = "median",
  dashboard = TRUE,
  titvRaw = FALSE
)

dev.off()



##############################################
#######Figure2D
##############################################

############################################################
## MutationalPatterns:
## Use cfDNA_mut to draw the 96 trinucleotide mutation spectrum
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(IRanges)
  library(BSgenome)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(MutationalPatterns)
  library(ggplot2)
})

############################################################
## 1. Set output directory
############################################################

outdir_sig <- file.path(outdir, "cfDNA_mut_MutationalPatterns")
dir.create(outdir_sig, recursive = TRUE, showWarnings = FALSE)

############################################################
## 2. Keep only SNV / SNP
############################################################

cfDNA_snv <- cfDNA_mut %>%
  dplyr::filter(
    Variant_Type == "SNP",
    !is.na(Chromosome),
    !is.na(Start_Position),
    !is.na(Reference_Allele),
    !is.na(Tumor_Seq_Allele2),
    Reference_Allele %in% c("A", "T", "C", "G"),
    Tumor_Seq_Allele2 %in% c("A", "T", "C", "G"),
    Reference_Allele != Tumor_Seq_Allele2
  ) %>%
  dplyr::mutate(
    Chromosome = as.character(Chromosome),
    Start_Position = as.numeric(Start_Position),
    End_Position = as.numeric(End_Position),
    sample_id_use = as.character(Tumor_Sample_Barcode)
  )

############################################################
## 3. Keep only primary chromosomes
############################################################

main_chr <- paste0("chr", c(1:22, "X", "Y"))

cfDNA_snv <- cfDNA_snv %>%
  dplyr::filter(Chromosome %in% main_chr)

cat("SNVs used:", nrow(cfDNA_snv), "\n")
cat("Samples used:", length(unique(cfDNA_snv$sample_id_use)), "\n")

############################################################
## 4. Convert to a GRanges list
############################################################

vcfs_list <- split(cfDNA_snv, cfDNA_snv$sample_id_use)

vcfs_gr <- lapply(vcfs_list, function(df) {
  
  gr <- GenomicRanges::GRanges(
    seqnames = df$Chromosome,
    ranges = IRanges::IRanges(
      start = df$Start_Position,
      end = df$Start_Position
    )
  )
  
  mcols(gr)$REF <- df$Reference_Allele
  mcols(gr)$ALT <- df$Tumor_Seq_Allele2
  
  ## Key step: set genome and seqlevel style
  GenomeInfoDb::seqlevelsStyle(gr) <- "UCSC"
  GenomeInfoDb::genome(gr) <- "hg38"
  
  return(gr)
})

############################################################
## 5. Build the 96 mutation matrix
############################################################

# Load required R packages
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)  # Assume the hg19 reference genome; change to hg38 or another genome if needed
library(MutationalPatterns)
library(dplyr)

#BiocManager::install("MutationalPatterns")

# Set the root directory, i.e., the path containing breast; adjust as needed
root_dir <- "/400T/cfDNA_resource/BAM/breast"  

# 1. Get all sample folder paths (first-level subfolders, e.g., RTCG0P0007-1-TWN1)
sample_folders <- list.dirs(
  path = root_dir, 
  full.names = TRUE, 
  recursive = FALSE
)

# 2. Iterate over sample folders and filter .filtered.vcf.gz files
vcf_files <- lapply(sample_folders, function(sample_folder) {
  # Find files ending with .filtered.vcf.gz in each sample folder
  filtered_vcf <- list.files(
    path = sample_folder, 
    pattern = "\\.filtered\\.vcf\\.gz$",  # Key step: match .filtered.vcf.gz
    full.names = TRUE
  )
  
  if (length(filtered_vcf) > 0) {
    return(filtered_vcf[1])  # If multiple files exist, take the first one (adjust as needed)
  } else {
    warning(paste("样本文件夹", basename(sample_folder), "内无 .filtered.vcf.gz 文件"))
    return(NULL)
  }
}) %>% unlist()  # Organize as a vector of file paths

# 3. Define sample names (extracted from folder names)
sample_names <- basename(sample_folders[-1])

# 4. Read VCF files as GRanges objects
# Specify the reference genome; hg19 is used here, change as needed, e.g., "BSgenome.Hsapiens.UCSC.hg38"
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"  
# Correct usage
vcfs <- read_vcfs_as_granges(
  vcf_files = vcf_files,          # Vector of VCF file paths (e.g., the 48 files you provided)
  sample_names = sample_names,    # Vector of sample names (one-to-one with vcf_files)
  genome = Hsapiens,              # BSgenome object (e.g., Hsapiens for hg38)
  group = "auto",                 # Chromosome grouping mode (default "auto+sex"; options include "auto")
  type = "snv",                   # Variant type (default "snv"; options include "indel")
  change_seqnames = TRUE,         # Whether to convert chromosome names (e.g., chr1 to 1)
  remove_duplicate_variants = TRUE # Whether to remove duplicate variants
)

# 5. Mutation type statistics and visualization (example workflow)
# Count mutation type occurrences
type_occurrences <- mut_type_occurrences(vcfs, ref_genome)  

# Draw the mutation spectrum plot
p1 <- plot_spectrum(type_occurrences)  
# Draw the mutation spectrum plot with confidence intervals
p2 <- plot_spectrum(type_occurrences, CT = TRUE)  
# Draw the mutation spectrum plot with hidden legend for panel assembly
p3 <- plot_spectrum(type_occurrences, CT = TRUE, legend = FALSE)  

# Load the plot assembly package
library(gridExtra)  
# Show assembled plots
grid.arrange(p1, p2, p3, ncol = 1, widths = c(3, 3, 1.75))  

# 6. Build the mutation matrix and visualize the 96 mutation spectrum
mut_mat <- mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)  


# 1. Calculate the row sums
row_sums <- rowSums(mut_mat)

# 2. Add row sums as a new column
mut_mat_with_sum <- cbind(mut_mat, sum = row_sums)

colnames(mut_mat_with_sum)[21]='breast_cancer'

# save(vcfs, file = '/400T/ckn/WES_code/breast_tissue/breast_cfDNA.Rdata')

# 1. Select the mutation matrix for a single sample (e.g., the first sample)
single_sample_mat <- mut_mat_with_sum[, 21, drop = FALSE]  # Preserve matrix dimensions

# 2. Draw the 96 mutation spectrum
plot_96_profile(single_sample_mat, condensed = TRUE)

# 3. Draw mutation spectra for multiple samples (e.g., the first 3 samples)
multi_sample_mat <- mut_mat[, 1:3]
plot_96_profile(multi_sample_mat, condensed = TRUE)
# Draw the 96 mutation signature spectrum (simplified display)
mut_mat_with_sum[,21]=as.numeric(mut_mat_with_sum[,21])
# Use drop = FALSE to preserve matrix dimensions, even when selecting one column
plot_96_profile(mut_mat_with_sum[, 21, drop = FALSE], condensed = TRUE)
# Draw sample-level mutation spectra with confidence intervals and legends
plot_spectrum(type_occurrences, by = "sample", CT = TRUE, legend = TRUE)  

# Optional: compare mutation spectra if grouping information such as tissue type is available
# Here all samples are assumed to belong to one group; replace tissue content as needed
tissue <- rep("cfDNA_breast", length(sample_names))  
plot_compare_profiles(type_occurrences, profile_names = tissue)  









