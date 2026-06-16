library(dplyr)
library(grid)
library(gtable)

############################################################
## 1. Prepare data
############################################################

table1_df <- clinical_df_all2 %>%
  distinct(clinical_id, .keep_all = TRUE)

n_total <- nrow(table1_df)

clean_value <- function(x) {
  x <- as.character(x)
  x[x %in% c("", "/", "NA", "N/A", "na", "Na", "NULL")] <- NA
  x
}

fmt_np <- function(n, denom = n_total) {
  ifelse(
    is.na(n) | n == 0,
    "",
    paste0(n, " (", sprintf("%.1f", n / denom * 100), "%)")
  )
}

make_rows <- function(data, var, label, levels_use = NULL, level_labels = NULL) {
  
  x <- clean_value(data[[var]])
  
  if (!is.null(levels_use)) {
    x <- factor(x, levels = levels_use)
  }
  
  tab <- table(x, useNA = "no")
  
  if (!is.null(levels_use)) {
    tab <- tab[levels_use]
    tab[is.na(tab)] <- 0
  }
  
  out <- data.frame(
    Characteristic = "",
    Level = names(tab),
    N_percent = fmt_np(as.numeric(tab)),
    stringsAsFactors = FALSE
  )
  
  if (!is.null(level_labels)) {
    out$Level <- level_labels[out$Level]
  }
  
  out$Characteristic[1] <- label
  out
}


############################################################
## 2. Build Table 1 content
############################################################

## Display pT as T1-T4
table1_df <- table1_df %>%
  mutate(
    pT_table = clean_value(pT_raw),
    pT_table = ifelse(!is.na(pT_table), paste0("T", pT_table), NA),
    
    Stage_table = clean_value(Stage_group),
    Stage_table = ifelse(Stage_table == "0", NA, Stage_table),
    
    ER_table = clean_value(ER_group),
    PR_table = clean_value(PR_group),
    HER2_table = clean_value(HER2_group),
    Ki67_table = clean_value(Ki67_status_group),
    M_table = clean_value(M_group),
    Subtype_table = clean_value(yunzhong_subtype),
    cfDNA_table = clean_value(cfDNA_group),
    TMB_table = clean_value(TMB_group)
  )

## Remove incomplete subtype categories such as Luminal if they should not be displayed
table1_df$Subtype_table[table1_df$Subtype_table == "Luminal"] <- NA

table_rows <- bind_rows(
  make_rows(
    table1_df,
    "Age_group",
    "Age (years)",
    levels_use = c("<60", ">=60")
  ),
  
  make_rows(
    table1_df,
    "pT_table",
    "Pathological tumor size (pT)",
    levels_use = c("T1", "T2", "T3", "T4")
  ),
  
  make_rows(
    table1_df,
    "cN_group",
    "Clinical nodal status",
    levels_use = c("negative", "positive"),
    level_labels = c(
      "negative" = "Negative",
      "positive" = "Positive"
    )
  ),
  
  make_rows(
    table1_df,
    "Stage_table",
    "Stage",
    levels_use = c("I", "II", "III", "IV")
  ),
  
  make_rows(
    table1_df,
    "ER_table",
    "ER status",
    levels_use = c("negative", "positive"),
    level_labels = c(
      "negative" = "Negative",
      "positive" = "Positive"
    )
  ),
  
  make_rows(
    table1_df,
    "PR_table",
    "PR status",
    levels_use = c("negative", "positive"),
    level_labels = c(
      "negative" = "Negative",
      "positive" = "Positive"
    )
  ),
  
  make_rows(
    table1_df,
    "HER2_table",
    "HER2 status",
    levels_use = c("negative", "equivocal", "positive"),
    level_labels = c(
      "negative" = "Negative",
      "equivocal" = "Equivocal",
      "positive" = "Positive"
    )
  ),
  
  make_rows(
    table1_df,
    "Ki67_table",
    "Ki67",
    levels_use = c("low", "high"),
    level_labels = c(
      "low" = "Low",
      "high" = "High"
    )
  ),
  
  make_rows(
    table1_df,
    "M_table",
    "Metastasis status",
    levels_use = c("M0", "M1")
  ),
  
  make_rows(
    table1_df,
    "Subtype_table",
    "Molecular subtype",
    levels_use = c(
      "Luminal A",
      "Luminal B",
      "HER2-enriched",
      "Triple Negative"
    )
  ),
  
  make_rows(
    table1_df,
    "cfDNA_table",
    "ctDNA group",
    levels_use = c("cfDNA-low", "cfDNA-high")
  ),
  
  make_rows(
    table1_df,
    "TMB_table",
    "TMB group",
    levels_use = c("TMB-low", "TMB-high")
  )
)

## Remove completely empty rows
table_rows <- table_rows %>%
  filter(N_percent != "")

#write.csv(
#  table_rows,
#  "E:/code/400T/TCGA/Figure/Figure1_Table1_clinicopathologic_characteristics.csv",
#  row.names = FALSE
#)


############################################################
## 3. Draw publication-style Table 1
############################################################

draw_table1 <- function(table_rows, n_total, out_pdf, out_png = NULL) {
  
  nr <- nrow(table_rows)
  
  ## Page parameters
  page_w <- 7.2
  page_h <- max(6.0, 1.15 + nr * 0.28)
  
  if (!is.null(out_png)) {
    png(out_png, width = page_w, height = page_h, units = "in", res = 600)
  } else {
    pdf(out_pdf, width = page_w, height = page_h)
  }
  
  grid.newpage()
  
  ## Title
  grid.text(
    paste0("Table 1 | Clinicopathologic patient characteristics (n = ", n_total, ")"),
    x = unit(0.02, "npc"),
    y = unit(0.97, "npc"),
    just = c("left", "top"),
    gp = gpar(fontsize = 16, fontface = "bold", col = "black")
  )
  
  ## Table area
  left <- 0.02
  right <- 0.98
  top <- 0.91
  row_h <- 0.032
  header_h <- 0.040
  
  col1_x <- 0.03
  col2_x <- 0.52
  col3_x <- 0.84
  
  table_top <- top
  table_bottom <- top - header_h - nr * row_h
  
  ## Header background
  grid.rect(
    x = unit((left + right) / 2, "npc"),
    y = unit(table_top - header_h / 2, "npc"),
    width = unit(right - left, "npc"),
    height = unit(header_h, "npc"),
    gp = gpar(fill = "#D9D8C8", col = NA)
  )
  
  ## Header text
  grid.text(
    "Clinicopathologic characteristics",
    x = unit(col1_x, "npc"),
    y = unit(table_top - header_h / 2, "npc"),
    just = c("left", "center"),
    gp = gpar(fontsize = 12, fontface = "bold")
  )
  
  grid.text(
    "N (%)",
    x = unit(col3_x, "npc"),
    y = unit(table_top - header_h / 2, "npc"),
    just = c("left", "center"),
    gp = gpar(fontsize = 12, fontface = "bold.italic")
  )
  
  ## Top line and header bottom line
  grid.lines(
    x = unit(c(left, right), "npc"),
    y = unit(c(table_top, table_top), "npc"),
    gp = gpar(lwd = 1.2, col = "black")
  )
  grid.lines(
    x = unit(c(left, right), "npc"),
    y = unit(c(table_top - header_h, table_top - header_h), "npc"),
    gp = gpar(lwd = 0.7, col = "black")
  )
  
  ## Row background and content
  y0 <- table_top - header_h
  
  for (i in seq_len(nr)) {
    
    y_mid <- y0 - (i - 0.5) * row_h
    y_line <- y0 - i * row_h
    
    ## Alternating light row background
    if (i %% 2 == 1) {
      grid.rect(
        x = unit((left + right) / 2, "npc"),
        y = unit(y_mid, "npc"),
        width = unit(right - left, "npc"),
        height = unit(row_h, "npc"),
        gp = gpar(fill = "#F4F4EF", col = NA)
      )
    }
    
    ## Add a bold line before each major section
    if (table_rows$Characteristic[i] != "") {
      grid.lines(
        x = unit(c(left, right), "npc"),
        y = unit(c(y0 - (i - 1) * row_h, y0 - (i - 1) * row_h), "npc"),
        gp = gpar(lwd = 0.9, col = "black")
      )
    }
    
    grid.text(
      table_rows$Characteristic[i],
      x = unit(col1_x, "npc"),
      y = unit(y_mid, "npc"),
      just = c("left", "center"),
      gp = gpar(fontsize = 11, col = "black")
    )
    
    grid.text(
      table_rows$Level[i],
      x = unit(col2_x, "npc"),
      y = unit(y_mid, "npc"),
      just = c("left", "center"),
      gp = gpar(fontsize = 11, col = "black")
    )
    
    grid.text(
      table_rows$N_percent[i],
      x = unit(col3_x, "npc"),
      y = unit(y_mid, "npc"),
      just = c("left", "center"),
      gp = gpar(fontsize = 11, col = "black")
    )
  }
  
  ## Bottom line
  grid.lines(
    x = unit(c(left, right), "npc"),
    y = unit(c(table_bottom, table_bottom), "npc"),
    gp = gpar(lwd = 1.2, col = "black")
  )
  
  dev.off()
  
  if (!is.null(out_png)) {
    pdf(out_pdf, width = page_w, height = page_h)
    draw_table1(table_rows, n_total, out_pdf = out_pdf, out_png = NULL)
  }
}


draw_table1(
  table_rows = table_rows,
  n_total = n_total,
  out_pdf = "E:/code/400T/TCGA/Figure/Figure1_Table1_clinicopathologic_characteristics.pdf",
  out_png = "E:/code/400T/TCGA/Figure/Figure1_Table1_clinicopathologic_characteristics.png"
)



