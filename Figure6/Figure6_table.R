## Generate a publication-style clinical summary table for selected BC GEO cohorts.
## Input:  E:/code/400T/GEO/KM_risk_score_results/*_clinical_extracted.csv
## Output: Desktop CSV/HTML, and DOCX if flextable + officer are installed.

base_dir <- "E:/code/400T/GEO/KM_risk_score_results"
out_dir <- "C:/Users/cknckn11/Desktop"

datasets <- list(
  GSE4922  = "GSE4922_GPL96_clinical_extracted.csv",
  GSE6532  = "GSE6532_GPL570_clinical_extracted.csv",
  GSE7390  = "GSE7390_clinical_extracted.csv",
  GSE11121 = "GSE11121_clinical_extracted.csv",
  GSE19783 = "GSE19783_GPL6480_clinical_extracted.csv",
  GSE31448 = "GSE31448_clinical_extracted.csv",
  GSE45255 = "GSE45255_clinical_extracted.csv"
)

read_one <- function(file) {
  path <- file.path(base_dir, file)
  if (!file.exists(path)) stop("Missing file: ", path)
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

dat <- lapply(datasets, read_one)
n_by_ds <- vapply(dat, nrow, integer(1))

clean_value <- function(x) {
  x <- as.character(x)
  x[x %in% c("", "NA", "N/A", "na", "--", "null", "NULL")] <- NA
  x
}

get_col <- function(df, col) {
  if (col %in% names(df)) clean_value(df[[col]]) else rep(NA_character_, nrow(df))
}

num_value <- function(x) {
  x <- clean_value(x)
  out <- rep(NA_real_, length(x))
  hit <- !is.na(x) & grepl("-?\\d+(\\.\\d+)?", x)
  out[hit] <- as.numeric(sub(".*?(-?\\d+(\\.\\d+)?).*", "\\1", x[hit]))
  out
}

format_npct <- function(n, denom) {
  if (is.na(n) || n == 0) return("")
  sprintf("%d (%.1f%%)", n, 100 * n / denom)
}

count_category <- function(x, category, denom) {
  format_npct(sum(x == category, na.rm = TRUE), denom)
}

age_group <- function(df, ds) {
  x <- switch(ds,
    GSE4922 = get_col(df, "age_at_diagnosis"),
    GSE6532 = get_col(df, "age"),
    GSE7390 = get_col(df, "age"),
    GSE31448 = get_col(df, "age_at_diagnosis"),
    GSE45255 = get_col(df, "patient_age"),
    rep(NA_character_, nrow(df))
  )
  y <- rep(NA_character_, length(x))
  y[grepl("<=\\s*50", x)] <- "Age <50"
  y[grepl(">\\s*50", x)] <- "Age >=50"
  z <- num_value(x)
  idx <- is.na(y) & !is.na(z)
  y[idx & z < 50] <- "Age <50"
  y[idx & z >= 50] <- "Age >=50"
  y
}

survival_group <- function(df, ds) {
  y <- rep(NA_character_, nrow(df))
  if (ds == "GSE7390") {
    e <- num_value(get_col(df, "e_os"))
    y[e == 0] <- "Living"; y[e == 1] <- "Dead"
  } else if (ds == "GSE19783") {
    s <- tolower(get_col(df, "death_status"))
    y[grepl("alive", s)] <- "Living"
    y[grepl("dead", s)] <- "Dead"
  }
  y
}

metastasis_group <- function(df, ds) {
  e <- switch(ds,
    GSE4922 = num_value(get_col(df, "dfs_event_0_censored_1_event_defined_as_any_type_of_recurrence_local_regional_or_distant_or_death_from_breast_cancer")),
    GSE6532 = num_value(get_col(df, "e_dmfs")),
    GSE7390 = num_value(get_col(df, "e_dmfs")),
    GSE11121 = num_value(get_col(df, "e_dmfs")),
    GSE31448 = num_value(get_col(df, "dfs_evt")),
    GSE45255 = num_value(get_col(df, "dmfs_event_defined_as_distant_metastasis_or_death_from_breast_cancer")),
    rep(NA_real_, nrow(df))
  )
  y <- rep(NA_character_, length(e))
  y[e == 0] <- "M0 / no distant event"
  y[e == 1] <- "M1 / distant event"
  y
}

grade_group <- function(df, ds) {
  x <- switch(ds,
    GSE4922 = get_col(df, "elston_ngs_histologic_grade"),
    GSE6532 = get_col(df, "grade"),
    GSE7390 = get_col(df, "grade"),
    GSE11121 = get_col(df, "grade"),
    GSE31448 = get_col(df, "sbr_grade"),
    GSE45255 = get_col(df, "histological_grade"),
    rep(NA_character_, nrow(df))
  )
  y <- rep(NA_character_, length(x))
  y[grepl("1", x)] <- "Grade 1"
  y[grepl("2", x)] <- "Grade 2"
  y[grepl("3", x)] <- "Grade 3"
  y
}

er_group <- function(df, ds) {
  x <- switch(ds,
    GSE4922 = get_col(df, "er_status"),
    GSE6532 = get_col(df, "er"),
    GSE7390 = get_col(df, "er"),
    GSE19783 = get_col(df, "estrogen_receptor_status"),
    GSE31448 = get_col(df, "er_ihc"),
    GSE45255 = get_col(df, "er_status"),
    rep(NA_character_, nrow(df))
  )
  s <- tolower(x)
  y <- rep(NA_character_, length(x))
  y[grepl("negative|er-|^0$", s)] <- "ER negative"
  y[grepl("positive|er\\+|^1$", s)] <- "ER positive"
  y
}

pr_group <- function(df, ds) {
  x <- switch(ds,
    GSE6532 = get_col(df, "pgr"),
    GSE31448 = get_col(df, "pr_ihc"),
    GSE45255 = get_col(df, "pgr_status"),
    rep(NA_character_, nrow(df))
  )
  s <- tolower(x)
  y <- rep(NA_character_, length(x))
  y[grepl("negative|pgr-|^0$", s)] <- "PR negative"
  y[grepl("positive|pgr\\+|^1$", s)] <- "PR positive"
  y
}

her2_group <- function(df, ds) {
  x <- switch(ds,
    GSE19783 = get_col(df, "her2_fish_status"),
    GSE31448 = get_col(df, "erbb2_ihc_status"),
    GSE45255 = get_col(df, "her2_status"),
    rep(NA_character_, nrow(df))
  )
  s <- tolower(x)
  y <- rep(NA_character_, length(x))
  y[grepl("negative|he-|neg|^0$", s)] <- "HER2 negative"
  y[grepl("positive|he\\+|pos|^1$", s)] <- "HER2 positive"
  y
}

tumor_stage_group <- function(df, ds) {
  y <- rep(NA_character_, nrow(df))
  cm <- rep(NA_real_, nrow(df))
  if (ds == "GSE4922") cm <- num_value(get_col(df, "tumor_size_mm")) / 10
  if (ds == "GSE6532") cm <- num_value(get_col(df, "size"))
  if (ds == "GSE7390") cm <- num_value(get_col(df, "size"))
  if (ds == "GSE11121") cm <- num_value(get_col(df, "size_in_cm"))
  if (ds == "GSE45255") cm <- num_value(get_col(df, "size_mm")) / 10
  if (ds == "GSE31448") {
    s <- tolower(get_col(df, "pt"))
    y[grepl("pt1", s)] <- "T1"; y[grepl("pt2", s)] <- "T2"
    y[grepl("pt3", s)] <- "T3"; y[grepl("pt4", s)] <- "T4"
  }
  idx <- is.na(y) & !is.na(cm)
  y[idx & cm <= 2] <- "T1"
  y[idx & cm > 2 & cm <= 5] <- "T2"
  y[idx & cm > 5] <- "T3"
  y
}

node_group <- function(df, ds) {
  x <- switch(ds,
    GSE4922 = get_col(df, "lymph_node_status"),
    GSE6532 = get_col(df, "node"),
    GSE7390 = get_col(df, "node"),
    GSE11121 = get_col(df, "node"),
    GSE31448 = get_col(df, "pn"),
    GSE45255 = get_col(df, "ln_status"),
    rep(NA_character_, nrow(df))
  )
  s <- tolower(x)
  y <- rep(NA_character_, length(x))
  y[grepl("ln-|negative|^0$", s)] <- "N0 / node negative"
  y[grepl("ln\\+|positive|pos|^1$", s)] <- "N1 / node positive"
  y[grepl("^2$", s)] <- "N2"
  y[grepl("^3$", s)] <- "N3"
  y
}

subtype_group <- function(df, ds) {
  x <- switch(ds,
    GSE19783 = get_col(df, "breast_cancer_subtype"),
    GSE31448 = get_col(df, "molecular_subtype"),
    rep(NA_character_, nrow(df))
  )
  s <- tolower(x)
  y <- rep(NA_character_, length(x))
  y[grepl("basal|tnbc", s)] <- "Basal-like / TNBC"
  y[grepl("erbb2|her2", s)] <- "HER2-enriched"
  y[grepl("lum a|luminala|luminal a", s)] <- "Luminal A"
  y[grepl("lum b|luminalb|luminal b", s)] <- "Luminal B"
  y
}

row_def <- data.frame(
  Characteristic = c(
    rep("Age (years)", 2), rep("Survival status", 2),
    rep("Metastasis status", 2), rep("Grade", 3),
    rep("ER status", 2), rep("PR status", 2), rep("HER2 status", 2),
    rep("Tumor stage", 4), rep("Lymph node stage", 4), rep("Subtype", 4)
  ),
  Category = c(
    "Age <50", "Age >=50",
    "Living", "Dead",
    "M0 / no distant event", "M1 / distant event",
    "Grade 1", "Grade 2", "Grade 3",
    "ER negative", "ER positive",
    "PR negative", "PR positive",
    "HER2 negative", "HER2 positive",
    "T1", "T2", "T3", "T4",
    "N0 / node negative", "N1 / node positive", "N2", "N3",
    "Basal-like / TNBC", "HER2-enriched", "Luminal A", "Luminal B"
  ),
  Type = c(
    rep("Age", 2), rep("Survival", 2), rep("Met", 2), rep("Grade", 3),
    rep("ER", 2), rep("PR", 2), rep("HER2", 2), rep("T", 4),
    rep("N", 4), rep("Subtype", 4)
  ),
  stringsAsFactors = FALSE
)

group_by_type <- list(
  Age = age_group, Survival = survival_group, Met = metastasis_group,
  Grade = grade_group, ER = er_group, PR = pr_group, HER2 = her2_group,
  T = tumor_stage_group, N = node_group, Subtype = subtype_group
)

tab <- row_def[, c("Characteristic", "Category")]
for (ds in names(dat)) {
  values <- character(nrow(row_def))
  for (i in seq_len(nrow(row_def))) {
    g <- group_by_type[[row_def$Type[i]]](dat[[ds]], ds)
    values[i] <- count_category(g, row_def$Category[i], n_by_ds[[ds]])
  }
  tab[[paste0(ds, "\n(n=", n_by_ds[[ds]], ")")]] <- values
}

csv_out <- file.path(out_dir, "selected_GEO_clinical_summary_R.csv")
html_out <- file.path(out_dir, "selected_GEO_clinical_summary_R.html")
docx_out <- file.path(out_dir, "selected_GEO_clinical_summary_R.docx")
pdf_out <- file.path(out_dir, "selected_GEO_clinical_summary_R.pdf")

write.csv(tab, csv_out, row.names = FALSE, fileEncoding = "UTF-8")

html_escape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x
}

headers <- paste0("<th>", gsub("\n", "<br>", html_escape(names(tab))), "</th>", collapse = "")
body <- character(nrow(tab))
for (i in seq_len(nrow(tab))) {
  cls <- if (i == 1 || tab$Characteristic[i] != tab$Characteristic[i - 1]) " class='group-start'" else ""
  cells <- paste0("<td>", html_escape(as.character(unlist(tab[i, ], use.names = FALSE))), "</td>", collapse = "")
  body[i] <- paste0("<tr", cls, ">", cells, "</tr>")
}

css <- "
<style>
body{font-family:'Times New Roman',serif;margin:24px}
table{border-collapse:collapse;font-size:12px}
caption{caption-side:top;text-align:left;font-weight:bold;font-size:16px;margin-bottom:10px}
th,td{padding:6px 9px;border-bottom:1px solid #ddd;text-align:center;white-space:nowrap}
th{border-top:3px solid #000;border-bottom:2px solid #000;background:#f7f7f7}
td:first-child,td:nth-child(2),th:first-child,th:nth-child(2){text-align:left}
tr.group-start td{border-top:2px solid #000}
</style>"

html <- paste0(
  "<!doctype html><html><head><meta charset='utf-8'>", css, "</head><body>",
  "<table><caption>Table 1<br><span style='font-weight:normal'>Summary of selected BC-related GEO expression datasets and corresponding clinical characteristics</span></caption>",
  "<thead><tr>", headers, "</tr></thead><tbody>",
  paste(body, collapse = "\n"),
  "</tbody></table></body></html>"
)
writeLines(html, html_out, useBytes = TRUE)

if (requireNamespace("flextable", quietly = TRUE) &&
    requireNamespace("officer", quietly = TRUE)) {
  ft <- flextable::flextable(tab)
  ft <- flextable::set_caption(
    ft,
    caption = "Table 1. Summary of selected BC-related GEO expression datasets and corresponding clinical characteristics"
  )
  ft <- flextable::theme_booktabs(ft)
  ft <- flextable::fontsize(ft, size = 8, part = "all")
  ft <- flextable::font(ft, fontname = "Times New Roman", part = "all")
  ft <- flextable::align(ft, align = "center", part = "all")
  ft <- flextable::align(ft, j = c("Characteristic", "Category"), align = "left", part = "all")
  ft <- flextable::autofit(ft)
  doc <- officer::read_docx()
  doc <- officer::body_add_flextable(doc, ft)
  print(doc, target = docx_out)
} else {
  message("Packages flextable/officer are not installed; DOCX was not generated.")
}

if (requireNamespace("gridExtra", quietly = TRUE) &&
    requireNamespace("grid", quietly = TRUE)) {
  pdf_tab <- tab
  names(pdf_tab) <- gsub("\n", "\n", names(pdf_tab), fixed = TRUE)
  tg <- gridExtra::tableGrob(
    pdf_tab,
    rows = NULL,
    theme = gridExtra::ttheme_minimal(
      base_size = 6.5,
      base_family = "Times",
      core = list(
        fg_params = list(hjust = 0.5, x = 0.5),
        padding = grid::unit(c(3, 3), "pt")
      ),
      colhead = list(
        fg_params = list(fontface = "bold", hjust = 0.5, x = 0.5),
        padding = grid::unit(c(4, 4), "pt")
      )
    )
  )
  pdf(pdf_out, width = 16, height = 9, family = "Times")
  grid::grid.newpage()
  grid::grid.text(
    "Table 1. Summary of selected BC-related GEO expression datasets and corresponding clinical characteristics",
    x = 0.01, y = 0.98, just = c("left", "top"),
    gp = grid::gpar(fontfamily = "Times", fontsize = 11, fontface = "bold")
  )
  grid::pushViewport(grid::viewport(y = 0.46, height = 0.88))
  grid::grid.draw(tg)
  grid::popViewport()
  dev.off()
} else {
  message("Package gridExtra is not installed; PDF was not generated.")
}

message("Saved: ", csv_out)
message("Saved: ", html_out)
if (file.exists(docx_out)) message("Saved: ", docx_out)
if (file.exists(pdf_out)) message("Saved: ", pdf_out)
