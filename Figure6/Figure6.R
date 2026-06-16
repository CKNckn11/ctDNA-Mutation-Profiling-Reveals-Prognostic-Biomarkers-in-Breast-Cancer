## KM curves for GEO cohorts using the 7-gene LASSO risk score
## No GEOquery, no hgu133plus2.db, no survminer.
## Input files are under E:/code/400T/GEO.

base_dir <- "E:/code/400T/GEO"
out_dir <- file.path(base_dir, "KM_risk_score_results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
old_pdf <- list.files(out_dir, pattern = "^KM_.*_risk_score\\.pdf$", full.names = TRUE)
if (length(old_pdf) > 0) unlink(old_pdf)
write_combined_sample_file <- FALSE

if (!requireNamespace("survival", quietly = TRUE)) {
  stop("The 'survival' package is required for KM/log-rank analysis.")
}
if (!requireNamespace("survminer", quietly = TRUE)) {
  stop("The 'survminer' package is required for ggsurvplot. Install it with: install.packages('survminer')")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("The 'ggplot2' package is required. Install it with: install.packages('ggplot2')")
}

## Coefficients from the LASSO result.
## Risk score = sum(coef_i * z-scored expression_i) within each GEO dataset.
risk_coef <- c(
  ARID1B =  1.16870213,
  JAK1   = -0.31956337,
  MUC4   =  0.31647725,
  CFB    = -0.24106792,
  PDCD1  = -0.22657677,
  NFKBIA = -0.22434813,
  FCGBP  = -0.07262415
)
risk_genes <- names(risk_coef)

gse_files <- c(
  # GSE2034 skipped: series matrix has relapse event but no usable follow-up time field.
  # GSE3494_GPL96 skipped: no survival endpoint fields in series matrix.
  # GSE3494_GPL97 skipped: no survival endpoint fields in series matrix.
  GSE1456_GPL96 = file.path(base_dir, "GSE1456", "GSE1456-GPL96_series_matrix.txt.gz"),
  GSE1456_GPL97 = file.path(base_dir, "GSE1456", "GSE1456-GPL97_series_matrix.txt.gz"),
  GSE4922_GPL96 = file.path(base_dir, "GSE4922", "GSE4922-GPL96_series_matrix.txt.gz"),
  GSE4922_GPL97 = file.path(base_dir, "GSE4922", "GSE4922-GPL97_series_matrix.txt.gz"),
  GSE7390 = file.path(base_dir, "GSE7390", "GSE7390_series_matrix.txt.gz"),
  GSE9195 = file.path(base_dir, "GSE9195", "GSE9195_series_matrix.txt.gz"),
  GSE16446 = file.path(base_dir, "GSE16446", "GSE16446_series_matrix.txt.gz"),
  GSE17907_GPL570 = file.path(base_dir, "GSE17907", "GSE17907-GPL570_series_matrix.txt.gz"),
  # GSE17907_GPL9128 skipped: GPL9128 has no NCBI annot.gz; SOFT annotation is ~3.75 GB.
  GSE19615 = file.path(base_dir, "GSE19615", "GSE19615_series_matrix.txt.gz"),
  GSE20194 = file.path(base_dir, "GSE20194", "GSE20194_series_matrix.txt.gz"),
  GSE22093 = file.path(base_dir, "GSE22093", "GSE22093_series_matrix.txt.gz"),
  GSE45255 = file.path(base_dir, "GSE45255", "GSE45255_series_matrix.txt.gz"),
  GSE20685 = file.path(base_dir, "GSE20685", "GSE20685_series_matrix.txt.gz"),
  GSE6532_GPL570 = file.path(base_dir, "GSE6532", "GSE6532-GPL570_series_matrix.txt.gz"),
  GSE6532_GPL96 = file.path(base_dir, "GSE6532", "GSE6532-GPL96_series_matrix.txt.gz"),
  GSE6532_GPL97 = file.path(base_dir, "GSE6532", "GSE6532-GPL97_series_matrix.txt.gz"),
  GSE11121 = file.path(base_dir, "GSE11121", "GSE11121_series_matrix.txt.gz"),
  GSE12093 = file.path(base_dir, "GSE12093", "GSE12093_series_matrix.txt.gz"),
  GSE2603 = file.path(base_dir, "GSE2603", "GSE2603_series_matrix.txt.gz"),
  GSE25066 = file.path(base_dir, "GSE25066", "GSE25066_series_matrix.txt.gz"),
  GSE25055 = file.path(base_dir, "GSE25055", "GSE25055_series_matrix.txt.gz"),
  GSE25065 = file.path(base_dir, "GSE25065", "GSE25065_series_matrix.txt.gz"),
  GSE42568 = file.path(base_dir, "GSE42568", "GSE42568_series_matrix.txt.gz"),
  GSE12276 = file.path(base_dir, "GSE12276", "GSE12276_series_matrix.txt.gz"),
  GSE103091 = file.path(base_dir, "GSE103091", "GSE103091_series_matrix.txt.gz"),
  GSE20711 = file.path(base_dir, "GSE20711", "GSE20711_series_matrix.txt.gz"),
  GSE21653 = file.path(base_dir, "GSE21653", "GSE21653_series_matrix.txt.gz"),
  GSE58812 = file.path(base_dir, "GSE58812", "GSE58812_series_matrix.txt.gz"),
  GSE88770 = file.path(base_dir, "GSE88770", "GSE88770_series_matrix.txt.gz"),
  GSE2990 = file.path(base_dir, "GSE2990", "GSE2990_series_matrix.txt.gz"),
  GSE9893 = file.path(base_dir, "GSE9893", "GSE9893_series_matrix.txt.gz"),
  GSE17705 = file.path(base_dir, "GSE17705", "GSE17705_series_matrix.txt.gz"),
  GSE26971 = file.path(base_dir, "GSE26971", "GSE26971_series_matrix.txt.gz"),
  GSE31448 = file.path(base_dir, "GSE31448", "GSE31448_series_matrix.txt.gz"),
  GSE48390 = file.path(base_dir, "GSE48390", "GSE48390_series_matrix.txt.gz"),
  GSE19783_GPL6480 = file.path(base_dir, "GSE19783", "GSE19783-GPL6480_series_matrix.txt.gz"),
  # GSE19783_GPL8227 skipped: GPL8227 is a miRNA platform and cannot map mRNA risk genes.
  GSE22219 = file.path(base_dir, "GSE22219", "GSE22219_series_matrix.txt.gz"),
  GSE22220_GPL6098 = file.path(base_dir, "GSE22220", "GSE22220-GPL6098_series_matrix.txt.gz"),
  GSE22220_GPL8178 = file.path(base_dir, "GSE22220", "GSE22220-GPL8178_series_matrix.txt.gz"),
  GSE24185 = file.path(base_dir, "GSE24185", "GSE24185_series_matrix.txt.gz"),
  GSE32646 = file.path(base_dir, "GSE32646", "GSE32646_series_matrix.txt.gz"),
  GSE58644 = file.path(base_dir, "GSE58644", "GSE58644_series_matrix.txt.gz")
)

stopifnot(all(file.exists(gse_files)))

strip_quote <- function(x) gsub('^"|"$', "", x)

num_value <- function(x) {
  x <- as.character(x)
  x[x %in% c("", "NA", "N/A", "na", "null", "NULL")] <- NA_character_
  pattern <- "-?\\d+(\\.\\d+)?"
  out <- rep(NA_real_, length(x))
  hit <- !is.na(x) & grepl(pattern, x)
  out[hit] <- suppressWarnings(as.numeric(sub(paste0(".*?(", pattern, ").*"), "\\1", x[hit])))
  out
}

event_value <- function(x) {
  z <- trimws(tolower(as.character(x)))
  out <- rep(NA_real_, length(z))
  out[z %in% c("1", "yes", "y", "true", "dead", "deceased", "event")] <- 1
  out[z %in% c("0", "no", "n", "false", "alive", "censored", "none")] <- 0
  need <- is.na(out)
  out[need] <- suppressWarnings(as.numeric(z[need]))
  need <- is.na(out)
  if (any(need)) out[need] <- num_value(z[need])
  out
}

convert_time_to_months <- function(time, unit) {
  unit <- tolower(as.character(unit))
  out <- as.numeric(time)
  out[unit == "years"] <- out[unit == "years"] * 12
  out[unit == "days"] <- out[unit == "days"] / 30.4375
  out
}

platform_bucket <- function(platform) {
  num <- as.integer(sub("^GPL", "", platform))
  if (is.na(num)) stop("Bad platform id: ", platform)
  if (num < 1000) return("GPLnnn")
  paste0("GPL", floor(num / 1000), "nnn")
}

get_platform_annotation <- function(platform) {
  annot_path <- file.path(base_dir, paste0(platform, ".annot.gz"))
  if (file.exists(annot_path) && file.info(annot_path)$size > 100000) {
    return(annot_path)
  }
  
  annot_url <- paste0(
    "https://ftp.ncbi.nlm.nih.gov/geo/platforms/",
    platform_bucket(platform), "/", platform, "/annot/", platform, ".annot.gz"
  )
  message("Downloading ", platform, " annotation: ", annot_url)
  ok <- tryCatch({
    utils::download.file(annot_url, annot_path, mode = "wb", quiet = FALSE)
    file.exists(annot_path) && file.info(annot_path)$size > 100000
  }, error = function(e) FALSE, warning = function(w) FALSE)
  if (ok) return(annot_path)
  
  soft_path <- file.path(base_dir, paste0(platform, "_family.soft.gz"))
  if (file.exists(soft_path) && file.info(soft_path)$size > 100000) {
    return(soft_path)
  }
  soft_url <- paste0(
    "https://ftp.ncbi.nlm.nih.gov/geo/platforms/",
    platform_bucket(platform), "/", platform, "/soft/", platform, "_family.soft.gz"
  )
  message("Annotation file not available for ", platform, "; trying SOFT: ", soft_url)
  utils::download.file(soft_url, soft_path, mode = "wb", quiet = FALSE)
  if (!file.exists(soft_path) || file.info(soft_path)$size < 100000) {
    stop("Cannot download annotation or SOFT file for ", platform,
         ". Tried: ", annot_url, " and ", soft_url)
  }
  soft_path
}

read_platform_gene_map <- function(platform, genes) {
  annot_path <- get_platform_annotation(platform)
  message("Using ", platform, " annotation: ", annot_path)
  con <- gzfile(annot_path, open = "rt")
  on.exit(close(con), add = TRUE)
  
  header <- NULL
  repeat {
    line <- readLines(con, n = 1, warn = FALSE)
    if (length(line) == 0) break
    if (strip_quote(strsplit(line, "\t", fixed = TRUE)[[1]][1]) == "ID") {
      header <- strip_quote(strsplit(line, "\t", fixed = TRUE)[[1]])
      break
    }
  }
  if (is.null(header)) stop("Cannot find ", platform, " annotation header.")
  
  id_col <- match("ID", header)
  symbol_col <- match("Gene Symbol", header)
  if (is.na(symbol_col)) symbol_col <- grep("gene.?symbol", header, ignore.case = TRUE)[1]
  if (is.na(symbol_col)) symbol_col <- match("Gene_Symbol", header)
  if (is.na(symbol_col)) symbol_col <- match("Symbol", header)
  if (is.na(id_col) || is.na(symbol_col)) {
    stop("Cannot identify ID/Gene Symbol/Gene_Symbol/Symbol columns in ", platform, " annotation.")
  }
  
  res_probe <- character()
  res_gene <- character()
  repeat {
    line <- readLines(con, n = 5000, warn = FALSE)
    if (length(line) == 0) break
    for (ln in line) {
      fields <- strip_quote(strsplit(ln, "\t", fixed = TRUE)[[1]])
      if (length(fields) < max(id_col, symbol_col)) next
      probe <- fields[id_col]
      sym_raw <- fields[symbol_col]
      if (is.na(sym_raw) || sym_raw == "" || sym_raw == "---") next
      syms <- unlist(strsplit(sym_raw, " /// |;|,"))
      syms <- trimws(syms)
      syms <- syms[syms %in% genes]
      if (length(syms) > 0) {
        res_probe <- c(res_probe, rep(probe, length(syms)))
        res_gene <- c(res_gene, syms)
      }
    }
  }
  gene_map <- unique(data.frame(probe_id = res_probe, gene = res_gene, stringsAsFactors = FALSE))
  if (nrow(gene_map) == 0) stop("No probes mapped to requested risk genes for ", platform, ".")
  write.csv(gene_map, file.path(out_dir, paste0(platform, "_probe_gene_map_risk10.csv")), row.names = FALSE)
  gene_map
}

read_series_platform <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  series_platform <- NA_character_
  sample_platform <- NA_character_
  repeat {
    line <- readLines(con, n = 1, warn = FALSE)
    if (length(line) == 0 || line == "!series_matrix_table_begin") break
    if (grepl("^!Series_platform_id", line)) {
      parts <- strip_quote(strsplit(line, "\t", fixed = TRUE)[[1]])
      if (length(parts) > 1) series_platform <- parts[2]
    }
    if (grepl("^!Sample_platform_id", line)) {
      parts <- strip_quote(strsplit(line, "\t", fixed = TRUE)[[1]])
      vals <- unique(parts[-1])
      vals <- vals[grepl("^GPL", vals)]
      if (length(vals) > 0) sample_platform <- vals[1]
    }
  }
  if (!is.na(sample_platform)) return(sample_platform)
  if (!is.na(series_platform)) return(series_platform)
  stop("Cannot detect platform for: ", path)
}

read_series_matrix_subset <- function(path, probe_ids) {
  message("Reading expression matrix: ", path)
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  
  header <- NULL
  repeat {
    line <- readLines(con, n = 1, warn = FALSE)
    if (length(line) == 0) break
    if (line == "!series_matrix_table_begin") {
      header_line <- readLines(con, n = 1, warn = FALSE)
      header <- strip_quote(strsplit(header_line, "\t", fixed = TRUE)[[1]])
      break
    }
  }
  if (is.null(header)) stop("Cannot find expression table in: ", path)
  samples <- header[-1]
  
  keep <- new.env(hash = TRUE, parent = emptyenv())
  for (p in probe_ids) keep[[p]] <- TRUE
  
  values <- list()
  repeat {
    lines <- readLines(con, n = 1000, warn = FALSE)
    if (length(lines) == 0) break
    for (ln in lines) {
      if (ln == "!series_matrix_table_end") break
      fields <- strip_quote(strsplit(ln, "\t", fixed = TRUE)[[1]])
      probe <- fields[1]
      if (!exists(probe, envir = keep, inherits = FALSE)) next
      vals <- suppressWarnings(as.numeric(fields[-1]))
      if (length(vals) == length(samples)) values[[probe]] <- vals
    }
    if (any(lines == "!series_matrix_table_end")) break
  }
  
  if (length(values) == 0) stop("No requested probe rows found in: ", path)
  mat <- do.call(rbind, values)
  colnames(mat) <- samples
  mat
}

normalize_key <- function(x) {
  x <- trimws(tolower(as.character(x)))
  x <- gsub("[^a-z0-9]+", "_", x)
  gsub("^_|_$", "", x)
}

read_series_clinical <- function(path, gse) {
  message("Reading clinical annotations: ", path)
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  
  records <- list()
  sample_count <- 0
  repeat {
    line <- readLines(con, n = 1, warn = FALSE)
    if (length(line) == 0 || line == "!series_matrix_table_begin") break
    if (!grepl("^!Sample_", line)) next
    
    parts <- strip_quote(strsplit(line, "\t", fixed = TRUE)[[1]])
    if (length(parts) < 2) next
    tag <- normalize_key(sub("^!Sample_", "", parts[1]))
    vals <- parts[-1]
    if (sample_count == 0) {
      sample_count <- length(vals)
      records <- replicate(sample_count, list(dataset = gse), simplify = FALSE)
    }
    if (length(vals) < sample_count) vals <- c(vals, rep(NA_character_, sample_count - length(vals)))
    vals <- vals[seq_len(sample_count)]
    
    for (i in seq_len(sample_count)) {
      val <- vals[i]
      if (tag == "characteristics_ch1" && grepl(":", val, fixed = TRUE)) {
        key <- normalize_key(sub(":.*$", "", val))
        v <- trimws(sub("^[^:]*:\\s*", "", val))
        records[[i]][[key]] <- v
      } else if (tag %in% c("geo_accession", "title", "source_name_ch1", "organism_ch1", "platform_id")) {
        records[[i]][[tag]] <- val
      }
    }
  }
  all_keys <- unique(unlist(lapply(records, names)))
  out <- as.data.frame(
    lapply(all_keys, function(k) vapply(records, function(r) if (!is.null(r[[k]])) r[[k]] else NA_character_, character(1))),
    stringsAsFactors = FALSE
  )
  names(out) <- all_keys
  out
}

get_col <- function(df, name) {
  if (name %in% names(df)) df[[name]] else rep(NA_character_, nrow(df))
}

make_endpoint_df <- function(clin, endpoint, time_col, event_col, unit) {
  data.frame(
    geo_accession = clin$geo_accession,
    endpoint = endpoint,
    time = num_value(get_col(clin, time_col)),
    event = event_value(get_col(clin, event_col)),
    time_unit = unit,
    stringsAsFactors = FALSE
  )
}

collapse_probe_to_gene <- function(probe_expr, gene_map) {
  gene_map <- gene_map[gene_map$probe_id %in% rownames(probe_expr), , drop = FALSE]
  out <- list()
  selected <- data.frame(gene = character(), probe_id = character(), stringsAsFactors = FALSE)
  for (g in unique(gene_map$gene)) {
    probes <- intersect(gene_map$probe_id[gene_map$gene == g], rownames(probe_expr))
    if (length(probes) == 1) {
      chosen <- probes
    } else if (length(probes) > 1) {
      iqr <- apply(probe_expr[probes, , drop = FALSE], 1, IQR, na.rm = TRUE)
      chosen <- probes[which.max(iqr)]
    } else {
      next
    }
    out[[g]] <- probe_expr[chosen, ]
    selected <- rbind(selected, data.frame(gene = g, probe_id = chosen, stringsAsFactors = FALSE))
  }
  mat <- do.call(rbind, out)
  rownames(mat) <- names(out)
  attr(mat, "selected_probes") <- selected
  mat
}

make_endpoint_table <- function(gse, clin) {
  gse_base <- sub("_GPL.*$", "", gse)
  if (gse_base == "GSE1456") {
    rbind(
      make_endpoint_df(clin, "OS", "surv_death", "death", "years"),
      make_endpoint_df(clin, "RFS", "surv_relapse", "relapse", "years"),
      make_endpoint_df(clin, "DSS", "surv_death", "death_bc", "years")
    )
  } else if (gse_base == "GSE4922") {
    make_endpoint_df(
      clin, "DFS", "dfs_time_yrs",
      "dfs_event_0_censored_1_event_defined_as_any_type_of_recurrence_local_regional_or_distant_or_death_from_breast_cancer",
      "years"
    )
  } else if (gse_base == "GSE17907") {
    make_endpoint_df(clin, "MFS", "mfsdel_month", "mfs", "months")
  } else if (gse_base == "GSE19615") {
    make_endpoint_df(clin, "DRFS", "distant_recurrence_free_survival_mo", "distant_recur_yn", "months")
  } else if (gse_base == "GSE103091") {
    rbind(
      make_endpoint_df(clin, "OS", "os_days", "death", "days"),
      make_endpoint_df(clin, "MFS", "mfs_days", "death", "days")
    )
  } else if (gse_base == "GSE25055" || gse_base == "GSE25065") {
    make_endpoint_df(clin, "DRFS", "drfs_even_time_years", "drfs_1_event_0_censored", "years")
  } else if (gse_base == "GSE2990") {
    make_endpoint_df(clin, "RFS", "time_rfs", "distant_rfs", "days")
  } else if (gse_base == "GSE9893") {
    make_endpoint_df(clin, "DMFS", "follow_up_period_months", "distant_metastases", "months")
  } else if (gse_base == "GSE17705") {
    make_endpoint_df(clin, "DMFS", "event_time_years", "distant_relapse_1_dr_0_censored", "years")
  } else if (gse_base == "GSE26971") {
    make_endpoint_df(clin, "MFS", "metastasis_free_interval_months", "metastasis_1", "months")
  } else if (gse_base == "GSE31448") {
    make_endpoint_df(clin, "DFS", "dfs_time_months", "dfs_evt", "months")
  } else if (gse_base == "GSE48390") {
    clin$event_drfs <- ifelse(grepl("disease-free", tolower(get_col(clin, "event"))), 0,
                              ifelse(grepl("recurrence|metastasis|mortality", tolower(get_col(clin, "event"))), 1, NA))
    make_endpoint_df(clin, "DRFS", "survival_time_year", "event_drfs", "years")
  } else if (gse_base == "GSE19783") {
    clin$bc_death_event <- ifelse(grepl("dead of bc", tolower(get_col(clin, "death_status"))), 1,
                                  ifelse(grepl("alive|dead of other", tolower(get_col(clin, "death_status"))), 0, NA))
    make_endpoint_df(clin, "BCSS", "disease_free_survival_time_months", "bc_death_event", "months")
  } else if (gse_base == "GSE22219" || gse_base == "GSE22220") {
    make_endpoint_df(clin, "DRFS", "distant_relapse_free_survival", "distant_relapse_event", "years")
  } else if (gse_base == "GSE58644") {
    make_endpoint_df(clin, "RFS", "time", "event", "months")
  } else
    if (gse_base == "GSE20685") {
      rbind(
        make_endpoint_df(clin, "OS", "follow_up_duration_years", "event_death", "years"),
        make_endpoint_df(clin, "MFS", "follow_up_duration_years", "event_metastasis", "years"),
        make_endpoint_df(clin, "RFS", "follow_up_duration_years", "regional_relapse", "years")
      )
    } else if (gse_base == "GSE20711") {
      rbind(
        make_endpoint_df(clin, "OS", "t_os", "e_os", "years"),
        make_endpoint_df(clin, "RFS", "t_rfs", "e_rfs", "years")
      )
    } else if (gse_base == "GSE21653") {
      make_endpoint_df(clin, "DFS", "dfs_time_months", "dfs_evt", "months")
    } else if (gse_base == "GSE58812") {
      rbind(
        make_endpoint_df(clin, "OS", "os_days", "death", "days"),
        make_endpoint_df(clin, "MFS", "mfs_days", "death", "days")
      )
    } else if (gse_base == "GSE88770") {
      rbind(
        make_endpoint_df(clin, "OS", "os_or_last_contact_years", "death", "years"),
        make_endpoint_df(clin, "DRFS", "drfs_or_last_contact_years", "drfs_event", "years")
      )
    } else if (gse_base == "GSE7390") {
      rbind(
        make_endpoint_df(clin, "OS", "t_os", "e_os", "days"),
        make_endpoint_df(clin, "RFS", "t_rfs", "e_rfs", "days"),
        make_endpoint_df(clin, "DMFS", "t_dmfs", "e_dmfs", "days")
      )
    } else if (gse_base == "GSE9195") {
      rbind(
        make_endpoint_df(clin, "RFS", "t_rfs", "e_rfs", "days"),
        make_endpoint_df(clin, "DMFS", "t_dmfs", "e_dmfs", "days")
      )
    } else if (gse_base == "GSE16446") {
      rbind(
        make_endpoint_df(clin, "OS", "os_time", "os_event", "days"),
        make_endpoint_df(clin, "DMFS", "dmfs_time", "dmfs_event", "days")
      )
    } else if (gse_base == "GSE45255") {
      rbind(
        make_endpoint_df(clin, "DFS", "dfs_time", "dfs_event_defined_as_any_type_of_recurrence_or_death_from_breast_cancer", "years"),
        make_endpoint_df(clin, "DMFS", "dmfs_time", "dmfs_event_defined_as_distant_metastasis_or_death_from_breast_cancer", "years"),
        make_endpoint_df(clin, "DSS", "dss_time", "dss_event_defined_as_death_from_breast_cancer", "years")
      )
    } else if (gse_base == "GSE6532") {
      rbind(
        make_endpoint_df(clin, "RFS", "t_rfs", "e_rfs", "days"),
        make_endpoint_df(clin, "DMFS", "t_dmfs", "e_dmfs", "days")
      )
    } else if (gse_base == "GSE11121") {
      make_endpoint_df(clin, "DMFS", "t_dmfs", "e_dmfs", "months")
    } else if (gse_base == "GSE12093") {
      make_endpoint_df(clin, "DFS", "dfs_time", "dfs_status", "months")
    } else if (gse_base == "GSE25066") {
      make_endpoint_df(clin, "DRFS", "drfs_even_time_years", "drfs_1_event_0_censored", "years")
    } else if (gse_base == "GSE42568") {
      rbind(
        make_endpoint_df(clin, "OS", "overall_survival_time_days", "overall_survival_event", "days"),
        make_endpoint_df(clin, "RFS", "relapse_free_survival_time_days", "relapse_free_survival_event", "days")
      )
    } else {
      data.frame(geo_accession = character(), endpoint = character(), time = numeric(),
                 event = numeric(), time_unit = character(), stringsAsFactors = FALSE)
    }
}

calc_risk_score <- function(gene_expr, coef_vec) {
  genes_use <- intersect(names(coef_vec), rownames(gene_expr))
  if (length(genes_use) < 3) {
    stop("Too few risk genes found in expression matrix: ", paste(genes_use, collapse = ", "))
  }
  z <- t(scale(t(gene_expr[genes_use, , drop = FALSE])))
  z[is.na(z)] <- 0
  score <- as.numeric(crossprod(coef_vec[genes_use], z))
  names(score) <- colnames(gene_expr)
  list(score = score, genes_used = genes_use)
}

km_one_endpoint <- function(score_df, endpoint_df, gse, endpoint) {
  ep <- endpoint_df[endpoint_df$endpoint == endpoint, ]
  ep$time <- convert_time_to_months(ep$time, ep$time_unit)
  ep$time_unit <- "months"
  dat <- merge(ep, score_df, by = "geo_accession", all = FALSE)
  dat <- dat[!is.na(dat$time) & !is.na(dat$event) & !is.na(dat$risk_score) & dat$time > 0, ]
  if (nrow(dat) < 20 || sum(dat$event == 1) < 5 || length(unique(dat$event)) < 2) {
    return(list(
      data = dat,
      stats = data.frame(dataset = gse, endpoint = endpoint, n = nrow(dat),
                         event = sum(dat$event == 1, na.rm = TRUE),
                         high_n = NA, low_n = NA, HR = NA, lower95 = NA,
                         upper95 = NA, cox_p = NA, logrank_p = NA,
                         median_cutoff = NA, cox_method = NA,
                         note = "too_few_events_or_invalid")
    ))
  }
  
  cutoff <- median(dat$risk_score, na.rm = TRUE)
  dat$risk_group <- ifelse(dat$risk_score > cutoff, "High", "Low")
  dat$risk_group <- factor(dat$risk_group, levels = c("Low", "High"))
  
  if (length(unique(dat$risk_group)) < 2) {
    return(list(
      data = dat,
      stats = data.frame(dataset = gse, endpoint = endpoint, n = nrow(dat),
                         event = sum(dat$event == 1, na.rm = TRUE),
                         high_n = sum(dat$risk_group == "High"),
                         low_n = sum(dat$risk_group == "Low"),
                         HR = NA, lower95 = NA, upper95 = NA, cox_p = NA,
                         logrank_p = NA, median_cutoff = cutoff,
                         cox_method = NA,
                         note = "one_group_only")
    ))
  }
  
  if (requireNamespace("coxphf", quietly = TRUE)) {
    fit_cox <- coxphf::coxphf(survival::Surv(time, event) ~ risk_group, data = dat)
    HR <- unname(exp(fit_cox$coefficients[1]))
    lower95 <- unname(exp(fit_cox$ci.lower[1]))
    upper95 <- unname(exp(fit_cox$ci.upper[1]))
    cox_p <- unname(fit_cox$prob[1])
    cox_method <- "Firth Cox"
  } else {
    fit_cox <- survival::coxph(survival::Surv(time, event) ~ risk_group, data = dat)
    s <- summary(fit_cox)
    ci <- s$conf.int[1, ]
    HR <- unname(ci["exp(coef)"])
    lower95 <- unname(ci["lower .95"])
    upper95 <- unname(ci["upper .95"])
    cox_p <- unname(s$coef[1, "Pr(>|z|)"])
    cox_method <- "Cox"
  }
  lr <- survival::survdiff(survival::Surv(time, event) ~ risk_group, data = dat)
  logrank_p <- 1 - pchisq(lr$chisq, df = length(lr$n) - 1)
  
  stats <- data.frame(
    dataset = gse,
    endpoint = endpoint,
    n = nrow(dat),
    event = sum(dat$event == 1, na.rm = TRUE),
    high_n = sum(dat$risk_group == "High"),
    low_n = sum(dat$risk_group == "Low"),
    HR = HR,
    lower95 = lower95,
    upper95 = upper95,
    cox_p = cox_p,
    logrank_p = logrank_p,
    median_cutoff = cutoff,
    cox_method = cox_method,
    note = "ok"
  )
  list(data = dat, stats = stats)
}

plot_km_base <- function(dat, stats, file_prefix) {
  if (nrow(dat) == 0 || !"risk_group" %in% names(dat) || length(unique(dat$risk_group)) < 2) return(invisible(NULL))
  
  fit <- survival::survfit(survival::Surv(time, event) ~ risk_group, data = dat)
  title_text <- paste0(stats$dataset, " ", stats$endpoint, " by 10-gene risk score")
  hr_text <- paste0(
    "Ref: Low risk\n",
    "n=", stats$n, ", events=", stats$event,
    "\nLow=", stats$low_n, ", High=", stats$high_n,
    "\n",
    "High risk: HR=", sprintf("%.2f", stats$HR),
    " (95% CI ", sprintf("%.2f", stats$lower95), "-",
    sprintf("%.2f", stats$upper95), ")",
    "\n", stats$cox_method, " P=",
    ifelse(stats$cox_p < 0.001, "<0.001", sprintf("%.3f", stats$cox_p)),
    "\nLog-rank P=",
    ifelse(stats$logrank_p < 0.001, "<0.001", sprintf("%.3f", stats$logrank_p))
  )
  
  p <- survminer::ggsurvplot(
    fit,
    data = dat,
    pval = TRUE,
    conf.int = FALSE,
    risk.table = TRUE,
    risk.table.height = 0.28,
    palette = c("#4DBBD5", "#E64B35"),
    xlab = paste0("Time (", unique(dat$time_unit)[1], ")"),
    ylab = "Survival probability",
    legend.title = "Risk group",
    legend.labs = c("Low risk", "High risk"),
    title = title_text
  )
  
  p$plot <- p$plot +
    ggplot2::theme_bw() +
    ggplot2::theme(
      aspect.ratio = 1,
      plot.title = ggplot2::element_text(hjust = 0.5, size = 11),
      axis.text = ggplot2::element_text(size = 9, color = "black"),
      axis.title = ggplot2::element_text(size = 10),
      legend.title = ggplot2::element_text(size = 9),
      legend.text = ggplot2::element_text(size = 8)
    ) +
    ggplot2::annotate(
      "text",
      x = max(dat$time, na.rm = TRUE) * 0.04,
      y = 0.08,
      label = hr_text,
      hjust = 0,
      vjust = 0,
      size = 2.0
    )
  
  p$table <- p$table +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 8, color = "black"),
      axis.text.y = ggplot2::element_text(size = 8, color = "black"),
      axis.title.x = ggplot2::element_text(size = 9),
      axis.title.y = ggplot2::element_text(size = 9)
    )
  
  pdf_file <- file.path(out_dir, paste0(file_prefix, ".pdf"))
  grDevices::pdf(pdf_file, width = 5.8, height = 7.0)
  print(p)
  grDevices::dev.off()
  pdf_file
}

append_csv <- function(x, file) {
  write.table(
    x,
    file = file,
    sep = ",",
    row.names = FALSE,
    col.names = !file.exists(file),
    append = file.exists(file),
    quote = TRUE
  )
}

summary_file <- file.path(out_dir, "KM_risk_score_summary.csv")
all_samples_file <- file.path(out_dir, "KM_risk_score_all_samples.csv")
manifest_file <- file.path(out_dir, "KM_PDF_manifest.csv")
process_log_file <- file.path(out_dir, "KM_processing_log.csv")
for (f in c(summary_file, all_samples_file, manifest_file, process_log_file)) {
  if (file.exists(f)) unlink(f)
}

gene_map_cache <- new.env(parent = emptyenv())

for (gse in names(gse_files)) {
  message("\n===== Processing ", gse, " =====")
  tryCatch({
    ## Read clinical annotations first. If no usable survival endpoint exists,
    ## skip expression loading to save memory and time.
    clin <- read_series_clinical(gse_files[[gse]], gse)
    write.csv(clin, file.path(out_dir, paste0(gse, "_clinical_extracted.csv")), row.names = FALSE)
    endpoints <- make_endpoint_table(gse, clin)
    
    if (nrow(endpoints) == 0 || length(unique(endpoints$endpoint)) == 0) {
      message("No recognized survival endpoint for ", gse, "; skip expression reading and KM.")
      append_csv(
        data.frame(dataset = gse, endpoint = NA, n = NA, event = NA,
                   note = "no_recognized_survival_endpoint", stringsAsFactors = FALSE),
        process_log_file
      )
    } else {
      message("Recognized endpoints for ", gse, ": ", paste(unique(endpoints$endpoint), collapse = ", "))
      
      platform <- read_series_platform(gse_files[[gse]])
      if (!exists(platform, envir = gene_map_cache, inherits = FALSE)) {
        assign(platform, read_platform_gene_map(platform, risk_genes), envir = gene_map_cache)
      }
      gene_map <- get(platform, envir = gene_map_cache, inherits = FALSE)
      
      probe_expr <- read_series_matrix_subset(gse_files[[gse]], unique(gene_map$probe_id))
      gene_expr <- collapse_probe_to_gene(probe_expr, gene_map)
      selected <- attr(gene_expr, "selected_probes")
      write.csv(selected, file.path(out_dir, paste0(gse, "_selected_probes_risk10.csv")), row.names = FALSE)
      
      risk <- calc_risk_score(gene_expr, risk_coef)
      score_df <- data.frame(
        geo_accession = names(risk$score),
        risk_score = as.numeric(risk$score),
        genes_used = paste(risk$genes_used, collapse = ";"),
        stringsAsFactors = FALSE
      )
      write.csv(score_df, file.path(out_dir, paste0(gse, "_risk_score_samples.csv")), row.names = FALSE)
      
      for (ep in unique(endpoints$endpoint)) {
        ans <- km_one_endpoint(score_df, endpoints, gse, ep)
        dat <- ans$data
        st <- ans$stats
        append_csv(st, summary_file)
        message("  ", ep, ": n=", st$n, ", events=", st$event, ", note=", st$note)
        append_csv(
          data.frame(dataset = gse, endpoint = ep, n = st$n,
                     event = st$event, note = st$note, stringsAsFactors = FALSE),
          process_log_file
        )
        if (nrow(dat) > 0) {
          dat$dataset <- gse
          dat$endpoint <- ep
          write.csv(dat, file.path(out_dir, paste0(gse, "_", ep, "_KM_input.csv")), row.names = FALSE)
          if (write_combined_sample_file) append_csv(dat, all_samples_file)
        }
        if (st$note == "ok") {
          pdf_file <- plot_km_base(dat, st, paste0("KM_", gse, "_", ep, "_risk_score"))
          append_csv(
            data.frame(
              dataset = gse,
              endpoint = ep,
              pdf_file = pdf_file,
              n = st$n,
              event = st$event,
              high_n = st$high_n,
              low_n = st$low_n,
              HR = st$HR,
              cox_p = st$cox_p,
              logrank_p = st$logrank_p,
              stringsAsFactors = FALSE
            ),
            manifest_file
          )
        }
        rm(ans, dat, st)
        gc(verbose = FALSE)
      }
    }
  }, error = function(e) {
    warning("Failed ", gse, ": ", conditionMessage(e))
    append_csv(
      data.frame(dataset = gse, endpoint = NA, n = NA, event = NA,
                 note = paste0("failed: ", conditionMessage(e)), stringsAsFactors = FALSE),
      process_log_file
    )
  })
  rm(list = intersect(c("probe_expr", "gene_expr", "selected", "risk", "score_df", "clin", "endpoints", "gene_map", "platform"), ls()))
  gc(verbose = FALSE)
  message("Finished ", gse, "; memory cleaned.")
}
if (file.exists(summary_file)) {
  km_stats <- read.csv(summary_file, stringsAsFactors = FALSE, check.names = FALSE)
  km_stats$FDR_logrank <- p.adjust(km_stats$logrank_p, method = "BH")
  write.csv(km_stats, summary_file, row.names = FALSE)
  rm(km_stats)
  gc(verbose = FALSE)
}

message("Done. KM risk-score results saved to: ", out_dir)













