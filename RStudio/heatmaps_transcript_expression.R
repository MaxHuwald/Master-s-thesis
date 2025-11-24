setwd("/project/sysviro/users/Max/analyses/Conditions/Results/Heatmaps")

library(dplyr)
library(edgeR)
library(pheatmap)
library(viridisLite)
# ---- Helper: extract condition name from filename ----
my_cols <- (inferno(100))

get_condition_name <- function(filepath) {
  filename <- basename(filepath)
  sub("\\..*", "", filename)  # take everything before first dot
}

# --- Helper to order rows by peak condition (and strength) ---
# primary_mat:   matrix used to detect peak condition per transcript (rows)
# col_order:     vector with the desired left->right order of conditions (must match column order)
# secondary_mat: optional matrix used to break ties and rank within each peak group by absolute strength
order_rows_by_peak <- function(primary_mat, col_order, secondary_mat = NULL) {
  stopifnot(all(colnames(primary_mat) == col_order[seq_along(colnames(primary_mat))]))
  peak_idx  <- max.col(primary_mat, ties.method = "first")
  peak_col  <- colnames(primary_mat)[peak_idx]
  col_rank  <- setNames(seq_along(col_order), col_order)
  primary_key <- as.integer(col_rank[peak_col])
  
  if (!is.null(secondary_mat)) {
    stopifnot(all(dim(secondary_mat) == dim(primary_mat)),
              all(colnames(secondary_mat) == colnames(primary_mat)),
              all(rownames(secondary_mat) == rownames(primary_mat)))
    peak_strength <- secondary_mat[cbind(seq_len(nrow(primary_mat)), peak_idx)]
  } else {
    peak_strength <- primary_mat[cbind(seq_len(nrow(primary_mat)), peak_idx)]
  }
  
  ord <- order(primary_key, -peak_strength, rownames(primary_mat))
  primary_mat[ord, , drop = FALSE]
}

# ---- Main function (absolute color scale; no row scaling) ----
generate_transcript_heatmap <- function(
    input_dir,
    output_prefix,
    file_pattern = "_counts.*\\.tsv$",
    desired_order = c("CHX", "2h", "4h DMSO", "4h STM2457", "6h", "8h"),
    top_n = 114,
    rename_map = c("4h_DMSO" = "4h DMSO", "4h_STM2457" = "4h STM2457"),
    width = 12, height = 12,
    order_by_peak = TRUE   # NEW: order rows by peak condition (kinetics)
) {
  message("Scanning: ", input_dir)
  files <- list.files(input_dir, pattern = file_pattern, full.names = TRUE)
  if (length(files) == 0) stop("No files matched in: ", input_dir)
  message("Found ", length(files), " count files.")
  
  # ---- Read single file ----
  read_counts <- function(f) {
    df <- read.table(f, sep = "\t", header = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
    if (ncol(df) == 4) {
      colnames(df) <- c("transcript_id", "length", "mapped_reads", "unmapped_reads")
      df <- df[, c("transcript_id", "mapped_reads")]
    } else if (ncol(df) == 2) {
      colnames(df) <- c("transcript_id", "mapped_reads")
    } else stop(paste("Unexpected number of columns in file:", f))
    cond <- get_condition_name(f)
    colnames(df)[2] <- cond
    df
  }
  
  # ---- Merge all files ----
  merged <- read_counts(files[1])
  if (length(files) > 1) {
    for (f in files[-1]) merged <- dplyr::full_join(merged, read_counts(f), by = "transcript_id")
  }
  merged <- merged %>% dplyr::filter(transcript_id != "*")
  
  # ---- Count matrix ----
  count_mat <- merged
  rownames(count_mat) <- count_mat$transcript_id
  count_mat <- count_mat[, -1, drop = FALSE]
  count_mat[] <- lapply(count_mat, function(x) as.numeric(x))
  count_mat[is.na(count_mat)] <- 0
  rownames(count_mat) <- merged$transcript_id
  
  # ---- Clean column names and reorder ----
  colnames(count_mat) <- dplyr::recode(colnames(count_mat), !!!rename_map, .default = colnames(count_mat))
  existing <- desired_order[desired_order %in% colnames(count_mat)]
  if (length(existing) == 0) stop("None of the desired_order columns found in data.")
  count_mat <- count_mat[, existing, drop = FALSE]
  
  # ---- Normalize ----
  cpm_mat <- edgeR::cpm(count_mat, log = FALSE)
  log_cpm <- log2(cpm_mat + 1)
  
  # ---- Top variable transcripts ----
  var_genes <- apply(log_cpm, 1, var)
  if (length(var_genes) < top_n) top_n <- length(var_genes)
  top_transcripts <- names(sort(var_genes, decreasing = TRUE))[1:top_n]
  heatmap_data <- log_cpm[top_transcripts, , drop = FALSE]
  
  # ---- Order rows: by peak column (kinetics), then by peak strength ----
  if (order_by_peak) {
    heatmap_data <- order_rows_by_peak(
      primary_mat   = heatmap_data,
      col_order     = colnames(heatmap_data),
      secondary_mat = heatmap_data  # absolute peak strength in log_cpm
    )
  } else {
    heatmap_data <- heatmap_data[order(rownames(heatmap_data)), , drop = FALSE]
  }
  
  # ---- Auto-title based on eU / non-eU ----
  if (grepl("non_eU", input_dir, ignore.case = TRUE) || grepl("non_eU", output_prefix, ignore.case = TRUE)) {
    dataset_label <- "non-eU-filtered"
  } else if (grepl("eU", input_dir, ignore.case = TRUE) || grepl("eU", output_prefix, ignore.case = TRUE)) {
    dataset_label <- "eU-filtered"
  } else {
    dataset_label <- "Total"
  }
  plot_title <- paste("HSV-2 Transcript Expression Across Conditions (", dataset_label, ")", sep = "")
  
  # ---- Manual relabeling for specific transcripts (rows) ----
  manual_map <- c(
    "mRNA.TU.12-TSS.35.1-0--96.52" = "ORF_UL36-1",
    "mRNA.TU.12-TSS.35.2-2--3.26" = "ORF_UL36-2",
    "mRNA.TU.21-TSS.53.1-0--9.27" = "ORF_RL2-5"
  )
  
  rn <- rownames(heatmap_data)
  for (old in names(manual_map)) {
    rn[rn == old] <- manual_map[[old]]
  }
  rownames(heatmap_data) <- rn
  
  # ---- Plot ----
  pdf(paste0(output_prefix, ".pdf"), width = width, height = height)
  pheatmap(
    heatmap_data,
    scale = "none",
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize_row = ifelse(nrow(heatmap_data) > 80, 6, 10),
    fontsize_col = 14,
    angle_col = 0,
    color = colorRampPalette(c("#FFFFCC", "#FFFFFF", "#762A83"))(100),
    main = plot_title
  )
  dev.off()
  
  message("Wrote: ", paste0(output_prefix, ".pdf"))
}

# ==========================================
# Row-scaled version: per-transcript normalization (0–1 across columns)
# ==========================================
generate_transcript_heatmap_rowscaled <- function(
    input_dir,
    output_prefix,
    file_pattern = "_counts.*\\.tsv$",
    desired_order = c("CHX", "2h", "4h DMSO", "4h STM2457", "6h", "8h"),
    top_n = 114,
    rename_map = c("4h_DMSO" = "4h DMSO", "4h_STM2457" = "4h STM2457"),
    width = 12, height = 12,
    order_by_peak = TRUE   # NEW: order rows by peak condition (kinetics)
) {
  message("Scanning: ", input_dir)
  files <- list.files(input_dir, pattern = file_pattern, full.names = TRUE)
  if (length(files) == 0) stop("No files matched in: ", input_dir)
  message("Found ", length(files), " count files.")
  
  # ---- Read single file ----
  read_counts <- function(f) {
    df <- read.table(f, sep = "\t", header = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
    if (ncol(df) == 4) {
      colnames(df) <- c("transcript_id", "length", "mapped_reads", "unmapped_reads")
      df <- df[, c("transcript_id", "mapped_reads")]
    } else if (ncol(df) == 2) {
      colnames(df) <- c("transcript_id", "mapped_reads")
    } else stop(paste("Unexpected number of columns in file:", f))
    cond <- get_condition_name(f)
    colnames(df)[2] <- cond
    df
  }
  
  # ---- Merge all files ----
  merged <- read_counts(files[1])
  if (length(files) > 1) {
    for (f in files[-1]) merged <- dplyr::full_join(merged, read_counts(f), by = "transcript_id")
  }
  merged <- merged %>% dplyr::filter(transcript_id != "*")
  
  # ---- Count matrix ----
  count_mat <- merged
  rownames(count_mat) <- count_mat$transcript_id
  count_mat <- count_mat[, -1, drop = FALSE]
  count_mat[] <- lapply(count_mat, function(x) as.numeric(x))
  count_mat[is.na(count_mat)] <- 0
  rownames(count_mat) <- merged$transcript_id
  
  # ---- Clean column names and reorder ----
  colnames(count_mat) <- dplyr::recode(colnames(count_mat), !!!rename_map, .default = colnames(count_mat))
  existing <- desired_order[desired_order %in% colnames(count_mat)]
  if (length(existing) == 0) stop("None of the desired_order columns found in data.")
  count_mat <- count_mat[, existing, drop = FALSE]
  
  # ---- Normalize ----
  cpm_mat <- edgeR::cpm(count_mat, log = FALSE)
  log_cpm <- log2(cpm_mat + 1)
  
  # ---- Top variable transcripts (on log_cpm) ----
  var_genes <- apply(log_cpm, 1, var)
  if (length(var_genes) < top_n) top_n <- length(var_genes)
  top_transcripts <- names(sort(var_genes, decreasing = TRUE))[1:top_n]
  heatmap_data <- log_cpm[top_transcripts, , drop = FALSE]
  
  # ---- Row-wise scaling to [0,1] per transcript ----
  row_min <- apply(heatmap_data, 1, min, na.rm = TRUE)
  row_max <- apply(heatmap_data, 1, max, na.rm = TRUE)
  row_rng <- row_max - row_min
  row_rng[row_rng == 0] <- 1  # avoid division by zero for flat rows
  heatmap_data_scaled <- sweep(sweep(heatmap_data, 1, row_min, "-"), 1, row_rng, "/")
  
  # ---- Order rows: find peak on scaled, rank within group by absolute peak strength from unscaled ----
  if (order_by_peak) {
    heatmap_data_scaled <- order_rows_by_peak(
      primary_mat   = heatmap_data_scaled,     # where is the peak (0..1)
      col_order     = colnames(heatmap_data_scaled),
      secondary_mat = heatmap_data             # rank by absolute peak in log_cpm
    )
  } else {
    heatmap_data_scaled <- heatmap_data_scaled[order(rownames(heatmap_data_scaled)), , drop = FALSE]
  }
  
  # ---- Auto-title based on eU / non-eU + row-scaled tag ----
  if (grepl("non_eU", input_dir, ignore.case = TRUE) || grepl("non_eU", output_prefix, ignore.case = TRUE)) {
    dataset_label <- "non-eU-filtered"
  } else if (grepl("eU", input_dir, ignore.case = TRUE) || grepl("eU", output_prefix, ignore.case = TRUE)) {
    dataset_label <- "eU-filtered"
  } else {
    dataset_label <- "Total"
  }
  plot_title <- paste0("HSV-2 transcript expression across conditions (", dataset_label, "; row-scaled per transcript)")
  
  # ---- Manual relabeling for specific transcripts (rows) ----
  manual_map <- c(
    "mRNA.TU.12-TSS.35.1-0--96.52" = "ORF_UL36-1",
    "mRNA.TU.12-TSS.35.2-2--3.26" = "ORF_UL36-2",
    "mRNA.TU.21-TSS.53.1-0--9.27" = "ORF_RL2-5"
  )
  
  rn <- rownames(heatmap_data_scaled)
  for (old in names(manual_map)) {
    rn[rn == old] <- manual_map[[old]]
  }
  rownames(heatmap_data_scaled) <- rn
  
  # ---- Plot ----
  pdf(paste0(output_prefix, "_rowscaled.pdf"), width = width, height = height)
  pheatmap::pheatmap(
    heatmap_data_scaled,
    scale = "none",
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize_row = ifelse(nrow(heatmap_data_scaled) > 80, 6, 10),
    fontsize_col = 14,
    angle_col = 0,
    color = colorRampPalette(c("#FFFFCC", "#FFFFFF", "#762A83"))(100),
    main = plot_title
  )
  dev.off()
  
  message("Wrote: ", paste0(output_prefix, "_rowscaled.pdf"))
}

# ==========================================
# Row-scaled version with FLIPPED axes:
# rows = conditions, columns = transcripts
# ==========================================
generate_transcript_heatmap_rowscaled_flipped <- function(
    input_dir,
    output_prefix,
    file_pattern = "_counts.*\\.tsv$",
    desired_order = c("CHX", "2h", "4h DMSO", "4h STM2457", "6h", "8h"),
    top_n = 114,
    rename_map = c("4h_DMSO" = "4h DMSO", "4h_STM2457" = "4h STM2457"),
    width = 20, height = 10,
    order_by_peak = TRUE   # keep the same kinetics-based ordering of transcripts
) {
  message("Scanning: ", input_dir)
  files <- list.files(input_dir, pattern = file_pattern, full.names = TRUE)
  if (length(files) == 0) stop("No files matched in: ", input_dir)
  message("Found ", length(files), " count files.")
  
  # ---- Read single file ----
  read_counts <- function(f) {
    df <- read.table(f, sep = "\t", header = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
    if (ncol(df) == 4) {
      colnames(df) <- c("transcript_id", "length", "mapped_reads", "unmapped_reads")
      df <- df[, c("transcript_id", "mapped_reads")]
    } else if (ncol(df) == 2) {
      colnames(df) <- c("transcript_id", "mapped_reads")
    } else stop(paste("Unexpected number of columns in file:", f))
    cond <- get_condition_name(f)
    colnames(df)[2] <- cond
    df
  }
  
  # ---- Merge all files ----
  merged <- read_counts(files[1])
  if (length(files) > 1) {
    for (f in files[-1]) merged <- dplyr::full_join(merged, read_counts(f), by = "transcript_id")
  }
  merged <- merged %>% dplyr::filter(transcript_id != "*")
  
  # ---- Count matrix ----
  count_mat <- merged
  rownames(count_mat) <- count_mat$transcript_id
  count_mat <- count_mat[, -1, drop = FALSE]
  count_mat[] <- lapply(count_mat, function(x) as.numeric(x))
  count_mat[is.na(count_mat)] <- 0
  rownames(count_mat) <- merged$transcript_id
  
  # ---- Clean column names and reorder ----
  colnames(count_mat) <- dplyr::recode(colnames(count_mat), !!!rename_map, .default = colnames(count_mat))
  existing <- desired_order[desired_order %in% colnames(count_mat)]
  if (length(existing) == 0) stop("None of the desired_order columns found in data.")
  count_mat <- count_mat[, existing, drop = FALSE]
  
  # ---- Normalize ----
  cpm_mat <- edgeR::cpm(count_mat, log = FALSE)
  log_cpm <- log2(cpm_mat + 1)
  
  # ---- Top variable transcripts (on log_cpm) ----
  var_genes <- apply(log_cpm, 1, var)
  if (length(var_genes) < top_n) top_n <- length(var_genes)
  top_transcripts <- names(sort(var_genes, decreasing = TRUE))[1:top_n]
  heatmap_data <- log_cpm[top_transcripts, , drop = FALSE]
  
  # ---- Row-wise scaling to [0,1] per transcript ----
  row_min <- apply(heatmap_data, 1, min, na.rm = TRUE)
  row_max <- apply(heatmap_data, 1, max, na.rm = TRUE)
  row_rng <- row_max - row_min
  row_rng[row_rng == 0] <- 1
  heatmap_data_scaled <- sweep(sweep(heatmap_data, 1, row_min, "-"), 1, row_rng, "/")
  
  # ---- Keep same kinetics-based ordering of transcripts BEFORE transpose ----
  if (order_by_peak) {
    # use scaled matrix to find the peak condition, but rank within group by absolute peak in log_cpm
    ordered_scaled <- order_rows_by_peak(
      primary_mat   = heatmap_data_scaled,
      col_order     = colnames(heatmap_data_scaled),
      secondary_mat = heatmap_data
    )
  } else {
    ordered_scaled <- heatmap_data_scaled[order(rownames(heatmap_data_scaled)), , drop = FALSE]
  }
  
  # ---- Flip axes: rows = conditions, columns = transcripts ----
  heatmap_data_scaled_t <- t(ordered_scaled)
  
  # ---- Auto-title ----
  if (grepl("non_eU", input_dir, ignore.case = TRUE) || grepl("non_eU", output_prefix, ignore.case = TRUE)) {
    dataset_label <- "non-eU-filtered"
  } else if (grepl("eU", input_dir, ignore.case = TRUE) || grepl("eU", output_prefix, ignore.case = TRUE)) {
    dataset_label <- "eU-filtered"
  } else {
    dataset_label <- "Total"
  }
  plot_title <- paste0("HSV-2 transcript expression across conditions (", dataset_label, "; row-scaled)")
  
  # ---- Manual relabeling for specific transcripts (columns) ----
  manual_map <- c(
    "mRNA.TU.12-TSS.35.1-0--96.52" = "ORF_UL36-1",
    "mRNA.TU.12-TSS.35.2-2--3.26" = "ORF_UL36-2",
    "mRNA.TU.21-TSS.53.1-0--9.27" = "ORF_RL2-5"
  )
  
  cn <- colnames(heatmap_data_scaled_t)
  for (old in names(manual_map)) {
    cn[cn == old] <- manual_map[[old]]
  }
  colnames(heatmap_data_scaled_t) <- cn
  
  
  # ---- Plot ----
  pdf(paste0(output_prefix, "_rowscaled_flipped.pdf"), width = width, height = height)
  pheatmap::pheatmap(
    heatmap_data_scaled_t,
    scale = "none",
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_rownames = TRUE,   # now shows condition names (rows)
    show_colnames = TRUE,   # now shows transcript IDs (columns)
    fontsize_row = 14,      # conditions: usually fewer; make larger
    fontsize_col = ifelse(ncol(heatmap_data_scaled_t) > 80, 11, 12),
    angle_col = 90,         # rotate transcript IDs for readability
    color = my_cols,
    main = plot_title
  )
  dev.off()
  
  message("Wrote: ", paste0(output_prefix, "_rowscaled_flipped.pdf"))
}

# ==========================================
# Example calls
# ==========================================

generate_transcript_heatmap(
  input_dir = "/project/sysviro/users/Max/analyses/Conditions/transcriptome/non_eU/counts",
  output_prefix = "HSV2_transcript_heatmap_non_eU",
  desired_order = c("CHX", "2h", "4h DMSO", "4h STM2457", "6h", "8h"),
  top_n = 114
)

generate_transcript_heatmap(
  input_dir = "/project/sysviro/users/Max/analyses/Conditions/transcriptome/eU/counts",
  output_prefix = "HSV2_transcript_heatmap_eU",
  desired_order = c("CHX", "2h", "4h DMSO", "4h STM2457", "6h", "8h"),
  top_n = 114
)

generate_transcript_heatmap(
  input_dir = "/project/sysviro/users/Max/analyses/Conditions/transcriptome/counts",
  output_prefix = "HSV2_transcript_heatmap_total",
  desired_order = c("CHX", "2h", "4h DMSO", "4h STM2457", "6h", "8h"),
  top_n = 114
)

generate_transcript_heatmap_rowscaled_flipped(
  input_dir = "/project/sysviro/users/Max/analyses/Conditions/transcriptome/non_eU/counts",
  output_prefix = "HSV2_transcript_heatmap_non_eU",
  desired_order = c("CHX", "2h", "4h DMSO", "4h STM2457", "6h", "8h"),
  top_n = 114
)

generate_transcript_heatmap_rowscaled_flipped(
  input_dir = "/project/sysviro/users/Max/analyses/Conditions/transcriptome/eU/counts",
  output_prefix = "HSV2_transcript_heatmap_eU",
  desired_order = c("CHX", "2h", "4h DMSO", "4h STM2457", "6h", "8h"),
  top_n = 114
)

generate_transcript_heatmap_rowscaled(
  input_dir = "/project/sysviro/users/Max/analyses/Conditions/transcriptome/counts",
  output_prefix = "HSV2_transcript_heatmap_total",
  desired_order = c("CHX", "2h", "4h DMSO", "4h STM2457", "6h", "8h"),
  top_n = 114
)

generate_transcript_heatmap_rowscaled_flipped(
  input_dir = "/project/sysviro/users/Max/analyses/Conditions/transcriptome/counts",
  output_prefix = "HSV2_transcript_heatmap_total",
  desired_order = c("CHX", "2h", "4h DMSO", "4h STM2457", "6h", "8h"),
  top_n = 114
)
