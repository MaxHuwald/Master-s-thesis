library(readr)
library(dplyr)
library(tidyr)
library(pheatmap)
library(viridisLite)

my_cols <- (viridis(100))

# --- 1) Folder ---
data_dir <- "//project/sysviro/users/Max/analyses/Conditions/halflives"

# --- 2) Files ---
files <- list.files(data_dir, pattern = "\\.halflives\\.csv$", full.names = TRUE)

# --- 3–5) Read, clean, add Condition, and bind ---
all_df <- lapply(files, function(f) {
  condition <- sub("\\.halflives\\.csv$", "", basename(f))
  
  raw <- readr::read_csv(f, show_col_types = FALSE)
  raw <- raw[, !grepl("^\\.\\.\\.|^$", names(raw)), drop = FALSE]
  
  raw %>%
    dplyr::select(transcript, percentage_modified, reads, pred_t5) %>%
    dplyr::mutate(Condition = condition)
}) %>%
  dplyr::bind_rows()

# --- 4) Rename conditions (cleanup) ---
all_df <- all_df %>%
  mutate(Condition = recode(Condition,
                            "4h_DMSO" = "4h DMSO",
                            "4h_STM2457" = "4h STM2457"))

# --- 5) Create wide matrix ---
mat_df <- all_df %>%
  dplyr::select(transcript, Condition, pred_t5) %>%
  tidyr::pivot_wider(names_from = Condition, values_from = pred_t5)

# --- 6) Convert to matrix and fill NAs ---
mat <- mat_df %>%
  tibble::column_to_rownames("transcript") %>%
  as.matrix()
mat[is.na(mat)] <- 0

# --- 7) Column order ---
desired_order <- c("2h", "4h DMSO", "4h STM2457", "6h", "8h")
present_order <- desired_order[desired_order %in% colnames(mat)]
mat <- mat[, present_order, drop = FALSE]

# --- 8) Optional manual transcript renames ---
manual_map <- c(
  "mRNA.TU.12-TSS.35.1-0--96.52" = "ORF_UL36-1",
  "mRNA.TU.12-TSS.35.2-2--3.26" = "ORF_UL36-2",
  "mRNA.TU.21-TSS.53.1-0--9.27" = "ORF_RL2-5"
)
rn <- rownames(mat)
for (old in names(manual_map)) {
  rn[rn == old] <- manual_map[[old]]
}
rownames(mat) <- rn

# --- 9) Row-wise min–max scaling (0–1 per transcript) ---
row_min <- apply(mat, 1, min, na.rm = TRUE)
row_max <- apply(mat, 1, max, na.rm = TRUE)
row_rng <- row_max - row_min
row_rng[row_rng == 0] <- 1  # avoid division by zero
mat_scaled <- sweep(sweep(mat, 1, row_min, "-"), 1, row_rng, "/")

# --- Flip heatmap orientation: make conditions rows, transcripts columns ---
mat_scaled_t <- t(mat_scaled)


# --- 10) Heatmap (row-wise scaled values) ---
pdf("transcript_half-lives_minmax_scaled.pdf", width = 20, height = 10)
pheatmap(
  mat_scaled_t,
  scale = "none",  # already scaled manually
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 14,
  fontsize_col = 11,
  angle_col = 90,
  color = my_cols,
  main = "Transcript half-lives (row-scaled per transcript)"
)
dev.off()
