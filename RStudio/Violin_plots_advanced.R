setwd("/project/sysviro/users/Max/analyses/Conditions/Results/Violinplots")
library(Rsamtools)
library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
# -------------------------------
# 1. Define directories
# -------------------------------
bam_dirs <- c(
  HSV2    = "/project/sysviro/users/Max/analyses/Conditions/trimmed",
  hg38    = "/project/sysviro/users/Max/analyses/Conditions/hg38"
)

RNAKINET_DIR <- "/project/sysviro/users/Max/analyses/RNAkinet/"

# -------------------------------
# 2. List BAM and CSV files
# -------------------------------
bam_files <- lapply(names(bam_dirs), function(cond) {
  files <- list.files(bam_dirs[[cond]], pattern = "\\.bam$", full.names = TRUE)
  data.frame(file_bam = files, Condition = cond, stringsAsFactors = FALSE)
}) %>% bind_rows()

csv_files <- list.files(RNAKINET_DIR, pattern = "\\.csv$", full.names = TRUE)

# -------------------------------
# 3. Function to extract sample name from BAM or CSV
# -------------------------------
get_sample_name <- function(fname) {
  # Use everything before first dot as sample
  base <- basename(fname)
  strsplit(base, "\\.")[[1]][1]
}

bam_files$Sample <- sapply(bam_files$file_bam, get_sample_name)
csv_samples <- sapply(csv_files, get_sample_name)
csv_files_df <- data.frame(file_csv = csv_files, Sample = csv_samples, stringsAsFactors = FALSE)

# -------------------------------
# 4. Merge BAM and CSV files by Sample
# -------------------------------
pairs_df <- merge(bam_files, csv_files_df, by = "Sample", all.x = TRUE)

# -------------------------------
# 5. Function to extract BAM tags and merge with CSV
# -------------------------------
extract_bam_data <- function(bam_file, csv_file, condition, sample, eu_cutoff = 0.9) {
  
  # --- Read BAM tags ---
  param <- ScanBamParam(what = c("qname"), tag = c("pt", "XB"))
  bam_tags <- scanBam(BamFile(bam_file), param = param)[[1]]
  
  n_reads <- length(bam_tags$qname)
  
  # Fill missing tags with NA if not present
  pt_tag <- if (!is.null(bam_tags$tag$pt)) bam_tags$tag$pt else rep(NA, n_reads)
  XB_tag <- if (!is.null(bam_tags$tag$XB)) bam_tags$tag$XB else rep(NA, n_reads)
  
  df_bam <- data.frame(
    read_id = bam_tags$qname,
    pt = pt_tag,
    XB = XB_tag,
    stringsAsFactors = FALSE
  )
  
  # --- Merge CSV if exists ---
  if (!is.na(csv_file) && file.exists(csv_file)) {
    df_csv <- fread(csv_file)
    setnames(df_csv, c("read_id","5eu_mod_score","5eu_modified_prediction"))
    
    df <- df_bam %>%
      left_join(df_csv, by = "read_id") %>%
      mutate(eU_label = ifelse(!is.na(`5eu_mod_score`) & `5eu_mod_score` >= eu_cutoff, "new mRNA", "old mRNA"))
    
  } else {
    df <- df_bam
    df$`5eu_mod_score` <- NA
    df$eU_label <- NA
  }
  
  df$Condition <- condition
  df$Sample <- sample
  return(df)
}

# -------------------------------
# 6. Apply to all BAMs
# -------------------------------
all_data <- mapply(
  extract_bam_data,
  pairs_df$file_bam,
  pairs_df$file_csv,
  pairs_df$Condition,
  pairs_df$Sample,
  SIMPLIFY = FALSE
) %>% bind_rows()

# -------------------------------
# 7. Example: plotting pt by sample and eU
# -------------------------------

all_data <- all_data %>%
  mutate(Condition = recode(Condition,
                         "HSV2" = "HSV-2"))

#rename 4h conditions
all_data <- all_data %>%
  mutate(Sample = recode(Sample,
                         "4h_DMSO" = "4h DMSO",
                         "4h_STM2457" = "4h STM2457"))

# Create a combined grouping variable
all_data$Group <- paste(all_data$Condition, all_data$eU_label)

#Create subsets:
all_data_new_mRNA <- all_data %>%
  filter(eU_label == "new mRNA", Sample != "4h", Sample !="CHX")

all_data_no_4h <- all_data %>%
  filter(Sample != "4h")
#remove eU CHX
all_data_no_4h <- all_data_no_4h %>%
  filter(!(Sample == "CHX" & eU_label == "new mRNA"))


# Make Label a factor in the desired order
all_data_new_mRNA$Sample <- factor(all_data_new_mRNA$Sample,
                          levels = c("CHX", "2h", "4h DMSO", "4h STM2457", "6h", "8h"))


# Make Label a factor in the desired order
all_data$Sample <- factor(all_data$Sample,
                         levels = c("CHX", "2h","4h", "4h DMSO", "4h STM2457", "6h", "8h"))

# Make Label a factor in the desired order
all_data_no_4h$Sample <- factor(all_data_no_4h$Sample,
                          levels = c("CHX", "2h", "4h DMSO", "4h STM2457", "6h", "8h"))

# Define desired order of Group levels:
desired_group_levels <- c(
  "hg38 old mRNA",
  "hg38 new mRNA",
  "HSV-2 old mRNA",
  "HSV-2 new mRNA"
)

# Convert Group to factor with this order
all_data_no_4h$Group <- factor(all_data_no_4h$Group, levels = desired_group_levels)


# Example: two dark colors from Oranges and Purples
condition_colors <- c("HSV-2" = brewer.pal(5, "Purples")[4],  # darkest purple
                      "hg38"    = brewer.pal(5, "Greens")[4])  # darkest orange

condition_colors1 <- c(
  "HSV-2 new mRNA" = brewer.pal(5, "Purples")[2],   # dark purple
  "HSV-2 old mRNA" = brewer.pal(5, "Purples")[5],   # lighter purple
  "hg38 new mRNA"    = brewer.pal(5, "Greens")[2],   # light green
  "hg38 old mRNA"    = brewer.pal(5, "Greens")[5]    # dark green
)

# Violin plot eU scores
p1 <- ggplot(all_data %>% filter(Sample != "4h"), aes(x = Sample, y = `5eu_mod_score`, fill = Condition)) +
  geom_violin(trim = FALSE, scale = "width",
              draw_quantiles = c(0.25, 0.5, 0.75),
              color = "black", linewidth = 0.5,
              position = position_dodge(width = 0.9), alpha = 0.8) +
  #geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.5,
  #             position = position_dodge(width = 0.9), color = "black") +
  scale_fill_manual(values = condition_colors) +
  labs(
    title = "Violin Plots of eU-labelling score by timepoint (width-scaled)",
    x = "Condition",
    y = "Score",
    fill = "Alignment"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),   # remove vertical major grid lines
    panel.grid.minor.x = element_blank(),   # remove vertical minor grid lines
    panel.grid.major.y = element_line(color = "grey80", linewidth = 0.4),
    panel.grid.minor.y = element_blank(),
    # ---- TEXT SIZE CONTROL ----
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text  = element_text(size = 16),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x  = element_text(size = 18),
    axis.text.y  = element_text(size = 18)) +
geom_hline(yintercept = 0.9, color = "red", linetype = "dotted", linewidth = .8)

# Violin plot pt scores
p2 <- ggplot(all_data, aes(x = Sample, y = pt, fill = Condition)) +
  geom_violin(trim = FALSE, scale = "width",
              draw_quantiles = c(0.25, 0.5, 0.75),
              color = "black", linewidth = 0.5,
              position = position_dodge(width = 0.9), alpha = 0.8) +
  #geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.5,
  #             position = position_dodge(width = 0.9), color = "black") +
  scale_fill_manual(values = condition_colors) +
  labs(
    title = "Violin Plots of polyA-tail length by timepoint",
    x = "Condition",
    y = "polyA-tail length",
    fill = "Alignment"
  ) +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text  = element_text(size = 12),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold"),
        axis.text.x  = element_text(size = 14, angle = 0, hjust = 0.5),
        axis.text.y  = element_text(size = 14))


# Violin plot
p3 <- ggplot(all_data_no_4h, aes(x = Sample, y = pt, fill = Group)) +
  geom_violin(trim = FALSE, scale = "width",
              draw_quantiles = c(0.25, 0.5, 0.75),
              color = "black", linewidth = 0.5,
              position = position_dodge(width = 0.9), alpha = 0.8) +
  scale_fill_manual(values = condition_colors1) +
  labs(
    title = "Violin Plots of polyA-tail length by Sample, Alignment, and mRNA age",
    x = "Condition",
    y = "polyA-tail length",
    fill = "Alignment + mRNA age"
  ) +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text  = element_text(size = 12),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold"),
        axis.text.x  = element_text(size = 14, angle = 0, hjust = 0.5),
        axis.text.y  = element_text(size = 14))

p4 <- ggplot(all_data %>% filter(Sample != "4h"), aes(x = Sample, y = `5eu_mod_score`, fill = Condition)) +
  geom_violin(trim = FALSE, scale = "count",
              draw_quantiles = c(0.25, 0.5, 0.75),
              color = "black", linewidth = 0.5,
              position = position_dodge(width = 0.9), alpha = 0.8) +
  #geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.5,
  #             position = position_dodge(width = 0.9), color = "black") +
  scale_fill_manual(values = condition_colors) +
  labs(
    title = "Violin Plots of eU-labelling score by timepoint (count-scaled)",
    x = "Condition",
    y = "Score",
    fill = "Alignment"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),   # remove vertical major grid lines
    panel.grid.minor.x = element_blank(),   # remove vertical minor grid lines
    panel.grid.major.y = element_line(color = "grey80", linewidth = 0.4),
    panel.grid.minor.y = element_blank(),
    # ---- TEXT SIZE CONTROL ----
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text  = element_text(size = 16),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x  = element_text(size = 18),
    axis.text.y  = element_text(size = 18)) +
  geom_hline(yintercept = 0.9, color = "red", linetype = "dotted", linewidth = .8)


ggsave("violin_eU_scores_total.pdf", plot = p1, width = 18, height = 5, dpi = 300)
ggsave("violin_eU_scores_total_count.pdf", plot = p4, width = 18, height = 5, dpi = 300)
ggsave("violin_polyA_total.pdf", plot = p2, width = 16, height = 6, dpi = 300)
ggsave("violin_pt_eU_condition_sample.pdf", plot = p3, width = 16, height = 6, dpi = 300)

combined_p14 <- p4 / p1  # "/" stacks vertically, "|" would stack horizontally
ggsave("violin_eU_scores_combined.png", plot = combined_p14,
       width = 18, height = 12, dpi = 300)

# -------------------------------
# 8. Separate violin plots for HSV2 and hg38
# -------------------------------

# Subset data
all_data_no_4h_HSV2 <- all_data_no_4h %>% filter(Condition == "HSV-2")
all_data_no_4h_hg38 <- all_data_no_4h %>% filter(Condition == "hg38")

# HSV2-only plot (p4)
p5 <- ggplot(all_data_no_4h_HSV2, aes(x = Sample, y = pt, fill = Group)) +
  geom_violin(trim = FALSE, scale = "width",
              draw_quantiles = c(0.25, 0.5, 0.75),
              color = "black", linewidth = 0.5,
              position = position_dodge(width = 0.9), alpha = 0.8) +
  scale_fill_manual(values = condition_colors1[c("HSV-2 old mRNA", "HSV-2 new mRNA")]) +
  labs(
    title = "HSV-2: poly(A)-tail length by condition and mRNA age (width-scaled)",
    x = "Condition",
    y = "poly(A)-tail length",
    fill = "mRNA age"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),   # remove vertical major grid lines
    panel.grid.minor.x = element_blank(),   # remove vertical minor grid lines
    panel.grid.major.y = element_line(color = "grey80", linewidth = 0.4),
    panel.grid.minor.y = element_blank(),
    # ---- TEXT SIZE CONTROL ----
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text  = element_text(size = 16),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x  = element_text(size = 18),
    axis.text.y  = element_text(size = 18))

# hg38-only plot (p5)
p6 <- ggplot(all_data_no_4h_hg38, aes(x = Sample, y = pt, fill = Group)) +
  geom_violin(trim = FALSE, scale = "width",
              draw_quantiles = c(0.25, 0.5, 0.75),
              color = "black", linewidth = 0.5,
              position = position_dodge(width = 0.9), alpha = 0.8) +
  scale_fill_manual(values = condition_colors1[c("hg38 old mRNA", "hg38 new mRNA")]) +
  labs(
    title = "hg38: poly(A)-tail length by condition and mRNA age (width-scaled)",
    x = "Condition",
    y = "poly(A)-tail length",
    fill = "mRNA age"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),   # remove vertical major grid lines
    panel.grid.minor.x = element_blank(),   # remove vertical minor grid lines
    panel.grid.major.y = element_line(color = "grey80", linewidth = 0.4),
    panel.grid.minor.y = element_blank(),
    # ---- TEXT SIZE CONTROL ----
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text  = element_text(size = 16),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x  = element_text(size = 18),
    axis.text.y  = element_text(size = 18))

# Save plots
combined_p56 <- p6 / p5  # "/" stacks vertically, "|" would stack horizontally
ggsave("violin_pt_eU_HSV2_hg38_combined.png", plot = combined_p56,
       width = 20, height = 10, dpi = 300)
