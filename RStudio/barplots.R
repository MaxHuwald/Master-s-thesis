setwd("/project/sysviro/users/Max/analyses/Conditions/Results/Barplots")
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(stringr)
library(readr)
library(patchwork)   # for combining plots
library(scales)

# file paths
file_wdx3 <- "/project/sysviro/users/Max/WarpDemuX/HSV2_333_ARPE19_WDX-3/Demultiplexed/barcode_predictions.all.csv"
file_wdx4 <- "/project/sysviro/users/Max/WarpDemuX/HSV2_333_ARPE19_WDX-4/Demultiplexed/barcode_predictions.all.csv"

# read them
df_wdx3 <- read.csv(file_wdx3)
df_wdx4 <- read.csv(file_wdx4)

# add extra column to label
df_wdx3$Experiment <- "Run 3"
df_wdx4$Experiment <- "Run 4"

# combine into one data frame
df_all <- rbind(df_wdx3, df_wdx4)

# remove header rows
df_all <- df_all %>%
  filter(read_id != "read_id")   # remove rows where read_id equals header

df_all <- df_all %>%
  mutate(cond = case_when(
    Experiment == "Run 3" & predicted_barcode == 4 ~ "2h",
    Experiment == "Run 3" & predicted_barcode == 5 ~ "CHX",
    Experiment == "Run 4" & predicted_barcode == 4 ~ "4h_DMSO",
    Experiment == "Run 4" & predicted_barcode == 5 ~ "4h_STM2457",
    Experiment == "Run 4" & predicted_barcode == 6 ~ "6h",
    Experiment == "Run 4" & predicted_barcode == 11 ~ "8h",
    Experiment == "Run 3" & predicted_barcode == -1 ~ "Run 3 unassigned",
    Experiment == "Run 3" & predicted_barcode == 6 ~ "Run 3 unassigned",
    Experiment == "Run 3" & predicted_barcode == 11 ~ "Run 3 unassigned",
    Experiment == "Run 4" & predicted_barcode == -1 ~ "Run 4 unassigned",
    TRUE ~ NA_character_
  ))

### --- Histogram section ---
df_hist <- df_all %>%
  # rename for plotting consistency
  mutate(cond = recode(cond,
                       "CHX" = "Run 3 WDX-5 (CHX)",
                       "2h" = "Run 3 WDX-4 (2h)",
                       "4h_DMSO" = "Run 4 WDX-4 (4h DMSO)",
                       "4h_STM2457" = "Run 4 WDX-5 (4h STM2457)",
                       "6h" = "Run 4 WDX-6 (6h)",
                       "8h" = "Run 4 WDX-11 (8h)")) %>%
  # replace remaining NA with experiment-specific unassigned if needed
  mutate(cond = ifelse(is.na(cond) & Experiment == "Run 3", "Run 3 unassigned",
                       ifelse(is.na(cond) & Experiment == "Run 4", "Run 4 unassigned", cond))) %>%
  # convert numeric and add percent column for plotting
  mutate(confidence_score = as.numeric(confidence_score),
         confidence_score_pct = confidence_score * 100) %>%
  # set factor levels to ensure consistent ordering & color mapping
  mutate(cond = factor(cond, levels = c("Run 3 WDX-5 (CHX)", "Run 3 WDX-4 (2h)", "Run 4 WDX-4 (4h DMSO)", "Run 4 WDX-5 (4h STM2457)", "Run 4 WDX-6 (6h)", "Run 4 WDX-11 (8h)",
                                        "Run 3 unassigned", "Run 4 unassigned")))

# now proceed with histogram plotting using df_hist (which contains ALL reads)
cond_values <- levels(df_hist$cond)
bin_breaks <- seq(0, 100, by = 1)

### ---barplots section--- ###
# keep only >= 0.95 confidence
df_all <- df_all %>%
  mutate(cond = ifelse(confidence_score <= 0.95 | 
                         cond %in% c("Run 3 unassigned", "Run 4 unassigned"),
                       "low confidence", cond))

### alignments ###
align_dir <- "/project/sysviro/users/Max/analyses/Conditions/read_ids_per_alignment"
eu_dir <- "/project/sysviro/users/Max/analyses/Conditions/read_ids_per_eU"

read_alignment_file <- function(filepath, alignment_type) {
  cond <- strsplit(basename(filepath), "\\.")[[1]][1]
  read_ids <- read_tsv(filepath, col_names = "read_id", show_col_types = FALSE) %>%
    mutate(alignment = alignment_type, cond = cond)
  return(read_ids)
}

read_eu_file <- function(filepath, eu_type) {
  cond <- strsplit(basename(filepath), "\\.")[[1]][1]
  read_ids <- read_tsv(filepath, col_names = "read_id", show_col_types = FALSE) %>%
    mutate(eU_label = eu_type, cond = cond)
  return(read_ids)
}

files_hg38 <- list.files(file.path(align_dir, "hg38"), full.names = TRUE)
files_hsv2 <- list.files(file.path(align_dir, "trimmed"), full.names = TRUE)
files_eu <- list.files(file.path(eu_dir, "trimmed"), full.names = TRUE)

df_hg38 <- bind_rows(lapply(files_hg38, read_alignment_file, alignment_type = "hg38"))
df_hsv2 <- bind_rows(lapply(files_hsv2, read_alignment_file, alignment_type = "HSV-2"))
df_eu <- bind_rows(lapply(files_eu, read_eu_file, eu_type = "eU"))

df_alignments <- bind_rows(df_hg38, df_hsv2)

# drop 4h combined condition
df_alignments <- df_alignments %>%
  filter(cond != "4h")
df_eu <- df_eu %>%
  filter(cond != "4h")

# drop duplicates
df_alignments <- df_alignments %>%
  group_by(read_id, cond) %>%
  summarise(alignment = ifelse(n_distinct(alignment) > 1, "unaligned", first(alignment)),
            .groups = "drop")

# join with df_all
df_all_extended <- df_all %>%
  left_join(df_alignments, by = c("read_id", "cond")) %>%
  mutate(alignment = ifelse(is.na(alignment), "unaligned", alignment))

df_all_final <- df_all_extended %>%
  left_join(df_eu, by = c("read_id", "cond")) %>%
  mutate(eU_label = ifelse(is.na(eU_label), "non-eU", eU_label))

# rename treatments
df_all_final <- df_all_final %>%
  mutate(cond = recode(cond,
                       "4h_DMSO" = "4h DMSO",
                       "4h_STM2457" = "4h STM2457"))

df_alignment_plot <- df_all_final %>%
  filter(cond != "low confidence")

df_eu_plot <- df_all_final %>%
  filter(cond != "low confidence", eU_label == "eU")

df_all_final <- df_all_final %>%
  mutate(cond = factor(cond, levels = c("low confidence", "CHX", "2h", "4h DMSO", "4h STM2457", "6h", "8h")))

df_alignment_plot <- df_alignment_plot %>%
  mutate(cond = factor(cond, levels = c("CHX", "2h", "4h DMSO", "4h STM2457", "6h", "8h")))

df_alignment_plot <- df_alignment_plot %>%
  mutate(alignment = factor(alignment, levels = c("unaligned", "HSV-2", "hg38")))

df_eu_plot <- df_eu_plot %>%
  mutate(cond = factor(cond, levels = c("CHX", "2h", "4h DMSO", "4h STM2457", "6h", "8h")))

df_eu_plot <- df_eu_plot %>%
  mutate(alignment = factor(alignment, levels = c("unaligned", "HSV-2", "hg38")))


### Colors ###
dark2_colors <- brewer.pal(n = 7, "Dark2")
custom_colors <- c(
  "low confidence" = "black",
  "2h" = dark2_colors[1],
  "CHX" = dark2_colors[2],
  "4h DMSO" = dark2_colors[3],
  "4h STM2457" = dark2_colors[4],
  "6h" = dark2_colors[5],
  "8h" = dark2_colors[6]
)

custom_colors_hist <- c(
  "Run 3 WDX-4 (2h)" = dark2_colors[1],
  "Run 3 WDX-5 (CHX)" = dark2_colors[2],
  "Run 4 WDX-4 (4h DMSO)" = dark2_colors[3],
  "Run 4 WDX-5 (4h STM2457)" = dark2_colors[4],
  "Run 4 WDX-6 (6h)" = dark2_colors[5],
  "Run 4 WDX-11 (8h)" = dark2_colors[6],
  "Run 3 unassigned" = "black",
  "Run 4 unassigned" = "black"
)
purples <- brewer.pal(5, "Purples")
greens <- brewer.pal(5, "Greens")

alignment_colors <- c(
  "unaligned" = "grey50",
  "HSV-2" = purples[4],
  "hg38" = greens[4]
)
alignment_colors_eU <- c(
  "unaligned" = "grey50",
  "HSV-2" = purples[2],
  "hg38" = greens[2]
)
alignment_colors_combined <- c(
  "Unaligned total" = "grey50",
  "Unaligned eU" = "grey50",
  "HSV-2 total" = purples[4],
  "hg38 total" = greens[4],
  "HSV-2 eU" = purples[2],
  "hg38 eU" = greens[2]
)

eu_colors <- c(
  "eU" = brewer.pal(5, "Blues")[3],
  "non-eU" = brewer.pal(5, "Oranges")[3]
)


plot_list <- list()
max_count <- 0

for (cond_var in cond_values) {
  # subset only rows for this condition
  subset_data <- subset(df_hist, cond == cond_var)
  if (nrow(subset_data) == 0) next
  
  # compute counts for this condition
  counts <- hist(subset_data$confidence_score_pct, breaks = bin_breaks, plot = FALSE)$counts
  max_count <- max(max_count, max(counts))
  
  # create the plot
  p <- ggplot(subset_data, aes(x = confidence_score_pct)) +
    geom_histogram(breaks = bin_breaks,
                   fill = custom_colors_hist[as.character(cond_var)],
                   color = "black") +
    geom_vline(xintercept = c(90, 98),
               linetype = "dashed", color = "grey30", linewidth = 0.6) +
    geom_vline(xintercept = 95,
               linetype = "dashed", color = "red", linewidth = 0.6) +
    labs(title = cond_var,
         x = "Confidence Score (%)",
         y = "Frequency") +
    xlim(0, 100) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid.major.y = element_blank(),   # remove major horizontal grid lines
      panel.grid.minor = element_blank(),     # remove minor grid lines
      panel.grid.major.x = element_blank(),   # remove vertical grid lines
      axis.text = element_text(color = "black"),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0)  # bold header and centered
    )
  
  plot_list[[as.character(cond_var)]] <- p
}

# decide what the alternate max should be
alt_max_count <- max_count * 0.02   # adjust as needed

# get plot order
n_plots <- length(plot_list)

# apply different y-limits
for (i in seq_along(plot_list)) {
  
  if (i > n_plots - 2) {
    # last two histograms
    plot_list[[i]] <- plot_list[[i]] +
      scale_y_continuous(
        limits = c(0, alt_max_count),
        labels = scales::scientific
      )
  } else {
    # all other histograms
    plot_list[[i]] <- plot_list[[i]] +
      scale_y_continuous(
        limits = c(0, max_count),
        labels = scales::scientific
      )
  }
}


# combine into a grid
combined_histograms <- wrap_plots(plot_list, ncol = 2)

# save
ggsave("confidence_histograms.png", combined_histograms,
       width = 16, height = 12, dpi = 300)



# Plot with manual colors
p <- ggplot(filter(df_all_final, !is.na(cond)), aes(x = Experiment, fill = cond)) +
  geom_bar(position = "stack", color = "black", width = 0.7) +
  scale_fill_manual(values = custom_colors) +  # <-- custom colors
  labs(
    x = "Experiment",
    y = "Number of reads",
    fill = "Condition"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.position = "right",
    legend.key.size = unit(0.6, "cm")
  )


p_alignment <- ggplot(df_alignment_plot, aes(x = cond, fill = alignment)) +
  geom_bar(position = "stack", color = "black", width = 0.7) +
  scale_fill_manual(values = alignment_colors) +
  labs(
    x = "Condition",
    y = "Number of reads",
    fill = "Alignment"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # <-- rotate labels
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.position = "right",
    legend.key.size = unit(0.6, "cm")
  )


p_eu_alignment <- ggplot(filter(df_eu_plot, cond != "CHX"), aes(x = cond, fill = alignment)) +
  geom_bar(position = "stack", color = "black", width = 0.7) +
  scale_fill_manual(values = alignment_colors_eU) +
  labs(
    x = "Condition",
    y = "Number of eU-labelled reads",
    fill = "Alignment"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),  # diagonal labels
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.position = "right",
    legend.key.size = unit(0.6, "cm")
  )

p_eu_status <- ggplot(filter(df_all_final, cond != "low confidence"), aes(x = cond, fill = eU_label)) +
  geom_bar(position = "stack", color = "black", width = 0.7) +
  scale_fill_manual(values = eu_colors) +
  labs(
    x = "Condition",
    y = "Number of reads",
    fill = "eU-label"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),  # diagonal labels
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.position = "right",
    legend.key.size = unit(0.6, "cm")
  )


#plot normalised

# Plot with manual colors
p1 <- ggplot(filter(df_all_final, !is.na(cond)), aes(x = Experiment, fill = cond)) +
  geom_bar(position = "fill", color = "black", width = 0.7) +
  scale_fill_manual(values = custom_colors) +  # <-- custom colors
  labs(
    x = "Experiment",
    y = "Number of reads (normalised)",
    fill = "Condition"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.position = "right",
    legend.key.size = unit(0.6, "cm")
  )


p1_alignment <- ggplot(df_alignment_plot, aes(x = cond, fill = alignment)) +
  geom_bar(position = "fill", color = "black", width = 0.7) +
  scale_fill_manual(values = alignment_colors) +
  labs(
    x = "Condition",
    y = "Number of reads (normalised)",
    fill = "Alignment"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # <-- rotate labels
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.position = "right",
    legend.key.size = unit(0.6, "cm")
  )


p1_eu_alignment <- ggplot(filter(df_eu_plot, cond != "CHX"), aes(x = cond, fill = alignment)) +
  geom_bar(position = "fill", color = "black", width = 0.7) +
  scale_fill_manual(values = alignment_colors_eU) +
  labs(
    x = "Condition",
    y = "Number of eU-labelled reads (normalised)",
    fill = "Alignment"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),  # diagonal labels
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.position = "right",
    legend.key.size = unit(0.6, "cm")
  )

p1_eu_status <- ggplot(filter(df_all_final, cond != "low confidence"), aes(x = cond, fill = eU_label)) +
  geom_bar(position = "fill", color = "black", width = 0.7) +
  scale_fill_manual(values = eu_colors) +
  labs(
    x = "Condition",
    y = "Number of reads (normalised)",
    fill = "eU-label"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),  # diagonal labels
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.position = "right",
    legend.key.size = unit(0.6, "cm")
  )


# Save at high resolution
ggsave("stacked_bar_publication.png", plot = p, width = 6, height = 4, dpi = 300)
ggsave("stacked_bar_publication.pdf", plot = p, width = 6, height = 4)
ggsave("stacked_bar_alignment.png", plot = p_alignment, width = 6, height = 4, dpi = 300)
ggsave("stacked_bar_alignment.pdf", plot = p_alignment, width = 6, height = 4)
ggsave("stacked_bar_alignment_eU.png", plot = p_eu_alignment, width = 6, height = 4, dpi = 300)
ggsave("stacked_bar_alignment_eU.pdf", plot = p_eu_alignment, width = 6, height = 4)
ggsave("stacked_bar_publication_proportional.png", plot = p1, width = 6, height = 4, dpi = 300)
ggsave("stacked_bar_publication_proportional.pdf", plot = p1, width = 6, height = 4)
ggsave("stacked_bar_alignment_proportional.png", plot = p1_alignment, width = 6, height = 4, dpi = 300)
ggsave("stacked_bar_alignment_proportional.pdf", plot = p1_alignment, width = 6, height = 4)
ggsave("stacked_bar_alignment_eU_proportional.png", plot = p1_eu_alignment, width = 6, height = 4, dpi = 300)
ggsave("stacked_bar_alignment_eU_proportional.pdf", plot = p1_eu_alignment, width = 6, height = 4)
ggsave("stacked_bar_eU_status.png", plot = p_eu_status, width = 6, height = 4, dpi = 300)
ggsave("stacked_bar_eU_status.pdf", plot = p_eu_status, width = 6, height = 4)
ggsave("stacked_bar_eU_status_proportional.png", plot = p1_eu_status, width = 6, height = 4, dpi = 300)
ggsave("stacked_bar_eU_status_proportional.pdf", plot = p1_eu_status, width = 6, height = 4)




###combine condition plots###
combined_stacked_barplots <- (
  (p_alignment | p_eu_status | p_eu_alignment ) /
    (p1_alignment | p1_eu_status | p1_eu_alignment )
) +
  plot_annotation(
    title = "Read alignment and eU-label distribution across conditions",
    subtitle = "Top: absolute read counts | Bottom: proportional composition",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5)
    )
  )

# Save the combined grid as PNG and PDF
ggsave("combined_stacked_barplots_reordered.png", combined_stacked_barplots,
       width = 18, height = 8, dpi = 300)
ggsave("combined_stacked_barplots_reordered.pdf", combined_stacked_barplots,
       width = 18, height = 8)


###combine condition plots normalised only###
combined_stacked_barplots_norm <- (
    (p1_alignment | p1_eu_status | p1_eu_alignment )
) +
  plot_annotation(
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5)
    )
  )

# Save the combined grid as PNG and PDF
ggsave("combined_stacked_barplots_reordered_norm.png", combined_stacked_barplots_norm,
       width = 18, height = 6, dpi = 300)
ggsave("combined_stacked_barplots_reordered_norm.pdf", combined_stacked_barplots_norm,
       width = 18, height = 6)
