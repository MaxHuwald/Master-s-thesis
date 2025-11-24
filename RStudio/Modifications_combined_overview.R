setwd("/project/sysviro/users/Max/analyses/Conditions/Results/Modifications/thesis")
library(dplyr)
library(readr)
library(stringr)
library(purrr)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(tidyr)
# Base directory
base_dir <- "/project/sysviro/users/Max/analyses/Modkit"
holden_file <- "/project/sysviro/HOLDEN/HG38/histograms/ARPE19_polyA_TSO_IVT.sup-allMods.trimAdapters.dorado.0.9.0_hg38_probabilities.tsv"
# Get all *_probabilities.tsv files
all_files <- list.files(base_dir, pattern = "_probabilities.tsv$", 
                        recursive = TRUE, full.names = TRUE)

# Function to parse filenames into metadata
parse_filename <- function(filepath) {
  folder <- basename(dirname(filepath))   # eU, non_eU, or Total
  fname <- basename(filepath)
  
  if (folder %in% c("eU", "non_eU")) {
    # Example: 4h.hg38.eU.0.9.filtered_probabilities.tsv
    condition <- str_extract(fname, "^[^.]+")  # everything before first dot
    alignment_label <- str_match(fname, "^[^.]+\\.([^.]+)")[,2]  # between first and second dot
    
  } else if (folder == "Total") {
    # Example 1: 2h.trimAdapters.dorado.1.1.1.filtered.aligned_hg38.primary.fwd_probabilities.tsv
    if (str_detect(fname, "aligned_hg38")) {
      condition <- str_extract(fname, "^[^.]+")  # everything before first dot
      alignment_label <- "hg38"
      
      # Example 2: 2h.trimAdapters.dorado.1.1.1.untrimmed.filtered.aligned.primary_probabilities.tsv
    } else if (str_detect(fname, "untrimmed")) {
      condition <- str_extract(fname, "^[^.]+")  # everything before first dot
      alignment_label <- "untrimmed"
    } else {
      condition <- NA
      alignment_label <- NA
    }
  } else {
    condition <- NA
    alignment_label <- NA
  }
  
  tibble(filepath = filepath,
         folder = folder,
         condition = condition,
         alignment_label = alignment_label)
}

# Parse all filenames
metadata <- map_dfr(all_files, parse_filename)

# Read and combine data
all_data <- metadata %>%
  mutate(data = map(filepath, ~ read_tsv(.x, show_col_types = FALSE))) %>%
  unnest(data)


### holden (IVT) data
holden_data <- read_tsv(holden_file, show_col_types = FALSE)

# Step 2 — Define the folder/alignment combinations
combos <- expand.grid(
  folder = c("Total", "eU"),
  alignment_label = c("hg38", "untrimmed"),
  stringsAsFactors = FALSE
)

# Step 3 — Duplicate data for each combo
holden_data_all <- bind_rows(
  lapply(seq_len(nrow(combos)), function(i) {
    holden_data %>%
      mutate(
        folder = combos$folder[i],
        alignment_label = combos$alignment_label[i],
        condition = "IVT"
      )
  })
)

all_data <- bind_rows(all_data, holden_data_all)

# Recode modification codes to human-readable names
all_data <- all_data %>%
  mutate(
    code = case_when(
      code == "a"       ~ "m6A",
      code == "17596"   ~ "Inosine",
      code == "69426"   ~ "2OmeA",
      code == "17802"   ~ "pseU",
      code == "19227"   ~ "2OmeU",
      code == "-" & primary_base == "A" ~ "A",
      code == "-" & primary_base == "T" ~ "U",
      TRUE ~ code
    )
  )
all_data <- all_data %>%
  mutate(condition = recode(condition,
                            "4h_DMSO" = "4h DMSO",
                            "4h_STM2457" = "4h STM2457"))


all_data <- all_data %>%
  mutate(alignment_label = recode(alignment_label, "untrimmed" = "HSV-2"))

all_data <- all_data %>%
  mutate(folder = recode(folder, "non_eU" = "non-eU"))
# Add probability bin center
all_data <- all_data %>%
  mutate(bin_center = (range_start + range_end) / 2)

# optional: filter out 4h cond.

#all_data <- all_data %>%
# filter(condition != "4h")

# order levels

condition_levels <- c("CHX", "2h", "4h DMSO", "4h STM2457", "4h", "6h", "8h", "IVT")
all_data$condition <- factor(all_data$condition, levels = condition_levels)
condition_levels_eU <- c("IVT", "2h", "4h DMSO", "4h STM2457", "6h", "8h")

# Colorscheminng (not implemented)

all_data <- all_data %>%
  mutate(
    code_pm = dplyr::case_when(
      code == "m6A" ~ "m^6*A",   # superscript 6
      TRUE          ~ code       # default: leave as-is
    )
  )


plot_mod_distribution_hist <- function(df, folder, alignment, cond, code,
                                       max_frac_global = NULL) {
  df_cond <- df %>%
    filter(folder == !!folder,
           alignment_label == !!alignment,
           condition == !!cond,
           code == !!code)
  
  if (nrow(df_cond) == 0) return(NULL)
  
  if (is.null(max_frac_global)) {
    df_max_frac <- df %>%
      filter(folder == !!folder,
             alignment_label == !!alignment,
             code == !!code)
    max_frac <- max(df_max_frac$frac, na.rm = TRUE)
  } else {
    max_frac <- max_frac_global
  }
  
  bin_width <- min(diff(sort(unique(df_cond$bin_center))))
  if (!is.finite(bin_width) || bin_width <= 0) bin_width <- 0.00390625
  
  # ---- NEW: special title for IVT ----
  if (as.character(cond) == "IVT") {
    title_expr <- "IVT"
  } else {
    code_pm_str  <- df_cond$code_pm[1]
    code_pm_expr <- parse(text = code_pm_str)[[1]]
    title_expr   <- bquote(.(folder) ~ .(alignment) ~ .(cond) ~ .(code_pm_expr))
  }
  # -----------------------------------
  
  ggplot(df_cond, aes(x = bin_center, y = frac)) +
    geom_col(width = bin_width, position = "dodge", color = "blue", fill = "blue") +
    geom_vline(xintercept = 0.98, linetype = "dotted", color = "black", size = 0.6) +  # vertical line at 0.98
    scale_fill_manual(values = "blue") +
    scale_y_continuous(labels = scales::scientific, limits = c(0, max_frac)) +
    labs(
      title = title_expr,
      x = "Modification probability",
      y = "Fraction"
    ) +
    coord_cartesian(xlim = c(0.95, 1.0)) +   # if you want the x-zoom
    theme_classic() +
    theme(
      plot.title  = element_text(size = 9, face = "bold", color = "black"),
      axis.title.x = element_text(size = 10, color = "black"),
      axis.title.y = element_text(size = 10, color = "black"),
      axis.text.x  = element_text(size = 10, color = "black"),
      axis.text.y  = element_text(size = 10, color = "black"),
      legend.position = "none"
    )
}



mod_code <- "pseU"  # or "pseU", etc.

# max fraction shared between eU + non-eU for EACH alignment
max_frac_hg38 <- all_data %>%
  filter(alignment_label == "hg38",
         folder %in% c("eU", "non-eU"),
         code == mod_code) %>%
  summarise(max_frac = max(frac, na.rm = TRUE)) %>%
  pull(max_frac)

max_frac_hsv2 <- all_data %>%
  filter(alignment_label == "HSV-2",
         folder %in% c("eU", "non-eU"),
         code == mod_code) %>%
  summarise(max_frac = max(frac, na.rm = TRUE)) %>%
  pull(max_frac)

## ---------- hg38 ----------

# non-eU hg38
plots_hg38_non_eU <- purrr::map(
  condition_levels[!condition_levels %in% c("4h", "IVT")],
  ~ plot_mod_distribution_hist(
    all_data %>% filter(condition != "4h"),
    folder   = "non-eU",
    alignment = "hg38",
    cond     = .x,
    code     = mod_code,
    max_frac_global = max_frac_hg38
  )
)

# eU hg38
plots_hg38_eU <- purrr::map(
  condition_levels_eU[!condition_levels_eU %in% c("4h", "CHX")],
  ~ plot_mod_distribution_hist(
    all_data %>%
      filter(!condition %in% c("4h", "CHX")) %>%
      mutate(condition = factor(condition, levels = condition_levels_eU)),
    folder   = "eU",
    alignment = "hg38",
    cond     = .x,
    code     = mod_code,
    max_frac_global = max_frac_hg38
  )
)

## ---------- HSV-2 ----------

# non-eU HSV-2
plots_hsv2_non_eU <- purrr::map(
  condition_levels[!condition_levels %in% c("4h", "IVT")],
  ~ plot_mod_distribution_hist(
    all_data %>% filter(condition != "4h"),
    folder   = "non-eU",
    alignment = "HSV-2",
    cond     = .x,
    code     = mod_code,
    max_frac_global = max_frac_hsv2
  )
)

# eU HSV-2
plots_hsv2_eU <- purrr::map(
  condition_levels_eU[!condition_levels_eU %in% c("4h", "CHX")],
  ~ plot_mod_distribution_hist(
    all_data %>%
      filter(!condition %in% c("4h", "CHX")) %>%
      mutate(condition = factor(condition, levels = condition_levels_eU)),
    folder   = "eU",
    alignment = "HSV-2",
    cond     = .x,
    code     = mod_code,
    max_frac_global = max_frac_hsv2
  )
)

# wrap into grids (same as before, or tweak ncol)
combined_hg38_non_eU <- wrap_plots(plots_hg38_non_eU, ncol = 2)
combined_hg38_eU     <- wrap_plots(plots_hg38_eU,     ncol = 2)
combined_hsv2_non_eU <- wrap_plots(plots_hsv2_non_eU, ncol = 2)
combined_hsv2_eU     <- wrap_plots(plots_hsv2_eU,     ncol = 2)


# wrap into grids (same as before, or tweak ncol)
combined_hg38_non_eU <- wrap_plots(plots_hg38_non_eU, ncol = 2)
combined_hg38_eU     <- wrap_plots(plots_hg38_eU,     ncol = 2)
combined_hsv2_non_eU <- wrap_plots(plots_hsv2_non_eU, ncol = 2)
combined_hsv2_eU     <- wrap_plots(plots_hsv2_eU,     ncol = 2)

big_plot <- (combined_hg38_non_eU | combined_hg38_eU) /
  (combined_hsv2_non_eU | combined_hsv2_eU)

# optional overall title
#big_plot <- big_plot + plot_annotation(
#  title = paste0("Modification probability distributions – ", mod_code)
#)

ggsave(paste0("combined_eU_non_eU_hg38_HSV2_", mod_code, ".png"),
       plot   = big_plot,
       width  = 9.5,   # adjust as you like
       height = 10,
       dpi    = 300,
       units  = "in")
ggsave(paste0("combined_hsv2_eU_", mod_code, ".png"), 
       plot = combined_hsv2_eU, width = 7, height = 8, dpi = 300, device = "png", units = "in")
