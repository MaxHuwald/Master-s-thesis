# ----------------------------
# Setup
# ----------------------------
setwd("/project/sysviro/users/Max/analyses/Conditions/Results/Coverageplots/normalised_2h")

library(Gviz)
library(tidyverse)
library(stringr)
library(dplyr)
library(RColorBrewer)
library(data.table)
library(grid)
library(txdbmaker)

# Genome/chromosome
myChr <- "HSV2-st333"
myStart <- 1
myEnd <- 155503
options(ucscChromosomeNames = FALSE)

# Path to Visualisation directory
base_dir <- "/project/sysviro/users/Max/analyses/Conditions/subsampled/Visualisation"

# ----------------------------
# List bedgraph files
# ----------------------------
bed_files <- list.files(base_dir, pattern = "\\.bedgraph$", full.names = TRUE)

# Gene models
modelsPos <- makeTxDbFromGFF("/project/sysviro/users/Max/Genome/HSV2-st333_Annotation_v4.forward.gff3")
modelsNeg <- makeTxDbFromGFF("/project/sysviro/users/Max/Genome/HSV2-st333_Annotation_v4.reverse.gff3")
rtrackFor <- GeneRegionTrack(modelsPos, genome = myChr, chromosome = myChr,
                             name = "Gene Model", col="black", fill="grey", stacking="squish",
                             shape="smallArrow", background.title = "transparent", showFeatureId = TRUE,
                             geneSymbol = TRUE, transcriptAnnotation = "symbol",just.group = "above", fontcolor.group = "black", cex.group = 1.04)
rtrackRev <- GeneRegionTrack(modelsNeg, genome = myChr, chromosome = myChr,
                             name = "Gene Model", col="black", fill="grey", stacking="squish",
                             shape="smallArrow", background.title = "transparent", showFeatureId = TRUE,
                             geneSymbol = TRUE, transcriptAnnotation = "symbol",just.group = "above", fontcolor.group = "black", cex.group = 1.04)
gtrack <- GenomeAxisTrack(col="black", cex = 1.3, fontcolor = "black") # genome axis track

# ----------------------------
# Parse metadata
# ----------------------------
file_info <- tibble(file = bed_files) %>%
  mutate(filename = basename(file)) %>%
  mutate(condition = str_extract(filename, "^[^\\.]+")) %>%
  mutate(condition = case_when(
    condition == "4h_DMSO"    ~ "4h DMSO",
    condition == "4h_STM2457" ~ "4h STM2457",
    TRUE ~ condition
  )) %>%
  mutate(
    eU_status = case_when(
      str_detect(filename, "\\.eU\\.") ~ "eU",
      str_detect(filename, "non_eU")   ~ "non-eU",
      TRUE                             ~ "total"
    ),
    direction = case_when(
      str_detect(filename, "\\.fwd\\.") ~ "fwd",
      str_detect(filename, "\\.rev\\.") ~ "rev",
      TRUE                              ~ NA_character_
    )
  )

# Remove erroneous CHX eU rows
file_info <- file_info %>%
  filter(!(condition == "CHX" & eU_status == "eU"))

# ----------------------------
# Color maps
# ----------------------------
color_map <- c(
  "eU"     = brewer.pal(5, "Purples")[2],
  "non-eU" = brewer.pal(5, "Purples")[5]
)

combined_color_map <- c(
  "4h DMSO non-eU"     = brewer.pal(5, "Blues")[5],
  "4h DMSO eU"         = brewer.pal(5, "Blues")[2],
  "4h STM2457 non-eU"  = brewer.pal(5, "Oranges")[5],
  "4h STM2457 eU"      = brewer.pal(5, "Oranges")[2]
)
# Replace underscores for consistency
names(combined_color_map) <- str_replace_all(names(combined_color_map), "_", " ")

combined_color_map2 <- c(
  "CHX eU"     = brewer.pal(5, "Blues")[2],
  "2h eU"      = brewer.pal(5, "Oranges")[2],
  "CHX non eU" = brewer.pal(5, "Blues")[5],
  "2h non eU"  = brewer.pal(5, "Oranges")[5]
)

# ----------------------------
# Add max_value column
# ----------------------------
file_info <- file_info %>%
  rowwise() %>%
  mutate(max_value = max(fread(file, col.names = c("chromosome","start","end","value"))$value, na.rm = TRUE)) %>%
  ungroup()

# ----------------------------
# Compute global normalisation factors from "total" tracks
# ----------------------------
total_max_fwd <- file_info %>%
  filter(eU_status == "total", direction == "fwd") %>%
  pull(max_value) %>% max(na.rm = TRUE)

total_max_rev <- file_info %>%
  filter(eU_status == "total", direction == "rev") %>%
  pull(max_value) %>% max(na.rm = TRUE)

message("Normalization factors - FWD: ", total_max_fwd, ", REV: ", total_max_rev)

# ----------------------------
# Create DataTracks (normalised)
# ----------------------------
dataTracks <- list()

for (i in seq_len(nrow(file_info))) {
  row <- file_info[i, ]
  track_name <- paste(row$condition, row$eU_status, row$direction, sep = " ")
  
  df <- fread(row$file, col.names = c("chromosome","start","end","value"))
  norm_factor <- if(row$direction == "fwd") total_max_fwd else total_max_rev
  df$value <- df$value / norm_factor
  
  dt <- DataTrack(
    range = df, type = "l",
    chromosome = myChr, genome = myChr,
    fill = color_map[row$eU_status], col = color_map[row$eU_status],
    col.axis = "black",
    background.title = "transparent",
    background.panel = "transparent",
    ylim = if(row$direction == "rev") c(1, 0) else c(0, 1)
  )
  
  dataTracks[[track_name]] <- dt
}

# ----------------------------
# Create DataTracks (RAW, no normalisation)
# ----------------------------
ylim_info_raw <- file_info %>%
  filter(eU_status != "total") %>%
  group_by(condition, direction) %>%
  summarise(shared_max = max(max_value, na.rm = TRUE), .groups = "drop")

dataTracks_raw <- list()

for (i in seq_len(nrow(file_info))) {
  row <- file_info[i, ]
  track_name <- paste(row$condition, row$eU_status, row$direction, sep = " ")
  
  df <- fread(row$file, col.names = c("chromosome","start","end","value"))
  
  if(row$eU_status != "total") {
    shared_max <- ylim_info_raw %>%
      filter(condition == row$condition, direction == row$direction) %>%
      pull(shared_max) %>% max(na.rm = TRUE)
  } else {
    shared_max <- row$max_value
  }
  
  dt <- DataTrack(
    range = df, type = "l",
    chromosome = myChr, genome = myChr,
    fill = color_map[row$eU_status], col = color_map[row$eU_status],
    col.axis = "black",
    background.title = "transparent",
    background.panel = "transparent",
    ylim = if(row$direction == "rev") c(shared_max, 0) else c(0, shared_max)
  )
  
  dataTracks_raw[[track_name]] <- dt
}

# ----------------------------
# Mega plot function with legend
# ----------------------------
mega_plot <- function(dataTracks, file_info, cond_order, from, to, chr, color_map) {
  
  # Forward overlays
  fwd_overlays <- lapply(cond_order, function(cond) {
    cond_tracks <- file_info %>%
      filter(condition == cond, eU_status != "total", direction == "fwd") %>%
      mutate(name = paste(condition, eU_status, direction, sep = " ")) %>%
      pull(name)
    # skip any missing tracks
    cond_tracks <- cond_tracks[cond_tracks %in% names(dataTracks)]
    
    OverlayTrack(
      trackList = lapply(cond_tracks, function(n) dataTracks[[n]]),
      name = cond,
      col.axis = "black",
      col.title = "black",
      background.title = "transparent",
      col.frame = "transparent"
    )
  })
  
  # Reverse overlays
  rev_overlays <- lapply(rev(cond_order), function(cond) {
    cond_tracks <- file_info %>%
      filter(condition == cond, eU_status != "total", direction == "rev") %>%
      mutate(name = paste(condition, eU_status, direction, sep = " ")) %>%
      pull(name)
    cond_tracks <- cond_tracks[cond_tracks %in% names(dataTracks)]
    
    OverlayTrack(
      trackList = lapply(cond_tracks, function(n) dataTracks[[n]]),
      name = cond,
      col.axis = "black",
      col.title = "black",
      background.title = "transparent",
      col.frame = "transparent"
    )
  })
  
  # Combine list
  all_tracks <- c(fwd_overlays, list(rtrackFor, gtrack, rtrackRev), rev_overlays)
  
  # Plot
  plotTracks(
    all_tracks,
    from = from, to = to, chromosome = chr,
    type = "l",
    cex.title = 1,
    cex.axis = 0.7,
    title.width = 1.5,
    track.padding = 0.1
  )
  
  # Legend
  legend_colors <- color_map[names(color_map)]
  labels <- names(legend_colors)
  legend_x <- 0.01
  legend_y <- 0.56
  line_spacing <- 0.03
  
  for (i in seq_along(labels)) {
    y_pos <- legend_y - 0.04 - (i - 1) * line_spacing
    grid.lines(x = unit(c(legend_x, legend_x + 0.03), "npc"),
               y = unit(rep(y_pos, 2), "npc"),
               gp = gpar(col = legend_colors[i], lwd = 4))
    grid.text(labels[i],
              x = unit(legend_x + 0.04, "npc"),
              y = unit(y_pos, "npc"),
              just = "left", gp = gpar(fontsize = 11))
  }
}







###Megaplot for combined conds:
mega_plot_combined <- function(file_info, cond_order, from, to, chr, color_map) {
  
  dataTracks <- list()
  
  # Create DataTracks inside the function
  for (i in seq_len(nrow(file_info))) {
    row <- file_info[i, ]
    
    # Skip tracks not in combined color map or 'total'
    color_key <- paste(row$condition, row$eU_status)
    color_key <- str_replace_all(color_key, "_", " ")
    if (!(color_key %in% names(color_map)) || row$eU_status == "total") next
    
    track_name <- paste(row$condition, row$eU_status, row$direction, sep = " ")
    
    df <- fread(row$file, col.names = c("chromosome","start","end","value"))
    norm_factor <- if(row$direction == "fwd") total_max_fwd else total_max_rev
    df$value <- df$value / norm_factor
    
    dt <- DataTrack(
      range = df, type = "l",
      chromosome = chr, genome = chr,
      fill = color_map[color_key], col = color_map[color_key],
      col.axis = "black",
      background.title = "transparent",
      background.panel = "transparent",
      ylim = if(row$direction == "rev") c(1, 0) else c(0, 1)
    )
    
    dataTracks[[track_name]] <- dt
  }
  
  # Function to create overlays safely
  create_overlay <- function(cond, direction = "fwd") {
    cond_tracks <- names(dataTracks)[str_detect(names(dataTracks), paste0("^", cond, ".*", direction))]
    if (length(cond_tracks) == 0) return(NULL)
    if (length(cond_tracks) == 1) return(dataTracks[[cond_tracks]])
    OverlayTrack(
      trackList = lapply(cond_tracks, function(n) dataTracks[[n]]),
      name = cond,
      col.axis = "black",
      col.title = "black",
      background.title = "transparent",
      col.frame = "transparent"
    )
  }
  
  fwd_overlays <- lapply(cond_order, create_overlay)
  rev_overlays <- lapply(rev(cond_order), function(cond) create_overlay(cond, direction = "rev"))
  
  # Remove NULLs
  fwd_overlays <- fwd_overlays[!sapply(fwd_overlays, is.null)]
  rev_overlays <- rev_overlays[!sapply(rev_overlays, is.null)]
  
  # Combine tracks
  all_tracks <- c(fwd_overlays, list(rtrackFor, gtrack, rtrackRev), rev_overlays)
  
  # Plot
  plotTracks(
    all_tracks,
    from = from, to = to, chromosome = chr,
    type = "l",
    cex.title = 1,
    cex.axis = 0.7,
    title.width = 1.5,
    track.padding = 0.1
  )
  
  # Legend
  legend_colors <- color_map[names(color_map)]
  labels <- names(legend_colors)
  legend_x <- 0.85
  legend_y <- 1
  line_spacing <- 0.03
  
  for (i in seq_along(labels)) {
    y_pos <- legend_y - 0.04 - (i - 1) * line_spacing
    grid.lines(x = unit(c(legend_x, legend_x + 0.03), "npc"),
               y = unit(rep(y_pos, 2), "npc"),
               gp = gpar(col = legend_colors[i], lwd = 2))
    grid.text(labels[i],
              x = unit(legend_x + 0.04, "npc"),
              y = unit(y_pos, "npc"),
              just = "left", gp = gpar(fontsize = 9))
  }
}

###megaplot for non_normalised combined
mega_plot_combined_raw <- function(file_info, cond_order, from, to, chr, color_map) {
  
  dataTracks <- list()
  
  # Determine shared max per condition and direction (excluding 'total')
  ylim_info <- file_info %>%
    filter(eU_status != "total") %>%
    group_by(condition, direction) %>%
    summarise(shared_max = max(max_value, na.rm = TRUE), .groups = "drop")
  
  # Create DataTracks using raw values
  for (i in seq_len(nrow(file_info))) {
    row <- file_info[i, ]
    
    color_key <- paste(row$condition, row$eU_status)
    color_key <- str_replace_all(color_key, "_", " ")
    if (!(color_key %in% names(color_map)) || row$eU_status == "total") next
    
    track_name <- paste(row$condition, row$eU_status, row$direction, sep = " ")
    
    df <- fread(row$file, col.names = c("chromosome","start","end","value"))
    
    # Determine y-limit for this track
    if(row$eU_status != "total") {
      shared_max <- ylim_info %>%
        filter(condition == row$condition, direction == row$direction) %>%
        pull(shared_max) %>% max(na.rm = TRUE)
    } else {
      shared_max <- row$max_value
    }
    
    dt <- DataTrack(
      range = df, type = "l",
      chromosome = chr, genome = chr,
      fill = color_map[color_key], col = color_map[color_key],
      col.axis = "black",
      background.title = "transparent",
      background.panel = "transparent",
      ylim = if(row$direction == "rev") c(shared_max, 0) else c(0, shared_max)
    )
    
    dataTracks[[track_name]] <- dt
  }
  
  # Helper function for overlays
  create_overlay <- function(cond, direction = "fwd") {
    cond_tracks <- names(dataTracks)[str_detect(names(dataTracks), paste0("^", cond, ".*", direction))]
    if (length(cond_tracks) == 0) return(NULL)
    if (length(cond_tracks) == 1) return(dataTracks[[cond_tracks]])
    OverlayTrack(
      trackList = lapply(cond_tracks, function(n) dataTracks[[n]]),
      name = cond,
      col.axis = "black",
      col.title = "black",
      background.title = "transparent",
      col.frame = "transparent"
    )
  }
  
  fwd_overlays <- lapply(cond_order, create_overlay)
  rev_overlays <- lapply(rev(cond_order), function(cond) create_overlay(cond, direction = "rev"))
  
  # Remove NULLs
  fwd_overlays <- fwd_overlays[!sapply(fwd_overlays, is.null)]
  rev_overlays <- rev_overlays[!sapply(rev_overlays, is.null)]
  
  # Combine tracks
  all_tracks <- c(fwd_overlays, list(rtrackFor, gtrack, rtrackRev), rev_overlays)
  
  # Plot
  plotTracks(
    all_tracks,
    from = from, to = to, chromosome = chr,
    type = "l",
    cex.title = 1,
    cex.axis = 0.7,
    title.width = 1.5,
    track.padding = 0.1
  )
  
  # Legend
  legend_colors <- color_map[names(color_map)]
  labels <- names(legend_colors)
  legend_x <- 0.2
  legend_y <- 0.7
  line_spacing <- 0.03
  
  for (i in seq_along(labels)) {
    y_pos <- legend_y - 0.04 - (i - 1) * line_spacing
    grid.lines(x = unit(c(legend_x, legend_x + 0.03), "npc"),
               y = unit(rep(y_pos, 2), "npc"),
               gp = gpar(col = legend_colors[i], lwd = 5))
    grid.text(labels[i],
              x = unit(legend_x + 0.04, "npc"),
              y = unit(y_pos, "npc"),
              just = "left", gp = gpar(fontsize = 14))
  }
}


mega_plot_combined_normalised_overlay <- function(file_info, cond_order, from, to, chr, color_map) {
  
  dataTracks <- list()
  
  # Create DataTracks inside the function
  for (i in seq_len(nrow(file_info))) {
    row <- file_info[i, ]
    
    # Skip tracks not in combined color map or 'total'
    color_key <- paste(row$condition, row$eU_status)
    color_key <- str_replace_all(color_key, "_", " ")
    if (!(color_key %in% names(color_map)) || row$eU_status == "total") next
    
    track_name <- paste(row$condition, row$eU_status, row$direction, sep = " ")
    
    df <- fread(row$file, col.names = c("chromosome","start","end","value"))
    norm_factor <- if(row$direction == "fwd") total_max_fwd else total_max_rev
    df$value <- df$value / norm_factor
    
    dt <- DataTrack(
      range = df, type = "l",
      chromosome = chr, genome = chr,
      fill = color_map[color_key], col = color_map[color_key],
      col.axis = "black",
      background.title = "transparent",
      background.panel = "transparent",
      ylim = if(row$direction == "rev") c(1, 0) else c(0, 1)
    )
    
    dataTracks[[track_name]] <- dt
  }
  
  # Combine all forward tracks into a single overlay
  fwd_tracks <- names(dataTracks)[str_detect(names(dataTracks), "fwd$")]
  forward_overlay <- OverlayTrack(
    trackList = lapply(fwd_tracks, function(n) dataTracks[[n]]),
    name = "Forward (all conditions)",
    col.axis = "black",
    col.title = "black",
    background.title = "transparent",
    col.frame = "transparent"
  )
  
  # Combine all reverse tracks into a single overlay
  rev_tracks <- names(dataTracks)[str_detect(names(dataTracks), "rev$")]
  reverse_overlay <- OverlayTrack(
    trackList = lapply(rev_tracks, function(n) dataTracks[[n]]),
    name = "Reverse (all conditions)",
    col.axis = "black",
    col.title = "black",
    background.title = "transparent",
    col.frame = "transparent"
  )
  
  # Combine tracks with gene models and genome axis
  all_tracks <- c(list(forward_overlay, rtrackFor, gtrack, rtrackRev, reverse_overlay))
  
  # Plot
  plotTracks(
    all_tracks,
    from = from, to = to, chromosome = chr,
    type = "l",
    cex.title = 1,
    cex.axis = 0.7,
    title.width = 1.5,
    track.padding = 0.1
  )
  
  # Legend
  legend_colors <- color_map[names(color_map)]
  labels <- names(legend_colors)
  legend_x <- 0.85
  legend_y <- 1
  line_spacing <- 0.03
  
  for (i in seq_along(labels)) {
    y_pos <- legend_y - 0.04 - (i - 1) * line_spacing
    grid.lines(x = unit(c(legend_x, legend_x + 0.03), "npc"),
               y = unit(rep(y_pos, 2), "npc"),
               gp = gpar(col = legend_colors[i], lwd = 2))
    grid.text(labels[i],
              x = unit(legend_x + 0.04, "npc"),
              y = unit(y_pos, "npc"),
              just = "left", gp = gpar(fontsize = 9))
  }
}


# ----------------------------
# Direction-specific mega plots
# ----------------------------

mega_plot_forward <- function(dataTracks, file_info, cond_order, from, to, chr, color_map) {
  
  # Forward overlays only
  fwd_overlays <- lapply(cond_order, function(cond) {
    cond_tracks <- file_info %>%
      filter(condition == cond, eU_status != "total", direction == "fwd") %>%
      mutate(name = paste(condition, eU_status, direction, sep = " ")) %>%
      pull(name)
    
    cond_tracks <- cond_tracks[cond_tracks %in% names(dataTracks)]
    if (length(cond_tracks) == 0) return(NULL)
    
    OverlayTrack(
      trackList = lapply(cond_tracks, function(n) dataTracks[[n]]),
      name = cond,
      col.axis = "black",
      col.title = "black",
      background.title = "transparent",
      col.frame = "transparent"
    )
  })
  
  fwd_overlays <- fwd_overlays[!sapply(fwd_overlays, is.null)]
  
  # Combine forward overlays + gene model + genome axis
  all_tracks <- c(fwd_overlays, list(rtrackFor, gtrack))
  
  # Plot
  plotTracks(
    all_tracks,
    from = from, to = to, chromosome = chr,
    type = "l",
    cex.title = 1,
    cex.axis = 0.7,
    title.width = 1.5,
    track.padding = 0.1
  )
  
  # Legend
  legend_colors <- color_map[names(color_map)]
  labels <- names(legend_colors)
  legend_x <- 0.2
  legend_y <- 0.16
  line_spacing <- 0.03
  
  for (i in seq_along(labels)) {
    y_pos <- legend_y - 0.04 - (i - 1) * line_spacing
    grid.lines(x = unit(c(legend_x, legend_x + 0.03), "npc"),
               y = unit(rep(y_pos, 2), "npc"),
               gp = gpar(col = legend_colors[i], lwd = 6))
    grid.text(labels[i],
              x = unit(legend_x + 0.04, "npc"),
              y = unit(y_pos, "npc"),
              just = "left", gp = gpar(fontsize = 16))
  }
}


mega_plot_reverse <- function(dataTracks, file_info, cond_order, from, to, chr, color_map) {
  
  # Reverse overlays only
  rev_overlays <- lapply(rev(cond_order), function(cond) {
    cond_tracks <- file_info %>%
      filter(condition == cond, eU_status != "total", direction == "rev") %>%
      mutate(name = paste(condition, eU_status, direction, sep = " ")) %>%
      pull(name)
    
    cond_tracks <- cond_tracks[cond_tracks %in% names(dataTracks)]
    if (length(cond_tracks) == 0) return(NULL)
    
    OverlayTrack(
      trackList = lapply(cond_tracks, function(n) dataTracks[[n]]),
      name = cond,
      col.axis = "black",
      col.title = "black",
      background.title = "transparent",
      col.frame = "transparent"
    )
  })
  
  rev_overlays <- rev_overlays[!sapply(rev_overlays, is.null)]
  
  # Combine: genome axis and gene model at top, then overlays
  all_tracks <- c(list(gtrack, rtrackRev), rev_overlays)
  
  # Plot
  plotTracks(
    all_tracks,
    from = from, to = to, chromosome = chr,
    type = "l",
    cex.title = 1,
    cex.axis = 0.7,
    title.width = 1.5,
    track.padding = 0.1
  )
  
  # Legend
  legend_colors <- color_map[names(color_map)]
  labels <- names(legend_colors)
  legend_x <- 0.1
  legend_y <- 0.8
  line_spacing <- 0.03
  
  for (i in seq_along(labels)) {
    y_pos <- legend_y - 0.04 - (i - 1) * line_spacing
    grid.lines(x = unit(c(legend_x, legend_x + 0.03), "npc"),
               y = unit(rep(y_pos, 2), "npc"),
               gp = gpar(col = legend_colors[i], lwd = 6))
    grid.text(labels[i],
              x = unit(legend_x + 0.04, "npc"),
              y = unit(y_pos, "npc"),
              just = "left", gp = gpar(fontsize = 16))
  }
}




# ----------------------------
# Control conditions to plot
# ----------------------------
unique_conditions <- intersect(c("CHX", "2h", "4h DMSO", "6h", "8h"),
                               unique(file_info$condition))

# ----------------------------
# Output PDFs for original plots
# ----------------------------
pdf("GViz_megaplot_normalised.pdf", width = 20, height = 10)
mega_plot(dataTracks, file_info, unique_conditions, from = myStart, to = myEnd, chr = myChr, color_map = color_map)
dev.off()

pdf("GViz_megaplot_raw.pdf", width = 20, height = 10)
mega_plot(dataTracks_raw, file_info, unique_conditions, from = myStart, to = myEnd, chr = myChr, color_map = color_map)
dev.off()

pdf("GViz_megaplot_forward.pdf", width = 20, height = 10)
mega_plot_forward(dataTracks_raw, file_info, unique_conditions, from = myStart, to = myEnd, chr = myChr, color_map = color_map)
dev.off()

pdf("GViz_megaplot_reverse.pdf", width = 20, height = 10)
mega_plot_reverse(dataTracks_raw, file_info, unique_conditions, from = myStart, to = myEnd, chr = myChr, color_map = color_map)
dev.off()



# ----------------------------
# Custom plot using combined_color_map
# ----------------------------
pdf("GViz_megaplot_combined.pdf", width = 18, height = 10)
mega_plot_combined_normalised_overlay(file_info, cond_order = c("4h DMSO", "4h STM2457"),
                   from = myStart, to = myEnd, chr = myChr,
                   color_map = combined_color_map)
dev.off()

pdf("GViz_megaplot_combined_raw.pdf", width = 18, height = 10)
mega_plot_combined_raw(file_info, cond_order = c("4h DMSO", "4h STM2457"),
                   from = myStart, to = myEnd, chr = myChr,
                   color_map = combined_color_map)
dev.off()







# Your legend data
legend_colors <- combined_color_map[names(combined_color_map)]
labels <- names(legend_colors)

# Open PNG device (adjust width/height as needed)
png("legend_only_4h.png", width = 1000, height = 500, res = 300)

# Start a new grid page
grid.newpage()

# Legend layout parameters
legend_x <- 0.05
legend_y <- 0.9
line_spacing <- 0.18

# Draw legend
for (i in seq_along(labels)) {
  y_pos <- legend_y - (i - 1) * line_spacing
  
  grid.lines(
    x = unit(c(legend_x, legend_x + 0.18), "npc"),
    y = unit(rep(y_pos, 2), "npc"),
    gp = gpar(col = legend_colors[i], lwd = 9)
  )
  
  grid.text(
    labels[i],
    x = unit(legend_x + 0.23, "npc"),
    y = unit(y_pos + 0.01, "npc"),
    just = "left",
    gp = gpar(fontsize = 18)
  )
}

dev.off()

