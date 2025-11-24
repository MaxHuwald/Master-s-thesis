# Set working directory for outputs
setwd("/project/sysviro/users/Max/analyses/Conditions/Results/Modifications_pileup/coverage")

# Load libraries
library(dplyr)
library(readr)
library(stringr)
library(purrr)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(tidyr)
library(Gviz)
library(GenomicRanges)
library(GenomicFeatures)  # for makeTxDbFromGFF
library(data.table)       # using data.table for fast reading

base_dir <- "/project/sysviro/users/Max/analyses/Modkit"
# Genome/chromosome
myChr   <- "HSV2-st333"
myStart <- 1
myEnd   <- 155503

options(ucscChromosomeNames = FALSE)

# Gene models (forward/reverse)
modelsPos <- makeTxDbFromGFF("/project/sysviro/users/Max/Genome/HSV2-st333_Annotation_v4.forward.gff3")
modelsNeg <- makeTxDbFromGFF("/project/sysviro/users/Max/Genome/HSV2-st333_Annotation_v4.reverse.gff3")

rtrackFor <- GeneRegionTrack(
  modelsPos,
  genome           = myChr,
  chromosome       = myChr,
  name             = "Gene Model",
  col              = "black",
  fill             = "grey",
  stacking         = "squish",
  shape            = "smallArrow",
  background.title = "transparent"
)

rtrackRev <- GeneRegionTrack(
  modelsNeg,
  genome           = myChr,
  chromosome       = myChr,
  name             = "Gene Model",
  col              = "black",
  fill             = "grey",
  stacking         = "squish",
  shape            = "smallArrow",
  background.title = "transparent"
)

# Genome axis track
gtrack <- GenomeAxisTrack(
  col       = "black",
  cex       = 1.5,
  labelPos  = "above",
  fontcolor = "black"
)

# Find all BEDs
all_files <- list.files(
  base_dir,
  pattern   = "\\.bed$",
  recursive = TRUE,
  full.names = TRUE
)

# Parse filenames
parse_filename <- function(filepath) {
  folder <- basename(dirname(filepath))          # expects "eU" or "non_eU" (or "Total")
  fname  <- basename(filepath)
  
  condition <- str_extract(fname, "^[^.]+")
  
  alignment_label <- dplyr::case_when(
    str_detect(fname, "\\buntrimmed\\b") ~ "untrimmed",
    str_detect(fname, "aligned_hg38")    ~ "hg38",
    TRUE                                 ~ NA_character_
  )
  
  motif <- dplyr::case_when(
    str_detect(fname, "\\.pseU\\.bed$") ~ "pseU",
    str_detect(fname, "\\.m6a\\.bed$")  ~ "m6A",
    TRUE                                ~ NA_character_
  )
  
  cutoff <- stringr::str_match(fname, "\\.(\\d+(?:\\.\\d+)?)\\.filtered")[, 2]
  cutoff <- suppressWarnings(as.numeric(cutoff))
  
  tibble(
    filepath, folder, condition, alignment_label, motif, cutoff
  )
}

metadata <- purrr::map_dfr(all_files, parse_filename) %>%
  # keep only the targets you want
  dplyr::filter(
    folder %in% c("eU", "non_eU"),
    alignment_label == "untrimmed",
    !is.na(motif),
    !(condition == "4h"),                    # drop raw "4h"
    !(condition == "CHX" & folder == "eU")   # drop CHX in eU
  )

# Define Modkit BED column names
modkit_cols <- c(
  "chrom", "start", "end", "mod_code_motif", "score", "strand",
  "compat_start", "compat_end", "color", "N_valid_cov", "pct_mod",
  "N_mod", "N_canonical", "N_other_mod", "N_delete", "N_fail",
  "N_diff", "N_nocall"
)

# Read each BED and keep only rows matching the filename motif
# - pseU files  -> keep mod_code_motif == "17802"
# - m6A  files  -> keep mod_code_motif == "a"
# Also coerce types and overwrite mod_code_motif with clean label.

all_data <- metadata %>%
  mutate(data = purrr::map2(filepath, motif, ~{
    dt <- data.table::fread(.x, header = FALSE, col.names = modkit_cols)
    
    # filter by motif indicated in the filename
    if (.y == "pseU") {
      dt <- dt[mod_code_motif == "17802"]
    } else if (.y == "m6A") {
      dt <- dt[mod_code_motif == "a"]
    } else {
      dt <- dt[0]  # unknown motif -> empty
    }
    
    # enforce numeric types where needed
    dt[, `:=`(
      start   = as.integer(start),
      end     = as.integer(end),
      pct_mod = as.numeric(pct_mod)
    )]
    
    # replace raw code with clean motif label
    dt[, mod_code_motif := .y]
    
    dt
  })) %>%
  tidyr::unnest(data)

# Filter for coverage and valid pct_mod
all_data1 <- all_data %>%
  filter(N_valid_cov >= 50, !is.na(pct_mod))

# Recode condition names
all_data1 <- all_data1 %>%
  mutate(condition = recode(
    condition,
    "4h_DMSO"     = "4h DMSO",
    "4h_STM2457"  = "4h STM2457"
  ))

# Recode folder names
all_data1 <- all_data1 %>%
  mutate(folder = recode(
    folder,
    "non_eU" = "old mRNA",
    "eU"     = "new mRNA"
  ))

# Recode alignment label
all_data1 <- all_data1 %>%
  mutate(alignment_label = recode(
    alignment_label,
    "untrimmed" = "HSV2"
  ))

# Convert percentage values to 0–1 range
all_data1 <- all_data1 %>%
  mutate(pct_mod = pct_mod / 100)


all_data2 <- all_data1 %>%
  filter(!folder == "new mRNA")

# Order factor levels
condition_levels <- c("CHX", "2h", "4h DMSO", "4h STM2457", "6h", "8h")
all_data1$condition <- factor(all_data1$condition, levels = condition_levels)
all_data2$condition <- factor(all_data2$condition, levels = condition_levels)
# Stoichiometry colors (old vs new mRNA)
stoich_colors_m6a <- c(
  "old mRNA" = brewer.pal(5, "Reds")[5],
  "new mRNA" = brewer.pal(5, "Reds")[2]
)

stoich_colors_pseU <- c(
  "old mRNA" = brewer.pal(5, "Blues")[5],
  "new mRNA" = brewer.pal(5, "Blues")[2]
)


# ---------- AUTOMATED PLOTTING (m6A / pseU) --------------------------------

# Build an OverlayTrack for a single condition, given motif & strand
.build_condition_overlay <- function(df, condition, motif, strand, chr, from, to, col_old, col_new, type = "l") {
  # subset to this condition/motif/strand/window
  dsub <- df %>%
    dplyr::filter(
      .data$condition == !!condition,
      .data$mod_code_motif == !!motif,
      .data$chrom == !!chr,
      .data$strand == !!strand,
      .data$start < !!to,
      .data$end   > !!from,
      .data$folder %in% c("old mRNA", "new mRNA")
    )
  
  if (nrow(dsub) == 0) return(NULL)
  
  d_old <- dsub %>% dplyr::filter(folder == "old mRNA")
  d_new <- dsub %>% dplyr::filter(folder == "new mRNA")
  
  mkGR <- function(dd) {
    if (nrow(dd) == 0) return(GenomicRanges::GRanges())
    GenomicRanges::GRanges(
      seqnames = dd$chrom,
      ranges   = IRanges::IRanges(start = as.integer(dd$start) + 1L, end = as.integer(dd$end)),
      strand   = dd$strand,
      score    = as.numeric(dd$pct_mod)  # 0–1 already
    )
  }
  
  mkDT <- function(gr, col) {
    if (length(gr) == 0) return(NULL)
    tr <- Gviz::DataTrack(
      range      = gr,
      genome     = chr,
      chromosome = chr,
      name       = as.character(condition),  # never NULL
      type       = type,
      ylim       = c(0, 1),
      na.rm      = TRUE
    )
    Gviz::displayPars(tr) <- list(
      fill           = col,
      col.histogram  = col,
      col            = col,      # line color for type="l"
      col.grid       = NA,
      col.background = NA,
      baseline       = NA,
      col.baseline   = NA,
      alpha          = 1
    )
    tr
  }
  
  tr_old <- mkDT(mkGR(d_old), col_old)
  tr_new <- mkDT(mkGR(d_new), col_new)
  
  # If both exist, overlay; if only one, return that track
  tracks <- Filter(Negate(is.null), list(tr_old, tr_new))
  if (length(tracks) == 0) return(NULL)
  if (length(tracks) == 1) return(tracks[[1]])
  
  Gviz::OverlayTrack(
    trackList        = tracks,
    name             = as.character(condition),
    background.title = "transparent"
  )
}

.plot_mod_series <- function(df, motif, strand = "+", chr = myChr,
                             from = myStart, to = myEnd,
                             conditions = levels(df$condition),
                             type = "l", pdf_file = NULL) {
  stopifnot(motif %in% c("m6A","pseU"))
  # palettes
  if (motif == "m6A") {
    col_old <- unname(stoich_colors_m6a["old mRNA"])
    col_new <- unname(stoich_colors_m6a["new mRNA"])
  } else {
    col_old <- unname(stoich_colors_pseU["old mRNA"])
    col_new <- unname(stoich_colors_pseU["new mRNA"])
  }
  
  # If reverse strand, flip the condition order so it becomes:
  # 8h, 6h, 4h STM2457, 4h DMSO, 2h, CHX (given your factor levels)
  conds <- if (strand == "-") rev(conditions) else conditions
  
  # Build overlays in the chosen order
  overlays <- lapply(conds, function(cond) {
    .build_condition_overlay(df, cond, motif, strand, chr, from, to, col_old, col_new, type = type)
  })
  overlays <- Filter(Negate(is.null), overlays)
  if (length(overlays) == 0) stop("No data to plot for given parameters.")
  
  # Pick the gene model by strand
  gene_track <- if (strand == "+") rtrackFor else rtrackRev
  
  # Compose track list:
  # - forward: conditions (top→bottom), then gene, then genome axis (as you had)
  # - reverse: genome axis, gene, then conditions (top→bottom)  ✅
  if (strand == "-") {
    tracks <- c(list(gtrack, gene_track), overlays)
    sizes  <- c(0.8, 1.2, rep(1, length(overlays)))
  } else {
    tracks <- c(overlays, list(gene_track, gtrack))
    sizes  <- c(rep(1, length(overlays)), 1.2, 0.8)
  }
  
  title_txt <- sprintf("%s stoichiometry (%s strand)", motif, ifelse(strand=="+","forward","reverse"))
  
  if (!is.null(pdf_file)) pdf(pdf_file, width = 20, height = 10)
  on.exit(if (!is.null(pdf_file)) dev.off(), add = TRUE)
  
  Gviz::plotTracks(
    tracks,
    sizes            = sizes,
    from             = from,
    to               = to,
    chromosome       = chr,
    background.title = "transparent",
    cex.axis         = 1,
    col.axis         = "black"
    
  )
}

.plot_mod_series_combined <- function(
    df, motif,
    chr        = myChr,
    from       = myStart,
    to         = myEnd,
    conditions = levels(df$condition),
    type       = "l",
    pdf_file   = NULL
) {
  stopifnot(motif %in% c("m6A", "pseU"))
  
  # Choose color palette
  if (motif == "m6A") {
    col_old <- unname(stoich_colors_m6a["old mRNA"])
    col_new <- unname(stoich_colors_m6a["new mRNA"])
  } else {
    col_old <- unname(stoich_colors_pseU["old mRNA"])
    col_new <- unname(stoich_colors_pseU["new mRNA"])
  }
  
  # Forward strand: CHX → 8h (your normal factor order)
  conds_pos <- conditions
  
  # Reverse strand: flipped, 8h → CHX (as in your strand-specific plot)
  conds_neg <- rev(conditions)
  
  # Build overlays
  overlays_pos <- lapply(conds_pos, function(cond) {
    .build_condition_overlay(
      df        = df,
      condition = cond,
      motif     = motif,
      strand    = "+",
      chr       = chr,
      from      = from,
      to        = to,
      col_old   = col_old,
      col_new   = col_new,
      type      = type
    )
  })
  overlays_pos <- Filter(Negate(is.null), overlays_pos)
  
  overlays_neg <- lapply(conds_neg, function(cond) {
    .build_condition_overlay(
      df        = df,
      condition = cond,
      motif     = motif,
      strand    = "-",
      chr       = chr,
      from      = from,
      to        = to,
      col_old   = col_old,
      col_new   = col_new,
      type      = type
    )
  })
  overlays_neg <- Filter(Negate(is.null), overlays_neg)
  
  if (length(overlays_pos) == 0 && length(overlays_neg) == 0) {
    stop("No data to plot for given parameters.")
  }
  
  # Build track list: forward series → forward gene model → axis → reverse gene model → reverse series
  tracks <- c(
    overlays_pos,
    if (length(overlays_pos) > 0) list(rtrackFor) else list(),
    list(gtrack),
    if (length(overlays_neg) > 0) list(rtrackRev) else list(),
    overlays_neg
  )
  
  # Relative sizes: 1 for data tracks, 1.2 for gene models, 0.8 for axis
  sizes <- c(
    rep(1, length(overlays_pos)),
    if (length(overlays_pos) > 0) 1.2 else NULL,
    0.8,
    if (length(overlays_neg) > 0) 1.2 else NULL,
    rep(1, length(overlays_neg))
  )
  
  if (!is.null(pdf_file)) grDevices::pdf(pdf_file, width = 20, height = 10)
  on.exit(if (!is.null(pdf_file)) grDevices::dev.off(), add = TRUE)
  
  Gviz::plotTracks(
    tracks,
    sizes            = sizes,
    from             = from,
    to               = to,
    chromosome       = chr,
    background.title = "transparent",
    cex.axis         = 1,
    col.axis         = "black"
  )
}


# Public, motif-specific wrappers (hard-coded colors)
plot_m6A_series <- function(strand = "+", chr = myChr, from = myStart, to = myEnd,
                            conditions = levels(all_data1$condition),
                            type = "l", pdf_file = NULL) {
  .plot_mod_series(all_data1, motif = "m6A", strand = strand, chr = chr, from = from, to = to,
                   conditions = conditions, type = type, pdf_file = pdf_file)
}

plot_pseU_series <- function(strand = "+", chr = myChr, from = myStart, to = myEnd,
                             conditions = levels(all_data2$condition),
                             type = "l", pdf_file = NULL) {
  .plot_mod_series(all_data2, motif = "pseU", strand = strand, chr = chr, from = from, to = to,
                   conditions = conditions, type = type, pdf_file = pdf_file)
}

plot_m6A_combined_series <- function(
    chr        = myChr,
    from       = myStart,
    to         = myEnd,
    conditions = levels(all_data1$condition),
    type       = "l",
    pdf_file   = NULL
) {
  .plot_mod_series_combined(
    df         = all_data1,
    motif      = "m6A",
    chr        = chr,
    from       = from,
    to         = to,
    conditions = conditions,
    type       = type,
    pdf_file   = pdf_file
  )
}

plot_pseU_combined_series <- function(
    chr        = myChr,
    from       = myStart,
    to         = myEnd,
    conditions = levels(all_data2$condition),
    type       = "l",
    pdf_file   = NULL
) {
  .plot_mod_series_combined(
    df         = all_data2,
    motif      = "pseU",
    chr        = chr,
    from       = from,
    to         = to,
    conditions = conditions,
    type       = type,
    pdf_file   = pdf_file
  )
}


# 1) pseU • forward strand
plot_pseU_series(
  strand = "+", chr = myChr, from = myStart, to = myEnd,
  type = "l",                      # or "histogram"
  pdf_file = "pseU_forward_series_non_eU.pdf"
)

# 2) pseU • reverse strand
plot_pseU_series(
  strand = "-", chr = myChr, from = myStart, to = myEnd,
  type = "l",
  pdf_file = "pseU_reverse_series_non_eU.pdf"
)

# 3) m6A • forward strand
plot_m6A_series(
  strand = "+", chr = myChr, from = myStart, to = myEnd,
  type = "l",
  pdf_file = "m6A_forward_series.pdf"
)

# 4) m6A • reverse strand
plot_m6A_series(
  strand = "-", chr = myChr, from = myStart, to = myEnd,
  type = "l",
  pdf_file = "m6A_reverse_series.pdf"
)

# 5) pseU • combined (forward + axis + reverse)
plot_pseU_combined_series(
  chr      = myChr,
  from     = myStart,
  to       = myEnd,
  type     = "l",
  pdf_file = "pseU_combined_series_non_eU.pdf"
)

# 6) m6A • combined (forward + axis + reverse)
plot_m6A_combined_series(
  chr      = myChr,
  from     = myStart,
  to       = myEnd,
  type     = "l",
  pdf_file = "m6A_combined_series.pdf"
)


export_combined_legends <- function(
    file = "stoichiometry_legends.pdf",
    width = 6, height = 4,
    layout = c("vertical","horizontal"),
    # styling
    swatch_len = 0.35, line_spacing = 0.28, lwd = 8, cex_text = 2,
    show_box = TRUE
) {
  layout <- match.arg(layout)
  
  # Local renderer (no titles)
  .render_legend <- function(color_map) {
    if (show_box) {
      grid::grid.rect(gp = grid::gpar(col = "grey70", fill = "white", lwd = 1))
    }
    
    y0 <- 0.80
    labels <- names(color_map)
    for (i in seq_along(labels)) {
      y <- y0 - (i - 1) * line_spacing
      
      grid::grid.lines(
        x = grid::unit(c(0.08, 0.08 + swatch_len), "npc"),
        y = grid::unit(rep(y, 2), "npc"),
        gp = grid::gpar(col = unname(color_map[i]), lwd = lwd, lineend = "round")
      )
      grid::grid.text(
        labels[i],
        x = grid::unit(0.08 + swatch_len + 0.05, "npc"),
        y = grid::unit(y, "npc"),
        just = "left",
        gp = grid::gpar(cex = cex_text)
      )
    }
  }
  
  grDevices::pdf(file, width = width, height = height)
  on.exit(grDevices::dev.off(), add = TRUE)
  grid::grid.newpage()
  
  if (layout == "vertical") {
    vp <- grid::viewport(layout = grid::grid.layout(
      nrow = 2, ncol = 1,
      heights = grid::unit.c(grid::unit(0.5, "npc"), grid::unit(0.5, "npc"))
    ))
    grid::pushViewport(vp)
    
    grid::pushViewport(grid::viewport(layout.pos.row = 1))
    .render_legend(stoich_colors_m6a)
    grid::popViewport()
    
    grid::pushViewport(grid::viewport(layout.pos.row = 2))
    .render_legend(stoich_colors_pseU)
    grid::popViewport()
    
  } else { # horizontal
    vp <- grid::viewport(layout = grid::grid.layout(
      nrow = 1, ncol = 2,
      widths = grid::unit.c(grid::unit(0.5, "npc"), grid::unit(0.5, "npc"))
    ))
    grid::pushViewport(vp)
    
    grid::pushViewport(grid::viewport(layout.pos.col = 1))
    .render_legend(stoich_colors_m6a)
    grid::popViewport()
    
    grid::pushViewport(grid::viewport(layout.pos.col = 2))
    .render_legend(stoich_colors_pseU)
    grid::popViewport()
  }
}

# --- Example usage ---
# Stacked (vertical)
export_combined_legends(
  file = "stoichiometry_legends_vertical.pdf",
  layout = "vertical",
  width = 4, height = 5,
  line_spacing = 0.23   # tighter spacing (default was 0.28)
)


# Side-by-side (horizontal)
export_combined_legends(
  file = "stoichiometry_legends_horizontal.pdf",
  layout = "horizontal", width = 7, height = 3.5
)
