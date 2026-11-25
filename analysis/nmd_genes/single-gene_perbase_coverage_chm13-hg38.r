##Use this scrip to generate a plot of the per-base coverage across a gene of interest

# --- packages ---
library(readr)
library(dplyr)
library(ggplot2)
library(patchwork)
library(grid)
library(IRanges)
library(shadowtext)
library(cowplot)

setwd("/Users/00104561/Library/CloudStorage/OneDrive-UWA/Research/Projects/AIM1_DATA-REANALYSIS/nmd_genes/")


# --- helpers ---------------------------------------------------------------

agg_cov <- function(df) {
  df %>%
    mutate(depth = as.numeric(depth)) %>%
    group_by(gene, chr, start, end) %>%
    summarise(
      n_samples    = n(),
      mean_depth   = mean(depth, na.rm = TRUE),
      median_depth = median(depth, na.rm = TRUE),
      q25          = quantile(depth, 0.25, na.rm = TRUE),
      q75          = quantile(depth, 0.75,  na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(start)
}

to_bed_style <- function(exons_df, cov_df, force = c("auto","gtf1","bed0")) {
  force <- match.arg(force)
  ex <- exons_df %>% mutate(start = as.integer(start), end = as.integer(end))
  if (force == "gtf1") { ex$start <- ex$start - 1L; return(ex) }
  if (force == "bed0") return(ex)
  cov_min <- suppressWarnings(min(cov_df$start, na.rm = TRUE))
  ex_min  <- suppressWarnings(min(ex$start,     na.rm = TRUE))
  if (is.finite(cov_min) && is.finite(ex_min) && (ex_min == cov_min + 1L)) {
    ex$start <- ex$start - 1L
  }
  ex
}

master_gene_model <- function(exons_df, x_min, x_max) {
  ex_ir     <- IRanges::IRanges(start = as.integer(exons_df$start),
                                end   = as.integer(exons_df$end))
  ex_merged <- IRanges::reduce(ex_ir)
  master_exons <- data.frame(
    start = as.integer(IRanges::start(ex_merged)),
    end   = as.integer(IRanges::end(ex_merged))
  ) |>
    dplyr::transmute(start = pmax(start, x_min),
                     end   = pmin(end,   x_max)) |>
    dplyr::filter(end > start) |>
    dplyr::arrange(start)

  master_introns <- if (nrow(master_exons) >= 2) {
    data.frame(start = master_exons$end[-nrow(master_exons)],
               end   = master_exons$start[-1])
  } else master_exons[0, c("start","end")]

  strand_master <- exons_df |>
    dplyr::count(strand) |>
    dplyr::arrange(dplyr::desc(n)) |>
    dplyr::slice_head(n = 1) |>
    dplyr::pull(strand)
  if (length(strand_master) == 0 || is.na(strand_master)) strand_master <- "+"

  list(exons = master_exons, introns = master_introns, strand = strand_master)
}

make_arrow_df <- function(df, strand = "+", step = NULL, pad = 4L) {
  if (nrow(df) == 0) return(df[0, ])
  out <- vector("list", nrow(df))
  for (i in seq_len(nrow(df))) {
    s <- df$start[i]; e <- df$end[i]
    if (is.na(s) || is.na(e) || e - s <= (pad * 2)) next
    local_step <- if (is.null(step)) max(round((e - s) / 20), 50L) else step
    if (strand == "+") { xs <- seq(s + pad, e - pad, by = local_step); xe <- pmin(xs + pad, e - 1L)
    } else              { xe <- seq(e - pad, s + pad, by = -local_step); xs <- pmax(xe - pad, s + 1L) }
    out[[i]] <- data.frame(x = xs, xend = xe)
  }
  dplyr::bind_rows(out)
}

plot_cov <- function(df, agg, xlim, y_cap, y_breaks, show_xticks = FALSE, title_suffix = "") {
  mean_depth_all <- mean(agg$mean_depth, na.rm = TRUE)

  # Average SD of depth across samples
  avg_sd_all <- df %>%
    dplyr::group_by(sample) %>%
    dplyr::summarise(sd_depth = sd(depth, na.rm = TRUE), .groups = "drop") %>%
    dplyr::summarise(avg_sd = mean(sd_depth, na.rm = TRUE)) %>%
    dplyr::pull(avg_sd)

  # positions for the labels (top-right)
  span <- diff(xlim)
  x_right <- xlim[2] - 0.01 * span  # small padding from the right edge

  ggplot(agg) +
    geom_rect(aes(xmin = start, xmax = end, ymin = q25, ymax = q75), alpha = 0.2) +
    geom_segment(aes(x = start, xend = end, y = mean_depth, yend = mean_depth), linewidth = 0.4) +
    # Global mean line across the span
    geom_hline(yintercept = mean_depth_all, colour = "red", linewidth = 0.4, linetype = "dashed") +
    # labels (no shadowtext)
    annotate("label",
             x = x_right, y = Inf, hjust = 1, vjust = 1.2,
             label = sprintf("Mean depth: %.1f×", mean_depth_all),
             fill = "white", alpha = 0.85, label.size = 0,
             colour = "red", fontface = "bold", size = 3.3) +
    annotate("label",
             x = x_right, y = Inf, hjust = 1, vjust = 2.6,
             label = sprintf("Mean SD: %.1f×", avg_sd_all),
             fill = "white", alpha = 0.85, label.size = 0,
             colour = "black", fontface = "bold", size = 3.3) +
    coord_cartesian(xlim = xlim, ylim = c(0, y_cap), expand = FALSE) +
    scale_y_continuous(breaks = y_breaks, limits = c(0, y_cap), expand = expansion()) +
    labs(
      title = NULL,
      x = if (show_xticks) paste0(unique(agg$chr), " position") else NULL,
      y = "Depth (mean, IQR band)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title   = element_text(face = "bold"),
      axis.text.x  = if (show_xticks) element_text() else element_blank(),
      axis.ticks.x = if (show_xticks) element_line() else element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin  = margin(2, 4, 2, 4)
    )
}

plot_master <- function(model, xlim, xlab = NULL) {
  arrows_df <- make_arrow_df(model$introns, strand = model$strand)
  ggplot() +
    { if (nrow(model$introns))
      geom_segment(data = model$introns,
                   aes(x = start, xend = end, y = 1, yend = 1), linewidth = 0.35) else NULL } +
    { if (nrow(arrows_df))
      geom_segment(data = arrows_df,
                   aes(x = x, xend = xend, y = 1, yend = 1),
                   arrow = arrow(length = unit(0.02, "npc")), linewidth = 0.3) else NULL } +
    { if (nrow(model$exons))
      geom_segment(data = model$exons,
                   aes(x = start, xend = end, y = 1, yend = 1), linewidth = 6) else NULL } +
    coord_cartesian(xlim = xlim, expand = FALSE) +
    scale_y_continuous(limits = c(0.5, 1.5)) +
    labs(x = xlab, y = NULL) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid   = element_blank(),
      plot.margin  = margin(2, 5, 0, 5)
    )
}

lims_from_union <- function(df_cov, ex_df) {
  xmin <- min(min(df_cov$start, na.rm=TRUE), min(ex_df$start, na.rm=TRUE))
  xmax <- max(max(df_cov$end,   na.rm=TRUE), max(ex_df$end,   na.rm=TRUE))
  c(xmin, xmax)
}
# --- load data -------------------------------------------------------------
gene <- "ACTN2"
df_hg38      <- read_tsv(paste0("./data/", gene, ".perbase_mosdepth_hg38.tsv"),  show_col_types = FALSE)
exons_hg38   <- read_tsv(paste0("./data/", gene, ".grch38.exons.tsv"),           show_col_types = FALSE)
df_chm13     <- read_tsv(paste0("./data/", gene, ".perbase_mosdepth_chm13.tsv"), show_col_types = FALSE)
exons_chm13  <- read_tsv(paste0("./data/", gene, ".chm13.exons.tsv"),            show_col_types = FALSE)

#run to exclude sample
exclude <- "D21-0091"
df_hg38  <- df_hg38  %>% filter(sample != exclude)
df_chm13 <- df_chm13 %>% filter(sample != exclude)

# Harmonise exon coords to BED-style against each coverage set
exons_hg38  <- to_bed_style(exons_hg38,  df_hg38,  force = "auto")
exons_chm13 <- to_bed_style(exons_chm13, df_chm13, force = "auto")

# --- per-reference aggregates & models ------------------------------------

agg_38   <- agg_cov(df_hg38)
agg_t2t  <- agg_cov(df_chm13)

# derive union x-lims (coverage ∪ exons), then equalise span across references
xlim_38_union  <- lims_from_union(df_hg38,  exons_hg38)
xlim_t2t_union <- lims_from_union(df_chm13, exons_chm13)
span_max <- max(diff(xlim_38_union), diff(xlim_t2t_union))

# anchor left (set to "center" if preferred)
xlim_38  <- c(xlim_38_union[1],  xlim_38_union[1]  + span_max)
xlim_t2t <- c(xlim_t2t_union[1], xlim_t2t_union[1] + span_max)

# build master models clipped to those x-lims
model_38  <- master_gene_model(exons_hg38,  xlim_38[1],  xlim_38[2])
model_t2t <- master_gene_model(exons_chm13, xlim_t2t[1], xlim_t2t[2])

# shared y-scale (use max of mean & upper IQR across both)
y_cap <- max(
  max(agg_38$mean_depth,  na.rm = TRUE),
  max(agg_t2t$mean_depth, na.rm = TRUE),
  max(agg_38$q75,         na.rm = TRUE),
  max(agg_t2t$q75,        na.rm = TRUE)
)
y_breaks <- pretty(c(0, y_cap), n = 5)

# --- plots ----------------------------------------------------------------

p38_cov  <- plot_cov(df_hg38,  agg_38,  xlim_38,  y_cap, y_breaks, show_xticks = FALSE, title_suffix = " (GRCh38)")
p38_master <- plot_master(model_38, xlim_38, xlab = NULL)

pT2T_cov <- plot_cov(df_chm13, agg_t2t, xlim_t2t, y_cap, y_breaks, show_xticks = FALSE, title_suffix = " (CHM13)")
pT2T_master <- plot_master(model_t2t, xlim_t2t, xlab = NULL)

# --- build two-row stacks for each reference (unchanged) ---
hg38_stack  <- patchwork::wrap_plots(p38_cov,  p38_master,  ncol = 1, heights = c(3, 1))
chm13_stack <- patchwork::wrap_plots(pT2T_cov, pT2T_master, ncol = 1, heights = c(3, 1))

# remove y-axis titles from individual panels
hg38_stack  <- hg38_stack  & theme(axis.title.y = element_blank())
chm13_stack <- chm13_stack & theme(axis.title.y = element_blank())

# combine the two stacks vertically
combined <- cowplot::plot_grid(
  hg38_stack, chm13_stack,
  ncol = 1, rel_heights = c(1, 1), align = "v"
)

# build final figure with spacing for labels - HPRC
final_plot <- cowplot::ggdraw() +
  # shift plots inward to make room for left y-axis label
  cowplot::draw_plot(combined, x = 0.03, y = 0, width = 0.88, height = 0.9) +
  # left-aligned main title and subtitle
  cowplot::draw_label(
    paste0(unique(df_hg38$gene), " coverage (Linear vs HPRC)"),
    x = 0.03, y = 0.99, hjust = 0, vjust = 1,
    fontface = "bold", size = 14
  ) +
  cowplot::draw_label(
    "Mean depth with IQR band; dashed red line = mean coverage across span",
    x = 0.03, y = 0.96, hjust = 0, vjust = 1, size = 10
  ) +
  # shared y-axis label
  cowplot::draw_label(
    "Coverage (mean, IQR band)",
    angle = 90, x = 0.02, y = 0.5, vjust = 0.5, size = 12, fontface = "bold"
  )


pdf(paste0("./", gene, ".grch38-linear-hprc_perbase_cov.pdf"), width = 10, height = 8, )
final_plot
dev.off()

ggplot2::ggsave(
  filename = sprintf("%s.grch38-linear-hprc_perbase_cov.cairo.pdf", gene),
  plot     = final_plot,
  device   = cairo_pdf,          # needs the Cairo package installed
  width    = 10, height = 8, units = "in",
  bg       = "white",
  onefile  = FALSE              # a single page
  #useDingbats = FALSE            # safer font embedding
)

# high-res, flattened raster export
ggplot2::ggsave(
  filename = sprintf("%s.grch38-linear-hprc_perbase_cov.png", gene),
  plot     = final_plot,
  width    = 10, height = 4, units = "in", dpi = 1000, bg = "white"
)
# (Optional) quick sanity check:
min(df_hg38$start);  min(exons_hg38$start)
min(df_chm13$start); min(exons_chm13$start)


# --- load data -------------------------------------------------------------
gene <- "COX6A1"

# both per-base coverages are hg38-aligned
df_linear <- read_tsv(paste0("../../../nmd_genes/", gene, ".perbase_mosdepth_linear_hg38.tsv"), show_col_types = FALSE)
df_hprc   <- read_tsv(paste0("../../../nmd_genes/", gene, ".perbase_mosdepth_hprc_hg38.tsv"),   show_col_types = FALSE)

# single exon set for hg38
exons_hg38 <- read_tsv(paste0("../../../nmd_genes/", gene, ".grch38.exons.tsv"), show_col_types = FALSE)

# --- optionally exclude a sample ------------------------------------------
exclude <- "D21-0091"
df_linear <- df_linear %>% dplyr::filter(sample != exclude)
df_hprc   <- df_hprc   %>% dplyr::filter(sample != exclude)

# Harmonise exon coords to BED-style against each coverage set (keeps them consistent on hg38)
exons_hg38_lin  <- to_bed_style(exons_hg38, df_linear, force = "auto")
exons_hg38_hprc <- to_bed_style(exons_hg38, df_hprc,   force = "auto")

# --- aggregates ------------------------------------------------------------
agg_linear <- agg_cov(df_linear)
agg_hprc   <- agg_cov(df_hprc)

# --- shared x-axis (same reference/gene model) -----------------------------
# take the union of (coverage ∪ exons) across both datasets, then use one shared span
xlim_lin_union  <- lims_from_union(df_linear, exons_hg38_lin)
xlim_hprc_union <- lims_from_union(df_hprc,   exons_hg38_hprc)

x_start <- min(xlim_lin_union[1],  xlim_hprc_union[1],  na.rm = TRUE)
x_end   <- max(xlim_lin_union[2],  xlim_hprc_union[2],  na.rm = TRUE)
xlim_shared <- c(x_start, x_end)

# build one master gene model on hg38 and clip to the shared limits
model_hg38 <- master_gene_model(exons_hg38, xlim_shared[1], xlim_shared[2])

# --- shared y-scale --------------------------------------------------------
y_cap <- max(
  max(agg_linear$mean_depth, na.rm = TRUE),
  max(agg_hprc$mean_depth,   na.rm = TRUE),
  max(agg_linear$q75,        na.rm = TRUE),
  max(agg_hprc$q75,          na.rm = TRUE)
)
y_breaks <- pretty(c(0, y_cap), n = 5)

# --- plots -----------------------------------------------------------------
# Suffixes clarify method; both are hg38
pLinear_cov   <- plot_cov(df_linear, agg_linear, xlim_shared, y_cap, y_breaks,
                          show_xticks = FALSE, title_suffix = " (Linear, GRCh38)")
pLinear_model <- plot_master(model_hg38, xlim_shared, xlab = NULL)

pHPRC_cov     <- plot_cov(df_hprc,   agg_hprc,   xlim_shared, y_cap, y_breaks,
                          show_xticks = FALSE, title_suffix = " (HPRC, GRCh38)")
pHPRC_model   <- plot_master(model_hg38, xlim_shared,
                             xlab = paste0(unique(df_linear$chr), " position"))

# two-row stacks (coverage + gene model) for each method
linear_stack <- patchwork::wrap_plots(pLinear_cov, pLinear_model, ncol = 1, heights = c(3, 1))
hprc_stack   <- patchwork::wrap_plots(pHPRC_cov,   pHPRC_model,   ncol = 1, heights = c(3, 1))

# remove per-panel y-axis titles; we’ll draw one shared label later
linear_stack <- linear_stack & theme(axis.title.y = element_blank())
hprc_stack   <- hprc_stack   & theme(axis.title.y = element_blank())

# stack methods vertically
combined <- cowplot::plot_grid(
  linear_stack, hprc_stack,
  ncol = 1, rel_heights = c(1, 1), align = "v"
)

# final figure with shared y-axis label and title
final_plot <- cowplot::ggdraw() +
  cowplot::draw_plot(combined, x = 0.03, y = 0, width = 0.88, height = 0.9) +
  cowplot::draw_label(
    paste0(unique(df_linear$gene), " coverage (Linear vs HPRC, GRCh38)"),
    x = 0.03, y = 0.99, hjust = 0, vjust = 1, fontface = "bold", size = 14
  ) +
  cowplot::draw_label(
    "Mean depth with IQR band; dashed red line = mean coverage across span",
    x = 0.03, y = 0.96, hjust = 0, vjust = 1, size = 10
  ) +
  cowplot::draw_label(
    "Depth (mean, IQR band)",
    angle = 90, x = 0.02, y = 0.5, vjust = 0.5, size = 12, fontface = "bold"
  )

pdf(paste0("./", gene, ".grch38_linear_vs_hprc_perbase_cov.pdf"), width = 10, height = 8)
final_plot
dev.off()

# (Optional sanity checks)
min(df_linear$start);  min(exons_hg38_lin$start)
min(df_hprc$start);    min(exons_hg38_hprc$start)
