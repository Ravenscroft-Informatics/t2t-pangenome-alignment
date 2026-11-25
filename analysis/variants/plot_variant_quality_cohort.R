# ----------------------------
# Variant quality score plots
# ----------------------------

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(patchwork)
library(forcats)

save_pdf <- function(plot, file, width = 12, height = 5) {
  ggsave(filename = file.path(out_dir, file), plot = plot, device = cairo_pdf,
         width = width, height = height, units = "in")
}

# ---- Paths / inputs ----
setwd("/Users/00104561/Library/CloudStorage/OneDrive-UWA/Research/Projects/AIM1_DATA-REANALYSIS/benchmark/variant_quality_score/")

source("/Users/00104561/Library/CloudStorage/OneDrive-UWA/Research/Projects/AIM1_DATA-REANALYSIS/benchmark/variant_quality_score/code/parse_bcftools_stats_vqc.R")

data_lines <- readLines("./data/hprc_bcftools_stats_vqc_cohort.tsv")

# Assay labels must match the order of rows after transpose (see sanity check below).
assay_labels <- c(
  "genome", "clinvar", "exome", "mask", "syntenic",
  "genome", "clinvar", "exome", "mask"
)
assay <- c("genome", "exome", "clinvar", "mask", "genome", "exome", "clinvar", "mask")

# ---- Build tidy dataframe ----
df <- create_quality_dataframe(data_lines)

# Transpose: rows become samples/assays, columns become quality score bins
df_t <- as.data.frame(t(df))

# Safety check: ensure assay_labels length matches transposed row count
stopifnot(nrow(df_t) == length(assay_labels))

df_t <- df_t %>%
  mutate(
    sample_id = sub("\\..*", "", rownames(df_t)),
    genome    = if_else(str_detect(rownames(df_t), "chm13|T2T"), "chm13", "grch38"),
    assay     = assay_labels
  ) %>%
  select(sample_id, genome, assay, everything())

rownames(df_t) <- NULL

# Treat 0 counts as missing (so they don't clutter plots if that's intended)
df_t[df_t == 0] <- NA

df_t$genome <- fct_recode(df_t$genome, "CHM13" = "chm13")
df_t$genome <- fct_recode(df_t$genome, "GRCh38" = "grch38")

# ---- Subsets used for plots ----
masked   <- filter(df_t, assay == "mask")
exome    <- filter(df_t, assay == "exome")
genome_w <- filter(df_t, assay == "genome")
clinvar  <- filter(df_t, assay == "clinvar")

# Your original logic: syntenic plus (grch38 genome) in the same panel
syntenic <- filter(df_t, assay == "syntenic" | (assay == "genome" & genome == "GRCh38"))

# ---- Helper: convert to long + plot ----
make_plot <- function(df_subset,
                      title,
                      y_max,
                      show_legend = FALSE,
                      x_break_by = 10) {
  
  plot_df <- df_subset %>%
    pivot_longer(
      cols = 4:ncol(df_subset),
      names_to = "quality_score",
      values_to = "value"
    ) %>%
    mutate(quality_score = as.numeric(quality_score))
  
  # Define x breaks based on this subset (fixes the plot_df bug in your script)
  x_min <- min(plot_df$quality_score, na.rm = TRUE)
  x_max <- max(plot_df$quality_score, na.rm = TRUE)
  
  p <- ggplot(plot_df, aes(x = quality_score, y = value, color = genome, group = genome)) +
    geom_line(linewidth = 1) +
    labs(
      title = title,
      x = "Quality score",
      y = "Number of variants",
      color = "Genome"
    ) +
    scale_x_continuous(breaks = seq(x_min, x_max, by = x_break_by)) +
    scale_y_continuous(limits = c(0, y_max)) +
    scale_colour_manual(values=c("pink", "skyblue")) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 1, vjust = 0.5),
      panel.grid.minor = element_blank()
    )
  
  if (!show_legend) {
    p <- p + theme(legend.position = "none")
  }
  
  p
}

# ---- Make plots ----
genome_plot  <- make_plot(genome_w, "Genome-wide Variants", 600000, show_legend = FALSE)
exome_plot   <- make_plot(exome,    "Exome Variants",        17000, show_legend = FALSE)
synteny_plot <- make_plot(syntenic, "Syntenic Variants",    600000, show_legend = FALSE)
mask_plot    <- make_plot(masked,   "Genome Mask Variants", 600000, show_legend = FALSE)
clinvar_plot <- make_plot(clinvar,  "ClinVar Variants",       7000, show_legend = TRUE)

# ---- Combine ----
plot_linear <- (genome_plot + exome_plot + synteny_plot + mask_plot + clinvar_plot) +
  plot_layout(ncol = 2)

# Output directory (created if missing)
out_dir <- "plots"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

save_pdf(plot_linear, "./linear_vqc_comparison.pdf", width = 10, height = 8)

#hprc 
plot_hprc <- (genome_plot + exome_plot + synteny_plot + mask_plot + clinvar_plot) +
  plot_layout(ncol = 2)

save_pdf(plot_hprc, "./hprc_vqc_comparison.pdf", width = 10, height = 8)
