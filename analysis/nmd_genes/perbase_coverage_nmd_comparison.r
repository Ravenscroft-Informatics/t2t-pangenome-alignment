# ==============================================================================
# Script: nmd_perbase_compare_chm13_vs_grch38.R
#
# Purpose:
#   Compare the coverage and coverage uniformity of NMD genes between GRCh38 and CHM13.
#
#   Given per-gene summary statistics (derived from per-base mosdepth outputs),
#   this script:
#     1) Reads GRCh38 and CHM13 per-gene per-sample summary tables.
#     2) Harmonises sample IDs (strips ".hg38"/".chm13" suffixes) and labels
#        references as GRCh38 and CHM13 (ordered: CHM13 first).
#     3) Computes averages across all genes:
#          - per sample per reference
#          - overall per reference
#     4) Computes reference-to-reference differences per gene per sample for:
#          - mean depth (CHM13 − GRCh38)
#          - SD depth   (CHM13 − GRCh38)
#     5) Summarises per-gene differences across samples and selects top-N genes
#        where CHM13 or GRCh38 has higher mean coverage and less uniform coverage.
#     6) Generates publication-ready PDF plots and exports tidy TSV outputs.
#
# Inputs (TSV; per-gene summaries derived from per-base mosdepth):
#   - ./data/all.perbase_mosdepth.summary.linear_hg38.tsv
#   - ./data/all.perbase_mosdepth.summary.linear_chm13.tsv
#
# Expected columns (minimum):
#   gene, sample,
#   min_depth, max_depth, range_depth, mean_depth, median_depth, sd_depth
#
# Outputs (written to ./plots/):
#   PDFs:
#     - ave_cov_all_nmd_genes.pdf
#     - ave_sd_all_nmd_genes.pdf
#     - scatter_sd_chm13_vs_grch38_persample.pdf
#     - scatter_sd_chm13_vs_grch38_mean.pdf
#     - scatter_mean_chm13_vs_grch38_mean.pdf
#     - scatter_mean_chm13_vs_grch38.pdf
#     - t2t_genes_higher_mean_depth.pdf
#     - grch38_genes_higher_mean_depth.pdf
#
#   TSVs:
#     - persample_summary_across_nmd_genes.tsv
#     - total_summary_by_reference.tsv
#     - compare_mean_depth_per_genexsample.tsv
#     - compare_sd_depth_per_genexsample.tsv
#     - gene_summary_mean_depth.tsv
#     - gene_summary_sd_depth.tsv
#
# Key assumptions:
#   - Differences are defined as CHM13 − GRCh38.
#   - Sample IDs must match across references after suffix stripping.
#   - Metrics are already per-gene summaries (this script does NOT use interval-
#     level mosdepth output).
#   - Genes are labelled in scatter plots only where |mean_diff| exceeds a
#     threshold (>2 for SD depth; >1 for mean depth).
#   - D21-0091 is excluded from all analyses.
#
# Author: Chiara Folland
# ==============================================================================

setwd("~/Library/CloudStorage/OneDrive-UWA/Research/Projects/AIM1_DATA-REANALYSIS/benchmark/batch1/nmd_genes/")

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(forcats)
  library(ggplot2)
  library(scales)
  library(ggrepel)
})

# ------------------------------- Config ---------------------------------------
path_hg38  <- "./data/all.perbase_mosdepth.summary.linear_hg38.tsv"
path_chm13 <- "./data/all.perbase_mosdepth.summary.linear_chm13.tsv"

options(pillar.sigfig = 6)

out_dir <- "plots"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

exclude_sample <- "D21-0091"
top_n <- 20

ref_levels <- c("CHM13", "GRCh38")
ref_colors <- c(CHM13 = "#f5cccd", GRCh38 = "#abc3eb")

# If you use aspect ratios in themes, set a default here:
asrt <- 1.2

# ---------------------------- Helper Functions --------------------------------
assert_readable <- function(p) {
  if (!file.exists(p)) stop("File not found: ", p, call. = FALSE)
}

save_pdf <- function(plot, filename, width = 7, height = 5) {
  ggsave(
    filename = file.path(out_dir, filename),
    plot = plot,
    device = cairo_pdf,
    width = width,
    height = height,
    units = "in"
  )
}

# Ensure expected columns exist
assert_cols <- function(df, cols) {
  missing <- setdiff(cols, names(df))
  if (length(missing)) {
    stop("Missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
}

# ------------------------------- Load Data ------------------------------------
assert_readable(path_hg38)
assert_readable(path_chm13)

df_hg38  <- read_tsv(path_hg38,  show_col_types = FALSE)
df_chm13 <- read_tsv(path_chm13, show_col_types = FALSE)

required_cols <- c(
  "gene", "sample",
  "min_depth", "max_depth", "range_depth", "mean_depth", "median_depth", "sd_depth"
)
assert_cols(df_hg38, required_cols)
assert_cols(df_chm13, required_cols)

# Harmonise sample IDs + label reference
df_hg38 <- df_hg38 %>%
  mutate(
    reference = "GRCh38",
    sample = str_remove(sample, "\\.hg38$")
  )

df_chm13 <- df_chm13 %>%
  mutate(
    reference = "CHM13",
    sample = str_remove(sample, "\\.chm13$")
  )

# Exclude sample (if present)
df_hg38  <- df_hg38  %>% filter(sample != exclude_sample)
df_chm13 <- df_chm13 %>% filter(sample != exclude_sample)

# Harmonise sample IDs + label reference
df_hg38 <- df_hg38 %>%
  mutate(
    reference = "GRCh38",
    sample = str_remove(sample, "\\.hg38$")
  )

df_chm13 <- df_chm13 %>%
  mutate(
    reference = "CHM13",
    sample = str_remove(sample, "\\.chm13$")
  )

df_combined <- bind_rows(df_hg38, df_chm13) %>%
  mutate(reference = factor(reference, levels = ref_levels)) 

# ----------------------- Averages Across All Genes ----------------------------
# Per sample per reference averages across genes
all_nmd_sample_summary <- df_combined %>%
  group_by(sample, reference) %>%
  summarise(
    mean_min_depth    = round(mean(min_depth, na.rm = TRUE), 2),
    mean_max_depth    = round(mean(max_depth, na.rm = TRUE), 2),
    mean_range_depth  = round(mean(range_depth, na.rm = TRUE), 2),
    mean_mean_depth   = round(mean(mean_depth, na.rm = TRUE), 2),
    mean_median_depth = round(mean(median_depth, na.rm = TRUE), 2),
    mean_sd_depth     = round(mean(sd_depth, na.rm = TRUE), 2),
    n_genes = n(),
    .groups = "drop"
  ) %>%
  arrange(sample, reference)

write_tsv(all_nmd_sample_summary, file.path(out_dir, "persample_summary_across_nmd_genes.tsv"))

# Overall per reference averages across all gene×sample rows
total_summary <- df_combined %>%
  group_by(reference) %>%
  summarise(
    mean_min_depth    = mean(min_depth, na.rm = TRUE),
    mean_max_depth    = mean(max_depth, na.rm = TRUE),
    mean_range_depth  = mean(range_depth, na.rm = TRUE),
    mean_mean_depth   = mean(mean_depth, na.rm = TRUE),
    mean_median_depth = mean(median_depth, na.rm = TRUE),
    mean_sd_depth     = mean(sd_depth, na.rm = TRUE),
    n_rows = n(),
    .groups = "drop"
  ) %>%
  arrange(reference)

write_tsv(total_summary, file.path(out_dir, "total_summary_by_reference.tsv"))

# --------------------- Reference-to-reference comparisons ---------------------
# Wide comparisons per gene×sample
wide_mean <- df_combined %>%
  select(gene, sample, reference, mean_depth) %>%
  pivot_wider(names_from = reference, values_from = mean_depth)

wide_sd <- df_combined %>%
  select(gene, sample, reference, sd_depth) %>%
  pivot_wider(names_from = reference, values_from = sd_depth)

# Sanity check: both references present for each row
# (If there are missing values, diffs will be NA; that's usually preferable to stopping.)
compare_mean <- wide_mean %>%
  mutate(
    diff_mean = CHM13 - GRCh38,
  )

compare_sd <- wide_sd %>%
  mutate(
    diff_sd = CHM13 - GRCh38,
  )

write_tsv(compare_mean, file.path(out_dir, "compare_mean_depth_per_genexsample.tsv"))
write_tsv(compare_sd,   file.path(out_dir, "compare_sd_depth_per_genexsample.tsv"))

# Per-gene summaries across samples (mean depth)
gene_summary_mean <- compare_mean %>%
  group_by(gene) %>%
  summarise(
    mean_GRCh38 = mean(GRCh38, na.rm = TRUE),
    mean_CHM13  = mean(CHM13,  na.rm = TRUE),
    mean_diff   = mean(diff_mean, na.rm = TRUE),
    median_diff = median(diff_mean, na.rm = TRUE),
    mean_pct_change = ifelse(is.finite(mean_GRCh38) & mean_GRCh38 != 0,
                             100 * (mean_CHM13 - mean_GRCh38) / mean_GRCh38, NA_real_),
    n_samples = sum(is.finite(GRCh38) & is.finite(CHM13)),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_diff))

write_tsv(gene_summary_mean, file.path(out_dir, "gene_summary_mean_depth.tsv"))

# Per-gene summaries across samples (SD depth)
gene_summary_sd <- compare_sd %>%
  group_by(gene) %>%
  summarise(
    mean_GRCh38 = mean(GRCh38, na.rm = TRUE),
    mean_CHM13  = mean(CHM13,  na.rm = TRUE),
    mean_diff   = mean(diff_sd, na.rm = TRUE),
    median_diff = median(diff_sd, na.rm = TRUE),
    mean_pct_change = ifelse(is.finite(mean_GRCh38) & mean_GRCh38 != 0,
                             100 * (mean_CHM13 - mean_GRCh38) / mean_GRCh38, NA_real_),
    n_samples = sum(is.finite(GRCh38) & is.finite(CHM13)),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_diff))

write_tsv(gene_summary_sd, file.path(out_dir, "gene_summary_sd_depth.tsv"))

# ------------------------------- Plots ----------------------------------------
# 1) Boxplot: per-sample mean coverage across NMD genes (mean of gene means)
p_ave_cov <- ggplot(all_nmd_sample_summary, aes(x = reference, y = mean_mean_depth, fill = reference)) +
  geom_boxplot(outlier.shape = 21, outlier.size = 2) +
  stat_summary(fun = mean, geom = "point", color = "black", size = 2.2, shape = 23) +
  scale_fill_manual(values = ref_colors) +
  scale_y_continuous(labels = label_number(accuracy = 0.01)) +
  labs(x = NULL, y = "Mean coverage across NMD genes") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none", aspect.ratio = asrt)

save_pdf(p_ave_cov, "ave_cov_all_nmd_genes.pdf", width = 4.5, height = 6)

# 2) Boxplot: per-sample mean SD across NMD genes
p_ave_sd <- ggplot(all_nmd_sample_summary, aes(x = reference, y = mean_sd_depth, fill = reference)) +
  geom_boxplot(outlier.shape = 21, outlier.size = 2) +
  stat_summary(fun = mean, geom = "point", color = "black", size = 2.2, shape = 23) +
  scale_fill_manual(values = ref_colors) +
  labs(x = NULL, y = "Mean SD(depth) across NMD genes") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none", aspect.ratio = asrt)

save_pdf(p_ave_sd, "ave_sd_all_nmd_genes.pdf", width = 4.5, height = 6)

# 3) Scatter: SD depth CHM13 vs GRCh38
p_scatter_sd_persample <- ggplot(compare_sd, aes(x = GRCh38, y = CHM13)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "red") +
  coord_cartesian(xlim = c(0, 400), ylim = c(0, 400)) +
  labs(
    x = "Per-sample SD depth (GRCh38)",
    y = "Per-sample SD depth (CHM13)",
    title = ""
  ) +
  theme_minimal(base_size = 13) + 
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)) +
  theme(panel.grid.major = element_line(color = "lightgray", size = 0.5), panel.grid.minor = element_line(color = "lightgray", size = 0.25))

save_pdf(p_scatter_sd_persample, "scatter_sd_chm13_vs_grch38_persample.pdf", width = 6, height = 6)

# 4) Scatter: SD depth CHM13 vs GRCh38 (zoomed)
p_scatter_sd_mean <- ggplot(gene_summary_sd, aes(x = mean_GRCh38, y = mean_CHM13)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "red") +
  geom_text_repel(data = . %>% filter(abs(mean_diff) > 2), aes(label = gene), size = 2.5, max.overlaps = 20) +
  coord_cartesian(xlim = c(5, 60), ylim = c(5, 60)) +
  labs(
    x = "Mean SD depth (GRCh38)",
    y = "Mean SD depth (CHM13)",
    title = ""
  ) +
  theme_minimal(base_size = 13) + 
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)) +
  theme(panel.grid.major = element_line(color = "lightgray", size = 0.5), panel.grid.minor = element_line(color = "lightgray", size = 0.25))


save_pdf(p_scatter_sd_mean, "scatter_sd_chm13_vs_grch38_mean.pdf", width = 6, height = 6)

# 5) Scatter: mean depth CHM13 vs GRCh38
p_scatter_mean_mean <- ggplot(gene_summary_mean, aes(x = mean_GRCh38, y = mean_CHM13)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "red") +
  geom_text_repel(data = . %>% filter(abs(mean_diff) > 1), aes(label = gene), size = 2.5, max.overlaps = 20) +
  coord_cartesian(xlim = c(25, 50), ylim = c(25, 50)) +
  labs(
    x = "Mean depth (GRCh38)",
    y = "Mean depth (CHM13)",
    title = ""
  ) +
  theme_minimal(base_size = 13) + 
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)) +
  theme(panel.grid.major = element_line(color = "lightgray", size = 0.5), panel.grid.minor = element_line(color = "lightgray", size = 0.25))

save_pdf(p_scatter_mean_mean, "scatter_mean_chm13_vs_grch38_mean.pdf", width = 6, height = 6)

p_scatter_mean <- ggplot(compare_mean, aes(x = GRCh38, y = CHM13)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "red") +
  labs(
    x = "Mean depth (GRCh38)",
    y = "Mean depth (CHM13)",
    title = "Per-gene mean depth: CHM13 vs GRCh38"
  ) +
  theme_minimal(base_size = 13)

save_pdf(p_scatter_mean, "scatter_mean_chm13_vs_grch38.pdf", width = 7, height = 7)

# 6) Top genes: higher mean coverage in CHM13
t2t_higher <- gene_summary_mean %>%
  filter(is.finite(mean_diff), mean_diff > 0) %>%
  slice_head(n = top_n) %>%
  select(gene)

plot_t2t <- compare_mean %>%
  semi_join(t2t_higher, by = "gene") %>%
  group_by(gene) %>%
  mutate(mean_diff_gene = mean(diff_mean, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(gene = fct_reorder(gene, mean_diff_gene, .desc = TRUE))

p_t2t <- ggplot(plot_t2t, aes(x = gene, y = diff_mean)) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", color = "black", size = 2.2, shape = 23) +
  labs(x = NULL, y = "Mean depth difference (CHM13 − GRCh38)") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_pdf(p_t2t, "t2t_genes_higher_mean_depth.pdf", width = 12, height = 5)

# 7) Top genes: higher mean coverage in GRCh38
hg38_higher <- gene_summary_mean %>%
  filter(is.finite(mean_diff), mean_diff < 0) %>%
  mutate(abs_diff = -mean_diff) %>%
  arrange(desc(abs_diff)) %>%
  slice_head(n = top_n) %>%
  select(gene)

plot_hg38 <- compare_mean %>%
  semi_join(hg38_higher, by = "gene") %>%
  mutate(diff_grch38 = -diff_mean) %>%
  group_by(gene) %>%
  mutate(mean_diff_gene = mean(diff_mean, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(gene = fct_reorder(gene, mean_diff_gene, .desc = FALSE))

p_hg38 <- ggplot(plot_hg38, aes(x = gene, y = diff_grch38)) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", color = "black", size = 2.2, shape = 23) +
  labs(x = NULL, y = "Mean depth difference (GRCh38 − CHM13)") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_pdf(p_hg38, "grch38_genes_higher_mean_depth.pdf", width = 12, height = 5)

message("Done. PDFs and TSVs written to: ", normalizePath(out_dir))

