#load libraries
library(ggplot2)
library(readr)
library(patchwork)
library(dplyr)
library(forcats)
library(scales)
library(ggh4x)

#set working dir
setwd("/Users/00104561/Library/CloudStorage/OneDrive-UWA/Research/Projects/AIM1_DATA-REANALYSIS/benchmark/variants/variant_stats/")

# -- data & basic recodes (yours) --
df <- as.data.frame(read_csv("./data/hprc-linear_variant_stats.csv")) %>%
  mutate(pass_pct = 100 * (Dysgu_PASS / Dysgu_SVs),
         DV_var    = DV_var    / 1e6,
         DV_SNPs   = DV_SNPs   / 1e6,
         DV_Indels = DV_Indels / 1e6,
         Dysgu_SVs = Dysgu_SVs / 1e6)

df$data   <- fct_recode(df$data, "Linear" = "linear", "HPRC" = "hprc")
df$genome <- fct_recode(df$genome, "CHM13" = "chm13", "GRCh38" = "grch38")

df_mendel <- as.data.frame(read_csv("./data/hprc-linear_mendel.csv")) %>%
  mutate(nLowQual = nLowQual / 1e6,
         Rate     = Rate * 100)

# -- ensure BOTH facetting vars have the same 2 levels and don't drop empty facets --
df$data        <- factor(df$data, levels = c("Linear","HPRC"))
df_mendel$Type <- factor(df_mendel$Type, levels = c("Linear","HPRC"))  # keeps 2 facet columns even if one is empty

common_bits <- list(
  theme_bw(base_size = 14),
  scale_fill_manual(values = c("#f5cccd","#abc3eb")),
  stat_summary(fun = mean, geom = "point", color = "black", size = 2.1, shape = 23),
  theme(legend.position = "none"),
  theme(axis.title.y = element_text(size = 10.5))
)

# ---- plots (unchanged except: facet drop=FALSE; make sure 'scales' is plural) ----
viol_plot <- ggplot(df_mendel, aes(Ref, Rate, fill = Ref)) +
  geom_boxplot(outlier.shape = 21, outlier.size = 2) +
  labs(y = "Mendelian Violation Rate (%)", x = NULL) +
  scale_y_continuous(labels = label_number(accuracy = 0.1)) +
  facet_wrap(~Type, nrow = 1, drop = FALSE) + common_bits

denovo_plot <- ggplot(df_mendel, aes(Ref, Denovo, fill = Ref)) +
  geom_boxplot(outlier.shape = 21, outlier.size = 2) +
  labs(y = "De Novo Variants", x = NULL) +
  scale_y_continuous(labels = label_number(accuracy = 0.1)) +
  facet_wrap(~Type, nrow = 1, drop = FALSE) + common_bits

lowq_plot <- ggplot(df_mendel, aes(Ref, nLowQual, fill = Ref)) +
  geom_boxplot(outlier.shape = 21, outlier.size = 2) +
  labs(y = "Low Quality Variants (Millions)", x = NULL) +
  scale_y_continuous(labels = label_number(accuracy = 0.1)) +
  facet_wrap(~Type, nrow = 1, drop = FALSE) + common_bits

dv_plot <- ggplot(df, aes(genome, DV_var, fill = genome)) +
  geom_boxplot(outlier.shape = 21, outlier.size = 2) +
  labs(y = "DeepVariant Calls (Millions)", x = NULL) +
  scale_y_continuous(labels = label_number(accuracy = 0.1)) +
  facet_wrap(~data, nrow = 1, drop = FALSE) + common_bits

dv_snps_plot <- ggplot(df, aes(genome, DV_SNPs, fill = genome)) +
  geom_boxplot(outlier.shape = 21, outlier.size = 2) +
  labs(y = "DeepVariant SNP Calls (Millions)", x = NULL) +
  scale_y_continuous(labels = label_number(accuracy = 0.1)) +
  facet_wrap(~data, nrow = 1, drop = FALSE) + common_bits

dv_indels_plot <- ggplot(df, aes(genome, DV_Indels, fill = genome)) +
  geom_boxplot(outlier.shape = 21, outlier.size = 2) +
  labs(y = "DeepVariant InDel Calls (Millions)", x = NULL) +
  scale_y_continuous(labels = label_number(accuracy = 0.1)) +
  facet_wrap(~data, nrow = 1, drop = FALSE) + common_bits

dv_tstv_plot <- ggplot(df, aes(genome, ts_tv, fill = genome)) +
  geom_boxplot(outlier.shape = 21, outlier.size = 2) +
  labs(y = "DeepVariant Ts/Tv", x = NULL) +
  scale_y_continuous(labels = label_number(accuracy = 0.1)) +
  facet_wrap(~data, nrow = 1, drop = FALSE) + common_bits

dysgu_svs_plot <- ggplot(df, aes(genome, Dysgu_SVs, fill = genome)) +
  geom_boxplot(outlier.shape = 21, outlier.size = 2) +
  labs(y = "Dysgu SV Calls (Millions)", x = NULL) +
  scale_y_continuous(labels = label_number(accuracy = 0.1)) +
  facet_wrap(~data, nrow = 1, drop = FALSE) + common_bits

dysgu_pass_plot <- ggplot(df, aes(genome, pass_pct, fill = genome)) +
  geom_boxplot(outlier.shape = 21, outlier.size = 2) +
  labs(y = "Proportion of Pass SVs (%)", x = NULL) +
  scale_y_continuous(labels = label_number(accuracy = 0.1)) +
  facet_wrap(~data, nrow = 1, drop = FALSE) + common_bits

# ---- assemble + fix panel sizes ----
stacked <- (dv_plot + dysgu_svs_plot + dysgu_pass_plot +
              viol_plot + denovo_plot + lowq_plot) +
  patchwork::plot_layout(ncol = 3, guides = "collect") &
  ggh4x::force_panelsizes(rows = unit(55, "mm"),
                          cols = unit(35, "mm")) &
  theme(
    plot.margin = margin(6, 6, 6, 6),
    panel.background = element_blank(),
    plot.background  = element_rect(fill = "transparent", colour = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# ---- SAVE (explicit size so nothing is clipped) ----

ggsave(
  filename = "./figures/hprc-linear_grch38-chm13_variants.pdf",
  plot     = stacked,
  device   = cairo_pdf,
  width    = 300,      # mm (increase if still tight)
  height   = 160,      # mm
  units    = "mm",
  bg       = "transparent",
  limitsize = FALSE
)


