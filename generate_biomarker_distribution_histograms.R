#!/usr/bin/env Rscript
# Histograms: raw biomarker concentrations (D0/D2/D4) and patient-level fold changes.
# Supports verification of manuscript Methods text on log-scale LMMs.

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(ggplot2)
  library(ggpubr)
})

args <- commandArgs(trailingOnly = FALSE)
f <- sub("^--file=", "", args[grep("^--file=", args)])
root <- if (length(f)) dirname(normalizePath(f)) else getwd()
setwd(root)
source(file.path(root, "glp_paths.R"))
source(file.path(root, "glp_figure_theme.R"))
source(file.path(root, "mixed_model_helpers.R"))

paths <- glp_ensure_dirs(glp_paths(root))
rev <- glp_ensure_dirs(glp_manuscript_revision_paths(root))
out_fig <- rev$figures_distribution
out_stats <- rev$stats

DAY_LEVELS <- c("D0", "D2", "D4")
DAY_LABELS <- c("D0" = "Day 0", "D2" = "Day 2", "D4" = "Day 4")
ASSAYS <- c("Active GLP-1", "C-Peptide", "Glucagon")

# Minimal theme: no in-figure titles/captions (captions saved separately for Word).
theme_hist_clean <- function() {
  theme_bw(base_size = 11) +
    theme(
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.caption = element_blank(),
      strip.text = element_text(size = 10, face = "bold"),
      strip.background = element_rect(fill = "grey92", colour = NA),
      panel.grid.minor = element_blank()
    )
}

skewness_simple <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 3) return(NA_real_)
  m <- mean(x)
  s <- sd(x)
  if (s == 0) return(0)
  mean(((x - m) / s)^3)
}

outcomes <- fread("cafpint.bms.outcomes.glp.csv") %>% as.data.frame()

assay_long <- outcomes %>%
  filter(Assay %in% ASSAYS) %>%
  transmute(
    cafpintid,
    biomarker = Assay,
    Day = factor(Day, levels = DAY_LEVELS),
    value_pM = Calc..Concentration
  ) %>%
  filter(is.finite(value_pM), value_pM > 0)

il6_long <- outcomes %>%
  distinct(cafpintid, Day, .keep_all = TRUE) %>%
  transmute(
    cafpintid,
    biomarker = "IL-6",
    Day = factor(Day, levels = DAY_LEVELS),
    value_pM = il6
  ) %>%
  filter(is.finite(value_pM), value_pM > 0)

biomarker_long <- bind_rows(assay_long, il6_long) %>%
  mutate(
    biomarker = factor(biomarker, levels = c(ASSAYS, "IL-6")),
    ln_pM = log(value_pM),
    log10_1p = log10(1 + value_pM)
  )

# Per-biomarker display limits (99th percentile) so Active GLP-1 is not crushed on one bar.
biomarker_limits <- biomarker_long %>%
  group_by(biomarker) %>%
  summarise(
    x_max = if (unique(biomarker) == "Active GLP-1") {
      as.numeric(quantile(value_pM, 0.95, na.rm = TRUE))
    } else {
      as.numeric(quantile(value_pM, 0.99, na.rm = TRUE))
    },
    binwidth = if (unique(biomarker) == "Active GLP-1") {
      0.05
    } else {
      as.numeric(quantile(value_pM, 0.99, na.rm = TRUE)) / 25
    },
    .groups = "drop"
  )

biomarker_plot <- biomarker_long %>%
  inner_join(biomarker_limits, by = "biomarker") %>%
  filter(value_pM <= x_max)

# --- Skewness table ---
skew_tbl <- biomarker_long %>%
  group_by(biomarker) %>%
  summarise(
    n = n(),
    skew_raw_pM = skewness_simple(value_pM),
    skew_ln_pM = skewness_simple(ln_pM),
    skew_log10_1p = skewness_simple(log10_1p),
    median_raw = median(value_pM),
    p95_raw = quantile(value_pM, 0.95),
    p99_raw = quantile(value_pM, 0.99),
    max_raw = max(value_pM),
    .groups = "drop"
  )
fwrite(skew_tbl, file.path(out_stats, "biomarker_skewness_by_scale.csv"))
fwrite(biomarker_limits, file.path(out_stats, "biomarker_histogram_axis_limits.csv"))

# --- Patient-level fold changes (D2/D0, D4/D0) ---
fold_df <- biomarker_long %>%
  select(cafpintid, biomarker, Day, value_pM) %>%
  pivot_wider(names_from = Day, values_from = value_pM) %>%
  filter(!is.na(D0), D0 > 0) %>%
  mutate(
    fold_D2_vs_D0 = D2 / D0,
    fold_D4_vs_D0 = D4 / D0,
    log2_fc_D2 = log2(1 + D2 / D0),
    log2_fc_D4 = log2(1 + D4 / D0)
  ) %>%
  pivot_longer(
    c(fold_D2_vs_D0, fold_D4_vs_D0),
    names_to = "contrast",
    values_to = "fold_vs_D0"
  ) %>%
  mutate(
    contrast = dplyr::recode(contrast,
      "fold_D2_vs_D0" = "D2 / D0",
      "fold_D4_vs_D0" = "D4 / D0"
    ),
    log2_fc = if_else(contrast == "D2 / D0", log2_fc_D2, log2_fc_D4)
  )

fold_summary <- fold_df %>%
  filter(is.finite(fold_vs_D0)) %>%
  group_by(biomarker, contrast) %>%
  summarise(
    n = n(),
    median_fold = median(fold_vs_D0),
    q25_fold = quantile(fold_vs_D0, 0.25),
    q75_fold = quantile(fold_vs_D0, 0.75),
    pct_fold_gt_1 = 100 * mean(fold_vs_D0 > 1),
    pct_fold_lt_1 = 100 * mean(fold_vs_D0 < 1),
    .groups = "drop"
  )
fwrite(fold_summary, file.path(out_stats, "biomarker_fold_change_summary.csv"))

# For transformed scales, compute binwidth per biomarker.
transform_limits <- biomarker_long %>%
  group_by(biomarker) %>%
  summarise(
    binwidth = max(diff(range(log10_1p, na.rm = TRUE)) / 25, 0.02),
    .groups = "drop"
  )

ln_limits <- biomarker_long %>%
  group_by(biomarker) %>%
  summarise(
    binwidth = max(diff(range(ln_pM, na.rm = TRUE)) / 25, 0.05),
    .groups = "drop"
  )

units_lab <- c(
  "Active GLP-1" = "Active GLP-1 (pM)",
  "C-Peptide" = "C-Peptide (pg/mL)",
  "Glucagon" = "Glucagon (pM)",
  "IL-6" = "IL-6 (pg/mL)"
)

hist_raw_xlab <- function(biomarker) {
  switch(
    biomarker,
    "Active GLP-1" = "Concentration (pM)",
    "C-Peptide" = "Concentration (pg/mL)",
    "Glucagon" = "Concentration (pM)",
    "IL-6" = "Concentration (pg/mL)",
    "Concentration"
  )
}

hist_log10_xlab <- function(biomarker) {
  paste0("log10(1 + ", assay_concentration_label(biomarker), ")")
}

plot_hist_by_biomarker <- function(data, x_var, xlab, fill = "#8DA0CB", use_binwidth_col = TRUE) {
  plots <- data %>%
    group_by(biomarker) %>%
    group_split(.keep = TRUE) %>%
    lapply(function(d) {
      biomarker_chr <- as.character(d$biomarker[1])
      x_vals <- d[[x_var]]
      bw <- if (use_binwidth_col && "binwidth" %in% names(d)) {
        unique(d$binwidth)[1]
      } else {
        max(diff(range(x_vals, na.rm = TRUE)) / 25, 0.02)
      }
      x_label <- if (is.function(xlab)) xlab(biomarker_chr) else xlab
      ggplot(d, aes(x = .data[[x_var]])) +
        geom_histogram(
          binwidth = bw,
          boundary = 0,
          fill = fill,
          colour = "white",
          linewidth = 0.25
        ) +
        facet_wrap(~Day, ncol = 3, labeller = labeller(Day = DAY_LABELS)) +
        labs(
          title = units_lab[biomarker_chr],
          x = x_label,
          y = "Count"
        ) +
        theme_hist_clean() +
        theme(plot.title = element_text(size = 10, face = "bold", hjust = 0))
    })
  ggarrange(plotlist = plots, ncol = 1, nrow = length(plots))
}

p_raw <- plot_hist_by_biomarker(
  biomarker_plot, "value_pM", hist_raw_xlab
)
save_glp_fig(p_raw, file.path(out_fig, "histogram_raw_by_day"), w = 12, h = 14)

p_log10 <- biomarker_long %>%
  inner_join(transform_limits, by = "biomarker") %>%
  { plot_hist_by_biomarker(., "log10_1p", hist_log10_xlab, fill = "#E6A800") }
save_glp_fig(p_log10, file.path(out_fig, "histogram_log10_1p_by_day"), w = 12, h = 14)

p_ln <- biomarker_long %>%
  inner_join(ln_limits, by = "biomarker") %>%
  plot_hist_by_biomarker("ln_pM", "ln(concentration)", fill = "#66C2A5")
save_glp_fig(p_ln, file.path(out_fig, "histogram_ln_by_day"), w = 12, h = 14)

# Active GLP-1 zoomed panel (full cohort, raw scale 0–99th pct)
glp_hi <- as.numeric(quantile(
  biomarker_long$value_pM[biomarker_long$biomarker == "Active GLP-1"],
  0.95, na.rm = TRUE
))
p_glp_zoom <- ggplot(
  biomarker_long %>%
    filter(biomarker == "Active GLP-1", value_pM <= glp_hi),
  aes(x = value_pM)
) +
  geom_histogram(binwidth = 0.05, boundary = 0, fill = "#8DA0CB", colour = "white") +
  facet_wrap(~Day, labeller = labeller(Day = DAY_LABELS)) +
  coord_cartesian(xlim = c(0, glp_hi * 1.05)) +
  labs(x = "Active GLP-1 (pM)", y = "Count") +
  theme_hist_clean()
save_glp_fig(p_glp_zoom, file.path(out_fig, "histogram_active_glp1_zoom"), w = 10, h = 4)

# Fold-change histograms (per biomarker x contrast panel)
fold_plot_df <- fold_df %>%
  filter(is.finite(fold_vs_D0), fold_vs_D0 > 0) %>%
  group_by(biomarker, contrast) %>%
  filter(fold_vs_D0 <= quantile(fold_vs_D0, 0.99, na.rm = TRUE)) %>%
  ungroup()

p_fold <- fold_plot_df %>%
  group_by(biomarker, contrast) %>%
  group_split() %>%
  lapply(function(d) {
    bw <- max(diff(range(d$fold_vs_D0, na.rm = TRUE)) / 20, 0.05)
    ggplot(d, aes(x = fold_vs_D0)) +
      geom_histogram(binwidth = bw, boundary = 0, fill = "#FC8D62", colour = "white") +
      geom_vline(xintercept = 1, linetype = "dashed", colour = "grey30") +
      labs(
        title = paste(d$biomarker[1], d$contrast[1]),
        x = "Fold change vs Day 0",
        y = "Count"
      ) +
      theme_hist_clean() +
      theme(plot.title = element_text(size = 9, face = "bold"))
  }) %>%
  ggarrange(plotlist = ., ncol = 2, nrow = 4, common.legend = FALSE)

save_glp_fig(p_fold, file.path(out_fig, "histogram_fold_change_vs_D0"), w = 11, h = 14)

p_log2fc <- fold_df %>%
  filter(is.finite(log2_fc)) %>%
  group_by(biomarker, contrast) %>%
  group_split() %>%
  lapply(function(d) {
    bw <- max(diff(range(d$log2_fc, na.rm = TRUE)) / 20, 0.05)
    ggplot(d, aes(x = log2_fc)) +
      geom_histogram(binwidth = bw, fill = "#A6D854", colour = "white") +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "grey30") +
      labs(
        title = paste(d$biomarker[1], d$contrast[1]),
        x = "log2(1 + Dx / D0)",
        y = "Count"
      ) +
      theme_hist_clean() +
      theme(plot.title = element_text(size = 9, face = "bold"))
  }) %>%
  ggarrange(plotlist = ., ncol = 2, nrow = 4)

save_glp_fig(p_log2fc, file.path(out_fig, "histogram_log2_fold_change_vs_D0"), w = 11, h = 14)

p_combo <- ggarrange(p_raw, p_fold, ncol = 1, nrow = 2, heights = c(1.15, 1))
save_glp_fig(p_combo, file.path(out_fig, "histogram_raw_and_fold_combined"), w = 12, h = 22)

# --- Figure captions for Word (not embedded in images) ---
captions <- c(
  "FIGURE CAPTIONS — Biomarker distribution checks (CAFPINT GLP cohort, n = 84)",
  "",
  "Figure S1. Raw biomarker concentrations by study day.",
  "Histograms of untransformed assay values at Day 0, Day 2, and Day 4. X-axes show assay-specific",
  "units: active GLP-1 and glucagon (pM); C-peptide and IL-6 (pg/mL). Values above the 99th",
  "percentile within each biomarker were excluded from display so that low-magnitude analytes",
  "(e.g., Active GLP-1) remain visible. Distributions are right-skewed for Active GLP-1,",
  "glucagon, and IL-6, motivating log-scale longitudinal modeling.",
  "",
  "Figure S2. Active GLP-1 concentrations (zoomed raw scale).",
  "Histograms of Active GLP-1 at each study day using a fixed bin width of 0.05 pM, with the",
  "display limited to the 95th percentile of observed values (~2 pM). Outliers above this range",
  "are omitted from display only. This panel shows the spread that is compressed when Active",
  "GLP-1 is shown on the same axis span as higher-magnitude biomarkers.",
  "",
  "Figure S3. log10(1 + concentration) by study day.",
  "Histograms of the log10(1 + concentration) transform applied in the primary linear mixed-effects",
  "models, with assay-specific units on each panel x-axis (C-peptide, pg/mL; active GLP-1 and",
  "glucagon, pM; IL-6, pg/mL). This is the scale used in the analysis code for longitudinal",
  "trajectory modeling.",
  "",
  "Figure S4. Natural-log concentration by study day.",
  "Histograms of ln(concentration) shown for comparison with natural-log modeling assumptions.",
  "",
  "Figure S5. Patient-level fold change versus Day 0.",
  "Histograms of the ratio of each patient's concentration at Day 2 or Day 4 to their Day 0 value",
  "(same biomarker). The dashed vertical line indicates no change (fold = 1). Values above the",
  "99th percentile within each biomarker and contrast were excluded for display. Active GLP-1,",
  "glucagon, and IL-6 fold changes remain right-skewed.",
  "",
  "Figure S6. Patient-level log2 fold change versus Day 0.",
  "Histograms of log2(1 + Dx/D0) for each biomarker at Day 2 and Day 4 relative to Day 0.",
  "The dashed vertical line indicates no change (log2 fold = 0).",
  "",
  "Figure S7. Combined distribution summary.",
  "Top: raw concentrations by biomarker and day (Figure S1). Bottom: patient-level fold changes",
  "versus Day 0 (Figure S5).",
  "",
  "FILE KEY (for matching exported filenames to captions):",
  "  histogram_raw_by_day.*          -> Figure S1",
  "  histogram_active_glp1_zoom.*  -> Figure S2",
  "  histogram_log10_1p_by_day.*   -> Figure S3",
  "  histogram_ln_by_day.*           -> Figure S4",
  "  histogram_fold_change_vs_D0.*   -> Figure S5",
  "  histogram_log2_fold_change_vs_D0.* -> Figure S6",
  "  histogram_raw_and_fold_combined.*  -> Figure S7"
)
writeLines(captions, file.path(out_stats, "figure_captions_distribution_histograms.txt"))

# Methods verification note (unchanged substance)
note <- c(
  "# Methods paragraph verification: log-scale LMM and fold-differences",
  "",
  "## Why Active GLP-1 looked like a single bar (fixed in this script)",
  "- Active GLP-1 concentrations are ~0.05–2 pM (median ~0.22 pM).",
  "- C-peptide and IL-6 range into thousands; wide automatic binning collapses GLP-1 into one bin.",
  "- Fix: per-biomarker 99th-percentile display limit + biomarker-specific bin width;",
  "  plus a dedicated zoomed Active GLP-1 panel (`histogram_active_glp1_zoom`).",
  "",
  "## Scale used in LMM code",
  "- Analysis code uses **log10(1 + concentration)**, not natural log.",
  "- See `biomarker_skewness_by_scale.csv` for skewness on raw, ln, and log10(1+x) scales."
)
writeLines(note, file.path(out_stats, "methods_log_lmm_verification.md"))

message("Saved histograms to: ", out_fig)
message("Saved figure captions: ", file.path(out_stats, "figure_captions_distribution_histograms.txt"))
sync_fixed_axis_figures(root)
message("Synced supplemental histograms to: ", rev$figures_fixed_axis)
