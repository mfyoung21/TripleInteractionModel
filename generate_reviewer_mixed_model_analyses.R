#!/usr/bin/env Rscript
# Reviewer-response analyses: longitudinal mixed models with 95% CIs.
# Produces model-based main figures + descriptive supplements + stats tables.
#
# Addresses:
#  - Figs 1-3: LMM trajectories (primary) vs Wilcoxon boxplots (descriptive)
#  - Fig 2: mortality x time interaction
#  - Figs 3 & 6: baseline-adjusted TGC models; hyper subgroup labeled
#  - Fig 5: proportional-change LMM with TGC x time within LCA
#  - Fig 4: exploratory mediation documentation
#  - Sensitivity: joint model for informative dropout (death)

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(ggplot2)
  library(ggpubr)
  library(patchwork)
  library(rstatix)
  library(mediation)
  library(quantreg)
  library(JM)
  library(survival)
  if (requireNamespace("survminer", quietly = TRUE)) library(survminer)
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
out_fig <- rev$figures_analysis
out_fig_ms <- rev$figures_manuscript
out_stats <- rev$stats
out_docs <- rev$docs

save_ms_fig <- function(p, stem, w, h, dpi = GLP_PUBLICATION_DPI) {
  save_reviewer_fig(p, file.path(out_fig_ms, stem), w = w, h = h, dpi = dpi)
}

ASSAYS <- c("Active GLP-1", "C-Peptide", "Glucagon")

# --- Build analysis dataset (same as regenerate_all_figures.R) ---
outcomes <- fread("cafpint.bms.outcomes.glp.csv") %>% as.data.frame()
cohort_v8 <- fread("GLP_regenerated/cafpintv8_glp_cohort_n84.csv") %>% as.data.frame()

gcg.df.med <- outcomes %>%
  dplyr::select(
    Sample, Assay, Calc..Concentration, Day, cafpintid,
    lcaprob_1011_dichot25, death90, insulin, il6,
    tgc = tgc.x, avginsdose, twbga
  ) %>%
  mutate(
    lca = factor(lcaprob_1011_dichot25, levels = c(0, 1), ordered = TRUE,
                 labels = c("Hypoinflamed", "Hyperinflamed")),
    tgc = factor(tgc, levels = c(1, 2), ordered = TRUE,
                 labels = c("TGC80-110", "TGC150-180")),
    Day = factor(Day, levels = DAY_LEVELS, ordered = TRUE),
    log10c = log10(1 + Calc..Concentration),
    death90 = as.integer(death90),
    death_f = factor(death90, levels = c(0, 1), labels = c("Survivor", "Nonsurvivor"))
  )

surv_meta <- cohort_v8 %>%
  distinct(cafpintid, .keep_all = TRUE) %>%
  transmute(
    cafpintid,
    death90_pt = as.integer(death90),
    deathday = pmax(as.numeric(deathday), 0.5, na.rm = TRUE)
  )

all_interactions <- list()
all_main_effects <- list()
all_emmeans <- list()
all_contrasts <- list()

record_model <- function(
  fig,
  analysis,
  assay,
  cohort,
  formula,
  fit,
  interactions,
  emm,
  contrasts = NULL,
  main_effects = NULL
) {
  if (!is.null(fit$error)) {
    message("[", fig, "] ", assay, " (", cohort, "): ", fit$error)
    return(invisible(NULL))
  }
  anova_tbl <- interactions %>%
    mutate(figure = fig, analysis = analysis, assay = assay, cohort = cohort, formula = deparse(formula))
  emm_tbl <- emm %>%
    mutate(figure = fig, analysis = analysis, assay = assay, cohort = cohort, formula = deparse(formula))
  all_interactions[[length(all_interactions) + 1]] <<- anova_tbl
  all_emmeans[[length(all_emmeans) + 1]] <<- emm_tbl
  if (!is.null(main_effects) && nrow(main_effects)) {
    main_effects <- main_effects %>%
      mutate(figure = fig, analysis = analysis, assay = assay, cohort = cohort, formula = deparse(formula))
    all_main_effects[[length(all_main_effects) + 1]] <<- main_effects
  }
  if (!is.null(contrasts) && nrow(contrasts)) {
    contrasts <- contrasts %>%
      mutate(figure = fig, analysis = analysis, assay = assay, cohort = cohort)
    all_contrasts[[length(all_contrasts) + 1]] <<- contrasts
  }
  invisible(list(fit = fit, emm = emm_tbl))
}

primary_main_term <- function(formula_rhs) {
  terms <- unique(unlist(strsplit(formula_rhs, " * ", fixed = TRUE)))
  terms <- terms[terms %in% c("lca", "tgc", "death_f", "Day")]
  if (length(terms) >= 1) terms[[1]] else NULL
}

human_model_term <- function(term, term_labels = NULL) {
  if (!is.null(term_labels) && !is.null(term_labels[[term]])) {
    return(term_labels[[term]])
  }
  switch(
    term,
    death_f = "90-day mortality",
    lca = "subphenotype",
    tgc = "treatment arm",
    Day = "study day",
    term
  )
}

wrap_fig_caption <- function(text, width = 88) {
  paste(strwrap(text, width = width), collapse = "\n")
}

run_trajectory_by_assay <- function(data, group_var, formula_rhs, fig, analysis, cohort_label,
                                    ylab = "log10(1 + concentration)",
                                    color_values = NULL, facet_var = NULL,
                                    shared_caption = NULL,
                                    legend_title = NULL,
                                    term_labels = NULL,
                                    subtitle_readable = FALSE,
                                    title_style = c("full", "assay_only"),
                                    collect_legend = FALSE,
                                    caption_width = 88,
                                    shared_title = NULL) {
  title_style <- match.arg(title_style)
  plots <- list()
  n_assay <- length(ASSAYS)
  for (i in seq_along(ASSAYS)) {
    assay <- ASSAYS[[i]]
    d <- data %>% filter(Assay == assay)
    form <- as.formula(paste("log10c ~", formula_rhs, "+ (1 | cafpintid)"))
    fit <- fit_lmm(d, form)
    int_terms <- unique(unlist(strsplit(formula_rhs, " * ", fixed = TRUE)))
    int_terms <- int_terms[int_terms %in% c("lca", "tgc", "death_f", "Day")]
    if (length(int_terms) >= 2) {
      int_term <- paste0(int_terms[1], ":", int_terms[2])
    } else {
      int_term <- "Day"
    }
    interactions <- anova_interactions(fit, int_term)
    main_term <- primary_main_term(formula_rhs)
    main_effects <- if (!is.null(main_term)) anova_main_effects(fit, main_term) else NULL
    emm <- emmeans_trajectory(fit, as.formula(paste("~", group_var, "| Day")))
    record_model(fig, analysis, assay, cohort_label, form, fit, interactions, emm, main_effects = main_effects)

    title <- if (title_style == "assay_only") {
      assay
    } else {
      paste0(assay, ": ", analysis, " (", cohort_label, ")")
    }

    if (subtitle_readable) {
      main_label <- human_model_term(main_term, term_labels)
      int_label <- if (length(int_terms) >= 2) {
        paste(
          human_model_term(int_terms[[1]], term_labels),
          "x",
          human_model_term(int_terms[[2]], term_labels)
        )
      } else {
        "interaction"
      }
      subtitle <- paste0(
        main_label, " main effect p = ", signif(main_effects$p_value[1], 3),
        "; ", int_label, " interaction p = ", signif(interactions$p_value[1], 3)
      )
    } else {
      subtitle <- NULL
    }

    plots[[assay]] <- plot_trajectory_lmm(
      emm, x_var = "Day", color_var = group_var,
      facet_var = facet_var, ylab = resolve_assay_ylab(assay, ylab), title = title, subtitle = subtitle,
      color_values = color_values,
      legend_title = legend_title,
      show_legend = !collect_legend || i == n_assay,
      xlab = if (collect_legend && i < n_assay) NULL else "Study day"
    )
  }

  if (collect_legend) {
    p_panels <- patchwork::wrap_plots(plots, ncol = 1)
    p <- p_panels +
      plot_layout(guides = "collect") &
      theme(legend.position = "bottom", legend.box = "horizontal")
    if (!is.null(shared_title)) {
      p <- p + patchwork::plot_annotation(
        title = shared_title,
        theme = theme(
          plot.title = element_text(face = "bold", hjust = 0.5, size = 15, margin = margin(b = 6))
        )
      )
    }
  } else {
    p <- ggarrange(plotlist = plots, ncol = 1, nrow = 3, align = "v")
    if (!is.null(shared_title)) {
      p <- p + patchwork::plot_annotation(
        title = shared_title,
        theme = theme(
          plot.title = element_text(face = "bold", hjust = 0.5, size = 15, margin = margin(b = 6))
        )
      )
    }
  }
  p
}

# ============================================================================
# FIGURES 1-3: LMM trajectories (primary) + descriptive boxplots (supplement)
# ============================================================================
message("=== Figures 1-3: LCA, mortality, TGC ===")

p_fig01_lmm <- run_trajectory_by_assay(
  gcg.df.med, "lca", "lca * Day",
  fig = "Fig1", analysis = "LCA x time LMM", cohort_label = "Full cohort (n=84)",
  color_values = GLP_LCA_COLORS
)
save_reviewer_fig(p_fig01_lmm, file.path(out_fig, "fig01_lca_trajectory_lmm"), w = 14, h = 16)

manuscript_lca_plots <- lapply(ASSAYS, function(assay) {
  d <- gcg.df.med %>% filter(Assay == assay)
  form <- log10c ~ lca * Day + (1 | cafpintid)
  fit <- fit_lmm(d, form)
  emm <- emmeans_trajectory(fit, ~ lca | Day)
  main_p <- anova_main_effects(fit, "lca")$p_value[1]
  int_p <- anova_interactions(fit, "lca:Day")$p_value[1]
  ylab_txt <- assay_log10_ylab(assay, "mean")
  plot_trajectory_lmm(
    emm,
    "Day",
    "lca",
    ylab = ylab_txt,
    title = assay,
    color_values = GLP_LCA_COLORS
  )
})
p_fig01_manuscript <- ggarrange(plotlist = manuscript_lca_plots, ncol = 1, nrow = 3, align = "v") +
  patchwork::plot_annotation(
    title = "Longitudinal trajectories by inflammatory subphenotype",
    theme = theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16)
    )
  )
save_ms_fig(p_fig01_manuscript, "fig01_manuscript_lca_pooled_lmm", w = 10, h = 14)

p_fig02_manuscript <- run_trajectory_by_assay(
  gcg.df.med, "death_f", "death_f * Day",
  fig = "Fig2", analysis = "Mortality x time LMM", cohort_label = "Full cohort (n=84)",
  ylab = "Model-based mean\nlog10(1 + concentration)",
  color_values = c("Survivor" = "#66C2A5", "Nonsurvivor" = "#FC8D62"),
  legend_title = "90-day mortality",
  subtitle_readable = TRUE,
  title_style = "assay_only",
  collect_legend = TRUE,
  shared_title = "Biomarker trajectories by 90-day mortality"
)
save_ms_fig(p_fig02_manuscript, "fig02_manuscript_mortality_lmm", w = 10, h = 13.5)

p_fig02_lmm <- run_trajectory_by_assay(
  gcg.df.med, "death_f", "death_f * Day",
  fig = "Fig2", analysis = "Mortality x time LMM", cohort_label = "Full cohort (n=84)",
  ylab = "Model-based mean\nlog10(1 + concentration)",
  color_values = c("Survivor" = "#66C2A5", "Nonsurvivor" = "#FC8D62"),
  legend_title = "90-day mortality",
  subtitle_readable = TRUE,
  title_style = "assay_only",
  collect_legend = TRUE
)
save_reviewer_fig(p_fig02_lmm, file.path(out_fig, "fig02_death90_trajectory_lmm"), w = 14, h = 16)

p_fig03_pooled_lmm <- run_trajectory_by_assay(
  gcg.df.med, "tgc", "tgc * Day",
  fig = "Fig3", analysis = "TGC x time LMM (unadjusted)", cohort_label = "Full cohort (n=84)",
  color_values = c("TGC80-110" = "#8DA0CB", "TGC150-180" = "#A6D854")
)
save_reviewer_fig(p_fig03_pooled_lmm, file.path(out_fig, "fig03_tgc_trajectory_lmm_pooled"), w = 14, h = 16)

# Fig 3 hyper subgroup (clarify pooled vs hyper)
gcg_hyper <- gcg.df.med %>% filter(lca == "Hyperinflamed")
p_fig03_hyper_lmm <- run_trajectory_by_assay(
  gcg_hyper, "tgc", "tgc * Day",
  fig = "Fig3", analysis = "TGC x time LMM (unadjusted)", cohort_label = "Hyperinflammatory subgroup (n=42)",
  color_values = c("TGC80-110" = "#8DA0CB", "TGC150-180" = "#A6D854")
)
save_reviewer_fig(p_fig03_hyper_lmm, file.path(out_fig, "fig03_tgc_trajectory_lmm_hyper"), w = 14, h = 16)

# Baseline-adjusted TGC models (post-baseline values)
gcg_adj <- gcg.df.med %>%
  group_by(Assay, cafpintid) %>%
  mutate(
    baseline = log10c[as.character(Day) == "D0"][1],
    change_from_baseline = log10c - baseline,
    post_baseline = as.character(Day) != "D0"
  ) %>%
  ungroup() %>%
  filter(post_baseline)

for (cohort_label in c("Full cohort (n=84)", "Hyperinflammatory subgroup (n=42)")) {
  d <- if (cohort_label == "Full cohort (n=84)") gcg_adj else gcg_adj %>% filter(lca == "Hyperinflamed")
  tag <- if (cohort_label == "Full cohort (n=84)") "pooled" else "hyper"
  plots <- list()
  for (assay in ASSAYS) {
    dd <- d %>% filter(Assay == assay)
    form <- log10c ~ tgc * Day + baseline + (1 | cafpintid)
    fit <- fit_lmm(dd, form)
    interactions <- anova_interactions(fit, "tgc:Day")
    main_effects <- anova_main_effects(fit, "tgc")
    emm <- emmeans_trajectory(fit, ~ tgc | Day)
    record_model("Fig3/6", "TGC x time LMM (baseline-adjusted)", assay, cohort_label, form, fit, interactions, emm, main_effects = main_effects)
    ctr <- compare_emmeans_contrasts(fit, ~ tgc | Day)
    record_model("Fig3/6", "TGC pairwise at each day (adjusted)", assay, cohort_label, form, fit, interactions, emm, ctr, main_effects = main_effects)
    plots[[assay]] <- plot_trajectory_lmm(
      emm, "Day", "tgc",
      ylab = assay_log10_ylab(assay, "baseline_adjusted"),
      title = paste0(assay, ": baseline-adjusted treatment (", cohort_label, ")"),
      color_values = c("TGC80-110" = "#8DA0CB", "TGC150-180" = "#A6D854")
    )
  }
  p_adj <- ggarrange(plotlist = plots, ncol = 1, nrow = 3, align = "v")
  save_reviewer_fig(
    p_adj,
    file.path(out_fig, paste0("fig03_tgc_baseline_adjusted_lmm_", tag)),
    w = 14, h = 16
  )
}

# Change-from-baseline TGC model (alternative reviewer suggestion)
gcg_chg <- gcg.df.med %>%
  group_by(Assay, cafpintid) %>%
  mutate(
    baseline = log10c[as.character(Day) == "D0"][1],
    change_from_baseline = log10c - baseline,
    post_baseline = as.character(Day) != "D0"
  ) %>%
  ungroup() %>%
  filter(post_baseline)

for (tag_info in list(
  list(data = gcg_chg, tag = "pooled", label = "Full cohort (n=84)"),
  list(data = gcg_chg %>% filter(lca == "Hyperinflamed"), tag = "hyper", label = "Hyperinflammatory subgroup (n=42)")
)) {
  plots <- list()
  for (assay in ASSAYS) {
    dd <- tag_info$data %>% filter(Assay == assay)
    form <- change_from_baseline ~ tgc * Day + (1 | cafpintid)
    fit <- fit_lmm(dd, form, y_var = "change_from_baseline")
    interactions <- anova_interactions(fit, "tgc:Day")
    emm <- emmeans_trajectory(fit, ~ tgc | Day, y_var = "change_from_baseline")
    record_model("Fig3/6", "Change from baseline LMM", assay, tag_info$label, form, fit, interactions, emm)
    plots[[assay]] <- plot_trajectory_lmm(
      emm, "Day", "tgc",
      ylab = assay_log10_ylab(assay, "change_from_baseline"),
      title = paste0(assay, ": change from baseline (", tag_info$label, ")"),
      color_values = c("TGC80-110" = "#8DA0CB", "TGC150-180" = "#A6D854")
    )
  }
  p_chg <- ggarrange(plotlist = plots, ncol = 1, nrow = 3, align = "v")
  save_reviewer_fig(
    p_chg,
    file.path(out_fig, paste0("fig03_tgc_change_from_baseline_lmm_", tag_info$tag)),
    w = 14, h = 16
  )
}

# Descriptive boxplots (supplement) — reuse logic from regenerate_all_figures.R
source_lines <- readLines("regenerate_all_figures.R")
# Minimal descriptive plots with clear labeling
descriptive_box <- function(data, fill_var, colors, labels, name, fname) {
  p <- ggplot(data, aes(x = Day, y = log10c, fill = .data[[fill_var]])) +
    geom_boxplot(position = position_dodge(0.75), outlier.shape = NA, alpha = 0.85) +
    geom_jitter(position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.75),
                alpha = 0.2, size = 0.5) +
    facet_wrap(~Assay, scales = "free_y") +
    scale_fill_manual(values = colors, labels = labels) +
    labs(
      title = name,
      x = "Study day",
      y = "log10(1 + concentration)",
      fill = name
    ) +
    theme_glp_fig()
  save_reviewer_fig(p, file.path(out_fig, fname), w = 14, h = 11)
}

descriptive_box(
  gcg.df.med, "lca",
  GLP_LCA_COLORS,
  GLP_LCA_LABELS, "Latent class",
  "fig01_lca_boxplot_descriptive"
)
descriptive_box(
  gcg.df.med, "death_f",
  c("Survivor" = "#66C2A5", "Nonsurvivor" = "#FC8D62"),
  c("Survivor", "Nonsurvivor"), "90-day mortality",
  "fig02_death90_boxplot_descriptive"
)
descriptive_box(
  gcg.df.med, "tgc",
  c("TGC80-110" = "#8DA0CB", "TGC150-180" = "#A6D854"),
  c("TGC80-110", "TGC150-180"), "Treatment arm",
  "fig03_tgc_boxplot_descriptive_pooled"
)
descriptive_box(
  gcg_hyper, "tgc",
  c("TGC80-110" = "#8DA0CB", "TGC150-180" = "#A6D854"),
  c("TGC80-110", "TGC150-180"), "Treatment arm",
  "fig03_tgc_boxplot_descriptive_hyper"
)

# ============================================================================
# FIGURE 5: Proportional trajectories (C-peptide, insulin, insdoseprop)
# ============================================================================
message("=== Figure 5: proportional-change LMM ===")

day_file <- if (file.exists("cafpintday_v4.csv")) "cafpintday_v4.csv" else stop("Need cafpintday_v4.csv")
cafpint.day <- fread(day_file) %>% as.data.frame()

assay_long <- outcomes %>%
  filter(Assay %in% c("C-Peptide", "Active GLP-1")) %>%
  transmute(cafpintid, Day = factor(Day, levels = DAY_LEVELS), name = Assay,
            value = Calc..Concentration, lcaprob_1011_dichot25, tgc = tgc.x)

insulin_long <- outcomes %>%
  distinct(cafpintid, Day, .keep_all = TRUE) %>%
  transmute(cafpintid, Day = factor(Day, levels = DAY_LEVELS), name = "insulin",
            value = insulin, lcaprob_1011_dichot25, tgc = tgc.x) %>%
  filter(!is.na(value))

gcg_long <- bind_rows(assay_long, insulin_long)

endo_ids <- gcg_long %>%
  filter(name %in% c("insulin", "C-Peptide")) %>%
  distinct(cafpintid, Day, name, .keep_all = TRUE) %>%
  pivot_wider(names_from = name, values_from = value) %>%
  filter(!is.na(`C-Peptide`), !is.na(insulin)) %>%
  pull(cafpintid) %>% unique()

gcg_prop <- gcg_long %>%
  filter(name %in% c("insulin", "C-Peptide"), cafpintid %in% endo_ids) %>%
  distinct(cafpintid, name, Day, .keep_all = TRUE) %>%
  pivot_wider(id_cols = c(cafpintid, name, lcaprob_1011_dichot25, tgc),
              names_from = Day, values_from = value) %>%
  filter(!is.na(D0), D0 > 0) %>%
  mutate(
    D2 = log10_fold_change(D0, D2),
    D4 = log10_fold_change(D0, D4),
    D0 = 0
  ) %>%
  pivot_longer(c(D0, D2, D4), names_to = "Day", values_to = "prop_value") %>%
  mutate(Day = factor(Day, levels = DAY_LEVELS))

caf.insdose <- cafpint.day %>%
  dplyr::select(cafpintid, day, insdose) %>%
  filter(day %in% c(0, 2, 4)) %>%
  pivot_wider(names_from = day, names_prefix = "insdose.", values_from = insdose)

ins_prop <- caf.insdose %>%
  filter(cafpintid %in% endo_ids) %>%
  filter(all(c("insdose.0", "insdose.2", "insdose.4") %in% names(.))) %>%
  mutate(
    insdoseprop.0 = 0,
    insdoseprop.2 = log10(insdose.2 + 1) - log10(insdose.0 + 1),
    insdoseprop.4 = log10(insdose.4 + 1) - log10(insdose.0 + 1)
  ) %>%
  dplyr::select(cafpintid, insdoseprop.0, insdoseprop.2, insdoseprop.4) %>%
  pivot_longer(-cafpintid, names_to = "day_col", values_to = "prop_value") %>%
  mutate(
    Day = factor(
      dplyr::recode(day_col, insdoseprop.0 = "D0", insdoseprop.2 = "D2", insdoseprop.4 = "D4"),
      levels = DAY_LEVELS
    ),
    name = "insdoseprop"
  ) %>%
  left_join(
    outcomes %>% distinct(cafpintid, .keep_all = TRUE) %>%
      transmute(cafpintid, lcaprob_1011_dichot25, tgc = tgc.x),
    by = "cafpintid"
  )

fig5_df <- bind_rows(gcg_prop, ins_prop) %>%
  filter(cafpintid %in% endo_ids, is.finite(prop_value)) %>%
  mutate(
    analyte = factor(
      name,
      levels = c("C-Peptide", "insulin", "insdoseprop"),
      labels = c("C-Peptide", "insulin", "insulin-dose proportion")
    ),
    lca = factor(lcaprob_1011_dichot25, levels = c(0, 1),
                 labels = c("Hypoinflamed", "Hyperinflamed")),
    tgc = factor(tgc, levels = c(1, 2), labels = c("TGC80-110", "TGC150-180")),
    Day = factor(as.character(Day), levels = DAY_LEVELS)
  )

fig5_plots <- list()
fig5_interactions <- list()
fig5_emm_all <- list()

for (an in levels(fig5_df$analyte)) {
  dd <- fig5_df %>% filter(analyte == an)
  for (lca_lvl in levels(dd$lca)) {
    sub <- dd %>% filter(lca == lca_lvl, as.character(Day) != "D0")
    form_sub <- prop_value ~ tgc * Day + (1 | cafpintid)
    fit_sub <- fit_lmm(sub, form_sub, y_var = "prop_value")
    int <- anova_interactions(fit_sub, "tgc:Day") %>%
      mutate(analyte = an, lca = lca_lvl)
    main <- anova_main_effects(fit_sub, "tgc") %>%
      mutate(analyte = an, lca = lca_lvl, term = "tgc (main)")
    fig5_interactions[[length(fig5_interactions) + 1]] <- bind_rows(int, main)
    emm <- emmeans_trajectory(fit_sub, ~ tgc | Day, y_var = "prop_value") %>%
      mutate(analyte = an, lca = lca_lvl)
    fig5_emm_all[[length(fig5_emm_all) + 1]] <- emm
  }
}

fig5_emm <- bind_rows(fig5_emm_all)
fig5_int <- bind_rows(fig5_interactions)

p_fig05_lmm <- ggplot(
  fig5_emm %>% filter(as.character(Day) != "D0"),
  aes(x = Day, y = estimate, color = tgc, group = tgc)
) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = tgc), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.2) +
  facet_grid(analyte ~ lca, scales = "free_y") +
  scale_color_manual(values = c("TGC80-110" = "#8DA0CB", "TGC150-180" = "#A6D854")) +
  scale_fill_manual(values = c("TGC80-110" = "#8DA0CB", "TGC150-180" = "#A6D854")) +
  scale_x_discrete(labels = DAY_LABELS[c("D2", "D4")]) +
  labs(
    title = "Log fold change from Day 0 — Hyperinflammatory group",
    x = "Study day",
    y = "Log fold change from Day 0",
    color = "Treatment arm",
    fill = "Treatment arm"
  ) +
  theme_glp_fig()

save_reviewer_fig(p_fig05_lmm, file.path(out_fig, "fig05_proportional_trajectory_lmm"), w = 12, h = 14)

# Descriptive median supplement (old Fig 5 style)
p_fig05_desc <- ggplot(
  fig5_df %>% filter(name %in% c("C-Peptide", "insulin", "insdoseprop")),
  aes(x = Day, y = prop_value, color = analyte, group = analyte)
) +
  stat_summary(fun = median, geom = "line", linewidth = 1.5) +
  stat_summary(fun = median, geom = "point") +
  facet_grid(lca ~ tgc) +
  labs(
    title = "Descriptive supplement (empirical medians, no CI)",
    y = "Proportional change from Day 0 (unitless)"
  ) +
  theme_glp_fig()
save_reviewer_fig(p_fig05_desc, file.path(out_fig, "fig05_proportional_trajectory_descriptive"), w = 12, h = 10)

# ============================================================================
# FIGURE 6: IL-6 & Active GLP-1 proportional change (hyper only, units labeled)
# ============================================================================
message("=== Figure 6: IL-6 / GLP-1 hyper subgroup ===")

gcg_long6 <- gcg.df.med %>%
  dplyr::select(cafpintid, Day, name = Assay, value = Calc..Concentration,
                lcaprob_1011_dichot25, tgc) %>%
  bind_rows(
    gcg.df.med %>%
      dplyr::select(cafpintid, Day, il6, insulin, lcaprob_1011_dichot25, tgc) %>%
      distinct(cafpintid, Day, .keep_all = TRUE) %>%
      pivot_longer(c(il6, insulin), names_to = "name", values_to = "value")
  )

gcg_prop6 <- gcg_long6 %>%
  pivot_wider(names_from = Day, values_from = value) %>%
  filter(!is.na(D0), D0 > 0) %>%
  mutate(
    baseline_log = log10(1 + D0),
    D0p = 0,
    D2p = log10_fold_change(D0, D2),
    D4p = log10_fold_change(D0, D4)
  ) %>%
  pivot_longer(c(D0p, D2p, D4p), names_to = "Day", values_to = "prop_value") %>%
  mutate(
    Day = factor(
      dplyr::recode(Day, D0p = "D0", D2p = "D2", D4p = "D4"),
      levels = DAY_LEVELS
    ),
    lca = factor(lcaprob_1011_dichot25, levels = c(0, 1),
                 labels = c("Hypoinflamed", "Hyperinflamed")),
    tgc = factor(tgc, levels = c("TGC80-110", "TGC150-180"))
  )

il6_glp1 <- gcg_prop6 %>%
  filter(name %in% c("il6", "Active GLP-1"), lca == "Hyperinflamed", is.finite(prop_value)) %>%
  mutate(
    marker = factor(name, levels = c("il6", "Active GLP-1"),
                    labels = c("IL-6", "Active GLP-1"))
  )

fig6_emm <- list()
fig6_int <- list()
for (mk in levels(il6_glp1$marker)) {
  dd <- il6_glp1 %>% filter(marker == mk, as.character(Day) != "D0")
  form_post <- prop_value ~ tgc * Day + (1 | cafpintid)
  fit_post <- fit_lmm(dd, form_post, y_var = "prop_value")
  fig6_int[[length(fig6_int) + 1]] <- bind_rows(
    anova_interactions(fit_post, "tgc:Day") %>%
      mutate(marker = mk, model = "post-baseline log fold change (D2/D4)"),
    anova_main_effects(fit_post, "tgc") %>%
      mutate(marker = mk, model = "post-baseline log fold change (D2/D4)", term = "tgc (main)")
  )
  fig6_emm[[length(fig6_emm) + 1]] <- emmeans_trajectory(fit_post, ~ tgc | Day, y_var = "prop_value") %>%
    mutate(marker = mk, model = "post-baseline log fold change (D2/D4)")

  dd_adj <- il6_glp1 %>% filter(marker == mk, as.character(Day) != "D0")
  form_adj <- prop_value ~ tgc * Day + baseline_log + (1 | cafpintid)
  fit_adj <- fit_lmm(dd_adj, form_adj, y_var = "prop_value")
  fig6_int[[length(fig6_int) + 1]] <- bind_rows(
    anova_interactions(fit_adj, "tgc:Day") %>%
      mutate(marker = mk, model = "baseline-adjusted log fold change (D2/D4)"),
    anova_main_effects(fit_adj, "tgc") %>%
      mutate(marker = mk, model = "baseline-adjusted log fold change (D2/D4)", term = "tgc (main)")
  )
  fig6_emm[[length(fig6_emm) + 1]] <- emmeans_trajectory(fit_adj, ~ tgc | Day, y_var = "prop_value") %>%
    mutate(marker = mk, model = "baseline-adjusted log fold change (D2/D4)")
}

fig6_emm_df <- bind_rows(fig6_emm)
fig6_int_df <- bind_rows(fig6_int)

fig6_plot_df <- fig6_emm_df %>%
  filter(model == "baseline-adjusted log fold change (D2/D4)", as.character(Day) != "D0")

p_fig06_lmm <- ggplot(
  fig6_plot_df,
  aes(x = Day, y = estimate, color = tgc, group = tgc)
) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = tgc), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.5) +
  facet_wrap(~marker, scales = "free_y", ncol = 1) +
  scale_color_manual(values = c("TGC80-110" = "#8DA0CB", "TGC150-180" = "#A6D854")) +
  scale_fill_manual(values = c("TGC80-110" = "#8DA0CB", "TGC150-180" = "#A6D854")) +
  scale_x_discrete(labels = DAY_LABELS[c("D2", "D4")]) +
  labs(
    title = "IL-6 and Active GLP-1 — hyperinflammatory subgroup only",
    x = "Study day",
    y = "Log fold change from Day 0",
    color = "Treatment arm",
    fill = "Treatment arm"
  ) +
  theme_glp_fig()
save_reviewer_fig(p_fig06_lmm, file.path(out_fig, "fig06_il6_glp1_hyper_baseline_adjusted_lmm"), w = 11, h = 10)

fig6_post_df <- fig6_emm_df %>%
  filter(model == "post-baseline log fold change (D2/D4)", as.character(Day) != "D0")

p_fig05_manuscript <- ggplot(fig6_post_df, aes(x = Day, y = estimate, color = tgc, group = tgc)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = tgc), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.5) +
  facet_wrap(~marker, scales = "free_y", ncol = 1) +
  scale_color_manual(values = c("TGC80-110" = "#8DA0CB", "TGC150-180" = "#A6D854")) +
  scale_fill_manual(values = c("TGC80-110" = "#8DA0CB", "TGC150-180" = "#A6D854")) +
  scale_x_discrete(labels = DAY_LABELS[c("D2", "D4")]) +
  labs(
    title = "Log fold change from Day 0 — hyperinflammatory subgroup",
    x = "Study day",
    y = "Log fold change from Day 0",
    color = "Treatment arm",
    fill = "Treatment arm"
  ) +
  theme_glp_fig()
save_ms_fig(p_fig05_manuscript, "fig05_manuscript_il6_glp1_hyper", w = 11, h = 10)

# Endo.prop (Figure 6 alternate panel) with units
endo.df <- gcg.df.med %>%
  filter(Assay %in% c("insulin", "C-Peptide")) %>%
  dplyr::select(-log10c) %>%
  pivot_wider(names_from = Assay, values_from = Calc..Concentration) %>%
  filter(!is.na(`C-Peptide`), !is.na(insulin), insulin > 0) %>%
  mutate(
    endo.prop = log10(`C-Peptide` / insulin),
    lca_lab = factor(lcaprob_1011_dichot25, levels = c(0, 1),
                     labels = c("Hypoinflamed", "Hyperinflamed")),
    tgc_lab = factor(tgc, levels = c("TGC80-110", "TGC150-180"))
  )

endo_emm <- list()
endo_int <- list()
for (lca_lvl in levels(endo.df$lca_lab)) {
  dd <- endo.df %>% filter(lca_lab == lca_lvl)
  form <- endo.prop ~ tgc_lab * Day + (1 | cafpintid)
  fit <- fit_lmm(dd, form, y_var = "endo.prop")
  endo_int[[length(endo_int) + 1]] <- anova_interactions(fit, "tgc_lab:Day") %>%
    mutate(lca = lca_lvl)
  endo_emm[[length(endo_emm) + 1]] <- emmeans_trajectory(fit, ~ tgc_lab | Day, y_var = "endo.prop") %>%
    mutate(lca = lca_lvl)
}
endo_emm_df <- bind_rows(endo_emm)

p_endo_lmm <- ggplot(endo_emm_df, aes(x = Day, y = estimate, color = tgc_lab, group = tgc_lab)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = tgc_lab), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.2) +
  facet_wrap(~lca) +
  scale_color_manual(values = c("TGC80-110" = "#8DA0CB", "TGC150-180" = "#A6D854")) +
  scale_fill_manual(values = c("TGC80-110" = "#8DA0CB", "TGC150-180" = "#A6D854")) +
  labs(
    title = "Endotyping proportion by treatment arm",
    x = "Study day",
    y = "log10(C-peptide / insulin)",
    color = "Treatment arm",
    fill = "Treatment arm"
  ) +
  theme_glp_fig()
save_reviewer_fig(p_endo_lmm, file.path(out_fig, "fig06_endo_prop_lmm"), w = 12, h = 8)

# ============================================================================
# FIGURE 4: Exploratory mediation documentation
# ============================================================================
message("=== Figure 4: exploratory mediation ===")

med_note <- c(
  "# Figure 4 — Exploratory mediation (hypothesis-generating)",
  "",
  "## Limitation (reviewer comment)",
  "Subphenotype (LCA) is defined at Day 0, while the pooled mediator uses mean GLP-1 across D0/D2/D4.",
  "This does not establish clean temporal separation between exposure and mediator.",
  "",
  "## Recommendation",
  "Label as **exploratory / hypothesis-generating** in the manuscript.",
  "Primary causal inference should not rely on this analysis.",
  "",
  "## Implemented sensitivity",
  "See `crosslagged_mediation/run_crosslagged_mediation.R` for visit-specific mediators (D0, D2, D4).",
  "",
  "## Pooled analysis reproduced for documentation"
)
writeLines(med_note, file.path(out_docs, "fig04_mediation_exploratory_note.md"))

med_df <- gcg.df.med %>%
  filter(Assay == "Active GLP-1", !is.na(log10c), !is.na(lca), !is.na(death90))
if (nrow(med_df) >= 30) {
  set.seed(42)
  m_pool <- mediate(
    rq(log10c ~ lca, data = med_df),
    glm(death90 ~ log10c + lca, family = binomial(), data = med_df),
    treat = "lca", mediator = "log10c", boot = TRUE, sims = 500
  )
  med_sum <- tibble(
    analysis = "Pooled visits (exploratory)",
    exposure = "LCA (hyper vs hypo, D0 assignment)",
    mediator = "log10(1+GLP-1) pooled across D0/D2/D4",
    outcome = "death90",
    ACME = m_pool$d0,
    ACME_lo = m_pool$d0.ci[1],
    ACME_hi = m_pool$d0.ci[2],
    ADE = m_pool$z0,
    ADE_lo = m_pool$z0.ci[1],
    ADE_hi = m_pool$z0.ci[2],
    TE = m_pool$tau.coef,
    TE_lo = m_pool$tau.ci[1],
    TE_hi = m_pool$tau.ci[2],
    label = "EXPLORATORY — not primary causal evidence"
  )
  fwrite(med_sum, file.path(out_stats, "fig04_mediation_exploratory_summary.csv"))
  save_med_plot(
    m_pool,
    file.path(out_fig, "fig04_mediation_exploratory"),
    "EXPLORATORY: Active GLP-1 mediation of LCA-associated mortality (pooled visits)"
  )
}

# ============================================================================
# Sensitivity: Joint model for informative dropout (death)
# ============================================================================
message("=== Sensitivity: joint model (Active GLP-1) ===")

jm_results <- list()
glp_long <- gcg.df.med %>%
  filter(Assay == "Active GLP-1") %>%
  left_join(surv_meta, by = "cafpintid") %>%
  filter(is.finite(log10c)) %>%
  mutate(
    day_num = as.numeric(factor(Day, levels = DAY_LEVELS)),
    death90_pt = coalesce(death90_pt, death90)
  )

if (n_distinct(glp_long$cafpintid) >= 20) {
  # Shared random intercept joint model: longitudinal GLP-1 + survival
  lme_fit <- lme(
    log10c ~ day_num + lca + tgc,
    random = ~ 1 | cafpintid,
    data = glp_long,
    na.action = na.omit,
    control = lmeControl(opt = "optim")
  )
  surv_df <- glp_long %>%
    distinct(cafpintid, .keep_all = TRUE) %>%
    mutate(
      event = death90_pt,
      surv_time = pmax(deathday, 1)
    )
  cox_fit <- coxph(Surv(surv_time, event) ~ lca + tgc, data = surv_df, x = TRUE)
  jm_fit <- tryCatch(
    jointModel(lme_fit, cox_fit, timeVar = "day_num", method = "weibull-PH-aGH"),
    error = function(e) e
  )
  if (inherits(jm_fit, "jointModel")) {
    jm_sum <- tibble(
      analysis = "Joint model sensitivity (Active GLP-1)",
      longitudinal = "log10(1+GLP-1) ~ day + LCA + TGC + random intercept",
      survival = "Cox PH: death90 ~ LCA + TGC (time = deathday)",
      association = "Weibull-PH-aGH",
      loglik = as.numeric(logLik(jm_fit)),
      note = "Sensitivity for informative dropout due to early death; not primary analysis."
    )
    fwrite(jm_sum, file.path(out_stats, "sensitivity_joint_model_active_glp1.csv"))
    jm_results <- jm_sum
  } else {
    fwrite(
      tibble(
        analysis = "Joint model",
        status = "failed",
        error = conditionMessage(jm_fit)
      ),
      file.path(out_stats, "sensitivity_joint_model_active_glp1.csv")
    )
  }
}

# ============================================================================
# Mortality bar chart + Kaplan-Meier (full CAF-PINT cohort, n = 293)
# ============================================================================
message("=== Mortality bar chart + Kaplan-Meier (n = 293) ===")

km_master_path <- "../CAFPINT Caitlin/causal_mediation_analyses/data/mediation_master.csv"
km_alluv_path <- "_zip_ref/cafpint_insulin/CAFPINT insulin/cafpint.lca.lngr.dropout.tgc.alluv.csv"

if (file.exists(km_master_path)) {
  km_alluv <- if (file.exists(km_alluv_path)) {
    fread(km_alluv_path) %>%
      as.data.frame() %>%
      distinct(cafpintid, .keep_all = TRUE) %>%
      transmute(cafpintid, deathday_alluv = deathday)
  } else {
    tibble(cafpintid = character(), deathday_alluv = numeric())
  }

  km_df <- fread(km_master_path) %>%
    as.data.frame() %>%
    rename(cafpintid = Individual) %>%
    left_join(km_alluv, by = "cafpintid") %>%
    filter(!is.na(lcaprob_1011_dichot25)) %>%
    mutate(
      lca_grp = factor(
        lcaprob_1011_dichot25,
        levels = c(0, 1),
        labels = c("Hypoinflammatory", "Hyperinflammatory")
      ),
      tgc_grp = factor(
        tgc,
        levels = c(1, 2),
        labels = c(GLP_TGC_IIT_LABEL, GLP_TGC_CIT_LABEL)
      ),
      event = as.integer(death90),
      time = if_else(
        event == 1L,
        pmin(pmax(coalesce(deathday_alluv, 90), 0.5), 90),
        90
      )
    )

  km_bar_df <- km_df %>%
    group_by(lca_grp) %>%
    summarise(
      n = n(),
      deaths = sum(event),
      pct = 100 * deaths / n,
      .groups = "drop"
    )
  fisher_p <- fisher.test(table(km_df$lca_grp, km_df$event))$p.value

  p_mort_bar <- ggplot(km_bar_df, aes(x = lca_grp, y = pct, fill = lca_grp)) +
    geom_col(width = 0.62, colour = "black", linewidth = 0.4) +
    geom_text(
      aes(label = sprintf("%d / %d\n(%.1f%%)", deaths, n, pct)),
      vjust = -0.35,
      size = 4.2,
      fontface = "bold"
    ) +
    scale_fill_manual(values = GLP_SUBPHENOTYPE_COLORS, guide = "none") +
    scale_y_continuous(
      limits = c(0, max(km_bar_df$pct) * 1.28),
      breaks = seq(0, 40, 10),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(
      title = "90-day mortality by baseline inflammatory subphenotype",
      subtitle = sprintf(
        "Full CAF-PINT cohort (n = %d); Fisher exact p %s",
        nrow(km_df),
        if (fisher_p < 0.001) "< 0.001" else paste("=", signif(fisher_p, 3))
      ),
      x = NULL,
      y = "90-day mortality (%)"
    ) +
    theme_glp_fig()

  save_ms_fig(p_mort_bar, "mortality90_bar_lca_n293", w = 7.5, h = 6.5, dpi = 300)
  fwrite(
    km_bar_df %>%
      mutate(
        cohort = "Full CAF-PINT (n = 293)",
        fisher_exact_p = fisher_p
      ),
    file.path(out_stats, "mortality90_bar_lca_n293.csv")
  )

  km_fit_lca <- survfit(Surv(time, event) ~ lca_grp, data = km_df)
  km_diff_lca <- survdiff(Surv(time, event) ~ lca_grp, data = km_df)
  km_p_lca <- 1 - pchisq(km_diff_lca$chisq, length(km_diff_lca$n) - 1)

  km_stats <- tibble(
    cohort = "Full CAF-PINT (mediation_master + alluv deathday)",
    n = nrow(km_df),
    n_hypo = sum(km_df$lca_grp == "Hypoinflammatory"),
    n_hyper = sum(km_df$lca_grp == "Hyperinflammatory"),
    deaths_hypo = sum(km_df$lca_grp == "Hypoinflammatory" & km_df$event == 1L),
    deaths_hyper = sum(km_df$lca_grp == "Hyperinflammatory" & km_df$event == 1L),
    logrank_chisq = unname(km_diff_lca$chisq),
    logrank_p = km_p_lca,
    note = "Outcome: 90-day mortality. Time = deathday for events (capped at 90); 90 for censored survivors."
  )
  fwrite(km_stats, file.path(out_stats, "km_mortality90_lca_full_cohort_n293.csv"))

  km_fit_tgc <- survfit(Surv(time, event) ~ tgc_grp, data = km_df)
  km_diff_tgc <- survdiff(Surv(time, event) ~ tgc_grp, data = km_df)
  km_p_tgc <- 1 - pchisq(km_diff_tgc$chisq, length(km_diff_tgc$n) - 1)
  fwrite(
    tibble(
      cohort = "Full CAF-PINT (n = 293)",
      n_iit = sum(km_df$tgc_grp == GLP_TGC_IIT_LABEL),
      n_cit = sum(km_df$tgc_grp == GLP_TGC_CIT_LABEL),
      deaths_iit = sum(km_df$tgc_grp == GLP_TGC_IIT_LABEL & km_df$event == 1L),
      deaths_cit = sum(km_df$tgc_grp == GLP_TGC_CIT_LABEL & km_df$event == 1L),
      logrank_chisq = unname(km_diff_tgc$chisq),
      logrank_p = km_p_tgc
    ),
    file.path(out_stats, "km_mortality90_tgc_full_cohort_n293.csv")
  )

  km_theme <- theme_glp_fig() +
    theme(
      legend.position.inside = c(0.78, 0.88),
      legend.position = "inside",
      legend.background = element_rect(fill = alpha("white", 0.92), colour = NA),
      legend.key.width = unit(1.6, "cm")
    )

  if (requireNamespace("survminer", quietly = TRUE)) {
    p_km_lca <- ggsurvplot(
      km_fit_lca,
      data = km_df,
      pval = TRUE,
      pval.coord = c(2, 0.12),
      risk.table = TRUE,
      risk.table.height = 0.28,
      palette = unname(GLP_SUBPHENOTYPE_COLORS[c("Hypoinflammatory", "Hyperinflammatory")]),
      linetype = "solid",
      size = 1.15,
      legend.title = "Subphenotype",
      legend.labs = c("Hypoinflammatory", "Hyperinflammatory"),
      xlab = "Days from randomization",
      ylab = "90-day survival probability",
      title = "90-day survival by baseline inflammatory subphenotype",
      subtitle = paste0(
        "Full CAF-PINT cohort (n = ", nrow(km_df), "); log-rank p = ", signif(km_p_lca, 3)
      ),
      ggtheme = km_theme,
      font.main = c(16, "bold"),
      font.submain = c(12, "plain"),
      font.legend = c(11, "plain")
    )
    save_ms_fig(p_km_lca$plot, "km_mortality90_lca_full_cohort_n293", w = 9, h = 7, dpi = 300)
    save_reviewer_fig(
      p_km_lca$table,
      file.path(out_fig_ms, "km_mortality90_lca_full_cohort_n293_risktable"),
      w = 9,
      h = 2.5,
      dpi = 300
    )

    p_km_tgc <- ggsurvplot(
      km_fit_tgc,
      data = km_df,
      pval = TRUE,
      pval.coord = c(2, 0.12),
      risk.table = TRUE,
      risk.table.height = 0.28,
      palette = c("black", "black"),
      linetype = c("solid", "dashed"),
      size = 1.15,
      legend.title = "Treatment arm",
      legend.labs = c(GLP_TGC_IIT_LABEL, GLP_TGC_CIT_LABEL),
      xlab = "Days from randomization",
      ylab = "90-day survival probability",
      title = "90-day survival by randomized treatment arm",
      subtitle = paste0(
        "Full CAF-PINT cohort (n = ", nrow(km_df), "); log-rank p = ", signif(km_p_tgc, 3),
        ". Solid = IIT; dashed = CIT."
      ),
      ggtheme = km_theme,
      font.main = c(16, "bold"),
      font.submain = c(12, "plain"),
      font.legend = c(11, "plain")
    )
    save_ms_fig(p_km_tgc$plot, "km_mortality90_tgc_full_cohort_n293", w = 9, h = 7, dpi = 300)
    save_reviewer_fig(
      p_km_tgc$table,
      file.path(out_fig_ms, "km_mortality90_tgc_full_cohort_n293_risktable"),
      w = 9,
      h = 2.5,
      dpi = 300
    )
  } else {
    p_km_base <- ggplot(
      tibble(
        time = km_fit_lca$time,
        surv = km_fit_lca$surv,
        strata = rep(names(km_fit_lca$strata), km_fit_lca$strata)
      ),
      aes(x = time, y = surv, color = strata)
    ) +
      geom_step(linewidth = 1.1) +
      scale_color_manual(
        values = GLP_SUBPHENOTYPE_COLORS[c("Hypoinflammatory", "Hyperinflammatory")],
        labels = c("Hypoinflammatory", "Hyperinflammatory")
      ) +
      labs(
        title = "90-day survival by baseline inflammatory subphenotype",
        subtitle = paste0("Full CAF-PINT cohort (n = ", nrow(km_df), "); log-rank p = ", signif(km_p_lca, 3)),
        x = "Days from randomization",
        y = "Survival probability",
        color = "Subphenotype"
      ) +
      theme_glp_fig()
    save_ms_fig(p_km_base, "km_mortality90_lca_full_cohort_n293", w = 9, h = 7, dpi = 300)
  }
}

# Manuscript supplemental Fig 9 alias (hyper baseline-adjusted trajectories)
if (file.exists(file.path(out_fig, "fig03_tgc_baseline_adjusted_lmm_hyper.pdf"))) {
  for (ext in c("pdf", "png", "tif")) {
    src <- file.path(out_fig, paste0("fig03_tgc_baseline_adjusted_lmm_hyper.", ext))
    if (file.exists(src)) {
      file.copy(
        src,
        file.path(out_fig_ms, paste0("supp09_manuscript_hyper_baseline_adjusted_lmm.", ext)),
        overwrite = TRUE
      )
    }
  }
}

# ============================================================================
# Export combined stats tables
# ============================================================================
interactions_df <- bind_rows(all_interactions)
main_effects_df <- if (length(all_main_effects)) bind_rows(all_main_effects) else tibble()
emmeans_df <- bind_rows(all_emmeans)

fwrite(interactions_df, file.path(out_stats, "mixed_model_interaction_tests.csv"))
if (nrow(main_effects_df)) {
  fwrite(main_effects_df, file.path(out_stats, "mixed_model_main_effect_tests.csv"))
}
fwrite(emmeans_df, file.path(out_stats, "mixed_model_emmeans_trajectories.csv"))
if (length(all_contrasts)) {
  contrasts_df <- bind_rows(all_contrasts)
  fwrite(contrasts_df, file.path(out_stats, "mixed_model_pairwise_contrasts.csv"))
}
fwrite(fig5_int, file.path(out_stats, "fig05_proportional_interaction_tests.csv"))
fwrite(fig5_emm, file.path(out_stats, "fig05_proportional_emmeans.csv"))
fwrite(fig6_int_df, file.path(out_stats, "fig06_il6_glp1_interaction_tests.csv"))
fwrite(fig6_emm_df, file.path(out_stats, "fig06_il6_glp1_emmeans.csv"))
fwrite(bind_rows(endo_int), file.path(out_stats, "fig06_endo_prop_interaction_tests.csv"))
fwrite(endo_emm_df, file.path(out_stats, "fig06_endo_prop_emmeans.csv"))

# Methods / Results labeling guide for manuscript
label_guide <- c(
  "# Analysis labeling guide (reviewer response)",
  "",
  "| Output | Role in manuscript |",
  "|--------|-------------------|",
  "| `figures/manuscript/fig01_manuscript_lca_pooled_lmm` | **Primary** Fig 1 (arm-pooled, log scale) |",
  "| `figures/manuscript/fig02_manuscript_mortality_lmm` | **Primary** Fig 2 |",
  "| `figures/manuscript/fig05_manuscript_il6_glp1_hyper` | **Primary** Fig 5 (hyper IL-6 / GLP-1 log fold change) |",
  "| `figures/manuscript/supp09_manuscript_hyper_baseline_adjusted_lmm` | **Supplemental Fig 9** alias |",
  "| `figures/manuscript/km_mortality90_lca_full_cohort_n293` | **New** Kaplan-Meier full cohort (n=293) |",
  "| `figures/analysis/fig01-03_*` | Detailed LMM / boxplot analysis figures |",
  "| `figures/analysis/fig05_proportional_trajectory_lmm` | **Supplemental Fig 8** (log fold change; D2/D4 only) |",
  "| `mixed_model_main_effect_tests.csv` | **Primary** treatment/LCA/mortality main effects |",
  "| `mixed_model_interaction_tests.csv` | **Secondary** group x time interactions |",
  "",
  "See `docs/REVIEWER_RESPONSE_MANUSCRIPT_EDITS.md` for suggested Results/caption text.",
  "",
  "## Primary contrasts",
  "- Fig 1: subphenotype main effect (arm-pooled)",
  "- Fig 2: mortality main effect",
  "- Fig 5 / 6: treatment main effect in baseline-adjusted log fold-change models",
  "- Report interactions as secondary only"
)
writeLines(label_guide, file.path(out_docs, "MANUSCRIPT_LABELING_GUIDE.md"))

manuscript_captions <- c(
  "FIG01_MANUSCRIPT_CAPTION.txt",
  "Fig 1. Longitudinal trajectories of active GLP-1, C-peptide, and glucagon by inflammatory subphenotype.",
  "Model-based mean log10(1 + concentration) +/- 1 SE from linear mixed-effects models pooling treatment arms",
  "(subphenotype x study day + random intercept; n = 84). Blue, hypoinflammatory; orange, hyperinflammatory.",
  "Reported subphenotype fold differences use this arm-pooled model.",
  "",
  "FIG05_MANUSCRIPT_CAPTION.txt",
  "Fig 5. Baseline-adjusted log fold change in IL-6 and active GLP-1 among hyperinflammatory patients by treatment arm (n = 42).",
  "Y-axis: log10(1 + Dx) - log10(1 + D0). Days 2 and 4 from mixed models adjusting for Day-0 baseline.",
  "Primary contrast: treatment main effect. Treatment x time interactions are secondary.",
  "",
  "KM_MANUSCRIPT_CAPTION.txt",
  "Kaplan-Meier curve of 90-day survival by baseline inflammatory subphenotype in the full CAF-PINT cohort (n = 293).",
  "P value from log-rank test comparing hypoinflammatory and hyperinflammatory groups."
)
writeLines(manuscript_captions, file.path(out_docs, "manuscript_figure_captions_draft.txt"))

manifest <- tibble(
  file = c(
    list.files(out_fig, pattern = "\\.(pdf|png)$", full.names = TRUE, recursive = TRUE),
    list.files(out_fig_ms, pattern = "\\.(pdf|png)$", full.names = TRUE, recursive = TRUE)
  )
) %>%
  mutate(category = "manuscript_revision_figure")
fwrite(manifest, file.path(out_stats, "manuscript_revision_figure_manifest.csv"))

sync_fixed_axis_figures(root)
message("  Fixed-axis bundle: ", rev$figures_fixed_axis)

message("\nDone. Outputs:")
message("  Manuscript figures: ", out_fig_ms)
message("  Analysis figures:   ", out_fig)
message("  Stats:              ", out_stats)
message("  Docs:               ", out_docs)
message("  See docs/MANUSCRIPT_LABELING_GUIDE.md for primary vs descriptive labeling.")
