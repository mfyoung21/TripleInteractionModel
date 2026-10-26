#!/usr/bin/env Rscript
# n = 84 GLP assay cohort — Daniela-style emmeans trajectories, all stratifications.
#
# Stratifications:
#   lca      — hypo vs hyper (visit * subclass)
#   tgc      — TGC 80-110 vs 150-180 (visit * tgc)
#   tgc_lca  — 4 groups: hypo/hyper x TGC arm (visit * subclass * tgc)
#
# Inference: raw p-values only (4 markers per panel; no FDR).
# Data: cafpint.bms.outcomes.glp.csv (natural log transform, Daniela template).

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(lme4)
  if (requireNamespace("lmerTest", quietly = TRUE)) library(lmerTest)
  library(emmeans)
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(car)
})

args <- commandArgs(trailingOnly = FALSE)
f <- sub("^--file=", "", args[grep("^--file=", args)])
script_dir <- if (length(f)) dirname(normalizePath(f)) else getwd()
repo_root <- normalizePath(file.path(script_dir, "..", ".."))
setwd(repo_root)

out_fig <- file.path(script_dir, "figures")
out_stats <- file.path(script_dir, "stats")
for (sub in c("", "lca", "tgc", "tgc_lca")) {
  dir.create(file.path(out_fig, sub), recursive = TRUE, showWarnings = FALSE)
}
dir.create(out_stats, recursive = TRUE, showWarnings = FALSE)

OUTCOMES_PATH <- "cafpint.bms.outcomes.glp.csv"
PRIOR_ANOVA <- file.path(script_dir, "..", "stats", "daniela_style_anova_all_terms.csv")

MARKERS <- tibble::tibble(
  marker_id = c("Active_GLP1", "C_Peptide", "Glucagon", "IL_6"),
  assay_name = c("Active GLP-1", "C-Peptide", "Glucagon", "IL-6"),
  slug = c("active_glp1", "c_peptide", "glucagon", "il6"),
  ylab = c(
    "Log(Active GLP-1 pM)",
    "Log(C-peptide pM)",
    "Log(Glucagon pM)",
    "Log(IL-6 pg/mL)"
  ),
  title_prefix = c("Active GLP-1", "C-peptide", "Glucagon", "IL-6")
)

STRATS <- list(
  lca = list(
    label = "By subphenotype (hypo vs hyper)",
    formula = value_log ~ visit * subclass + (1 | cafpintid),
    interaction_term = "visit:subclass",
    emm_spec = ~ visit | subclass
  ),
  tgc = list(
    label = "By TGC arm",
    formula = value_log ~ visit * tgc + (1 | cafpintid),
    interaction_term = "visit:tgc",
    emm_spec = ~ visit | tgc
  ),
  tgc_lca = list(
    label = "By subphenotype x TGC arm",
    formula = value_log ~ visit * subclass * tgc + (1 | cafpintid),
    interaction_term = "visit:subclass:tgc",
    emm_spec = ~ visit * tgc | subclass
  )
)

log_transform_safe <- function(x) {
  if (any(x == 0, na.rm = TRUE)) log(x + 1) else log(x)
}

visit_recode <- function(day_chr) {
  dplyr::recode(as.character(day_chr), D0 = "0", D2 = "2", D4 = "4", .default = NA_character_)
}

build_glp_cohort_long <- function(outcomes) {
  classes <- outcomes %>%
    distinct(cafpintid, .keep_all = TRUE) %>%
    transmute(
      cafpintid,
      subclass = lcaprob_1011_dichot25,
      tgc = tgc.x
    )

  assay_long <- outcomes %>%
    filter(Assay %in% c("Active GLP-1", "C-Peptide", "Glucagon")) %>%
    transmute(
      cafpintid,
      assay = Assay,
      visit = as.integer(visit_recode(Day)),
      value_raw = Calc..Concentration
    )

  il6_long <- outcomes %>%
    distinct(cafpintid, Day, .keep_all = TRUE) %>%
    transmute(
      cafpintid,
      assay = "IL-6",
      visit = as.integer(visit_recode(Day)),
      value_raw = il6
    )

  bind_rows(assay_long, il6_long) %>%
    inner_join(classes, by = "cafpintid") %>%
    filter(visit %in% c(0, 2, 4), !is.na(value_raw), is.finite(value_raw)) %>%
    mutate(
      visit = factor(visit, levels = c(0, 2, 4)),
      tgc = factor(tgc, levels = c(1, 2)),
      subclass = factor(subclass, levels = c(0, 1)),
      value_log = log_transform_safe(value_raw)
    )
}

fit_model <- function(data, strat) {
  lmer(
    strat$formula,
    data = data,
    REML = TRUE,
    control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
  )
}

anova_term_p <- function(model, term) {
  a <- car::Anova(model, type = 3)
  pcol <- intersect(c("Pr(>Chisq)", "Pr(>F)"), colnames(a))[1]
  if (!term %in% rownames(a)) return(NA_real_)
  as.numeric(a[term, pcol])
}

full_anova <- function(model) {
  a <- car::Anova(model, type = 3)
  pcol <- intersect(c("Pr(>Chisq)", "Pr(>F)"), colnames(a))[1]
  tibble(
    term = rownames(a),
    Chisq = a[, "Chisq"],
    Df = a[, "Df"],
    p_value = a[, pcol]
  )
}

PLOT_SIZE_IN <- 9
PLOT_TITLE_SIZE <- 16
PLOT_AXIS_TITLE_SIZE <- 14
PLOT_AXIS_TEXT_SIZE <- 12
PLOT_LEGEND_SIZE <- 11
PLOT_ANNOT_SIZE <- 5.5

format_p_label <- function(p) {
  if (!is.finite(p)) return("NA")
  format(signif(p, 3), scientific = FALSE)
}

p_annotation_label <- function(raw_p) {
  paste0("Interaction p = ", format_p_label(raw_p))
}

plot_title <- function(title) {
  paste0("Plasma ", title, " Levels Over Time")
}

daniela_theme <- function(...) {
  theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = PLOT_TITLE_SIZE, hjust = 0.5),
      axis.title = element_text(face = "bold", size = PLOT_AXIS_TITLE_SIZE),
      axis.text = element_text(face = "bold", size = PLOT_AXIS_TEXT_SIZE),
      legend.position = "bottom",
      legend.text = element_text(face = "bold", size = PLOT_LEGEND_SIZE),
      aspect.ratio = 1,
      plot.margin = margin(14, 44, 10, 10),
      ...
    )
}

add_p_annotate <- function(y_top, raw_p) {
  annotate(
    "text",
    x = 4,
    y = y_top * 1.07,
    hjust = 1,
    vjust = 0,
    label = p_annotation_label(raw_p),
    size = PLOT_ANNOT_SIZE,
    fontface = "bold",
    lineheight = 0.95
  )
}

make_plot_lca <- function(data, title, ylab, raw_p) {
  strat <- STRATS$lca
  model <- fit_model(data, strat)
  emm_df <- as.data.frame(emmeans(model, strat$emm_spec))

  y_top <- max(emm_df$upper.CL, na.rm = TRUE)
  ggplot(emm_df, aes(x = visit, y = emmean, group = subclass, color = subclass)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.5) +
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2, linewidth = 0.8) +
    scale_color_manual(
      name = NULL,
      values = c("0" = "blue", "1" = "red"),
      labels = c("Hypoinflamed", "Hyperinflamed")
    ) +
    add_p_annotate(y_top, raw_p) +
    labs(x = "Day", y = ylab, title = plot_title(title)) +
    coord_cartesian(clip = "off") +
    daniela_theme()
}

make_plot_tgc <- function(data, title, ylab, raw_p) {
  strat <- STRATS$tgc
  model <- fit_model(data, strat)
  emm_df <- as.data.frame(emmeans(model, strat$emm_spec))

  y_top <- max(emm_df$upper.CL, na.rm = TRUE)
  ggplot(emm_df, aes(x = visit, y = emmean, group = tgc, color = tgc, linetype = tgc)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.5) +
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2, linewidth = 0.8) +
    scale_color_manual(
      name = NULL, values = c("1" = "blue", "2" = "red"),
      labels = c(expression(paste("Treated with ", TGC[80 - 110])),
                 expression(paste("Treated with ", TGC[150 - 180])))
    ) +
    scale_linetype_manual(
      name = NULL, values = c("1" = "solid", "2" = "22"),
      labels = c(expression(paste("Treated with ", TGC[80 - 110])),
                 expression(paste("Treated with ", TGC[150 - 180])))
    ) +
    add_p_annotate(y_top, raw_p) +
    labs(x = "Day", y = ylab, title = plot_title(title)) +
    coord_cartesian(clip = "off") +
    daniela_theme()
}

make_plot_tgc_lca <- function(data, title, ylab, raw_p) {
  strat <- STRATS$tgc_lca
  model <- fit_model(data, strat)
  emm_df <- as.data.frame(emmeans(model, strat$emm_spec)) %>%
    mutate(group = interaction(subclass, tgc, sep = "_", lex.order = FALSE))

  group_levels <- c("1_1", "1_2", "0_1", "0_2")
  group_labels <- list(
    "1_1" = expression(paste("Hyperinflamed & Treated with ", TGC[80 - 110])),
    "1_2" = expression(paste("Hyperinflamed & Treated with ", TGC[150 - 180])),
    "0_1" = expression(paste("Hypoinflamed & Treated with ", TGC[80 - 110])),
    "0_2" = expression(paste("Hypoinflamed & Treated with ", TGC[150 - 180]))
  )

  y_top <- max(emm_df$upper.CL, na.rm = TRUE)
  ggplot(emm_df, aes(x = visit, y = emmean, group = group, color = group, linetype = group)) +
    geom_line(position = position_dodge(0.2), linewidth = 0.9) +
    geom_point(position = position_dodge(0.2), size = 2.5) +
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2,
                  position = position_dodge(0.2), linewidth = 0.8) +
    scale_color_manual(name = NULL, values = c("1_1" = "red", "1_2" = "red", "0_1" = "blue", "0_2" = "blue"),
                       breaks = group_levels, labels = group_labels) +
    scale_linetype_manual(name = NULL, values = c("1_1" = "solid", "1_2" = "22", "0_1" = "solid", "0_2" = "22"),
                          breaks = group_levels, labels = group_labels) +
    add_p_annotate(y_top, raw_p) +
    labs(x = "Day", y = ylab, title = plot_title(title)) +
    coord_cartesian(clip = "off") +
    daniela_theme(legend.key.width = unit(1.8, "cm")) +
    guides(
      color = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(shape = NA)),
      linetype = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(shape = NA))
    )
}

plot_fn <- list(lca = make_plot_lca, tgc = make_plot_tgc, tgc_lca = make_plot_tgc_lca)

save_plot <- function(p, path_no_ext, w = PLOT_SIZE_IN, h = PLOT_SIZE_IN) {
  ggsave(paste0(path_no_ext, ".pdf"), p, width = w, height = h)
  ggsave(paste0(path_no_ext, ".png"), p, width = w, height = h, dpi = 220)
}

# --- Load and audit ---
stopifnot(file.exists(OUTCOMES_PATH))
outcomes <- fread(OUTCOMES_PATH) %>% as.data.frame()
cohort_ids <- sort(unique(outcomes$cafpintid))
stopifnot(length(cohort_ids) == 84)

all_long <- build_glp_cohort_long(outcomes)
marker_data <- split(all_long, all_long$assay)

pt_meta <- outcomes %>%
  distinct(cafpintid, .keep_all = TRUE) %>%
  transmute(cafpintid, subclass = lcaprob_1011_dichot25, tgc = tgc.x)

audit_rows <- MARKERS %>%
  rowwise() %>%
  mutate(
    n_patients = n_distinct(marker_data[[assay_name]]$cafpintid),
    n_obs = nrow(marker_data[[assay_name]]),
    n_visit0 = sum(marker_data[[assay_name]]$visit == 0),
    n_visit2 = sum(marker_data[[assay_name]]$visit == 2),
    n_visit4 = sum(marker_data[[assay_name]]$visit == 4),
    value_raw_min = min(marker_data[[assay_name]]$value_raw, na.rm = TRUE),
    value_raw_max = max(marker_data[[assay_name]]$value_raw, na.rm = TRUE)
  ) %>%
  ungroup()

cohort_audit <- tibble(
  metric = c(
    "n_patients", "hypo_n", "hyper_n", "tgc80_n", "tgc150_n",
    "hypo_tgc80", "hypo_tgc150", "hyper_tgc80", "hyper_tgc150"
  ),
  value = c(
    length(cohort_ids),
    sum(pt_meta$subclass == 0), sum(pt_meta$subclass == 1),
    sum(pt_meta$tgc == 1), sum(pt_meta$tgc == 2),
    sum(pt_meta$subclass == 0 & pt_meta$tgc == 1),
    sum(pt_meta$subclass == 0 & pt_meta$tgc == 2),
    sum(pt_meta$subclass == 1 & pt_meta$tgc == 1),
    sum(pt_meta$subclass == 1 & pt_meta$tgc == 2)
  )
)

fwrite(all_long, file.path(out_stats, "n84_analysis_long.csv"))
fwrite(audit_rows, file.path(out_stats, "data_audit_by_marker.csv"))
fwrite(cohort_audit, file.path(out_stats, "cohort_audit.csv"))

# --- Fit all models ---
interaction_rows <- list()
anova_rows <- list()

for (strat_key in names(STRATS)) {
  strat <- STRATS[[strat_key]]

  for (i in seq_len(nrow(MARKERS))) {
    meta <- MARKERS[i, ]
    d <- marker_data[[meta$assay_name]]
    model <- fit_model(d, strat)
    raw_p <- anova_term_p(model, strat$interaction_term)

    anova_rows[[length(anova_rows) + 1]] <- full_anova(model) %>%
      mutate(
        marker = meta$marker_id,
        assay = meta$assay_name,
        stratification = strat_key,
        n_patients = n_distinct(d$cafpintid),
        n_obs = nrow(d),
        interaction_term = strat$interaction_term
      )

    interaction_rows[[length(interaction_rows) + 1]] <- tibble(
      marker = meta$marker_id,
      assay = meta$assay_name,
      stratification = strat_key,
      interaction_term = strat$interaction_term,
      n_patients = n_distinct(d$cafpintid),
      n_obs = nrow(d),
      p_value_raw = raw_p
    )
  }
}

interaction_df <- bind_rows(interaction_rows)
anova_df <- bind_rows(anova_rows)

fwrite(interaction_df, file.path(out_stats, "interaction_raw_all_strats.csv"))
fwrite(anova_df, file.path(out_stats, "anova_all_terms_all_strats.csv"))

# --- Cross-check vs prior daniela folder run ---
if (file.exists(PRIOR_ANOVA)) {
  prior <- fread(PRIOR_ANOVA)

  verify <- interaction_df %>%
    left_join(
      prior %>%
        filter(term %in% c("visit:tgc", "visit:subclass", "visit:subclass:tgc")) %>%
        dplyr::select(marker, stratification, term, p_value_prior = p_value),
      by = c("marker", "stratification", "interaction_term" = "term")
    ) %>%
    mutate(
      abs_diff = abs(p_value_raw - p_value_prior),
      match_ok = is.finite(abs_diff) & abs_diff < 1e-5
    )

  fwrite(verify, file.path(out_stats, "verification_vs_daniela_style_anova.csv"))
  n_bad <- sum(!verify$match_ok, na.rm = TRUE)
  if (n_bad > 0) {
    warning("Verification: ", n_bad, " interaction p-values differ from prior daniela_style_anova_all_terms.csv")
  } else {
    message("Verification passed: interaction p-values match prior daniela_style_anova_all_terms.csv")
  }
}

# --- Plots ---
for (strat_key in names(STRATS)) {
  plots <- list()
  strat_dir <- file.path(out_fig, strat_key)

  for (i in seq_len(nrow(MARKERS))) {
    meta <- MARKERS[i, ]
    d <- marker_data[[meta$assay_name]]
    stat_row <- interaction_df %>%
      filter(stratification == strat_key, marker == meta$marker_id)

    p <- plot_fn[[strat_key]](
      d, meta$title_prefix, meta$ylab,
      stat_row$p_value_raw
    )
    plots[[meta$marker_id]] <- p
    save_plot(p, file.path(strat_dir, paste0(meta$slug, "_", strat_key)))
  }

  panel <- (plots[["Active_GLP1"]] | plots[["C_Peptide"]]) /
    (plots[["Glucagon"]] | plots[["IL_6"]]) +
    plot_layout(guides = "collect") &
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold", size = PLOT_TITLE_SIZE),
      axis.title = element_text(face = "bold", size = PLOT_AXIS_TITLE_SIZE),
      axis.text = element_text(face = "bold", size = PLOT_AXIS_TEXT_SIZE),
      legend.text = element_text(face = "bold", size = PLOT_LEGEND_SIZE)
    )

  save_plot(
    panel, file.path(out_fig, paste0("panel_", strat_key, "_n84")),
    w = PLOT_SIZE_IN * 2, h = PLOT_SIZE_IN * 2
  )
}

sig <- interaction_df %>%
  filter(p_value_raw < 0.05) %>%
  arrange(stratification, p_value_raw)
fwrite(sig, file.path(out_stats, "significant_interactions_raw05.csv"))

message("\n=== n = 84 GLP stratified analysis complete ===")
message("Cohort: ", length(cohort_ids), " patients | hypo/hyper: ",
        cohort_audit$value[2], "/", cohort_audit$value[3])
print(interaction_df %>%
        dplyr::select(stratification, marker, p_value_raw) %>%
        arrange(stratification, marker))
message("\nOutputs: ", script_dir)
