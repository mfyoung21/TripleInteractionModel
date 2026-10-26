# n = 84 GLP incretin trajectories — figures with table-matching p-value blocks (raw p only)
# TGC panel:     visit:tgc + visit
# LCA panel:     visit:subclass + subclass
# TGC×LCA panel: visit:subclass:tgc + visit:subclass + visit:tgc + subclass
# Note: no FDR — only 4 markers per panel (< 10 tests).

library(dplyr)
library(tidyr)
library(lme4)
if (requireNamespace("lmerTest", quietly = TRUE)) {
  library(lmerTest)
}
library(emmeans)
library(ggplot2)
library(patchwork)
library(data.table)
library(car)

args <- commandArgs(trailingOnly = FALSE)
f <- sub("^--file=", "", args[grep("^--file=", args)])
script_dir <- if (length(f)) dirname(normalizePath(f)) else getwd()
repo_root <- normalizePath(file.path(script_dir, "..", ".."))
setwd(repo_root)

out_fig <- file.path(script_dir, "figures")
out_fig_ind <- file.path(out_fig, "individual")
out_stats <- file.path(script_dir, "stats")
dir.create(file.path(out_fig_ind, "tgc_lca"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_fig_ind, "tgc"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_fig_ind, "lca"), recursive = TRUE, showWarnings = FALSE)
dir.create(out_stats, recursive = TRUE, showWarnings = FALSE)

marker_slug <- c(
  Active_GLP1 = "active_glp1",
  C_Peptide = "c_peptide",
  Glucagon = "glucagon",
  IL_6 = "il6"
)

outcomes <- fread("cafpint.bms.outcomes.glp.csv") %>% as.data.frame()

classes <- outcomes %>%
  distinct(cafpintid, .keep_all = TRUE) %>%
  transmute(cafpintid, subclass = lcaprob_1011_dichot25, tgc = tgc.x)

biomarkers <- outcomes %>%
  filter(Assay %in% c("Active GLP-1", "C-Peptide", "Glucagon")) %>%
  transmute(
    cafpintid,
    visit = as.integer(dplyr::recode(as.character(Day), D0 = "0", D2 = "2", D4 = "4")),
    tgc = tgc.x,
    marker = Assay,
    value = Calc..Concentration
  ) %>%
  bind_rows(
    outcomes %>%
      distinct(cafpintid, Day, .keep_all = TRUE) %>%
      transmute(
        cafpintid,
        visit = as.integer(dplyr::recode(as.character(Day), D0 = "0", D2 = "2", D4 = "4")),
        tgc = tgc.x,
        marker = "IL-6",
        value = il6
      )
  ) %>%
  pivot_wider(names_from = marker, values_from = value) %>%
  rename(
    Active_GLP1 = `Active GLP-1`,
    C_Peptide = `C-Peptide`,
    IL_6 = `IL-6`
  ) %>%
  merge(classes[, c("cafpintid", "subclass")], by = "cafpintid") %>%
  filter(visit %in% c(0, 2, 4))

biomarkers$visit <- factor(biomarkers$visit, levels = c(0, 2, 4))
biomarkers$tgc <- factor(biomarkers$tgc, levels = c(1, 2))
biomarkers$subclass <- factor(biomarkers$subclass, levels = c(0, 1))

log_transform_safe <- function(x) {
  if (any(x == 0, na.rm = TRUE)) log(x + 1) else log(x)
}

for (m in names(marker_slug)) {
  biomarkers[[m]] <- log_transform_safe(biomarkers[[m]])
}

y_labels <- c(
  Active_GLP1 = "Log(Active GLP-1 pM)",
  C_Peptide = "Log(C-peptide pM)",
  Glucagon = "Log(Glucagon pM)",
  IL_6 = "Log(IL-6 pg/mL)"
)

marker_titles <- c(
  Active_GLP1 = "Active GLP-1",
  C_Peptide = "C-peptide",
  Glucagon = "Glucagon",
  IL_6 = "IL-6"
)

PLOT_SIZE_IN <- 9
PLOT_TITLE_SIZE <- 16
PLOT_AXIS_TITLE_SIZE <- 14
PLOT_AXIS_TEXT_SIZE <- 12
PLOT_LEGEND_SIZE <- 11
PLOT_ANNOT_SIZE <- 3.6

format_p_label <- function(p) {
  if (!is.finite(p)) return("NA")
  if (p < 0.001) return("< 0.001")
  format(signif(p, 3), scientific = FALSE)
}

term_p <- function(model, term) {
  a <- car::Anova(model, type = 3)
  pcol <- intersect(c("Pr(>Chisq)", "Pr(>F)"), colnames(a))[1]
  as.numeric(a[term, pcol])
}

format_p_line <- function(term_label, p) {
  paste0(term_label, ": ", format_p_label(p))
}

p_block_tgc <- function(model) {
  paste0(
    format_p_line("visit:tgc", term_p(model, "visit:tgc")), "\n",
    format_p_line("visit", term_p(model, "visit"))
  )
}

p_block_lca <- function(model) {
  paste0(
    format_p_line("visit:subclass", term_p(model, "visit:subclass")), "\n",
    format_p_line("subclass", term_p(model, "subclass"))
  )
}

p_block_tgc_lca <- function(model) {
  paste0(
    format_p_line("visit:subclass:tgc", term_p(model, "visit:subclass:tgc")), "\n",
    format_p_line("visit:subclass", term_p(model, "visit:subclass")), "\n",
    format_p_line("visit:tgc", term_p(model, "visit:tgc")), "\n",
    format_p_line("subclass", term_p(model, "subclass"))
  )
}

plot_title <- function(marker) {
  paste0("Plasma ", marker_titles[[marker]], " Levels Over Time")
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
      plot.margin = margin(12, 12, 10, 10),
      ...
    )
}

add_stats_to_plot <- function(p, label, n_lines) {
  top_expand <- if (n_lines >= 4) 0.48 else 0.30
  stats_df <- data.frame(x = 4, y = Inf, label = label, stringsAsFactors = FALSE)

  p +
    scale_y_continuous(expand = expansion(mult = c(0.06, top_expand))) +
    geom_label(
      data = stats_df,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 1,
      vjust = 1.08,
      size = PLOT_ANNOT_SIZE,
      fontface = "bold",
      label.size = 0.2,
      label.padding = unit(0.3, "lines"),
      fill = alpha("white", 0.94),
      color = "black",
      lineheight = 0.85
    )
}

make_marker_plot_tgc <- function(data, marker, y_lab) {
  form <- as.formula(paste(marker, "~ visit * tgc + (1 | cafpintid)"))
  model <- lmer(form, data = data)
  emm_df <- as.data.frame(emmeans(model, ~ visit | tgc))
  p_label <- p_block_tgc(model)

  p <- ggplot(emm_df, aes(x = visit, y = emmean, group = tgc, color = tgc, linetype = tgc)) +
    geom_line(position = position_dodge(0.2), linewidth = 0.9) +
    geom_point(position = position_dodge(0.2), size = 2.5) +
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                  width = 0.2, position = position_dodge(0.2), linewidth = 0.8) +
    scale_color_manual(name = NULL, values = c("1" = "blue", "2" = "red"),
                       labels = c("TGC 80-110", "TGC 150-180")) +
    scale_linetype_manual(name = NULL, values = c("1" = "solid", "2" = "22"),
                          labels = c("TGC 80-110", "TGC 150-180")) +
    labs(x = "Day", y = y_lab, title = plot_title(marker)) +
    daniela_theme(legend.title = element_blank())

  add_stats_to_plot(p, p_label, 2)
}

make_marker_plot_lca <- function(data, marker, y_lab) {
  form <- as.formula(paste(marker, "~ visit * subclass + (1 | cafpintid)"))
  model <- lmer(form, data = data)
  emm_df <- as.data.frame(emmeans(model, ~ visit | subclass))
  p_label <- p_block_lca(model)

  p <- ggplot(emm_df, aes(x = visit, y = emmean, group = subclass, color = subclass)) +
    geom_line(position = position_dodge(0.2), linewidth = 0.9) +
    geom_point(position = position_dodge(0.2), size = 2.5) +
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                  width = 0.2, position = position_dodge(0.2), linewidth = 0.8) +
    scale_color_manual(name = NULL, values = c("0" = "blue", "1" = "red"),
                       labels = c("Hypoinflamed", "Hyperinflamed")) +
    labs(x = "Day", y = y_lab, title = plot_title(marker)) +
    daniela_theme(legend.title = element_blank())

  add_stats_to_plot(p, p_label, 2)
}

make_marker_plot_tgc_lca <- function(data, marker, y_lab) {
  form <- as.formula(paste(marker, "~ visit * subclass * tgc + (1 | cafpintid)"))
  model <- lmer(form, data = data)
  emm_df <- as.data.frame(emmeans(model, ~ visit * tgc | subclass)) %>%
    mutate(group = interaction(subclass, tgc, sep = "_"))

  group_levels <- c("0_1", "0_2", "1_1", "1_2")
  group_labels <- list(
    "0_1" = expression(paste("Hypoinflamed, ", plain(TGC)[80 - 110])),
    "0_2" = expression(paste("Hypoinflamed, ", plain(TGC)[150 - 180])),
    "1_1" = expression(paste("Hyperinflamed, ", plain(TGC)[80 - 110])),
    "1_2" = expression(paste("Hyperinflamed, ", plain(TGC)[150 - 180]))
  )
  color_values <- c("0_1" = "blue", "0_2" = "blue", "1_1" = "red", "1_2" = "red")
  linetype_values <- c("0_1" = "solid", "0_2" = "22", "1_1" = "solid", "1_2" = "22")
  p_label <- p_block_tgc_lca(model)

  p <- ggplot(emm_df, aes(x = visit, y = emmean, group = group,
                     color = group, linetype = group)) +
    geom_line(position = position_dodge(0.2), linewidth = 0.9) +
    geom_point(position = position_dodge(0.2), size = 2.5) +
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                  width = 0.2, position = position_dodge(0.2), linewidth = 0.8) +
    scale_color_manual(name = NULL, values = color_values,
                       breaks = group_levels, labels = group_labels) +
    scale_linetype_manual(name = NULL, values = linetype_values,
                          breaks = group_levels, labels = group_labels) +
    labs(x = "Day", y = y_lab, title = plot_title(marker)) +
    daniela_theme(
      legend.title = element_blank(),
      legend.key.width = unit(1.8, "cm")
    ) +
    guides(
      color = guide_legend(override.aes = list(shape = NA)),
      linetype = guide_legend(override.aes = list(shape = NA))
    )

  add_stats_to_plot(p, p_label, 4)
}

save_plot <- function(p, fname, w = PLOT_SIZE_IN, h = PLOT_SIZE_IN) {
  ggsave(file.path(out_fig, paste0(fname, ".pdf")), p, width = w, height = h)
  ggsave(file.path(out_fig, paste0(fname, ".png")), p, width = w, height = h, dpi = 220)
}

collect_annotation_rows <- function(strat, marker, model) {
  terms <- switch(strat,
    tgc = c("visit:tgc", "visit"),
    lca = c("visit:subclass", "subclass"),
    tgc_lca = c("visit:subclass:tgc", "visit:subclass", "visit:tgc", "subclass")
  )
  lapply(terms, function(tr) {
    data.frame(
      stratification = strat,
      marker = marker,
      term = tr,
      p_value_raw = term_p(model, tr),
      stringsAsFactors = FALSE
    )
  }) %>% bind_rows()
}

panel_theme <- theme(
  legend.position = "bottom",
  plot.title = element_text(face = "bold", size = PLOT_TITLE_SIZE),
  axis.title = element_text(face = "bold", size = PLOT_AXIS_TITLE_SIZE),
  axis.text = element_text(face = "bold", size = PLOT_AXIS_TEXT_SIZE),
  legend.text = element_text(face = "bold", size = PLOT_LEGEND_SIZE)
)

annotation_rows <- list()
plots_tgc <- lapply(names(y_labels), function(m) {
  form <- as.formula(paste(m, "~ visit * tgc + (1|cafpintid)"))
  model <- lmer(form, data = biomarkers)
  annotation_rows[[paste0("tgc_", m)]] <<-
    collect_annotation_rows("tgc", m, model)
  p <- make_marker_plot_tgc(biomarkers, m, y_labels[[m]])
  save_plot(p, file.path("individual", "tgc", paste0(marker_slug[[m]], "_tgc_table_annotated")))
  p
})
final_tgc <- wrap_plots(plots_tgc, ncol = 2) +
  plot_layout(guides = "collect") &
  panel_theme &
  plot_annotation(
    title = "Panel A: TGC pooled (n = 84)",
    subtitle = "Do intensive vs standard insulin trajectories differ, averaging hypo + hyper?",
    theme = theme(
      plot.title = element_text(face = "bold", size = 17, hjust = 0.5),
      plot.subtitle = element_text(face = "italic", size = 12, hjust = 0.5)
    )
  )
save_plot(final_tgc, "panel_tgc_table_annotated", w = PLOT_SIZE_IN * 2, h = PLOT_SIZE_IN * 2 + 0.8)

plots_lca <- lapply(names(y_labels), function(m) {
  form <- as.formula(paste(m, "~ visit * subclass + (1|cafpintid)"))
  model <- lmer(form, data = biomarkers)
  annotation_rows[[paste0("lca_", m)]] <<-
    collect_annotation_rows("lca", m, model)
  p <- make_marker_plot_lca(biomarkers, m, y_labels[[m]])
  save_plot(p, file.path("individual", "lca", paste0(marker_slug[[m]], "_lca_table_annotated")))
  p
})
final_lca <- wrap_plots(plots_lca, ncol = 2) +
  plot_layout(guides = "collect") &
  panel_theme &
  plot_annotation(
    title = "Panel B: LCA pooled (n = 84)",
    subtitle = "Do hypo vs hyper subphenotype trajectories differ, averaging both TGC arms?",
    theme = theme(
      plot.title = element_text(face = "bold", size = 17, hjust = 0.5),
      plot.subtitle = element_text(face = "italic", size = 12, hjust = 0.5)
    )
  )
save_plot(final_lca, "panel_lca_table_annotated", w = PLOT_SIZE_IN * 2, h = PLOT_SIZE_IN * 2 + 0.8)

plots_tgc_lca <- lapply(names(y_labels), function(m) {
  form <- as.formula(paste(m, "~ visit * subclass * tgc + (1|cafpintid)"))
  model <- lmer(form, data = biomarkers)
  annotation_rows[[paste0("tgc_lca_", m)]] <<-
    collect_annotation_rows("tgc_lca", m, model)
  p <- make_marker_plot_tgc_lca(biomarkers, m, y_labels[[m]])
  save_plot(p, file.path("individual", "tgc_lca", paste0(marker_slug[[m]], "_tgc_lca_table_annotated")))
  p
})
final_tgc_lca <- wrap_plots(plots_tgc_lca, ncol = 2) +
  plot_layout(guides = "collect") &
  panel_theme &
  plot_annotation(
    title = "Panel C: TGC x LCA factorial (n = 84)",
    subtitle = "Does treatment effect differ by subphenotype over time? (4 groups, 3-way interaction)",
    theme = theme(
      plot.title = element_text(face = "bold", size = 17, hjust = 0.5),
      plot.subtitle = element_text(face = "italic", size = 12, hjust = 0.5)
    )
  )
save_plot(final_tgc_lca, "panel_tgc_lca_table_annotated", w = PLOT_SIZE_IN * 2, h = PLOT_SIZE_IN * 2 + 0.8)

annot_out <- bind_rows(annotation_rows)
fwrite(annot_out, file.path(out_stats, "figure_annotation_pvalues.csv"))

message("Done. Table-annotated figures (raw p only): ", out_fig)
message("Annotation CSV: ", file.path(out_stats, "figure_annotation_pvalues.csv"))
