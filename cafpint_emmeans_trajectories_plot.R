# Biomarkers over time by LCA subclass and TGC arm
# Mixed-effects model with emmeans
# ---------------------------------------------------------------------------
# Faithful adaptation of Daniela Markovic's blue-journal IL-6 / IL-8 template.
# Only data paths and marker names are changed for the GLP manuscript.
# ---------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(lme4)
if (requireNamespace("lmerTest", quietly = TRUE)) {
  library(lmerTest)
}
library(emmeans)
library(ggplot2)
library(grid)
library(patchwork)
library(data.table)
library(car)

# --- Paths ---
args <- commandArgs(trailingOnly = FALSE)
f <- sub("^--file=", "", args[grep("^--file=", args)])
script_dir <- if (length(f)) dirname(normalizePath(f)) else getwd()
repo_root <- normalizePath(file.path(script_dir, "..", ".."))
setwd(repo_root)

out_fig <- file.path(script_dir, "figures")
out_fig_ind <- file.path(out_fig, "individual")
out_stats <- file.path(script_dir, "stats")
dir.create(out_fig, recursive = TRUE, showWarnings = FALSE)
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

# --- Load GLP data ---
outcomes <- fread("cafpint.bms.outcomes.glp.csv") %>% as.data.frame()

classes <- outcomes %>%
  distinct(cafpintid, .keep_all = TRUE) %>%
  transmute(
    cafpintid,
    subclass = lcaprob_1011_dichot25,
    tgc = tgc.x
  )

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
  pivot_wider(names_from = marker, values_from = value)

# Column names for plotting loop
marker_cols <- c(
  "Active GLP-1" = "Active_GLP1",
  "C-Peptide" = "C_Peptide",
  "Glucagon" = "Glucagon",
  "IL-6" = "IL_6"
)

biomarkers <- biomarkers %>%
  rename(
    Active_GLP1 = `Active GLP-1`,
    C_Peptide = `C-Peptide`,
    IL_6 = `IL-6`
  ) %>%
  merge(classes[, c("cafpintid", "subclass")], by = "cafpintid")

biomarkers <- biomarkers %>%
  filter(visit %in% c(0, 2, 4))

biomarkers$visit    <- factor(biomarkers$visit, levels = c(0, 2, 4))
biomarkers$tgc      <- factor(biomarkers$tgc, levels = c(1, 2))
biomarkers$subclass <- factor(biomarkers$subclass, levels = c(0, 1))

# Log transform 
log_transform_safe <- function(x) {
  if (any(x == 0, na.rm = TRUE)) log(x + 1) else log(x)
}

for (col in unname(marker_cols)) {
  biomarkers[[col]] <- log_transform_safe(biomarkers[[col]])
}

y_labels <- c(
  Active_GLP1 = "Log(Active GLP-1 pM)",
  C_Peptide   = "Log(C-peptide pM)",
  Glucagon    = "Log(Glucagon pM)",
  IL_6        = "Log(IL-6 pg/mL)"
)

PLOT_SIZE_IN <- 9
PLOT_TITLE_SIZE <- 16
PLOT_AXIS_TITLE_SIZE <- 14
PLOT_AXIS_TEXT_SIZE <- 12
PLOT_LEGEND_SIZE <- 11
PLOT_ANNOT_SIZE <- 5.5

marker_titles <- c(
  Active_GLP1 = "Active GLP-1",
  C_Peptide   = "C-peptide",
  Glucagon    = "Glucagon",
  IL_6        = "IL-6"
)

format_p_label <- function(p) {
  if (!is.finite(p)) return("NA")
  format(signif(p, 3), scientific = FALSE)
}

p_annotation_label <- function(raw_p) {
  paste0("Interaction p = ", format_p_label(raw_p))
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

interaction_p <- function(model, term) {
  a <- car::Anova(model, type = 3)
  pcol <- intersect(c("Pr(>Chisq)", "Pr(>F)"), colnames(a))[1]
  as.numeric(a[term, pcol])
}

interaction_p_for_marker <- function(data, marker, strat, interaction_term) {
  form <- switch(strat,
    tgc = as.formula(paste(marker, "~ visit * tgc + (1|cafpintid)")),
    lca = as.formula(paste(marker, "~ visit * subclass + (1|cafpintid)")),
    tgc_lca = as.formula(paste(marker, "~ visit * subclass * tgc + (1|cafpintid)"))
  )
  interaction_p(lmer(form, data = data), interaction_term)
}

# --- Plot function (combined: subclass x tgc) ---
make_marker_plot_tgc_lca <- function(data, marker, y_lab, raw_p) {
  form  <- as.formula(paste(marker, "~ visit * subclass * tgc + (1 | cafpintid)"))
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
  color_values    <- c("0_1" = "blue", "0_2" = "blue", "1_1" = "red", "1_2" = "red")
  linetype_values <- c("0_1" = "solid", "0_2" = "22", "1_1" = "solid", "1_2" = "22")

  y_top <- max(emm_df$upper.CL, na.rm = TRUE)

  ggplot(emm_df, aes(x = visit, y = emmean, group = group,
                     color = group, linetype = group)) +
    geom_line(position = position_dodge(0.2), linewidth = 0.9) +
    geom_point(position = position_dodge(0.2), size = 2.5) +
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                  width = 0.2, position = position_dodge(0.2), linewidth = 0.8) +
    scale_color_manual(name = NULL, values = color_values,
                       breaks = group_levels, labels = group_labels) +
    scale_linetype_manual(name = NULL, values = linetype_values,
                          breaks = group_levels, labels = group_labels) +
    add_p_annotate(y_top, raw_p) +
    labs(x = "Day", y = y_lab, title = plot_title(marker)) +
    coord_cartesian(clip = "off") +
    daniela_theme(
      legend.title = element_blank(),
      legend.key.width = unit(1.8, "cm")
    ) +
    guides(
      color    = guide_legend(override.aes = list(shape = NA)),
      linetype = guide_legend(override.aes = list(shape = NA))
    )
}

# --- Similar ggplot style: TGC only (full cohort) ---
make_marker_plot_tgc <- function(data, marker, y_lab, raw_p) {
  form  <- as.formula(paste(marker, "~ visit * tgc + (1 | cafpintid)"))
  model <- lmer(form, data = data)
  emm_df <- as.data.frame(emmeans(model, ~ visit | tgc))
  y_top <- max(emm_df$upper.CL, na.rm = TRUE)

  ggplot(emm_df, aes(x = visit, y = emmean, group = tgc, color = tgc, linetype = tgc)) +
    geom_line(position = position_dodge(0.2), linewidth = 0.9) +
    geom_point(position = position_dodge(0.2), size = 2.5) +
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                  width = 0.2, position = position_dodge(0.2), linewidth = 0.8) +
    scale_color_manual(name = NULL, values = c("1" = "blue", "2" = "red"),
                       labels = c("TGC 80-110", "TGC 150-180")) +
    scale_linetype_manual(name = NULL, values = c("1" = "solid", "2" = "22"),
                          labels = c("TGC 80-110", "TGC 150-180")) +
    add_p_annotate(y_top, raw_p) +
    labs(x = "Day", y = y_lab, title = plot_title(marker)) +
    coord_cartesian(clip = "off") +
    daniela_theme(legend.title = element_blank())
}

# --- Similar ggplot style: LCA only (full cohort) ---
make_marker_plot_lca <- function(data, marker, y_lab, raw_p) {
  form  <- as.formula(paste(marker, "~ visit * subclass + (1 | cafpintid)"))
  model <- lmer(form, data = data)
  emm_df <- as.data.frame(emmeans(model, ~ visit | subclass))
  y_top <- max(emm_df$upper.CL, na.rm = TRUE)

  ggplot(emm_df, aes(x = visit, y = emmean, group = subclass, color = subclass)) +
    geom_line(position = position_dodge(0.2), linewidth = 0.9) +
    geom_point(position = position_dodge(0.2), size = 2.5) +
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                  width = 0.2, position = position_dodge(0.2), linewidth = 0.8) +
    scale_color_manual(name = NULL, values = c("0" = "blue", "1" = "red"),
                       labels = c("Hypoinflamed", "Hyperinflamed")) +
    add_p_annotate(y_top, raw_p) +
    labs(x = "Day", y = y_lab, title = plot_title(marker)) +
    coord_cartesian(clip = "off") +
    daniela_theme(legend.title = element_blank())
}

save_plot <- function(p, fname, w = PLOT_SIZE_IN, h = PLOT_SIZE_IN) {
  ggsave(file.path(out_fig, paste0(fname, ".pdf")), p, width = w, height = h)
  ggsave(file.path(out_fig, paste0(fname, ".png")), p, width = w, height = h, dpi = 220)
}

interaction_rows <- list()
for (strat in c("tgc", "lca", "tgc_lca")) {
  int_term <- switch(strat,
    tgc = "visit:tgc",
    lca = "visit:subclass",
    tgc_lca = "visit:subclass:tgc"
  )
  for (m in names(y_labels)) {
    interaction_rows[[length(interaction_rows) + 1]] <- data.frame(
      stratification = strat,
      marker = m,
      interaction_term = int_term,
      p_value_raw = interaction_p_for_marker(biomarkers, m, strat, int_term),
      stringsAsFactors = FALSE
    )
  }
}
fwrite(
  bind_rows(interaction_rows),
  file.path(out_stats, "interaction_raw_all_strats.csv")
)

# --- Combined TGC x LCA (4 lines, one panel) ---
plots_tgc_lca <- lapply(names(y_labels), function(m) {
  raw_p <- interaction_p_for_marker(biomarkers, m, "tgc_lca", "visit:subclass:tgc")
  p <- make_marker_plot_tgc_lca(biomarkers, m, y_labels[[m]], raw_p)
  save_plot(p, file.path("individual", "tgc_lca", paste0(marker_slug[[m]], "_tgc_lca_daniela_style")))
  p
})
names(plots_tgc_lca) <- names(y_labels)

panel_theme <- theme(
  legend.position = "bottom",
  plot.title = element_text(face = "bold", size = PLOT_TITLE_SIZE),
  axis.title = element_text(face = "bold", size = PLOT_AXIS_TITLE_SIZE),
  axis.text = element_text(face = "bold", size = PLOT_AXIS_TEXT_SIZE),
  legend.text = element_text(face = "bold", size = PLOT_LEGEND_SIZE)
)

final_tgc_lca <- wrap_plots(plots_tgc_lca, ncol = 2) +
  plot_layout(guides = "collect") &
  panel_theme
print(final_tgc_lca)
save_plot(final_tgc_lca, "panel_tgc_lca_", w = PLOT_SIZE_IN * 2, h = PLOT_SIZE_IN * 2)

# --- TGC only ---
plots_tgc <- lapply(names(y_labels), function(m) {
  raw_p <- interaction_p_for_marker(biomarkers, m, "tgc", "visit:tgc")
  p <- make_marker_plot_tgc(biomarkers, m, y_labels[[m]], raw_p)
  save_plot(p, file.path("individual", "tgc", paste0(marker_slug[[m]], "_tgc_daniela_style")))
  p
})
final_tgc <- wrap_plots(plots_tgc, ncol = 2) + plot_layout(guides = "collect") & panel_theme
save_plot(final_tgc, "panel_tgc_", w = PLOT_SIZE_IN * 2, h = PLOT_SIZE_IN * 2)

# --- LCA only ---
plots_lca <- lapply(names(y_labels), function(m) {
  raw_p <- interaction_p_for_marker(biomarkers, m, "lca", "visit:subclass")
  p <- make_marker_plot_lca(biomarkers, m, y_labels[[m]], raw_p)
  save_plot(p, file.path("individual", "lca", paste0(marker_slug[[m]], "_lca_daniela_style")))
  p
})
final_lca <- wrap_plots(plots_lca, ncol = 2) + plot_layout(guides = "collect") & panel_theme
save_plot(final_lca, "panel_lca_", w = PLOT_SIZE_IN * 2, h = PLOT_SIZE_IN * 2)

# --- Export full ANOVA tables ---
anova_rows <- list()
for (strat in c("tgc", "lca", "tgc_lca")) {
  for (m in names(y_labels)) {
    form <- switch(strat,
      tgc = as.formula(paste(m, "~ visit * tgc + (1|cafpintid)")),
      lca = as.formula(paste(m, "~ visit * subclass + (1|cafpintid)")),
      tgc_lca = as.formula(paste(m, "~ visit * subclass * tgc + (1|cafpintid)"))
    )
    fit <- lmer(form, data = biomarkers)
    a <- car::Anova(fit, type = 3)
    anova_rows[[length(anova_rows) + 1]] <- data.frame(
      marker = m,
      stratification = strat,
      term = rownames(a),
      Chisq = a[, "Chisq"],
      Df = a[, "Df"],
      p_value = a[, "Pr(>Chisq)"],
      stringsAsFactors = FALSE
    )
  }
}
anova_out <- bind_rows(anova_rows)
fwrite(anova_out, file.path(out_stats, "_anova_all_terms.csv"))

message("Done. Combined panels: ", out_fig)
message("Individual plots: ", out_fig_ind, " (tgc_lca/, tgc/, lca/)")
message("Stats: ", out_stats)
