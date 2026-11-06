# ============================================================
# IPCW DIAGNOSTICS (ORIGINAL cox_data, NO BOOTSTRAP)
# ============================================================
# Requires (from 1_Cox_dataset.R): cox_data, fm_cn, followup_period.
# IPCW come only from fit_cox_all(), computed for all children at their own
# surv_time (event time, censoring time, or end of follow-up).

library(dplyr)
library(survival)
library(broom)
library(data.table)
library(ggplot2)

# ============================================================
# COX IPCW FUNCTION (ONE PROTOCOL)
# ============================================================
# Censoring Cox model as in 2_Cox_IPCW.R; returns the IPCW (wt) of every child
# in d, with covariates, sorted by surv_time.
# d = cox_data of one protocol

fit_cox_all <- function(d) {
  
  fml <- fm_cn
  # Country_ID can't be a covariate in single-country analyses
  if (n_distinct(d$Country_ID) < 2) fml <- update(fml, . ~ . - Country_ID)
  
  coxph.censor <- coxph(fml, data = d, model = TRUE)
  
  broom::augment(coxph.censor, newdata = d, type.predict = "expected") %>%
    dplyr::mutate(prob = exp(-.fitted),
                  wt = 1 / prob) %>%
    arrange(surv_time)
}

# ------------------------------------------------------------
# Weights from fit_cox_all() on the original data, by protocol
# ------------------------------------------------------------
wt_data <- cox_data %>%
  split(.$protocol) %>%
  lapply(fit_cox_all) %>%
  bind_rows(.id = "protocol")

# ------------------------------------------------------------
# 1. Overall distribution of IPCW by protocol
# ------------------------------------------------------------
# Children still uncensored and outcome-free at Day 730 have surv_time = 730,
# so their fit_cox_all() weight is the weight at Day 730.
wt_summary <- wt_data %>%
  filter(censored == 0) %>%
  group_by(protocol) %>%
  summarise(
    n      = n(),
    ESS    = sum(wt)^2 / sum(wt^2),
    Mean   = mean(wt),
    SD     = sd(wt),
    Min    = min(wt),
    P1     = quantile(wt, 0.01),
    P25    = quantile(wt, 0.25),
    Median = median(wt),
    P75    = quantile(wt, 0.75),
    P99    = quantile(wt, 0.99),
    Max    = max(wt),
    .groups = "drop"
  )

wt_summary_730 <- wt_summary %>%
  mutate(across(where(is.numeric) & !n, ~ round(.x, 2)))

wt_out_dir <- file.path("../Rota-TTE/suppl_outputs", paste(country_choice, collapse = "_"))
dir.create(wt_out_dir, recursive = TRUE, showWarnings = FALSE)
fwrite(wt_summary_730, file.path(wt_out_dir, paste0("wt_summary_", outcome_choice, ".csv")))

# ------------------------------------------------------------
# 2. Density of IPCW by protocol (uncensored children)
# ------------------------------------------------------------
wt_density_plot <- wt_data %>%
  filter(censored == 0) %>%
  group_by(protocol) %>%
  ggplot(aes(x = wt, colour = protocol, fill = protocol)) +
  geom_density(linewidth = 1.2, alpha = 0.2) +
  # Same palette as 3_plot_figs_first_outcome.R: turquoise (Delayed), coral (Recommended)
  scale_colour_manual(values = c(Delayed = "#00CED1", Recommended = "#F8766D")) +
  scale_fill_manual(values = c(Delayed = "#00CED1", Recommended = "#F8766D")) +
  labs(x = "IPCW", y = "Density",
       title = if (outcome_choice == "eia") "EIA Infection" else "EIA Diarrhea") +
  # IPCW >= 1 by construction, so mark 1 as the reference break
  scale_x_continuous(limits = c(0, 15), breaks = c(1, 5, 10, 15)) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    legend.position = "bottom",
    legend.text = element_text(size = 12),
    legend.title = element_blank(),
    axis.text = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 18, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20)
  )

wt_density_plot
ggsave(file.path(wt_out_dir, paste0("wt_density_", outcome_choice, ".png")),
       wt_density_plot, width = 8, height = 6, units = "in", dpi = 300)
