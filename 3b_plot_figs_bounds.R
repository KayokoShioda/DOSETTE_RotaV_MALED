# ============================================================
# FIGURES: IPCW ESTIMATE WITH NONPARAMETRIC BOUNDS
# ============================================================
# Derived from 3_plot_figs_first_outcome.R. Edits are marked ## NPB-CHANGE n.
#   solid line = median IPCW estimate from the main analysis (df_result_*.csv in dir_name)
#   ribbon     = bounds-based 95% interval (df_bounds_*.csv in dir_name_bounds):
#                2.5th bootstrap percentile of the lower bound to
#                97.5th bootstrap percentile of the upper bound
# Same layout, colours and ribbon style as the original figures; only the
# interval is replaced. Bootstrap percentile intervals of the IPCW estimate
# are kept in the Day-730 table (lowerCI / upperCI) but are NOT drawn.
#
# Needs: dir_name (IPCW run), dir_name_bounds (bounds run), outcome_choice,
# country_choice, followup_period and dose windows in the environment.

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)

## NPB-CHANGE 1: two input directories (IPCW results and bounds results); figures go with the bounds run
out_dir        <- file.path("../Rota-TTE/outputs", dir_name)
out_dir_bounds <- file.path("../Rota-TTE/outputs", dir_name_bounds)
fig_dir        <- file.path("../Rota-TTE/figures", dir_name_bounds)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
outcome_label <- if (outcome_choice == "eia") "EIA infection" else "EIA diarrhea"
fig_name <- if (outcome_choice == "eia") "figs_eia_infections" else "figs_eia_diarrhea"
# Panel letters: infection a/c/e, diarrhea b/d/f
panel_letters <- if (outcome_choice == "eia") c("a", "c", "e") else c("b", "d", "f")
# Exception: BR-only EIA infection uses "b" for the RD panel
if (outcome_choice == "eia" && identical(country_choice, "BR")) panel_letters[2] <- "b"

df_result_1 <- fread(file.path(out_dir, "df_result_1.csv"))
df_result_2 <- fread(file.path(out_dir, "df_result_2.csv"))

## NPB-CHANGE 2: read the bounds
df_bounds_1 <- fread(file.path(out_dir_bounds, "df_bounds_1.csv"))
df_bounds_2 <- fread(file.path(out_dir_bounds, "df_bounds_2.csv"))

# Risk ratio (Delayed / Recommended) per day within each bootstrap sample.
# If both risks are 0 (no events yet in either arm), the ratio is 0/0, so RR = 1.
df_rr_boot <- df_result_1 %>%
  dplyr::select(boot_id, surv_time_censoring, rec_risk = cumrisk) %>%
  inner_join(df_result_2 %>% dplyr::select(boot_id, surv_time_censoring, del_risk = cumrisk),
             by = c("boot_id", "surv_time_censoring")) %>%
  mutate(rr = if_else(del_risk == 0 & rec_risk == 0, 1, del_risk / rec_risk))

# Risk difference (Delayed - Recommended) per day within each bootstrap sample.
df_rd_boot <- df_result_1 %>%
  dplyr::select(boot_id, surv_time_censoring, rec_risk = cumrisk) %>%
  inner_join(df_result_2 %>% dplyr::select(boot_id, surv_time_censoring, del_risk = cumrisk),
             by = c("boot_id", "surv_time_censoring")) %>%
  mutate(rd = del_risk - rec_risk)

## NPB-CHANGE 3: RR and RD bounds per bootstrap (Delayed = exposed, Recommended = reference)
##   RR: [lo_del / hi_rec, hi_del / lo_rec]   RD: [lo_del - hi_rec, hi_del - lo_rec]
df_rr_rd_bounds <- df_bounds_1 %>%
  dplyr::select(boot_id, surv_time_censoring, rec_lo = lo, rec_hi = hi) %>%
  inner_join(df_bounds_2 %>% dplyr::select(boot_id, surv_time_censoring, del_lo = lo, del_hi = hi),
             by = c("boot_id", "surv_time_censoring")) %>%
  # rec_hi / rec_lo == 0: ratio is unbounded; Inf (not NA) so quantile() keeps these draws.
  # 0/0 (del_lo == rec_hi == 0, or del_hi == rec_lo == 0) is set to 1, as for the
  # IPCW RR in 3_plot_figs_first_outcome.R
  mutate(rr_lo = case_when(del_lo == 0 & rec_hi == 0 ~ 1,
                           rec_hi > 0 ~ del_lo / rec_hi,
                           TRUE ~ Inf),
         rr_hi = case_when(del_hi == 0 & rec_lo == 0 ~ 1,
                           rec_lo > 0 ~ del_hi / rec_lo,
                           TRUE ~ Inf),
         rd_lo = del_lo - rec_hi,
         rd_hi = del_hi - rec_lo)

fwrite(df_rr_rd_bounds, file.path(out_dir_bounds, "df_rr_rd_bounds.csv"))

# ============================================================
# 32. SUMMARIZE AND PLOT
# ============================================================

# Cumulative risk: median IPCW estimate across bootstrap samples (CI columns kept for the table only)
res_all <- bind_rows(df_result_1, df_result_2) %>%
  group_by(surv_time_censoring, protocol) %>%
  summarise(median  = median(cumrisk),
            lowerCI = quantile(cumrisk, probs = 0.025),
            upperCI = quantile(cumrisk, probs = 0.975),
            .groups = "drop")

## NPB-CHANGE 4: bounds across bootstraps, joined to the IPCW medians
##   lo / hi         = median lower / upper bound (point-estimate bounds)
##   lo_q025 / hi_q975 = 95% bounds-based interval drawn as the ribbon
bounds_all <- bind_rows(df_bounds_1, df_bounds_2) %>%
  group_by(surv_time_censoring, protocol) %>%
  # Percentiles first: summarise() evaluates sequentially, so once lo = median(lo)
  # is defined, a later quantile(lo, ...) would see the scalar, not the column.
  summarise(lo_q025 = quantile(lo, probs = 0.025),
            hi_q975 = quantile(hi, probs = 0.975),
            lo      = median(lo),
            hi      = median(hi),
            .groups = "drop")

res_all <- res_all %>% left_join(bounds_all, by = c("surv_time_censoring", "protocol"))
res_all$protocol <- as.factor(res_all$protocol)

# Dose window boundaries
rec_window_days <- c(rec_dose_1_start, rec_dose_1_end, rec_dose_2_start, rec_dose_2_end)
del_window_days <- c(del_dose_1_start, del_dose_1_end, del_dose_2_start, del_dose_2_end)

# Dark green (Delayed) and coral (Recommended)
protocol_cols <- c(Delayed = "#00CED1", Recommended = "#F8766D")

## NPB-CHANGE 5: ribbon spans the bounds-based 95% interval instead of the IPCW CI
plot_risk <- res_all %>%
  ggplot(aes(x = surv_time_censoring, y = median, colour = protocol)) +
  geom_vline(xintercept = del_window_days, color = protocol_cols[["Delayed"]],
             linetype = "dashed", linewidth = 0.8, alpha = 0.6) +
  geom_vline(xintercept = rec_window_days, color = protocol_cols[["Recommended"]],
             linetype = "dashed", linewidth = 0.8, alpha = 0.6) +
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = lo_q025, ymax = hi_q975, fill = protocol), linetype = 2, alpha = 0.2) +
  scale_colour_manual(values = protocol_cols) +
  scale_fill_manual(values = protocol_cols) +
  labs(x = "Days Since Birth", y = "Cumulative Risk") +
  scale_x_continuous(limits = c(0, followup_period), breaks = seq(0, followup_period, 50)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    legend.position = "bottom",
    legend.text = element_text(size = 12),
    legend.title = element_blank(),
    axis.text = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 18, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20)
  ) +
  ggtitle(paste0(panel_letters[1], ". ", outcome_label))

plot_risk

# Risk ratio (Delayed / Recommended): median IPCW estimate across bootstrap samples
rr_all <- df_rr_boot %>%
  group_by(surv_time_censoring) %>%
  summarise(median  = median(rr),
            lowerCI = quantile(rr, probs = 0.025),
            upperCI = quantile(rr, probs = 0.975),
            .groups = "drop")

## NPB-CHANGE 6: RR / RD bounds across bootstraps, joined to the IPCW medians
##   *_lo / *_hi         = median bounds (point-estimate bounds)
##   *_lo_q025 / *_hi_q975 = 95% bounds-based interval drawn as the ribbon
rr_rd_bounds_all <- df_rr_rd_bounds %>%
  group_by(surv_time_censoring) %>%
  # Percentiles first (see note above on sequential evaluation in summarise)
  summarise(rr_lo_q025 = quantile(rr_lo, probs = 0.025, na.rm = TRUE),
            rr_hi_q975 = quantile(rr_hi, probs = 0.975, na.rm = TRUE),
            rd_lo_q025 = quantile(rd_lo, probs = 0.025, na.rm = TRUE),
            rd_hi_q975 = quantile(rd_hi, probs = 0.975, na.rm = TRUE),
            rr_lo      = median(rr_lo, na.rm = TRUE),
            rr_hi      = median(rr_hi, na.rm = TRUE),
            rd_lo      = median(rd_lo, na.rm = TRUE),
            rd_hi      = median(rd_hi, na.rm = TRUE),
            .groups = "drop")

rr_all <- rr_all %>%
  left_join(rr_rd_bounds_all %>% dplyr::select(surv_time_censoring, rr_lo, rr_hi, rr_lo_q025, rr_hi_q975),
            by = "surv_time_censoring")

## NPB-CHANGE 7: ribbon spans the bounds-based 95% interval instead of the IPCW CI
plot_rr <- rr_all %>%
  ggplot(aes(x = surv_time_censoring, y = median)) +
  geom_vline(xintercept = del_window_days, color = protocol_cols[["Delayed"]],
             linetype = "dashed", linewidth = 0.8, alpha = 0.6) +
  geom_vline(xintercept = rec_window_days, color = protocol_cols[["Recommended"]],
             linetype = "dashed", linewidth = 0.8, alpha = 0.6) +
  geom_line(color = "#1f77b4", linewidth = 1.2) +
  geom_ribbon(aes(ymin = rr_lo_q025, ymax = rr_hi_q975), fill = "#1f77b4", alpha = 0.2, colour = NA) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", linewidth = 0.8) +
  labs(x = "Days Since Birth", y = "Risk Ratio \n (Delayed / Recommended)") +
  scale_x_continuous(limits = c(0, followup_period), breaks = seq(0, followup_period, 50)) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    axis.text = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 18, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20)
  ) +
  ggtitle(paste0(panel_letters[3], ". ", outcome_label))

plot_rr

# Risk difference (Delayed - Recommended): median IPCW estimate across bootstrap samples
rd_all <- df_rd_boot %>%
  group_by(surv_time_censoring) %>%
  summarise(median  = median(rd),
            lowerCI = quantile(rd, probs = 0.025),
            upperCI = quantile(rd, probs = 0.975),
            .groups = "drop") %>%
  left_join(rr_rd_bounds_all %>% dplyr::select(surv_time_censoring, rd_lo, rd_hi, rd_lo_q025, rd_hi_q975),
            by = "surv_time_censoring")                                   ## NPB-CHANGE 8

## NPB-CHANGE 9: ribbon spans the bounds-based 95% interval instead of the IPCW CI
plot_rd <- rd_all %>%
  ggplot(aes(x = surv_time_censoring, y = median)) +
  geom_vline(xintercept = del_window_days, color = protocol_cols[["Delayed"]],
             linetype = "dashed", linewidth = 0.8, alpha = 0.6) +
  geom_vline(xintercept = rec_window_days, color = protocol_cols[["Recommended"]],
             linetype = "dashed", linewidth = 0.8, alpha = 0.6) +
  geom_line(color = "#1f77b4", linewidth = 1.2) +
  geom_ribbon(aes(ymin = rd_lo_q025, ymax = rd_hi_q975), fill = "#1f77b4", alpha = 0.2, colour = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
  labs(x = "Days Since Birth", y = "Risk Difference \n (Delayed - Recommended)") +
  scale_x_continuous(limits = c(0, followup_period), breaks = seq(0, followup_period, 50)) +
  scale_y_continuous(labels = scales::percent) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    axis.text = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 18, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20)
  ) +
  ggtitle(paste0(panel_letters[2], ". ", outcome_label))

plot_rd

library(patchwork)
figs_risk_rd    <- plot_risk + plot_rd
figs_risk_rr_rd <- plot_risk + plot_rd + plot_rr

## NPB-CHANGE 10: file names carry a _bounds suffix
ggsave(file.path(fig_dir, paste0(fig_name, "_bounds_risk_rd.png")),
       figs_risk_rd, width = 16, height = 6, units = "in", dpi = 300)

ggsave(file.path(fig_dir, paste0(fig_name, "_bounds_risk_rr_rd.png")),
       figs_risk_rr_rd, width = 24, height = 6, units = "in", dpi = 300)

## NPB-CHANGE 11: Day-730 table has IPCW median + IPCW 95% CI + median bounds
##   + bounds-based 95% interval (bound_lo_q025 / bound_hi_q975, the plotted ribbon)
res_rr_rd <- bind_rows(
  res_all %>% mutate(measure = paste0("Cumrisk_", as.character(protocol))) %>%
    dplyr::select(measure, surv_time_censoring, median, lowerCI, upperCI,
                  bound_lo = lo, bound_hi = hi, bound_lo_q025 = lo_q025, bound_hi_q975 = hi_q975),
  rr_all %>% mutate(measure = "RR") %>%
    dplyr::select(measure, surv_time_censoring, median, lowerCI, upperCI,
                  bound_lo = rr_lo, bound_hi = rr_hi, bound_lo_q025 = rr_lo_q025, bound_hi_q975 = rr_hi_q975),
  rd_all %>% mutate(measure = "RD") %>%
    dplyr::select(measure, surv_time_censoring, median, lowerCI, upperCI,
                  bound_lo = rd_lo, bound_hi = rd_hi, bound_lo_q025 = rd_lo_q025, bound_hi_q975 = rd_hi_q975)
) %>%
  filter(surv_time_censoring == followup_period)

fwrite(res_rr_rd, file.path(out_dir_bounds, "res_rr_rd_bounds.csv"))
