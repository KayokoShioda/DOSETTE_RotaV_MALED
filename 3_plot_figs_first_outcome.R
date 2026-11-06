library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)

# dir_name, outcome_choice, followup_period and dose windows must exist in the environment
out_dir <- file.path("../Rota-TTE/outputs", dir_name)
fig_dir <- file.path("../Rota-TTE/figures", dir_name)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
outcome_label <- switch(outcome_choice, eia = "EIA infection",
                        severe_diarrhea = "Severe diarrhea", "EIA diarrhea")
fig_name <- switch(outcome_choice, eia = "figs_eia_infections",
                   severe_diarrhea = "figs_severe_diarrhea", "figs_eia_diarrhea")
# Panel letters: infection a/c/e, diarrhea b/d/f, none otherwise
panel_letters <- switch(outcome_choice,
                        eia = c("a. ", "c. ", "e. "),
                        eia_diarrhea = c("b. ", "d. ", "f. "),
                        c("", "", ""))
# Exception: BR-only EIA infection uses "b" for the RD panel
if (outcome_choice == "eia" && identical(country_choice, "BR")) panel_letters[2] <- "b. "

df_result_1 <- fread(file.path(out_dir, "df_result_1.csv"))
df_result_2 <- fread(file.path(out_dir, "df_result_2.csv"))

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

fwrite(df_rr_boot, file.path(out_dir, "df_rr_boot.csv"))
fwrite(df_rd_boot, file.path(out_dir, "df_rd_boot.csv"))

# ============================================================
# 32. SUMMARIZE AND PLOT
# ============================================================

# Cumulative risk: median and 95% CI across bootstrap samples
res_all <- bind_rows(df_result_1, df_result_2) %>%
  group_by(surv_time_censoring, protocol) %>%
  summarise(median  = median(cumrisk),
            lowerCI = quantile(cumrisk, probs = 0.025),
            upperCI = quantile(cumrisk, probs = 0.975),
            .groups = "drop")

res_all$protocol <- as.factor(res_all$protocol)

# Dose window boundaries (defined in 1_Cox_dataset.R / 4_PersonDayDataset_first_outcome.R)
rec_window_days <- c(rec_dose_1_start, rec_dose_1_end, rec_dose_2_start, rec_dose_2_end)
del_window_days <- c(del_dose_1_start, del_dose_1_end, del_dose_2_start, del_dose_2_end)

# Dark turquoise (Delayed) and coral (Recommended)
protocol_cols <- c(Delayed = "#00CED1", Recommended = "#F8766D")

plot_risk <- res_all %>%
  ggplot(aes(x = surv_time_censoring, y = median, colour = protocol)) +
  geom_vline(xintercept = del_window_days, color = protocol_cols[["Delayed"]],
             linetype = "dashed", linewidth = 0.8, alpha = 0.6) +
  geom_vline(xintercept = rec_window_days, color = protocol_cols[["Recommended"]],
             linetype = "dashed", linewidth = 0.8, alpha = 0.6) +
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = lowerCI, ymax = upperCI, fill = protocol), linetype = 2, alpha = 0.2) +
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
  # ggtitle("a. EIA infection")
  ggtitle(paste0(panel_letters[1], outcome_label))

plot_risk

# Risk ratio (Delayed / Recommended): median and 95% CI across bootstrap samples
rr_all <- df_rr_boot %>%
  group_by(surv_time_censoring) %>%
  summarise(median  = median(rr),
            lowerCI = quantile(rr, probs = 0.025),
            upperCI = quantile(rr, probs = 0.975),
            .groups = "drop")

plot_rr <- rr_all %>%
  ggplot(aes(x = surv_time_censoring, y = median)) +
  geom_vline(xintercept = del_window_days, color = protocol_cols[["Delayed"]],
             linetype = "dashed", linewidth = 0.8, alpha = 0.6) +
  geom_vline(xintercept = rec_window_days, color = protocol_cols[["Recommended"]],
             linetype = "dashed", linewidth = 0.8, alpha = 0.6) +
  geom_line(color = "#1f77b4", linewidth = 1.2) +
  geom_ribbon(aes(ymin = lowerCI, ymax = upperCI), fill = "#1f77b4", alpha = 0.2, colour = NA) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", linewidth = 0.8) +
  labs(x = "Days Since Birth", y = "Risk Ratio \n (Delayed / Recommended)") +
  scale_x_continuous(limits = c(0, followup_period), breaks = seq(0, followup_period, 50)) +
  # coord_cartesian zooms without dropping data, so the ribbon isn't cut where the CI exceeds 5
  coord_cartesian(ylim = c(0, 5)) +
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
  # ggtitle("c. EIA infection")
  ggtitle(paste0(panel_letters[3], outcome_label))

plot_rr

# Risk difference (Delayed - Recommended): median and 95% CI across bootstrap samples
rd_all <- df_rd_boot %>%
  group_by(surv_time_censoring) %>%
  summarise(median  = median(rd),
            lowerCI = quantile(rd, probs = 0.025),
            upperCI = quantile(rd, probs = 0.975),
            .groups = "drop")

plot_rd <- rd_all %>%
  ggplot(aes(x = surv_time_censoring, y = median)) +
  geom_vline(xintercept = del_window_days, color = protocol_cols[["Delayed"]],
             linetype = "dashed", linewidth = 0.8, alpha = 0.6) +
  geom_vline(xintercept = rec_window_days, color = protocol_cols[["Recommended"]],
             linetype = "dashed", linewidth = 0.8, alpha = 0.6) +
  geom_line(color = "#1f77b4", linewidth = 1.2) +
  geom_ribbon(aes(ymin = lowerCI, ymax = upperCI), fill = "#1f77b4", alpha = 0.2, colour = NA) +
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
  ggtitle(paste0(panel_letters[2], outcome_label))

plot_rd

library(patchwork)
figs_risk_rd    <- plot_risk + plot_rd
figs_risk_rr_rd <- plot_risk + plot_rd + plot_rr

ggsave(file.path(fig_dir, paste0(fig_name, "_risk_rd.png")),
       figs_risk_rd,
       width = 16,
       height = 6,
       units = "in",
       dpi = 300)

ggsave(file.path(fig_dir, paste0(fig_name, "_rr.png")),
       plot_rr + ggtitle(outcome_label),
       width = 8,
       height = 6,
       units = "in",
       dpi = 300)

ggsave(file.path(fig_dir, paste0(fig_name, "_risk_rr_rd.png")),
       figs_risk_rr_rd,
       width = 24,
       height = 6,
       units = "in",
       dpi = 300)

# Save the cumulative risk per protocol, RR and RD at Day 730
res_rr_rd <- bind_rows(
  res_all %>% mutate(measure = paste0("Cumrisk_", as.character(protocol))),
  rr_all %>% mutate(measure = "RR"),
  rd_all %>% mutate(measure = "RD")
) %>%
  filter(surv_time_censoring == followup_period) %>%
  dplyr::select(measure, surv_time_censoring, median, lowerCI, upperCI)

fwrite(res_rr_rd, file.path(out_dir, "res_rr_rd.csv"))
