# ============================================================
# SEVERE DIARRHEA OUTCOME
# (run after 1_Cox_dataset.R, before 2_Cox_IPCW.R)
# Severe diarrhea: score >= 6 on the first day of a diarrhea episode (dint3 == 1)
# ============================================================

maled_diarrheaseverity <- read.csv(
  "../Rota-TTE/misc/Original Data/maled_diarrheaseverity.csv"
) %>%
  filter(substr(Pid, 1, 2) %in% country_choice)

# Check: among children in vax_inf_cov_data, any stooltype == "D1" <=> in maled_diarrheaseverity
no_severity_pids <- setdiff(unique(vax_inf_cov_data$Pid), maled_diarrheaseverity$Pid)

d1_check <- original_cov %>%
  filter(pid %in% vax_inf_cov_data$Pid) %>%
  group_by(pid) %>%
  summarise(any_D1 = any(stooltype == "D1", na.rm = TRUE)) %>%
  mutate(in_severity = pid %in% maled_diarrheaseverity$Pid)
print(table(any_D1 = d1_check$any_D1, in_severity = d1_check$in_severity))
print(sprintf("Any D1 but not in maled_diarrheaseverity: %d; no D1 but in maled_diarrheaseverity: %d",
              sum(d1_check$any_D1 & !d1_check$in_severity),
              sum(!d1_check$any_D1 & d1_check$in_severity)))

# Children in vax_inf_cov_data but not in maled_diarrheaseverity had no diarrhea
# (67 of 346 for c("BR", "PE")); kept, with no severe diarrhea event
print(sprintf("Pids in vax_inf_cov_data not in maled_diarrheaseverity (no diarrhea, kept), %d out of %d",
              length(no_severity_pids), n_distinct(vax_inf_cov_data$Pid)))

# ============================================================
# FIRST SEVERE DIARRHEA AGE
# ============================================================

severe_diarrhea_days <- maled_diarrheaseverity %>%
  filter(score >= 6, dint3 == 1, age >= 0) %>%
  group_by(Pid) %>%
  summarise(severe_diarrhea_age = min(age), .groups = "drop")

outcome_age <- "severe_diarrhea_age"

# Person-level data without infection_age / diarrhea_age
# Children with no data in maled_diarrheaseverity: severe_diarrhea = 0, severe_diarrhea_age = NA
person_data <- vax_inf_cov_data %>%
  mutate(sex = cafsex,
         severe_diarrhea_age = if_else(
           Pid %in% no_severity_pids, NA_real_,
           as.numeric(severe_diarrhea_days$severe_diarrhea_age[
             match(Pid, severe_diarrhea_days$Pid)])),
         severe_diarrhea = if_else(Pid %in% no_severity_pids, 0, 
                                   as.numeric(!is.na(severe_diarrhea_age))))


# ============================================================
# CLONE, CENSOR, AND DEFINE OUTCOME (same as 1_Cox_dataset.R)
# ============================================================

cox_data <- bind_rows(
  create_clone(person_data, "Recommended"),
  create_clone(person_data, "Delayed")
) %>%
  mutate(
    outcome  = as.numeric(coalesce(.data[[outcome_age]] < censoring_age, FALSE)),
    censored = as.numeric(censoring_age < followup_period & outcome == 0),
    surv_time = if_else(outcome == 1, .data[[outcome_age]], censoring_age),
    surv_time_censoring = surv_time,
    protocol = factor(protocol, levels = c("Recommended", "Delayed"))
  ) %>%
  dplyr::select(Pid, Country_ID, sex, enrollhaz, enrollwaz, wamiscore,
                protocol, censoring_age, outcome, censored,
                surv_time, surv_time_censoring)
