# ============================================================
# ROTAVIRUS TARGET TRIAL EMULATION
# Cox REGRESSION DATASET
# (censored at first EIA-confirmed infection or
#  first EIA-confirmed diarrhea, set by outcome_choice)
# ============================================================


# ============================================================
# 1. PACKAGES
# ============================================================

library(tidyverse)
library(survival)
library(broom)
library(gtsummary)

# Outcome option: "eia" or "eia_diarrhea" (outcome_choice set in master.R)
outcome_age <- if (outcome_choice == "eia") "infection_age" else "diarrhea_age"

# Country option: "BR", "PE", or c("BR", "PE") (country_choice set in master.R)


# ============================================================
# 2. LOAD DATA
# ============================================================

original_vax <- read.csv(
  "../Rota-TTE/misc/Original Data/vax_and_titer_wdates.csv"
) %>%
  filter(Country_ID %in% country_choice)

original_cov <- read_csv(
  "../Rota-TTE/misc/Original Data/Tac_Sep2018.csv"
) %>% subset(country_id %in% country_choice)

dropout <- read.csv(
  "../Rota-TTE/misc/Original Data/dropout.csv"
) %>%
  filter(substr(Pid, 1, 2) %in% country_choice)


# ============================================================
# 3. FIRST EIA-CONFIRMED INFECTION / DIARRHEA DATES
# ============================================================
# EIA-confirmed infection: rota_micro == 1
# EIA-confirmed diarrhea: rota_micro == 1 & stooltype == "D1"

get_dates <- function(data) {
  data %>%
    mutate(date = as.Date(date, format = "%m/%d/%Y")) %>%
    filter(!is.na(date)) %>%
    distinct(pid, date)
}

eia_dates <- original_cov %>%
  filter(rota_micro == 1) %>%
  get_dates()

eia_diarrhea_dates <- original_cov %>%
  filter(rota_micro == 1, stooltype == "D1") %>%
  get_dates()


# ============================================================
# 4. BASELINE PARTICIPANT DATA
# ============================================================

baseline_cov <- original_cov %>%
  rename(Pid = pid) %>%
  arrange(Pid, date) %>%
  group_by(Pid) %>%
  summarise(
    Country_ID = first(country_id),
    cafsex = first(cafsex),
    dob = first(dob),
    enrollhaz = first(enrollhaz),
    enrollwaz = first(enrollwaz),
    wamiscore = first(wamiscore),
    .groups = "drop"
  )


# ============================================================
# 5. MERGE VACCINATION + BASELINE + DROPOUT
# ============================================================

# Each `Pid` has two rows; de-duplicate on the dose dates/ages.
vax_doses <- original_vax %>%
  distinct(Pid, rotadate1, rotaaged1, rotadate2, rotaaged2)

# stopifnot(
#   "Pid has >1 distinct dose record in original_vax" = !anyDuplicated(vax_doses$Pid),
#   "Duplicate Pid in dropout" = !anyDuplicated(dropout$Pid)
# ) # checked

vax_inf_cov_data <- vax_doses %>%
  # Pids missing from original_cov (100) or original_vax (12) are dropped by inner_join.
  inner_join(baseline_cov, by = "Pid") %>%
  mutate(birth_date = as.Date(dob, format = "%m/%d/%Y")) %>%
  inner_join(dropout, by = "Pid") %>%
  # 1 kid is dropped here because rotaaged1 == 0, a presumed data-entry error
  filter(is.na(rotaaged1) | rotaaged1 != 0)


# ============================================================
# 6. TARGET TRIAL VACCINATION WINDOWS
# ============================================================

# Recommended
rec_dose_1_start <- 7 * 7 # 7 weeks
rec_dose_1_end   <- 9 * 7 - 1 # <9 weeks

rec_dose_2_start <- 16 * 7 # 16 weeks
rec_dose_2_end   <- 18 * 7 - 1 # <18 weeks

# Delayed
del_dose_1_start <- 9 * 7 # 9 weeks
del_dose_1_end   <- 15 * 7 - 1 # <15 weeks

del_dose_2_start <- 18 * 7 # 18 weeks
del_dose_2_end   <- 25 * 7 # 25 weeks

# Administrative follow-up
followup_period <- 730

# Formula for censoring weight
fm_cn <- as.formula("Surv(surv_time_censoring, censored) ~ sex + enrollhaz + enrollwaz + wamiscore + Country_ID")


# ============================================================
# 7. CREATE CLONE FUNCTION
# ============================================================

create_clone <- function(
    data,
    protocol_number
) {
  
  if (protocol_number == "Recommended") {
    
    # -----------------------------
    # Recommended protocol
    # -----------------------------
    
    clone <- data %>%
      mutate(
        
        protocol = "Recommended",
        
        censoring_age = case_when(
          
          # Criteria 1: dropped out before Dose 1 ended and did not receive dose 1
          dropage < followup_period &
            (dropage < rec_dose_1_end & is.na(rotaaged1)) ~ dropage,
          
          # Criteria 2: didn't receive dose 1
          is.na(rotaaged1) ~ rec_dose_1_end,
          
          # Criteria 3: Got first dose before protocol
          rotaaged1 < rec_dose_1_start ~ rotaaged1,
          
          # Criteria 4: Got first dose after protocol
          rotaaged1 > rec_dose_1_end ~ rec_dose_1_end,
          
          # Criteria 5: Dropped out after dose 1 but before dose 2
          dropage < followup_period &
            (dropage >= rotaaged1) &
            (dropage <= rotaaged2 | is.na(rotaaged2)) &
            (rotaaged1 >= rec_dose_1_start &
               rotaaged1 <= rec_dose_1_end) &
            (dropage < rec_dose_2_end) ~ dropage,
          
          # Criteria 6: Did not receive dose 2
          is.na(rotaaged2) ~ rec_dose_2_end,
          
          # Criteria 7: Second dose within 28 days of first dose
          (rotaaged2 - rotaaged1 < 28) ~ rotaaged2,
          
          # Criteria 8: Got second dose before protocol
          rotaaged2 < rec_dose_2_start ~ rotaaged2,
          
          # Criteria 9: Got second dose after protocol
          rotaaged2 > rec_dose_2_end ~ rec_dose_2_end,
          
          # Criteria 10: Dropped out after dose 2
          dropage < followup_period &
            (dropage > rotaaged2) &
            (rotaaged1 >= rec_dose_1_start &
               rotaaged1 <= rec_dose_1_end) &
            (rotaaged2 >= rec_dose_2_start &
               rotaaged2 <= rec_dose_2_end) ~ dropage,
          
          # Criteria 11: Followed both dosing protocols
          (rotaaged1 >= rec_dose_1_start) &
            (rotaaged1 <= rec_dose_1_end) &
            (rotaaged2 >= rec_dose_2_start) &
            (rotaaged2 <= rec_dose_2_end) &
            (rotaaged2 - rotaaged1 >= 28) ~ followup_period
        )
      )
    
  } else {
    
    # -----------------------------
    # Delayed protocol
    # -----------------------------
    
    clone <- data %>%
      mutate(
        
        protocol = "Delayed",
        
        censoring_age = case_when(
          
          # Criteria 1: dropped out before Dose 1 ended and did not receive dose 1
          dropage < followup_period &
            (dropage < del_dose_1_end & is.na(rotaaged1)) ~ dropage,
          
          # Criteria 2: didn't receive dose 1
          is.na(rotaaged1) ~ del_dose_1_end,
          
          # Criteria 3: Got first dose before protocol
          rotaaged1 < del_dose_1_start ~ rotaaged1,
          
          # Criteria 4: Got first dose after protocol
          rotaaged1 > del_dose_1_end ~ del_dose_1_end,
          
          # Criteria 5: Dropped out after dose 1 but before dose 2
          dropage < followup_period &
            (dropage >= rotaaged1) &
            (dropage <= rotaaged2 | is.na(rotaaged2)) &
            (rotaaged1 >= del_dose_1_start &
               rotaaged1 <= del_dose_1_end) &
            (dropage < del_dose_2_end) ~ dropage,
          
          # Criteria 6: Did not receive dose 2
          is.na(rotaaged2) ~ del_dose_2_end,
          
          # Criteria 7: Second dose within 28 days of first dose
          (rotaaged2 - rotaaged1 < 28) ~ rotaaged2,
          
          # Criteria 8: Got second dose before protocol
          rotaaged2 < del_dose_2_start ~ rotaaged2,
          
          # Criteria 9: Got second dose after protocol
          rotaaged2 > del_dose_2_end ~ del_dose_2_end,
          
          # Criteria 10: Dropped out after dose 2
          dropage < followup_period &
            (dropage > rotaaged2) &
            (rotaaged1 >= del_dose_1_start &
               rotaaged1 <= del_dose_1_end) &
            (rotaaged2 >= del_dose_2_start &
               rotaaged2 <= del_dose_2_end) ~ dropage,
          
          # Criteria 11: Followed both dosing protocols
          (rotaaged1 >= del_dose_1_start) &
            (rotaaged1 <= del_dose_1_end) &
            (rotaaged2 >= del_dose_2_start) &
            (rotaaged2 <= del_dose_2_end) &
            (rotaaged2 - rotaaged1 >= 28) ~ followup_period
        )
      )
  }
  
  return(clone)
}


# ============================================================
# 8. FIRST EIA INFECTION AND DIARRHEA AGES
# ============================================================

first_age <- function(dates, age_name) {
  dates %>%
    rename(Pid = pid) %>%
    inner_join(baseline_cov %>% dplyr::select(Pid, dob), by = "Pid") %>%
    mutate(age = as.numeric(date - as.Date(dob, format = "%m/%d/%Y"))) %>%
    filter(!is.na(age), age >= 0) %>%
    group_by(Pid) %>%
    summarise(age = min(age), .groups = "drop") %>%
    rename(!!age_name := age)
}

infection_days <- first_age(eia_dates, "infection_age")
diarrhea_days  <- first_age(eia_diarrhea_dates, "diarrhea_age")


# ============================================================
# DESCRIPTIVE STATISTICS (BEFORE CLONING AND CENSORING)
# ============================================================
# One row per participant from vax_inf_cov_data.
# Timing categories among two-dose recipients: "Recommended" = both doses in the
# recommended windows, "Delayed" = both doses in the delayed windows,
# "Other" = any other combination (kept visible rather than assumed).

desc_data <- vax_inf_cov_data |>
  distinct(Pid, .keep_all = TRUE) |>
  mutate(
    sex = factor(cafsex),
    age_dose1 = rotaaged1,
    n_doses = factor(
      case_when(
        !is.na(rotaaged1) & !is.na(rotaaged2) ~ "Two",
        !is.na(rotaaged1) | !is.na(rotaaged2) ~ "One",
        TRUE ~ "Zero"
      ),
      levels = c("Two", "One", "Zero")
    ),
    # Non-two-dose recipients shown as an explicit "Missing" level so they count in the denominator
    timing = case_when(
      n_doses != "Two" ~ "Missing",
      between(rotaaged1, rec_dose_1_start, rec_dose_1_end) &
        between(rotaaged2, rec_dose_2_start, rec_dose_2_end) &
        (rotaaged2 - rotaaged1) >= 28 ~ "Recommended timing",
      between(rotaaged1, del_dose_1_start, del_dose_1_end) &
        between(rotaaged2, del_dose_2_start, del_dose_2_end) &
        (rotaaged2 - rotaaged1) >= 28 ~ "Delayed timing",
      TRUE ~ "Other timing"
    ) |>
      factor(levels = c("Recommended timing", "Delayed timing", "Other timing", "Missing")),
    eia_infection = Pid %in% infection_days$Pid[infection_days$infection_age < followup_period],
    eia_diarrhea  = Pid %in% diarrhea_days$Pid[diarrhea_days$diarrhea_age < followup_period]
  )

desc_tbl <- desc_data |>
  dplyr::select(sex, age_dose1, enrollhaz, enrollwaz, wamiscore,
                n_doses, timing, eia_infection, eia_diarrhea) |>
  tbl_summary(
    type = list(c(age_dose1, enrollhaz, enrollwaz, wamiscore) ~ "continuous",
                c(eia_infection, eia_diarrhea) ~ "dichotomous"),
    statistic = list(all_continuous() ~ "{mean} ({sd})",
                     all_categorical() ~ "{n} ({p}%)"),
    digits = list(all_continuous() ~ 1,
                  all_categorical() ~ c(0, 1)),
    value = list(eia_infection ~ TRUE, eia_diarrhea ~ TRUE),
    label = list(
      sex ~ "Sex",
      age_dose1 ~ "Age at first dose (days), among those with a first dose",
      enrollhaz ~ "Height at birth z-score",
      enrollwaz ~ "Weight at birth z-score",
      wamiscore ~ "WAMI index",
      n_doses ~ "Number of rotavirus vaccine doses",
      timing ~ "Timing of two-dose rotavirus vaccine",
      eia_infection ~ "At least one EIA-confirmed infection during two-year follow-up",
      eia_diarrhea ~ "At least one EIA-confirmed diarrhea during two-year follow-up"
    ),
    missing_text = "Missing"
  )
desc_tbl

# Save descriptive table (EIA only; folder depends on country, not outcome)
if (outcome_choice == "eia") {
  desc_dir <- file.path("../Rota-TTE/suppl_outputs", paste(country_choice, collapse = "_"))
  dir.create(desc_dir, recursive = TRUE, showWarnings = FALSE)
  gt::gtsave(gtsummary::as_gt(desc_tbl), file.path(desc_dir, "desc_table.html"))
}


# ============================================================
# 9. PERSON-LEVEL DATASET (ONE ROW PER PARTICIPANT)
# ============================================================

person_data <- vax_inf_cov_data %>%
  left_join(infection_days, by = "Pid") %>%
  left_join(diarrhea_days, by = "Pid") %>%
  mutate(sex = cafsex)


# ============================================================
# 10. CLONE, CENSOR, AND DEFINE OUTCOME
# ============================================================
# One row per participant per protocol.
# outcome = 1 if event occurs before protocol censoring_age (strict <: an event on the
# censoring day, including window-end days, counts as censored, by design).
# censored = 1 if follow-up ends at protocol deviation/dropout
# (< followup_period) without an outcome.
# surv_time = surv_time_censoring = age at outcome, or censoring_age.

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

# Check
# sum(!is.na(person_data$infection_age) &  
#       person_data$Country_ID %in% c("BR")  &                   
#       person_data$rotaaged1 >= rec_dose_1_start & person_data$rotaaged1 <= rec_dose_1_end &     
#       person_data$rotaaged2 >= rec_dose_2_start & person_data$rotaaged2 <= rec_dose_2_end & 
#       (person_data$rotaaged2 - person_data$rotaaged1 >= 28) , na.rm = T)
# 5