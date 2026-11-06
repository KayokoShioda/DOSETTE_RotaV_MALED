# ============================================================
# ROTAVIRUS TARGET TRIAL EMULATION
# POOLED LOGISTIC REGRESSION DATASET
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
library(ggplot2)
library(dplyr)

# Outcome option: "eia" or "eia_diarrhea"
# outcome_choice <- "eia_diarrhea" # now set in master.R
outcome_age <- if (outcome_choice == "eia") "infection_age" else "diarrhea_age"

# Country option: "BR", "PE", or c("BR", "PE")
# country_choice <- c("BR", "PE") # now set in master.R


# ============================================================
# 2. LOAD VACCINATION DATA
# ============================================================

original_vax <- read.csv(
  "../Rota-TTE/misc/Original Data/vax_and_titer_wdates.csv"
) %>%
  filter(
    Country_ID %in% country_choice
  )


# ============================================================
# 3. LOAD INFECTION + COVARIATE DATA
# ============================================================

original_cov <- read_csv(
  "../Rota-TTE/misc/Original Data/Tac_Sep2018.csv"
) %>% subset(country_id %in% country_choice)


# ============================================================
# 4. ALL EIA-CONFIRMED ROTAVIRUS INFECTIONS
# ============================================================
#EIA-confirmed infection: rota_micro == 1
#EIA-confirmed dirrhea: rota_micro == 1 & subset_cov$stooltype == "D1"
#Keep ALL infection dates for each participant.

eia_dates <- original_cov %>%
  filter(
    rota_micro == 1
  ) %>%
  mutate(
    date = as.Date(
      date,
      format = "%m/%d/%Y"
    )
  ) %>%
  filter(
    !is.na(date)
  ) %>%
  arrange(
    pid,
    date
  ) %>%
  distinct(
    pid,
    date
  )


# ============================================================
# 5. ALL EIA-CONFIRMED ROTAVIRUS DIARRHEA
# ============================================================
#EIA-confirmed rotavirus diarrhea: rota_micro == 1 AND stooltype == "D1"
#Keep ALL diarrhea dates for each participant.

eia_diarrhea_dates <- original_cov %>%
  filter(
    rota_micro == 1,
    stooltype == "D1"
  ) %>%
  mutate(
    date = as.Date(
      date,
      format = "%m/%d/%Y"
    )
  ) %>%
  filter(
    !is.na(date)
  ) %>%
  arrange(
    pid,
    date
  ) %>%
  distinct(
    pid,
    date
  )


# ============================================================
# 6. CREATE BASELINE PARTICIPANT DATA
# ============================================================

baseline_cov <- original_cov %>%
  rename(
    Pid = pid
  ) %>%
  arrange(
    Pid,
    date
  ) %>%
  group_by(
    Pid
  ) %>%
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
# 7. CREATE TIME-VARYING COVARIATES
# ============================================================
#zwei, zlen, zwfl, score, and ebf

tv_covariates <- original_cov %>%
  rename(
    Pid = pid
  ) %>%
  mutate(
    
    # Convert measurement date
    date = as.Date(
      date,
      format = "%m/%d/%Y"
    ),
    
    # Convert date of birth
    dob = as.Date(
      dob,
      format = "%m/%d/%Y"
    ),
    
    birth_date = dob,
    
    age_day = as.numeric(
      date - birth_date
    )
    
  ) %>%
  
  filter(
    age_day >= 0
  ) %>%
  
  dplyr::select(
    Pid,
    age_day,
    zwei,
    zlen,
    zwfl,
    score,
    ebf
  ) %>%
  
  arrange(
    Pid,
    age_day
  )

# KJ: add fever
maled_fever <- read.csv("../Rota-TTE/misc/Original Data/maled_fever.csv")

fever <- maled_fever %>% dplyr::select(Pid, age_day = age, saffev) %>%
  filter(substr(Pid, 1, 2) %in% country_choice)

# One row per Pid-day on both sides, so the joins below cannot duplicate person-days
stopifnot(
  "Duplicate Pid-age_day in tv_covariates" = !anyDuplicated(tv_covariates[c("Pid", "age_day")]),
  "Duplicate Pid-age_day in fever" = !anyDuplicated(fever[c("Pid", "age_day")])
)

# 0 kids are dropped here: every Pid/age_day row of tv_covariates has a match in fever (not vice versa),
# so inner_join keeps all rows and kids
tv_covariates <- tv_covariates %>%
  inner_join(fever %>% dplyr::select(Pid, age_day, saffev), 
             by=c("Pid", "age_day"))

# ============================================================
# 8. MERGE VACCINATION + BASELINE INFORMATION
# ============================================================

stopifnot(
  "Pid has >1 distinct dose record in original_vax" =
    !anyDuplicated(distinct(original_vax, Pid, rotadate1, rotaaged1, rotadate2, rotaaged2)$Pid)
)
vax_inf_cov_data <- original_vax %>%
  # KJ: Keep rotavirus vaccine data only.
  # Each `Pid` has two rows with different `VRFSID` and `datecol` values,
  # so de-duplicate on the unique rotadate1, rotaaged1, rotadate2, and rotaaged2.
  # Note there are 100 individuals without covariate data that are dropped here.
  # 100 kids are dropped here because they are in original_vax but not in original_cov;
  # the 12 kids in original_cov but not in original_vax are also not carried forward
  # because inner_join keeps only Pids present in both.
  distinct(Pid, rotadate1, rotaaged1, rotadate2, rotaaged2) %>%
  inner_join(
    baseline_cov,
    by = "Pid"
  )

# ============================================================
# 9. CREATE BIRTH DATE
# ============================================================

vax_inf_cov_data <- vax_inf_cov_data %>%
  mutate(
    
    birth_date = as.Date(
      dob,
      format = "%m/%d/%Y"
    )
    
  )


# ============================================================
# 10. LOAD DROPOUT DATA
# ============================================================

dropout <- read.csv(
  "../Rota-TTE/misc/Original Data/dropout.csv"
) %>%
  filter(substr(Pid, 1, 2) %in% country_choice)


# ============================================================
# 11. MERGE DROPOUT DATA
# ============================================================

stopifnot("Duplicate Pid in dropout" = !anyDuplicated(dropout$Pid))

# 0 kids are dropped here: every Pid in vax_inf_cov_data has a match in dropout,
# so inner_join keeps all kids
vax_inf_cov_data <- vax_inf_cov_data %>%
  inner_join(
    dropout,
    by = "Pid"
  )


# ============================================================
# 12. REMOVE DOSE 1 RECORDED ON DAY 0
# ============================================================

# 1 kid is dropped here because rotaaged1 == 0, a presumed data-entry error
vax_inf_cov_data <- vax_inf_cov_data %>%
  filter(
    is.na(rotaaged1) |
      rotaaged1 != 0
  )

# ============================================================
# 13. TARGET TRIAL VACCINATION WINDOWS
# ============================================================

#Recommended 
rec_dose_1_start <- 7 * 7 # 7 weeks
rec_dose_1_end   <- 9 * 7 - 1 # <9 weeks

rec_dose_2_start <- 16 * 7 # 16 weeks
rec_dose_2_end   <- 18 * 7 -1 # <18 weeks


#Delayed 
del_dose_1_start <- 9 * 7 # 9 weeks
del_dose_1_end   <- 15 * 7 - 1 # <15 weeks

del_dose_2_start <- 18 * 7 # 18 weeks 
del_dose_2_end   <- 25 * 7 # 25 weeks


#Administrative follow-up
followup_period <- 730


# ============================================================
# 14. CREATE CLONE FUNCTION
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
  
  
  # ----------------------------------------------------------
  # Follow-up ends at protocol censoring
  # ----------------------------------------------------------
  
  clone <- clone %>%
    mutate(
      
      age_end = censoring_age
      
    )
  
  
  return(clone)
}



# ============================================================
# 15. CREATE FIRST EIA INFECTION AND DIARRHEA AGES
# ============================================================
# Convert the first EIA-confirmed infection and EIA-confirmed diarrhea date into age in days.

infection_days <- eia_dates %>%
  rename(
    Pid = pid
  ) %>%
  left_join(
    baseline_cov %>%
      dplyr::select(
        Pid,
        dob
      ),
    by = "Pid"
  ) %>%
  mutate(
    
    birth_date = as.Date(
      dob,
      format = "%m/%d/%Y"
    ),
    
    infection_age = as.numeric(
      date - birth_date
    )
    
  ) %>%
  filter(
    !is.na(infection_age),
    infection_age >= 0
  ) %>%
  dplyr::select(
    Pid,
    infection_age
  ) %>%
  distinct() %>% 
  group_by(Pid) %>%
  summarise(infection_age=min(infection_age)) %>%
  ungroup()

diarrhea_days <- eia_diarrhea_dates %>%
  rename(
    Pid = pid
  ) %>%
  left_join(
    baseline_cov %>%
      dplyr::select(
        Pid,
        dob
      ),
    by = "Pid"
  ) %>%
  mutate(
    
    birth_date = as.Date(
      dob,
      format = "%m/%d/%Y"
    ),
    
    diarrhea_age = as.numeric(
      date - birth_date
    )
    
  ) %>%
  filter(
    !is.na(diarrhea_age),
    diarrhea_age >= 0
  ) %>%
  dplyr::select(
    Pid,
    diarrhea_age
  ) %>%
  distinct() %>% 
  group_by(Pid) %>%
  summarise(diarrhea_age=min(diarrhea_age)) %>%
  ungroup()


# ============================================================
# 16. START PERSON-LEVEL DATASET AND MERGE EIA INFECTION AND DIARRHEA
# ============================================================
# One row per participant (not yet cloned).

person_day <- vax_inf_cov_data %>%
  left_join(
    infection_days %>%
      mutate(
        eia_event = 1
      ),
    by = c(
      "Pid"
    )
  ) %>%
  left_join(
    diarrhea_days %>%
      mutate(
        eia_diarrhea_event = 1
      ),
    by = c(
      "Pid"
    )
  ) %>%
  mutate(
    
    eia_event = ifelse(
      is.na(eia_event),
      0,
      eia_event
    ),
    
    eia_diarrhea_event = ifelse(
      is.na(eia_diarrhea_event),
      0,
      eia_diarrhea_event
    )
    
  )


# ============================================================
# 17. CREATE PERSON-DAY DATASET (UNCLONED)
# ============================================================
# TIME ZERO = BIRTH
# Follow-up here continues until:
#   administrative end at day 730
#   outcome
# Protocol-specific censoring (deviation / dropout) is applied after cloning.

person_day <- person_day %>%
  rowwise() %>%
  mutate(
    
    age_day = list(
      seq(
        1,
        min(
          followup_period,
          .data[[outcome_age]]
          , na.rm = T),
        by = 1
      )
    )
    
  ) %>%
  ungroup() %>%
  unnest(age_day)


# ============================================================
# 18. MERGE TIME-VARYING COVARIATES
# ============================================================

person_day <- person_day %>%
  left_join(
    tv_covariates,
    by = c(
      "Pid",
      "age_day"
    )
  ) %>%
  group_by(
    Pid
  ) %>%
  arrange(
    age_day,
    .by_group = TRUE
  ) %>%
  fill(
    zwei,
    zlen,
    zwfl,
    score,
    ebf,
    saffev,
    .direction = "down"
  ) %>%
  ungroup()


# ============================================================
# 19. EXCLUSIVE BREASTFEEDING
# ============================================================

# For initial days where ebf is NA, use the child's first available ebf
person_day <- person_day %>%
  group_by(
    Pid
  ) %>%
  mutate(
    
    exclusive_breastfeeding =
      ifelse(
        
        age_day <= coalesce(ebf, first(na.omit(ebf))),
        
        1,
        0
      )
    
  ) %>%
  ungroup()


# ============================================================
# 21. CALENDAR DATE
# ============================================================
# age_day = time since birth
#
# calendar_date = actual calendar date for that person-day

person_day <- person_day %>%
  mutate(
    
    calendar_date =
      birth_date + age_day
    
  )


# ============================================================
# 22. CALENDAR MONTH
# ============================================================
# Calendar month is time-varying because it changes
# as the child moves through follow-up.
#
# 01 = January
# 02 = February
# ...
# 12 = December

person_day <- person_day %>%
  mutate(
    
    calendar_month =
      factor(
        format(
          calendar_date,
          "%m"
        ),
        levels = sprintf(
          "%02d",
          1:12
        )
      )
    
  )


# ============================================================
# 23. SUBSET VARIABLES, CHECK MISSINGNESS, AND FILL LEADING NAs
# ============================================================

# 1. Select variables (infection_age is kept temporarily to derive
#    the censoring and outcome indicators after cloning)
person_day_subset <- person_day %>% dplyr::select(Pid, 
                                                  Country_ID, 
                                                  wamiscore, 
                                                  sex = cafsex, 
                                                  enrollhaz, 
                                                  enrollwaz, 
                                                  zwei, 
                                                  zlen, 
                                                  zwfl, 
                                                  saffev,
                                                  exclusive_breastfeeding, 
                                                  calendar_month,
                                                  all_of(outcome_age),
                                                  surv_time_censoring = age_day) 

# 2. Check that all covariates have the same NA count within each Pid,
# and, for every covariate, all NAs fall in the first n rows 

na_by_pid <- person_day_subset %>% 
  group_by(Pid) %>%
  summarise(across(c(sex, enrollhaz, enrollwaz, zwei, zlen, zwfl, saffev, exclusive_breastfeeding),
                   ~ sum(is.na(.x))),
            .groups = "drop")

covs <- c("zwei", "zlen", "zwfl", "saffev")

na_check <- person_day_subset %>%
  arrange(Pid, surv_time_censoring) %>%
  group_by(Pid) %>%
  summarise(
    same_n_na  = n_distinct(sapply(pick(all_of(covs)), function(x) sum(is.na(x)))) == 1,
    leading_na = all(sapply(pick(all_of(covs)), function(x) all(is.na(x[seq_len(sum(is.na(x)))])))),
    .groups = "drop")

# Within each Pid, leading NAs (rows before the first non-NA value) take that first non-NA value. 
# saffev: all leading NAs = 0.
person_day_subset <- person_day_subset %>%
  arrange(Pid, surv_time_censoring) %>%
  group_by(Pid) %>%
  mutate(across(all_of(setdiff(covs, "saffev")), ~ {
    i <- which(!is.na(.x))[1]
    if (is.na(i)) .x else replace(.x, seq_len(i - 1), .x[i])
  }),
  saffev = {
    i <- which(!is.na(saffev))[1]
    if (is.na(i)) saffev else replace(saffev, seq_len(i - 1), 0)
  }) %>%
  ungroup() %>%
  # time_bin does not depend on the bootstrap sample, so compute it once
  mutate(time_bin = factor(ceiling(surv_time_censoring / 7)))


# ============================================================
# 24. CREATE RECOMMENDED CLONE
# ============================================================

clone_recommended <- create_clone(
  vax_inf_cov_data,
  "Recommended"
)


# ============================================================
# 25. CREATE DELAYED CLONE
# ============================================================

clone_delayed <- create_clone(
  vax_inf_cov_data,
  "Delayed"
)


# ============================================================
# 26. COMBINE CLONES
# ============================================================

clones <- bind_rows(
  clone_recommended,
  clone_delayed
)


# ============================================================
# 27. APPLY CLONES TO PERSON-DAY DATASET
# ============================================================
# Each person-day row is duplicated once per protocol, then follow-up is
# truncated at the protocol-specific censoring age. Rows with NA censoring_age
# are kept, matching the previous na.rm = TRUE behaviour in the expansion.

person_day_subset <- person_day_subset %>%
  inner_join(
    clones %>% dplyr::select(Pid, protocol, censoring_age),
    by = "Pid",
    relationship = "many-to-many"
  ) %>%
  filter(
    is.na(censoring_age) | surv_time_censoring <= censoring_age
  ) %>%
  arrange(Pid, protocol, surv_time_censoring)


# ============================================================
# 28. CENSORING INDICATOR
# ============================================================
# censored = 1 on the last day of follow-up for anyone
# without an infection before censoring, EXCLUDING day-730 completers
# (those who reach the administrative end of follow-up get censored = 0)

person_day_subset <- person_day_subset %>%
  group_by(Pid, protocol) %>%
  mutate(
    censored = as.numeric(
      surv_time_censoring == max(surv_time_censoring) &
        surv_time_censoring < followup_period &
        !coalesce(.data[[outcome_age]] < censoring_age, FALSE)
    )
  ) %>%
  ungroup()


# ============================================================
# 29. OUTCOME: INFECTION BEFORE CENSORING
# ============================================================
# outcome = 1 only on the event day (surv_time_censoring == outcome age) of a
# child-protocol whose first outcome occurs before censoring_age (strict <: an event
# on the censoring day, including window-end days, counts as censored, by design)

person_day_subset <- person_day_subset %>%
  mutate(
    outcome = as.numeric(
      coalesce(.data[[outcome_age]] < censoring_age &
                 surv_time_censoring == .data[[outcome_age]], FALSE)
    )
  ) %>%
  dplyr::select(Pid, 
                Country_ID, 
                wamiscore, 
                sex, 
                enrollhaz, 
                enrollwaz, 
                protocol, 
                censoring_age, 
                zwei, 
                zlen, 
                zwfl, 
                saffev,
                exclusive_breastfeeding, 
                censored, 
                outcome,
                calendar_month,
                surv_time_censoring,
                time_bin)
