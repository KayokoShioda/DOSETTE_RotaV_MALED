# ============================================================
# 30. POOLED LOGISTIC FUNCTION (ONE PROTOCOL, ONE SAMPLE)
# ============================================================
# d       = person-day data of one protocol, with a unique child id (id)
# n_denom = number of children in d (denominator of the cumulative risk)

fit_pooled <- function(d, n_denom) {
  
  fml <- censored ~ as.factor(time_bin) + sex + enrollhaz + enrollwaz + zwei + zlen + zwfl + exclusive_breastfeeding + saffev + calendar_month +
    wamiscore + Country_ID
  # Country_ID can't be a covariate in single-country analyses
  if (n_distinct(d$Country_ID) < 2) fml <- update(fml, . ~ . - Country_ID)
  
  fit <- glm(
    fml,
    data = d,
    family = binomial()
  )
  
  d$hazard_hat <- predict(fit, newdata = d, type = "response")
  d$hazard_hat <- ifelse(is.na(d$hazard_hat), 0, d$hazard_hat)
  
  d <- d %>%
    arrange(id, surv_time_censoring) %>%
    group_by(id) %>%
    mutate(prob = cumprod(1 - hazard_hat)) %>%
    ungroup()
  
  # outcome = 1 only on the event day, so this is one row per case
  cases <- d %>%
    filter(outcome == 1)
  
  data.table::setorder(cases, surv_time_censoring)
  cases$risk <- cumsum(1 / cases$prob) / n_denom
  
  cases[cases$surv_time_censoring <= followup_period, ] %>%
    group_by(surv_time_censoring) %>%
    summarise(cumrisk = max(risk, na.rm = TRUE), .groups = "drop") %>%
    complete(surv_time_censoring = 1:followup_period) %>%
    fill(cumrisk, .direction = "down") %>%
    # FLAG: days before the first event have no cumrisk after fill(), so they are NA.
    # They are set to 0 (no events yet means cumulative risk = 0), so every
    # bootstrap sample has a value on every day and the RR can be computed.
    mutate(cumrisk = replace_na(cumrisk, 0))
}


# ============================================================
# 31. BOOTSTRAP (PARALLELIZED)
# ============================================================

library(parallel)

nsims <- 500
seed  <- 123
ncores <- max(1, detectCores() - 1)
pids  <- unique(person_day_subset$Pid)

# One bootstrap iteration: resample Pids, fit both protocols, return both risk curves
boot_once <- function(i) {
  
  # Reproducible, independent draw per iteration
  set.seed(seed + i)
  
  # Resample Pids with replacement. Each draw gets a new id, so duplicated children
  # are treated as separate children. Both protocol copies of a Pid are kept together.
  boot_ids <- tibble(Pid = sample(pids, replace = TRUE)) %>%
    mutate(id = row_number())
  
  bootm <- boot_ids %>%
    inner_join(person_day_subset, by = "Pid", relationship = "many-to-many")
  
  # Recommended
  boot_rec <- bootm %>% filter(protocol == "Recommended")
  res_rec  <- fit_pooled(boot_rec, n_distinct(boot_rec$id)) %>%
    mutate(protocol = "Recommended", boot_id = i)
  
  # Delayed
  boot_del <- bootm %>% filter(protocol == "Delayed")
  res_del  <- fit_pooled(boot_del, n_distinct(boot_del$id)) %>%
    mutate(protocol = "Delayed", boot_id = i)
  
  list(res_rec = res_rec, res_del = res_del)
}

# Parallel setup
cl <- makeCluster(ncores)

# Export objects the workers need
clusterExport(cl, c("person_day_subset", "pids", "fit_pooled",
                    "followup_period", "seed"), environment())

# Load packages in workers
clusterEvalQ(cl, {
  library(dplyr)
  library(tidyr)
  library(data.table)
  NULL
})

# Run bootstrap iterations in parallel
boot_res <- parLapply(cl, 1:nsims, boot_once)

stopCluster(cl)

# Reassemble results
df_result_1 <- bind_rows(lapply(boot_res, `[[`, "res_rec"))
df_result_2 <- bind_rows(lapply(boot_res, `[[`, "res_del"))

library(dplyr)
library(tidyr)
library(data.table)

# dir_name is set in master.R
out_dir <- file.path("../Rota-TTE/outputs", dir_name)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Covariates of the pooled logistic censoring model in fit_pooled() (keep in sync)
check_missing_covariates(
  person_day_subset,
  c("time_bin", "sex", "enrollhaz", "enrollwaz", "zwei", "zlen", "zwfl",
    "exclusive_breastfeeding", "saffev", "calendar_month", "wamiscore", "Country_ID"),
  out_dir, "person_day_subset"
)

fwrite(df_result_1, file.path(out_dir, "df_result_1.csv"))
fwrite(df_result_2, file.path(out_dir, "df_result_2.csv"))
