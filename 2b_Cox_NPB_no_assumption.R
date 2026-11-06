# ============================================================
# NONPARAMETRIC BOUNDS: NO-ASSUMPTION VERSION
# ============================================================
# Sensitivity analysis to the IPCW assumption. All censoring other than
# administrative censoring at followup_period (protocol deviation AND
# drop-out) is allowed to be arbitrarily informative:
#   lower bound: censored children never have the outcome
#   upper bound: censored children have the outcome at the time of censoring
# No censoring model is fitted.
#
# Run after 1_Cox_dataset.R (needs cox_data, person_data, followup_period,
# dir_name_bounds). Same bootstrap structure as 2_Cox_IPCW.R.
# Output (per bootstrap x day x protocol): df_bounds_1.csv (Recommended),
# df_bounds_2.csv (Delayed), columns surv_time_censoring, lo, hi, protocol, boot_id.

library(dplyr)
library(tidyr)
library(parallel)
library(data.table)

# ============================================================
# BOUNDS FUNCTION (ONE PROTOCOL, ONE SAMPLE)
# ============================================================
# d       = cox_data of one protocol, with a unique child id (id)
# n_denom = number of children in d

bound_arm <- function(d, n_denom) {
  admin <- d$outcome == 0 & d$censoring_age >= followup_period
  t <- seq(1, followup_period, 1)
  lo <- sapply(t, function(x) sum(d$outcome == 1 & d$surv_time <= x, na.rm = TRUE)) / n_denom
  hi <- sapply(t, function(x) sum(d$surv_time <= x & !admin, na.rm = TRUE)) / n_denom
  data.frame(surv_time_censoring = t, lo = lo, hi = pmin(hi, 1))
}


# ============================================================
# BOOTSTRAP (PARALLELIZED)  -- same resampling as 2_Cox_IPCW.R
# ============================================================

nsims <- 500
seed  <- 123
ncores <- max(1, detectCores() - 1)
pids  <- unique(cox_data$Pid)

boot_once <- function(i) {

  set.seed(seed + i)

  boot_ids <- tibble(Pid = sample(pids, replace = TRUE)) %>%
    mutate(id = row_number())

  bootm <- boot_ids %>%
    inner_join(cox_data, by = "Pid", relationship = "many-to-many")

  boot_rec <- bootm %>% filter(protocol == "Recommended")
  res_rec  <- bound_arm(boot_rec, n_distinct(boot_rec$id)) %>%
    mutate(protocol = "Recommended", boot_id = i)

  boot_del <- bootm %>% filter(protocol == "Delayed")
  res_del  <- bound_arm(boot_del, n_distinct(boot_del$id)) %>%
    mutate(protocol = "Delayed", boot_id = i)

  list(res_rec = res_rec, res_del = res_del)
}

cl <- makeCluster(ncores)

clusterExport(cl, c("cox_data", "pids", "bound_arm", "followup_period", "seed"),
              environment())

clusterEvalQ(cl, {
  library(dplyr)
  library(tidyr)
  NULL
})

boot_res <- parLapply(cl, 1:nsims, boot_once)

stopCluster(cl)

df_bounds_1 <- bind_rows(lapply(boot_res, `[[`, "res_rec"))
df_bounds_2 <- bind_rows(lapply(boot_res, `[[`, "res_del"))

# dir_name_bounds is set in master.R
out_dir <- file.path("../Rota-TTE/outputs", dir_name_bounds)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fwrite(df_bounds_1, file.path(out_dir, "df_bounds_1.csv"))
fwrite(df_bounds_2, file.path(out_dir, "df_bounds_2.csv"))
