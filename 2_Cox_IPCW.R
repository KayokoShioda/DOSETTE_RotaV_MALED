# Open data frame
df_result_1 <- data.frame()
df_result_2 <- data.frame()
df_rr_boot <- data.frame()

# ============================================================
# COX IPCW FUNCTION (ONE PROTOCOL, ONE SAMPLE)
# ============================================================
# d       = cox_data of one protocol, with a unique child id (id)
# n_denom = number of children in d (denominator of the cumulative risk)

fit_cox <- function(d, n_denom) {
  
  fml <- fm_cn
  # Country_ID can't be a covariate in single-country analyses
  if (n_distinct(d$Country_ID) < 2) fml <- update(fml, . ~ . - Country_ID)
  
  # model = TRUE so predict() works inside functions/workers
  coxph.censor <- coxph(fml, data = d, model = TRUE)
  # Calculate weighted outcome over study period
  cases <- broom::augment(coxph.censor,
                          newdata = d[d$outcome == 1, ],
                          type.predict = "expected") %>%
    dplyr::mutate(prob = exp(-.fitted),
                  pid = d[d$outcome == 1, ]$Pid)
  cases <- cases[order(cases$surv_time), ]
  
  cases$wt <- 1 / cases$prob
  # Create time data
  cases <- cases %>%
    subset(., select = c(pid, wt, surv_time))
  cases <- data.frame(surv_time = seq(1, followup_period, 1)) %>%
    left_join(., cases, by = "surv_time")
  
  # Compute cumsum
  cases$risk <- cumsum(replace_na(cases$wt, 0)) / n_denom
  # Extract the maximum risk (i.e., cumulative risk) on each event day
  res <- cases[which(cases$surv_time <= followup_period), ] %>%
    group_by(surv_time) %>%
    summarise(cumrisk = max(risk, na.rm = TRUE))
  
  # surv_time_censoring is the column name expected by 3_plot_figs_first_outcome.R
  res %>%
    arrange(surv_time) %>%
    rename(surv_time_censoring = surv_time) %>%
    as.data.frame()
}


# ============================================================
# BOOTSTRAP (PARALLELIZED)
# ============================================================

library(parallel)

nsims <- 500
seed  <- 123
ncores <- max(1, detectCores() - 1)
pids  <- unique(cox_data$Pid)

# One bootstrap iteration: resample Pids, fit both protocols, return both risk curves
boot_once <- function(i) {
  
  set.seed(seed + i)
  
  # Each draw gets a new id, so duplicated children are treated as separate children.
  # Both protocol copies of a Pid are kept together.
  boot_ids <- tibble(Pid = sample(pids, replace = TRUE)) %>%
    mutate(id = row_number())
  
  bootm <- boot_ids %>%
    inner_join(cox_data, by = "Pid", relationship = "many-to-many")
  
  boot_rec <- bootm %>% filter(protocol == "Recommended")
  res_rec  <- fit_cox(boot_rec, n_distinct(boot_rec$id)) %>%
    mutate(protocol = "Recommended", boot_id = i)
  
  boot_del <- bootm %>% filter(protocol == "Delayed")
  res_del  <- fit_cox(boot_del, n_distinct(boot_del$id)) %>%
    mutate(protocol = "Delayed", boot_id = i)
  
  list(res_rec = res_rec, res_del = res_del)
}

cl <- makeCluster(ncores)

clusterExport(cl, c("cox_data", "pids", "fit_cox", "fm_cn",
                    "followup_period", "seed"), environment())

clusterEvalQ(cl, {
  library(dplyr)
  library(tidyr)
  library(survival)
  library(broom)
  NULL
})

boot_res <- parLapply(cl, 1:nsims, boot_once)

stopCluster(cl)

# Reassemble results
df_result_1 <- bind_rows(lapply(boot_res, `[[`, "res_rec"))
df_result_2 <- bind_rows(lapply(boot_res, `[[`, "res_del"))

library(data.table)

# dir_name is set in master.R
out_dir <- file.path("../Rota-TTE/outputs", dir_name)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Covariates of the censoring model
check_missing_covariates(cox_data, all.vars(fm_cn[[3]]), out_dir, "cox_data")

fwrite(df_result_1, file.path(out_dir, "df_result_1.csv"))
fwrite(df_result_2, file.path(out_dir, "df_result_2.csv"))
