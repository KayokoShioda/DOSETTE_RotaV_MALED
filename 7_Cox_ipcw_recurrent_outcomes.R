# ============================================================
# 30. COX IPCW FUNCTION (ONE PROTOCOL, ONE SAMPLE)
# ============================================================
# cd = cox_data_recurrent_outcomes of one protocol and sample (one row per child, with a
#      unique child id), used both to fit the censoring model (fm_cn, from
#      6_Cox_dataset_recurrent_outcomes.R) and to compute each child's IPC weight

# Cumulative-count column follows outcome_choice (set in master.R)
cum_t_var <- if (outcome_choice == "eia") "cum_eia_t" else "cum_eia_diarrhea_t"

fit_cox <- function(cd) {
  
  fml <- fm_cn
  # Country_ID can't be a covariate in single-country analyses
  if (n_distinct(cd$Country_ID) < 2) fml <- update(fml, . ~ . - Country_ID)
  
  fit <- coxph(fml, data = cd, model = TRUE)
  
  # Weight at the end of follow-up: 1 / K(t_end | x) = 1 / exp(-H(t_end | x))
  cd <- broom::augment(fit,
                       newdata = cd,
                       type.predict = "expected") %>%
    dplyr::mutate(prob = exp(-.fitted))
  
  cd$w <- 1/cd$prob

  cd
}

# ============================================================
# 31. BOOTSTRAP (PARALLELIZED)
# ============================================================

library(parallel)

nsims <- 500
seed  <- 123
ncores <- max(1, detectCores() - 1)
pids  <- unique(cox_data_recurrent_outcomes$Pid)

# One bootstrap iteration: resample Pids, fit both protocols, return the IRR
boot_once <- function(i) {
  
  # Reproducible, independent draw per iteration
  set.seed(seed + i)
  
  # Resample Pids with replacement. Each draw gets a new id, so duplicated children
  # are treated as separate children. Both protocol copies of a Pid are kept together.
  boot_ids <- tibble(Pid = sample(pids, replace = TRUE)) %>%
    mutate(id = row_number())
  
  bootc <- boot_ids %>%
    inner_join(cox_data_recurrent_outcomes, by = "Pid", relationship = "many-to-many")
  
  res_rec <- fit_cox(filter(bootc, protocol == "Recommended")) %>% mutate(boot_id = i)
  res_del <- fit_cox(filter(bootc, protocol == "Delayed")) %>% mutate(boot_id = i)
  
  tryCatch(fit_irr(rbind(res_rec, res_del)), error = function(e) c(irr_nb = NA_real_))
}

# Weighted rate ratio (Delayed vs Recommended) by naive NB + IPW. CIs come from the
# bootstrap, which accounts for within-child correlation.
fit_irr <- function(d) {
  # Naive NB + IPW: one row per child-protocol, adherent children only (never censored).
  # Outcome = total episodes, offset = total follow-up days, weight = IPC weight at the
  # end of follow-up, 1 / K(t_end | x) = 1 / exp(-H(t_end | x)) from the Cox censoring model.
  nb_d <- d %>%
    filter(censored == 0) %>%
    transmute(
      ttl_events = .data[[cum_t_var]],
      ttl_time   = surv_time_censoring,
      w,
      protocol   = factor(protocol, levels = c("Recommended", "Delayed"))
    )
  
  irr_nb <- tryCatch({
    nb_model <- MASS::glm.nb(ttl_events ~ protocol + offset(log(ttl_time)),
                             data = nb_d,
                             weights = w)
    unname(exp(coef(nb_model)["protocolDelayed"]))
  }, error = function(e) NA_real_)
  
  c(irr_nb = irr_nb)
}

# Parallel bootstrap (PSOCK cluster works on all platforms)
cl <- makeCluster(ncores)
clusterExport(cl, c("cox_data_recurrent_outcomes", "fm_cn", "pids",
                    "fit_cox", "fit_irr", "seed", "cum_t_var"))
clusterEvalQ(cl, {
  library(dplyr)
  library(survival)
  NULL
})
# nsims x 1 matrix: column irr_nb (naive NB + IPW)
boot_irr <- do.call(rbind, parLapply(cl, seq_len(nsims), boot_once))
stopCluster(cl)

library(dplyr)
library(data.table)

# dir_name is set in master.R
out_dir <- file.path("../Rota-TTE/outputs", dir_name)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Covariates of the censoring model
check_missing_covariates(cox_data_recurrent_outcomes, all.vars(fm_cn[[3]]), out_dir,
                         "cox_data_recurrent_outcomes")

fwrite(data.table(boot_id = seq_len(nsims), boot_irr),
       file.path(out_dir, "boot_irr.csv"))

# Percentile 95% CI (computed on the log scale; equivalent for percentiles)
# One row per model (irr_nb)
summarise_irr <- function(m) {
  irr_ci <- exp(quantile(log(boot_irr[, m]), probs = c(0.025, 0.975), na.rm = TRUE))
  
  tibble(
    model         = m,
    boot_median   = median(boot_irr[, m], na.rm = TRUE),
    lowerCI       = irr_ci[[1]],
    upperCI       = irr_ci[[2]],
    n_boot_failed = sum(is.na(boot_irr[, m]))
  )
}
irr_summary <- bind_rows(lapply(colnames(boot_irr), summarise_irr))
irr_summary

fwrite(irr_summary, file.path(out_dir, "irr_summary.csv"))
