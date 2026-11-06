# Output/figure directory name built from outcome_choice and country_choice
make_dir_name <- function(suffix, nsims = 500) {
  paste0("20261105_",
         if (outcome_choice == "eia") "EIA_infections" else "EIA_diarrhea", "_",
         paste(country_choice, collapse = "_"), "_",
         nsims, "sims", suffix)
}

# Write a warning file to out_dir if any row of `data` has an NA in `covs`;
# remove a stale warning file from an earlier run if there are none.
check_missing_covariates <- function(data, covs, out_dir, data_name) {
  file <- file.path(out_dir, paste0("WARNING_missing_covariates_", data_name, ".txt"))
  miss <- dplyr::filter(data, dplyr::if_any(dplyr::all_of(covs), is.na))
  if (nrow(miss) == 0) {
    if (file.exists(file)) file.remove(file)
    return(invisible(FALSE))
  }
  n_na <- colSums(is.na(data[covs]))
  msg <- c(
    sprintf("WARNING: %d of %d rows (%d of %d Pids) in %s have a missing covariate.",
            nrow(miss), nrow(data), dplyr::n_distinct(miss$Pid),
            dplyr::n_distinct(data$Pid), data_name),
    "",
    "Rows with NA, by covariate:",
    paste0("  ", names(n_na), ": ", n_na),
    "",
    "Affected Pids:",
    paste0("  ", sort(unique(miss$Pid)))
  )
  writeLines(msg, file)
  warning(msg[1], " See ", file, call. = FALSE)
  invisible(TRUE)
}

outcomes  <- c("eia", "eia_diarrhea")
countries <- list(c("BR", "PE")) # list("BR", "PE", c("BR", "PE"))

# 1a. Cox-based clone-censor-weight analyses (first occurrence)
for (outcome_choice in outcomes) {
  for (country_choice in countries) {
    source("../Rota-TTE/1_Cox_dataset.R")
    dir_name <- make_dir_name("_Cox")
    source("../Rota-TTE/2_Cox_IPCW.R")
    source("../Rota-TTE/3_plot_figs_first_outcome.R")
  }
}

# 1a-2. Cox-based clone-censor-weight analysis, severe diarrhea (first occurrence)
outcome_choice <- "severe_diarrhea"
for (country_choice in countries) {
  source("../Rota-TTE/1_Cox_dataset.R")
  source("../Rota-TTE/severe_diarrhea.R")
  dir_name <- paste0("20261105_Severe_diarrhea_",
                     paste(country_choice, collapse = "_"), "_500sims_Cox")
  source("../Rota-TTE/2_Cox_IPCW.R")
  source("../Rota-TTE/3_plot_figs_first_outcome.R")
}

# 1b. Nonparametric bounds analysis
for (outcome_choice in outcomes) {
  for (country_choice in countries) {
    source("../Rota-TTE/1_Cox_dataset.R")
    dir_name <- make_dir_name("_Cox")                 # existing IPCW results (df_result_*.csv)
    for (b in c("no_assumption")) {
      dir_name_bounds <- make_dir_name(paste0("_NPB_", b))
      source(paste0("../Rota-TTE/2b_Cox_NPB_", b, ".R"))
      source("../Rota-TTE/3b_plot_figs_bounds.R")
    }
  }
}

# 2. Pooled logistic regression (first occurrence)
for (outcome_choice in outcomes) {
  for (country_choice in countries) {
    source("../Rota-TTE/4_PersonDayDataset_first_outcome.R")
    dir_name <- make_dir_name("_plr_first_outcome")
    source("../Rota-TTE/5_Plr_first_outcome.R")
    source("../Rota-TTE/3_plot_figs_first_outcome.R")
  }
}

# 3. Cox-based IPCW + negative binomial regression (recurrent outcomes)
for (outcome_choice in outcomes) {
  for (country_choice in countries) {
    source("../Rota-TTE/6_Cox_dataset_recurrent_outcomes.R")
    dir_name <- make_dir_name("_cox_recurrent_outcomes", nsims = 500)
    source("../Rota-TTE/7_Cox_ipcw_recurrent_outcomes.R")
  }
}

# 4. Weight diagnostics
for (outcome_choice in outcomes) {
  for (country_choice in countries) {
    source("../Rota-TTE/1_Cox_dataset.R")
    source("../Rota-TTE/weights_diagnostics.R")
  }
}
