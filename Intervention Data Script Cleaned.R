## ============================================================
## INTERVENTION ARM — CLEAN + BASELINE + SURVIVAL DATASET (365d)
## GitHub-ready standalone script
## ============================================================

library(readr)
library(dplyr)
library(stringr)
library(lubridate)
library(tidyr)
library(tibble)
library(gtsummary)
library(gt)
library(survival)
library(survminer)

## ------------------------------------------------------------
## Configuration
## ------------------------------------------------------------

DATA_PATH <- "data/raw/intervention.csv"

FOLLOWUP_DAYS <- 365
ANALYSIS_YEAR <- 2025

BASE_EVENT <- "lawatan_pertama_mi_arm_1"
HEP_EVENT  <- "status_rawatan_hep_arm_1"

EXCLUDED_IDS <- c(3, 73, 94)

## ------------------------------------------------------------
## Helper functions
## ------------------------------------------------------------

as_date_any <- function(x) {
  if (inherits(x, "Date")) return(x)
  if (inherits(x, c("POSIXct", "POSIXt"))) return(as.Date(x))
  
  x_chr <- str_squish(as.character(x))
  x_chr[x_chr %in% c("x", "X", "", "NA", "999")] <- NA_character_
  
  dt <- suppressWarnings(parse_date_time(
    x_chr,
    orders = c(
      "ymd HMS", "dmy HMS", "mdy HMS",
      "ymd", "dmy", "mdy"
    ),
    tz = "UTC"
  ))
  
  out <- as.Date(dt)
  
  x_num <- suppressWarnings(as.numeric(x_chr))
  idx_serial <- is.na(out) & !is.na(x_num) & x_num > 20000 & x_num < 70000
  out[idx_serial] <- as.Date(x_num[idx_serial], origin = "1899-12-30")
  
  out
}

yes_no_factor <- function(x) {
  factor(x, levels = c("No", "Yes"))
}

summarise_missing <- function(df, vars = names(df)) {
  df %>%
    summarise(across(
      all_of(intersect(vars, names(df))),
      ~ sum(is.na(.)),
      .names = "{.col}"
    )) %>%
    pivot_longer(
      cols = everything(),
      names_to = "variable",
      values_to = "n_missing"
    ) %>%
    mutate(
      n_total = nrow(df),
      pct_missing = round(100 * n_missing / n_total, 1)
    ) %>%
    arrange(desc(n_missing), variable)
}

create_baseline_table <- function(df, cohort_label = "Intervention") {
  baseline_vars <- intersect(c(
    "Sex", "Nationality", "Race", "Education", "Employment",
    "Marital status", "Housing", "HIV co-infection",
    "Hx of Incarceration", "Current Methadone",
    "Opioids", "ATS", "Benzos", "Cannabis", "Others",
    "Recent drug use", "Hx of injecting drug use"
  ), names(df))
  
  df %>%
    select(Age, all_of(baseline_vars)) %>%
    tbl_summary(
      statistic = list(
        all_continuous() ~ "{mean} ({sd})",
        all_categorical() ~ "{n} ({p}%)"
      ),
      missing = "no"
    ) %>%
    modify_spanning_header(
      all_stat_cols() ~ paste0("**", cohort_label, " (N = ", nrow(df), ")**")
    ) %>%
    bold_labels()
}

create_cascade_table <- function(df) {
  n_base <- nrow(df)
  
  tibble(
    stage = c(
      "Baseline extract rows",
      "RTK positive",
      "HCV RNA detected among RTK positive",
      "DAA treatment started among detected",
      "DAA completed among started"
    ),
    count = c(
      n_base,
      sum(df$`RTK Antibody Result` == "Positive", na.rm = TRUE),
      sum(
        df$`RTK Antibody Result` == "Positive" &
          df$`HCV RNA Result` == "Detected",
        na.rm = TRUE
      ),
      sum(
        df$`RTK Antibody Result` == "Positive" &
          df$`HCV RNA Result` == "Detected" &
          df$`DAA Treatment Started` == "Yes",
        na.rm = TRUE
      ),
      sum(
        df$`RTK Antibody Result` == "Positive" &
          df$`HCV RNA Result` == "Detected" &
          df$`DAA Treatment Started` == "Yes" &
          df$`DAA Completed` == "Yes",
        na.rm = TRUE
      )
    )
  ) %>%
    mutate(percentage = round(100 * count / first(count), 1))
}

create_tti_dataset <- function(df, followup_days = 365) {
  tti_raw <- df %>%
    filter(
      `RTK Antibody Result` == "Positive",
      `HCV RNA Result` == "Detected"
    ) %>%
    transmute(
      record_id,
      group = "Intervention",
      rtk_date = `Date RTK Antibody`,
      daa_date = `Date Start DAA Treatment`,
      Age,
      Sex,
      Employment,
      Housing,
      `Current Methadone`,
      `Hx of Incarceration`,
      `Recent drug use`,
      `Hx of injecting drug use`
    ) %>%
    filter(!is.na(rtk_date)) %>%
    mutate(
      time_raw = as.numeric(daa_date - rtk_date),
      interval_status = case_when(
        is.na(daa_date) ~ "no_daa_date",
        is.na(time_raw) ~ "invalid_time_na",
        time_raw < 0 ~ "invalid_negative",
        TRUE ~ "ok"
      )
    )
  
  invalid_intervals <- tti_raw %>%
    filter(interval_status %in% c("invalid_time_na", "invalid_negative"))
  
  if (nrow(invalid_intervals) > 0) {
    stop("Invalid intervention-arm treatment intervals detected. Check intervention_outputs$invalid_intervals.")
  }
  
  tti_raw %>%
    mutate(
      event = if_else(
        !is.na(daa_date) & !is.na(time_raw) & time_raw <= followup_days,
        1L,
        0L
      ),
      time_days = if_else(event == 1L, time_raw, as.numeric(followup_days))
    ) %>%
    select(-time_raw, -interval_status)
}

summarise_tti <- function(tti_df, followup_days = 365) {
  tti_df %>%
    summarise(
      risk_set_n = n(),
      events = sum(event == 1L, na.rm = TRUE),
      censored = sum(event == 0L, na.rm = TRUE),
      never_started = sum(is.na(daa_date)),
      late_started = sum(
        !is.na(daa_date) &
          as.numeric(daa_date - rtk_date) > followup_days,
        na.rm = TRUE
      )
    )
}

create_censoring_table <- function(tti_df, followup_days = 365) {
  tti_df %>%
    mutate(
      time_from_rtk = as.numeric(daa_date - rtk_date),
      censor_status = case_when(
        event == 1L ~ "Started DAA <= 365 days",
        event == 0L & is.na(daa_date) ~ "Never started DAA",
        event == 0L & !is.na(daa_date) & time_from_rtk > followup_days ~
          "Started DAA > 365 days",
        TRUE ~ "Other"
      )
    ) %>%
    count(censor_status, name = "n") %>%
    mutate(
      percentage = round(100 * n / sum(n), 1)
    )
}

create_crude_median_table <- function(tti_df) {
  tti_df %>%
    filter(event == 1L, !is.na(time_days)) %>%
    summarise(
      n_initiators = n(),
      median_days = median(time_days, na.rm = TRUE),
      p25_days = quantile(time_days, 0.25, na.rm = TRUE),
      p75_days = quantile(time_days, 0.75, na.rm = TRUE),
      min_days = min(time_days, na.rm = TRUE),
      max_days = max(time_days, na.rm = TRUE)
    )
}

create_km_median_table <- function(km_fit, tti_df) {
  km_med <- survminer::surv_median(km_fit)
  
  tibble(
    median_days = as.numeric(km_med$median),
    lower_95_ci = as.numeric(km_med$lower),
    upper_95_ci = as.numeric(km_med$upper),
    risk_set_n = nrow(tti_df),
    events = sum(tti_df$event == 1L, na.rm = TRUE),
    censored = sum(tti_df$event == 0L, na.rm = TRUE)
  )
}

create_km_summary_table <- function(km_fit, times = c(30, 90, 180, 365)) {
  summary(km_fit, times = times) |>
    with(
      tibble(
        time_days = time,
        survival = surv,
        lower_95_ci = lower,
        upper_95_ci = upper
      )
    )
}

## ------------------------------------------------------------
## Import REDCap export
## ------------------------------------------------------------

intervention_raw <- read_csv(
  DATA_PATH,
  show_col_types = FALSE
)

## ------------------------------------------------------------
## Define eligible cohort
## ------------------------------------------------------------

eligible_ids <- intervention_raw %>%
  filter(
    redcap_event_name == BASE_EVENT,
    site_allo == 1
  ) %>%
  pull(record_id) %>%
  unique()

intervention_filtered <- intervention_raw %>%
  filter(
    record_id %in% eligible_ids,
    redcap_event_name %in% c(BASE_EVENT, HEP_EVENT),
    !record_id %in% EXCLUDED_IDS
  )

baseline_n <- intervention_filtered %>%
  filter(redcap_event_name == BASE_EVENT) %>%
  distinct(record_id) %>%
  nrow()

## ------------------------------------------------------------
## Split baseline and treatment-status events
## ------------------------------------------------------------

intervention_baseline <- intervention_filtered %>%
  filter(redcap_event_name == BASE_EVENT) %>%
  select(
    record_id,
    
    date_rtk,
    hcv_status,
    rna_date,
    rna_status,
    daa_prescribe,
    daa_not_prescribed_reason___1,
    
    dob,
    gender,
    nat,
    eth,
    edu,
    employment_status,
    mar_status,
    housing,
    housing_oth,
    hiv_pos,
    
    prison_ever,
    past_mth_prison,
    
    current_mmt,
    mmt_rx_ever,
    
    assist_ever_opioids,
    assist_ever_ats,
    assist_ever_cocaine,
    assist_ever_sedatives,
    assist_ever_cannabis,
    assist_ever_inhalants,
    assist_ever_hallucinogens,
    assist_ever_other,
    assist_inject_ever,
    
    assist_freq_cannabis,
    assist_freq_cocaine,
    assist_freq_ats,
    assist_freq_inhalants,
    assist_freq_sedatives,
    assist_freq_hallucinogens,
    assist_freq_opioids,
    assist_freq_other
  ) %>%
  rename(YOB = dob)

intervention_treatment <- intervention_filtered %>%
  filter(redcap_event_name == HEP_EVENT) %>%
  select(
    record_id,
    rna_date_treat = rna_date,
    rna_status_treat = rna_status,
    daa_prescribe_treat = daa_prescribe,
    daa_not_prescribed_reason_treat = daa_not_prescribed_reason___1,
    daa_dispensed,
    daa_treatment_status,
    daa_not_completed_reason___1
  )

intervention_merged <- intervention_baseline %>%
  left_join(intervention_treatment, by = "record_id") %>%
  mutate(
    rna_date_final = coalesce(rna_date_treat, rna_date),
    rna_status_final = coalesce(rna_status_treat, rna_status),
    daa_prescribe_final = coalesce(daa_prescribe_treat, daa_prescribe),
    daa_not_prescribed_reason_final = coalesce(
      daa_not_prescribed_reason_treat,
      daa_not_prescribed_reason___1
    )
  )

stopifnot(nrow(intervention_merged) == baseline_n)

## ------------------------------------------------------------
## Harmonise to control-arm variable names
## ------------------------------------------------------------

intervention_clean <- intervention_merged %>%
  mutate(
    `Date RTK Antibody` = as_date_any(date_rtk),
    `Date HCV RNA` = as_date_any(rna_date_final),
    `Date Start DAA Treatment` = as_date_any(daa_dispensed),
    
    `RTK Antibody Result` = case_when(
      as.character(hcv_status) == "1" ~ "Positive",
      as.character(hcv_status) == "0" ~ "Negative",
      TRUE ~ NA_character_
    ),
    
    `HCV RNA Result` = case_when(
      as.character(rna_status_final) == "1" ~ "Detected",
      as.character(rna_status_final) == "0" ~ "Not Detected",
      TRUE ~ NA_character_
    ),
    
    `DAA Treatment Started` = case_when(
      as.character(daa_prescribe_final) == "1" ~ "Yes",
      as.character(daa_prescribe_final) == "0" ~ "No",
      TRUE ~ NA_character_
    ),
    
    `DAA Completed` = case_when(
      as.character(daa_treatment_status) == "1" ~ "Yes",
      as.character(daa_treatment_status) == "0" ~ "No",
      TRUE ~ NA_character_
    ),
    
    `Reason Not Started` = as.character(daa_not_prescribed_reason_final),
    `Reason Not Completed` = as.character(daa_not_completed_reason___1),
    
    YOB = suppressWarnings(as.integer(YOB)),
    Age = if_else(!is.na(YOB), as.numeric(ANALYSIS_YEAR - YOB), NA_real_)
  )

## ------------------------------------------------------------
## Demographics and baseline covariates
## ------------------------------------------------------------

intervention_clean <- intervention_clean %>%
  mutate(
    Sex = case_when(
      as.character(gender) == "0" ~ "M",
      as.character(gender) == "1" ~ "F",
      TRUE ~ NA_character_
    ),
    
    Nationality = case_when(
      as.character(nat) == "0" ~ "Malaysian",
      as.character(nat) == "1" ~ "Other",
      TRUE ~ NA_character_
    ),
    
    Race_raw = case_when(
      eth == 1 ~ "Malay",
      eth == 2 ~ "Chinese",
      eth == 3 ~ "Indian",
      eth == 4 ~ "Others",
      TRUE ~ NA_character_
    ),
    
    Race = case_when(
      Race_raw == "Malay" ~ "Malay",
      Race_raw %in% c("Chinese", "Indian", "Others") ~ "Non-Malay",
      TRUE ~ NA_character_
    ),
    
    Education_raw = case_when(
      edu == 0 ~ "No Formal education",
      edu %in% c(1, 2) ~ "Primary school",
      edu %in% c(3, 4) ~ "Secondary school",
      edu %in% c(5, 6) ~ "Tertiary education",
      TRUE ~ NA_character_
    ),
    
    Education = case_when(
      Education_raw %in% c("No Formal education", "Primary school") ~ "Primary or less",
      Education_raw %in% c("Secondary school", "Tertiary education") ~ "Secondary or more",
      TRUE ~ NA_character_
    ),
    
    Employment_raw = case_when(
      employment_status == 0 ~ "Unemployed",
      employment_status == 1 ~ "Employed, part-time",
      employment_status == 2 ~ "Employed, full time",
      TRUE ~ NA_character_
    ),
    
    Employment = case_when(
      Employment_raw %in% c("Employed, full time", "Employed, part-time") ~ "Employed",
      Employment_raw == "Unemployed" ~ "Not employed",
      TRUE ~ NA_character_
    ),
    
    `Marital status_raw` = case_when(
      mar_status %in% c(1, 3) ~ "Married",
      mar_status %in% c(2, 4) ~ "Divorced",
      mar_status %in% c(5, 6) ~ "Single",
      TRUE ~ NA_character_
    ),
    
    `Marital status` = case_when(
      `Marital status_raw` == "Married" ~ "Married",
      `Marital status_raw` %in% c("Single", "Divorced") ~ "Not married",
      TRUE ~ NA_character_
    ),
    
    housing_other_lower = str_to_lower(str_squish(coalesce(housing_oth, ""))),
    
    Housing = case_when(
      housing %in% 1:9 ~ "Stable",
      housing == 10 ~ "Unstable",
      housing == 11 & str_detect(housing_other_lower, "tempat kerja") ~ "Unstable",
      housing == 11 ~ "Stable",
      TRUE ~ NA_character_
    ),
    
    `HIV co-infection` = case_when(
      as.character(hiv_pos) == "1" ~ "Yes",
      as.character(hiv_pos) == "0" ~ "No",
      TRUE ~ NA_character_
    )
  ) %>%
  select(-housing_other_lower)

## ------------------------------------------------------------
## Incarceration, methadone, and drug-use variables
## ------------------------------------------------------------

intervention_clean <- intervention_clean %>%
  mutate(
    incarceration_status = case_when(
      prison_ever == 1 & past_mth_prison == 1 ~ "Recent incarceration",
      prison_ever == 1 & (past_mth_prison == 0 | is.na(past_mth_prison)) ~ "Past incarceration",
      prison_ever == 0 ~ "Never incarcerated",
      TRUE ~ NA_character_
    ),
    
    incarceration_status = factor(
      incarceration_status,
      levels = c(
        "Never incarcerated",
        "Past incarceration",
        "Recent incarceration"
      )
    ),
    
    `Hx of Incarceration` = case_when(
      prison_ever == 1 ~ "Yes",
      prison_ever == 0 ~ "No",
      TRUE ~ NA_character_
    ),
    
    `Current Methadone` = case_when(
      mmt_rx_ever == 1 & current_mmt == 1 ~ "Yes",
      mmt_rx_ever == 1 & current_mmt == 0 ~ "No",
      mmt_rx_ever == 0 ~ "No",
      TRUE ~ NA_character_
    ),
    
    Opioids = case_when(
      as.character(assist_ever_opioids) == "1" ~ "Yes",
      as.character(assist_ever_opioids) == "0" ~ "No",
      TRUE ~ NA_character_
    ),
    
    ATS = case_when(
      as.character(assist_ever_ats) == "1" ~ "Yes",
      as.character(assist_ever_ats) == "0" ~ "No",
      TRUE ~ NA_character_
    ),
    
    Benzos = case_when(
      as.character(assist_ever_sedatives) == "1" ~ "Yes",
      as.character(assist_ever_sedatives) == "0" ~ "No",
      TRUE ~ NA_character_
    ),
    
    Cannabis = case_when(
      as.character(assist_ever_cannabis) == "1" ~ "Yes",
      as.character(assist_ever_cannabis) == "0" ~ "No",
      TRUE ~ NA_character_
    ),
    
    Others = case_when(
      assist_ever_cocaine == 1 |
        assist_ever_inhalants == 1 |
        assist_ever_hallucinogens == 1 |
        assist_ever_other == 1 ~ "Yes",
      assist_ever_cocaine == 0 &
        assist_ever_inhalants == 0 &
        assist_ever_hallucinogens == 0 &
        assist_ever_other == 0 ~ "No",
      TRUE ~ NA_character_
    ),
    
    `Hx of injecting drug use` = case_when(
      as.character(assist_inject_ever) == "1" ~ "Yes",
      as.character(assist_inject_ever) == "0" ~ "No",
      TRUE ~ NA_character_
    ),
    
    max_freq = pmax(
      assist_freq_cannabis,
      assist_freq_cocaine,
      assist_freq_ats,
      assist_freq_inhalants,
      assist_freq_sedatives,
      assist_freq_hallucinogens,
      assist_freq_opioids,
      assist_freq_other,
      na.rm = TRUE
    ),
    
    max_freq = if_else(is.infinite(max_freq), NA_real_, as.numeric(max_freq)),
    
    `Recent drug use` = case_when(
      max_freq %in% c(2, 3, 4, 6) ~ "Yes",
      max_freq == 0 ~ "No",
      TRUE ~ NA_character_
    )
  )

## ------------------------------------------------------------
## Final analytic intervention dataset
## ------------------------------------------------------------

intervention_clean <- intervention_clean %>%
  transmute(
    record_id,
    
    Age,
    Sex = factor(Sex, levels = c("M", "F")),
    Nationality = factor(Nationality, levels = c("Malaysian", "Other")),
    Race = factor(Race, levels = c("Malay", "Non-Malay")),
    Education = factor(Education, levels = c("Primary or less", "Secondary or more")),
    Employment = factor(Employment, levels = c("Employed", "Not employed")),
    `Marital status` = factor(`Marital status`, levels = c("Married", "Not married")),
    Housing = factor(Housing, levels = c("Stable", "Unstable")),
    `HIV co-infection` = yes_no_factor(`HIV co-infection`),
    
    incarceration_status,
    `Hx of Incarceration` = yes_no_factor(`Hx of Incarceration`),
    `Current Methadone` = yes_no_factor(`Current Methadone`),
    
    Opioids = yes_no_factor(Opioids),
    ATS = yes_no_factor(ATS),
    Benzos = yes_no_factor(Benzos),
    Cannabis = yes_no_factor(Cannabis),
    Others = yes_no_factor(Others),
    `Recent drug use` = yes_no_factor(`Recent drug use`),
    `Hx of injecting drug use` = yes_no_factor(`Hx of injecting drug use`),
    
    Race_raw,
    Education_raw,
    Employment_raw,
    `Marital status_raw`,
    
    `Date RTK Antibody`,
    `RTK Antibody Result`,
    `Date HCV RNA`,
    `HCV RNA Result`,
    `Date Start DAA Treatment`,
    `DAA Treatment Started`,
    `DAA Completed`,
    `Reason Not Started`,
    `Reason Not Completed`
  )

## ------------------------------------------------------------
## Quality checks
## ------------------------------------------------------------

analysis_vars <- c(
  "Age", "Sex", "Nationality", "Race", "Education", "Employment",
  "Marital status", "Housing", "HIV co-infection",
  "Hx of Incarceration", "Current Methadone",
  "Opioids", "ATS", "Benzos", "Cannabis", "Others",
  "Recent drug use", "Hx of injecting drug use"
)

intervention_missing <- summarise_missing(
  intervention_clean,
  analysis_vars
)

intervention_start_prescription_check <- intervention_clean %>%
  summarise(
    started_without_prescription = sum(
      `DAA Treatment Started` == "No" &
        !is.na(`Date Start DAA Treatment`),
      na.rm = TRUE
    )
  )

## ------------------------------------------------------------
## Baseline table
## ------------------------------------------------------------

intervention_baseline_table <- create_baseline_table(
  intervention_clean,
  cohort_label = "Intervention"
)

## ------------------------------------------------------------
## Cascade table
## ------------------------------------------------------------

intervention_cascade <- create_cascade_table(
  intervention_clean
)

## ------------------------------------------------------------
## Survival dataset
## ------------------------------------------------------------

tti_intervention <- create_tti_dataset(
  intervention_clean,
  followup_days = FOLLOWUP_DAYS
)

intervention_tti_summary <- summarise_tti(
  tti_intervention,
  followup_days = FOLLOWUP_DAYS
)

intervention_censoring <- create_censoring_table(
  tti_intervention,
  followup_days = FOLLOWUP_DAYS
)

## ------------------------------------------------------------
## Kaplan-Meier analysis
## ------------------------------------------------------------

km_intervention <- survfit(
  Surv(time_days, event) ~ 1,
  data = tti_intervention
)

intervention_crude_median <- create_crude_median_table(
  tti_intervention
)

intervention_km_median <- create_km_median_table(
  km_fit = km_intervention,
  tti_df = tti_intervention
)

intervention_km_summary <- create_km_summary_table(
  km_fit = km_intervention,
  times = c(30, 90, 180, 365)
)

km_plot_intervention <- ggsurvplot(
  km_intervention,
  data = tti_intervention,
  risk.table = TRUE,
  conf.int = TRUE,
  xlab = "Days from RTK antibody date",
  ylab = "Proportion not yet started DAA",
  break.time.by = 30,
  xlim = c(0, FOLLOWUP_DAYS)
)

## ------------------------------------------------------------
## Cox model
## ------------------------------------------------------------

cox_intervention <- coxph(
  Surv(time_days, event) ~
    Age +
    Sex +
    Employment +
    Housing +
    `Current Methadone` +
    `Hx of Incarceration` +
    `Recent drug use` +
    `Hx of injecting drug use`,
  data = tti_intervention
)

cox_intervention_ph <- cox.zph(cox_intervention)

cox_required_vars <- c(
  "time_days", "event", "Age", "Sex", "Employment",
  "Housing", "Current Methadone", "Hx of Incarceration",
  "Recent drug use", "Hx of injecting drug use"
)

cox_complete_case_summary <- tibble(
  risk_set_n = nrow(tti_intervention),
  cox_complete_cases = sum(complete.cases(tti_intervention[, cox_required_vars])),
  rows_dropped_for_missing_covariates =
    risk_set_n - cox_complete_cases
)

cox_missing_by_variable <- tti_intervention %>%
  summarise(across(
    all_of(cox_required_vars),
    ~ sum(is.na(.)),
    .names = "{.col}"
  )) %>%
  pivot_longer(
    cols = everything(),
    names_to = "variable",
    values_to = "n_missing"
  ) %>%
  arrange(desc(n_missing), variable)

## ------------------------------------------------------------
## Final output object
## ------------------------------------------------------------

intervention_outputs <- list(
  clean_data = intervention_clean,
  missingness = intervention_missing,
  start_prescription_check = intervention_start_prescription_check,
  baseline_table = intervention_baseline_table,
  cascade = intervention_cascade,
  tti_dataset = tti_intervention,
  tti_summary = intervention_tti_summary,
  censoring = intervention_censoring,
  crude_median = intervention_crude_median,
  km_fit = km_intervention,
  km_median = intervention_km_median,
  km_summary = intervention_km_summary,
  km_plot = km_plot_intervention,
  cox_model = cox_intervention,
  cox_ph = cox_intervention_ph,
  cox_complete_case_summary = cox_complete_case_summary,
  cox_missing_by_variable = cox_missing_by_variable
)

## Minimal console output
intervention_outputs$tti_summary
intervention_outputs$cascade
intervention_outputs$start_prescription_check