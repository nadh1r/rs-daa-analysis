## ============================================================
## CONTROL DATA — CLEAN + BASELINE + SURVIVAL DATASET (365d)
## ============================================================

library(readxl)
library(dplyr)
library(stringr)
library(lubridate)
library(tidyr)
library(tibble)
library(gtsummary)
library(gt)
library(survival)
library(survminer)

FOLLOWUP_DAYS <- 365

## ------------------------------------------------------------
## Helper functions
## ------------------------------------------------------------

has_col <- function(df, nm) {
  nm %in% names(df)
}

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

standardise_yes_no <- function(x) {
  case_when(
    x %in% c("Yes", "Y", "YES", "yes") ~ "Yes",
    x %in% c("No", "N", "NO", "no") ~ "No",
    TRUE ~ x
  )
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

check_control_rules <- function(df) {
  violations <- list()
  
  if (has_col(df, "DAA Treatment Started")) {
    
    if (has_col(df, "Reason Not Started")) {
      violations$no_reason_when_not_started <- df %>%
        filter(
          `DAA Treatment Started` == "No",
          is.na(`Reason Not Started`) | `Reason Not Started` %in% c("x", "")
        )
      
      violations$reason_present_when_started <- df %>%
        filter(
          `DAA Treatment Started` == "Yes",
          !is.na(`Reason Not Started`),
          `Reason Not Started` != "x"
        )
    }
    
    if (has_col(df, "Date Start DAA Treatment")) {
      violations$date_present_when_not_started <- df %>%
        filter(
          `DAA Treatment Started` == "No",
          !is.na(`Date Start DAA Treatment`),
          `Date Start DAA Treatment` != "x"
        )
      
      violations$date_missing_when_started <- df %>%
        filter(
          `DAA Treatment Started` == "Yes",
          is.na(`Date Start DAA Treatment`) |
            `Date Start DAA Treatment` %in% c("x", "", "999")
        )
    }
    
    if (has_col(df, "DAA Completed")) {
      violations$completed_present_when_not_started <- df %>%
        filter(
          `DAA Treatment Started` == "No",
          !is.na(`DAA Completed`),
          `DAA Completed` != "x"
        )
      
      violations$completed_invalid_when_started <- df %>%
        filter(
          `DAA Treatment Started` == "Yes",
          !is.na(`DAA Completed`),
          !`DAA Completed` %in% c("Yes", "No", "x")
        )
    }
  }
  
  if (has_col(df, "DAA Completed") && has_col(df, "Reason Not Completed")) {
    violations$completed_yes_but_reason_present <- df %>%
      filter(
        `DAA Completed` == "Yes",
        !is.na(`Reason Not Completed`),
        `Reason Not Completed` != "x"
      )
    
    violations$completed_no_but_reason_missing <- df %>%
      filter(
        `DAA Completed` == "No",
        is.na(`Reason Not Completed`) | `Reason Not Completed` %in% c("x", "")
      )
  }
  
  tibble(
    rule = names(violations),
    n_violations = vapply(violations, nrow, integer(1))
  )
}

clean_control_data <- function(path, sheet = 1) {
  raw <- read_excel(path, sheet = sheet)
  
  rules_df <- raw %>%
    filter(!is.na(Name)) %>%
    mutate(across(where(is.character), str_squish))
  
  yn_vars <- intersect(c(
    "DAA Treatment Started", "DAA Completed",
    "HIV co-infection", "Hx of Incarceration",
    "Current Methadone", "Recent drug use", "Hx of injecting drug use",
    "Opioids", "ATS", "Benzos", "Cannabis", "Others"
  ), names(rules_df))
  
  date_vars <- intersect(c(
    "Date RTK Antibody",
    "Date HCV RNA",
    "Date Start DAA Treatment"
  ), names(rules_df))
  
  clean_df <- rules_df %>%
    mutate(across(where(is.character), ~ ifelse(.x %in% c("999", ""), NA_character_, .x))) %>%
    mutate(across(all_of(yn_vars), standardise_yes_no)) %>%
    mutate(across(all_of(date_vars), as_date_any)) %>%
    mutate(
      No = suppressWarnings(as.integer(No)),
      Age = suppressWarnings(as.numeric(Age)),
      Year = suppressWarnings(as.integer(Year))
    )
  
  list(
    raw = raw,
    rules = rules_df,
    clean = clean_df,
    rule_violations = check_control_rules(rules_df)
  )
}

collapse_control_covariates <- function(df) {
  df %>%
    mutate(
      Race_raw = Race,
      Race = str_squish(as.character(Race)),
      Race = case_when(
        Race == "Malay" ~ "Malay",
        is.na(Race) ~ NA_character_,
        TRUE ~ "Non-Malay"
      ),
      Race = factor(Race, levels = c("Malay", "Non-Malay")),
      
      Education_raw = Education,
      Education = str_squish(as.character(Education)),
      Education = case_when(
        Education %in% c(
          "No Formal education", "No formal education",
          "No Formal Education", "Primary school", "Primary"
        ) ~ "Primary or less",
        Education %in% c(
          "Secondary school", "Secondary",
          "Tertiary education", "Tertiary"
        ) ~ "Secondary or more",
        is.na(Education) ~ NA_character_,
        TRUE ~ Education
      ),
      Education = factor(Education, levels = c("Primary or less", "Secondary or more")),
      
      Employment_raw = Employment,
      Employment = str_squish(as.character(Employment)),
      Employment = case_when(
        Employment %in% c("Employed, full time", "Employed, part-time") ~ "Employed",
        Employment == "Unemployed" ~ "Not employed",
        is.na(Employment) ~ NA_character_,
        TRUE ~ Employment
      ),
      Employment = factor(Employment, levels = c("Employed", "Not employed")),
      
      `Marital status_raw` = `Marital status`,
      `Marital status` = str_squish(as.character(`Marital status`)),
      `Marital status` = case_when(
        `Marital status` == "Married" ~ "Married",
        `Marital status` %in% c("Single", "Divorced") ~ "Not married",
        is.na(`Marital status`) ~ NA_character_,
        TRUE ~ `Marital status`
      ),
      `Marital status` = factor(`Marital status`, levels = c("Married", "Not married"))
    )
}

create_baseline_table <- function(df, cohort_label = "Control") {
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
  tti <- df %>%
    filter(
      `RTK Antibody Result` == "Positive",
      `HCV RNA Result` == "Detected"
    ) %>%
    transmute(
      No,
      group = "Control",
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
    mutate(time_raw = as.numeric(daa_date - rtk_date))
  
  invalid_dates <- tti %>%
    filter(!is.na(time_raw), time_raw < 0)
  
  if (nrow(invalid_dates) > 0) {
    stop("Found DAA start date earlier than RTK antibody date.")
  }
  
  tti %>%
    mutate(
      event = if_else(
        !is.na(daa_date) & !is.na(time_raw) & time_raw <= followup_days,
        1L,
        0L
      ),
      time_days = if_else(event == 1L, time_raw, as.numeric(followup_days))
    ) %>%
    select(-time_raw)
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
      censor_status = case_when(
        event == 1L ~ "Started DAA ≤ 365 days",
        event == 0L & is.na(daa_date) ~ "Never started DAA",
        event == 0L &
          !is.na(daa_date) &
          as.numeric(daa_date - rtk_date) > followup_days ~ "Started DAA > 365 days",
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
## Main control-arm pipeline
## ------------------------------------------------------------

control_import <- clean_control_data(
  path = "data/raw/control.xlsx",
  sheet = 1
)

control_clean <- control_import$clean %>%
  collapse_control_covariates()

analysis_vars <- c(
  "Age", "Sex", "Nationality", "Race", "Education", "Employment",
  "Marital status", "Housing", "HIV co-infection",
  "Hx of Incarceration", "Current Methadone",
  "Opioids", "ATS", "Benzos", "Cannabis", "Others",
  "Recent drug use", "Hx of injecting drug use"
)

control_missing <- summarise_missing(control_clean, analysis_vars)

control_baseline_table <- create_baseline_table(
  control_clean,
  cohort_label = "Control"
)

control_cascade <- create_cascade_table(control_clean)

tti_control <- create_tti_dataset(
  control_clean,
  followup_days = FOLLOWUP_DAYS
)

control_tti_summary <- summarise_tti(
  tti_control,
  followup_days = FOLLOWUP_DAYS
)

control_censoring <- create_censoring_table(
  tti_control,
  followup_days = FOLLOWUP_DAYS
)

km_control <- survfit(
  Surv(time_days, event) ~ 1,
  data = tti_control
)

control_crude_median <- create_crude_median_table(tti_control)

control_km_median <- create_km_median_table(
  km_fit = km_control,
  tti_df = tti_control
)

control_km_summary <- create_km_summary_table(
  km_fit = km_control,
  times = c(30, 90, 180, 365)
)

km_plot_control <- ggsurvplot(
  km_control,
  data = tti_control,
  risk.table = TRUE,
  conf.int = TRUE,
  xlab = "Days from RTK antibody date",
  ylab = "Proportion not yet started DAA",
  break.time.by = 30,
  xlim = c(0, FOLLOWUP_DAYS)
)

cox_control <- coxph(
  Surv(time_days, event) ~
    Age +
    Sex +
    Employment +
    Housing +
    `Current Methadone` +
    `Hx of Incarceration` +
    `Recent drug use` +
    `Hx of injecting drug use`,
  data = tti_control
)

cox_control_ph <- cox.zph(cox_control)

control_outputs <- list(
  rule_violations = control_import$rule_violations,
  missingness = control_missing,
  baseline_table = control_baseline_table,
  cascade = control_cascade,
  tti_dataset = tti_control,
  tti_summary = control_tti_summary,
  censoring = control_censoring,
  crude_median = control_crude_median,
  km_fit = km_control,
  km_median = control_km_median,
  km_summary = control_km_summary,
  km_plot = km_plot_control,
  cox_model = cox_control,
  cox_ph = cox_control_ph
)

control_outputs$tti_summary
control_outputs$cascade
control_outputs$rule_violations

