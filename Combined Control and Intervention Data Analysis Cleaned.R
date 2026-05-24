## ============================================================
## COMBINED ANALYSIS — BASELINE + CASCADE + TTI + COX
## GitHub-ready script: no local exports, no DOCX/HTML/PNG output
## Requires:
##   - control84
##   - int_clean
##   - as_date_any()
## ============================================================

library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(purrr)
library(rlang)
library(gtsummary)
library(survival)
library(survminer)
library(splines)

FOLLOWUP_DAYS <- 365

stopifnot(exists("control84"))
stopifnot(exists("int_clean"))
stopifnot(exists("as_date_any"), is.function(as_date_any))

## ------------------------------------------------------------
## Helper functions
## ------------------------------------------------------------

lock_common_levels <- function(df) {
  df %>%
    mutate(
      group = factor(group, levels = c("Control", "Intervention")),
      Age = suppressWarnings(as.numeric(Age)),
      
      Sex = factor(Sex, levels = c("M", "F")),
      Nationality = factor(Nationality, levels = c("Malaysian", "Other")),
      Race = factor(Race, levels = c("Malay", "Non-Malay")),
      Education = factor(Education, levels = c("Primary or less", "Secondary or more")),
      Employment = factor(Employment, levels = c("Employed", "Not employed")),
      `Marital status` = factor(`Marital status`, levels = c("Married", "Not married")),
      Housing = factor(Housing, levels = c("Stable", "Unstable")),
      
      across(
        c(
          `HIV co-infection`, `Hx of Incarceration`, `Current Methadone`,
          Opioids, ATS, Benzos, Cannabis, Others,
          `Recent drug use`, `Hx of injecting drug use`
        ),
        ~ factor(.x, levels = c("No", "Yes"))
      )
    )
}

make_baseline_dataset <- function(control_df, intervention_df) {
  bind_rows(
    control_df %>% mutate(group = "Control"),
    intervention_df %>% mutate(group = "Intervention")
  ) %>%
    select(
      group,
      Age,
      Sex,
      Nationality,
      Race,
      Education,
      Employment,
      `Marital status`,
      Housing,
      `HIV co-infection`,
      `Hx of Incarceration`,
      `Current Methadone`,
      Opioids,
      ATS,
      Benzos,
      Cannabis,
      Others,
      `Recent drug use`,
      `Hx of injecting drug use`
    ) %>%
    lock_common_levels()
}

create_combined_baseline_table <- function(df) {
  df %>%
    tbl_summary(
      by = group,
      statistic = list(
        all_continuous() ~ "{mean} ({sd})",
        all_categorical() ~ "{n} ({p}%)"
      ),
      missing = "no"
    ) %>%
    add_p(
      include = -Nationality,
      test = list(
        all_continuous() ~ "wilcox.test",
        all_categorical() ~ "chisq.test",
        Sex ~ "fisher.test",
        Race ~ "fisher.test",
        Education ~ "fisher.test",
        Opioids ~ "fisher.test"
      ),
      pvalue_fun = ~ style_pvalue(.x, digits = 3)
    ) %>%
    bold_labels()
}

cascade_counts <- function(df, group_label) {
  df %>%
    summarise(
      group = group_label,
      baseline = n(),
      rtk_pos = sum(`RTK Antibody Result` == "Positive", na.rm = TRUE),
      rna_det = sum(
        `RTK Antibody Result` == "Positive" &
          `HCV RNA Result` == "Detected",
        na.rm = TRUE
      ),
      daa_started = sum(
        `RTK Antibody Result` == "Positive" &
          `HCV RNA Result` == "Detected" &
          `DAA Treatment Started` == "Yes",
        na.rm = TRUE
      ),
      daa_completed = sum(
        `RTK Antibody Result` == "Positive" &
          `HCV RNA Result` == "Detected" &
          `DAA Treatment Started` == "Yes" &
          `DAA Completed` == "Yes",
        na.rm = TRUE
      )
    )
}

create_combined_cascade <- function(control_df, intervention_df) {
  bind_rows(
    cascade_counts(control_df, "Control"),
    cascade_counts(intervention_df, "Intervention")
  ) %>%
    mutate(
      group = factor(group, levels = c("Control", "Intervention")),
      rtk_pos_pct_baseline = round(100 * rtk_pos / baseline, 1),
      rna_det_pct_rtk_pos = if_else(rtk_pos > 0, round(100 * rna_det / rtk_pos, 1), NA_real_),
      daa_started_pct_rna_det = if_else(rna_det > 0, round(100 * daa_started / rna_det, 1), NA_real_),
      daa_completed_pct_started = if_else(daa_started > 0, round(100 * daa_completed / daa_started, 1), NA_real_),
      daa_completed_pct_baseline = round(100 * daa_completed / baseline, 1)
    )
}

build_tti_365 <- function(df, id_col, group_label, followup_days = 365) {
  tmp <- df %>%
    filter(
      `RTK Antibody Result` == "Positive",
      `HCV RNA Result` == "Detected"
    ) %>%
    transmute(
      id = .data[[id_col]],
      group = group_label,
      rtk_date = as_date_any(`Date RTK Antibody`),
      daa_date = as_date_any(`Date Start DAA Treatment`),
      
      Age,
      Sex,
      Nationality,
      Race,
      Education,
      Employment,
      `Marital status`,
      Housing,
      `HIV co-infection`,
      `Hx of Incarceration`,
      `Current Methadone`,
      Opioids,
      ATS,
      Benzos,
      Cannabis,
      Others,
      `Recent drug use`,
      `Hx of injecting drug use`
    ) %>%
    filter(!is.na(rtk_date)) %>%
    mutate(
      time_raw = as.numeric(daa_date - rtk_date),
      invalid_negative = !is.na(time_raw) & time_raw < 0
    )
  
  invalid_intervals <- tmp %>%
    filter(invalid_negative) %>%
    select(id, group, rtk_date, daa_date, time_raw)
  
  if (nrow(invalid_intervals) > 0) {
    stop("Negative time intervals detected. Inspect invalid_intervals before analysis.")
  }
  
  tmp %>%
    mutate(
      event = if_else(!is.na(time_raw) & time_raw <= followup_days, 1L, 0L),
      time_days = if_else(event == 1L, time_raw, as.numeric(followup_days)),
      never_flag = is.na(daa_date),
      late_flag = !is.na(time_raw) & time_raw > followup_days
    ) %>%
    select(-time_raw, -invalid_negative) %>%
    lock_common_levels()
}

create_risk_flow <- function(control_df, intervention_df) {
  bind_rows(
    control_df %>% mutate(group = "Control"),
    intervention_df %>% mutate(group = "Intervention")
  ) %>%
    group_by(group) %>%
    summarise(
      baseline = n(),
      rtk_pos = sum(`RTK Antibody Result` == "Positive", na.rm = TRUE),
      rna_detected = sum(
        `RTK Antibody Result` == "Positive" &
          `HCV RNA Result` == "Detected",
        na.rm = TRUE
      ),
      rna_detected_valid_rtk_date = sum(
        `RTK Antibody Result` == "Positive" &
          `HCV RNA Result` == "Detected" &
          !is.na(as_date_any(`Date RTK Antibody`)),
        na.rm = TRUE
      ),
      excluded_missing_rtk_date = rna_detected - rna_detected_valid_rtk_date,
      .groups = "drop"
    ) %>%
    mutate(group = factor(group, levels = c("Control", "Intervention")))
}

create_censoring_breakdown <- function(df) {
  df %>%
    mutate(
      started_le365 = event == 1L,
      never_started = event == 0L & never_flag,
      late_started = event == 0L & late_flag,
      censored_total = event == 0L
    ) %>%
    group_by(group) %>%
    summarise(
      started_le365 = sum(started_le365, na.rm = TRUE),
      never_started = sum(never_started, na.rm = TRUE),
      late_started = sum(late_started, na.rm = TRUE),
      total_censored = sum(censored_total, na.rm = TRUE),
      total_risk_set = n(),
      .groups = "drop"
    )
}

create_crude_median_table <- function(df) {
  df %>%
    filter(event == 1L) %>%
    group_by(group) %>%
    summarise(
      n_initiators = n(),
      median_days = median(time_days, na.rm = TRUE),
      p25_days = quantile(time_days, 0.25, na.rm = TRUE),
      p75_days = quantile(time_days, 0.75, na.rm = TRUE),
      min_days = min(time_days, na.rm = TRUE),
      max_days = max(time_days, na.rm = TRUE),
      .groups = "drop"
    )
}

audit_model_missing <- function(df, vars) {
  tibble(
    variable = vars,
    n_missing = map_int(vars, ~ sum(is.na(df[[.x]]))),
    n_total = nrow(df),
    pct_missing = round(100 * n_missing / n_total, 1)
  ) %>%
    arrange(desc(n_missing), variable)
}

sep_check <- function(df, variable) {
  vq <- sym(variable)
  
  out <- df %>%
    filter(!is.na(!!vq), !is.na(event)) %>%
    count(level = !!vq, event, name = "n") %>%
    pivot_wider(
      names_from = event,
      values_from = n,
      values_fill = 0
    )
  
  if (!"0" %in% names(out)) out[["0"]] <- 0L
  if (!"1" %in% names(out)) out[["1"]] <- 0L
  
  out %>%
    transmute(
      variable = variable,
      level,
      events = as.integer(.data[["1"]]),
      nonevents = as.integer(.data[["0"]]),
      sep_flag = events == 0 | nonevents == 0
    ) %>%
    arrange(desc(sep_flag), level)
}

check_all_separation <- function(df, variables) {
  map_dfr(variables, ~ sep_check(df, .x)) %>%
    group_by(variable) %>%
    mutate(any_sep = any(sep_flag)) %>%
    ungroup()
}

cox_to_tbl <- function(fit, estimate_label = "HR") {
  sm <- summary(fit)
  
  out <- tibble(
    term = rownames(sm$coefficients),
    estimate = sm$coefficients[, "exp(coef)"],
    ci_low = sm$conf.int[, "lower .95"],
    ci_high = sm$conf.int[, "upper .95"],
    p_value = sm$coefficients[, "Pr(>|z|)"]
  )
  
  names(out)[names(out) == "estimate"] <- estimate_label
  out
}

extract_group_effect <- function(fit, term = "groupIntervention") {
  cox_to_tbl(fit) %>%
    filter(.data$term == term)
}

zph_to_tbl <- function(zph_obj) {
  as.data.frame(zph_obj$table) %>%
    rownames_to_column("term") %>%
    as_tibble() %>%
    rename(p_value = p)
}

compare_aic_bic <- function(...) {
  models <- list(...)
  tibble(
    model = names(models),
    AIC = map_dbl(models, AIC),
    BIC = map_dbl(models, BIC)
  )
}

## ------------------------------------------------------------
## 1) Combined baseline table
## ------------------------------------------------------------

baseline_tbl <- make_baseline_dataset(
  control_df = control84,
  intervention_df = int_clean
)

tbl1 <- create_combined_baseline_table(baseline_tbl)

## ------------------------------------------------------------
## 2) Cascade
## ------------------------------------------------------------

cmb_cascade <- create_combined_cascade(
  control_df = control84,
  intervention_df = int_clean
)

## ------------------------------------------------------------
## 3) Survival dataset
## ------------------------------------------------------------

cmb_tti_365 <- bind_rows(
  build_tti_365(control84, "No", "Control", FOLLOWUP_DAYS),
  build_tti_365(int_clean, "record_id", "Intervention", FOLLOWUP_DAYS)
)

tti_definitions <- c(
  "Risk set: RTK antibody positive, HCV RNA detected, and valid RTK date.",
  "Event: DAA treatment started or dispensed within 365 days from RTK date.",
  "Administrative censoring: 365 days.",
  "Censoring subtypes: never started and late started after 365 days.",
  "Negative time intervals are treated as data errors and stop the analysis."
)

risk_flow <- create_risk_flow(control84, int_clean)

censor_breakdown <- create_censoring_breakdown(cmb_tti_365)

crude_median_initiators <- create_crude_median_table(cmb_tti_365)

event_table <- cmb_tti_365 %>%
  count(group, event, name = "n") %>%
  group_by(group) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  ungroup()

## ------------------------------------------------------------
## 4) Missingness and separation checks
## ------------------------------------------------------------

vars_model <- c(
  "Age",
  "Sex",
  "Employment",
  "Housing",
  "HIV co-infection",
  "Current Methadone",
  "Hx of Incarceration",
  "Recent drug use",
  "Hx of injecting drug use"
)

model_missing <- audit_model_missing(cmb_tti_365, vars_model)

vars_cat_check <- c(
  "Sex",
  "Race",
  "Education",
  "Employment",
  "Marital status",
  "Housing",
  "HIV co-infection",
  "Hx of Incarceration",
  "Current Methadone",
  "Opioids",
  "ATS",
  "Benzos",
  "Cannabis",
  "Others",
  "Recent drug use",
  "Hx of injecting drug use"
) %>%
  intersect(names(cmb_tti_365))

sep_cmb <- check_all_separation(cmb_tti_365, vars_cat_check)

drop_vars_sep_cmb <- sep_cmb %>%
  group_by(variable) %>%
  summarise(any_sep = any(sep_flag), .groups = "drop") %>%
  filter(any_sep) %>%
  pull(variable)

## ------------------------------------------------------------
## 5) Kaplan-Meier
## ------------------------------------------------------------

cmb_fit_km <- survfit(
  Surv(time_days, event) ~ group,
  data = cmb_tti_365
)

cmb_km_plot <- ggsurvplot(
  cmb_fit_km,
  data = cmb_tti_365,
  conf.int = TRUE,
  risk.table = TRUE,
  xlab = "Days from RTK antibody date",
  ylab = "Proportion not yet started DAA",
  break.time.by = 30,
  xlim = c(0, FOLLOWUP_DAYS),
  surv.median.line = "hv",
  legend.title = "Study group",
  legend.labs = c("Historical control", "Intervention")
)

logrank_test <- survdiff(
  Surv(time_days, event) ~ group,
  data = cmb_tti_365
)

logrank_p <- 1 - pchisq(
  logrank_test$chisq,
  df = length(logrank_test$n) - 1
)

km_median_tbl <- survminer::surv_median(cmb_fit_km) %>%
  transmute(
    group = str_remove(strata, "^group="),
    median_days = median,
    lower_95_ci = lower,
    upper_95_ci = upper
  )

## ------------------------------------------------------------
## 6) Cox model formulas
## ------------------------------------------------------------

f_group_only <- Surv(time_days, event) ~ group

f_adj <- Surv(time_days, event) ~
  group +
  Age +
  Employment +
  Housing +
  `Marital status` +
  `HIV co-infection` +
  `Current Methadone` +
  `Hx of Incarceration` +
  `Recent drug use` +
  `Hx of injecting drug use`

f_min <- Surv(time_days, event) ~
  group +
  Age +
  Housing +
  `Current Methadone` +
  `Hx of Incarceration` +
  `Hx of injecting drug use`

f_full <- Surv(time_days, event) ~
  group +
  Age +
  Race +
  Education +
  Employment +
  `Marital status` +
  Housing +
  `HIV co-infection` +
  `Hx of Incarceration` +
  `Current Methadone` +
  ATS +
  Benzos +
  Cannabis +
  Others +
  `Recent drug use` +
  `Hx of injecting drug use`

f_adj_strat <- Surv(time_days, event) ~
  group +
  Age +
  Employment +
  Housing +
  `Marital status` +
  `HIV co-infection` +
  `Current Methadone` +
  `Hx of Incarceration` +
  `Hx of injecting drug use` +
  strata(`Recent drug use`)

f_age_spline_strat <- Surv(time_days, event) ~
  group +
  ns(Age, df = 3) +
  Employment +
  Housing +
  `Marital status` +
  `HIV co-infection` +
  `Current Methadone` +
  `Hx of Incarceration` +
  `Hx of injecting drug use` +
  strata(`Recent drug use`)

## ------------------------------------------------------------
## 7) Fit Cox models
## ------------------------------------------------------------

m_group_only <- coxph(
  f_group_only,
  data = cmb_tti_365,
  na.action = na.exclude
)

m_adj <- coxph(
  f_adj,
  data = cmb_tti_365,
  na.action = na.exclude
)

m_min <- coxph(
  f_min,
  data = cmb_tti_365,
  na.action = na.exclude
)

m_full <- coxph(
  f_full,
  data = cmb_tti_365,
  na.action = na.exclude
)

m_adj_strat <- coxph(
  f_adj_strat,
  data = cmb_tti_365,
  na.action = na.exclude
)

m_age_spline_strat <- coxph(
  f_age_spline_strat,
  data = cmb_tti_365,
  na.action = na.exclude
)

tbl_group_only <- cox_to_tbl(m_group_only, "HR")
tbl_adj <- cox_to_tbl(m_adj, "aHR")
tbl_min <- cox_to_tbl(m_min, "HR")
tbl_full <- cox_to_tbl(m_full, "HR")
tbl_adj_strat <- cox_to_tbl(m_adj_strat, "aHR")
tbl_age_spline_strat <- cox_to_tbl(m_age_spline_strat, "aHR")

## ------------------------------------------------------------
## 8) PH diagnostics
## ------------------------------------------------------------

ph_adj <- cox.zph(m_adj)
ph_min <- cox.zph(m_min)
ph_full <- cox.zph(m_full)
ph_strat <- cox.zph(m_adj_strat)

ph_tables <- list(
  adjusted = zph_to_tbl(ph_adj),
  minimal = zph_to_tbl(ph_min),
  full = zph_to_tbl(ph_full),
  stratified_adjusted = zph_to_tbl(ph_strat)
)

model_comparison <- compare_aic_bic(
  minimal = m_min,
  adjusted = m_adj,
  full = m_full,
  stratified_adjusted = m_adj_strat
)

lrt_model_comparison <- anova(
  m_min,
  m_adj,
  m_full,
  test = "LRT"
)

age_spline_lrt <- anova(
  m_adj_strat,
  m_age_spline_strat,
  test = "LRT"
)

## ------------------------------------------------------------
## 9) Influence diagnostics for stratified adjusted model
## ------------------------------------------------------------

DFBETA_THRESHOLD <- 0.5

dfbetas_strat <- residuals(
  m_adj_strat,
  type = "dfbetas"
)

influence_tbl <- cmb_tti_365 %>%
  mutate(
    row_id = row_number(),
    max_abs_dfbetas = apply(abs(dfbetas_strat), 1, max, na.rm = TRUE),
    influential = max_abs_dfbetas > DFBETA_THRESHOLD
  ) %>%
  filter(influential) %>%
  select(
    row_id,
    id,
    group,
    time_days,
    event,
    max_abs_dfbetas,
    influential
  ) %>%
  arrange(desc(max_abs_dfbetas))

influential_rows <- influence_tbl$row_id

m_adj_strat_sens <- if (length(influential_rows) > 0) {
  coxph(
    f_adj_strat,
    data = cmb_tti_365 %>% slice(-influential_rows),
    na.action = na.exclude
  )
} else {
  m_adj_strat
}

sensitivity_group_effect <- bind_rows(
  extract_group_effect(m_adj_strat) %>%
    mutate(model = "Main stratified adjusted"),
  extract_group_effect(m_adj_strat_sens) %>%
    mutate(model = "Sensitivity excluding influential rows")
) %>%
  select(model, everything())

## ------------------------------------------------------------
## 10) Residual summaries
## ------------------------------------------------------------

martingale_residuals <- residuals(
  m_adj_strat,
  type = "martingale"
)

deviance_residuals <- residuals(
  m_adj_strat,
  type = "deviance"
)

residual_summary <- list(
  martingale = summary(martingale_residuals),
  deviance = summary(deviance_residuals)
)

## ------------------------------------------------------------
## 11) Final output object
## ------------------------------------------------------------

combined_outputs <- list(
  definitions = tti_definitions,
  
  baseline_dataset = baseline_tbl,
  baseline_table = tbl1,
  
  cascade = cmb_cascade,
  risk_flow = risk_flow,
  
  tti_dataset = cmb_tti_365,
  event_table = event_table,
  censoring_breakdown = censor_breakdown,
  crude_median_initiators = crude_median_initiators,
  
  missingness = model_missing,
  separation = sep_cmb,
  variables_with_separation = drop_vars_sep_cmb,
  
  km_fit = cmb_fit_km,
  km_plot = cmb_km_plot,
  logrank_test = logrank_test,
  logrank_p = logrank_p,
  km_median = km_median_tbl,
  
  formulas = list(
    group_only = f_group_only,
    adjusted = f_adj,
    minimal = f_min,
    full = f_full,
    stratified_adjusted = f_adj_strat,
    age_spline_stratified = f_age_spline_strat
  ),
  
  models = list(
    group_only = m_group_only,
    adjusted = m_adj,
    minimal = m_min,
    full = m_full,
    stratified_adjusted = m_adj_strat,
    age_spline_stratified = m_age_spline_strat,
    stratified_sensitivity = m_adj_strat_sens
  ),
  
  model_tables = list(
    group_only = tbl_group_only,
    adjusted = tbl_adj,
    minimal = tbl_min,
    full = tbl_full,
    stratified_adjusted = tbl_adj_strat,
    age_spline_stratified = tbl_age_spline_strat
  ),
  
  ph_tests = list(
    adjusted = ph_adj,
    minimal = ph_min,
    full = ph_full,
    stratified_adjusted = ph_strat
  ),
  
  ph_tables = ph_tables,
  model_comparison = model_comparison,
  lrt_model_comparison = lrt_model_comparison,
  age_spline_lrt = age_spline_lrt,
  
  influence = list(
    dfbeta_threshold = DFBETA_THRESHOLD,
    influence_table = influence_tbl,
    sensitivity_group_effect = sensitivity_group_effect
  ),
  
  residuals = list(
    martingale = martingale_residuals,
    deviance = deviance_residuals,
    summary = residual_summary
  )
)

## ------------------------------------------------------------
## Minimal console output
## ------------------------------------------------------------

combined_outputs$event_table
combined_outputs$censoring_breakdown
combined_outputs$km_median
combined_outputs$model_tables$stratified_adjusted
combined_outputs$influence$sensitivity_group_effect