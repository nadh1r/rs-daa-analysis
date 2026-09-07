# ==============================================================================
# RS-DAA THESIS: REPRODUCIBLE STATISTICAL ANALYSIS
# Primary mITT cohort excludes cirrhosis requiring specialist referral
# R version used for the thesis: 4.5.0
# ==============================================================================

# IMPORTANT COHORT RULE
# The HCV-antibody-positive population is used for baseline and cascade summaries.
# The primary mITT population is constructed once below and excludes cirrhosis.
# Every time-to-event, 30-day, RMST and IPTW analysis uses that locked mITT object.

# ---- 0. User settings ---------------------------------------------------------
CONTROL_FILE <- "path/to/control_masterlist.xlsx"
INTERVENTION_FILE <- "path/to/redcap_intervention_export.csv"
OUTPUT_DIR <- "RSDAA_thesis_outputs"

FOLLOWUP_DAYS <- 365L
RAPID_DAYS <- c(7L, 14L, 30L)
BOOTSTRAP_REPS <- 2000L
BOOTSTRAP_SEED <- 20260907L

required_packages <- c(
  "readxl", "readr", "dplyr", "tidyr", "stringr", "lubridate",
  "survival", "survRM2", "ggplot2", "scales", "broom"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Install the following packages before running the script: ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(lubridate)
  library(survival)
  library(ggplot2)
})

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
stopifnot(file.exists(CONTROL_FILE), file.exists(INTERVENTION_FILE))

# ---- 1. Reusable functions ----------------------------------------------------
as_date_any <- function(x) {
  if (inherits(x, "Date")) return(x)
  if (inherits(x, c("POSIXct", "POSIXt"))) return(as.Date(x))
  z <- str_squish(as.character(x))
  z[z %in% c("x", "X", "", "NA", "999")] <- NA_character_
  parsed <- suppressWarnings(parse_date_time(
    z,
    orders = c("ymd HMS", "dmy HMS", "mdy HMS", "ymd", "dmy", "mdy"),
    tz = "UTC"
  ))
  out <- as.Date(parsed)
  serial <- suppressWarnings(as.numeric(z))
  use_serial <- is.na(out) & !is.na(serial) & serial > 20000 & serial < 70000
  out[use_serial] <- as.Date(serial[use_serial], origin = "1899-12-30")
  out
}

yes_no <- function(x) {
  z <- str_to_lower(str_squish(as.character(x)))
  case_when(
    z %in% c("yes", "y", "1") ~ "Yes",
    z %in% c("no", "n", "0") ~ "No",
    TRUE ~ NA_character_
  )
}

assert_columns <- function(data, columns, label) {
  absent <- setdiff(columns, names(data))
  if (length(absent)) {
    stop(label, " is missing required columns: ", paste(absent, collapse = ", "))
  }
  invisible(TRUE)
}

format_p <- function(p) ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))

extract_cox_group <- function(model, label) {
  s <- summary(model)
  term <- "groupIntervention"
  stopifnot(term %in% rownames(s$coefficients))
  tibble(
    model = label,
    HR = unname(s$coefficients[term, "exp(coef)"]),
    lower_95_CI = unname(s$conf.int[term, "lower .95"]),
    upper_95_CI = unname(s$conf.int[term, "upper .95"]),
    p_value = unname(s$coefficients[term, "Pr(>|z|)"])
  )
}

weighted_smd_continuous <- function(x, treatment, weights = rep(1, length(x))) {
  calc <- function(g) {
    keep <- treatment == g
    m <- weighted.mean(x[keep], weights[keep])
    v <- weighted.mean((x[keep] - m)^2, weights[keep])
    c(mean = m, variance = v)
  }
  a <- calc(1L); b <- calc(0L)
  unname((a["mean"] - b["mean"]) / sqrt((a["variance"] + b["variance"]) / 2))
}

weighted_smd_binary <- function(x, treatment, level,
                                weights = rep(1, length(x))) {
  z <- as.numeric(x == level)
  p1 <- weighted.mean(z[treatment == 1L], weights[treatment == 1L])
  p0 <- weighted.mean(z[treatment == 0L], weights[treatment == 0L])
  pooled <- (p1 + p0) / 2
  if (pooled %in% c(0, 1)) return(NA_real_)
  (p1 - p0) / sqrt(pooled * (1 - pooled))
}

# Explicit Nelson-Aalen calculation. survfit supplies event times, risk sets and
# event counts; H(t) is calculated as cumulative sum d(t)/n(t).
nelson_aalen_by_group <- function(data) {
  levels(data$group) |>
    lapply(function(g) {
      one <- droplevels(filter(data, group == g))
      fit <- survfit(Surv(time_days, event) ~ 1, data = one, conf.type = "log")
      keep <- fit$n.event > 0
      h <- cumsum(fit$n.event[keep] / fit$n.risk[keep])
      var_h <- cumsum(fit$n.event[keep] / fit$n.risk[keep]^2)
      se_h <- sqrt(var_h)
      lower_h <- pmax(0, h - qnorm(0.975) * se_h)
      upper_h <- h + qnorm(0.975) * se_h
      tibble(
        group = g,
        time = fit$time[keep],
        n_risk = fit$n.risk[keep],
        n_event = fit$n.event[keep],
        cumulative_events = cumsum(fit$n.event[keep]),
        cumulative_hazard = h,
        cumulative_incidence = 1 - exp(-h),
        lower_95_CI = 1 - exp(-lower_h),
        upper_95_CI = 1 - exp(-upper_h)
      )
    }) |>
    bind_rows()
}

na_median <- function(data) {
  est <- nelson_aalen_by_group(data)
  est |>
    group_by(group) |>
    summarise(
      median_days = if (any(cumulative_incidence >= 0.5))
        min(time[cumulative_incidence >= 0.5]) else NA_real_,
      .groups = "drop"
    )
}

# ---- 2. Import and clean historical-control data ------------------------------
control_raw <- readxl::read_excel(CONTROL_FILE, sheet = 1)
control_raw <- control_raw |>
  filter(!is.na(Name)) |>
  mutate(across(where(is.character), str_squish))

control_required <- c(
  "No", "Name", "Age", "Sex", "Nationality", "Race", "Education",
  "Employment", "Marital status", "Housing", "HIV co-infection",
  "Hx of Incarceration", "Current Methadone", "Opioids", "ATS", "Benzos",
  "Cannabis", "Others", "Recent drug use", "Hx of injecting drug use",
  "Date RTK Antibody", "RTK Antibody Result", "Date HCV RNA",
  "HCV RNA Result", "Date Start DAA Treatment", "DAA Treatment Started",
  "DAA Completed", "Reason Not Started"
)
assert_columns(control_raw, control_required, "Control file")

control <- control_raw |>
  mutate(
    across(where(is.character), ~ na_if(na_if(.x, "999"), "")),
    across(
      all_of(c(
        "HIV co-infection", "Hx of Incarceration", "Current Methadone",
        "Opioids", "ATS", "Benzos", "Cannabis", "Others",
        "Recent drug use", "Hx of injecting drug use",
        "DAA Treatment Started", "DAA Completed"
      )),
      yes_no
    ),
    across(
      all_of(c("Date RTK Antibody", "Date HCV RNA", "Date Start DAA Treatment")),
      as_date_any
    ),
    Age = suppressWarnings(as.numeric(Age)),
    Race = if_else(Race == "Malay", "Malay", "Non-Malay", missing = NA_character_),
    Education = case_when(
      Education %in% c("No Formal education", "No formal education",
                       "No Formal Education", "Primary school", "Primary") ~
        "Primary or less",
      Education %in% c("Secondary school", "Secondary", "Tertiary education",
                       "Tertiary") ~ "Secondary or more",
      TRUE ~ NA_character_
    ),
    Employment = case_when(
      Employment %in% c("Employed, full time", "Employed, part-time") ~ "Employed",
      Employment == "Unemployed" ~ "Not employed",
      TRUE ~ NA_character_
    ),
    `Marital status` = case_when(
      `Marital status` == "Married" ~ "Married",
      `Marital status` %in% c("Single", "Divorced") ~ "Not married",
      TRUE ~ NA_character_
    ),
    Cirrhosis = if_else(
      `HCV RNA Result` == "Detected" &
        str_detect(coalesce(as.character(`Reason Not Started`), ""),
                   regex("cirrhosis", ignore_case = TRUE)),
      "Yes", "No"
    )
  )

# ---- 3. Import and clean prospective intervention data ------------------------
intervention_raw <- readr::read_csv(INTERVENTION_FILE, show_col_types = FALSE)
intervention_required <- c(
  "record_id", "redcap_event_name", "site_allo", "date_rtk", "hcv_status",
  "rna_date", "rna_status", "daa_prescribe", "daa_not_prescribed_reason___1",
  "daa_not_prescribed_reason___5", "dob", "gender", "nat", "eth", "edu",
  "employment_status", "mar_status", "housing", "housing_oth", "hiv_pos",
  "prison_ever", "past_mth_prison", "current_mmt", "mmt_rx_ever",
  "assist_ever_opioids", "assist_ever_ats", "assist_ever_cocaine",
  "assist_ever_sedatives", "assist_ever_cannabis", "assist_ever_inhalants",
  "assist_ever_hallucinogens", "assist_ever_other", "assist_inject_ever",
  "assist_freq_cannabis", "assist_freq_cocaine", "assist_freq_ats",
  "assist_freq_inhalants", "assist_freq_sedatives", "assist_freq_hallucinogens",
  "assist_freq_opioids", "assist_freq_other", "daa_dispensed",
  "daa_treatment_status", "daa_not_completed_reason___1"
)
assert_columns(intervention_raw, intervention_required, "Intervention file")

BASE_EVENT <- "lawatan_pertama_mi_arm_1"
TREATMENT_EVENT <- "status_rawatan_hep_arm_1"

eligible_ids <- intervention_raw |>
  filter(redcap_event_name == BASE_EVENT, site_allo == 1) |>
  distinct(record_id) |>
  pull(record_id)

intervention_filtered <- intervention_raw |>
  filter(
    record_id %in% eligible_ids,
    redcap_event_name %in% c(BASE_EVENT, TREATMENT_EVENT),
    !record_id %in% c(3, 73, 94)
  )

intervention_base <- intervention_filtered |>
  filter(redcap_event_name == BASE_EVENT) |>
  select(
    record_id, date_rtk, hcv_status, rna_date, rna_status, daa_prescribe,
    daa_not_prescribed_reason___1, daa_not_prescribed_reason___5,
    dob, gender, nat, eth, edu, employment_status, mar_status, housing,
    housing_oth, hiv_pos, prison_ever, past_mth_prison, current_mmt, mmt_rx_ever,
    starts_with("assist_ever_"), assist_inject_ever, starts_with("assist_freq_")
  ) |>
  rename(YOB = dob)

intervention_treatment <- intervention_filtered |>
  filter(redcap_event_name == TREATMENT_EVENT) |>
  select(
    record_id,
    rna_date_2 = rna_date, rna_status_2 = rna_status,
    daa_prescribe_2 = daa_prescribe,
    daa_reason_1_2 = daa_not_prescribed_reason___1,
    daa_reason_5_2 = daa_not_prescribed_reason___5,
    daa_dispensed, daa_treatment_status, daa_not_completed_reason___1
  )

intervention_joined <- intervention_base |>
  left_join(intervention_treatment, by = "record_id") |>
  mutate(
    rna_date_final = coalesce(rna_date_2, rna_date),
    rna_status_final = coalesce(rna_status_2, rna_status),
    daa_prescribe_final = coalesce(daa_prescribe_2, daa_prescribe),
    daa_reason_1_final = coalesce(daa_reason_1_2, daa_not_prescribed_reason___1),
    cirrhosis_checkbox_final = as.integer(
      coalesce(as.character(daa_not_prescribed_reason___5) == "1", FALSE) |
        coalesce(as.character(daa_reason_5_2) == "1", FALSE)
    )
  )

stopifnot(!anyDuplicated(intervention_joined$record_id))

intervention <- intervention_joined |>
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
    `Reason Not Started` = as.character(daa_reason_1_final),
    Cirrhosis = if_else(
      as.character(rna_status_final) == "1" & cirrhosis_checkbox_final == 1L,
      "Yes", "No"
    ),
    Age = if_else(!is.na(suppressWarnings(as.integer(YOB))),
                  2025 - suppressWarnings(as.integer(YOB)), NA_integer_),
    Sex = case_when(as.character(gender) == "0" ~ "M",
                    as.character(gender) == "1" ~ "F"),
    Nationality = case_when(as.character(nat) == "0" ~ "Malaysian",
                            as.character(nat) == "1" ~ "Other"),
    Race = if_else(eth == 1, "Malay", "Non-Malay", missing = NA_character_),
    Education = case_when(
      edu %in% 0:2 ~ "Primary or less",
      edu %in% 3:6 ~ "Secondary or more"
    ),
    Employment = case_when(
      employment_status %in% 1:2 ~ "Employed",
      employment_status == 0 ~ "Not employed"
    ),
    `Marital status` = case_when(
      mar_status %in% c(1, 3) ~ "Married",
      mar_status %in% c(2, 4, 5, 6) ~ "Not married"
    ),
    housing_text = str_to_lower(str_squish(coalesce(housing_oth, ""))),
    Housing = case_when(
      housing %in% 1:9 ~ "Stable",
      housing == 10 ~ "Unstable",
      housing == 11 & str_detect(housing_text, "tempat kerja") ~ "Unstable",
      housing == 11 ~ "Stable"
    ),
    `HIV co-infection` = yes_no(hiv_pos),
    `Hx of Incarceration` = yes_no(prison_ever),
    `Current Methadone` = case_when(
      mmt_rx_ever == 1 & current_mmt == 1 ~ "Yes",
      mmt_rx_ever %in% c(0, 1) ~ "No"
    ),
    Opioids = yes_no(assist_ever_opioids),
    ATS = yes_no(assist_ever_ats),
    Benzos = yes_no(assist_ever_sedatives),
    Cannabis = yes_no(assist_ever_cannabis),
    Others = case_when(
      assist_ever_cocaine == 1 | assist_ever_inhalants == 1 |
        assist_ever_hallucinogens == 1 | assist_ever_other == 1 ~ "Yes",
      assist_ever_cocaine == 0 & assist_ever_inhalants == 0 &
        assist_ever_hallucinogens == 0 & assist_ever_other == 0 ~ "No"
    ),
    `Hx of injecting drug use` = yes_no(assist_inject_ever),
    max_recent_frequency = pmax(
      assist_freq_cannabis, assist_freq_cocaine, assist_freq_ats,
      assist_freq_inhalants, assist_freq_sedatives, assist_freq_hallucinogens,
      assist_freq_opioids, assist_freq_other, na.rm = TRUE
    ),
    max_recent_frequency = if_else(
      is.infinite(max_recent_frequency), NA_real_, max_recent_frequency
    ),
    `Recent drug use` = case_when(
      max_recent_frequency %in% c(2, 3, 4, 6) ~ "Yes",
      max_recent_frequency == 0 ~ "No"
    )
  ) |>
  select(
    record_id, Age, Sex, Nationality, Race, Education, Employment,
    `Marital status`, Housing, `HIV co-infection`, `Hx of Incarceration`,
    `Current Methadone`, Opioids, ATS, Benzos, Cannabis, Others,
    `Recent drug use`, `Hx of injecting drug use`, Cirrhosis,
    `Date RTK Antibody`, `RTK Antibody Result`, `Date HCV RNA`, `HCV RNA Result`,
    `Date Start DAA Treatment`, `DAA Treatment Started`, `DAA Completed`,
    `Reason Not Started`
  )

# ---- 4. Harmonise factor levels -----------------------------------------------
factorise <- function(data) {
  data |>
    mutate(
      Sex = factor(Sex, c("M", "F")),
      Nationality = factor(Nationality, c("Malaysian", "Other")),
      Race = factor(Race, c("Malay", "Non-Malay")),
      Education = factor(Education, c("Primary or less", "Secondary or more")),
      Employment = factor(Employment, c("Employed", "Not employed")),
      `Marital status` = factor(`Marital status`, c("Married", "Not married")),
      Housing = factor(Housing, c("Stable", "Unstable")),
      across(
        all_of(c(
          "HIV co-infection", "Hx of Incarceration", "Current Methadone",
          "Opioids", "ATS", "Benzos", "Cannabis", "Others",
          "Recent drug use", "Hx of injecting drug use", "Cirrhosis"
        )),
        ~ factor(.x, c("No", "Yes"))
      )
    )
}
control <- factorise(control)
intervention <- factorise(intervention)

# ---- 5. Data-quality and cirrhosis audit --------------------------------------
stopifnot(!anyDuplicated(control$No), !anyDuplicated(intervention$record_id))
stopifnot(!any(is.na(control$Cirrhosis)), !any(is.na(intervention$Cirrhosis)))
stopifnot(!any(intervention$Cirrhosis == "Yes" &
                intervention$`HCV RNA Result` != "Detected"))

cirrhosis_audit <- bind_rows(
  control |>
    filter(Cirrhosis == "Yes") |>
    transmute(id = No, group = "Control", `HCV RNA Result`,
              `Reason Not Started`, `Date Start DAA Treatment`),
  intervention |>
    filter(Cirrhosis == "Yes") |>
    transmute(id = record_id, group = "Intervention", `HCV RNA Result`,
              `Reason Not Started`, `Date Start DAA Treatment`)
)
readr::write_csv(cirrhosis_audit, file.path(OUTPUT_DIR, "cirrhosis_exclusion_audit.csv"))

# ---- 6. HCV-antibody-positive descriptive population -------------------------
# Cirrhosis is intentionally not excluded here: this is the broader descriptive
# population specified in the thesis methods.
antibody_positive <- bind_rows(
  control |> mutate(id = No, group = "Control"),
  intervention |> mutate(id = record_id, group = "Intervention")
) |>
  filter(`RTK Antibody Result` == "Positive") |>
  mutate(group = factor(group, c("Control", "Intervention")))

baseline_variables <- c(
  "Age", "Sex", "Nationality", "Race", "Education", "Employment",
  "Marital status", "Housing", "HIV co-infection", "Hx of Incarceration",
  "Current Methadone", "Opioids", "ATS", "Benzos", "Cannabis", "Others",
  "Recent drug use", "Hx of injecting drug use"
)

baseline_missingness <- antibody_positive |>
  summarise(across(all_of(baseline_variables), ~ sum(is.na(.x)))) |>
  pivot_longer(everything(), names_to = "variable", values_to = "n_missing")
readr::write_csv(baseline_missingness,
                 file.path(OUTPUT_DIR, "baseline_missingness.csv"))

# A long-format descriptive table preserves exact values and can be formatted in
# the thesis document without tying the analysis to a particular table package.
baseline_continuous <- antibody_positive |>
  group_by(group) |>
  summarise(
    variable = "Age", level = NA_character_, n = sum(!is.na(Age)),
    estimate = sprintf("%.1f (%.1f)", mean(Age, na.rm = TRUE), sd(Age, na.rm = TRUE)),
    p_value = wilcox.test(Age ~ group, data = antibody_positive)$p.value,
    .groups = "drop"
  )

categorical_variables <- setdiff(baseline_variables, "Age")
baseline_categorical <- lapply(categorical_variables, function(v) {
  tab <- table(antibody_positive[[v]], antibody_positive$group, useNA = "no")
  expected <- suppressWarnings(chisq.test(tab)$expected)
  p <- if (any(expected < 5)) fisher.test(tab)$p.value else chisq.test(tab)$p.value
  antibody_positive |>
    filter(!is.na(.data[[v]])) |>
    count(group, level = .data[[v]], name = "n") |>
    group_by(group) |>
    mutate(
      variable = v,
      estimate = sprintf("%d (%.1f%%)", n, 100 * n / sum(n)),
      p_value = p
    ) |>
    ungroup() |>
    transmute(group, variable, level = as.character(level), n, estimate, p_value)
}) |>
  bind_rows()

baseline_table <- bind_rows(baseline_continuous, baseline_categorical)
readr::write_csv(baseline_table, file.path(OUTPUT_DIR, "baseline_table.csv"))

# Conditional denominators follow the immediately preceding eligible stage.
cascade <- antibody_positive |>
  group_by(group) |>
  summarise(
    antibody_positive_n = n(),
    RNA_tested_n = sum(!is.na(`HCV RNA Result`)),
    RNA_detected_n = sum(`HCV RNA Result` == "Detected", na.rm = TRUE),
    DAA_started_n = sum(`HCV RNA Result` == "Detected" &
                          !is.na(`Date Start DAA Treatment`)),
    DAA_completed_n = sum(`HCV RNA Result` == "Detected" &
                            !is.na(`Date Start DAA Treatment`) &
                            `DAA Completed` == "Yes", na.rm = TRUE),
    RNA_tested_pct = 100 * RNA_tested_n / antibody_positive_n,
    RNA_detected_pct = 100 * RNA_detected_n / RNA_tested_n,
    DAA_started_pct = 100 * DAA_started_n / RNA_detected_n,
    DAA_completed_pct = 100 * DAA_completed_n / DAA_started_n,
    .groups = "drop"
  )
readr::write_csv(cascade, file.path(OUTPUT_DIR, "care_cascade.csv"))

# ---- 7. Construct and lock the primary mITT cohort ----------------------------
build_mitt <- function(data, id_column, group_label) {
  data |>
    filter(
      `RTK Antibody Result` == "Positive",
      `HCV RNA Result` == "Detected",
      Cirrhosis == "No"                         # PRIMARY EXCLUSION
    ) |>
    transmute(
      id = .data[[id_column]], group = group_label,
      rtk_date = as_date_any(`Date RTK Antibody`),
      daa_date = as_date_any(`Date Start DAA Treatment`),
      Cirrhosis, Age, Employment, Housing, `Marital status`,
      `HIV co-infection`, `Current Methadone`, `Hx of Incarceration`,
      `Recent drug use`, `Hx of injecting drug use`
    ) |>
    filter(!is.na(rtk_date)) |>
    mutate(
      raw_time = as.numeric(daa_date - rtk_date),
      invalid_negative_interval = !is.na(raw_time) & raw_time < 0
    ) |>
    filter(!invalid_negative_interval) |>
    mutate(
      event = as.integer(!is.na(raw_time) & raw_time <= FOLLOWUP_DAYS),
      time_days = if_else(event == 1L, raw_time, as.numeric(FOLLOWUP_DAYS)),
      never_started = is.na(daa_date),
      started_after_365 = !is.na(raw_time) & raw_time > FOLLOWUP_DAYS
    ) |>
    select(-raw_time, -invalid_negative_interval)
}

mitt <- bind_rows(
  build_mitt(control, "No", "Control"),
  build_mitt(intervention, "record_id", "Intervention")
) |>
  mutate(
    group = factor(group, c("Control", "Intervention")),
    Cirrhosis = factor(Cirrhosis, c("No", "Yes"))
  ) |>
  droplevels()

# Non-negotiable proof that cirrhosis cannot enter any subsequent analysis.
stopifnot(nrow(mitt) > 0L)
stopifnot(all(as.character(mitt$Cirrhosis) == "No"))
excluded_keys <- cirrhosis_audit |>
  transmute(id = as.character(id), group)
included_excluded_overlap <- mitt |>
  transmute(id = as.character(id), group = as.character(group)) |>
  inner_join(excluded_keys, by = c("id", "group"))
stopifnot(nrow(included_excluded_overlap) == 0L)

mitt_signature <- mitt |>
  count(group, event, name = "n") |>
  arrange(group, event)
readr::write_csv(mitt_signature, file.path(OUTPUT_DIR, "mitt_cohort_lock.csv"))
saveRDS(mitt, file.path(OUTPUT_DIR, "primary_mitt_cirrhosis_excluded.rds"))

adjustment_variables <- c(
  "Age", "Employment", "Housing", "Marital status", "HIV co-infection",
  "Current Methadone", "Hx of Incarceration", "Recent drug use",
  "Hx of injecting drug use"
)
model_data <- mitt |>
  drop_na(time_days, event, group, all_of(adjustment_variables))
stopifnot(all(as.character(model_data$Cirrhosis) == "No"))

# ---- 8. Nelson-Aalen cumulative-incidence analysis ----------------------------
na_estimates <- nelson_aalen_by_group(mitt)
readr::write_csv(na_estimates,
                 file.path(OUTPUT_DIR, "nelson_aalen_cumulative_incidence.csv"))

na_figure <- ggplot(na_estimates,
                    aes(time, cumulative_incidence, colour = group, fill = group)) +
  geom_step(linewidth = 0.8) +
  geom_ribbon(aes(ymin = lower_95_CI, ymax = upper_95_CI),
              alpha = 0.15, colour = NA) +
  coord_cartesian(xlim = c(0, FOLLOWUP_DAYS), ylim = c(0, 1)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Days from HCV antibody-positive date",
       y = "Cumulative probability of DAA initiation",
       colour = "Study group", fill = "Study group") +
  theme_classic()
ggsave(file.path(OUTPUT_DIR, "nelson_aalen_cumulative_incidence.png"),
       na_figure, width = 7.5, height = 5.5, dpi = 300)

na_medians <- na_median(mitt)
set.seed(BOOTSTRAP_SEED)
bootstrap_medians <- lapply(seq_len(BOOTSTRAP_REPS), function(b) {
  mitt |>
    group_by(group) |>
    group_modify(~ slice_sample(.x, n = nrow(.x), replace = TRUE)) |>
    ungroup() |>
    na_median() |>
    mutate(replicate = b)
}) |>
  bind_rows()

na_median_results <- na_medians |>
  left_join(
    bootstrap_medians |>
      group_by(group) |>
      summarise(
        bootstrap_valid_n = sum(!is.na(median_days)),
        lower_95_CI = quantile(median_days, 0.025, na.rm = TRUE),
        upper_95_CI = quantile(median_days, 0.975, na.rm = TRUE),
        .groups = "drop"
      ),
    by = "group"
  )
readr::write_csv(na_median_results,
                 file.path(OUTPUT_DIR, "nelson_aalen_median_bootstrap.csv"))

logrank <- survdiff(Surv(time_days, event) ~ group, data = mitt)
logrank_result <- tibble(
  chi_square = unname(logrank$chisq),
  df = length(logrank$n) - 1L,
  p_value = pchisq(logrank$chisq, length(logrank$n) - 1L, lower.tail = FALSE)
)
readr::write_csv(logrank_result, file.path(OUTPUT_DIR, "logrank_test.csv"))

# ---- 9. Primary Cox models and proportional-hazards diagnostics ---------------
cox_formula <- Surv(time_days, event) ~
  group + Age + Employment + Housing + `Marital status` +
  `HIV co-infection` + `Current Methadone` + `Hx of Incarceration` +
  `Recent drug use` + `Hx of injecting drug use`

cox_crude <- coxph(Surv(time_days, event) ~ group, data = model_data,
                   ties = "efron", x = TRUE)
cox_adjusted <- coxph(cox_formula, data = model_data,
                      ties = "efron", x = TRUE, na.action = na.exclude)

cox_results <- bind_rows(
  extract_cox_group(cox_crude, "Crude"),
  extract_cox_group(cox_adjusted, "Adjusted")
)
readr::write_csv(cox_results, file.path(OUTPUT_DIR, "cox_primary_results.csv"))

ph_adjusted <- cox.zph(cox_adjusted)
ph_table <- as.data.frame(ph_adjusted$table) |>
  tibble::rownames_to_column("term")
readr::write_csv(ph_table, file.path(OUTPUT_DIR, "cox_PH_tests_365_day.csv"))

png(file.path(OUTPUT_DIR, "cox_schoenfeld_plots_365_day.png"),
    width = 3000, height = 2400, res = 300)
plot(ph_adjusted)
dev.off()

# Events, person-time and Poisson confidence intervals.
incidence_rates <- mitt |>
  group_by(group) |>
  summarise(events = sum(event), person_days = sum(time_days),
            person_years = person_days / 365.25, .groups = "drop") |>
  rowwise() |>
  mutate(
    rate_per_100_person_years = 100 * events / person_years,
    lower_95_CI = 100 * poisson.test(events, T = person_years)$conf.int[1],
    upper_95_CI = 100 * poisson.test(events, T = person_years)$conf.int[2]
  ) |>
  ungroup()
readr::write_csv(incidence_rates, file.path(OUTPUT_DIR, "incidence_rates.csv"))

# ---- 10. Thirty-day adjusted Cox model ---------------------------------------
model_data_30 <- model_data |>
  mutate(
    event_30 = as.integer(event == 1L & time_days <= 30),
    time_30 = pmin(time_days, 30)
  )
stopifnot(all(as.character(model_data_30$Cirrhosis) == "No"))

cox_30_formula <- Surv(time_30, event_30) ~
  group + Age + Employment + Housing + `Marital status` +
  `HIV co-infection` + `Current Methadone` + `Hx of Incarceration` +
  `Recent drug use` + `Hx of injecting drug use`
cox_adjusted_30 <- coxph(cox_30_formula, data = model_data_30,
                         ties = "efron", x = TRUE, na.action = na.exclude)
cox_30_result <- extract_cox_group(cox_adjusted_30, "Adjusted: days 0-30")
readr::write_csv(cox_30_result,
                 file.path(OUTPUT_DIR, "cox_adjusted_30_day.csv"))

ph_30 <- cox.zph(cox_adjusted_30)
readr::write_csv(
  as.data.frame(ph_30$table) |> tibble::rownames_to_column("term"),
  file.path(OUTPUT_DIR, "cox_PH_tests_30_day.csv")
)
png(file.path(OUTPUT_DIR, "cox_schoenfeld_plots_30_day.png"),
    width = 3000, height = 2400, res = 300)
plot(ph_30)
dev.off()

# ---- 11. Restricted mean survival time to day 365 -----------------------------
arm <- as.integer(model_data$group == "Intervention")
rmst_covariates <- model.matrix(
  ~ Age + Employment + Housing + `Marital status` + `HIV co-infection` +
    `Current Methadone` + `Hx of Incarceration` + `Recent drug use` +
    `Hx of injecting drug use`,
  data = model_data
)[, -1, drop = FALSE]
varying <- apply(rmst_covariates, 2, function(x) length(unique(x)) > 1L)
rmst_covariates <- rmst_covariates[, varying, drop = FALSE]

rmst_unadjusted <- survRM2::rmst2(
  model_data$time_days, model_data$event, arm, tau = FOLLOWUP_DAYS
)
rmst_adjusted <- survRM2::rmst2(
  model_data$time_days, model_data$event, arm, tau = FOLLOWUP_DAYS,
  covariates = rmst_covariates
)
saveRDS(list(unadjusted = rmst_unadjusted, adjusted = rmst_adjusted),
        file.path(OUTPUT_DIR, "rmst_365_day_results.rds"))

capture.output(rmst_unadjusted, rmst_adjusted,
               file = file.path(OUTPUT_DIR, "rmst_365_day_results.txt"))

# ---- 12. RS-DAA implementation outcomes --------------------------------------
rsdaa_mitt <- filter(mitt, group == "Intervention")
rapid_start <- lapply(RAPID_DAYS, function(day) {
  successes <- sum(rsdaa_mitt$event == 1L & rsdaa_mitt$time_days <= day)
  test <- binom.test(successes, nrow(rsdaa_mitt))
  tibble(
    day = day, numerator = successes, denominator = nrow(rsdaa_mitt),
    proportion = unname(test$estimate), lower_95_CI = test$conf.int[1],
    upper_95_CI = test$conf.int[2]
  )
}) |>
  bind_rows()
readr::write_csv(rapid_start, file.path(OUTPUT_DIR, "rapid_start_feasibility.csv"))

# Penetration and fidelity source counts were manually verified by reconciling
# the RS-DAA masterlist, Panda treatment record and laboratory records. Only
# aggregate counts are reproduced here; identifiable source records are not
# included in this repository.
implementation_source_flow <- tibble(
  source = c("RS-DAA masterlist", "Panda treatment record", "Laboratory records"),
  identified_n = c(84L, 7L, 3L)
)
stopifnot(sum(implementation_source_flow$identified_n) == 94L)

penetration <- tibble(
  identified_n = 94L,
  ineligible_n = 2L,
  eligible_and_offered_n = 92L,
  enrolled_n = 84L,
  conventional_SOC_n = 8L
) |>
  mutate(penetration = enrolled_n / eligible_and_offered_n)
stopifnot(
  penetration$identified_n - penetration$ineligible_n ==
    penetration$eligible_and_offered_n,
  penetration$eligible_and_offered_n - penetration$enrolled_n ==
    penetration$conventional_SOC_n
)

fidelity <- tibble(
  indicator = c(
    "Bundled pretreatment assessment",
    "Community blood collection",
    "Telemedicine consultation",
    "CHW-supported medication delivery"
  ),
  numerator = c(84L, 73L, 36L, 33L),
  denominator = c(84L, 84L, 54L, 36L)
) |>
  mutate(proportion = numerator / denominator)

readr::write_csv(implementation_source_flow,
                 file.path(OUTPUT_DIR, "implementation_source_flow.csv"))
readr::write_csv(penetration,
                 file.path(OUTPUT_DIR, "implementation_penetration.csv"))
readr::write_csv(fidelity,
                 file.path(OUTPUT_DIR, "implementation_fidelity.csv"))

# ---- 13. IPTW sensitivity analysis -------------------------------------------
iptw_data <- model_data |>
  mutate(treatment = as.integer(group == "Intervention"))

ps_formula <- treatment ~
  Age + Employment + Housing + `Marital status` + `HIV co-infection` +
  `Current Methadone` + `Hx of Incarceration` + `Recent drug use` +
  `Hx of injecting drug use`
ps_model <- glm(ps_formula, data = iptw_data, family = binomial())
iptw_data <- iptw_data |>
  mutate(
    propensity_score = predict(ps_model, type = "response"),
    marginal_treatment = mean(treatment),
    stabilised_weight = if_else(
      treatment == 1L,
      marginal_treatment / propensity_score,
      (1 - marginal_treatment) / (1 - propensity_score)
    )
  )
stopifnot(all(iptw_data$propensity_score > 0 & iptw_data$propensity_score < 1))

weight_limits <- quantile(iptw_data$stabilised_weight, c(0.01, 0.99))
iptw_data <- iptw_data |>
  mutate(weight = pmin(pmax(stabilised_weight, weight_limits[1]), weight_limits[2]))
effective_sample_size <- with(iptw_data, sum(weight)^2 / sum(weight^2))

balance <- bind_rows(
  tibble(
    variable = "Age", level = NA_character_,
    SMD_unweighted = weighted_smd_continuous(iptw_data$Age,
                                             iptw_data$treatment),
    SMD_weighted = weighted_smd_continuous(iptw_data$Age,
                                           iptw_data$treatment,
                                           iptw_data$weight)
  ),
  lapply(setdiff(adjustment_variables, "Age"), function(v) {
    level <- levels(iptw_data[[v]])[2]
    tibble(
      variable = v, level = level,
      SMD_unweighted = weighted_smd_binary(
        iptw_data[[v]], iptw_data$treatment, level
      ),
      SMD_weighted = weighted_smd_binary(
        iptw_data[[v]], iptw_data$treatment, level, iptw_data$weight
      )
    )
  }) |> bind_rows()
) |>
  mutate(absolute_SMD_weighted = abs(SMD_weighted),
         balanced_below_0_10 = absolute_SMD_weighted < 0.10)
readr::write_csv(balance, file.path(OUTPUT_DIR, "IPTW_covariate_balance.csv"))

ps_plot <- ggplot(iptw_data, aes(propensity_score, fill = group)) +
  geom_histogram(bins = 20, position = "identity", alpha = 0.5) +
  labs(x = "Propensity score", y = "Count", fill = "Study group") +
  theme_classic()
ggsave(file.path(OUTPUT_DIR, "IPTW_propensity_score_overlap.png"),
       ps_plot, width = 7, height = 5, dpi = 300)

cox_iptw <- coxph(
  Surv(time_days, event) ~ group, data = iptw_data, weights = weight,
  robust = TRUE, ties = "efron", x = TRUE
)
cox_iptw_doubly_adjusted <- coxph(
  cox_formula, data = iptw_data, weights = weight,
  robust = TRUE, ties = "efron", x = TRUE
)
iptw_results <- bind_rows(
  extract_cox_group(cox_iptw, "IPTW weighted"),
  extract_cox_group(cox_iptw_doubly_adjusted, "IPTW doubly adjusted")
)
readr::write_csv(iptw_results, file.path(OUTPUT_DIR, "IPTW_cox_results.csv"))
readr::write_csv(
  tibble(
    n = nrow(iptw_data), effective_sample_size = effective_sample_size,
    lower_truncation = weight_limits[1], upper_truncation = weight_limits[2],
    min_weight = min(iptw_data$weight), max_weight = max(iptw_data$weight)
  ),
  file.path(OUTPUT_DIR, "IPTW_weight_summary.csv")
)

# ---- 14. Final reproducibility locks and record -------------------------------
# Re-run these checks at the end, after every model has been fitted.
stopifnot(all(as.character(mitt$Cirrhosis) == "No"))
stopifnot(all(as.character(model_data$Cirrhosis) == "No"))
stopifnot(all(as.character(model_data_30$Cirrhosis) == "No"))
stopifnot(all(as.character(iptw_data$Cirrhosis) == "No"))
stopifnot(identical(
  mitt_signature,
  mitt |> count(group, event, name = "n") |> arrange(group, event)
))

writeLines(capture.output(sessionInfo()), file.path(OUTPUT_DIR, "sessionInfo.txt"))
writeLines(
  c(
    "Primary mITT definition: antibody positive, RNA detected, valid antibody-test date,",
    "eligible for clinic initiation, and NO cirrhosis requiring specialist referral.",
    paste("mITT N =", nrow(mitt)),
    paste("Complete-case model N =", nrow(model_data)),
    paste("Bootstrap seed =", BOOTSTRAP_SEED),
    paste("Bootstrap replicates =", BOOTSTRAP_REPS)
  ),
  file.path(OUTPUT_DIR, "analysis_readme.txt")
)

cat("Analysis completed. Outputs saved in:", normalizePath(OUTPUT_DIR), "\n")
