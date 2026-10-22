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
CONTROL_FILE <- Sys.getenv("RSDAA_CONTROL_FILE", "path/to/control_masterlist.xlsx")
INTERVENTION_FILE <- Sys.getenv("RSDAA_INTERVENTION_FILE", "path/to/redcap_intervention_export.xlsx")
IMPLEMENTATION_FILE <- Sys.getenv("RSDAA_IMPLEMENTATION_FILE", "path/to/implementation_masterlist.xlsx")
OUTPUT_DIR <- Sys.getenv("RSDAA_OUTPUT_DIR", "RSDAA_thesis_outputs")
CONTROL_BASELINE_SHEET <- 1L
CONTROL_SVR_SHEET <- 2L
# Outputs contain participant-level audits: keep the output directory private.

FOLLOWUP_DAYS <- 365L
RAPID_DAYS <- c(7L, 14L, 30L)
BOOTSTRAP_REPS <- 5000L
BOOTSTRAP_SEED <- 20260817L

required_packages <- c(
  "readxl", "readr", "dplyr", "tidyr", "stringr", "lubridate",
  "survival", "survRM2", "ggplot2", "scales", "survminer"
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
 if (inherits(x, "POSIXt")) return(as.Date(x))
 z <- trimws(as.character(x)); z[z %in% c("", "x", "X", "999", "NA")] <- NA
 out <- as.Date(rep(NA_character_, length(z)))
 iso <- !is.na(z) & grepl("^[0-9]{4}-", z)
 out[iso] <- as.Date(substr(z[iso], 1, 10), format = "%Y-%m-%d")
 serial <- suppressWarnings(as.numeric(z))
 excel <- !is.na(serial) & serial > 20000 & serial < 70000
 out[excel] <- as.Date(serial[excel], origin = "1899-12-30")
 rest <- !is.na(z) & !iso & !excel
 out[rest] <- as.Date(z[rest], format = "%d/%m/%Y")
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

# Decimal half-up formatting matches the thesis baseline table (e.g. 48.55 -> 48.6).
format_age <- function(x) {
 r1 <- function(v) floor(v*10 + 0.5 + 1e-10)/10
 sprintf("%.1f (%.1f)", r1(mean(x, na.rm = TRUE)), r1(sd(x, na.rm = TRUE)))
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
control_raw <- readxl::read_excel(CONTROL_FILE, sheet = CONTROL_BASELINE_SHEET)
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

control_raw <- control_raw |> mutate(source_id = sprintf("C%03d", row_number()))
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
intervention_raw <- if (tolower(tools::file_ext(INTERVENTION_FILE)) %in% c("xlsx", "xls")) {
  readxl::read_excel(INTERVENTION_FILE)
} else readr::read_csv(INTERVENTION_FILE, show_col_types = FALSE)
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
stopifnot(!anyDuplicated(control$source_id), !anyDuplicated(intervention$record_id))
# Historical No repeats across calendar years; it is not a unique cohort key.
stopifnot(!any(is.na(control$Cirrhosis)), !any(is.na(intervention$Cirrhosis)))
stopifnot(!any(intervention$Cirrhosis == "Yes" &
                intervention$`HCV RNA Result` != "Detected"))

cirrhosis_audit <- bind_rows(
  control |>
    filter(Cirrhosis == "Yes") |>
    transmute(id = source_id, group = "Control", `HCV RNA Result`,
              `Reason Not Started`, `Date Start DAA Treatment`),
  intervention |>
    filter(Cirrhosis == "Yes") |>
    transmute(id = as.character(record_id), group = "Intervention", `HCV RNA Result`,
              `Reason Not Started`, `Date Start DAA Treatment`)
)
readr::write_csv(cirrhosis_audit, file.path(OUTPUT_DIR, "cirrhosis_exclusion_audit.csv"))

# ---- 6. HCV-antibody-positive descriptive population -------------------------
# Cirrhosis is intentionally not excluded here: this is the broader descriptive
# population specified in the thesis methods.
antibody_positive <- bind_rows(
  control |> mutate(id = source_id, group = "Control"),
  intervention |> mutate(id = as.character(record_id), group = "Intervention")
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
    estimate = format_age(Age),
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

# Full descriptive cascade; historical SVR fields come from worksheet 2 only.
# Retained early tests are the study operational classification, not strict SVR12.
SVR_OUTPUT_DIR <- file.path(OUTPUT_DIR, "SVR_audits")
dir.create(SVR_OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
yn_svr <- function(x) { z <- yes_no(x); ifelse(is.na(z), NA, z == "Yes") }
## CONTROL: match source rows to the already-cleaned 84-person cohort.
cs <- readxl::read_excel(CONTROL_FILE, sheet = CONTROL_SVR_SHEET)
key_svr <- function(x) gsub("[^0-9A-Za-z]", "", as.character(x))
stopifnot(!anyDuplicated(key_svr(cs$IC)), !anyDuplicated(key_svr(control$IC)))
match_svr <- match(key_svr(control$IC), key_svr(cs$IC))
stopifnot(!anyNA(match_svr), nrow(control) == 84L)
cs <- cs[match_svr, ]
control_svr <- control %>% transmute(
 group = "Historical SOC", source_row = match_svr + 1L,
 id = as.character(No), antibody = `RTK Antibody Result` == "Positive",
 rna_tested = `HCV RNA Result` %in% c("Detected", "Not Detected"),
 rna_detected = `HCV RNA Result` == "Detected",
 referred = as.character(Cirrhosis) == "Yes",
 start_date = as_date_any(`Date Start DAA Treatment`),
 completed = yn_svr(`DAA Completed`),
 tested = yn_svr(cs$`SVR tested`),
 test_date = as_date_any(cs$`Date SVR tested`),
 cured = yn_svr(cs$`SVR achieved`)
) %>% mutate(
 ## Historical completion is ESTIMATED, not a verified actual date.
 completion_date = start_date + 84L,
 completion_basis = "Estimated: start plus 84 days"
)
## User-verified source corrections, applied here without editing the workbook.
control_svr <- control_svr %>% mutate(
 test_date = case_when(
  id == "12" & start_date == as.Date("2023-08-05") ~ as.Date("2024-01-30"),
  id == "29" & start_date == as.Date("2024-05-03") ~ as.Date("2024-10-17"),
  TRUE ~ test_date
 )
)

svr_cascade_audit <- control_svr %>% mutate(
 across(c(antibody, rna_tested, rna_detected, referred, completed, tested),
        ~ coalesce(.x, FALSE)),
 eligible = antibody & rna_detected & !referred,
 initiated = eligible & !is.na(start_date),
 completed_stage = initiated & completed,
 days_post_completion = as.numeric(test_date - completion_date),
 weeks_post_completion = days_post_completion / 7,
 tested_completer = completed_stage & tested,
 assessed = tested_completer & !is.na(days_post_completion) &
            ((days_post_completion >= 84 & days_post_completion <= 168) |
             (group == "Historical SOC" & days_post_completion == 83) |
             (group == "RS-DAA" & days_post_completion %in% 75:81)),
 achieved = assessed & coalesce(cured, FALSE),
 timing = case_when(
  !tested_completer ~ "Not tested among completers",
  is.na(days_post_completion) ~ "Missing date: review",
  days_post_completion < 0 ~ "Before completion: review",
  days_post_completion < 84 & assessed ~ "Early: retained",
  days_post_completion < 84 ~ "Other early test: review / excluded",
  days_post_completion <= 168 ~ "12–24 weeks: retained",
  TRUE ~ "Beyond 24 weeks: excluded"
 )
)
make_cascade_svr <- function(d) {
 n <- c(sum(d$antibody), sum(d$antibody & d$rna_tested),
        sum(d$antibody & d$rna_detected), sum(d$eligible),
        sum(d$initiated), sum(d$completed_stage), sum(d$assessed), sum(d$achieved))
 denom <- c(n[1], n[1], n[2], n[3], n[4], n[5], n[6], n[7])
 tibble(group = unique(d$group), Stage = c("HCV antibody positive", "HCV RNA tested",
  "HCV RNA detected", "Eligible for clinic-based DAA initiation",
  "DAA treatment initiated", "DAA treatment completed",
  "SVR12 assessed (study operational classification)",
  "SVR12 achieved (study operational classification)"),
  n = n, N = denom, percentage = ifelse(denom > 0, 100*n/denom, NA_real_),
  Result = sprintf("%d/%d (%.1f%%)", n, N, percentage))
}
final_cascade_svr <- bind_rows(lapply(split(svr_cascade_audit,
                                          svr_cascade_audit$group), make_cascade_svr))
svr_timing_summary <- svr_cascade_audit %>% filter(tested_completer) %>%
 group_by(group, timing) %>% summarise(n = n(), cured = sum(cured %in% TRUE), .groups = "drop")
## Validation only: these expected counts never set the calculated results.
expected_svr <- c(84,84,72,60,57,47,20,19)
if (!identical(as.integer(final_cascade_svr$n), as.integer(expected_svr)))
 warning("Calculated cascade differs from thesis counts. Review the source and date audit.")
print(final_cascade_svr, n = Inf)
print(svr_timing_summary)
write.csv(final_cascade_svr, file.path(SVR_OUTPUT_DIR, "Control_Cascade_With_SVR.csv"), row.names = FALSE)
write.csv(svr_cascade_audit %>% select(group, id, source_row, start_date,
 completion_date, completion_basis, test_date, days_post_completion,
 weeks_post_completion, tested, cured, assessed, achieved, timing),
 file.path(SVR_OUTPUT_DIR, "Control_SVR_Date_Audit.csv"), row.names = FALSE, na = "")
write.csv(svr_timing_summary, file.path(SVR_OUTPUT_DIR, "Control_SVR_Timing_Summary.csv"), row.names = FALSE)

control_cascade <- final_cascade_svr
## INTERVENTION: retain treatment-status rows only for cleaned cohort IDs.
is <- intervention_raw %>%
 filter(redcap_event_name == "status_rawatan_hep_arm_1",
        record_id %in% intervention$record_id)
stopifnot(!anyDuplicated(is$record_id), nrow(is) == 84L)
intervention_svr <- intervention %>% transmute(
 record_id, id = as.character(record_id),
 antibody = `RTK Antibody Result` == "Positive",
 referred = as.character(Cirrhosis) == "Yes"
) %>% left_join(is %>% transmute(
 record_id, rna_tested = yn_svr(rna_tested),
 rna_detected = as.character(rna_status) == "1",
 prescribed = yn_svr(daa_prescribe),
 start_date = as_date_any(daa_dispensed),
 completed = yn_svr(daa_treatment_status),
 completion_date = as_date_any(daa_completed),
 tested = yn_svr(svr_test), test_date = as_date_any(svr_date),
 cured = yn_svr(svr_cure)
), by = "record_id") %>% mutate(
 group = "RS-DAA", source_row = NA_integer_,
 completion_basis = "Recorded completion date"
) %>% select(-record_id)

svr_cascade_audit <- intervention_svr %>% mutate(
 across(c(antibody, rna_tested, rna_detected, referred, completed, tested),
        ~ coalesce(.x, FALSE)),
 eligible = antibody & rna_detected & !referred,
 initiated = eligible & !is.na(start_date),
 completed_stage = initiated & completed,
 days_post_completion = as.numeric(test_date - completion_date),
 weeks_post_completion = days_post_completion / 7,
 tested_completer = completed_stage & tested,
 assessed = tested_completer & !is.na(days_post_completion) &
            ((days_post_completion >= 84 & days_post_completion <= 168) |
             (group == "Historical SOC" & days_post_completion == 83) |
             (group == "RS-DAA" & days_post_completion %in% 75:81)),
 achieved = assessed & coalesce(cured, FALSE),
 timing = case_when(
  !tested_completer ~ "Not tested among completers",
  is.na(days_post_completion) ~ "Missing date: review",
  days_post_completion < 0 ~ "Before completion: review",
  days_post_completion < 84 & assessed ~ "Early: retained",
  days_post_completion < 84 ~ "Other early test: review / excluded",
  days_post_completion <= 168 ~ "12–24 weeks: retained",
  TRUE ~ "Beyond 24 weeks: excluded"
 )
)
make_cascade_svr <- function(d) {
 n <- c(sum(d$antibody), sum(d$antibody & d$rna_tested),
        sum(d$antibody & d$rna_detected), sum(d$eligible),
        sum(d$initiated), sum(d$completed_stage), sum(d$assessed), sum(d$achieved))
 denom <- c(n[1], n[1], n[2], n[3], n[4], n[5], n[6], n[7])
 tibble(group = unique(d$group), Stage = c("HCV antibody positive", "HCV RNA tested",
  "HCV RNA detected", "Eligible for clinic-based DAA initiation",
  "DAA treatment initiated", "DAA treatment completed",
  "SVR12 assessed (study operational classification)",
  "SVR12 achieved (study operational classification)"),
  n = n, N = denom, percentage = ifelse(denom > 0, 100*n/denom, NA_real_),
  Result = sprintf("%d/%d (%.1f%%)", n, N, percentage))
}
final_cascade_svr <- bind_rows(lapply(split(svr_cascade_audit,
                                          svr_cascade_audit$group), make_cascade_svr))
svr_timing_summary <- svr_cascade_audit %>% filter(tested_completer) %>%
 group_by(group, timing) %>% summarise(n = n(), cured = sum(cured %in% TRUE), .groups = "drop")
## Validation only: these expected counts never set the calculated results.
expected_svr <- c(84,84,66,58,54,47,34,31)
if (!identical(as.integer(final_cascade_svr$n), as.integer(expected_svr)))
 warning("Calculated cascade differs from thesis counts. Review the source and date audit.")
print(final_cascade_svr, n = Inf)
print(svr_timing_summary)
write.csv(final_cascade_svr, file.path(SVR_OUTPUT_DIR, "Intervention_Cascade_With_SVR.csv"), row.names = FALSE)
write.csv(svr_cascade_audit %>% select(group, id, source_row, start_date,
 completion_date, completion_basis, test_date, days_post_completion,
 weeks_post_completion, tested, cured, assessed, achieved, timing),
 file.path(SVR_OUTPUT_DIR, "Intervention_SVR_Date_Audit.csv"), row.names = FALSE, na = "")
write.csv(svr_timing_summary, file.path(SVR_OUTPUT_DIR, "Intervention_SVR_Timing_Summary.csv"), row.names = FALSE)

intervention_cascade <- final_cascade_svr
cascade <- bind_rows(control_cascade, intervention_cascade)
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
      id = as.character(.data[[id_column]]), group = group_label,
      rtk_date = as_date_any(`Date RTK Antibody`),
      daa_date = as_date_any(`Date Start DAA Treatment`),
      Cirrhosis, across(all_of(baseline_variables)),
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
  build_mitt(control, "source_id", "Control"),
  build_mitt(intervention, "record_id", "Intervention")
) |>
  mutate(
    group = factor(group, c("Control", "Intervention")),
    Cirrhosis = factor(Cirrhosis, c("No", "Yes"))
  ) |>
  droplevels()

# Non-negotiable proof that cirrhosis cannot enter any subsequent analysis.
stopifnot(nrow(mitt) == 118L, sum(mitt$event) == 106L,
          sum(mitt$group == "Control") == 60L,
          sum(mitt$group == "Intervention") == 58L,
          sum(mitt$event[mitt$group == "Control"]) == 52L,
          sum(mitt$event[mitt$group == "Intervention"]) == 54L)
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

# Second baseline table: final RNA-detected, clinic-eligible 118-person cohort.
# Fisher exact tests for categorical variables, matching the parallel script.
baseline_118 <- bind_rows(lapply(baseline_variables, function(v) {
 if (v == "Age") {
  p <- wilcox.test(Age ~ group, data = mitt)$p.value
  return(mitt |> group_by(group) |> summarise(variable = v,
   level = NA_character_, n = sum(!is.na(Age)),
   estimate = format_age(Age), p_value = p,
   .groups = "drop"))
 }
 tab <- table(mitt[[v]], mitt$group)
 tab <- tab[rowSums(tab) > 0, , drop = FALSE]
 p <- if (nrow(tab) > 1) fisher.test(tab)$p.value else NA_real_
 mitt |> filter(!is.na(.data[[v]])) |> count(group, level = .data[[v]], name = "n") |>
  group_by(group) |> mutate(variable = v,
   estimate = sprintf("%d (%.1f%%)", n, 100*n/sum(n)), p_value = p) |>
  ungroup() |> mutate(level = as.character(level))
}))
readr::write_csv(baseline_118, file.path(OUTPUT_DIR, "baseline_table_118.csv"))
stopifnot(nrow(model_data) == 118L)

# ---- 8. Nelson-Aalen cumulative-incidence analysis ----------------------------
na_estimates <- nelson_aalen_by_group(mitt)
readr::write_csv(na_estimates,
                 file.path(OUTPUT_DIR, "nelson_aalen_cumulative_incidence.csv"))

# stype=2/ctype=1: exponential of explicit Nelson–Aalen hazard, not KM.
# Log-survival intervals equal normal cumulative-hazard intervals transformed
# to incidence (lower hazard bound truncated at zero). Bootstrap is median-only.
na_curve_fit_nc <- survfit(
  Surv(time_days, event) ~ group,
  data = mitt,
  stype = 2, ctype = 1,
  conf.int = 0.95, conf.type = "log"
)
# stype=2 gives exp(-H_NA); fun="event" gives 1-exp(-H_NA).
# log-survival confidence limits correspond to normal limits on H_NA,
# with a nonnegative lower hazard bound; transformation reverses survival bounds.
na_curve_plot_nc <- survminer::ggsurvplot(
  na_curve_fit_nc, data = mitt, fun = "event",
  conf.int = TRUE, conf.int.style = "ribbon",
  risk.table = "nrisk_cumevents",
  risk.table.title = "Number at risk (cumulative events)",
  risk.table.fontsize = 4, risk.table.height = 0.23,
  risk.table.y.text.col = TRUE,
  censor = TRUE, pval = FALSE,
  xlab = "Days from HCV antibody-positive date",
  ylab = "Cumulative incidence of DAA initiation",
  xlim = c(0, 365), ylim = c(0, 1), break.time.by = 30,
  legend.title = "Study group",
  legend.labs = c("Historical control", "RS-DAA intervention"),
  palette = c("#F8766D", "#00BFC4"),
  ggtheme = theme_classic(base_size = 13)
)
na_curve_plot_nc$plot <- na_curve_plot_nc$plot + theme(
  legend.position = "top",
  axis.title = element_text(size = 15),
  axis.text = element_text(size = 12)
)
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
OUT_TRANSFORMED_NA_NC <- file.path(
  OUTPUT_DIR,
  "Figure_Transformed_Nelson_Aalen_Cumulative_Probability.png"
)
png(OUT_TRANSFORMED_NA_NC, width = 4200, height = 3000, res = 300)
tryCatch(print(na_curve_plot_nc, newpage = FALSE), finally = dev.off())
cat("Figure saved:", normalizePath(OUT_TRANSFORMED_NA_NC), "\n")

na_medians <- na_median(mitt)
set.seed(BOOTSTRAP_SEED)
bootstrap_medians <- bind_rows(lapply(levels(mitt$group), function(g) {
 one <- filter(mitt, group == g) |> droplevels()
 values <- replicate(BOOTSTRAP_REPS, {
  sampled <- one[sample.int(nrow(one), nrow(one), replace = TRUE), , drop = FALSE]
  na_median(sampled)$median_days
 })
 tibble(group = g, median_days = as.numeric(values), replicate = seq_len(BOOTSTRAP_REPS))
}))

readr::write_csv(bootstrap_medians, file.path(OUTPUT_DIR, "nelson_aalen_bootstrap_replicates.csv"))

na_median_results <- na_medians |>
  left_join(
    bootstrap_medians |>
      group_by(group) |>
      summarise(
        bootstrap_repetitions = BOOTSTRAP_REPS,
        bootstrap_valid_n = sum(!is.na(median_days)),
        excluded_repetitions = sum(is.na(median_days)),
        median_reached_pct = 100*mean(!is.na(median_days)),
        lower_95_CI = quantile(median_days, 0.025, na.rm = TRUE, type = 6),
        upper_95_CI = quantile(median_days, 0.975, na.rm = TRUE, type = 6),
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
    rate_per_person_year = events / person_years,
    lower_95_CI = poisson.test(events, T = person_years)$conf.int[1],
    upper_95_CI = poisson.test(events, T = person_years)$conf.int[2]
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

readr::write_csv(implementation_source_flow, file.path(OUTPUT_DIR, "implementation_source_flow.csv"))
readr::write_csv(penetration, file.path(OUTPUT_DIR, "implementation_penetration.csv"))
IMPLEMENTATION_OUTPUT_DIR <- file.path(OUTPUT_DIR, "Implementation")
EXCLUDE_IDS <- c(3L, 73L)
EXCLUDE_EXAMPLE_99 <- TRUE
EXPECTED_N <- 84L
raw <- readxl::read_excel(IMPLEMENTATION_FILE, sheet = "KK PANDAMARAN")
required <- c("Redcap No", "Blood taken where", "Mysejahtera/Whatsapp VC/Physical", "Meds")
stopifnot(all(required %in% names(raw)))
d <- data.frame(
 redcap_id = suppressWarnings(as.integer(as.character(raw[["Redcap No"]]))),
 blood_location = as.character(raw[["Blood taken where"]]),
 consultation_mode = as.character(raw[["Mysejahtera/Whatsapp VC/Physical"]]),
 medication_mode = as.character(raw[["Meds"]]), stringsAsFactors = FALSE
)
d <- d[!is.na(d$redcap_id), ]
stopifnot(!anyDuplicated(d$redcap_id))
for (nm in c("blood_location", "consultation_mode", "medication_mode")) {
 d[[nm]] <- tolower(trimws(d[[nm]]));d[[nm]][d[[nm]] == ""] <- NA_character_
}
d <- d[!d$redcap_id %in% EXCLUDE_IDS, ]
if(EXCLUDE_EXAMPLE_99) d <- d[d$redcap_id != 99L, ]
if(nrow(d) != EXPECTED_N) warning(sprintf("Cohort contains %d, not %d. Review exclusions; counts have not been forced.", nrow(d),EXPECTED_N))
known <- function(x, allowed) all(is.na(x) | x %in% allowed)
stopifnot(known(d$blood_location,c("community","clinic")),
 known(d$consultation_mode,c("mysj","mysejahtera","wa","whatsapp","whatsapp vc","physical")),
 known(d$medication_mode,c("chw","postal","self pick-up")))
d$telemedicine <- d$consultation_mode %in% c("mysj","mysejahtera","wa","whatsapp","whatsapp vc")
d$physical <- d$consultation_mode %in% "physical"
d$chw_supported_delivery <- d$medication_mode %in% c("chw", "postal")
n_review <- sum(d$telemedicine | d$physical)
n_tele <- sum(d$telemedicine)
make_row <- function(indicator, n, denom) data.frame(Indicator=indicator,n=n,N=denom,
 Percent=if(denom>0)round(100*n/denom,1) else NA_real_)
## Bundled assessment: user-verified aggregate, not inferred from blood location.
BUNDLED_ASSESSMENT_COMPLETED <- 84L
stopifnot(BUNDLED_ASSESSMENT_COMPLETED <= nrow(d))
results <- do.call(rbind,list(
 make_row("Bundled pretreatment assessment (user-verified)", BUNDLED_ASSESSMENT_COMPLETED, nrow(d)),
 make_row("Community blood collection",sum(d$blood_location %in% "community"),nrow(d)),
 make_row("Clinic blood collection",sum(d$blood_location %in% "clinic"),nrow(d)),
 make_row("Blood location missing",sum(is.na(d$blood_location)),nrow(d)),
 make_row("Telemedicine among recorded consultations",n_tele,n_review),
 make_row("Physical review among recorded consultations",sum(d$physical),n_review),
 make_row("Consultation mode missing",sum(is.na(d$consultation_mode)),nrow(d)),
 make_row("Direct CHW delivery among telemedicine participants",sum(d$telemedicine & d$medication_mode %in% "chw"),n_tele),
 make_row("Postal delivery among telemedicine participants",sum(d$telemedicine & d$medication_mode %in% "postal"),n_tele),
 make_row("Self-pickup among telemedicine participants",sum(d$telemedicine & d$medication_mode %in% "self pick-up"),n_tele),
 make_row("CHW-supported delivery (including postal) among telemedicine participants",sum(d$telemedicine & d$medication_mode %in% c("chw","postal")),n_tele)
))
## Study classification: postal delivery counts as CHW-supported delivery.
## The 54 recorded consultations are not independently proof of initiation.
## Link REDCap/dispensing data if a denominator defined by initiation is required.
audit <- d[,c("redcap_id","blood_location","consultation_mode","medication_mode","telemedicine","chw_supported_delivery")]
dir.create(IMPLEMENTATION_OUTPUT_DIR, recursive=TRUE,showWarnings=FALSE)
write.csv(results,file.path(IMPLEMENTATION_OUTPUT_DIR,"Implementation_Fidelity.csv"),row.names=FALSE,na="")
write.csv(audit,file.path(IMPLEMENTATION_OUTPUT_DIR,"Implementation_Participant_Audit.csv"),row.names=FALSE,na="")
cat("Retained participants:",nrow(d),"\n");print(results,row.names=FALSE)

# Match implementation participants to the clinical cohort; do not infer the
# clinical denominator solely from availability of a consultation entry.
stopifnot(setequal(d$redcap_id, intervention$record_id))
stopifnot(setequal(d$redcap_id[d$telemedicine | d$physical],
 intervention$record_id[!is.na(intervention$`Date Start DAA Treatment`)]))

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

# Supplementary Cox specifications from the final cirrhosis-excluded script.
cox_minimal <- coxph(Surv(time_days, event) ~ group + Age + Housing +
 `Current Methadone` + `Hx of Incarceration` + `Hx of injecting drug use`,
 data = mitt, ties = "efron", x = TRUE)
cox_full <- coxph(update(cox_formula, . ~ . + Race + Education + ATS + Benzos + Cannabis + Others),
 data = mitt, ties = "efron", x = TRUE)
readr::write_csv(bind_rows(extract_cox_group(cox_minimal, "Minimal"),
 extract_cox_group(cox_full, "Full")),
 file.path(OUTPUT_DIR, "cox_supplementary_results.csv"))
rates_30 <- model_data_30 |> group_by(group) |>
 summarise(events = sum(event_30), person_days = sum(time_30),
 person_years = person_days/365.25, .groups = "drop") |> rowwise() |>
 mutate(rate_per_person_year = events/person_years,
 lower_95_CI = poisson.test(events, T = person_years)$conf.int[1],
 upper_95_CI = poisson.test(events, T = person_years)$conf.int[2]) |> ungroup()
readr::write_csv(rates_30, file.path(OUTPUT_DIR, "incidence_rates_30_day.csv"))
all_models <- list(crude = cox_crude, primary = cox_adjusted, minimal = cox_minimal,
 full = cox_full, day30 = cox_adjusted_30,
 IPTW = cox_iptw, IPTW_adjusted = cox_iptw_doubly_adjusted)
all_coefficients <- bind_rows(lapply(names(all_models), function(nm) {
 su <- summary(all_models[[nm]])
 tibble(model = nm, term = rownames(su$coefficients),
 HR = su$conf.int[, "exp(coef)"], lower_95_CI = su$conf.int[, "lower .95"],
 upper_95_CI = su$conf.int[, "upper .95"], p_value = su$coefficients[, "Pr(>|z|)"])
}))
readr::write_csv(all_coefficients, file.path(OUTPUT_DIR, "cox_all_coefficients.csv"))
model_counts <- bind_rows(lapply(names(all_models), function(nm) {
 m <- all_models[[nm]]
 tibble(model = nm, participants = m$n, events = m$nevent,
 formula = paste(deparse(formula(m)), collapse = " "))
}))
readr::write_csv(model_counts, file.path(OUTPUT_DIR, "model_sample_sizes_and_formulas.csv"))
stopifnot(all(model_counts$participants == 118L),
 all(model_counts$events == ifelse(model_counts$model == "day30", 59L, 106L)))

# Appendix weighted Nelson–Aalen cumulative hazard: deliberately untransformed.
w_na <- survfit(Surv(time_days, event) ~ group, data = iptw_data,
 weights = weight, stype = 2, ctype = 1)
png(file.path(OUTPUT_DIR, "IPTW_weighted_Nelson_Aalen_cumulative_hazard.png"),
 width = 2700, height = 1800, res = 300)
plot(w_na, fun = "cumhaz", col = c("#F8766D", "#00BFC4"), conf.int = FALSE,
 xlim = c(0, 365), xlab = "Days from HCV antibody-positive date",
 ylab = "Cumulative hazard of DAA treatment initiation",
 main = "IPTW-weighted Nelson–Aalen cumulative-hazard curve")
legend("bottomright", c("Historical control", "RS-DAA intervention"),
 col = c("#F8766D", "#00BFC4"), lty = 1, bty = "n")
dev.off()
weight_distribution <- function(w, label) tibble(stage = label, n = length(w),
 min = min(w), p01 = quantile(w, .01), p05 = quantile(w, .05), median = median(w),
 mean = mean(w), p95 = quantile(w, .95), p99 = quantile(w, .99), max = max(w),
 effective_sample_size = sum(w)^2/sum(w^2))
readr::write_csv(bind_rows(weight_distribution(iptw_data$stabilised_weight, "Before truncation"),
 weight_distribution(iptw_data$weight, "After truncation")),
 file.path(OUTPUT_DIR, "IPTW_weight_distributions.csv"))
readr::write_csv(iptw_data |> group_by(group) |> summarise(n = n(),
 min_ps = min(propensity_score), p01_ps = quantile(propensity_score, .01),
 p05_ps = quantile(propensity_score, .05), median_ps = median(propensity_score),
 p95_ps = quantile(propensity_score, .95), p99_ps = quantile(propensity_score, .99),
 max_ps = max(propensity_score)), file.path(OUTPUT_DIR, "IPTW_propensity_score_distribution.csv"))

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
