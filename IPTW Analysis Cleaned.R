## ============================================================
## IPTW SENSITIVITY ANALYSIS — RS-DAA VS HISTORICAL CONTROL
## GitHub-ready script: no local exports, no DOCX/PNG output
## Requires:
##   - cmb_tti_365
## ============================================================

library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(survival)
library(broom)
library(ggplot2)

stopifnot(exists("cmb_tti_365"))
stopifnot(all(c("id", "group", "time_days", "event") %in% names(cmb_tti_365)))

## ------------------------------------------------------------
## Configuration
## ------------------------------------------------------------

WEIGHT_TRUNC_PROBS <- c(0.01, 0.99)

iptw_vars <- c(
  "Age",
  "Employment",
  "Housing",
  "Marital status",
  "Current Methadone",
  "Hx of Incarceration",
  "Recent drug use",
  "Hx of injecting drug use"
)

iptw_vars <- intersect(iptw_vars, names(cmb_tti_365))

stopifnot(length(iptw_vars) > 0)

## ------------------------------------------------------------
## Helper functions
## ------------------------------------------------------------

lock_iptw_levels <- function(df) {
  df %>%
    mutate(
      group = factor(group, levels = c("Control", "Intervention")),
      treat = as.integer(group == "Intervention"),
      Age = suppressWarnings(as.numeric(Age)),
      
      Employment = factor(Employment, levels = c("Employed", "Not employed")),
      Housing = factor(Housing, levels = c("Stable", "Unstable")),
      `Marital status` = factor(`Marital status`, levels = c("Married", "Not married")),
      `Current Methadone` = factor(`Current Methadone`, levels = c("No", "Yes")),
      `Hx of Incarceration` = factor(`Hx of Incarceration`, levels = c("No", "Yes")),
      `Recent drug use` = factor(`Recent drug use`, levels = c("No", "Yes")),
      `Hx of injecting drug use` = factor(`Hx of injecting drug use`, levels = c("No", "Yes"))
    )
}

summarise_ps <- function(df) {
  df %>%
    group_by(group) %>%
    summarise(
      n = n(),
      min_ps = min(ps, na.rm = TRUE),
      p01_ps = quantile(ps, 0.01, na.rm = TRUE),
      p05_ps = quantile(ps, 0.05, na.rm = TRUE),
      median_ps = median(ps, na.rm = TRUE),
      p95_ps = quantile(ps, 0.95, na.rm = TRUE),
      p99_ps = quantile(ps, 0.99, na.rm = TRUE),
      max_ps = max(ps, na.rm = TRUE),
      .groups = "drop"
    )
}

summarise_weights <- function(w) {
  tibble(
    n = length(w),
    min_w = min(w, na.rm = TRUE),
    p01_w = quantile(w, 0.01, na.rm = TRUE),
    p05_w = quantile(w, 0.05, na.rm = TRUE),
    median_w = median(w, na.rm = TRUE),
    mean_w = mean(w, na.rm = TRUE),
    p95_w = quantile(w, 0.95, na.rm = TRUE),
    p99_w = quantile(w, 0.99, na.rm = TRUE),
    max_w = max(w, na.rm = TRUE)
  )
}

smd_cont <- function(x, treat, w = NULL) {
  if (is.null(w)) w <- rep(1, length(x))
  
  x1 <- x[treat == 1]
  x0 <- x[treat == 0]
  w1 <- w[treat == 1]
  w0 <- w[treat == 0]
  
  m1 <- weighted.mean(x1, w1, na.rm = TRUE)
  m0 <- weighted.mean(x0, w0, na.rm = TRUE)
  
  v1 <- weighted.mean((x1 - m1)^2, w1, na.rm = TRUE)
  v0 <- weighted.mean((x0 - m0)^2, w0, na.rm = TRUE)
  
  (m1 - m0) / sqrt((v1 + v0) / 2)
}

smd_bin <- function(x, treat, level, w = NULL) {
  if (is.null(w)) w <- rep(1, length(x))
  
  z <- as.numeric(x == level)
  
  p1 <- weighted.mean(z[treat == 1], w[treat == 1], na.rm = TRUE)
  p0 <- weighted.mean(z[treat == 0], w[treat == 0], na.rm = TRUE)
  
  p <- (p1 + p0) / 2
  
  if (p == 0 || p == 1) return(NA_real_)
  
  (p1 - p0) / sqrt(p * (1 - p))
}

create_balance_table <- function(df, vars, weight_var = "sw_trunc") {
  bind_rows(
    tibble(
      variable = "Age",
      level = NA_character_,
      smd_unweighted = smd_cont(df$Age, df$treat),
      smd_weighted = smd_cont(df$Age, df$treat, df[[weight_var]])
    ),
    bind_rows(lapply(setdiff(vars, "Age"), function(v) {
      lv <- levels(df[[v]])
      
      bind_rows(lapply(lv[-1], function(l) {
        tibble(
          variable = v,
          level = l,
          smd_unweighted = smd_bin(df[[v]], df$treat, l),
          smd_weighted = smd_bin(df[[v]], df$treat, l, df[[weight_var]])
        )
      }))
    }))
  ) %>%
    mutate(
      abs_smd_unweighted = abs(smd_unweighted),
      abs_smd_weighted = abs(smd_weighted),
      balanced_after_weighting = abs_smd_weighted < 0.10
    )
}

cox_extract <- function(fit, model_name) {
  sm <- summary(fit)
  
  tibble(
    model = model_name,
    variable = rownames(sm$coefficients),
    HR = sm$coefficients[, "exp(coef)"],
    ci_low = sm$conf.int[, "lower .95"],
    ci_high = sm$conf.int[, "upper .95"],
    p_value = sm$coefficients[, "Pr(>|z|)"]
  )
}

## ------------------------------------------------------------
## 1) Define IPTW analysis dataset
## ------------------------------------------------------------

iptw_dat <- cmb_tti_365 %>%
  select(
    id,
    group,
    time_days,
    event,
    all_of(iptw_vars)
  ) %>%
  lock_iptw_levels() %>%
  drop_na(time_days, event, group, all_of(iptw_vars))

iptw_sample_summary <- tibble(
  n = nrow(iptw_dat),
  events = sum(iptw_dat$event == 1L),
  censored = sum(iptw_dat$event == 0L)
)

iptw_event_table <- iptw_dat %>%
  count(group, event, name = "n") %>%
  group_by(group) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  ungroup()

## ------------------------------------------------------------
## 2) Propensity score model
## Probability of intervention-arm membership
## ------------------------------------------------------------

ps_model <- glm(
  treat ~
    Age +
    Employment +
    Housing +
    `Marital status` +
    `Current Methadone` +
    `Hx of Incarceration` +
    `Recent drug use` +
    `Hx of injecting drug use`,
  data = iptw_dat,
  family = binomial()
)

iptw_dat <- iptw_dat %>%
  mutate(
    ps = predict(ps_model, type = "response")
  )

ps_summary <- summarise_ps(iptw_dat)

ps_plot <- ggplot(
  iptw_dat,
  aes(x = ps, fill = group)
) +
  geom_histogram(alpha = 0.5, position = "identity", bins = 20) +
  labs(
    title = "Propensity score distribution by group",
    x = "Propensity score",
    y = "Count"
  )

## ------------------------------------------------------------
## 3) Stabilised IPTW weights
## ------------------------------------------------------------

p_treat <- mean(iptw_dat$treat == 1)

iptw_dat <- iptw_dat %>%
  mutate(
    sw = case_when(
      treat == 1L ~ p_treat / ps,
      treat == 0L ~ (1 - p_treat) / (1 - ps),
      TRUE ~ NA_real_
    )
  )

weight_summary <- summarise_weights(iptw_dat$sw)

w_bounds <- quantile(
  iptw_dat$sw,
  probs = WEIGHT_TRUNC_PROBS,
  na.rm = TRUE
)

iptw_dat <- iptw_dat %>%
  mutate(
    sw_trunc = pmin(pmax(sw, w_bounds[[1]]), w_bounds[[2]])
  )

weight_summary_truncated <- summarise_weights(iptw_dat$sw_trunc)

## ------------------------------------------------------------
## 4) Balance diagnostics
## ------------------------------------------------------------

balance_tbl <- create_balance_table(
  df = iptw_dat,
  vars = iptw_vars,
  weight_var = "sw_trunc"
)

## ------------------------------------------------------------
## 5) IPTW-weighted Cox models
## ------------------------------------------------------------

m_iptw <- coxph(
  Surv(time_days, event) ~ group,
  data = iptw_dat,
  weights = sw_trunc,
  robust = TRUE,
  ties = "efron"
)

m_iptw_strat <- coxph(
  Surv(time_days, event) ~ group + strata(`Recent drug use`),
  data = iptw_dat,
  weights = sw_trunc,
  robust = TRUE,
  ties = "efron"
)

m_iptw_dr <- coxph(
  Surv(time_days, event) ~
    group +
    Age +
    Employment +
    Housing +
    `Marital status` +
    `Current Methadone` +
    `Hx of Incarceration` +
    `Hx of injecting drug use` +
    strata(`Recent drug use`),
  data = iptw_dat,
  weights = sw_trunc,
  robust = TRUE,
  ties = "efron"
)

iptw_results <- bind_rows(
  cox_extract(m_iptw, "IPTW weighted Cox"),
  cox_extract(m_iptw_strat, "IPTW weighted Cox stratified by recent drug use"),
  cox_extract(m_iptw_dr, "Doubly adjusted IPTW weighted Cox")
) %>%
  filter(variable == "groupIntervention")

## ------------------------------------------------------------
## 6) Weighted Kaplan-Meier curve
## ------------------------------------------------------------

km_iptw <- survfit(
  Surv(time_days, event) ~ group,
  data = iptw_dat,
  weights = sw_trunc
)

## ------------------------------------------------------------
## Final output object
## ------------------------------------------------------------

iptw_outputs <- list(
  data = iptw_dat,
  
  variables = iptw_vars,
  truncation_probs = WEIGHT_TRUNC_PROBS,
  truncation_bounds = w_bounds,
  
  sample_summary = iptw_sample_summary,
  event_table = iptw_event_table,
  
  propensity_score = list(
    model = ps_model,
    summary = ps_summary,
    plot = ps_plot
  ),
  
  weights = list(
    before_truncation = weight_summary,
    after_truncation = weight_summary_truncated
  ),
  
  balance = balance_tbl,
  
  models = list(
    weighted_cox = m_iptw,
    weighted_cox_stratified = m_iptw_strat,
    doubly_adjusted_weighted_cox = m_iptw_dr
  ),
  
  results = iptw_results,
  
  weighted_km = km_iptw
)

## Minimal console output
iptw_outputs$sample_summary
iptw_outputs$propensity_score$summary
iptw_outputs$weights$before_truncation
iptw_outputs$weights$after_truncation
iptw_outputs$balance
iptw_outputs$results