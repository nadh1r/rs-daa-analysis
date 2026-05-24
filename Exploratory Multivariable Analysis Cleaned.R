## ============================================================
## OBJECTIVE 3 — EXPLORATORY ASSOCIATIONS WITH TTI
## GitHub-ready script: no local exports, no DOCX output
## Input:
##   - cmb_tti_365
## ============================================================

library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(rlang)
library(survival)
library(broom)

stopifnot(exists("cmb_tti_365"))
stopifnot(all(c("time_days", "event", "group") %in% names(cmb_tti_365)))

## ------------------------------------------------------------
## Configuration
## ------------------------------------------------------------

SCREENING_P_THRESHOLD <- 0.20

candidate_vars <- c(
  "Age",
  "Race",
  "Education",
  "Employment",
  "Marital status",
  "Housing",
  "HIV co-infection",
  "Hx of Incarceration",
  "Current Methadone",
  "ATS",
  "Benzos",
  "Cannabis",
  "Others",
  "Recent drug use",
  "Hx of injecting drug use"
)

candidate_vars <- intersect(candidate_vars, names(cmb_tti_365))

if (length(candidate_vars) == 0) {
  stop("No candidate variables found in cmb_tti_365.")
}

## ------------------------------------------------------------
## Helper functions
## ------------------------------------------------------------

backtick_var <- function(x) {
  paste0("`", x, "`")
}

lock_objective3_levels <- function(df) {
  df %>%
    mutate(
      group = factor(group, levels = c("Control", "Intervention")),
      
      event = case_when(
        is.numeric(event) ~ as.integer(event),
        is.logical(event) ~ as.integer(event),
        as.character(event) %in% c("1", "Yes", "Started", "Event") ~ 1L,
        as.character(event) %in% c("0", "No", "Not started", "No event") ~ 0L,
        TRUE ~ NA_integer_
      ),
      
      Age = suppressWarnings(as.numeric(Age)),
      
      Race = factor(Race, levels = c("Malay", "Non-Malay")),
      Education = factor(Education, levels = c("Primary or less", "Secondary or more")),
      Employment = factor(Employment, levels = c("Employed", "Not employed")),
      `Marital status` = factor(`Marital status`, levels = c("Married", "Not married")),
      Housing = factor(Housing, levels = c("Stable", "Unstable")),
      
      `HIV co-infection` = factor(`HIV co-infection`, levels = c("No", "Yes")),
      `Hx of Incarceration` = factor(`Hx of Incarceration`, levels = c("No", "Yes")),
      `Current Methadone` = factor(`Current Methadone`, levels = c("No", "Yes")),
      ATS = factor(ATS, levels = c("No", "Yes")),
      Benzos = factor(Benzos, levels = c("No", "Yes")),
      Cannabis = factor(Cannabis, levels = c("No", "Yes")),
      Others = factor(Others, levels = c("No", "Yes")),
      `Recent drug use` = factor(`Recent drug use`, levels = c("No", "Yes")),
      `Hx of injecting drug use` = factor(`Hx of injecting drug use`, levels = c("No", "Yes"))
    )
}

get_lrt_p <- function(anova_df) {
  p_col <- grep(
    "P\\(>\\|Chi\\|\\)|Pr\\(>\\|Chi\\|\\)",
    names(anova_df),
    value = TRUE
  )
  
  if (length(p_col) == 0) {
    return(NA_real_)
  }
  
  suppressWarnings(as.numeric(anova_df[[p_col[1]]][nrow(anova_df)]))
}

cox_to_tbl <- function(fit) {
  broom::tidy(
    fit,
    exponentiate = TRUE,
    conf.int = TRUE
  ) %>%
    transmute(
      term,
      HR = estimate,
      ci_low = conf.low,
      ci_high = conf.high,
      p_value = p.value
    )
}

fit_univariable_screen <- function(df, variable) {
  dat_cc <- df %>%
    select(time_days, event, group, all_of(variable)) %>%
    drop_na()
  
  if (nrow(dat_cc) < 10) {
    return(list(
      screen = tibble(
        variable = variable,
        p_lrt = NA_real_,
        n_used = nrow(dat_cc),
        events_used = sum(dat_cc$event == 1L),
        status = "too_few_rows"
      ),
      terms = tibble()
    ))
  }
  
  if (length(unique(dat_cc[[variable]])) < 2) {
    return(list(
      screen = tibble(
        variable = variable,
        p_lrt = NA_real_,
        n_used = nrow(dat_cc),
        events_used = sum(dat_cc$event == 1L),
        status = "no_variation"
      ),
      terms = tibble()
    ))
  }
  
  f_base <- Surv(time_days, event) ~ group
  f_full <- as.formula(
    paste0("Surv(time_days, event) ~ group + ", backtick_var(variable))
  )
  
  m_base <- coxph(f_base, data = dat_cc, ties = "efron")
  m_full <- coxph(f_full, data = dat_cc, ties = "efron")
  
  lrt <- anova(m_base, m_full, test = "LRT")
  p_lrt <- get_lrt_p(as.data.frame(lrt))
  
  terms <- cox_to_tbl(m_full) %>%
    filter(term != "groupIntervention") %>%
    mutate(variable = variable) %>%
    select(variable, term, HR, ci_low, ci_high, p_value)
  
  list(
    screen = tibble(
      variable = variable,
      p_lrt = p_lrt,
      n_used = m_full$n,
      events_used = m_full$nevent,
      status = "ok"
    ),
    terms = terms
  )
}

is_categorical_var <- function(x) {
  is.factor(x) || is.character(x) || is.logical(x)
}

separation_check <- function(df, variable) {
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

check_shortlisted_separation <- function(df, variables) {
  map_dfr(variables, ~ separation_check(df, .x))
}

create_event_support_table <- function(df, variables) {
  map_dfr(variables, function(variable) {
    vq <- sym(variable)
    
    df %>%
      filter(!is.na(!!vq), !is.na(event)) %>%
      count(
        variable = variable,
        level = !!vq,
        event,
        name = "n"
      ) %>%
      mutate(
        event_status = if_else(event == 1L, "Event", "No event")
      ) %>%
      select(variable, level, event_status, n)
  })
}

build_mv_formula <- function(variables) {
  rhs <- if (length(variables) == 0) {
    "group"
  } else {
    paste(c("group", backtick_var(variables)), collapse = " + ")
  }
  
  as.formula(paste("Surv(time_days, event) ~", rhs))
}

fit_exploratory_mv <- function(df, variables) {
  formula <- build_mv_formula(variables)
  
  dat_mv <- df %>%
    select(time_days, event, group, all_of(variables)) %>%
    drop_na()
  
  fit <- coxph(
    formula,
    data = dat_mv,
    ties = "efron",
    na.action = na.exclude
  )
  
  list(
    formula = formula,
    data = dat_mv,
    model = fit,
    table = cox_to_tbl(fit),
    model_n = tibble(
      n_total_risk_set = nrow(df),
      n_used_complete_case = nrow(dat_mv),
      n_dropped_missing = nrow(df) - nrow(dat_mv),
      events_used = sum(dat_mv$event == 1L),
      censored_used = sum(dat_mv$event == 0L)
    )
  )
}

zph_to_tbl <- function(zph_obj) {
  as.data.frame(zph_obj$table) %>%
    rownames_to_column("term") %>%
    as_tibble() %>%
    rename(p_value = p)
}

## ------------------------------------------------------------
## 0) Standardise input
## ------------------------------------------------------------

dat0 <- cmb_tti_365 %>%
  lock_objective3_levels()

stopifnot(all(na.omit(unique(dat0$event)) %in% c(0L, 1L)))

## ------------------------------------------------------------
## 1) Univariable screening using global LRT
## Each model: Surv(time_days, event) ~ group + X
## Compared against: Surv(time_days, event) ~ group
## ------------------------------------------------------------

uni_results <- map(
  candidate_vars,
  safely(~ fit_univariable_screen(dat0, .x))
)

uni_fail <- map2_dfr(candidate_vars, uni_results, function(variable, result) {
  tibble(
    variable = variable,
    ok = is.null(result$error),
    error = if (!is.null(result$error)) conditionMessage(result$error) else NA_character_
  )
}) %>%
  filter(!ok)

uni_ok <- uni_results %>%
  map("result") %>%
  compact()

uni_screen <- bind_rows(map(uni_ok, "screen")) %>%
  arrange(p_lrt)

uni_terms <- bind_rows(map(uni_ok, "terms")) %>%
  arrange(p_value)

shortlisted_vars <- uni_screen %>%
  filter(!is.na(p_lrt), p_lrt < SCREENING_P_THRESHOLD) %>%
  pull(variable)

## ------------------------------------------------------------
## 2) Feasibility check among shortlisted categorical variables
## ------------------------------------------------------------

shortlisted_cat_vars <- shortlisted_vars[
  map_lgl(shortlisted_vars, ~ is_categorical_var(dat0[[.x]]))
]

event_support_tbl <- create_event_support_table(
  dat0,
  shortlisted_cat_vars
)

sep_tbl <- check_shortlisted_separation(
  dat0,
  shortlisted_cat_vars
)

drop_vars_sep <- sep_tbl %>%
  group_by(variable) %>%
  summarise(any_sep = any(sep_flag), .groups = "drop") %>%
  filter(any_sep) %>%
  pull(variable)

## ------------------------------------------------------------
## 3) Exploratory multivariable Cox model
## Includes group + shortlisted variables without separation
## ------------------------------------------------------------

mv_vars <- setdiff(shortlisted_vars, drop_vars_sep)

mv_results <- fit_exploratory_mv(
  df = dat0,
  variables = mv_vars
)

cox_mv <- mv_results$model
mv_tbl <- mv_results$table
mv_model_n <- mv_results$model_n

## ------------------------------------------------------------
## 4) PH check
## ------------------------------------------------------------

ph_mv <- cox.zph(cox_mv)
ph_mv_tbl <- zph_to_tbl(ph_mv)

## ------------------------------------------------------------
## Final output object
## ------------------------------------------------------------

objective3_outputs <- list(
  input_data = dat0,
  
  candidate_variables = candidate_vars,
  screening_threshold = SCREENING_P_THRESHOLD,
  
  univariable = list(
    failures = uni_fail,
    screen = uni_screen,
    term_estimates = uni_terms,
    shortlisted_variables = shortlisted_vars
  ),
  
  feasibility = list(
    shortlisted_categorical_variables = shortlisted_cat_vars,
    event_support = event_support_tbl,
    separation = sep_tbl,
    variables_excluded_for_separation = drop_vars_sep
  ),
  
  multivariable = list(
    variables_entered = mv_vars,
    formula = mv_results$formula,
    model_data = mv_results$data,
    model = cox_mv,
    model_n = mv_model_n,
    table = mv_tbl
  ),
  
  proportional_hazards = list(
    test = ph_mv,
    table = ph_mv_tbl
  )
)

## Minimal console output
objective3_outputs$univariable$screen
objective3_outputs$feasibility$variables_excluded_for_separation
objective3_outputs$multivariable$model_n
objective3_outputs$multivariable$table
objective3_outputs$proportional_hazards$table