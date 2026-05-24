## ============================================================
## RMST SENSITIVITY ANALYSIS FOR NON-PROPORTIONAL HAZARDS
## GitHub-ready script: no local exports, no DOCX output
## Requires:
##   - cmb_tti_365
##   - FOLLOWUP_DAYS
## ============================================================

library(dplyr)
library(tibble)
library(survRM2)

stopifnot(exists("cmb_tti_365"))
stopifnot(exists("FOLLOWUP_DAYS"))
stopifnot(all(c("time_days", "event", "group") %in% names(cmb_tti_365)))

tti_complete <- cmb_tti_365 %>%
  mutate(
    group = factor(group, levels = c("Control", "Intervention")),
    group_binary = as.integer(group == "Intervention")
  )

rmst_res <- survRM2::rmst2(
  time = tti_complete$time_days,
  status = tti_complete$event,
  arm = tti_complete$group_binary,
  tau = FOLLOWUP_DAYS
)

rmst_control <- as.data.frame(rmst_res$RMST.arm0$result) %>%
  tibble::rownames_to_column("measure") %>%
  as_tibble() %>%
  mutate(group = "Control")

rmst_intervention <- as.data.frame(rmst_res$RMST.arm1$result) %>%
  tibble::rownames_to_column("measure") %>%
  as_tibble() %>%
  mutate(group = "Intervention")

rmst_by_arm <- bind_rows(
  rmst_control,
  rmst_intervention
) %>%
  select(group, everything())

rmst_contrast <- as.data.frame(rmst_res$unadjusted.result) %>%
  tibble::rownames_to_column("contrast") %>%
  as_tibble()

rmst_outputs <- list(
  tau = FOLLOWUP_DAYS,
  model = rmst_res,
  by_arm = rmst_by_arm,
  contrast = rmst_contrast
)

rmst_outputs$by_arm
rmst_outputs$contrast