# RS-DAA thesis: reproducible analytical workflow

`RSDAA_Reviewer_Analysis_v2.R` consolidates the analyses retained for the thesis comparison of rapid-start DAA treatment with historical standard of care. Exploratory analyses that were not retained are not part of this workflow. Raw participant data are supplied separately under the applicable access arrangements.

## Inputs and execution

Three local input files are required:

1. Historical control Excel workbook: worksheet 1 provides baseline and treatment data; worksheet 2 provides the reconciled SVR fields, matched to the same participants using the source IC field. Baseline covariates are not taken from worksheet 2.
2. Intervention REDCap export, in Excel or CSV format, containing baseline and treatment-status events, including completion and SVR dates/status.
3. Implementation Excel masterlist, worksheet `KK PANDAMARAN`, containing `Redcap No`, `Blood taken where`, `Mysejahtera/Whatsapp VC/Physical`, and `Meds`.

Set `CONTROL_FILE`, `INTERVENTION_FILE`, `IMPLEMENTATION_FILE` and `OUTPUT_DIR` at the beginning of the script, or supply the environment variables `RSDAA_CONTROL_FILE`, `RSDAA_INTERVENTION_FILE`, `RSDAA_IMPLEMENTATION_FILE` and `RSDAA_OUTPUT_DIR`. Run the entire script in a clean R session. Do not run isolated sections against objects left by an earlier analysis.

Install the required packages first:

```r
install.packages(c("readxl", "readr", "dplyr", "tidyr", "stringr", "lubridate",
                   "survival", "survRM2", "ggplot2", "scales", "survminer"))
```

The thesis reported R 4.5.0. Each execution writes `sessionInfo.txt` with the versions actually used.

## Populations and descriptive analyses

The antibody-positive descriptive population comprises 168 participants: 84 historical SOC and 84 RS-DAA. It is used for the broader baseline table and care cascade. A second baseline table describes the final clinic-eligible, RNA-detected analysis population.

The primary mITT cohort comprises 118 participants: 60 historical SOC and 58 RS-DAA. Twelve historical and eight RS-DAA participants with cirrhosis requiring specialist referral are excluded from this cohort. All time-to-event analyses use this cohort. Checks verify 52 historical and 54 RS-DAA initiation events within 365 days, and verify each fitted Cox model's actual participant and event counts. The 30-day analysis has 59 events. These checks validate calculated results; they do not set model estimates.

Baseline summaries include age mean/SD and categorical frequencies/percentages. Age comparisons use Wilcoxon rank-sum tests. The broader table uses chi-squared or Fisher exact tests according to expected cell counts; the 118-person table uses Fisher exact tests for categorical variables, matching the parallel analysis. Baseline comparisons are descriptive and do not select adjustment covariates. There is no multiple imputation. Complete-case model counts are exported explicitly.

## Care cascade and SVR classification

The cascade reports antibody positivity, RNA testing, RNA detection, clinic eligibility, DAA initiation, completion, SVR assessment and achievement with stage-specific denominators. Descriptive initiation includes any recorded initiation, including historical initiations after 365 days; it is distinct from the primary 365-day event definition.

The script applies the two source-date corrections confirmed during reconciliation, without modifying the workbook: historical source No 12 with initiation on 5 August 2023 has an SVR test date of 30 January 2024; source No 29 with initiation on 3 May 2024 has an SVR test date of 17 October 2024. Historical source numbers repeat between years, so cohort identity is not established by No alone.

Tests beyond 168 days after completion are excluded: 16 historical SOC and two RS-DAA. Four historical tests at 83 days and nine RS-DAA tests at 75–81 days are retained under the study's operational classification and flagged separately. These early tests must not be described as satisfying a strict minimum of 12 weeks after completion. Results without an eligible assessment are not classified as virological failure. The expected operational cascade is 20 assessed/19 achieved historically and 34 assessed/31 achieved in RS-DAA; discrepancies trigger review rather than replacement with target counts.

## Nelson–Aalen analysis and primary figure

Cumulative hazard is estimated by summing event counts divided by the number at risk immediately before each event time. The main figure transforms that hazard to cumulative incidence using one minus the exponential of negative cumulative hazard. This is the transformed Nelson–Aalen estimate, rather than the complement of Kaplan–Meier.

The main figure includes pointwise 95% confidence bands, the risk/cumulative-event table and the full 365-day horizon. It uses `survfit(stype = 2, ctype = 1, conf.type = "log")`, followed by the event transformation. These confidence limits correspond to normal limits on the Nelson–Aalen cumulative hazard transformed to the incidence scale, with the lower hazard bound constrained to zero. They are not bootstrap intervals or simultaneous confidence bands.

Median initiation time is the first event time at which transformed incidence reaches 0.50. Percentile confidence intervals use 5,000 participant-level resamples with replacement within each study group, retaining group sample sizes, seed 20260817 and quantile type 6. The output reports valid medians, excluded repetitions and the percentage reaching the median. Missing bootstrap medians are disclosed rather than silently omitted. The log-rank test compares the overall initiation-time distributions.

## Regression and sensitivity analyses

The primary adjusted Cox model includes study group, age, employment, housing, marital status, HIV co-infection, current methadone, incarceration history, recent drug use and injecting history. Efron handling of ties is used. Scaled Schoenfeld-residual tests and plots assess proportional hazards. The 365-day hazard ratios are interpreted in light of the observed non-proportional hazards and do not establish a constant effect over follow-up.

The same adjusted specification is fitted after administrative censoring at day 30. RMST analyses through day 365 include unadjusted and covariate-adjusted estimates of time remaining without initiation. Initiation rates and exact Poisson confidence limits are reported per person-year.
Supplementary Cox models are explicitly defined:

- Minimal: group, age, housing, current methadone, incarceration history and injecting history.
- Full: the primary specification plus race, education, ATS, benzodiazepines, cannabis and other drug use.

Actual sample sizes, event counts and formulas are exported for every Cox model. Historical appendix statements referring to a 138-person cohort should not be used to describe these 118-person fits.

IPTW uses the nine primary adjustment variables in a logistic propensity-score model. Stabilised weights are truncated at the pooled 1st and 99th percentiles. Outputs include propensity-score distributions, overlap histogram, weights before/after truncation, effective sample sizes, covariate balance, and robust weighted Cox estimates an additionally adjusted weighted Cox model. Incomplete overlap and any residual covariate imbalance should be acknowledged. The appendix weighted Nelson–Aalen figure shows untransformed cumulative hazard, distinct from the main incidence figure.

## Implementation outcomes

Penetration uses externally reconciled aggregate counts: 94 identified, two ineligible, 92 offered, 84 enrolled and eight following conventional SOC. Those aggregates are not inferred from the 84-person implementation worksheet. Bundled pretreatment assessment, 84/84, is a user-verified aggregate because that component is not a separate worksheet field.

Other fidelity measures are calculated from the implementation worksheet after excluding REDCap numbers 3 and 73 (invalid participants) and example row 99. MySejahtera is grouped with WhatsApp as telemedicine. The source yields community blood collection 73/84 and telemedicine 36/54. Among the 36 telemedicine participants, direct CHW delivery is 32, postal delivery one, and self-pickup three. Postal delivery is included within CHW-supported delivery under the study classification, giving 33/36 CHW-supported deliveries; the source-mode breakdown remains in the audit. The script checks the participant linkage for the 54 consultation records against clinical initiation records. No composite fidelity score is calculated.

Seven-day feasibility and its exact binomial confidence interval are calculated among the 58 clinic-eligible RS-DAA participants; 14- and 30-day proportions are also exported.

## Outputs and confidentiality

Outputs include both baseline tables, the full cascade, SVR timing audits, cohort checks, the main transformed figure, median bootstrap results, log-rank results, primary/supplementary/30-day models and diagnostics, rate and RMST results, implementation summaries and IPTW results/figure.
