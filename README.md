# RS-DAA Analysis

This repository contains the reproducible R analysis used in the main thesis evaluating rapid-start direct-acting antiviral treatment (RS-DAA) against the historical standard-of-care cohort.

## Main analysis script

`RSDAA_thesis_appendix_analysis.R` is the definitive analysis script for the thesis appendix. It replaces the earlier separate control, intervention, combined, exploratory, RMST and IPTW scripts.

The script performs:

- data cleaning and harmonisation for the prospective RS-DAA and historical standard-of-care cohorts;
- baseline descriptive analysis;
- HCV care-cascade analysis;
- Nelson–Aalen cumulative-incidence estimation;
- bootstrap confidence intervals for Nelson–Aalen-derived median initiation times;
- log-rank comparison;
- crude and adjusted Cox regression;
- proportional-hazards diagnostics using scaled Schoenfeld residuals;
- adjusted Cox analysis restricted to the first 30 days;
- unadjusted and covariate-adjusted restricted mean survival time analysis; and
- inverse probability of treatment weighting sensitivity analysis.

## Analysis populations

The HCV-antibody-positive descriptive population is used for baseline characteristics and the care cascade.

The primary modified intention-to-treat (mITT) cohort comprises participants who:

1. had a positive HCV antibody result;
2. had detectable HCV RNA;
3. had a valid HCV antibody-test date; and
4. remained eligible to initiate DAA through the participating clinic pathway.

Participants with cirrhosis requiring specialist referral are excluded from the primary mITT cohort. The script constructs this cohort once and includes checks that stop execution if a participant identified as having cirrhosis enters any time-to-event, 30-day, RMST or IPTW analysis.

## Missing data

Multiple imputation was not performed. Multivariable analyses use complete cases for the variables required by the relevant model. Missing HCV antibody-test dates and treatment-initiation dates are not imputed.

## Implementation-outcome data

Implementation outcomes were summarised descriptively using outcome-specific eligible populations. The source counts for penetration and component-specific fidelity were manually derived and verified through reconciliation of the RS-DAA masterlist, the Klinik Kesihatan Pandamaran (Panda) treatment record and laboratory records. The repository contains only the verified aggregate counts and does not contain identifiable source records.

Penetration was calculated as 84 enrolled participants among 92 eligible HCV antibody-positive PWUD who were offered RS-DAA (91.3%). The 94 individuals initially identified comprised 84 recorded in the RS-DAA masterlist, seven identified from the Panda treatment record and three from laboratory records. Two individuals did not meet the eligibility criteria, leaving 92 eligible individuals; eight of those eligible individuals followed the conventional standard-of-care pathway.

Component-specific fidelity was calculated using the applicable eligible population for each component:

- bundled pretreatment assessment: 84 of 84 enrolled participants;
- community blood collection: 73 of 84 enrolled participants;
- telemedicine consultation: 36 of 54 DAA initiators; and
- CHW-supported medication delivery: 33 of 36 telemedicine initiators.

No composite fidelity score was calculated because the components occurred at different pathway stages and did not share a common denominator. Seven-day feasibility and its exact binomial 95% confidence interval were calculated in R among the 58 participants eligible for clinic-based DAA initiation. The corresponding 14-day proportion was descriptive.

## Software requirements

The thesis analysis was conducted using R version 4.5.0. The script checks for the following packages:

- `readxl`
- `readr`
- `dplyr`
- `tidyr`
- `stringr`
- `lubridate`
- `survival`
- `survRM2`
- `ggplot2`
- `scales`
- `broom`

Install any missing packages before running the analysis:

```r
install.packages(c(
  "readxl", "readr", "dplyr", "tidyr", "stringr", "lubridate",
  "survival", "survRM2", "ggplot2", "scales", "broom"
))
```

## Running the analysis

1. Download or clone this repository.
2. Open `RSDAA_thesis_appendix_analysis.R` in RStudio.
3. At the beginning of the script, replace the two placeholder input paths with the locations of the control Excel workbook and intervention REDCap CSV export:

```r
CONTROL_FILE <- "path/to/control_masterlist.xlsx"
INTERVENTION_FILE <- "path/to/redcap_intervention_export.csv"
```

4. Run the script from the beginning in a clean R session.
5. Review the files produced in the `RSDAA_thesis_outputs` directory.

The script uses a fixed random seed and 2,000 participant-level bootstrap replicates so that the bootstrap confidence intervals can be reproduced.

## Reproducibility outputs

The output directory contains analysis tables, figures and supporting audit files, including:

- the cirrhosis-exclusion audit;
- the locked mITT cohort counts;
- baseline missingness;
- Nelson–Aalen estimates and bootstrap medians;
- Cox regression estimates and proportional-hazards diagnostics;
- 30-day Cox results;
- RMST results;
- implementation outcomes;
- IPTW balance and weighted estimates; and
- `sessionInfo.txt`, documenting the R and package versions used.

## Data confidentiality

The raw control and intervention datasets are not included in this repository. Do not upload identifiable participant data, names, record-level audit files or analysis outputs containing identifiers to a public repository. The `.gitignore` file should be configured to exclude raw data and generated output directories before any data are placed inside the local repository.

## Repository status

Earlier analysis scripts are superseded by `RSDAA_thesis_appendix_analysis.R`. Git retains the earlier versions in the repository history even after they are removed from the current file list.
