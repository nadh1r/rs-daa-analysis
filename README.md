# rs-daa-analysis

Reproducible R scripts for the RS-DAA quasi-experimental study evaluating rapid-start direct-acting antiviral (DAA) treatment initiation among people who use drugs (PWUD) in Malaysia.

## Repository contents

| File | Purpose |
|---|---|
| Control Data Script Cleaned.R | Cleaning and analysis pipeline for historical control cohort |
| Intervention Data Script Cleaned.R | Cleaning and analysis pipeline for intervention cohort |
| Combined Control and Intervention Data Analysis Cleaned.R | Main combined survival analysis |
| Exploratory Multivariable Analysis Cleaned.R | Exploratory Cox regression analyses |
| IPTW Analysis Cleaned.R | IPTW sensitivity analysis |
| RMST Sensitivity Analysis Cleaned.R | RMST sensitivity analysis for non-proportional hazards |

## Main analyses

Primary analyses included:

- Kaplan–Meier survival analysis
- Log-rank testing
- Cox proportional hazards regression
- Stratified Cox regression
- IPTW sensitivity analysis
- RMST sensitivity analysis

## Software and packages

Analyses were conducted in R using packages including:

- survival
- survminer
- dplyr
- broom
- flextable
- officer
- survRM2

## Reproducibility

Raw participant-level data are not publicly shared due to confidentiality and ethics restrictions.

This repository provides reproducible analytical workflows and statistical scripts used in the study.

## Citation

If using or adapting these scripts, please cite the associated RS-DAA thesis and related publications.
