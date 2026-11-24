# UK Biobank breast cancer PRS analysis code
Code for: Integrating Polygenic Risk Scores and Clinical-Reproductive Factors
Improves Breast Cancer Risk Stratification (Breast Cancer Research submission).

## Data
UK Biobank individual-level data are NOT included. Access via UK Biobank
application (https://www.ukbiobank.ac.uk); this study used Application 362917.
Fields used are listed in the manuscript Methods.

## Software
R 4.4.3; packages: survival 3.8-3, rstpm2 1.7.1, data.table 1.18.2.1,
dplyr 1.1.4, readr 2.1.6, gtsummary 2.5.0, msm, car, ggplot2, scales.
Full environment: sessionInfo.txt.

## Run order
1. main_analysis_final.R   # HRs, RERI, NRI/IDI/DCA, calibration,
                           # reference-profile absolute risk, marginal risks
2. sensitivity_analysis.R  # menarche-origin flexible parametric model
3. table1.R                # Table 1


## Notes
- Bootstrap loops use nboot = 1000 in the final version.
