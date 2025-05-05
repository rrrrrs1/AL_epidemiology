# AL_epidemiology

## Acute Leukaemia Epidemiology Analysis in China

This repository contains R scripts for analysing acute leukaemia epidemiology data in China, with a focus on incidence, mortality, and survival outcomes across different subtypes and age groups.

## Project Overview

This analysis examines epidemiological patterns and survival outcomes of acute leukaemia in China, stratified by:
- Leukaemia subtypes (non-APL AML, APL, ALL, and other forms)
- Age groups (child vs. adult populations, with specific analyses for 0-14 and 0-19 age ranges)
- Time periods (2016-2018 vs. 2019-2020)

The project includes comprehensive statistical analyses of incidence rates, mortality rates, mortality-to-incidence ratios (MIR), and survival outcomes with various prognostic factors.

## Repository Structure

### Core Files
- `source.R`: Source file containing common functions, variable definitions, and data preprocessing steps used across all analysis scripts

### Figure Scripts
R scripts that generate visualisations:
- `figure1.R`: Age-specific incidence and mortality rates by acute leukaemia subtype with BMJ-style formatting
- `figure2.R`: Survival analysis plots including:
  - Overall survival by acute leukaemia subtype across different age groups
  - Survival analysis stratified by age within each leukaemia subtype
  - Treatment analysis (allo-HSCT and TKI therapy effects)

### Table Scripts
R scripts that generate statistical tables:
- `table1.R`: Incidence, mortality rates and MIR by acute leukaemia subtype and sex, with separate analyses for different paediatric age thresholds (0-14 and 0-19)
- `table2.R`: Overall survival rates by acute leukaemia subtype and age group with follow-up statistics
- `table3.R`: Cox regression analysis including:
  - Univariate and multivariate analysis of survival by time period
  - Proportional hazards testing
  - Analysis of clinical factors including treatment modalities

## Getting Started

### Prerequisites

This project requires R with the following packages:

#### Core data handling and visualization
```R
library(dplyr)
library(tidyr)
library(ggplot2)
library(openxlsx)
library(svglite)
library(RColorBrewer)
library(ggsci)
library(scales)
library(patchwork)
library(showtext)
```

#### Survival analysis
```R
library(survival)
library(survminer)
library(ggsurvfit)
```

#### Other packages
```R
library(sf)
library(epitools)
library(epiR)
library(reshape2)
library(lubridate)
library(MetBrewer)
library(mice)
```

### Data Requirements

The analysis uses several data files (not included in the repository):
- Leukemia incidence and mortality data by age, sex, and subtype
- Patient-level data for survival analysis
- Demographic and clinical data for Cox regression models

### Running the Analysis

1. Set your working directory in each script to match your local file structure
2. Ensure all required data files are in the correct locations
3. Run the scripts in the following order:
   - Table generation scripts (`table1.R`, `table2.R`, `table3.R`)
   - Figure generation scripts (`figure1.R`, `figure2.R`)

## Analysis Components

### Epidemiological Analysis
- Age-standardised incidence rates (ASIR) by acute leukaemia subtype, sex, and age group
- Age-standardised mortality rates (ASMR) by acute leukaemia subtype, sex, and age group
- Crude rates and estimated annual case counts for China
- Mortality-to-incidence ratios (MIR) with confidence intervals
- Proportion calculations of different acute leukaemia subtypes
- Comparative analyses using two different paediatric age thresholds (0-14 and 0-19)

### Survival Analysis
- Overall survival (OS) by acute leukaemia subtype (AML non-M3, APL, ALL, other acute leukaemias)
- Age-stratified survival analyses (paediatric vs. adult outcomes)
- Effect of treatment modalities:
  - Allogeneic haematopoietic stem cell transplantation (allo-HSCT)
  - Tyrosine kinase inhibitor (TKI) therapy in Ph+ ALL
  - Combined effects of transplant and TKI generation
- Temporal trends in survival comparing 2016-2018 vs. 2019-2020 periods

### Statistical Methods
- Kaplan-Meier survival analysis with log-rank tests and visualisation
- Restricted mean survival time analysis
- Univariate and multivariate Cox proportional hazards regression
- Proportional hazards assumption testing
- Multiple imputation for handling missing clinical data
- Categorisation of laboratory values based on clinical thresholds
- Bonferroni correction for multiple comparisons

## Output

The scripts generate:
- Publication-quality tables in Excel format:
  - Epidemiological tables with incidence, mortality, and MIR statistics
  - Survival rate tables with confidence intervals
  - Cox regression hazard ratio tables with significance testing
- BMJ-style figures in both PDF and SVG formats:
  - Age-specific incidence and mortality rate curves by acute leukaemia subtype and sex
  - Kaplan-Meier survival curves with risk tables and statistical annotations
  - Treatment effect visualisation with stratification by age and therapy type
- Statistical summaries and model diagnostic outputs

## Contributors

This project was developed by researchers from the China Acute Leukaemia Epidemiology, Clinical and Multi-omics Data Consortium.

## License

This project is licensed under the MIT License - see the LICENSE file for details.
