# AL_epidemiology

## Leukemia Epidemiology Analysis in China

This repository contains R scripts for analyzing leukemia epidemiology data in China, with a focus on incidence, mortality, and survival outcomes across different subtypes and age groups.

## Project Overview

This analysis examines epidemiological patterns and survival outcomes of leukemia in China, stratified by:
- Leukemia subtypes (AML, APL, ALL, and other forms)
- Age groups (child vs. adult populations)
- Time periods (2016-2018 vs. 2019-2020)

The project includes comprehensive statistical analyses of incidence rates, mortality rates, mortality-to-incidence ratios (MIR), and survival outcomes with various prognostic factors.

## Repository Structure

### Figure Scripts
R scripts that generate visualizations:
- `figure1.R`: Age-specific incidence and mortality rates by leukemia subtype
- `figure2.R`: Survival analysis by leukemia subtype, age group, and treatment

### Table Scripts
R scripts that generate statistical tables:
- `table1.R`: Incidence, mortality rates and MIR by leukemia subtype and sex
- `table2.R`: Overall survival rates by leukemia subtype and age group
- `table3.R`: Cox regression analysis of survival by time period and clinical factors

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
- Age-standardized incidence and mortality rates
- Crude rates and case counts by leukemia subtype
- Mortality-to-incidence ratios
- Proportion of different leukemia subtypes

### Survival Analysis
- Overall survival (OS) by leukemia subtype and age group
- Effect of treatment modalities (allo-HSCT, TKI therapy)
- Temporal trends in survival (2016-2018 vs. 2019-2020)

### Statistical Methods
- Kaplan-Meier survival analysis with log-rank tests
- Univariate and multivariate Cox proportional hazards regression
- Multiple imputation for handling missing data

## Output

The scripts generate:
- Publication-quality tables (Excel format)
- BMJ-style figures (PDF/SVG format)
- Statistical summaries and model outputs

## Contributors

This project was developed by researchers from the China Acute Leukaemia Epidemiology, Clinical and Multi-omics Data Consortium.

## License

This project is licensed under the MIT License - see the LICENSE file for details.
