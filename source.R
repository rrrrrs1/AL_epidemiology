# BMJ Survival Analysis Code: Leukaemia Subtypes Survival Analysis
# Repository: 
# ------------------------------------------------------------------------------

# Software Information
# 
# R version 4.4.2 (2024-10-31 ucrt)
# Platform: x86_64-w64-mingw32/x64 running under Windows Server 2022 x64 (build 20348)
# Locale: Chinese (Simplified)_China.utf8

# 
# Core packages: tidyverse 2.0.0, dplyr 1.1.4, readr 2.1.5, tidyr 1.3.1, stringr 1.5.1, purrr 1.0.2
# Survival analysis: survival 3.7-0, survminer 0.5.0, ggsurvfit 1.1.0, survRM2 1.0-4
# Data handling: openxlsx 4.2.7.1, broom 1.0.7, readxl 1.4.3
# Visualization: ggpubr 0.6.0, ggrepel 0.9.6, RColorBrewer 1.1-3, scales 1.3.0, ggplot2 3.5.1
# Additional packages: ggh4x 0.2.8, magrittr 2.0.3, ggpp 0.5.8-1, patchwork 1.3.0, MetBrewer 0.2.0

setwd("C:/Data/rms/data/2412_YW")
rm(list=ls())

# ------------------------------------------------------------------------------
# 1. Load required packages
# ------------------------------------------------------------------------------

# Core data manipulation and visualization
library(broom)        # Convert statistical results to tidy data frames
library(survival)     # Core survival analysis functions
library(survminer)    # Enhanced survival plots
library(tableone)     # Create "Table 1" summary statistics
library(tidyverse)    # Data manipulation and visualization
library(ggpubr)       # Publication-ready plots
library(ggrepel)      # Text labels that avoid overlapping
library(stringr)      # String manipulation
library(plotrix)      # Additional plotting functions
library(RColorBrewer) # Colour palettes
library(openxlsx)     # Excel file handling
library(scales)       # Scale functions for visualisation
library(gridExtra)    # Arrange multiple plots
library(glue)         # String interpolation
library(purrr)        # Functional programming tools
library(epitools)     # Epidemiological tools
library(MetBrewer)    # Art-inspired colour palettes
library(ggsurvfit)    # More survival plotting options
library(patchwork)    # Combine multiple plots
library(ggpp)         # Position geoms for ggplot2
library(magrittr)     # Pipe operators
library(ggh4x)        # Extended ggplot2 functions
library(survRM2)      # Restricted mean survival time
library(gtsummary)    # Summary tables for regression models

# Set fonts
windowsFonts(serif='Times New Roman', sans='Arial', mono='Courier New')

# ------------------------------------------------------------------------------
# 2. Load data files
# ------------------------------------------------------------------------------

# Load case and treatment data
load(file="241101case4.Rdata")         # 2016-2020 patient data
load(file="241101treat.Rdata")         # Treatment information
load(file="241101treat_all.Rdata")     # Complete treatment data
load(file="241101hqms_time1.Rdata")    # Hospital quality metrics
load(file="241101gene.Rdata")          # Genetic data
load(file="241101key.Rdata")           # Key variables
load(file="241101AML_case_MICM1.Rdata") # AML-specific case data
load(file="241101ALL_case_MICM1.Rdata") # ALL-specific case data

# ------------------------------------------------------------------------------
# 3. Create age categories
# ------------------------------------------------------------------------------

# Add 60-year age threshold category
case <- mutate(case, age60 = if_else(diag_age >= 60, "60周岁及以上", "60周岁以下"))
AML_select <- mutate(AML_select, age60 = if_else(diag_age >= 60, "60周岁及以上", "60周岁以下"))
ALL_select <- mutate(ALL_select, age60 = if_else(diag_age >= 60, "60周岁及以上", "60周岁以下"))

# ------------------------------------------------------------------------------
# 4. Process induction therapy data
# ------------------------------------------------------------------------------

# Extract induction therapy data
quancheng <- treat_all %>% 
  filter(event_type == "诱导治疗") %>%  
  group_by(ID) %>%
  summarise(
    event_name_new1 = paste0(event_name, collapse = " @ "),
    event_content_new1 = paste0(event_content, collapse = " @ "),
    event_remarks_new1 = paste0(event_remarks, collapse = " @ ")
  )

# Process TKI therapy data
tki <- quancheng %>% 
  filter(
    str_detect(event_name_new1, "TKI") |
      str_detect(event_content_new1, "(?i)伊马替尼|达沙替尼|氟马替尼|尼洛替尼|奥雷巴替尼|泊那替尼|普纳替尼|格列卫|TKI|格尼可|豪森昕福|依尼舒|达希纳|耐立克") |
      str_detect(event_remarks_new1, "(?i)伊马替尼|达沙替尼|氟马替尼|尼洛替尼|奥雷巴替尼|泊那替尼|普纳替尼|格列卫|TKI|格尼可|豪森昕福|依尼舒|达希纳|耐立克")
  ) %>% 
  mutate(check_all = "TKI") %>% 
  mutate(TKI_type = case_when(
    grepl("奥雷巴替尼|耐立克", event_content_new1) ~ "Olverematinib",
    grepl("奥雷巴替尼|耐立克", event_remarks_new1) ~ "Olverematinib",
    grepl("泊那替尼|普纳替尼", event_content_new1) ~ "Ponatinib",
    grepl("泊那替尼|普纳替尼", event_remarks_new1) ~ "Ponatinib",
    grepl("氟马替尼|豪森昕福", event_content_new1) ~ "Flumatinib",
    grepl("氟马替尼|豪森昕福", event_remarks_new1) ~ "Flumatinib",
    grepl("尼洛替尼|达希纳", event_content_new1) ~ "Nilotinib",
    grepl("尼洛替尼|达希纳", event_remarks_new1) ~ "Nilotinib",
    grepl("达沙替尼|依尼舒", event_content_new1) ~ "Dasatinib",
    grepl("达沙替尼|依尼舒", event_remarks_new1) ~ "Dasatinib",
    grepl("格列卫|格尼可|伊马替尼", event_content_new1) ~ "Imatinib",
    grepl("格列卫|格尼可|伊马替尼", event_remarks_new1) ~ "Imatinib",
    TRUE ~ "Other TKI"
  ))

# Categorize TKI by generation
tki1 <- quancheng %>% 
  filter(
    str_detect(event_name_new1, "TKI") |
      str_detect(event_content_new1, "(?i)伊马替尼|达沙替尼|氟马替尼|尼洛替尼|奥雷巴替尼|泊那替尼|普纳替尼|格列卫|TKI|格尼可|豪森昕福|依尼舒|达希纳|耐立克") |
      str_detect(event_remarks_new1, "(?i)伊马替尼|达沙替尼|氟马替尼|尼洛替尼|奥雷巴替尼|泊那替尼|普纳替尼|格列卫|TKI|格尼可|豪森昕福|依尼舒|达希纳|耐立克")
  ) %>% 
  mutate(check_all = "TKI") %>% 
  mutate(TKI_type_general = case_when(
    grepl("(?i)奥雷巴替尼|泊那替尼|普纳替尼|耐立克", event_content_new1) ~ "Third-generation TKI",
    grepl("(?i)奥雷巴替尼|泊那替尼|普纳替尼|耐立克", event_remarks_new1) ~ "Third-generation TKI",
    grepl("(?i)达沙替尼|氟马替尼|尼洛替尼|豪森昕福|依尼舒|达希纳", event_content_new1) ~ "Second-generation TKI",
    grepl("(?i)达沙替尼|氟马替尼|尼洛替尼|豪森昕福|依尼舒|达希纳", event_remarks_new1) ~ "Second-generation TKI",
    grepl("格列卫|格尼可|伊马替尼", event_content_new1) ~ "First-generation TKI",
    grepl("格列卫|格尼可|伊马替尼", event_remarks_new1) ~ "First-generation TKI",
    TRUE ~ "Unknown TKI"
  ))

# Join TKI data to ALL_select dataset
ALL_select <- ALL_select %>% 
  left_join(tki[, c("ID", "TKI_type")], by = "ID") %>%
  left_join(tki1[, c("ID", "TKI_type_general")], by = "ID")

# ------------------------------------------------------------------------------
# 5. Define colour palettes
# ------------------------------------------------------------------------------

# Define colour palettes for visualisation
mycol3 <- brewer.pal(11, "RdYlBu")
cols <- c(brewer.pal(6, "Set3"), "grey") 
cols_two <- c("#67B6C3", "#FFAA91")
cols1 <- c(brewer.pal(12, "Set3"), "grey")
cols2 <- c(cols_two, rev(cols1))
col_surv <- c("#BC3C29FF", "#0072B5FF", "#E18727FF", "#20854EFF", "#7876B1FF", 
              "#6F99ADFF", "#EE4C97FF", "#FB9A99", "#53647F", "#aec08b", 
              "#767676FF", "#AB9FD1", "#370335FF", "#FDBF6F", "#842F91", 
              "#FF7F00", "#CAB2D6", "#ED402B", "#F48326", "#FD6AB0", 
              "#FFD711", "#256B71", "#2DC9BA", "#7570B3")
col_paper <- c('#ADB7D6', '#C5DEF4', '#BFD6D9', '#E9ADAF', '#A2BBD0', '#8299CE', '#E3CBCD', '#DBDCDD')

# Load additional colour palette from MetBrewer
paper_col <- rev(met.brewer('VanGogh2')) 
paper_col1 <- rgb(t(col2rgb(paper_col)/255), alpha = 0.8)

# ------------------------------------------------------------------------------
# 6. Define variable lists and factor levels
# ------------------------------------------------------------------------------

# AML FAB classification
vars_AML_fab <- c("AML/M3", "AML/M0", "AML/M4", "AML/M6", "AML/M5", "髓细胞肉瘤", 
                  "AML/M1", "AML/M2", "AML/未分型", "AML/M7")  

# ALL immunophenotype classification
vars_ALL_mianyi <- c("B-ALL/pro-B", "B-ALL/common-B", "B-ALL/mature-B", "B-ALL/pre-B", 
                     "B-ALL/未分型", "T-ALL", "ALL/未分型")

# T-ALL molecular markers
vars_tall_mol <- c("del(1p)/SIL::TAL1", "t(10;11)/CALM::AF10", "HOX11 expression", 
                   "HOX11L2 expression", "NOTCH1 mutations", "SET-NUP214", "BCR::ABL1 fusion")

# WHO 2008/2016 AML classification
Vars_WHO2_AML <- c("PML::RARA fusion", "RUNX1::RUNX1T1 fusion", "CBFB::MYH11 fusion", 
                   "DEK::NUP214 fusion", "RBM15::MRTFA fusion", "BCR::ABL1 fusion",
                   "KMT2A rearrangement", "MECOM rearrangement", "NUP98 rearrangement",
                   "NPM1 mutation", "biCEBPA mutation", "myelodysplasia-related(AML-MR)",
                   "other defined genetic alterations", "AML,defined by differentiation")

# WHO 2016 AML classification
Vars_WHO2016_AML <- c("PML::RARA fusion", "RUNX1::RUNX1T1 fusion", "CBFB::MYH11 fusion", 
                      "DEK::NUP214 fusion", "RBM15::MRTFA fusion", "BCR::ABL1 fusion",
                      "MLLT3::KMT2A fusion", "GATA2,MECOM", "NPM1 mutation", 
                      "biCEBPA mutation", "RUNX1 mutation",
                      "myelodysplasia-related changes(AML-MRC)", "t-AML", "NOS")

# WHO AML classification
Vars_WHO_AML <- c("PML::RARA fusion", "RUNX1::RUNX1T1 fusion", "CBFB::MYH11 fusion", 
                  "DEK::NUP214 fusion", "RBM15::MRTFA fusion", "BCR::ABL1 fusion",
                  "KMT2A rearrangement", "MECOM rearrangement", "NUP98 rearrangement",
                  "NPM1 mutation", "biCEBPA mutation", "myelodysplasia-related(AML-MR)",
                  "AML,defined by differentiation")

# AML FAB simplified
vars_AML_fab1 <- c("AML/M0", "AML/M1", "AML/M2", "AML/M3", "AML/M4", "AML/M5", 
                   "AML/M6", "AML/M7", "AML/未分型")

# WHO B-ALL classification
vars_WHO_ALL <- c("B-ALL伴高二倍体(51-65)", "B-ALL伴亚二倍体(30-39)", "B-ALL伴iamp21", "B-ALL伴BCR::ABL1 fusion",
                  "B-ALL伴BCR::ABL1-like", "B-ALL伴KMT2A rearrangement", "B-ALL伴ETV6::RUNX1 fusion",
                  "B-ALL伴TCF3::PBX1 fusion", "B-ALL伴IGH相关融合", "B-ALL伴TCF3::HLF fusion", "NOS")

# WHO B-ALL extended classification
vars_WHO2_ALL <- c("B-ALL伴高二倍体(51-65)", "B-ALL伴亚二倍体(30-39)", "B-ALL伴iamp21", "B-ALL伴BCR::ABL1 fusion",
                   "B-ALL伴BCR::ABL1-like", "B-ALL伴KMT2A rearrangement", "B-ALL伴ETV6::RUNX1 fusion",
                   "B-ALL伴TCF3::PBX1 fusion", "B-ALL伴IGH相关融合", "B-ALL伴TCF3::HLF fusion", 
                   "B-ALL伴other defined genetic abnormalities", "NOS")

# Additional variable lists for analyses
vars_allo_type <- c("单倍体", "脐血", "非血缘", "同胞全相合", "不详")
vars_sct_type <- c("Allo-SCT", "Auto-SCT", "类型不详")
vars_diag_year <- c("2016", "2017", "2018", "2019", "2020", "2021", "2022")
vars_diag_year1 <- c("2016年", "2017年", "2018年", "2019年", "2020年", "2021年", "2022年")
vars_ELN <- c("Favorable", "Intermediate", "Adverse")
vars_yizhi <- c('Allo-SCT', "未移植")
vars_Sanz <- c("Low", "Intermediate", "High")
vars_Sanz2 <- c("标危", "高危")
vars_al_subtype <- c("AML(non-M3)", "APL(M3)", "B-ALL", "T-ALL")
vars_al_type1 <- c("AML(non-M3)", "APL(M3)", "ALL", "Other AL")
vars_all_subtype <- c("B-ALL/pro-B", "B-ALL/common-B", "B-ALL/pre-B", "B-ALL/mature-B", "T-ALL", "ALL/未分型")
vars_al_type <- c("AML", "ALL", "MPAL", "other")
vars_aml_type <- c("AML(non-M3)", "APL(M3)")
vars_all_type <- c("B-ALL", "T-ALL", "ALL/unknown")
vars_al_type <- c("AML(non-M3)", "APL(M3)", "B-ALL", "T-ALL", "MPAL", "ALL/未分型", 
                  "AL/未分型")
vars_age_group <- c("1.儿童和青少年(15-18岁)", "2.青年(19-44岁)", "3.中年(45-59岁)", "4.年轻老年(60-74岁)", "5.老老年(≥75岁)")
vars_age_group_short <- c("15-18岁", "19-44岁", "45-59岁", "60-74岁", "≥75岁")
vars_age_group1 <- c("<1", "1-14", "15-18", "19-44", "45-59", "60-74", "75+")
vars_age_group2 <- c("0-14", "15-44", "45-59", "60-74", "75+")
vars_allage <- c("child", "adult", "old")
vars_age_group60 <- c("60周岁以下", "60周岁及以上")
vars_age_group60_e <- c("<60", "≥60")
vars_age_group10 <- c('10-19', '20-29', '30-39', '40-49', '50-59', '60-69', '70-79', '80-89', '90-99', "+100")
vars_age_group10_e <- c('15-19', '20-29', '30-39', '40-49', '50-59', '60-69', '70-79', '80+')
vars_age_group5 <- c('0-4', '5-9', '10-14', '15-19', '20-24', '25-29', '30-34', '35-39',
                     '40-44', '45-49', '50-54', '55-59', '60-64', '65-69', '70-74', '75-79',
                     '80-84', '85-89', '90-94', '95-99', '100+')
vars_WHO1_ALL <- c(vars_WHO_ALL[-11], vars_ALL_mianyi)
vars_gender <- c("男", "女")
vars_gender_e <- c("male", "female")

# ------------------------------------------------------------------------------
# 7. Define utility functions
# ------------------------------------------------------------------------------

# Function to create survival plots with risk tables
survival_plot4 <- function(surv, case, ylab, legend.title, legend.labs, palette) {
  ggsurv <- ggsurvplot(
    surv, data = case, pval = TRUE,
    risk.table = TRUE, tables.height = 0.3,
    risk.table.fontsize = 4,
    risk.table.title = NULL,
    tables.theme = theme_cleantable(),
    xlim = c(0, 60), ylim = c(0, 1),
    break.x.by = 12, break.y.by = 0.1, 
    surv.median.line = "hv", 
    palette = palette,
    xlab = "Months", ylab = ylab,
    legend.title = "Group",
    censor.size = 2, censor.shape = 124,
    legend = "top",
    legend.labs = legend.labs,
    font.x = c(14, "bold", "black"),
    font.y = c(14, "bold", "black"),
    font.tickslab = c(12, "bold", "black")
  )
  return(ggsurv)
}

# Function to extract and format survival times and rates
OS_select_func<-function(surv){
  res <- summary(surv,times=c(1,12,24,36,48,60),extend = TRUE)
  survtable <- as.data.frame(res[c("strata","time","n.risk","n.event","surv","lower","upper")])
  select <- survtable 
  select[,c("surv","lower","upper")]<-round(select[,c("surv","lower","upper")]*100,1)
  select<-select %>% mutate(OS1=paste0(surv,"(",lower ,"-", upper,")"))
  select1<-select[c("strata","time","OS1")]
  select1<-select1 %>% pivot_wider(names_from = time,values_from = OS1)
  test<-as.data.frame(summary(surv)$table[,c("records","median","0.95LCL" ,"0.95UCL")])
  test$strata=rownames(test)
  select1<-left_join(test,select1,by="strata")
  select1$median_OS=paste0(round(select1$median,1),"(",round(select1$`0.95LCL`,1) ,"-", round(select1$`0.95UCL`,1),")")
  select1$one<-str_remove_all(select1$`12`,"\\(.*")
  select1$three<-str_remove_all(select1$`36`,"\\(.*")
  select1$five<-str_remove_all(select1$`60`,"\\(.*")
  assign("select1",select1)
  return(select1)}
OS_select_func1<-function(surv){
  res <- summary(surv,times=c(1,12,24,36,48,60),extend = TRUE)
  survtable <- as.data.frame(res[c("time","n.risk","n.event","surv","lower","upper")])
  select <- survtable 
  select[,c("surv","lower","upper")]<-round(select[,c("surv","lower","upper")]*100,1)
  select<-select %>% mutate(OS1=paste0(surv,"(",lower ,"-", upper,")"))
  select1<-select[c("time","OS1")]
  select1<-select1 %>% pivot_wider(names_from = time,values_from = OS1)
  test<-summary(surv)$table[c("records","median","0.95LCL" ,"0.95UCL")]
  test<-t(as.data.frame(test))
  select1<-cbind(test,select1)
  select1$median_OS=paste0(round(select1$median,1),"(",round(select1$`0.95LCL`,1) ,"-", round(select1$`0.95UCL`,1),")")
  select1$one<-str_remove_all(select1$`12`,"\\(.*")
  select1$three<-str_remove_all(select1$`36`,"\\(.*")
  select1$five<-str_remove_all(select1$`60`,"\\(.*")
  assign("select1",select1,envir = .GlobalEnv)
  return(select1)}

# Function to perform Cox regression analysis
perform_cox_analysis <- function(data, time_var, death_var, covariates) {
  # Ensure the covariates are in the correct format
  covariate_formula <- paste(covariates, collapse = " + ")
  
  # Fit univariate Cox models
  univariate_results <- lapply(covariates, function(var) {
    model <- coxph(as.formula(paste("Surv(", time_var, ", ", death_var, ") ~", var)), data = data)
    tidy_results <- tidy(model)  # Use broom to tidy the results
    tidy_results <- tidy_results %>%
      mutate(
        HR = exp(estimate),  # Calculate HR
        conf.low = exp(estimate - 1.96 * std.error),  # Lower CI
        conf.high = exp(estimate + 1.96 * std.error)  # Upper CI
      )
    
    # Perform proportional hazards test
    ph_test <- cox.zph(model)
    ph_test_results <- ph_test
    
    # Add PH test results to the tidy results
    tidy_results <- tidy_results %>%
      mutate(
        ph_test_statistic = ph_test_results$table[1],
        ph_test_p_value = ph_test_results$table[5]
      )
    
    return(tidy_results)
  })
  
  # Combine univariate results
  univariate_results_df <- bind_rows(univariate_results, .id = "covariate")
  
  # Fit multivariate Cox model
  multivariate_model <- coxph(
    as.formula(paste("Surv(", time_var, ", ", death_var, ") ~", covariate_formula)), 
    data = data
  )
  
  multivariate_summary <- tidy(multivariate_model) %>%
    mutate(
      HR = exp(estimate),  # Calculate HR
      conf.low = exp(estimate - 1.96 * std.error),  # Lower CI
      conf.high = exp(estimate + 1.96 * std.error)  # Upper CI
    )
  
  # Perform PH test for multivariate model
  ph_test <- cox.zph(multivariate_model)
  ph_test_results1 <- ph_test$table
  p <- ggcoxzph(ph_test)
  
  # Return results
  return(list(
    univariate = univariate_results_df, 
    multivariate = multivariate_summary, 
    ph_test = ph_test_results1
  ))
}

# Function to calculate counts and proportions, removing NA values
preprocess_countsprop <- function(indat, x) {
  indat <- indat[!is.na(indat[, x]), ]
  pt_counts <- as.data.frame(table(indat[, x]))
  pt_counts <- pt_counts %>% 
    mutate(prop = paste0(round(100 * pt_counts$Freq / nrow(indat), 2), "%"))
  return(pt_counts)
}

# Function to calculate counts and proportions, removing NA values and zero counts
preprocess_countsprop1 <- function(indat, x) {
  indat <- indat[!is.na(indat[, x]), ]
  pt_counts <- as.data.frame(table(indat[, x]))
  pt_counts <- pt_counts[pt_counts$Freq != 0, ]
  pt_counts <- pt_counts %>% 
    mutate(prop = paste0(round(100 * pt_counts$Freq / nrow(indat), 2), "%"))
  return(pt_counts)
}

# Function to calculate proportions by group
preprocess_counts3 <- function(indat, x, y) {
  indat <- indat[!is.na(indat[, x]), ]
  pt_counts <- as.data.frame(table(indat[, x], indat[, y]))
  total <- as.data.frame(table(indat[, x]))
  row.names(total) <- total[, 1]
  pt_counts$total <- total[pt_counts[, 1], 2]
  pt_counts <- pt_counts %>% mutate(prop = round(Freq / total, 3))
  return(pt_counts)
}

# Function to calculate proportions by group, removing zero counts
preprocess_counts4 <- function(indat, x, y) {
  indat <- indat[!is.na(indat[, y]), ]
  indat <- indat[!is.na(indat[, x]), ]
  pt_counts <- as.data.frame(table(indat[, x], indat[, y]))
  total <- as.data.frame(table(indat[, x]))
  row.names(total) <- total[, 1]
  pt_counts$total <- total[pt_counts[, 1], 2]
  pt_counts <- pt_counts[!pt_counts$Freq == 0, ]
  pt_counts <- pt_counts %>% mutate(prop = round(Freq / total, 2))
  return(pt_counts)
}

# Bar plot functions for visualisation
barplot3 <- function(indat, x, y, z, vars, a, b) {
  options(digits = 2)
  fztb <- preprocess_counts3(indat, x, y)
  fztb <- fztb %>% filter(Var2 == z)
  fztb$Var1 <- factor(fztb$Var1, levels = vars)
  ggplot(fztb, aes(Var1, prop)) +
    geom_bar(stat = "identity", width = 0.8, fill = mycol3[4], col = "white") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 8)) +
    scale_y_continuous(labels = scales::percent) +
    xlab(a) + ylab(b) +
    geom_text(stat = "identity", label = paste0(round(fztb$prop * 100, 2), "%"),
              vjust = -0.5, size = 4, colour = "black")
}

barplot4 <- function(indat, x, y, z, vars, a, b) {
  options(digits = 2)
  fztb <- preprocess_counts3(indat, x, y)
  fztb <- fztb %>% filter(Var2 == z)
  fztb$Var1 <- factor(fztb$Var1, levels = vars)
  ggplot(fztb, aes(Var1, prop)) +
    geom_bar(stat = "identity", width = 0.8, fill = mycol3[4], col = "white") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 8)) +
    scale_y_continuous(labels = scales::percent) +
    xlab(a) + ylab(b) +
    geom_text(stat = "identity", 
              label = paste0(fztb$Freq, "例(", round(fztb$prop * 100, 2), "%", ")", sep = ""),
              vjust = -0.5, size = 3, colour = "black")
}

barplot5 <- function(indat, x, y, z, vars, a, b) {
  options(digits = 2)
  fztb <- preprocess_counts3(indat, x, y)
  fztb <- fztb %>% filter(Var2 == z)
  fztb$Var1 <- factor(fztb$Var1, levels = vars)
  ggplot(fztb, aes(prop, Var1)) +
    geom_bar(stat = "identity", width = 0.8, fill = mycol3[4], col = "white") +
    theme(axis.text.y = element_text(size = 8)) +
    scale_x_continuous(labels = scales::percent) +
    ylab(a) + xlab(b) +
    geom_text(stat = "identity", 
              label = paste0(fztb$Freq, "例(", round(fztb$prop * 100, 2), "%", ")", sep = ""),
              hjust = 1.2, size = 3, colour = "black")
}

# Stacked bar chart function for binary outcomes
barplot_stack1 <- function(indat, x, y, vars, labelx, labely, mycol, legend_title, legend_labels) {
  z <- preprocess_counts3(indat, x, y)
  z$Var1 <- factor(z$Var1, levels = vars)
  p <- ggplot(z, mapping = aes(x = Var1, y = Freq, fill = Var2)) +
    guides(fill = guide_legend(title = legend_title)) +
    geom_col(position = "fill") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 13)) + 
    scale_fill_manual(values = mycol, labels = legend_labels) +
    xlab(labelx) + ylab(labely) + 
    scale_y_continuous(labels = scales::percent_format()) +
    geom_text(aes(label = paste(round(prop * 100, 2), "%", sep = "")), 
              position = position_fill(vjust = .5), size = 4)
  return(p)
}

# Stacked bar chart function without labels
barplot_stack <- function(indat, x, y, vars, labelx, labely, mycol, legend_title) {
  z <- preprocess_counts4(indat, x, y)
  z$Var1 <- factor(z$Var1, levels = vars)
  p <- ggplot(z, mapping = aes(x = Var1, y = prop, fill = Var2)) +
    guides(fill = guide_legend(title = legend_title)) +
    geom_col(position = "fill") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 13)) + 
    scale_fill_manual(values = mycol) +
    xlab(labelx) + ylab(labely) + 
    scale_y_continuous(labels = scales::percent_format()) +
    geom_text(aes(label = paste(Freq, " (", round(prop * 100, 2), "%", ")", sep = "")), 
              position = position_fill(vjust = .5), size = 3)
  return(p)
}

# ------------------------------------------------------------------------------
# 8. Update survival dates for censoring after 31 Dec 2023
# ------------------------------------------------------------------------------

# Create new follow-up variables with consistent censoring date
case <- case %>% mutate(last_date4 = last_date, siwang_event4 = siwang_event)

# Check how many events are after the censoring date
nrow(case[!is.na(case$siwang_event) & case$siwang_event > as.Date("2023-12-31"), ])

# Update events after censoring date to be censored at 31 Dec 2023
case$last_date4[!is.na(case$siwang_event) & case$siwang_event > as.Date("2023-12-31")] <- "2023-12-31"
case$siwang_event4[!is.na(case$siwang_event) & case$siwang_event > as.Date("2023-12-31")] <- NA
case$last_date4[!is.na(case$last_date4) & case$last_date4 > as.Date("2023-12-31")] <- "2023-12-31"

# Create new survival variables
case <- case %>% mutate(OS_status4 = NA, early_death4 = NA, follow_status4 = NA)
case$OS_status4 <- ifelse(is.na(case$siwang_event4), 0, 1)

# Calculate overall survival time in days
case <- case %>% 
  mutate(OS4 = if_else(
    case$OS_status4 == 1,
    as.numeric(difftime(as.Date(siwang_event4), as.Date(final_diag_date), units = "days")),
    as.numeric(difftime(as.Date(last_date4), as.Date(final_diag_date), units = "days"))
  ))

# Set follow-up status and early death indicator
case$follow_status4 <- ifelse(is.na(case$siwang_event4), 1, 0)
case$early_death4 <- ifelse(case$OS_status4 == 1 & case$OS4 <= 30, 1, 0)

# Preserve original variables
case$OS_raw <- case$OS
case$OS_status_raw <- case$OS_status
case$follow_status_raw <- case$follow_status

# Update primary survival variables
case$OS <- case$OS4
case$OS_status <- case$OS_status4
case$follow_status <- case$follow_status4

# Check completeness of early death variables
table(is.na(case$early_death))
table(is.na(case$early_death1))

# ------------------------------------------------------------------------------
# 9. Add paediatric data and create combined dataset
# ------------------------------------------------------------------------------

# Read paediatric data
child <- read.xlsx('240307death -leukemia.xlsx')
child1 <- read.xlsx('C:\\Users\\data_analysis\\Desktop\\Import\\2025-01-01_儿童医学中心匹配CDC数据/death -leukemia的副本.xlsx')

# Prepare adult data for merging
child <- child %>% select(-is_sfz)
allage <- case[, names(child)]

# Combine adult and paediatric data
allage <- rbind(allage, child)
table(allage$age_group)

# Create age group categories
allage <- allage %>% mutate(group = ifelse(age_group %in% c("<1", "1-14"), "child", "adult"))

# Further categorize age groups
allage$age_group1 <- as.factor(allage$age_group)
allage$age_group2 <- as.character(allage$age_group1)
allage <- allage %>% 
  mutate(age_group2 = case_when(
    age_group1 %in% c("<1", "1-14") ~ '0-14',
    age_group1 %in% c("15-19", "20-39") ~ '15-39',
    TRUE ~ age_group2
  ))

# Create elderly category
allage <- allage %>% 
  mutate(group1 = if_else(age_group1 %in% c('60-74', '75+'), 'old', group))

# Set factor levels for age groups
allage$age_group1 <- factor(allage$age_group1, levels = vars_age_group1)

# Create time period variable
allage <- allage %>% 
  mutate(trend = case_when(
    diag_year %in% c('2019', '2020') ~ '2019-2020',
    TRUE ~ '2016-2018'
  )) 

# Remove cases with missing survival data
case <- case %>% filter(!is.na(OS)) # 3 newly registered cases excluded
allage <- allage %>% filter(!is.na(OS))