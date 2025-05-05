# Clear environment
rm(list=ls())

# Load required packages for data analysis
library(dplyr)          # Data manipulation
library(tidyr)          # Data reshaping
library(ggplot2)        # Visualization
library(openxlsx)       # Excel file operations
library(sf)             # Spatial data handling
library(RColorBrewer)   # Color palettes
library(epitools)       # Epidemiological tools
library(scales)         # Scale transformations
library(epiR)           # Epidemiological analyses
library(reshape2)       # Data reshaping
library(lubridate)      # Date handling
library(MetBrewer)      # Color palettes
library(showtext)       # Font rendering
library(tidyverse)      # Collection of data science packages

# Enable text rendering
showtext.auto(enable = T)

# Set working directory
setwd('/Users/yinwei/Desktop/paper数据分析汇总/paper-流行病学')

# Define variables and colors for visualization
vars_cause = c("non-APL AML", "APL", 
               'ALL', 'other AL',
               'CML', 'CLL', 'CMML', 'other CL',
               'other L',
               'leukemia')
vars_al <- c('acute leukemia', vars_cause[1:4])
vars_age <- c("0", "1", "5", "10", "15", "20", "25", "30", "35", "40", "45", "50", "55", "60", "65", "70", "75", "80", "85")

# Color palettes
colors = rev(met.brewer("VanGogh2"))
colors1 = rgb(t(col2rgb(colors) / 255), alpha = 0.8)
color1 = rgb(t(col2rgb(colors) / 255), alpha = 0.8)

# Function to add cause rows to dataframe
add_cause_rows <- function(df) {
  unique_causes <- unique(df$cause)
  result <- data.frame()
  
  for (cause_group in unique_causes) {
    # Extract data for the current cause group
    group_data <- df %>% filter(cause == cause_group)
    # Create new row with cause information
    cause_row <- as.data.frame(matrix(NA, nrow = 1, ncol = ncol(df)))
    colnames(cause_row) <- colnames(df)
    cause_row[1, 1] <- cause_group
    # Bind cause row and group data to result
    result <- bind_rows(result, cause_row, group_data)
  }
  # Remove cause column
  result <- result %>% select(-cause)
  return(result)
}

###################### Table 1 ######################
########Import data for 0-14 age group########
NCC <- read.xlsx('241005NCC_CR_ASR_AL_level6(儿校)-0-14.xlsx')

# Import dictionaries
doi_dic <- read.xlsx('/Users/yinwei/Desktop/paper数据分析汇总/paper-流行病学/241001NCC字典.xlsx', sheet = 'doi')
code2_dic <- read.xlsx('/Users/yinwei/Desktop/paper数据分析汇总/paper-流行病学/241001NCC字典.xlsx', sheet = 'code2') %>% 
  select(code2, province)
kind1_dic <- read.xlsx('/Users/yinwei/Desktop/paper数据分析汇总/paper-流行病学/241001NCC字典.xlsx', sheet = 'kind1')
kind2_dic <- data.frame(
  SDI_type = c("High SDI", "Middle high SDI", "Middle SDI", "Middle low SDI", "Low SDI", "Very low SDI", 'China'),
  kind2 = 6:0
)
sex_dic <- read.xlsx('/Users/yinwei/Desktop/paper数据分析汇总/paper-流行病学/241001NCC字典.xlsx', sheet = 'sex')

# Function to join dictionary mappings
add_mappings <- function(data) {
  data <- left_join(data, doi_dic)
  data <- left_join(data, code2_dic)
  data <- left_join(data, kind1_dic)
  data <- left_join(data, kind2_dic)
  data <- left_join(data, sex_dic)
  return(data)
}

########Calculate age-standardized incidence and mortality rates########
subtype_rate <- c()
for (i in unique(NCC$doi)) {
  for(h in unique(NCC$kind2)) {
    for(j in unique(NCC$kind1)) {
      for(g in unique(NCC$sex)) {
        for(k in c('Allages', 'Child', 'Adult')) {
          a = NCC %>% filter(doi==i) %>% 
            filter(cause %in% vars_al) %>% 
            filter(code2==0) %>% 
            filter(sex==g) %>% 
            filter(kind2==h) %>% 
            filter(kind1==j) %>% 
            filter(age_G %in% k)
          
          if (nrow(a) < 1) {
            next
          }
          
          a <- a %>% arrange(type2)
          total = sum(a$case[!a$cause=='acute leukemia'])
          a$prop <- round(a$case/total, 4)
          a$prop1 <- paste0(a$case, '(', round(a$case/total, 4)*100, '%)')
          a$prop2 <- paste0(round(a$est, 2), '(', round(a$case/total, 4)*100, '%)')
          subtype_rate <- rbind(subtype_rate, a)
        }
      }
    }
  }
}

subtype_rate1 <- add_mappings(subtype_rate)
subtype_rate2 <- subtype_rate1 %>% 
  filter(kind1_cn=='合计' & SDI_type=='China') %>%
  select(doi_en, kind1_en, SDI_type, sex_en, cause, age_G, rate:prop2)

########Calculate Mortality-to-Incidence Ratio (MIR)########
subtype_rate <- c()
for(h in unique(NCC$kind2)) {
  for(j in unique(NCC$kind1)) {
    for(g in unique(NCC$sex)) {
      for(k in c('Allages', 'Child', 'Adult')) {
        a = NCC %>% 
          filter(cause %in% vars_al) %>% 
          filter(code2==0) %>% 
          filter(sex==g) %>% 
          filter(kind2==h) %>% 
          filter(kind1==j) %>% 
          filter(age_G %in% k)
        
        a$cause <- factor(a$cause, levels=vars_al)
        a <- a %>% arrange(cause)
        
        if (nrow(a) < 1) {
          next
        }
        
        # Calculate MIR
        NCC_ratio <- a %>% select(doi, adj.rate, cause)
        names(NCC_ratio)[2] <- 'ASR'
        NCC_ratio$ASR <- as.numeric(NCC_ratio$ASR)
        NCC_ratio <- NCC_ratio %>% pivot_wider(names_from = doi, values_from = ASR)
        names(NCC_ratio)[1:3] <- c('Type', 'Incidence', 'Mortality')
        NCC_ratio$MIR <- NCC_ratio$Mortality / NCC_ratio$Incidence
        NCC_ratio$Mortality <- round(NCC_ratio$Mortality, 2)
        NCC_ratio$Incidence <- round(NCC_ratio$Incidence, 2)
        NCC_ratio$MIR1 <- paste0(round(NCC_ratio$MIR, 3)*100, '%')
        NCC_ratio$kind2 = h
        NCC_ratio$kind1 = j
        NCC_ratio$sex = g
        NCC_ratio$age_G = k
        subtype_rate <- rbind(subtype_rate, NCC_ratio)
      }
    }
  }
}

data <- subtype_rate
data <- left_join(data, kind1_dic)
data <- left_join(data, kind2_dic)
data <- left_join(data, sex_dic)
data1<-data %>% filter(kind1_cn=='合计'&SDI_type=='China')
data1<-data1 %>% select(-kind1,-kind2,-sex,-sex_cn,-kind1_cn )

# Combine incidence, mortality and MIR data
subtype_rate2 <- subtype_rate2 %>% select(doi_en, sex_en, cause, age_G, ASR)
data1 <- data1 %>% select(sex_en, Type, age_G, MIR)
names(data1)[2] <- 'cause'
names(data1)[4] <- 'ASR'
data1$ASR <- round(data1$ASR, 2)
data1$doi_en <- 'MIR'
addall <- rbind(subtype_rate2, data1)

# Reshape data for presentation
addall <- addall %>% pivot_wider(names_from = c(sex_en, doi_en), values_from = ASR)

my_vars <- function() {
  c(any_of(c("age_G", "cause")), starts_with("All"), starts_with('Male'), starts_with("Female"))
}
addall <- dplyr::select(addall, my_vars())
addall$cause <- factor(addall$cause, levels=vars_al)
addall <- addall %>% arrange(cause)
addall <- add_cause_rows(addall)

# Rename columns
names(addall) <- c(
  'age_G', 'All-sex_ASIR_inc', 'All-sex_ASMR_mor', 'ALL-sex_MIR',
  'Male_ASIR_inc', 'Male_ASMR_mor', 'Male_MIR',
  'Female_ASIR_inc', 'Female_ASMR_mor', 'Female_MIR'
)

addall_ASR = addall
addall_ASR$age_G[addall_ASR$age_G=='Child'] = "Child (0-14)"
addall_ASR$age_G[addall_ASR$age_G=='Adult'] = "Adult (≥15)"

############Import estimated cases data############
estimated_case <- read.xlsx('241005NCC_CR_ASR_AL_level6(儿校)_estimtaed_cases-0-14.xlsx')
estimated_case <- estimated_case %>% 
  filter(!is.na(estimated_case)) %>% 
  select(doi, kind1, kind2, sex, cause, age_G, estimated_case)

# Process case count and crude rates
subtype_rate <- c()
for (i in unique(NCC$doi)) {
  for(h in unique(NCC$kind2)) {
    for(j in unique(NCC$kind1)) {
      for(g in unique(NCC$sex)) {
        for(k in c('Allages', 'Child', 'Adult')) {
          a = NCC %>% 
            filter(doi==i) %>% 
            filter(cause %in% vars_al) %>% 
            filter(code2==0) %>% 
            filter(sex==g) %>% 
            filter(kind2==h) %>% 
            filter(kind1==j) %>% 
            filter(age_G %in% k)
          
          if (nrow(a) < 1) {
            next
          }
          
          a <- a %>% arrange(type2)
          total = sum(a$case[!a$cause=='acute leukemia'])
          a$prop <- round(a$case/total, 4)
          a$prop1 <- paste0(a$case, '(', round(a$case/total, 4)*100, '%)')
          a$prop2 <- paste0(round(a$est, 2), '(', round(a$case/total, 4)*100, '%)')
          subtype_rate <- rbind(subtype_rate, a)
        }
      }
    }
  }
}

subtype_rate1 <- add_mappings(subtype_rate)
subtype_rate2 <- subtype_rate1 %>% 
  filter(kind1_cn=='合计' & SDI_type=='China') %>%
  left_join(estimated_case) %>%
  ungroup() %>%
  select(doi_en, sex_en, cause, age_G, prop1, estimated_case, rate_CI)

# Create the final tables
addall <- subtype_rate2 %>% 
  pivot_wider(names_from = c(sex_en), values_from = c(prop1, estimated_case, rate_CI)) %>%
  select(
    doi_en, age_G, cause,
    `prop1_All-sex`, `estimated_case_All-sex`, `rate_CI_All-sex`,
    prop1_Male, estimated_case_Male, rate_CI_Male,
    prop1_Female, estimated_case_Female, rate_CI_Female
  )

addall$cause <- factor(addall$cause, levels=vars_al)
addall <- addall %>% arrange(cause)

names(addall) <- c(
  "doi_en", 'age_G', 'cause', 'All-sex_case', 'All-sex_estimated_case', 'All-sex_crude_rate',
  'Male_case', 'Male_estimated_case', 'Male_crude_rate',
  'Female_case', 'Female_estimated_case', 'Female_crude_rate'
)

# Split into incidence and mortality tables
addall_inc <- addall %>% 
  filter(doi_en %in% 'Incidence') %>% 
  select(-doi_en) %>%
  add_cause_rows()

addall_mor <- addall %>% 
  filter(doi_en=='Mortality') %>% 
  select(-doi_en) %>%
  add_cause_rows()

# Convert all columns to character and replace NA with empty string
addall_inc <- addall_inc %>% 
  mutate(across(everything(), as.character)) %>%
  mutate_all(~replace_na(., ""))

addall_mor <- addall_mor %>% 
  mutate(across(everything(), as.character)) %>%
  mutate_all(~replace_na(., ""))

# Update age group labels
addall_inc$age_G[addall_inc$age_G=='Child'] = "Child (0-14)"
addall_inc$age_G[addall_inc$age_G=='Adult'] = "Adult (≥15)"
addall_mor$age_G[addall_mor$age_G=='Child'] = "Child (0-14)"
addall_mor$age_G[addall_mor$age_G=='Adult'] = "Adult (≥15)"

# Save first set of tables
addall_ASR1 <- addall_ASR
addall_inc1 <- addall_inc
addall_mor1 <- addall_mor

# Create combined tables for incidence
result_inc1 <- bind_cols(
  addall_inc1,
  addall_ASR1 %>% select(ends_with('inc'))
) %>% select(age_G, starts_with("All"), starts_with("Male"), starts_with("Female"))

# Create combined tables for mortality
result_mor1 <- bind_cols(
  addall_mor1,
  addall_ASR1 %>% select(ends_with('mor'))
) %>% select(age_G, starts_with("All"), starts_with("Male"), starts_with("Female"))

###################### Repeat for 0-19 age group ######################
# Import data for 0-19 age group
NCC <- read.xlsx('241005NCC_CR_ASR_AL_level6(儿校)-0-19.xlsx')

# Import dictionaries (same as before)
doi_dic <- read.xlsx('/Users/yinwei/Desktop/paper数据分析汇总/paper-流行病学/241001NCC字典.xlsx', sheet = 'doi')
code2_dic <- read.xlsx('/Users/yinwei/Desktop/paper数据分析汇总/paper-流行病学/241001NCC字典.xlsx', sheet = 'code2') %>% 
  select(code2, province)
kind1_dic <- read.xlsx('/Users/yinwei/Desktop/paper数据分析汇总/paper-流行病学/241001NCC字典.xlsx', sheet = 'kind1')
kind2_dic <- data.frame(
  SDI_type = c("High SDI", "Middle high SDI", "Middle SDI", "Middle low SDI", "Low SDI", "Very low SDI", 'China'),
  kind2 = 6:0
)
sex_dic <- read.xlsx('/Users/yinwei/Desktop/paper数据分析汇总/paper-流行病学/241001NCC字典.xlsx', sheet = 'sex')

# Add mapping function (same as before)
add_mappings <- function(data) {
  data <- left_join(data, doi_dic)
  data <- left_join(data, code2_dic)
  data <- left_join(data, kind1_dic)
  data <- left_join(data, kind2_dic)
  data <- left_join(data, sex_dic)
  return(data)
}

########Calculate age-standardized incidence and mortality rates for 0-19 data########
subtype_rate <- c()
for (i in unique(NCC$doi)) {
  for(h in unique(NCC$kind2)) {
    for(j in unique(NCC$kind1)) {
      for(g in unique(NCC$sex)) {
        for(k in c('Allages', 'Child', 'Adult')) {
          a = NCC %>% filter(doi==i) %>% 
            filter(cause %in% vars_al) %>% 
            filter(code2==0) %>% 
            filter(sex==g) %>% 
            filter(kind2==h) %>% 
            filter(kind1==j) %>% 
            filter(age_G %in% k) 
          if (nrow(a) < 1) {
            next
          }
          a <- a %>% arrange(type2)
          total = sum(a$case[!a$cause=='acute leukemia'])
          a$prop <- round(a$case/total, 4)
          a$prop1 <- paste0(a$case, '(', round(a$case/total, 4)*100, '%)')
          a$prop2 <- paste0(round(a$est, 2), '(', round(a$case/total, 4)*100, '%)')
          subtype_rate <- rbind(subtype_rate, a)
        }
      }
    }
  }
}

subtype_rate1 <- add_mappings(subtype_rate)
subtype_rate2 <- subtype_rate1 %>% 
  filter(kind1_cn=='合计' & SDI_type=='China') %>%
  select(doi_en, kind1_en, SDI_type, sex_en, cause, age_G, rate:prop2)

########Calculate MIR for 0-19 data########
subtype_rate <- c()
for(h in unique(NCC$kind2)) {
  for(j in unique(NCC$kind1)) {
    for(g in unique(NCC$sex)) {
      for(k in c('Allages', 'Child', 'Adult')) {
        a = NCC %>% 
          filter(cause %in% vars_al) %>% 
          filter(code2==0) %>% 
          filter(sex==g) %>% 
          filter(kind2==h) %>% 
          filter(kind1==j) %>% 
          filter(age_G %in% k) 
        
        a$cause <- factor(a$cause, levels=vars_al)
        a <- a %>% arrange(cause)
        
        if (nrow(a) < 1) {
          next
        }
        
        NCC_ratio <- a %>% select(doi, adj.rate, cause) 
        names(NCC_ratio)[2] <- 'ASR'
        NCC_ratio$ASR <- as.numeric(NCC_ratio$ASR)
        NCC_ratio <- NCC_ratio %>% pivot_wider(names_from = doi, values_from = ASR)
        names(NCC_ratio)[1:3] <- c('Type', 'Incidence', 'Mortality')
        NCC_ratio$MIR <- NCC_ratio$Mortality / NCC_ratio$Incidence
        NCC_ratio$Mortality <- round(NCC_ratio$Mortality, 2)
        NCC_ratio$Incidence <- round(NCC_ratio$Incidence, 2)
        NCC_ratio$MIR1 <- paste0(round(NCC_ratio$MIR, 3)*100, '%')
        NCC_ratio$kind2 = h
        NCC_ratio$kind1 = j
        NCC_ratio$sex = g
        NCC_ratio$age_G = k
        subtype_rate <- rbind(subtype_rate, NCC_ratio)
      }
    }
  }
}

data <- subtype_rate
data <- left_join(data, kind1_dic)
data <- left_join(data, kind2_dic)
data <- left_join(data, sex_dic)
data1<-data %>% filter(kind1_cn=='合计'&SDI_type=='China')
data1<-data1 %>% select(-kind1,-kind2,-sex,-sex_cn,-kind1_cn )

# Combine incidence, mortality and MIR data for 0-19
subtype_rate2 <- subtype_rate2 %>% select(doi_en, sex_en, cause, age_G, ASR)
data1 <- data1 %>% select(sex_en, Type, age_G, MIR)
names(data1)[2] <- 'cause'
names(data1)[4] <- 'ASR'
data1$ASR <- round(data1$ASR, 2)
data1$doi_en <- 'MIR'
addall <- rbind(subtype_rate2, data1)

addall <- addall %>% pivot_wider(names_from = c(sex_en, doi_en), values_from = ASR)
addall <- dplyr::select(addall, my_vars())
addall$cause <- factor(addall$cause, levels=vars_al)
addall <- addall %>% arrange(cause)
addall <- add_cause_rows(addall)

names(addall) <- c(
  'age_G', 'All-sex_ASIR_inc', 'All-sex_ASMR_mor', 'ALL-sex_MIR',
  'Male_ASIR_inc', 'Male_ASMR_mor', 'Male_MIR',
  'Female_ASIR_inc', 'Female_ASMR_mor', 'Female_MIR'
)

addall_ASR = addall
addall_ASR$age_G[addall_ASR$age_G=='Child'] = "Child (0-19)"
addall_ASR$age_G[addall_ASR$age_G=='Adult'] = "Adult (≥20)"

############Import estimated cases data for 0-19############
estimated_case <- read.xlsx('241005NCC_CR_ASR_AL_level6(儿校)_estimtaed_cases-0-19.xlsx')
estimated_case <- estimated_case %>% 
  filter(!is.na(estimated_case)) %>% 
  select(doi, kind1, kind2, sex, cause, age_G, estimated_case)

# Process case count and crude rates for 0-19
subtype_rate <- c()
for (i in unique(NCC$doi)) {
  for(h in unique(NCC$kind2)) {
    for(j in unique(NCC$kind1)) {
      for(g in unique(NCC$sex)) {
        for(k in c('Allages', 'Child', 'Adult')) {
          a = NCC %>% 
            filter(doi==i) %>% 
            filter(cause %in% vars_al) %>% 
            filter(code2==0) %>% 
            filter(sex==g) %>% 
            filter(kind2==h) %>% 
            filter(kind1==j) %>% 
            filter(age_G %in% k) 
          
          if (nrow(a) < 1) {
            next
          }
          
          a <- a %>% arrange(type2)
          total = sum(a$case[!a$cause=='acute leukemia'])
          a$prop <- round(a$case/total, 4)
          a$prop1 <- paste0(a$case, '(', round(a$case/total, 4)*100, '%)')
          a$prop2 <- paste0(round(a$est, 2), '(', round(a$case/total, 4)*100, '%)')
          subtype_rate <- rbind(subtype_rate, a)
        }
      }
    }
  }
}

subtype_rate1 <- add_mappings(subtype_rate)
subtype_rate2 <- subtype_rate1 %>% 
  filter(kind1_cn=='合计' & SDI_type=='China') %>%
  left_join(estimated_case) %>%
  ungroup() %>%
  select(doi_en, sex_en, cause, age_G, prop1, estimated_case, rate_CI)

# Create the final tables for 0-19
addall <- subtype_rate2 %>% 
  pivot_wider(names_from = c(sex_en), values_from = c(prop1, estimated_case, rate_CI)) %>%
  select(
    doi_en, age_G, cause,
    `prop1_All-sex`, `estimated_case_All-sex`, `rate_CI_All-sex`,
    prop1_Male, estimated_case_Male, rate_CI_Male,
    prop1_Female, estimated_case_Female, rate_CI_Female
  )

addall$cause <- factor(addall$cause, levels=vars_al)
addall <- addall %>% arrange(cause)

names(addall) <- c(
  "doi_en", 'age_G', 'cause', 'All-sex_case', 'All-sex_estimated_case', 'All-sex_crude_rate',
  'Male_case', 'Male_estimated_case', 'Male_crude_rate',
  'Female_case', 'Female_estimated_case', 'Female_crude_rate'
)

# Split into incidence and mortality tables for 0-19
addall_inc <- addall %>% 
  filter(doi_en %in% 'Incidence') %>% 
  select(-doi_en) %>%
  add_cause_rows()

addall_mor <- addall %>% 
  filter(doi_en=='Mortality') %>% 
  select(-doi_en) %>%
  add_cause_rows()

# Convert all columns to character and replace NA with empty string
addall_inc <- addall_inc %>% 
  mutate(across(everything(), as.character)) %>%
  mutate_all(~replace_na(., ""))

addall_mor <- addall_mor %>% 
  mutate(across(everything(), as.character)) %>%
  mutate_all(~replace_na(., ""))

# Update age group labels for 0-19
addall_inc$age_G[addall_inc$age_G=='Child'] = "Child (0-19)"
addall_inc$age_G[addall_inc$age_G=='Adult'] = "Adult (≥20)"
addall_mor$age_G[addall_mor$age_G=='Child'] = "Child (0-19)"
addall_mor$age_G[addall_mor$age_G=='Adult'] = "Adult (≥20)"

# Create combined tables for incidence for 0-19
result_inc <- bind_cols(
  addall_inc,
  addall_ASR %>% select(ends_with('inc'))
) %>% select(age_G, starts_with("All"), starts_with("Male"), starts_with("Female"))

# Combine the results from both age group definitions (0-14 and 0-19)
result_inc_all <- rbind(result_inc, result_inc1)
result_inc_all <- result_inc_all[c(1,2,23,24,3,4,
                                   1+4,2+4,23+4,24+4,3+4,4+4,
                                   1+8,2+8,23+8,24+8,3+8,4+8,
                                   1+12,2+12,23+12,24+12,3+12,4+12,
                                   1+16,2+16,23+16,24+16,3+16,4+16),]

# Create combined tables for mortality for 0-19
result_mor <- bind_cols(
  addall_mor,
  addall_ASR %>% select(ends_with('mor'))
) %>% select(age_G, starts_with("All"), starts_with("Male"), starts_with("Female"))

# Combine the mortality results from both age group definitions
result_mor_all <- rbind(result_mor, result_mor1)
result_mor_all <- result_mor_all[c(1,2,23,24,3,4,
                                   1+4,2+4,23+4,24+4,3+4,4+4,
                                   1+8,2+8,23+8,24+8,3+8,4+8,
                                   1+12,2+12,23+12,24+12,3+12,4+12,
                                   1+16,2+16,23+16,24+16,3+16,4+16),]

# Export tables to Excel
write.xlsx(result_inc_all, 'bmj_output/bmj_table1_sex_inc.xlsx')
write.xlsx(result_mor_all, 'bmj_output/bmj_table1_sex_mor.xlsx')

