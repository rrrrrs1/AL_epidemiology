# Set working directory and load required libraries
setwd("C:\\Data\\rms\\data\\2412_YW")
library('mice')
source("source.R", encoding = "UTF-8")
#-----------------------------------------------------------------------------
# Table 3: Univariate Cox Regression Analysis
#-----------------------------------------------------------------------------

# Define covariates for univariate analysis
covariates <- c("trend")

# Fit survival model stratified by leukemia subtype, age group and time period
surv <- survfit(Surv(OS / 30, OS_status) ~ al_subtype_group1 + age_group1 + trend, allage)
median <- OS_select_func(surv)

# Initialize results dataframe
result <- data.frame()

# Loop through age groups and leukemia subtypes
for (a in unique(allage$age_group1)) {
  for (c in vars_al_type1) {
    # Filter data for current age group and leukemia subtype
    case1 <- allage %>% filter(
      age_group1 %in% a &
        al_subtype_group1 %in% c
    )
    case1 <- as.data.frame(case1)
    
    # Skip analysis if sample size is too small
    if (nrow(case1) < 10) {
      next
    } else {
      # Perform pairwise log-rank test for time periods
      surv_p <- pairwise_survdiff(
        Surv(OS / 30, OS_status) ~ trend,
        p.adjust.method = "bonferroni", 
        data = case1
      )
    }
    
    # Format p-values
    surv_p <- as.data.frame(as.matrix(
      ifelse(round(surv_p$p.value, 3) < 0.001, "<0.001", round(surv_p$p.value, 3))
    ))
    surv_p[is.na(surv_p)] <- "-"
    
    # Perform univariate Cox regression analysis
    sct_all <- perform_cox_analysis(case1, "OS", "OS_status", covariates)
    sct_all <- as.data.frame(sct_all$univariate)
    sct_all$strata <- paste0("al_subtype_group1=", c, ",age_group1=", a, ",trend=2019-2020")
    sct_all$log_rank_p <- surv_p$`2016-2018`
    result <- rbind(result, sct_all)
  }
}

# Format hazard ratios with confidence intervals
result$HR_CI <- paste0(round(result$HR, 2), "(", round(result$conf.low, 2), " to ", round(result$conf.high, 2), ")")

# Select relevant columns and format p-values
result1 <- result %>% select(strata, HR_CI, p.value, log_rank_p)
result1 <- result1 %>%
  mutate(
    cox_p_value = case_when(
      p.value < 0.001 ~ "p<0.001",
      p.value < 0.05 ~ "p<0.05",
      p.value >= 0.05 ~ paste0("p=", round(p.value, 2))
    )
  )

# Clean string formatting and join datasets
median$strata <- str_remove_all(median$strata, " ")
result1$strata <- str_remove_all(result1$strata, " ")
median <- left_join(median, result1)
median <- median %>% select(5, 1, 12, 6:11, 16, 18)

# Split and rename data by time periods
median0 <- median %>% filter(str_detect(strata, "2016"))
names(median0) <- paste0("trend2016_2018_", names(median0))
median1 <- median %>% filter(str_detect(strata, "2020"))
names(median1) <- paste0("trend2019_2020_", names(median1))
median <- cbind(median0[1:9], median1)

# Save univariate analysis results
write.xlsx(median, "bmj_output/bmj_table3_1.xlsx")

#-----------------------------------------------------------------------------
# Table 3: Multivariate Cox Regression Analysis with Hospital Level
#-----------------------------------------------------------------------------

# Define covariates for multivariate analysis
covariates <- c("trend", 'gender', 'diag_level_6')

# Fit survival model
surv <- survfit(Surv(OS / 30, OS_status) ~ al_subtype_group1 + age_group1 + trend, allage)
median <- OS_select_func(surv)

# Initialize results dataframe
result <- data.frame()

# Loop through age groups and leukemia subtypes
for (a in unique(allage$age_group1)) {
  for (c in vars_al_type1) {
    # Filter data for current age group and leukemia subtype
    case1 <- allage %>% filter(
      age_group1 %in% a &
        al_subtype_group1 %in% c
    )
    case1 <- as.data.frame(case1)
    
    # Skip analysis if sample size is too small
    if (nrow(case1) < 10) {
      next
    } else {
      # Perform pairwise log-rank test
      surv_p <- pairwise_survdiff(
        Surv(OS / 30, OS_status) ~ trend,
        p.adjust.method = "bonferroni", 
        data = case1
      )
    }
    
    # Format p-values
    surv_p <- as.data.frame(as.matrix(
      ifelse(round(surv_p$p.value, 3) < 0.001, "<0.001", round(surv_p$p.value, 3))
    ))
    surv_p[is.na(surv_p)] <- "-"
    
    # Perform multivariate Cox regression
    sct_all <- perform_cox_analysis(case1, "OS", "OS_status", covariates)
    sct_all <- as.data.frame(sct_all$multivariate)[1, ]
    sct_all$strata <- paste0("al_subtype_group1=", c, ",age_group1=", a, ",trend=2019-2020")
    sct_all$log_rank_p <- surv_p$`2016-2018`
    result <- rbind(result, sct_all)
  }
}

# Format hazard ratios with confidence intervals
result$HR_CI <- paste0(round(result$HR, 2), "(", round(result$conf.low, 2), " to ", round(result$conf.high, 2), ")")

# Select relevant columns and format p-values
result1 <- result %>% select(strata, HR_CI, p.value, log_rank_p)
result1 <- result1 %>%
  mutate(
    cox_p_value = case_when(
      p.value < 0.001 ~ "p<0.001",
      p.value < 0.05 ~ "p<0.05",
      p.value >= 0.05 ~ paste0("p=", round(p.value, 2))
    )
  )

# Clean string formatting and join datasets
median$strata <- str_remove_all(median$strata, " ")
result1$strata <- str_remove_all(result1$strata, " ")
median <- left_join(median, result1)
median <- median %>% select(5, 1, 12, 6:11, 16, 18)

# Split and rename data by time periods
median0 <- median %>% filter(str_detect(strata, "2016"))
names(median0) <- paste0("trend2016_2018_", names(median0))
median1 <- median %>% filter(str_detect(strata, "2020"))
names(median1) <- paste0("trend2019_2020_", names(median1))
median <- cbind(median0[1:9], median1)

# Save multivariate analysis results
write.xlsx(median, "bmj_output/bmj_table3_multicox.xlsx")

#-----------------------------------------------------------------------------
# Test Proportional Hazards Assumption
#-----------------------------------------------------------------------------

# Initialize results dataframe
result <- data.frame()
covariates <- c("trend")

# Loop through age groups and leukemia subtypes
for (a in unique(allage$age_group1)) {
  for (c in vars_al_type1) {
    # Filter data for current age group and leukemia subtype
    case1 <- allage %>% filter(
      al_subtype_group1 %in% c &
        age_group1 %in% a
    )
    case1 <- as.data.frame(case1)
    
    # Skip analysis if sample size is too small
    if (nrow(case1) < 10) {
      next
    } else {
      # Test proportional hazards assumption
      sct_all <- perform_cox_analysis(case1, "OS", "OS_status", covariates)
      cox_fit <- as.data.frame(sct_all$ph_test)
      cox_fit$strata <- paste0(a, c)
      result <- rbind(result, cox_fit)
    }
  }
}

# Evaluate proportional hazards assumption violations
result
table(result$ph_test_p_value < 0.05)

#-----------------------------------------------------------------------------
# Data Preprocessing
#-----------------------------------------------------------------------------

# Convert gender labels from Chinese to English
case$gender <- factor(case$gender, levels = c('男', '女'), labels = c('Male', 'Female'))
allage$gender <- factor(allage$gender, levels = c('男', '女'), labels = c('Male', 'Female'))

# Create age groups
case$age_group10 <- cut(case$diag_age,
                        breaks = c(9, 19, 29, 39, 49, 59, 69, Inf),
                        labels = c(vars_age_group10[1:6], "70+")
)

# Convert stem cell transplant status
case$allo_sct <- as.character(case$allo_sct)
case <- case %>% mutate(allo_sct = case_when(
  allo_sct %in% "1" ~ "Yes",
  allo_sct %in% "0" ~ "No",
  TRUE ~ allo_sct
))
case$allo_sct <- factor(case$allo_sct, c("Yes", "No"))

# Select variables for analysis
vars <- c('caseid', 'hospital', 'OS', 'OS_status', 'age_group1', 'gender', 'trend',
          'x_wbc', 'x_hb', 'x_plt', 'xintai_gs', 'al_subtype_group1', "allo_sct", 'WHO2_short')
case1 <- case %>% select(vars)

# Calculate missing data percentage
miss_percent <- colMeans(is.na(case1)) * 100
miss_n <- colSums(is.na(case1))
miss_table <- data.frame(
  Variable = names(miss_percent),
  Missing_n = miss_n,
  Missing_percent = miss_percent
)
print(miss_table)

# Perform multiple imputation
perform_imputation <- function(case1, m=5, maxit=50, seed=123) {
  imp <- mice(case1, m=m, maxit=maxit, seed=seed, printFlag=TRUE)
  complete_datasets <- lapply(1:m, function(i) complete(imp, i))
  return(list(imp=imp, complete_datasets=complete_datasets))
}

# Run imputation
imputation_result <- perform_imputation(case1)

# Select imputed dataset and create categorized version
new_case1 <- imputation_result$complete_datasets[[3]]
new_case1_cat <- new_case1

# Categorize white blood cell count
new_case1_cat$x_wbc <- cut(new_case1_cat$x_wbc,
                           breaks = c(0, 20, 100, Inf),
                           labels = c('<20', '20-100', '100+'),
                           right = FALSE)

# Define Anaemia based on gender-specific hemoglobin cutoffs
new_case1_cat <- new_case1_cat %>% mutate(
  x_hb = case_when(
    x_hb < 110 & gender == 'Female' ~ 'Anaemia',
    x_hb < 120 & gender == 'Male' ~ 'Anaemia',
    TRUE ~ 'Normal'
  )
)

# Categorize platelet count
new_case1_cat$x_plt <- cut(new_case1_cat$x_plt,
                           breaks = c(0, 40, Inf),
                           labels = c('<40', '40+'),
                           right = FALSE)

# Categorize bone marrow blast percentage
new_case1_cat$xintai_gs <- cut(new_case1_cat$xintai_gs,
                               breaks = c(0, 20, 50, 80, 100),
                               labels = c('<20%', '20-50%', '50-80%', '80+%'),
                               right = FALSE)

# Add diagnosis age to imputed datasets
new_case1 <- left_join(new_case1, case[, c('caseid', 'diag_age')])
new_case1_cat <- left_join(new_case1_cat, case[, c('caseid', 'diag_age')])

# Save imputed data
saveRDS(new_case1_cat, '250403_chabu_new_case1_cat.rds')

# Load preprocessed data
new_case1_cat <- readRDS(file = '250403_chabu_new_case1_cat.rds')
new_case1 <- readRDS(file = '250403_chabu_new_case1.rds')

# Set factor levels
new_case1_cat$x_hb <- factor(new_case1_cat$x_hb, levels = c('Normal', 'Anaemia'))
new_case1_cat$allo_sct <- factor(new_case1_cat$allo_sct, levels = c('No', 'Yes'))

# Convert variables to factors
covariates <- c('trend', 'age_group1', 'al_subtype_group1', 'gender',
                'x_plt', 'x_wbc', 'x_hb', 'xintai_gs', 'allo_sct')
new_case1_cat[covariates] <- lapply(new_case1_cat[covariates], factor)

covariates <- c('trend', 'age_group1', 'al_subtype_group1', 'gender', 'group')
allage[covariates] <- lapply(allage[covariates], factor)

covariates <- c('trend', 'age_group1', 'al_subtype_group1', 'gender', 'allo_sct')
new_case1[covariates] <- lapply(new_case1[covariates], factor)

#-----------------------------------------------------------------------------
# Table 3: Multivariate Cox Regression Analysis (without transplant)
#-----------------------------------------------------------------------------

# Define covariates for multivariate analysis
covariates <- c('trend', 'gender', 'x_plt', 'x_wbc', 'x_hb', 'xintai_gs')

# Fit survival model
surv <- survfit(Surv(OS / 30, OS_status) ~ al_subtype_group1 + age_group1 + trend, new_case1_cat)
median <- OS_select_func(surv)

# Initialize results dataframe
result <- data.frame()

# Loop through age groups and leukemia subtypes
for (a in unique(new_case1_cat$age_group1)) {
  for (c in vars_al_type1) {
    # Filter data for current age group and leukemia subtype
    case1 <- new_case1_cat %>% filter(
      age_group1 %in% a &
        al_subtype_group1 %in% c
    )
    case1 <- as.data.frame(case1)
    
    # Skip analysis if sample size is too small
    if (nrow(case1) < 10) {
      next
    } else {
      # Perform pairwise log-rank test
      surv_p <- pairwise_survdiff(
        Surv(OS / 30, OS_status) ~ trend,
        p.adjust.method = "bonferroni", 
        data = case1
      )
    }
    
    # Format p-values
    surv_p <- as.data.frame(as.matrix(
      ifelse(round(surv_p$p.value, 3) < 0.001, "<0.001", round(surv_p$p.value, 3))
    ))
    surv_p[is.na(surv_p)] <- "-"
    
    # Perform multivariate Cox regression
    sct_all <- perform_cox_analysis(case1, "OS", "OS_status", covariates)
    sct_all <- as.data.frame(sct_all$multivariate)[1, ]
    sct_all$strata <- paste0("al_subtype_group1=", c, ",age_group1=", a, ",trend=2019-2020")
    sct_all$log_rank_p <- surv_p$`2016-2018`
    result <- rbind(result, sct_all)
  }
}

# Format hazard ratios with confidence intervals
result$HR_CI <- paste0(round(result$HR, 2), "(", round(result$conf.low, 2), " to ", round(result$conf.high, 2), ")")

# Select relevant columns and format p-values
result1 <- result %>% select(strata, HR_CI, p.value, log_rank_p)
result1 <- result1 %>%
  mutate(
    cox_p_value = case_when(
      p.value < 0.001 ~ "p<0.001",
      p.value < 0.05 ~ "p<0.05",
      p.value >= 0.05 ~ paste0("p=", round(p.value, 2))
    )
  )

# Clean string formatting and join datasets
median$strata <- str_remove_all(median$strata, " ")
result1$strata <- str_remove_all(result1$strata, " ")
median <- left_join(median, result1)
median <- median %>% select(5, 1, 12, 6:11, 16, 18)

# Split and rename data by time periods
median0 <- median %>% filter(str_detect(strata, "2016"))
names(median0) <- paste0("trend2016_2018_", names(median0))
median1 <- median %>% filter(str_detect(strata, "2020"))
names(median1) <- paste0("trend2019_2020_", names(median1))
median <- cbind(median0[1:9], median1)

# Save multivariate analysis results (without transplant)
write.xlsx(median, "bmj_output/bmj_table3_multicox_withoutyizhi.xlsx")

#-----------------------------------------------------------------------------
# Table 3: Multivariate Cox Regression Analysis (with transplant)
#-----------------------------------------------------------------------------

# Add transplant to covariates for multivariate analysis
covariates <- c('trend', 'gender', 'x_plt', 'x_wbc', 'x_hb', 'xintai_gs', 'allo_sct')

# Fit survival model with selected leukemia subtypes (subtypes 1 and 3)
surv <- survfit(
  Surv(OS / 30, OS_status) ~ al_subtype_group1 + age_group1 + trend, 
  new_case1_cat[new_case1_cat$al_subtype_group1 %in% vars_al_type1[c(1, 3)], ]
)
median <- OS_select_func(surv)

# Initialize results dataframe
result <- data.frame()

# Loop through age groups and selected leukemia subtypes
for (a in unique(new_case1_cat$age_group1)) {
  for (c in vars_al_type1[c(1, 3)]) {
    # Filter data for current age group and leukemia subtype
    case1 <- new_case1_cat %>% filter(
      age_group1 %in% a &
        al_subtype_group1 %in% c
    )
    case1 <- as.data.frame(case1)
    
    # Skip analysis if sample size is too small
    if (nrow(case1) < 10) {
      next
    } else {
      # Perform pairwise log-rank test
      surv_p <- pairwise_survdiff(
        Surv(OS / 30, OS_status) ~ trend,
        p.adjust.method = "bonferroni", 
        data = case1
      )
    }
    
    # Format p-values
    surv_p <- as.data.frame(as.matrix(
      ifelse(round(surv_p$p.value, 3) < 0.001, "<0.001", round(surv_p$p.value, 3))
    ))
    surv_p[is.na(surv_p)] <- "-"
    
    # Perform multivariate Cox regression with transplant
    sct_all <- perform_cox_analysis(case1, "OS", "OS_status", covariates)
    sct_all <- as.data.frame(sct_all$multivariate)[1, ]
    sct_all$strata <- paste0("al_subtype_group1=", c, ",age_group1=", a, ",trend=2019-2020")
    sct_all$log_rank_p <- surv_p$`2016-2018`
    result <- rbind(result, sct_all)
  }
}

# Format hazard ratios with confidence intervals
result$HR_CI <- paste0(round(result$HR, 2), "(", round(result$conf.low, 2), " to ", round(result$conf.high, 2), ")")

# Select relevant columns and format p-values
result1 <- result %>% select(strata, HR_CI, p.value, log_rank_p)
result1 <- result1 %>%
  mutate(
    cox_p_value = case_when(
      p.value < 0.001 ~ "p<0.001",
      p.value < 0.05 ~ "p<0.05",
      p.value >= 0.05 ~ paste0("p=", round(p.value, 2))
    )
  )

# Clean string formatting and join datasets
median$strata <- str_remove_all(median$strata, " ")
result1$strata <- str_remove_all(result1$strata, " ")
median <- left_join(median, result1)
median <- median %>% select(5, 1, 12, 6:11, 16, 18)

# Split and rename data by time periods
median0 <- median %>% filter(str_detect(strata, "2016"))
names(median0) <- paste0("trend2016_2018_", names(median0))
median1 <- median %>% filter(str_detect(strata, "2020"))
names(median1) <- paste0("trend2019_2020_", names(median1))
median <- cbind(median0[1:9], median1)

# Save multivariate analysis results (with transplant)
write.xlsx(median, "bmj_output/bmj_table3_multicox_withyizhi.xlsx")