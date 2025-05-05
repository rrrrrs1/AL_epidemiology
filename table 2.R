# Set working directory and source required functions
setwd("C:\\Data\\rms\\data\\2412_YW")
source("21_source_paper_sup_changeage.R", encoding = "UTF-8")

# Function to calculate and format survival statistics at specific time points
OS_select_func <- function(surv) {
  # Calculate survival at specific time points (1, 12, 24, 36, 48, 60 months)
  res <- summary(surv, times = c(1, 12, 24, 36, 48, 60), extend = TRUE)
  survtable <- as.data.frame(res[c("strata", "time", "n.risk", "n.event", "surv", "lower", "upper")])
  
  # Format survival estimates with confidence intervals
  select <- survtable
  select[, c("surv", "lower", "upper")] <- round(select[, c("surv", "lower", "upper")] * 100, 1)
  select <- select %>% mutate(OS1 = paste0(surv, "(", lower, " to ", upper, ")"))
  
  # Pivot data to wide format
  select1 <- select[c("strata", "time", "OS1")]
  select1 <- select1 %>% pivot_wider(names_from = time, values_from = OS1)
  
  # Extract median survival statistics
  test <- as.data.frame(summary(surv)$table[, c("records", "median", "0.95LCL", "0.95UCL")])
  test$strata <- rownames(test)
  
  # Join data and format median OS
  select1 <- left_join(test, select1, by = "strata")
  select1$median_OS <- paste0(round(select1$median, 1), "(", round(select1$`0.95LCL`, 1), " to ", round(select1$`0.95UCL`, 1), ")")
  
  # Extract key survival percentages without confidence intervals
  select1$one <- str_remove_all(select1$`12`, "\\(.*")
  select1$three <- str_remove_all(select1$`36`, "\\(.*")
  select1$five <- str_remove_all(select1$`60`, "\\(.*")
  
  return(select1)
}

# Function to calculate overall survival statistics (without stratification)
OS_select_func1 <- function(surv) {
  # Calculate survival at specific time points
  res <- summary(surv, times = c(1, 12, 24, 36, 48, 60), extend = TRUE)
  survtable <- as.data.frame(res[c("time", "n.risk", "n.event", "surv", "lower", "upper")])
  
  # Format survival estimates with confidence intervals
  select <- survtable
  select[, c("surv", "lower", "upper")] <- round(select[, c("surv", "lower", "upper")] * 100, 1)
  select <- select %>% mutate(OS1 = paste0(surv, "(", lower, " to ", upper, ")"))
  
  # Pivot data to wide format
  select1 <- select[c("time", "OS1")]
  select1 <- select1 %>% pivot_wider(names_from = time, values_from = OS1)
  
  # Extract and transform overall survival statistics
  test <- summary(surv)$table[c("records", "median", "0.95LCL", "0.95UCL")]
  test <- t(as.data.frame(test))
  select1 <- cbind(test, select1)
  
  # Format median OS and extract key percentages
  select1$median_OS <- paste0(round(select1$median, 1), "(", round(select1$`0.95LCL`, 1), " to ", round(select1$`0.95UCL`, 1), ")")
  select1$one <- str_remove_all(select1$`12`, "\\(.*")
  select1$three <- str_remove_all(select1$`36`, "\\(.*")
  select1$five <- str_remove_all(select1$`60`, "\\(.*")
  
  # Add strata name for consistency with other results
  select1$strata <- "allages"
  
  return(select1)
}

# Main analysis - Calculate survival rates by group and subtype
case1 <- allage  # Use the full dataset

# 1. Calculate survival stratified by both group and subtype
surv <- survfit(Surv(OS / 30, OS_status) ~ group + al_subtype_group1, case1)
survival1 <- OS_select_func(surv)

# 2. Calculate survival stratified by subtype only
surv <- survfit(Surv(OS / 30, OS_status) ~ al_subtype_group1, case1)
survival2 <- OS_select_func(surv)

# 3. Calculate survival stratified by group only
surv <- survfit(Surv(OS / 30, OS_status) ~ group, case1)
survival3 <- OS_select_func(surv)

# 4. Calculate overall survival (no stratification)
surv <- survfit(Surv(OS / 30, OS_status) ~ NULL, case1)
survival4 <- OS_select_func1(surv)

# Combine all survival results
survival <- rbind(survival1, survival2, survival3, survival4)
write.csv(survival, "paper/AL_survival_rates.csv")

# Calculate median follow-up time and demographics
# Stratified by group and subtype
surv <- survfit(Surv(OS / 30, follow_status) ~ group + al_subtype_group1, case1)
age1 <- case1 %>% 
  group_by(group, al_subtype_group1) %>% 
  summarise(age = paste0(median(diag_age), "(", quantile(diag_age, 0.25), '-', quantile(diag_age, 0.75), ')'))
survival1 <- summary(surv)$table[, c("records", "events", "median", "0.95LCL", "0.95UCL")]

# Stratified by subtype only
surv <- survfit(Surv(OS / 30, follow_status) ~ al_subtype_group1, case1)
age2 <- case1 %>% 
  group_by(al_subtype_group1) %>% 
  summarise(age = paste0(median(diag_age), "(", quantile(diag_age, 0.25), '-', quantile(diag_age, 0.75), ')'))
survival2 <- summary(surv)$table[, c("records", "events", "median", "0.95LCL", "0.95UCL")]

# Stratified by group only
surv <- survfit(Surv(OS / 30, follow_status) ~ group, case1)
age3 <- case1 %>% 
  group_by(group) %>% 
  summarise(age = paste0(median(diag_age), "(", quantile(diag_age, 0.25), '-', quantile(diag_age, 0.75), ')'))
survival3 <- summary(surv)$table[, c("records", "events", "median", "0.95LCL", "0.95UCL")]

# Overall statistics (no stratification)
surv <- survfit(Surv(OS / 30, follow_status) ~ NULL, case1)
age4 <- case1 %>% 
  summarise(age = paste0(median(diag_age), "(", quantile(diag_age, 0.25), '-', quantile(diag_age, 0.75), ')'))
survival4 <- summary(surv)$table[c("records", "events", "median", "0.95LCL", "0.95UCL")]

# Combine all follow-up results
survival <- rbind(survival1, survival2, survival3, survival4)
age_all <- rbind(age1, age2, age3, age4)
survival <- as.data.frame(survival)
survival$age <- age_all$age
write.csv(survival, "paper/AL_median_followup.csv")

# Combine survival rates and follow-up data into a final table
t1 <- read.csv("paper/AL_survival_rates.csv")
t2 <- read.csv("paper/AL_median_followup.csv")

# Format median follow-up time
t2$median_follow_up <- paste0(round(t2$median, 1), "(", round(t2$X0.95LCL, 1), " to ", round(t2$X0.95UCL, 1), ")")

# Select and join relevant columns
t2 <- t2 %>% select(records, age, median_follow_up)
t2 <- left_join(t2, t1)
t2 <- t2 %>% select(strata, records, age, median_follow_up, median_OS, X1:X60)

# Reorder rows for better presentation
t2 <- t2[c(15, 9:12, 14, 5:8, 13, 1:4), ]

# Replace "NA(NA to NA)" with "not reached" for more clarity
t2$median_OS[t2$median_OS %in% "NA(NA to NA)"] <- "not reached"

# Save the final table to Excel
write.xlsx(t2, "bmj_output/bmj_table2.xlsx")