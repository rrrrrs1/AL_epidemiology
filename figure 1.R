# Clear environment
rm(list=ls())

###### 01 Load required packages ######
# Core data manipulation and visualization
library(dplyr)        # Data manipulation
library(tidyr)        # Data tidying
library(ggplot2)      # Data visualization
library(openxlsx)     # Excel file handling
library(svglite)      # SVG output
library(RColorBrewer) # Color palettes
library(ggsci)        # Scientific journal color palettes
library(scales)       # Scale functions for visualization
library(patchwork)    # Combine multiple plots
library(showtext)     # Font handling

# Enable automatic showtext functionality
showtext.auto(enable = TRUE)

# Add Times New Roman font
font_add("Times New Roman", "/Library/Fonts/Times New Roman.ttf")
font_add("Times New Roman Bold", "/Library/Fonts/Times New Roman Bold.ttf")
showtext_auto()

###### 02 Set working directory and define variables ######
setwd('/Users/yinwei/Desktop/paper数据分析汇总/paper-流行病学/')

# Define disease subtypes
vars_cause = c("non-APL AML", "APL",
               'ALL', 'other AL',
               'CML', 'CLL', 'CMML', 'other CL',
               'other L',
               'leukemia')
vars_al <- c('acute leukemia', vars_cause[1:4])

# Define age groups
vars_age = c("0", "1", "5", "10", "15", "20", "25", "30", "35", "40", 
             "45", "50", "55", "60", "65", "70", "75", "80", "85")

# Define age group labels
var_age_group2 <- c("0", "1-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", 
                    "35-39", "40-44", "45-49", "50-54", "55-59", "60-64", "65-69", 
                    "70-74", "75-79", "80-84", "≥85")

# Define JCO color palette
var_color7 = pal_jco(palette="default", alpha=1)(7)

###### 03 Load and prepare data ######
# Read input data
NCC <- read.xlsx('241005NCC_CR_ASR_AL_level6(儿校)-0-19.xlsx')
doi_dic <- read.xlsx('241001NCC字典.xlsx', sheet='doi')
code2_dic <- read.xlsx('241001NCC字典.xlsx', sheet='code2')
code2_dic <- code2_dic %>% select(code2, province)
kind1_dic <- read.xlsx('241001NCC字典.xlsx', sheet='kind1')
kind2_dic <- data.frame(
  SDI_type = c("High SDI", "Middle high SDI", "Middle SDI", "Middle low SDI", "Low SDI", "Very low SDI", 'China'),
  kind2 = 6:0
)
sex_dic <- read.xlsx('241001NCC字典.xlsx', sheet='sex')

# Function to add dictionary mappings to data
add_mappings <- function(data) {
  data <- left_join(data, doi_dic)
  data <- left_join(data, code2_dic)
  data <- left_join(data, kind1_dic)
  data <- left_join(data, kind2_dic)
  data <- left_join(data, sex_dic)
  return(data)
}

###### 04 Leukemia subtype analysis by age group for incidence and mortality - BMJ style ######
# Use Lancet Oncology standard color palette
lancet_colors <- pal_lancet("lanonc")(9)[c(1,3,5,7)]

# ========== Data Processing ==========
# 1. Process incidence data (doi = 1)
NCC_inc <- NCC %>%
  filter(doi == '1' & kind2 == '0' & kind1 == '0' & code2 == 0) %>%
  filter(cause %in% vars_al) %>%
  filter(age_G %in% vars_age) %>%
  mutate(
    sex_new = case_when(
      sex == 0 ~ 'Both',
      sex == 1 ~ 'Male',
      sex == 2 ~ 'Female'
    ),
    data_type = "Incidence"  # Mark as incidence data
  )

NCC_inc$rate <- round(NCC_inc$rate, 2)
inc_data <- NCC_inc %>% 
  select(cause, age_G, sex_new, rate, lwr, upr, data_type) %>%
  filter(!cause == 'acute leukemia')

# 2. Process mortality data (doi = 2)
NCC_mort <- NCC %>%
  filter(doi == '2' & kind2 == '0' & kind1 == '0' & code2 == 0) %>%
  filter(cause %in% vars_al) %>%
  filter(age_G %in% vars_age) %>%
  mutate(
    sex_new = case_when(
      sex == 0 ~ 'Both',
      sex == 1 ~ 'Male',
      sex == 2 ~ 'Female'
    ),
    data_type = "Mortality"  # Mark as mortality data
  )

NCC_mort$rate <- round(NCC_mort$rate, 2)
mort_data <- NCC_mort %>% 
  select(cause, age_G, sex_new, rate, lwr, upr, data_type) %>%
  filter(!cause == 'acute leukemia')

# 3. Combine data
combined_data <- bind_rows(inc_data, mort_data)

# 4. Recode age group labels
combined_data <- combined_data %>%
  mutate(
    age_G = factor(age_G, levels = vars_age),
    age_G = fct_recode(
      age_G,
      "0" = "0",
      "1-4" = "1",
      "5-9" = "5",
      "10-14" = "10",
      "15-19" = "15",
      "20-24" = "20",
      "25-29" = "25",
      "30-34" = "30",
      "35-39" = "35",
      "40-44" = "40",
      "45-49" = "45",
      "50-54" = "50",
      "55-59" = "55",
      "60-64" = "60",
      "65-69" = "65",
      "70-74" = "70",
      "75-79" = "75",
      "80-84" = "80",
      "85+" = "85"
    )
  )

# 5. Recode cause names and set factor levels
combined_data$cause[combined_data$cause=='other AL'] = 'Other AL'
combined_data$cause <- factor(combined_data$cause, levels = c(vars_al[2:4], "Other AL")) 
combined_data$data_type <- factor(combined_data$data_type, levels = c("Incidence", "Mortality"))
combined_data$sex_new <- factor(combined_data$sex_new, c('Both', 'Male', 'Female'))

# 6. Sort data by data type, age group, and cause
combined_data <- combined_data %>% arrange(data_type, age_G, cause)

###### 05 Create BMJ-style visualization ######
# Create optimized BMJ-style plot with confidence intervals
p_combined <- ggplot() +
  # Confidence interval ribbons
  geom_ribbon(
    data = combined_data,
    aes(x = age_G, y = rate, ymin = lwr, ymax = upr, fill = cause, group = cause),
    alpha = 0.15,
    color = NA
  ) +
  # Line plot
  geom_line(
    data = combined_data,
    aes(x = age_G, y = rate, color = cause, group = cause, linetype = cause),
    size = 0.8
  ) +
  # Data points
  geom_point(
    data = combined_data,
    aes(x = age_G, y = rate, color = cause, group = cause, shape = cause),
    size = 1.5,
    stroke = 0.7,
    fill = "white"
  ) +
  # Y-axis settings
  scale_y_continuous(
    name = "Crude rate (per 100 000)",
    expand = expansion(mult = c(0, 0.1)),
    breaks = scales::pretty_breaks(n = 5),
    limits = c(0, NA),  # Ensure Y-axis starts from 0
    labels = function(x) gsub(",", " ", format(x, big.mark = ",", scientific = FALSE))
  ) +
  # Create facets by data type and sex
  facet_grid(data_type ~ sex_new, scales = "free_y") +
  # Apply lancet oncology color palette
  scale_color_manual(values = lancet_colors) +
  scale_fill_manual(values = lancet_colors) +
  # Line types
  scale_linetype_manual(values = c("solid", "dashed", "dotdash", "longdash")) +
  # Point shapes
  scale_shape_manual(values = c(21, 22, 23, 24)) +
  # Labels
  labs(
    x = "Age group (years)"
  ) +
  # Theme settings for BMJ style
  theme_minimal() +
  theme(
    # Base font settings
    text = element_text(family = "Times New Roman Bold", color = "black", face = "bold"),
    
    # Axis settings
    axis.title = element_text(size = 14, face = "bold", family = "Times New Roman Bold"),
    axis.title.y = element_text(margin = margin(r = 10), family = "Times New Roman Bold", face = "bold"),
    axis.title.x = element_text(margin = margin(t = 10), family = "Times New Roman Bold", face = "bold"),
    axis.text = element_text(size = 8, family = "Times New Roman Bold", face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, family = "Times New Roman Bold", face = "bold"),
    axis.ticks = element_line(size = 0.3, color = "black"),
    axis.ticks.length = unit(0.15, "cm"),
    axis.line = element_line(size = 0.3, color = "black"),
    
    # Panel settings
    panel.grid.major.y = element_line(color = "gray92", size = 0.2),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.4),
    panel.background = element_rect(fill = "white", color = NA),
    panel.spacing = unit(1, "lines"),
    
    # Facet settings
    strip.background = element_rect(fill = "gray95", color = "black", size = 0.5),
    strip.text = element_text(face = "bold", size = 14, margin = margin(5, 0, 5, 0), family = "Times New Roman Bold"),
    strip.text.x = element_text(face = "bold", size = 14, margin = margin(5, 0, 5, 0), family = "Times New Roman Bold"),
    strip.text.y = element_text(face = "bold", size = 14, margin = margin(0, 5, 0, 5), angle = 270, family = "Times New Roman Bold"),
    
    # Legend settings
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.box = "horizontal",
    legend.background = element_rect(fill = "white", color = NA),
    legend.key = element_rect(fill = "white", color = NA),
    legend.key.size = unit(0.7, "cm"),
    legend.text = element_text(size = 12, family = "Times New Roman", face = "bold"),
    legend.margin = margin(t = 10, b = 5),
    legend.spacing.x = unit(0.3, "cm"),
    
    # Plot margins
    plot.margin = margin(t = 15, r = 15, b = 15, l = 15)
  )

# Display the plot
print(p_combined)

# Save the plot as PDF
ggsave("bmj_output/figure1.pdf", p_combined, width = 14, height = 9, scale = 0.85, device = cairo_pdf)
