# Set working directory and load source file
setwd("C:\\Data\\rms\\data\\2412_YW")
source("21_source_paper_sup_changeage.R", encoding = "UTF-8")

#------------------------------------------------------------
# Function definitions for survival plot generation
#------------------------------------------------------------

# Standard survival plot function with risk table
# Parameters:
# - surv: Survival fit object
# - title: Plot title
# - legend.labs: Labels for legend items
# - palette: Color palette for plotting
#------------------------------------------------------------
survival_plot4 <- function(surv, title, legend.labs, palette) {
  # Create ggsurvfit plot with customized elements
  p <- surv %>%
    ggsurvfit(linewidth = 0.8) +
    # Add risk table below the plot
    add_risktable(
      risktable_height = 0.25,
      risktable_stats = c("{n.risk}"),
      stats_label = list(n.risk = "Number at risk"),
      size = 4.5,
      plot.margin = margin(r = 0.1, l = 0.1, unit = "cm"),
      theme = list(
        theme_risktable_default(axis.text.y.size = 10, plot.title.size = 12),
        theme(
          plot.title = element_text(face = "bold", size = 12, family = 'serif'),
          text = element_text(family = 'serif'),
          axis.title.x = element_text(face = "bold", family = 'serif'),
          axis.title.y = element_text(family = 'serif'),
          axis.text = element_text(family = 'serif')
        )
      )
    ) +
    # Add symbols and markers
    ggsurvfit::add_risktable_strata_symbol(symbol = "\U25CF", size = 20) +
    add_censor_mark(shape = 1, size = 2, stroke = 1) +
    # Add p-value annotation
    add_pvalue(
      location = "annotation", x = 5, y = 0.1, hjust = 0,
      size = 4, caption = "log-rank {p.value}"
    ) +
    # Set labels
    labs(
      title = title,
      x = "Time from diagnosis (months)",
      y = "Overall Survival"
    ) +
    # Set axis scales
    scale_x_continuous(breaks = seq(0, 60, 12), limits = c(0, 60), expand = c(0.06, 0.06)) +
    scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1), labels = paste0(seq(0, 100, 10), "%"), expand = c(0.03, 0)) +
    # Set colors
    scale_color_manual(values = palette, labels = legend.labs) +
    scale_fill_manual(values = palette, labels = legend.labs) +
    # Apply theme and additional elements
    theme_classic() +
    guides(color = guide_legend(nrow = 1)) +
    add_quantile(y_value = 0.5) +
    theme(
      axis.text = element_text(size = 14, color = "black", family = 'serif'),
      axis.title.y = element_text(size = 14, color = "black", family = 'serif', vjust = 0.5),
      axis.title.x = element_text(size = 14, color = "black", family = 'serif'),
      panel.grid = element_blank(),
      plot.margin = margin(r = 0.1, l = 0.1, unit = "cm"),
      legend.text = element_text(size = 14, color = "black", hjust = 0, family = 'serif'),
      legend.background = element_blank(),
      legend.position = "bottom",
      plot.title = element_text(size = 14, color = "black", face = 'bold', hjust = 0, family = 'serif')
    )
  
  # Build final plot
  p <- ggsurvfit_build(p)
  return(p)
}

# Alternative survival plot function with multi-column legend
# Parameters same as survival_plot4, but configured for more legend items
#------------------------------------------------------------
survival_plot5 <- function(surv, title, legend.labs, palette) {
  # Create ggsurvfit plot with customized elements
  p <- surv %>%
    ggsurvfit(linewidth = 0.8) +
    # Add risk table
    add_risktable(
      risktable_height = 0.2,
      risktable_stats = c("{n.risk}"),
      stats_label = list(n.risk = "Number at risk"),
      size = 4.5,
      plot.margin = margin(r = 0.1, l = 0.1, unit = "cm"),
      theme = list(
        theme_risktable_default(axis.text.y.size = 10, plot.title.size = 12),
        theme(
          plot.title = element_text(face = "bold", size = 12, family = 'serif'),
          text = element_text(family = 'serif'),
          axis.title.x = element_text(face = "bold", family = 'serif'),
          axis.title.y = element_text(family = 'serif'),
          axis.text = element_text(family = 'serif')
        )
      )
    ) +
    # Add symbols and markers
    ggsurvfit::add_risktable_strata_symbol(symbol = "\U25CF", size = 20) +
    add_censor_mark(shape = 1, size = 2, stroke = 1) +
    # Add p-value annotation
    add_pvalue(
      location = "annotation", x = 5, y = 0.1, hjust = 0,
      size = 4, caption = "log-rank {p.value}"
    ) +
    # Set labels
    labs(
      title = title,
      x = "Time from diagnosis (months)",
      y = "Overall Survival"
    ) +
    # Set axis scales
    scale_x_continuous(breaks = seq(0, 60, 12), limits = c(0, 60), expand = c(0.06, 0.06)) +
    scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1), labels = paste0(seq(0, 100, 10), "%"), expand = c(0.03, 0)) +
    # Set colors
    scale_color_manual(values = palette, labels = legend.labs) +
    scale_fill_manual(values = palette, labels = legend.labs) +
    # Apply theme and additional elements
    theme_classic() +
    guides(color = guide_legend(ncol = 2)) +  # Two-column legend
    add_quantile(y_value = 0.5) +
    theme(
      axis.text = element_text(size = 14, color = "black", family = 'serif'),
      axis.title.y = element_text(size = 14, color = "black", vjust = 0.5, family = 'serif'),
      axis.title.x = element_text(size = 14, color = "black", family = 'serif'),
      panel.grid = element_blank(),
      plot.margin = margin(r = 0.1, l = 0.1, unit = "cm"),
      legend.text = element_text(size = 14, color = "black", hjust = 0, family = 'serif'),
      legend.background = element_blank(),
      legend.position = "bottom",
      plot.title = element_text(size = 14, color = "black", face = 'bold', hjust = 0, family = 'serif')
    )
  
  # Build final plot
  p <- ggsurvfit_build(p)
  return(p)
}

#------------------------------------------------------------
# Create survival plots for different age groups and subtypes
#------------------------------------------------------------

# P1: Survival analysis by leukemia subtype across different age groups
#------------------------------------------------------------

# Factor the subtype variable
allage$al_subtype_group1 <- factor(allage$al_subtype_group1, vars_al_type1)

# Plot 1A: All ages by subtype
case1 <- allage
surv <- survfit2(Surv(OS / 30, OS_status) ~ al_subtype_group1, case1)
p11 <- survival_plot4(surv, "All ages, by subtype", vars_al_type1, 
                      palette = paper_col[c(1, 8, 6, 4, 5, 3, 7)])

# Plot 1B: Children by subtype
case1 <- allage %>% filter(group %in% "child")
surv <- survfit2(Surv(OS / 30, OS_status) ~ al_subtype_group1, case1)
p12 <- survival_plot4(surv, "Children, by subtype", vars_al_type1, 
                      palette = paper_col[c(1, 8, 6, 4, 5, 3, 7)])

# Plot 1C: Adults by subtype
case1 <- allage %>% filter(!group %in% "child")
surv <- survfit2(Surv(OS / 30, OS_status) ~ al_subtype_group1, case1)
p13 <- survival_plot4(surv, "Adults, by subtype", vars_al_type1, 
                      palette = paper_col[c(1, 8, 6, 4, 5, 3, 7)])

# P2: Survival analysis by age group within each leukemia subtype
#------------------------------------------------------------
plot_list <- list()
OS_add <- data.frame()

# Generate plots for major leukemia subtypes
for (type in vars_al_type1[1:3]) {
  case1 <- allage %>% filter(al_subtype_group1 == type)
  case1$age_group1 <- as.character(case1$age_group1)
  surv <- survfit2(Surv(OS / 30, OS_status) ~ age_group1 + al_subtype_group1, case1)
  p <- survival_plot4(surv, paste0(type, ", by age group"), 
                      unique(case1$age_group1), 
                      palette = paper_col[c(1, 8, 6, 4, 5, 3, 7)])
  plot_list[[type]] <- p
  median <- OS_select_func(surv)
  OS_add <- rbind(OS_add, median)
}

# Combine plots horizontally
p2 <- wrap_plots(plot_list[1:3], nrow = 1)

# P3: Survival analysis by transplantation status and TKI treatment
#------------------------------------------------------------

# Create age group variable
case$age_group10 <- cut(case$diag_age,
                        breaks = c(9, 19, 29, 39, 49, 59, 69, Inf),
                        labels = c(vars_age_group10[1:6], "70+")
)

# Format transplant variable
case$allo_sct <- as.character(case$allo_sct)
case <- case %>% mutate(allo_sct = case_when(
  allo_sct %in% "1" ~ "Yes",
  allo_sct %in% "0" ~ "No",
  TRUE ~ allo_sct
))
case$allo_sct <- factor(case$allo_sct, c("Yes", "No"))
case <- case %>% arrange(allo_sct)

# Plot 3A: Acute leukemia by allogenic HSCT status
case1 <- case %>% filter(al_subtype_group1 %in% vars_al_type1[c(1, 3)] & 
                           !early_death90 %in% 1 & 
                           age_group10 %in% vars_age_group10[1:6])
surv <- survfit2(Surv(OS / 30, OS_status) ~ allo_sct, case1)
vars <- c("with allo-HSCT", "without allo-HSCT")
p31 <- survival_plot4(surv, "Acute leukaemia, by allo-HSCT", vars, 
                      palette = paper_col[c(1, 8, 2:7)])

# Plot 3B: Acute leukemia with allo-HSCT by age group
case1 <- case %>% filter(al_subtype_group1 %in% vars_al_type1[c(1, 3)] & 
                           !early_death90 %in% 1 & 
                           age_group10 %in% vars_age_group10[1:5] & 
                           allo_sct %in% "Yes")
surv <- survfit2(Surv(OS / 30, OS_status) ~ age_group10, case1)
p32 <- survival_plot4(surv, "Acute leukaemia with allo-HSCT, by age group", 
                      c("15-19", vars_age_group10[2:5]), 
                      palette = paper_col[c(1, 8, 6, 4, 5, 3, 7)])

# Plot 3C: Ph+ ALL by TKI generation and transplantation status
case1 <- case %>% filter(al_subtype_group1 == "ALL")
case1 <- left_join(case1, ALL_select[, c("caseid", "TKI_type", "TKI_type_general", "BCR-ABL1")])
case1 <- case1 %>% filter(!is.na(`BCR-ABL1`))

# Prepare data for TKI analysis
case2 <- case1 %>%
  filter(!TKI_type_general %in% c("Unknown TKI", 'Third-generation TKI')) %>%
  filter(!early_death90 %in% 1)

# Create survival analysis
surv <- survfit2(Surv(OS / 30, OS_status) ~ TKI_type_general + allo_sct, case2)
p33 <- survival_plot5(surv, "Ph+ ALL, by TKIs generation and allo-HSCT",
                      c(
                        "First TKIs,with allo-HSCT",
                        "First TKIs,without allo-HSCT",
                        "Second TKIs,with allo-HSCT",
                        "Second TKIs,without allo-HSCT"
                      ),
                      palette = paper_col[c(1, 8, 6, 4, 5, 3, 7)]
)

# Statistical comparison
surv <- pairwise_survdiff(Surv(OS/30, OS_status) ~ TKI_type_general + allo_sct,
                          p.adjust.method = "bonferroni", data = case2)
surv_p <- as.data.frame(as.matrix(ifelse(round(surv$p.value, 3) < 0.001, "<0.001",
                                         round(surv$p.value, 3))))
surv_p[is.na(surv_p)] <- "-"

#------------------------------------------------------------
# Combine all plots and save
#------------------------------------------------------------

# Arrange all plots into final figure
combined_plot <- ggarrange(
  ggarrange(p11, p12, p13,
            nrow = 1,
            labels = c("A"),
            font.label = list(size = 14, face = "bold")
  ),
  ggarrange(p2,
            nrow = 1,
            labels = c("B"),
            font.label = list(size = 14, face = "bold")
  ),
  ggarrange(p31, p32, p33,
            nrow = 1,
            labels = c("C", 'D', 'E'),
            font.label = list(size = 14, face = "bold")
  ),
  nrow = 3, heights = c(0.5, 0.55, 0.5)
)

# Save final figure in multiple formats
ggsave("bmj_output_changeage/figure2.pdf", 
       combined_plot, width = 26, height = 25, scale = 0.8, device = cairo_pdf)

ggsave("bmj_output_changeage/figure2.svg", 
       combined_plot, width = 60, height = 60, scale = 0.8, units = "cm")

