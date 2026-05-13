# 1. THE FOUNDATION (Data Cleaning & Wrangling)
library(tidyverse) 

# 2. THE CALCULATORS (Converting raw data to Effect Sizes)
library(esc)

# 3. THE ENGINE (The actual Meta-Analysis math)
library(meta)      # If you don't have this, run install.packages("meta")
library(metafor)   # Often used alongside 'meta' for advanced stats

# 4. THE SPECIALISTS (Outliers and specific diagnostics)
library(dmetar)

library(dplyr)

library(grid)

library(gridExtra)

install.packages("metamedian")
library(metamedian)


library(dplyr)

install.packages("estmeansd") # Run this only once
library(estmeansd)


df_final <- complete_raw_data_for_SGLT2i_meta_analysis %>%
  slice(1:60) %>% 
  rowwise() %>%
  mutate(
    # --- INTERVENTION GROUP (E) ---
    needs_calc_e = is.na(mean_e) && !is.na(median_e) && !is.na(iqr_low_e),
    shift_e = if (needs_calc_e && (median_e <= 0 || iqr_low_e <= 0)) 10 else 0,
    
    mean_e = if (needs_calc_e) {
      res_e = estmeansd::bc.mean.sd(q1 = iqr_low_e + shift_e, med = median_e + shift_e, q3 = iqr_high_e + shift_e, n = n_e)
      res_e$est.mean - shift_e
    } else { mean_e },
    
    sd_e = if (needs_calc_e) {
      estmeansd::bc.mean.sd(q1 = iqr_low_e + shift_e, med = median_e + shift_e, q3 = iqr_high_e + shift_e, n = n_e)$est.sd
    } else { sd_e },
    
    # --- CONTROL GROUP (C) ---
    needs_calc_c = is.na(mean_c) && !is.na(median_c) && !is.na(iqr_low_c),
    shift_c = if (needs_calc_c && (median_c <= 0 || iqr_low_c <= 0)) 10 else 0,
    
    mean_c = if (needs_calc_c) {
      res_c = estmeansd::bc.mean.sd(q1 = iqr_low_c + shift_c, med = median_c + shift_c, q3 = iqr_high_c + shift_c, n = n_c)
      res_c$est.mean - shift_c
    } else { mean_c },
    
    sd_c = if (needs_calc_c) {
      estmeansd::bc.mean.sd(q1 = iqr_low_c + shift_c, med = median_c + shift_c, q3 = iqr_high_c + shift_c, n = n_c)$est.sd
    } else { sd_c }
  ) %>%
  ungroup() %>%
  # Removing helper columns
  select(-starts_with("needs_calc"), -starts_with("shift"))


# View your complete dataset
head(df_final)

df_wide <- df_final %>%
  filter(value_type %in% c("Baseline", "Endpoint")) %>%
  select(study_id, outcome, groups, timepoint_label,
         value_type, mean_e, sd_e, n_e) %>%
  pivot_wider(
    names_from = value_type,
    values_from = c(mean_e, sd_e, n_e)
  )

r_value <- 0.7   # recommended for nerve conduction
df_wide <- df_wide %>%
  mutate(
    mean_change_calc = mean_e_Endpoint - mean_e_Baseline,
    sd_change_calc = sqrt(
      sd_e_Baseline^2 +
        sd_e_Endpoint^2 -
        2 * r_value * sd_e_Baseline * sd_e_Endpoint
    )
  )

df_final <- df_final %>%
  left_join(
    df_wide %>%
      select(study_id, outcome, groups, timepoint_label,
             mean_change_calc, sd_change_calc),
    by = c("study_id", "outcome", "groups", "timepoint_label")
  ) %>%
  mutate(
    mean_e = ifelse(
      value_type == "Change" & is.na(mean_e),
      mean_change_calc,
      mean_e
    ),
    sd_e = ifelse(
      value_type == "Change" & is.na(sd_e),
      sd_change_calc,
      sd_e
    )
  ) %>%
  select(-mean_change_calc, -sd_change_calc)


df_wide_c <- df_final %>%
  filter(value_type %in% c("Baseline", "Endpoint")) %>%
  select(study_id, outcome, groups, timepoint_label,
         value_type, mean_c, sd_c, n_c) %>%
  pivot_wider(
    names_from = value_type,
    values_from = c(mean_c, sd_c, n_c)
  )

df_wide_c <- df_wide_c %>%
  mutate(
    mean_change_calc_c = mean_c_Endpoint - mean_c_Baseline,
    sd_change_calc_c = sqrt(
      sd_c_Baseline^2 +
        sd_c_Endpoint^2 -
        2 * r_value * sd_c_Baseline * sd_c_Endpoint
    )
  )

df_final <- df_final %>%
  left_join(
    df_wide_c %>%
      select(study_id, outcome, groups, timepoint_label,
             mean_change_calc_c, sd_change_calc_c),
    by = c("study_id", "outcome", "groups", "timepoint_label")
  ) %>%
  mutate(
    mean_c = ifelse(
      value_type == "Change" & is.na(mean_c),
      mean_change_calc_c,
      mean_c
    ),
    sd_c = ifelse(
      value_type == "Change" & is.na(sd_c),
      sd_change_calc_c,
      sd_c
    )
  ) %>%
  select(-mean_change_calc_c, -sd_change_calc_c)


df_change <- df_final %>%
  filter(value_type == "Change")

outcome_list <- split(df_change, df_change$outcome)

df_change[17, "outcome"] <- "Peroneal_nerve_latency"

run_my_meta <- function(data_subset, fixed, random, inverted = FALSE) {
  outcome_name <- unique(data_subset$outcome)
  
  if (fixed == TRUE) {
    model_label <- "_Fixed"
  } else {
    model_label <- "_Random"
  }
  
  # 3. Add the Inversion label
  if (inverted == TRUE) {
    inv_label <- "_Inverted_"
  } else {
    inv_label <- "_Standard_"
  }
  plot_name <- paste0("Forest", model_label, inv_label, outcome_name)
  
  # 1. Run the meta-analysis
  m <- meta::metacont(n.e = n_e,
                      mean.e = mean_e,
                      sd.e = sd_e,
                      n.c = n_c,
                      mean.c = mean_c,
                      sd.c = sd_c,
                      studlab = study_id,
                      data = data_subset,
                      sm = "MD",
                      fixed = fixed,
                      random = random,
                      method.tau = "REML")
  
  if (inverted == TRUE) { # Invert study-level effects
    old_lower <- m$lower
    old_upper <- m$upper
    
    m$TE    <- -m$TE
    m$lower <- -old_upper
    m$upper <- -old_lower
    
    if (fixed == TRUE && random == FALSE) { 
      # Invert COMMON (fixed) effects pooled estimate
      # Note: m$TE.common is used for the fixed effect result
      old_lower_c <- m$lower.common
      old_upper_c <- m$upper.common
      
      m$TE.common    <- -m$TE.common
      m$lower.common <- -old_upper_c
      m$upper.common <- -old_lower_c}
    
    if (fixed == FALSE && random == TRUE) {# Invert random-effects pooled estimate
      old_lower_r <- m$lower.random
      old_upper_r <- m$upper.random
      
      m$TE.random    <- -m$TE.random
      m$lower.random <- -old_upper_r
      m$upper.random <- -old_lower_r
    }
  }
  
  # 2. Generate the Forest Plot
  meta::forest(m, 
               layout = "RevMan5",
               common = fixed,
               random = random,
               # Forces the "Test for overall effect" for the fixed model
               test.overall.common = fixed,
               test.overall.random = random,
               # Increase horizontal gap between columns (try 4mm or 5mm)
               colgap = "6mm", 
               # Adjust vertical spacing between rows (1.5 is more airy)
               spacing = 1,
               # Ensure the plot doesn't shrink too much
               plotwidth = "8 cm",
               # FORCE Mean and SD specifically
               digits.mean = 1,
               digits.sd = 1,
               
               digits.weight = 1,    
               digits.I2 = 0,        
               digits.pval = 2,
               # Change the column headers for the two groups
               label.e = "SGLT2i",
               label.c = "Standard care",
               label.left = "Standard care",
               label.right = "SGLT2i"
  )
  grid.text(
    plot_name,
    y = unit(0.97, "npc"),
    gp = gpar(fontsize = 14, fontface = "bold")
  )
  
  return(m) # Saves the math results in a list
}

run_my_meta(outcome_list$CNFL, fixed = TRUE, random = FALSE, inverted = FALSE)

leave_one_out_plot <- function(data_set, fixed, random, inverted) {
  meta_plot <- run_my_meta(data_set, fixed, random, inverted)
  leave_one_out <- metainf(meta_plot)
  
  forest(leave_one_out,
         layout = "RevMan5")}

leave_one_out_plot(outcome_list$Sural_nerve_amplitude, fixed = FALSE, random = TRUE, inverted = FALSE)

