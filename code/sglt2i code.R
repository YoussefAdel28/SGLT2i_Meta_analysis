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

install.packages("estmeansd") # Run this only once
library(estmeansd)

# this section created df_final which calculated mean and sd from median and IQR range
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


df_wide <- df_final %>%
  filter(value_type %in% c("Baseline", "Endpoint")) %>%
  select(study_id, outcome, groups, timepoint_label,
         value_type, mean_e, sd_e, n_e) %>%
  pivot_wider(
    names_from = value_type,
    values_from = c(mean_e, sd_e, n_e)
  )

r_value <- 0   # recommended for nerve conduction
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

df_change[17, "outcome"] <- "Peroneal_nerve_latency"

outcome_list <- split(df_change, df_change$outcome)
outcome_list_rightlimb <- split(df_change_r_0_right_limb_haggar_data, df_change_r_0_right_limb_haggar_data$outcome)
outcome_list_leftlimb <- split(df_change_r_0_left_limb_haggar_data, df_change_r_0_left_limb_haggar_data$outcome)

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
               label.c = "Control",
               label.left = "Control",
               label.right = "SGLT2i"
  )
  grid.text(
    plot_name,
    y = unit(0.97, "npc"),
    gp = gpar(fontsize = 14, fontface = "bold")
  )
  
  return(m) # Saves the math results in a list
}


leave_one_out_plot <- function(data_set, fixed, random, inverted) {
  meta_plot <- run_my_meta(data_set, fixed, random, inverted)
  leave_one_out <- metainf(meta_plot)
  
  forest(leave_one_out,
         layout = "RevMan5")}


run_and_save_meta_tiff <- function(data_subset, fixed, random, inverted = FALSE) {
  outcome_name <- unique(data_subset$outcome)
  
  # Create a filename (e.g., "Forest_fixed/random_CMS_6wks.tif")
  
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
  
  filename <- paste0("Forest", model_label, inv_label, outcome_name, ".tif") 
  
  
  tiff(filename, 
       width = 12, height = 8, 
       units = "in", 
       res = 600,            # 600 DPI is superior to the standard 300
       compression = "lzw")  # Lossless compression
  
  
  run_my_meta(data_subset, fixed, random, inverted)
  
  dev.off()
}

df_change_r_0_right_without_ishibashi <- df_change_r_0_right_limb_haggar_data %>% 
  filter(study_id != "Ishibashi 2022")
outcome_list_rightlimb_without_ishibashi <- split(df_change_r_0_right_without_ishibashi, df_change_r_0_right_without_ishibashi$outcome)

run_my_meta(outcome_list_rightlimb_without_ishibashi$Sural_nerve_velocity, fixed = TRUE, random = FALSE, inverted = FALSE)
leave_one_out_plot(outcome_list_leftlimb$Sural_nerve_amplitude, fixed = FALSE, random = TRUE, inverted = FALSE)
run_and_save_meta_tiff(outcome_list_rightlimb$Sural_nerve_velocity, fixed = TRUE, random = FALSE, inverted = FALSE)

library(grid)
library(gridExtra)


# 7. CORE FUNCTIONS
# ================================================================

# ----------------------------------------------------------------
# 7A. FOREST PLOT GROB GENERATOR
# Runs meta-analysis and captures forest plot as a grob object
# plot_label is now optional (NULL = no label handling here)
# ----------------------------------------------------------------
create_forest_grob <- function(data_subset, fixed, random,
                               inverted = FALSE) {
  grob <- grid::grid.grabExpr({
    
    # Tight margins to eliminate captured whitespace
    par(mar = c(0.1, 0.1, 0.1, 0.1))
    
    # Run meta-analysis
    m <- meta::metacont(
      n.e        = n_e,
      mean.e     = mean_e,
      sd.e       = sd_e,
      n.c        = n_c,
      mean.c     = mean_c,
      sd.c       = sd_c,
      studlab    = study_id,
      data       = data_subset,
      sm         = "MD",
      fixed      = fixed,
      random     = random,
      method.tau = "REML"
    )
    
    # Apply inversion if outcome direction needs flipping
    # (e.g. latency where reduction = improvement)
    if (inverted) {
      old_lower <- m$lower
      old_upper <- m$upper
      m$TE      <- -m$TE
      m$lower   <- -old_upper
      m$upper   <- -old_lower
      
      if (fixed && !random) {
        old_lc         <- m$lower.common
        old_uc         <- m$upper.common
        m$TE.common    <- -m$TE.common
        m$lower.common <- -old_uc
        m$upper.common <- -old_lc
      }
      if (!fixed && random) {
        old_lr         <- m$lower.random
        old_ur         <- m$upper.random
        m$TE.random    <- -m$TE.random
        m$lower.random <- -old_ur
        m$upper.random <- -old_lr
      }
    }
    
    # Draw forest plot
    meta::forest(
      m,
      layout              = "RevMan5",
      common              = fixed,
      random              = random,
      test.overall.common = fixed,
      test.overall.random = random,
      colgap              = "6mm",
      spacing             = 1,
      plotwidth           = "8cm",
      digits.mean         = 1,
      digits.sd           = 1,
      digits.weight       = 1,
      digits.I2           = 0,
      digits.pval         = 2,
      label.e             = "SGLT2i",
      label.c             = "Control",
      label.left          = "Control",
      label.right         = "SGLT2i"
    )
    
  }, wrap = TRUE)
  
  return(grob)
}

# ----------------------------------------------------------------
# 7B. LABEL STRIP GENERATOR
# Creates a thin grob containing just the bold outcome title
# left_offset aligns the title with the Study column of the plot
# Adjust left_offset if title is not aligned with study names:
#   too far left  → increase (e.g. 0.10, 0.12)
#   too far right → decrease (e.g. 0.06, 0.04)
# ----------------------------------------------------------------
make_label_strip <- function(label_text, left_offset = 0.085, font_size = 11) {
  
  grid::grid.grabExpr({
    grid.rect(gp = gpar(fill = "white", col = NA))
    grid.text(
      label = label_text,
      x     = unit(left_offset, "npc"),
      y     = unit(0.5, "npc"),
      just  = c("left", "center"),
      gp    = gpar(
        fontsize  = font_size,   # <-- controlled by argument now
        fontface  = "bold",
        col       = "black"
      )
    )
  }, wrap = TRUE)
}

# ----------------------------------------------------------------
# 7C. COMBINED FIGURE SAVER
# Accepts a named list where:
#   name  = label text for that outcome
#   value = forest plot grob for that outcome
# Automatically interleaves label strips and forest plots
# Heights: adjust strip_height to control label strip thickness
#          adjust plot_height to control forest plot height
# ----------------------------------------------------------------
save_combined_forest_tiff <- function(named_grob_list,
                                      filename,
                                      width        = 14,
                                      height       = NULL,
                                      res          = 600,
                                      strip_height = 0.25,
                                      plot_height  = 4.0,
                                      left_offset  = 0.085,
                                      font_size    = 11) {  # <-- added here
  
  n_plots <- length(named_grob_list)
  
  if (is.null(height)) {
    height <- n_plots * (plot_height + strip_height) + 0.5
  }
  
  all_grobs   <- list()
  all_heights <- c()
  
  for (i in seq_along(named_grob_list)) {
    
    label_text  <- names(named_grob_list)[i]
    forest_grob <- named_grob_list[[i]]
    
    label_strip <- make_label_strip(
      label_text,
      left_offset = left_offset,
      font_size   = font_size    # <-- passed through here
    )
    
    all_grobs   <- c(all_grobs, list(label_strip), list(forest_grob))
    all_heights <- c(all_heights, strip_height, plot_height)
  }
  
  tiff(
    filename,
    width       = width,
    height      = height,
    units       = "in",
    res         = res,
    compression = "lzw"
  )
  
  gridExtra::grid.arrange(
    grobs   = all_grobs,
    ncol    = 1,
    heights = all_heights
  )
  
  dev.off()
  message("Saved successfully: ", filename)
}

# ----------------------------------------------------------------

# ================================================================
# 8. GENERATE ALL FOREST PLOT GROBS
# ================================================================

# --- SURAL NERVE ---
grob_sural_velocity <- create_forest_grob(
  data_subset = outcome_list_rightlimb$Sural_nerve_velocity,
  fixed       = TRUE,
  random      = FALSE,
  inverted    = FALSE
)

grob_sural_amplitude <- create_forest_grob(
  data_subset = outcome_list_rightlimb$Sural_nerve_amplitude,
  fixed       = FALSE,
  random      = TRUE,
  inverted    = FALSE
)

# --- PERONEAL NERVE ---
grob_peroneal_latency <- create_forest_grob(
  data_subset = outcome_list_rightlimb$Peroneal_nerve_latency,
  fixed       = TRUE,
  random      = FALSE,
  inverted    = TRUE    # latency: reduction = improvement so inverted
)

grob_peroneal_velocity <- create_forest_grob(
  data_subset = outcome_list_rightlimb$Peroneal_nerve_velocity,
  fixed       = FALSE,
  random      = TRUE,
  inverted    = FALSE
)

grob_peroneal_amplitude <- create_forest_grob(
  data_subset = outcome_list_rightlimb$Peroneal_nerve_amplitude,
  fixed       = FALSE,
  random      = TRUE,
  inverted    = FALSE
)

# --- MDA ---
grob_mda <- create_forest_grob(
  data_subset = outcome_list$MDA,
  fixed       = FALSE,
  random      = TRUE,
  inverted    = TRUE
)

# ================================================================
# 9. SAVE COMBINED FIGURES
# ================================================================

# --- FIGURE 1: Sural Nerve (velocity + amplitude combined) ---
save_combined_forest_tiff(
  named_grob_list = list(
    "A. Sural Nerve Sensory Conduction Velocity" = grob_sural_velocity,
    "B. Sural Nerve Sensory Amplitude"           = grob_sural_amplitude
  ),
  filename     = "Figure1_Sural_Nerve.tif",
  width        = 14,
  res          = 600,
  strip_height = 0.35,   # increase strip height slightly when font is larger
  plot_height  = 2.3,
  left_offset  = 0.10,
  font_size    = 14      # <-- change this number to whatever size you want
)

# --- FIGURE 2: Peroneal Nerve (latency + velocity + amplitude) ---
save_combined_forest_tiff(
  named_grob_list = list(
    "A. Peroneal Nerve Motor Conduction Velocity" = grob_peroneal_velocity,
    "B. Peroneal Nerve Motor Amplitude"           = grob_peroneal_amplitude,
    "C. Peroneal Nerve Motor Latency"             = grob_peroneal_latency
  
  ),
  filename     = "Figure2_Peroneal_Nerve.tif",
  width        = 14,
  res          = 600,
  strip_height = 0.35,
  plot_height  = 2.3,
  left_offset  = 0.10,
  font_size    = 14
)

# --- FIGURE 3: MDA (standalone) ---
save_combined_forest_tiff(
  named_grob_list = list(
    "Malondialdehyde (MDA)" = grob_mda
  ),
  filename     = "Figure3_MDA.tif",
  width        = 14,
  res          = 600,
  strip_height = 0.35,
  plot_height  = 2,
  left_offset  = 0.10,
  font_size    = 14
)

# ================================================================
# CORNEAL CONFOCAL MICROSCOPY FOREST PLOTS
# ================================================================

# ----------------------------------------------------------------
# Mapping: outcome name → full title
# ----------------------------------------------------------------
corneal_outcomes <- list(
  CNFD = "A. Corneal Nerve Fibre Density (CNFD)",
  CNFL = "B. Corneal Nerve Fibre Length (CNFL)",
  CNBD = "C. Corneal Nerve Branch Density (CNBD)"
)

# ----------------------------------------------------------------
# Mapping: outcome name → model type
# TRUE = fixed effects, FALSE = random effects
# ----------------------------------------------------------------
fixed_effects_model <- c(
  CNFD = FALSE,
  CNFL = TRUE,    # <-- only CNFL uses fixed effects
  CNBD = FALSE
)

# ----------------------------------------------------------------
# Loop: generate one forest grob per corneal outcome
# ----------------------------------------------------------------
corneal_grobs <- list()

for (acronym in names(corneal_outcomes)) {
  
  full_title  <- corneal_outcomes[[acronym]]
  use_fixed   <- fixed_effects_model[[acronym]]
  use_random  <- !use_fixed
  
  if (!acronym %in% names(outcome_list_rightlimb)) {
    warning(paste("Outcome not found in dataset:", acronym, "— skipping"))
    next
  }
  
  corneal_grobs[[full_title]] <- create_forest_grob(
    data_subset = outcome_list_rightlimb[[acronym]],
    fixed       = use_fixed,
    random      = use_random,
    inverted    = FALSE
  )
}

# ----------------------------------------------------------------
# Save combined figure
# ----------------------------------------------------------------
save_combined_forest_tiff(
  named_grob_list = corneal_grobs,
  filename        = "Figure_Corneal_Outcomes.tif",
  width           = 14,
  res             = 600,
  strip_height    = 0.35,
  plot_height     = 2.3,
  left_offset     = 0.10,
  font_size       = 14
)


library(robvis)

rob2 <- rob_traffic_light(ROB2_raw_results, "ROB2")


library(ggplot2)

rob2_final <- rob2 + 
  theme(strip.text.y.left=element_text(angle = 0))

ggsave(
  filename    = "C:/Users/youss/Documents/2_Academics_Career/3_INTERNSHIP/7_Research/SR&MA/SGLT2i and peripheral neuropathy/Meta-analysis/SGLT2i_Meta_analysis/RoB2_Traffic_Light_Final.tif",
  plot        = rob2_final,
  width       = 10,
  height      = 7,    # slightly taller to accommodate horizontal labels
  dpi         = 600,
  device      = "tiff",
  compression = "lzw"
)