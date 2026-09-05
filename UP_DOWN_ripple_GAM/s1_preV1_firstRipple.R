rm(list = ls())

# --- 1. Load Required Libraries ---
packages <- c("mgcv", "gratia", "ggplot2", "dplyr", "tidyr")
for (p in packages) {
  if (!require(p, character.only = TRUE)) install.packages(p)
}

library(mgcv)
library(gratia)
library(ggplot2)
library(dplyr)
library(tidyr)
# --- 2. Load and Prepare Data ---
message("Loading data...")
dat <- read.csv("C:/Users/masah/Documents/GitHub/VR_NPX_analysis/UP_DOWN_ripple_GAM/UP_DOWN_info_GAM.csv")

dat$SessionID <- as.factor(dat$SessionID)
dat$AnimalID <- as.factor(dat$AnimalID)

dat <- dat %>%
  mutate(
    lastRippleNormalisedUP = case_when(
      lastRippleNormalisedUP < 0 ~ 0,
      lastRippleNormalisedUP > 1 ~ 1,
      TRUE ~ lastRippleNormalisedUP
    ),
    lastRippleNormalisedUP_z = as.numeric(scale(lastRippleNormalisedUP))
  )

# 1. Calculate Signed Geometric Mean (Updated names with HPC/V1 suffixes)
dat$geo_coherence <- sign(dat$lastRippleV1 * dat$lastRippleHPC) * 
  sqrt(abs(dat$lastRippleV1 * dat$lastRippleHPC))

dat$geo_coherenceLate <- sign(dat$lateUPV1 * dat$lateUPHPC) * 
  sqrt(abs(dat$lateUPV1 * dat$lateUPHPC))

dat$geo_coherenceRippleLate <- sign(dat$lateUPV1 * dat$lastRippleHPC) * 
  sqrt(abs(dat$lateUPV1 * dat$lastRippleHPC))

# 2. Calculate Signed Geometric Mean for PRE (Updated to use lastRippleHPCPRE)
dat$geo_coherencePRE <- sign(dat$lastRippleV1PRE * dat$lastRippleHPCPRE) * 
  sqrt(abs(dat$lastRippleV1PRE * dat$lastRippleHPCPRE))

# 1. Calculate Signed Geometric Mean (Updated names with HPC/V1 suffixes)
dat$geo_coherence1 <- sign(dat$firstRippleV1 * dat$firstRippleHPC) * 
  sqrt(abs(dat$firstRippleV1 * dat$firstRippleHPC))

# 2. Calculate Signed Geometric Mean for PRE (Updated to use lastRippleHPCPRE)
dat$geo_coherencePRE1 <- sign(dat$firstRippleV1PRE * dat$firstRippleHPC) * 
  sqrt(abs(dat$lastRippleV1PRE * dat$firstRippleHPC))

dat$geo_coherenceRippleNext <- sign(dat$lastRippleHPC * dat$nextUPV1) * 
  sqrt(abs(dat$lastRippleHPC * dat$nextUPV1))

dat$geo_coherenceNext <- sign(dat$lateUPHPC * dat$nextUPV1) * 
  sqrt(abs(dat$lateUPHPC * dat$nextUPV1))

dat$geo_coherenceNextV1 <- sign(dat$lateUPV1 * dat$nextUPV1) * 
  sqrt(abs(dat$lateUPV1 * dat$nextUPV1))


dat_Ripple <- dat %>%
  filter(
    # Keep if the value is <= 1s OR if the value is NaN/NA
    (RippleCounts>0) 
    #(UPDuration_Match >= 0.1    | is.na(UPDuration_Match)) &
    #(UPDuration_NonMatch >= 0.1 | is.na(UPDuration_NonMatch)) 
  )


dat_NoRipple <- dat %>%
  filter(
    # Keep if the value is <= 1s OR if the value is NaN/NA
    (is.na(RippleCounts)) 
    #(UPDuration_Match >= 0.1    | is.na(UPDuration_Match)) &
    #(UPDuration_NonMatch >= 0.1 | is.na(UPDuration_NonMatch)) 
  )

dat_lateUP_Ripple <- dat %>%
  filter(
    # Keep if the value is <= 1s OR if the value is NaN/NA
    #(UPDuration > 0.2    | is.na(UPDuration)) &
    (TimefromLastRipple < 0.1    | is.na(TimefromLastRipple)) 
    #(UPDuration_Match >= 0.1    | is.na(UPDuration_Match)) &
    #(UPDuration_NonMatch >= 0.1 | is.na(UPDuration_NonMatch)) 
  )

dat_lateUP_NoRipple <- dat %>%
  filter(
    # Keep if the value is <= 1s OR if the value is NaN/NA
    #(UPDuration > 0.2    | is.na(UPDuration)) &
    (TimefromLastRipple > 0.1    | is.na(TimefromLastRipple)) 
    #(UPDuration_Match >= 0.1    | is.na(UPDuration_Match)) &
    #(UPDuration_NonMatch >= 0.1 | is.na(UPDuration_NonMatch)) 
  )


dat_Time <- dat %>%
  filter(
    # Keep if the value is <= 1s OR if the value is NaN/NA
    (UPDuration > 0.2    | is.na(UPDuration)) 
    #(UPDuration_Match >= 0.1    | is.na(UPDuration_Match)) &
    #(UPDuration_NonMatch >= 0.1 | is.na(UPDuration_NonMatch)) 
  )


dat_Ripples <- dat %>%
  filter(
    # Keep if the value is <= 1s OR if the value is NaN/NA
    (RippleCounts>1) 
    #(UPDuration_Match >= 0.1    | is.na(UPDuration_Match)) &
    #(UPDuration_NonMatch >= 0.1 | is.na(UPDuration_NonMatch)) 
  )



# 3. Create Z-scores for all numeric columns
dat_clean <- dat %>%
  mutate(across(where(is.numeric), 
                ~ as.numeric(scale(.)), 
                .names = "{.col}_z"))

dat_Time <- dat_Time %>%
  mutate(across(where(is.numeric), 
                ~ as.numeric(scale(.)), 
                .names = "{.col}_z"))

dat_lateUP_Ripple <- dat_lateUP_Ripple %>%
  mutate(across(where(is.numeric), 
                ~ as.numeric(scale(.)), 
                .names = "{.col}_z"))
#
dat_lateUP_NoRipple <- dat_lateUP_NoRipple %>%
  mutate(across(where(is.numeric), 
                ~ as.numeric(scale(.)), 
                .names = "{.col}_z"))

dat_NoRipple <- dat_NoRipple %>%
  mutate(across(where(is.numeric), 
                ~ as.numeric(scale(.)), 
                .names = "{.col}_z"))

dat_Ripple <- dat_Ripple %>%
  mutate(across(where(is.numeric), 
                ~ as.numeric(scale(.)), 
                .names = "{.col}_z"))

dat_Ripples <- dat_Ripples %>%
  mutate(across(where(is.numeric), 
                ~ as.numeric(scale(.)), 
                .names = "{.col}_z"))

# --- 3. Clean Data ---
z_thresh <- 3.5  # include <99.9th centiles of data



my_folder <- "C:/Users/masah/Documents/GitHub/VR_NPX_analysis/UP_DOWN_ripple_GAM/preV1_HC_firstRipples"

# 2. Check if it exists; if not, create it
if (!dir.exists(my_folder)) {
  dir.create(my_folder, recursive = TRUE)
}

# 3. Change the working directory to that folder
setwd(my_folder)


########
######## Pre-ripple V1 predicts first ripple HC content
########

# List of the new Z-score columns to filter (Updated to include all relevant ripple predictors)
z_cols <- c(
  # "lastRippleV1_z", "lastRippleV1PRE_z",
  # "lastRippleHPC_z", "lastRippleHPCPRE_z",
  # "nextUPV1_z"
  "firstRippleV1_z", "firstRippleV1PRE_z",
  "firstRippleHPC_z", "firstRippleHPCPRE_z"
)

dat_clean1 <- dat_clean %>%
  filter(
    if_all(all_of(z_cols), ~ abs(.) < z_thresh | is.na(.))
  )

mdl_V1PRE <- bam(firstRippleHPC_z ~ 
                   # 1. Surviving Main Power Effects
                   # s(firstRippleHPC_z, k = 5) +
                   s(firstRippleV1PRE_z, k = 5) +
                   
                   # 4. Control Term
                   s(AnimalID, bs = "re") +
                   s(SessionID, bs = "re"), 
                 
                 data = dat_clean1, 
                 method = "fREML", 
                 discrete = TRUE, 
                 nthreads = 4)

message("\n--- FINAL MODEL SUMMARY ---")
print(summary(mdl_V1PRE))


### V1_pre
# Calculate scaling factors
raw_breaks <- c(-2,-1,0,1,2)
z_breaks <- sapply(raw_breaks, function(val) {
  dat_clean$firstRippleV1PRE_z[which.min(abs(dat_clean$firstRippleV1PRE - val))]
})

cairo_pdf("pre V1 predicts first ripple HC bias.pdf", width = 4.3, height = 4.3)
p_preV1_raw <- draw(mdl_V1PRE, select = "s(firstRippleV1PRE_z)", residuals = FALSE, rug = FALSE) + 
  theme_bw(base_family = "Arial") + 
  theme(aspect.ratio = 1) +
  scale_x_continuous(breaks = z_breaks, labels = raw_breaks) +
  labs(
    title = "V1 pre-ripple bias predicting first ripple HC bias", 
    x = "pre-ripple V1 bias", 
    y = "Partial Effect"
  )
# dev.new(noRStudioGD = TRUE)

print(p_preV1_raw)
dev.off()




# ==============================================================================
# --- MULTI-METRIC EFFECT SIZE CALCULATIONS ---
# ==============================================================================
message("\nCalculating 4 Effect Size Metrics (This may take a minute)...")

# EXACT formula components to rebuild the models robustly
formula_terms <- c(
  # "s(HC_post_z, k = 5)",
  "s(firstRippleV1PRE_z, k = 5)"
)


# EXACT labels output by summary() and smooth_estimates()
smooth_labels <- c(
  # "s(HC_post_z)",
  "s(firstRippleV1PRE_z)"
)

library(parallel)
library(pbapply)
library(dplyr)
library(stringr)
library(tidyr)

B <- 1000  # Set to 1000 for your final analysis run
RE_TERMS <- c("s(SessionID, bs = 're')", "s(AnimalID, bs = 're')")

message(sprintf("\nLaunching %d Case Bootstrap Replicates...", B))

run_one_bootstrap <- function(rep_id, original_data, formula_terms, smooth_labels) {
  
  # 1. Resample data with replacement
  boot_data <- original_data[sample(nrow(original_data), replace = TRUE), ]
  
  # 2. Fit Full Model
  full_form <- as.formula(paste(
    "firstRippleHPC_z ~",
    paste(c(formula_terms, RE_TERMS), collapse = " + ")
  ))
  
  mdl_f <- mgcv::bam(full_form, data = boot_data, method = "fREML", discrete = TRUE)
  
  f_sum <- summary(mdl_f)
  full_dev <- f_sum$dev.expl * 100
  res_df <- mdl_f$df.residual
  sum_tab <- as.data.frame(f_sum$s.table)
  sum_tab$Term <- rownames(sum_tab)
  
  # Pre-allocate container for this run's metrics
  run_res <- data.frame(Term = smooth_labels, Partial_Deviance = NA, Eta_Sq_Partial = NA, Amplitude = NA, RMS = NA, Rep = rep_id)
  
  # 3. Process every smooth term
  for(j in 1:length(smooth_labels)) {
    # A. Partial Deviance (Drop-One refit)
    act_terms <- c(formula_terms[-j], RE_TERMS)
    r_form <- as.formula(paste("firstRippleHPC_z ~", paste(act_terms, collapse = " + ")))
    mdl_r <- mgcv::bam(r_form, data = boot_data, method = "fREML", discrete = TRUE)
    run_res$Partial_Deviance[j] <- full_dev - (summary(mdl_r)$dev.expl * 100)
    
    # B. Partial Eta-Squared
    t_row <- sum_tab[sum_tab$Term == smooth_labels[j], ]
    if(nrow(t_row) == 1) {
      run_res$Eta_Sq_Partial[j] <- (t_row$F * t_row$edf) / ((t_row$F * t_row$edf) + res_df)
    }
    
    # C & D. Amplitude and RMS via completely self-contained custom grid
    vars_in_smooth <- all.vars(as.formula(paste("~", smooth_labels[j])))
    vars_in_smooth <- vars_in_smooth[!vars_in_smooth %in% c("s", "ti", "bs", "k")]
    
    # Build a clean evaluation sequence for variables in this smooth
    slice_args <- lapply(vars_in_smooth, function(v) {
      seq(min(boot_data[[v]], na.rm = TRUE), max(boot_data[[v]], na.rm = TRUE), length.out = 50)
    })
    names(slice_args) <- vars_in_smooth
    
    # Use base R expand.grid to keep data structure strictly local
    grid_clean <- do.call(expand.grid, slice_args)
    
    # Add dummy/mean columns for any other variables that predict() might structurally look for
    all_model_vars <- all.vars(full_form)[-1] 
    missing_vars <- setdiff(all_model_vars, colnames(grid_clean))
    
    for(mv in missing_vars) {
      if(mv == "SessionID") {
        grid_clean[[mv]] <- boot_data$SessionID[1] 
      } else if(is.numeric(boot_data[[mv]])) {
        grid_clean[[mv]] <- mean(boot_data[[mv]], na.rm = TRUE)
      } else {
        grid_clean[[mv]] <- boot_data[[mv]][1]
      }
    }
    
    # Isolate target smooth prediction using the lpmatrix
    Xp <- predict(mdl_f, newdata = grid_clean, type = "lpmatrix")
    smooth_cols <- grep(smooth_labels[j], colnames(Xp), fixed = TRUE)
    
    Xp_isolated <- matrix(0, nrow = nrow(Xp), ncol = ncol(Xp))
    Xp_isolated[, smooth_cols] <- Xp[, smooth_cols]
    
    fit_isolated <- Xp_isolated %*% coef(mdl_f)
    run_res$Amplitude[j] <- max(fit_isolated) - min(fit_isolated)
    run_res$RMS[j]       <- sqrt(mean(fit_isolated^2))
  }
  
  return(run_res)
}

# Run loop with progress bar
set.seed(42)
boot_results_list <- pblapply(1:B, run_one_bootstrap, 
                              original_data = dat_clean1, 
                              formula_terms = formula_terms, 
                              smooth_labels = smooth_labels)

# Clean out any NULLs if any structural issues happen
boot_results_list <- boot_results_list[!sapply(boot_results_list, is.null)]

# Combine results safely
boot_df <- do.call(rbind, boot_results_list)

# 4. Compute Median and 95% CIs from the Bootstrap Distribution
final_dashboard_data <- boot_df %>%
  dplyr::group_by(Term) %>%
  dplyr::summarise(
    dplyr::across(c(Partial_Deviance, Eta_Sq_Partial, Amplitude, RMS),
                  list(Val = ~median(.x, na.rm = TRUE),
                       Lwr = ~quantile(.x, probs = 0.025, na.rm = TRUE),
                       Upr = ~quantile(.x, probs = 0.975, na.rm = TRUE)),
                  .names = "{.col}__{.fn}")
  )

# ==============================================================================
# --- 3E. SAVE RAW BOOTSTRAP DISTRIBUTIONS (ALL REPLICATES) ---
# ==============================================================================
message("Saving raw bootstrap iterations for distribution archives...")

# Format the column headers cleanly before exporting
raw_iterations_clean <- boot_df %>%
  select(Rep, Term, Partial_Deviance, Eta_Sq_Partial, Amplitude, RMS) %>%
  rename(
    Replicate        = Rep,
    Deviance_Value   = Partial_Deviance,
    Eta_Sq_Value     = Eta_Sq_Partial,
    Amplitude_Value  = Amplitude,
    RMS_Value        = RMS
  ) %>%
  arrange(Term, Replicate)

# Export the raw simulation matrix to a CSV file
write.csv(raw_iterations_clean, "GAM_model_raw_bootstrap_iterations.csv", row.names = FALSE)

message("Raw bootstrap iteration file saved successfully: 'GAM_model_raw_bootstrap_iterations.csv'")


library(tidyverse)

# ==============================================================================
# --- 4. GRAPH BOOTSTRAPPED EFFECT SIZE DASHBOARD (BAR PLOT EDITION) ---
# ==============================================================================

# 1. Pivot the Median Estimates
eff_val <- final_dashboard_data %>%
  select(Term, ends_with("__Val")) %>%
  rename_with(~str_remove(., "__Val"), -Term) %>%
  pivot_longer(cols = -Term, names_to = "Metric", values_to = "Value")

# 2. Pivot the Lower Confidence Bounds
eff_lwr <- final_dashboard_data %>%
  select(Term, ends_with("__Lwr")) %>%
  rename_with(~str_remove(., "__Lwr"), -Term) %>%
  pivot_longer(cols = -Term, names_to = "Metric", values_to = "Lower")

# 3. Pivot the Upper Confidence Bounds
eff_upr <- final_dashboard_data %>%
  select(Term, ends_with("__Upr")) %>%
  rename_with(~str_remove(., "__Upr"), -Term) %>%
  pivot_longer(cols = -Term, names_to = "Metric", values_to = "Upper")

# 4. Merge everything cleanly and apply factor labels
plot_data <- eff_val %>%
  left_join(eff_lwr, by = c("Term", "Metric")) %>%
  left_join(eff_upr, by = c("Term", "Metric")) %>%
  mutate(Metric = factor(Metric, 
                         levels = c("Partial_Deviance", "Eta_Sq_Partial", "Amplitude", "RMS"),
                         labels = c("Deviance Explained (%)", "Partial Eta-Squared", 
                                    "Peak-to-Trough Amplitude", "RMS Effect")))

# 5. Construct the Custom Bar Plot
p_bars_with_ci <- ggplot(plot_data, aes(x = reorder(Term, Value), y = Value, fill = Term)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60", alpha = 0.5) +
  geom_bar(stat = "identity", show.legend = FALSE, alpha = 0.85, width = 0.75) +
  # Asymmetric error bars marking the lower and upper bounds cleanly
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.25, color = "black", size = 0.6) +
  coord_flip() +
  facet_wrap(~Metric, scales = "free_x", ncol = 2) +
  scale_fill_viridis_d(option = "mako", direction = -1) +
  theme_bw() +
  labs(title = "GAMM Effect Size & Variance Metrics",
       subtitle = "Bars denote bootstrap medians; error bars represent 95% non-parametric CIs.",
       x = NULL, y = NULL) +
  theme(
    strip.text = element_text(face = "bold", size = 11, colour = "black"),
    strip.background = element_rect(fill = "gray95", color = "gray80"),
    axis.text.y = element_text(face = "plain", size = 10, color = "black"),
    axis.text.x = element_text(size = 9, color = "black"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "gray95"),
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10, color = "gray40"),
    panel.spacing = unit(1.2, "lines")
  )

# Render to RStudio console
print(p_bars_with_ci)

# 6. Save vector graphic output
cairo_pdf("Model_Effect_Sizes_With_CI.pdf", width = 9.5, height = 6.5)
print(p_bars_with_ci)
dev.off()



# ==============================================================================
# --- 5. SAVE COMPREHENSIVE MODEL OUTPUT WITH 95% CIs ---
# ==============================================================================
library(dplyr)
library(tidyr)
library(stringr)

# 1. Format the basic model summary parameters
model_stats <- as.data.frame(summary(mdl_V1PRE)$s.table) %>%
  mutate(Term = rownames(.))

# 2. Flatten the nested stats into individual columns
flat_bootstrap_results <- final_dashboard_data %>%
  pivot_longer(cols = -Term, names_to = "Combined", values_to = "Value") %>%
  separate(Combined, into = c("Metric", "Stat"), sep = "__") %>%
  mutate(New_Col_Name = paste0(Metric, "_", Stat)) %>%
  select(-Metric, -Stat) %>%
  pivot_wider(names_from = New_Col_Name, values_from = Value)

# 3. Join the parametric summary table with the resampled metrics
combined_results <- model_stats %>%
  left_join(flat_bootstrap_results, by = "Term") %>%
  select(
    Term, edf, Ref.df, F, `p-value`,
    Partial_Deviance_Val, Partial_Deviance_Lwr, Partial_Deviance_Upr,
    Eta_Sq_Partial_Val, Eta_Sq_Partial_Lwr, Eta_Sq_Partial_Upr,
    Amplitude_Val, Amplitude_Lwr, Amplitude_Upr,
    RMS_Val, RMS_Lwr, RMS_Upr
  ) %>%
  rename(
    p_value              = `p-value`,
    Deviance_Median      = Partial_Deviance_Val,
    Deviance_CI_Lower    = Partial_Deviance_Lwr,
    Deviance_CI_Upper    = Partial_Deviance_Upr,
    Eta_Sq_Median        = Eta_Sq_Partial_Val,
    Eta_Sq_CI_Lower      = Eta_Sq_Partial_Lwr,
    Eta_Sq_CI_Upper      = Eta_Sq_Partial_Upr,
    Amplitude_Median     = Amplitude_Val,
    Amplitude_CI_Lower   = Amplitude_Lwr,
    Amplitude_CI_Upper   = Amplitude_Upr,
    RMS_Median           = RMS_Val,
    RMS_CI_Lower         = RMS_Lwr,
    RMS_CI_Upper         = RMS_Upr
  )

# 4. Export the comprehensive data table
write.csv(combined_results, "GAM_model_CI_output.csv", row.names = FALSE)

message("Master spreadsheet file saved successfully: 'GAM_model_CI_output.csv'")
message("\nAnalysis Pipeline Complete!")


