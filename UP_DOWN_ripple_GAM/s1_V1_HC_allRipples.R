rm(list = ls())

# ==============================================================================
# pre-ripple V1 + ripple HC -> ripple V1
# ==============================================================================

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
#dat <- read.csv("C:/Users/masah/Documents/GitHub/VR_NPX_analysis/UP_DOWN_ripple_GAM/pre_post_normalised_UP_ripple_20_120ms.csv")
#dat <- read.csv("C:/Users/masah/Documents/GitHub/VR_NPX_analysis/UP_DOWN_ripple_GAM/pre_post_normalised_UP_ripple_-150_-50ms.csv")
dat <- read.csv("C:/Users/masah/Documents/GitHub/VR_NPX_analysis/UP_DOWN_ripple_GAM/pre_post_normalised_UP_ripple.csv")

dat$SessionID <- as.factor(dat$SessionID)
dat$AnimalID <- as.factor(dat$AnimalID)

dat$geo_coherence <- sign(dat$V1_post * dat$HC_post) * 
  sqrt(abs(dat$V1_post * dat$HC_post))

# 2. Calculate Signed Geometric Mean for V1-PRE and HPC
dat$geo_coherencePRE <- sign(dat$V1_pre * dat$HC_post) * 
  sqrt(abs(dat$V1_pre * dat$HC_post))

dat_clean <- dat %>%
  mutate(across(where(is.numeric), 
                ~ as.numeric(scale(.)), 
                .names = "{.col}_z"))

# Check the results
head(dat_clean)

# Verify the new columns exist
colnames(dat_clean)
# 
# message(sprintf("Trimming extreme outliers beyond +/- %s Z-scores for %d variables...", z_thresh, length(z_cols)))
# 
# --- 3. Clean Data ---
# z_thresh <- 2.6  # include <99th centiles of data
z_thresh <- 3.5  # include <99th centiles of data

# List of the new Z-score columns to filter
z_cols <- c(
  "V1_pre_z", "V1_post_z", "HC_post_z" , "HC_pre_z",
  "geo_coherencePRE_z","geo_coherence_z"
)
# # List of the new Z-score columns to filter
# z_cols <- c(
#   "V1_pre_z", "V1_post_z", "HC_post_z" , "HC_pre_z",
#   "ripple_power_z"
# )
message(sprintf("Trimming extreme outliers beyond +/- %s Z-scores for %d variables...", z_thresh, length(z_cols)))

dat_clean <- dat_clean %>%
  filter(
    (NormalisedUP_L >=0    | NormalisedUP_R >=0),
    if_all(all_of(z_cols), ~ abs(.) < z_thresh | is.na(.))
  )
# dat_clean <- dat_clean %>%
#   filter(
#     # (NormalisedUP_L >=0    | NormalisedUP_R >=0), 
#     if_all(all_of(z_cols), ~ abs(.) < z_thresh | is.na(.))
#   )


my_folder <- "C:/Users/masah/Documents/GitHub/VR_NPX_analysis/UP_DOWN_ripple_GAM/V1_HC_allRipples"

# 2. Check if it exists; if not, create it
if (!dir.exists(my_folder)) {
  dir.create(my_folder, recursive = TRUE)
}

# 3. Change the working directory to that folder
setwd(my_folder)

# 
# # Fit the model where HC bias is predicted by the interaction of Pre-V1 and Time
# mdl_pre_post <- bam(geo_coherence_z ~ 
#                       s(ripple_power_z,k = 5) +
#                       # s(V1_pre_z, k = 5) +
#                       # ti(ripple_power_z, HC_post_z, k = c(5,5)) +
#                       # te(V1_pre_z, HC_post_z, k = c(5,5)) +
#                       s(AnimalID, bs = "re")+
#                       s(SessionID, bs = "re"), 
#                     data = dat_clean, method = "fREML", discrete = TRUE)
# 
# print(summary(mdl_pre_post))
# 
# 
# ### last ripple power -> last ripple coherence
# # Calculate scaling factors
# raw_breaks <- c(5,10,15,20)
# z_breaks <- sapply(raw_breaks, function(val) {
#   dat_clean$ripple_power_z[which.min(abs(dat_clean$ripple_power - val))]
# })
# 
# # cairo_pdf("lateUPHPC_nextUP", width = 4.3, height = 4.3)
# p_lag_raw <- draw(mdl_pre_post, select = "s(ripple_power_z)", residuals = FALSE, rug = FALSE) + 
#   theme_bw(base_family = "Arial") + 
#   theme(aspect.ratio = 1) +
#   scale_x_continuous(breaks = z_breaks, labels = raw_breaks) +
#   coord_cartesian(xlim = c(-1.2,3)) + 
#   labs(
#     title = "Ripple power on ripple coherence", 
#     x = "Ripple power", 
#     y = "Partial Effect"
#   )
# dev.new(noRStudioGD = TRUE)
# print(p_lag_raw)

#####
library(mgcv)
library(ggplot2)
library(dplyr)
# Fit the model where HC bias is predicted by the interaction of Pre-V1 and Time
mdl_pre_post <- bam(V1_post_z ~ 
                      s(HC_post_z,k = 5) +
                      s(V1_pre_z, k = 5) +
                      # ti(V1_pre_z, HC_post_z, k = c(5,5)) +
                      # te(V1_pre_z, HC_post_z, k = c(5,5)) +
                      s(AnimalID, bs = "re")+
                      s(SessionID, bs = "re"), 
                    data = dat_clean, method = "fREML", discrete = TRUE)

print(summary(mdl_pre_post))

# ---------------------------------------------------------
# Plot 1: Main Effect of HC_post_z
# ---------------------------------------------------------
### HC post
# Calculate scaling factors
cairo_pdf("HC ripple bias.pdf", width = 4.3, height = 4.3)
raw_breaks <- c(-2,-1,0,1,2)
z_breaks <- sapply(raw_breaks, function(val) {
  dat_clean$HC_post_z[which.min(abs(dat_clean$HC_post - val))]
})

# cairo_pdf("lateUPHPC_nextUP", width = 4.3, height = 4.3)
p_HC_raw <- draw(mdl_pre_post, select = "s(HC_post_z)", residuals = FALSE, rug = FALSE) + 
  theme_bw(base_family = "Arial") + 
  theme(aspect.ratio = 1) +
  scale_x_continuous(breaks = z_breaks, labels = raw_breaks) +
  labs(
    title = "HC ripple bias predicting post V1 track bias", 
    x = "HC bias", 
    y = "Partial Effect"
  )
# dev.new(noRStudioGD = TRUE)
print(p_HC_raw)
dev.off()

# ---------------------------------------------------------
# Plot 2: Main Effect of V1_pre_z
# ---------------------------------------------------------
# sm_v1 <- smooth_estimates(mdl_pre_post, smooth = "s(V1_pre_z)", n = 200) %>%
#   mutate(
#     .lower_ci = .estimate - (1.96 * .se),
#     .upper_ci = .estimate + (1.96 * .se)
#   )
# 
# p_v1 <- ggplot(sm_v1, aes(x = V1_pre_z, y = .estimate)) +
#   geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), alpha = 0.2, fill = "firebrick") +
#   geom_line(linewidth = 1, color = "firebrick") +
#   theme_bw() + 
#   labs(title = "Isolated Main Effect: V1 Pre (Z)", 
#        x = "V1_pre_z", y = "Partial Effect")
# 
# print(p_v1)


### V1_pre
# Calculate scaling factors
cairo_pdf("V1 pre ripple bias.pdf", width = 4.3, height = 4.3)
raw_breaks <- c(-2,-1,0,1,2)
z_breaks <- sapply(raw_breaks, function(val) {
  dat_clean$V1_pre_z[which.min(abs(dat_clean$V1_pre - val))]
})

# cairo_pdf("lateUPHPC_nextUP", width = 4.3, height = 4.3)
p_preV1_raw <- draw(mdl_pre_post, select = "s(V1_pre_z)", residuals = FALSE, rug = FALSE) + 
  theme_bw(base_family = "Arial") + 
  theme(aspect.ratio = 1) +
  scale_x_continuous(breaks = z_breaks, labels = raw_breaks) +
  labs(
    title = "V1 pre-ripple bias predicting post V1 bias", 
    x = "pre-ripple V1 bias", 
    y = "Partial Effect"
  )
# dev.new(noRStudioGD = TRUE)

print(p_preV1_raw)
dev.off()
# 
# # ---------------------------------------------------------
# # Plot 3: Pure Interaction Tensor (ti) 
# # ---------------------------------------------------------
# sm_2d <- smooth_estimates(mdl_pre_post, smooth = "ti(V1_pre_z,HC_post_z)", n = 100)
# 
# p_interaction <- ggplot(sm_2d, aes(x = V1_pre_z, y = HC_post_z, fill = .estimate)) +
#   geom_tile() + 
#   geom_contour(aes(z = .estimate), color = "black", alpha = 0.2) + 
#   scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "firebrick", midpoint = 0, name = "Effect") +
#   theme_minimal() +
#   theme(aspect.ratio = 1) +
#   labs(title = "Pure Tensor Interaction (ti Term)",
#        subtitle = "Variance unique to the combination",
#        x = "V1_pre_z", y = "HC_post_z")
# 
# print(p_interaction)

# 
# # ==============================================================================
# # --- RECREATING THE TOTAL LANDSCAPE ---
# # ==============================================================================
# 
# 1. Create a clean grid across the observed range of your data
# grid_range_v1 <- seq(min(dat_clean$V1_pre_z, na.rm=TRUE), max(dat_clean$V1_pre_z, na.rm=TRUE), length.out = 100)
# grid_range_hc <- seq(min(dat_clean$HC_post_z, na.rm=TRUE), max(dat_clean$HC_post_z, na.rm=TRUE), length.out = 100)
cairo_pdf("V1 pre and HC combined.pdf", width = 4.3, height = 3)

grid_range_v1 <- seq(-2.5,2.5, length.out = 50)
grid_range_hc <- seq(-2.5,2.5, length.out = 50)

pred_grid_total <- expand.grid(
  V1_pre_z  = grid_range_v1,
  HC_post_z = grid_range_hc,
  AnimalID  = dat_clean$AnimalID[1], # Held constant (ignored by type="terms" if excluded, but safe)
  SessionID  = dat_clean$SessionID[1]
)

# 2. Extract specific term components matrix
term_preds <- predict(mdl_pre_post, newdata = pred_grid_total, type = "terms")

# 3. Sum only your targets of interest
pred_grid_total$Reconstructed_Effect <-
  term_preds[, "s(HC_post_z)"] +
  term_preds[, "s(V1_pre_z)"]
  # term_preds[, "ti(V1_pre_z,HC_post_z)"]

# 4. Plot Combined Surface
p_combined <- ggplot(pred_grid_total, aes(x = V1_pre_z, y = HC_post_z, fill = Reconstructed_Effect)) +
  geom_tile() +
  geom_contour(aes(z = Reconstructed_Effect), color = "black", alpha = 0.2) +
  scale_fill_gradient2(high = "blue", mid = "white", low = "firebrick", midpoint = 0, name = "Total\nEffect") +
  theme_minimal() +
  # theme(aspect.ratio = 1) +
  labs(
    title = "Total Combined Effect Surface",
    subtitle = "Reconstructed: s(HC) + s(V1)",
    x = "V1_pre_z",
    y = "HC_post_z"
  )
# dev.new(noRStudioGD = TRUE)

print(p_combined)
dev.off()



# ==============================================================================
# --- MULTI-METRIC EFFECT SIZE CALCULATIONS ---
# ==============================================================================
message("\nCalculating 4 Effect Size Metrics (This may take a minute)...")

# EXACT formula components to rebuild the models robustly
formula_terms <- c(
  "s(HC_post_z, k = 5)",
  "s(V1_pre_z, k = 5)"
)


# EXACT labels output by summary() and smooth_estimates()
smooth_labels <- c(
  "s(HC_post_z)",
  "s(V1_pre_z)"
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
    "V1_post_z ~",
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
    r_form <- as.formula(paste("V1_post_z ~", paste(act_terms, collapse = " + ")))
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
                              original_data = dat_clean, 
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
model_stats <- as.data.frame(summary(mdl_pre_post)$s.table) %>%
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





