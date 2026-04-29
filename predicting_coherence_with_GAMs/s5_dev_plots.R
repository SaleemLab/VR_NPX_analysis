# ==============================================================================
# V1-HC COHERENCE: STAGE 5 - FINAL MODEL & DIAGNOSTICS
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
dat <- read.csv("D:/Code/MasaCoherenceDataAnalysis/v1_hc_data.csv")
dat$SessionID <- as.factor(dat$SessionID)

# --- 3. Clean Data ---
z_thresh <- 2.56
message(sprintf("Trimming extreme outliers beyond +/- %s Z-scores...", z_thresh))
dat_clean <- dat %>%
  filter(
    RipplePower_Z > -z_thresh & RipplePower_Z < z_thresh,
    SpindlePower_Match_Z > -z_thresh & SpindlePower_Match_Z < z_thresh
  ) %>%
  # Standardize the existing Coherence column from the CSV
  mutate(
    Event_Coherence_Post_Z = as.numeric(scale(Event_Coherence_Post))
  )

# ==============================================================================
# --- TASK 5.1: FIT THE MINIMAL ADEQUATE MODEL ---
# ==============================================================================
message("\nFitting the Final Minimal Adequate Model...")

mdl_final <- bam(Event_Coherence_Post ~ 
                   # 1. Surviving Main Power Effects
                   s(RipplePower_Z, k = 5) + 
                   s(SpindlePower_Match_Z, k = 5) + 
                   
                   # 2. The Baseline Phase Landscape (Synergy)
                   te(SOPhase_Match, SOPhase_NonMatch, bs = "cc", k = 8) + 
                   
                   # 3. The Winning Interaction (Local Ripple Gate)
                   ti(RipplePower_Z, SOPhase_Match, bs = c("tp", "cc"), k = c(5, 8)) + 
                   
                   # 4. Control Term
                   s(SessionID, bs = "re"), 
                 
                 data = dat_clean, 
                 method = "fREML", 
                 discrete = TRUE, 
                 nthreads = 4)

message("\n--- FINAL MODEL SUMMARY ---")
print(summary(mdl_final))

# ==============================================================================
# --- TASK 5.2: FORMAL DIAGNOSTICS (GAM.CHECK) ---
# ==============================================================================
message("\n--- RUNNING DIAGNOSTICS (gam.check) ---")
#gam.check(mdl_final)

# ==============================================================================
# --- STAGE 5 VISUALIZATIONS: PUBLICATION FIGURES ---
# ==============================================================================
message("\nGenerating Final Publication Figures...")

# 1. Main Effect: Ripple Power
dev.new(noRStudioGD = TRUE)
p_ripple <- draw(mdl_final, select = "s(RipplePower_Z)", residuals = FALSE, rug = FALSE) + 
  theme_bw() + labs(title = "Main Effect: Ripple Power", x = "Ripple Power (Z)", y = "Partial Effect")
print(p_ripple)

# 2. Main Effect: Match Spindle Power
dev.new(noRStudioGD = TRUE)
p_spindle <- draw(mdl_final, select = "s(SpindlePower_Match_Z)", residuals = FALSE, rug = FALSE) + 
  theme_bw() + labs(title = "Main Effect: Match Spindle Power", x = "Match Spindle Power (Z)", y = "Partial Effect")
print(p_spindle)

# 3. The 2D Joint Phase Landscape
message("Rendering Joint Phase Landscape...")
sm_2d <- smooth_estimates(mdl_final, smooth = "te(SOPhase_Match,SOPhase_NonMatch)", n = 100) %>%
  mutate(
    .lower_ci = .estimate - (1.96 * .se),
    .upper_ci = .estimate + (1.96 * .se)
  )

p_joint <- sm_2d %>%
  mutate(
    SOPhase_Match = ifelse(SOPhase_Match < 0, SOPhase_Match + 2*pi, SOPhase_Match),
    SOPhase_NonMatch = ifelse(SOPhase_NonMatch < 0, SOPhase_NonMatch + 2*pi, SOPhase_NonMatch)
  ) %>%
  arrange(SOPhase_Match, SOPhase_NonMatch) %>%
  ggplot(aes(x = SOPhase_Match, y = SOPhase_NonMatch, fill = .estimate)) +
  geom_tile(width = (2*pi)/99, height = (2*pi)/99) + 
  geom_contour(aes(z = .estimate), color = "black", alpha = 0.2) + 
  scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "firebrick", midpoint = 0, name = "Effect") +
  scale_x_continuous(breaks = c(0, pi, 2*pi), labels = c("0", expression(pi), expression(2*pi))) +
  scale_y_continuous(breaks = c(0, pi, 2*pi), labels = c("0", expression(pi), expression(2*pi))) +
  theme_minimal() +
  labs(title = "Baseline Inter-hemispheric Phase Synergy",
       x = "Match SO Phase (Rad)", y = "Non-Match SO Phase (Rad)")

dev.new(noRStudioGD = TRUE)
print(p_joint)

# 3b. Marginal Effect: Match SO Phase
p_marginal_match <- sm_2d %>%
  mutate(SOPhase_Match = ifelse(SOPhase_Match < 0, SOPhase_Match + 2*pi, SOPhase_Match)) %>%
  group_by(SOPhase_Match) %>%
  summarize(
    marginal_effect = mean(.estimate), 
    marginal_lower = mean(.lower_ci),
    marginal_upper = mean(.upper_ci),
    .groups = 'drop'
  ) %>%
  ggplot(aes(x = SOPhase_Match, y = marginal_effect)) +
  geom_ribbon(aes(ymin = marginal_lower, ymax = marginal_upper), alpha = 0.2, fill = "dodgerblue") +
  geom_line(linewidth = 1, color = "dodgerblue") +
  scale_x_continuous(breaks = c(0, pi, 2*pi), labels = c("0", expression(pi), expression(2*pi))) +
  theme_bw() +
  labs(title = "Marginal Effect: Match SO Phase", x = "Match SO Phase (Rad)", y = "Average Partial Effect")

dev.new(noRStudioGD = TRUE)
print(p_marginal_match)

# 3c. Marginal Effect: Non-Match SO Phase
p_marginal_nonmatch <- sm_2d %>%
  mutate(SOPhase_NonMatch = ifelse(SOPhase_NonMatch < 0, SOPhase_NonMatch + 2*pi, SOPhase_NonMatch)) %>%
  group_by(SOPhase_NonMatch) %>%
  summarize(
    marginal_effect = mean(.estimate), 
    marginal_lower = mean(.lower_ci),
    marginal_upper = mean(.upper_ci),
    .groups = 'drop'
  ) %>%
  ggplot(aes(x = SOPhase_NonMatch, y = marginal_effect)) +
  geom_ribbon(aes(ymin = marginal_lower, ymax = marginal_upper), alpha = 0.2, fill = "firebrick") +
  geom_line(linewidth = 1, color = "firebrick") +
  scale_x_continuous(breaks = c(0, pi, 2*pi), labels = c("0", expression(pi), expression(2*pi))) +
  theme_bw() +
  labs(title = "Marginal Effect: Non-Match SO Phase", x = "Non-Match SO Phase (Rad)", y = "Average Partial Effect")

dev.new(noRStudioGD = TRUE)
print(p_marginal_nonmatch)

# 4. The Local Ripple Gate (Interaction) - SWAPPED AXES
message("Rendering Ripple Gate Interaction (Phase on X, Power on Y)...")

# Force a dense grid to capture full data range
sm_gate <- smooth_estimates(mdl_final, smooth = "ti(RipplePower_Z,SOPhase_Match)", n = 100)

# Re-calculate tile dimensions for swapped orientation
rp_range <- range(sm_gate$RipplePower_Z)
n_unique_rp <- length(unique(sm_gate$RipplePower_Z))
n_unique_ph <- length(unique(sm_gate$SOPhase_Match))
rp_height <- (rp_range[2] - rp_range[1]) / (n_unique_rp - 1)
ph_width <- (2*pi) / (n_unique_ph - 1)

p_gate <- sm_gate %>%
  mutate(SOPhase_Match = ifelse(SOPhase_Match < 0, SOPhase_Match + 2*pi, SOPhase_Match)) %>%
  arrange(SOPhase_Match, RipplePower_Z) %>%
  ggplot(aes(x = SOPhase_Match, y = RipplePower_Z, fill = .estimate)) +
  geom_tile(width = ph_width, height = rp_height) + 
  geom_contour(aes(z = .estimate), color = "black", alpha = 0.2) +
  scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "firebrick", midpoint = 0, name = "Interaction") +
  scale_x_continuous(breaks = c(0, pi, 2*pi), labels = c("0", expression(pi), expression(2*pi))) +
  theme_bw() + 
  labs(title = "The Ripple Gate",
       subtitle = "Coherence boost depends on Local SO Phase",
       x = "Match SO Phase (Rad)", y = "Ripple Power (Z)")

dev.new(noRStudioGD = TRUE)
print(p_gate)

# ==============================================================================
# --- 3. MULTI-METRIC EFFECT SIZE CALCULATIONS ---
# ==============================================================================
message("\nCalculating 4 Effect Size Metrics (This may take a minute)...")

# EXACT formula components to rebuild the models robustly
formula_terms <- c(
  "s(RipplePower_Z, k = 5)", 
  "s(SpindlePower_Match_Z, k = 5)", 
  "te(SOPhase_Match, SOPhase_NonMatch, bs = 'cc', k = 8)", 
  "ti(RipplePower_Z, SOPhase_Match, bs = c('tp', 'cc'), k = c(5, 8))"
)

# EXACT labels output by summary() and smooth_estimates()
smooth_labels <- c(
  "s(RipplePower_Z)", 
  "s(SpindlePower_Match_Z)", 
  "te(SOPhase_Match,SOPhase_NonMatch)", 
  "ti(RipplePower_Z,SOPhase_Match)"
)

# Initialize results dataframe
effect_sizes <- data.frame(Term = smooth_labels, Partial_Deviance = NA, 
                           Eta_Sq_Partial = NA, Amplitude = NA, RMS = NA)

# 3A. Deviance Partitioning (Drop-One)
full_dev <- summary(mdl_final)$dev.expl * 100

message("  -> Running Drop-One Deviance Models:")
for(i in 1:length(formula_terms)) {
  message(sprintf("     Fitting reduced model without: %s", smooth_labels[i]))
  
  # Rebuild formula safely by taking all terms EXCEPT the current one, plus SessionID
  active_terms <- c(formula_terms[-i], "s(SessionID, bs = 're')")
  red_form <- as.formula(paste("Event_Coherence_Post ~", paste(active_terms, collapse = " + ")))
  
  # Fit reduced model
  mdl_red <- bam(red_form, data = dat_clean, method = "fREML", discrete = TRUE, nthreads = 4)
  
  # Calculate Dev Lost
  effect_sizes$Partial_Deviance[i] <- full_dev - (summary(mdl_red)$dev.expl * 100)
}

# 3B. Partial Eta-Squared from F-statistics
sum_tab <- as.data.frame(summary(mdl_final)$s.table)
sum_tab$Term <- rownames(sum_tab)
res_df <- mdl_final$df.residual

for(i in 1:length(smooth_labels)) {
  term_row <- sum_tab[sum_tab$Term == smooth_labels[i], ]
  if(nrow(term_row) == 1) {
    # Eta_p^2 = (F * edf) / (F * edf + df_residual)
    eta_sq <- (term_row$F * term_row$edf) / ((term_row$F * term_row$edf) + res_df)
    effect_sizes$Eta_Sq_Partial[i] <- eta_sq
  }
}

# 3C & 3D. Amplitude (Peak-to-Trough) and RMS (Average Impact)
for(i in 1:length(smooth_labels)) {
  sm <- smooth_estimates(mdl_final, smooth = smooth_labels[i], n = 100)
  effect_sizes$Amplitude[i] <- max(sm$.estimate) - min(sm$.estimate)
  effect_sizes$RMS[i] <- sqrt(mean(sm$.estimate^2))
}

# Print Results cleanly to console
message("\n--- COMPREHENSIVE EFFECT SIZES (RAW UNITS) ---")
print(effect_sizes %>% arrange(desc(Partial_Deviance)))

# ==============================================================================
# --- 4. EFFECT SIZE DASHBOARD VISUALIZATION ---
# ==============================================================================
# Reshape data for faceted plotting
eff_long <- effect_sizes %>%
  pivot_longer(cols = -Term, names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = factor(Metric, 
                         levels = c("Partial_Deviance", "Eta_Sq_Partial", "Amplitude", "RMS"),
                         labels = c("Deviance Explained (%)", "Partial Eta-Squared", 
                                    "Peak-to-Trough Amplitude", "RMS Effect")))

p_dashboard <- ggplot(eff_long, aes(x = reorder(Term, Value), y = Value, fill = Term)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  coord_flip() +
  facet_wrap(~Metric, scales = "free_x", ncol = 2) +
  scale_fill_viridis_d(option = "mako") +
  theme_bw() +
  labs(title = "Effect Size Dashboard",
       subtitle = "Comparing alternative metrics for biological importance",
       x = NULL, y = NULL) +
  theme(strip.text = element_text(face = "bold", size = 11))

dev.new(noRStudioGD = TRUE); print(p_dashboard)

message("\nAnalysis Pipeline Complete!")