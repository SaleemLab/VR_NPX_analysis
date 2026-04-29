# ==============================================================================
# V1-HC COHERENCE: FORMAL LME vs. GAM COMPARISON
# ==============================================================================
# Objective: Prove that the relationships driving V1-HC coherence are strictly
# non-linear by pitting a Linear Mixed Model against the Generalized Additive Model.

packages <- c("mgcv", "gratia", "ggplot2", "dplyr")
for (p in packages) if (!require(p, character.only = TRUE)) install.packages(p)

library(mgcv); library(gratia); library(ggplot2); library(dplyr)

# --- 1. Data Preparation ---
message("Loading and preparing data...")
dat <- read.csv("D:/Code/MasaCoherenceDataAnalysis/v1_hc_data.csv")
dat$SessionID <- as.factor(dat$SessionID)

z_thresh <- 2.56
dat_clean <- dat %>%
  filter(RipplePower_Z > -z_thresh & RipplePower_Z < z_thresh,
         SpindlePower_Match_Z > -z_thresh & SpindlePower_Match_Z < z_thresh)

# ==============================================================================
# --- 2. FIT THE MODELS ---
# ==============================================================================

# --- Model A: The Linear Mixed-Effects Model (LME) ---
# We use bam() with linear terms to ensure the underlying likelihood estimator
# is identical to the GAM, allowing for a perfectly valid AIC/ANOVA comparison.
message("\nFitting Model A: Linear Mixed-Effects Model (LME)...")
mdl_lme <- bam(Event_Coherence_Post ~ 
                 # Linear Main Effects
                 RipplePower_Z + 
                 SpindlePower_Match_Z + 
                 SOPhase_Match + 
                 SOPhase_NonMatch + 
                 # Linear Interaction (Analogous to the ti() Gate)
                 RipplePower_Z:SOPhase_Match + 
                 # Random Effect (identical across models)
                 s(SessionID, bs = "re"), 
               data = dat_clean, method = "fREML", discrete = TRUE, nthreads = 4)

# --- Model B: The Generalized Additive Model (GAM) ---
message("Fitting Model B: Generalized Additive Model (GAM)...")
mdl_gam <- bam(Event_Coherence_Post ~ 
                 # Non-linear Main Effects
                 s(RipplePower_Z, k = 5) + 
                 s(SpindlePower_Match_Z, k = 5) + 
                 # Non-linear Synergy & Gating
                 te(SOPhase_Match, SOPhase_NonMatch, bs = "cc", k = 8) + 
                 ti(RipplePower_Z, SOPhase_Match, bs = c("tp", "cc"), k = c(5, 8)) + 
                 # Random Effect
                 s(SessionID, bs = "re"), 
               data = dat_clean, method = "fREML", discrete = TRUE, nthreads = 4)


# ==============================================================================
# --- 3. FORMAL COMPARISON ---
# ==============================================================================

message("\n=========================================================")
message("--- RESULTS: LME vs. GAM ---")
message("=========================================================")

# 3A. ANOVA (Chi-Square Test)
message("\n1. ANOVA (Chi-Square Test for Non-Linearity)")
message("Tests if the added complexity of the GAM explains significantly more variance.")
print(anova(mdl_lme, mdl_gam, test = "Chisq"))

# 3B. AIC Comparison
message("\n2. Akaike Information Criterion (AIC)")
message("Lower is better. A Delta AIC > 10 provides overwhelming evidence for the winning model.")
aic_results <- AIC(mdl_lme, mdl_gam)
aic_results$Model <- c("Linear Mixed Model (LME)", "Generalized Additive Model (GAM)")
aic_results$Delta_AIC <- aic_results$AIC - min(aic_results$AIC)

# Reorder columns for clean printing
aic_results <- aic_results %>% select(Model, df, AIC, Delta_AIC) %>% arrange(AIC)
print(aic_results)

# 3C. Deviance Explained
message("\n3. Total Variance (Deviance) Explained")
dev_lme <- summary(mdl_lme)$dev.expl * 100
dev_gam <- summary(mdl_gam)$dev.expl * 100
message(sprintf("LME: %.2f%%", dev_lme))
message(sprintf("GAM: %.2f%%", dev_gam))
message(sprintf("The GAM captures %.2f%% MORE of the true biological variance.", (dev_gam - dev_lme)))

# ==============================================================================
# --- 4. VISUALIZATION: WHY THE LINEAR MODEL FAILS ---
# ==============================================================================
# A quick plot to show how the linear model fundamentally misinterprets the circular phase.

message("\nGenerating comparison plot...")

# Predict across the phase range holding other variables at 0 (mean)
pred_data <- data.frame(
  RipplePower_Z = 0,
  SpindlePower_Match_Z = 0,
  SOPhase_NonMatch = pi,
  SOPhase_Match = seq(0, 2*pi, length.out = 100),
  SessionID = dat_clean$SessionID[1] # Dummy session for prediction
)

# Get predictions (excluding the random effect for the plot)
pred_data$LME_Prediction <- predict(mdl_lme, newdata = pred_data, exclude = "s(SessionID)", type = "response")
pred_data$GAM_Prediction <- predict(mdl_gam, newdata = pred_data, exclude = "s(SessionID)", type = "response")

# Plot
p_compare <- ggplot(pred_data, aes(x = SOPhase_Match)) +
  geom_line(aes(y = LME_Prediction, color = "Linear Model"), linewidth = 1.2, linetype = "dashed") +
  geom_line(aes(y = GAM_Prediction, color = "GAM"), linewidth = 1.2) +
  scale_color_manual(values = c("Linear Model" = "darkgray", "GAM" = "dodgerblue3")) +
  scale_x_continuous(breaks = c(0, pi, 2*pi), labels = c("0", expression(pi), expression(2*pi))) +
  theme_bw() +
  labs(title = "Why we use GAMs for Neural Oscillations",
       subtitle = "Linear models force a straight line through circular phase data",
       x = "SO Phase (Rad)", y = "Predicted Coherence",
       color = "Model Fit") +
  theme(legend.position = "top")

dev.new(noRStudioGD = TRUE)
print(p_compare)

message("\nComparison script complete! Check your console for the final statistics.")