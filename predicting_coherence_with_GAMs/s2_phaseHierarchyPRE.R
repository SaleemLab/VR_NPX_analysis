# ==============================================================================
# V1-HC COHERENCE: STAGE 2 - SO PHASE HIERARCHY
# ==============================================================================

# --- 1. Load Required Libraries ---
packages <- c("mgcv", "gratia", "ggplot2", "dplyr")
for (p in packages) {
  if (!require(p, character.only = TRUE)) install.packages(p)
}

library(mgcv)
library(gratia)
library(ggplot2)
library(dplyr)

# --- 2. Load and Prepare Data ---
message("Loading data...")
dat <- read.csv("D:/Code/MasaCoherenceDataAnalysis/v1_hc_data_geo2.csv")

dat$SessionID <- as.factor(dat$SessionID)

# --- 3. Clean Data ---
z_thresh <- 2.56  

message(sprintf("Trimming extreme outliers beyond +/- %s Z-scores...", z_thresh))
dat_clean <- dat %>%
  filter(
    RipplePower_Z > -z_thresh & RipplePower_Z < z_thresh,
    SpindlePower_Match_Z > -z_thresh & SpindlePower_Match_Z < z_thresh,
    SpindlePower_NonMatch_Z > -z_thresh & SpindlePower_NonMatch_Z < z_thresh
  )

# ==============================================================================
# --- TASK 2.1: FIT THE HIERARCHY MODELS ---
# ==============================================================================
# Note: We have permanently dropped s(SpindlePower_NonMatch_Z) and s(AnimalID) 
# based on the Stage 1 Redundancy Audit.

message("\nFitting Model A (Unilateral: Match Phase Only)...")
mdl_A_unilateral <- bam(Event_Coherence_Pre_GeoMean ~ 
                          s(RipplePower_Z, k = 5) + 
                          s(SpindlePower_Match_Z, k = 5) + 
                          s(SOPhase_Match, bs = "cc", k = 8) + 
                          s(SessionID, bs = "re"), 
                        data = dat_clean, method = "fREML", discrete = TRUE, nthreads = 4)

message("Fitting Model B (Additive: Match + Non-Match Phase)...")
mdl_B_additive <- bam(Event_Coherence_Pre_GeoMean ~ 
                        s(RipplePower_Z, k = 5) + 
                        s(SpindlePower_Match_Z, k = 5) + 
                        s(SOPhase_Match, bs = "cc", k = 8) + 
                        s(SOPhase_NonMatch, bs = "cc", k = 8) + 
                        s(SessionID, bs = "re"), 
                      data = dat_clean, method = "fREML", discrete = TRUE, nthreads = 4)

message("Fitting Model C (Joint: 2D Phase Synergy)...")
mdl_C_joint <- bam(Event_Coherence_Pre_GeoMean ~ 
                     s(RipplePower_Z, k = 5) + 
                     s(SpindlePower_Match_Z, k = 5) + 
                     te(SOPhase_Match, SOPhase_NonMatch, bs = "cc", k = 8) + 
                     s(SessionID, bs = "re"), 
                   data = dat_clean, method = "fREML", discrete = TRUE, nthreads = 4)

print(summary(mdl_A_unilateral))
print(summary(mdl_B_additive))
print(summary(mdl_C_joint))
# ==============================================================================
# --- TASK 2.2: MODEL COMPARISON ---
# ==============================================================================
message("\n--- MODEL COMPARISON: ANOVA ---")
# Tests if each step-up in complexity explains significantly more variance
print(anova(mdl_A_unilateral, mdl_B_additive, mdl_C_joint, test = "Chisq"))

message("\n--- MODEL COMPARISON: AIC ---")
# Lower AIC is better. A difference (Delta AIC) > 10 is considered overwhelming evidence.
aic_results <- AIC(mdl_A_unilateral, mdl_B_additive, mdl_C_joint)
aic_results$Delta_AIC <- aic_results$AIC - min(aic_results$AIC)
print(aic_results)

# ==============================================================================
# --- STAGE 2 VISUALIZATIONS: THE JOINT LANDSCAPE ---
# ==============================================================================
# We will proactively plot Model C (The Joint Landscape) as it's the strongest hypothesis.
# If Model B or A wins the AIC test, we can revise this!

message("\nGenerating 2D Phase Landscape for Model C (Shifted to center \u03c0)...")

sm_2d <- smooth_estimates(mdl_C_joint, smooth = "te(SOPhase_Match,SOPhase_NonMatch)")

p_joint <- sm_2d %>%
  mutate(
    SOPhase_Match = ifelse(SOPhase_Match < 0, SOPhase_Match + 2*pi, SOPhase_Match),
    SOPhase_NonMatch = ifelse(SOPhase_NonMatch < 0, SOPhase_NonMatch + 2*pi, SOPhase_NonMatch)
  ) %>%
  arrange(SOPhase_Match, SOPhase_NonMatch) %>%
  ggplot(aes(x = SOPhase_Match, y = SOPhase_NonMatch, fill = .estimate)) +
  geom_tile(width = (2*pi)/99, height = (2*pi)/99) + 
  scale_fill_viridis_c(option = "D", name = "Effect") +
  scale_x_continuous(breaks = c(0, pi, 2*pi), labels = c("0", expression(pi), expression(2*pi))) +
  scale_y_continuous(breaks = c(0, pi, 2*pi), labels = c("0", expression(pi), expression(2*pi))) +
  theme_minimal() +
  labs(title = "Stage 2: Joint Phase Landscape (Model C)",
       subtitle = "If Model C wins AIC, this is your baseline phase rule",
       x = "Match Phase (Rad)", y = "Non-Match Phase (Rad)")

dev.new(noRStudioGD = TRUE)
print(p_joint)

message("Stage 2 execution complete. Check ANOVA and Delta_AIC outputs to declare a winner!")