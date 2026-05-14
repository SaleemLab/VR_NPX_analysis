# ==============================================================================
# V1-HC COHERENCE: STAGE 3.2 - RIPPLE x SPINDLE SYNERGY
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
dat <- read.csv("C:/Users/masah/Documents/GitHub/VR_NPX_analysis/predicting_coherence_with_GAMs/v1_hc_data_geo3.csv")
dat$SessionID <- as.factor(dat$SessionID)

# --- 3. Clean Data ---
z_thresh <- 2.56  

message(sprintf("Trimming extreme outliers beyond +/- %s Z-scores...", z_thresh))
dat_clean <- dat %>%
  filter(
    RipplePower_Z > -z_thresh & RipplePower_Z < z_thresh,
    SpindlePower_Match_Z > -z_thresh & SpindlePower_Match_Z < z_thresh
    # Note: SpindlePower_NonMatch_Z is officially removed from all data filtering and modeling
  )

# ==============================================================================
# --- TASK 3.2: FIT THE RIPPLE X SPINDLE SYNERGY MODELS ---
# ==============================================================================
message("\nFitting Baseline Model (Main Effects + Joint Phase)...")
# Our fully verified baseline from Stages 1 and 2
mdl_base <- bam(Event_Coherence_Pre_GeoMean ~ 
                  s(RipplePower_Z, k = 5) + 
                  s(SpindlePower_Match_Z, k = 5) + 
                  te(SOPhase_Match, SOPhase_NonMatch, bs = "cc", k = 8) + 
                  s(SessionID, bs = "re"), 
                data = dat_clean, method = "fREML", discrete = TRUE, nthreads = 4)

message("Fitting Synergy Model (+ Ripple/Spindle Interaction)...")
mdl_synergy <- bam(Event_Coherence_Pre_GeoMean ~ 
                     s(RipplePower_Z, k = 5) + 
                     s(SpindlePower_Match_Z, k = 5) + 
                     te(SOPhase_Match, SOPhase_NonMatch, bs = "cc", k = 8) + 
                     # The new Cross-Oscillator Interaction Term:
                     ti(RipplePower_Z, SpindlePower_Match_Z, k = c(5, 5)) + 
                     s(SessionID, bs = "re"), 
                   data = dat_clean, method = "fREML", discrete = TRUE, nthreads = 4)

print(summary(mdl_base))
print(summary(mdl_synergy))
# ==============================================================================
# --- TASK 3.3: MODEL COMPARISON ---
# ==============================================================================
message("\n--- MODEL COMPARISON: SUMMARY ---")
# Let's see if the interaction term itself is significant
print(summary(mdl_synergy))

message("\n--- MODEL COMPARISON: ANOVA ---")
print(anova(mdl_base, mdl_synergy, test = "Chisq"))

message("\n--- MODEL COMPARISON: AIC ---")
aic_results <- AIC(mdl_base, mdl_synergy)
aic_results$Delta_AIC <- aic_results$AIC - min(aic_results$AIC)
print(aic_results)

# ==============================================================================
# --- STAGE 3.2 VISUALIZATIONS ---
# ==============================================================================
message("\nGenerating Interaction Heatmap...")

# Visualizing the potential "Super-Coherence" state
dev.new(noRStudioGD = TRUE)
p_interaction <- draw(mdl_synergy, select = "ti(RipplePower_Z,SpindlePower_Match_Z)", rug = FALSE) + 
  scale_fill_viridis_c(option = "plasma", name = "Interaction\nEffect") +
  theme_bw() + 
  labs(title = "Stage 3.2: Ripple x Spindle Synergy",
       subtitle = "Does high ripple power amplify the coherence boost from a spindle?",
       x = "Ripple Power (Z)", y = "Match Spindle Power (Z)")

print(p_interaction)

message("Stage 3.2 complete. Check the p-value of the ti() term and the Delta AIC.")