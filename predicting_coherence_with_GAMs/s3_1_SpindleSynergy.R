# ==============================================================================
# V1-HC COHERENCE: STAGE 3.1 - SPINDLE COORDINATION
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
dat <- read.csv("D:/Code/MasaCoherenceDataAnalysis/v1_hc_data_prob.csv")
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
# --- TASK 3.1: FIT THE SPINDLE SYNERGY MODELS ---
# ==============================================================================
message("\nFitting Baseline Model (Main Effects + Joint Phase)...")
# Note: We include the NonMatch Spindle main effect here to obey the principle of marginality.
mdl_base <- bam(Event_Coherence_Post ~ 
                  s(RipplePower_Z, k = 5) + 
                  s(SpindlePower_Match_Z, k = 5) + 
                  s(SpindlePower_NonMatch_Z, k = 5) + 
                  te(SOPhase_Match, SOPhase_NonMatch, bs = "cc", k = 8) + 
                  s(SessionID, bs = "re"), 
                data = dat_clean, method = "fREML", discrete = TRUE, nthreads = 4)

message("Fitting Synergy Model (+ Match/NonMatch Spindle Interaction)...")
mdl_synergy <- bam(Event_Coherence_Post ~ 
                     s(RipplePower_Z, k = 5) + 
                     s(SpindlePower_Match_Z, k = 5) + 
                     s(SpindlePower_NonMatch_Z, k = 5) + 
                     te(SOPhase_Match, SOPhase_NonMatch, bs = "cc", k = 8) + 
                     # The new Interaction Term:
                     ti(SpindlePower_Match_Z, SpindlePower_NonMatch_Z, k = c(5, 5)) + 
                     s(SessionID, bs = "re"), 
                   data = dat_clean, method = "fREML", discrete = TRUE, nthreads = 4)

# ==============================================================================
# --- TASK 3.2: MODEL COMPARISON ---
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
# --- STAGE 3.1 VISUALIZATIONS ---
# ==============================================================================
message("\nGenerating Interaction Heatmap...")

# We only plot the interaction if you want to see what it looks like.
dev.new(noRStudioGD = TRUE)
p_interaction <- draw(mdl_synergy, select = "ti(SpindlePower_Match_Z,SpindlePower_NonMatch_Z)", rug = FALSE) + 
  scale_fill_viridis_c(option = "plasma", name = "Interaction\nEffect") +
  theme_bw() + 
  labs(title = "Stage 3.1: Spindle Coordination",
       subtitle = "Does coherence spike only when BOTH spindles are active?",
       x = "Match Spindle Power (Z)", y = "Non-Match Spindle Power (Z)")

print(p_interaction)

message("Stage 3.1 complete. Check the p-value of the ti() term and the Delta AIC.")