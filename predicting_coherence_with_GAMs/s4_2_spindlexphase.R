# ==============================================================================
# V1-HC COHERENCE: STAGE 4 - SPINDLE x PHASE GATING
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
dat <- read.csv("D:/Code/MasaCoherenceDataAnalysis/v1_hc_data.csv")
dat$SessionID <- as.factor(dat$SessionID)

# --- 3. Clean Data ---
z_thresh <- 2.56  

message(sprintf("Trimming extreme outliers beyond +/- %s Z-scores...", z_thresh))
dat_clean <- dat %>%
  filter(
    RipplePower_Z > -z_thresh & RipplePower_Z < z_thresh,
    SpindlePower_Match_Z > -z_thresh & SpindlePower_Match_Z < z_thresh
  )

# ==============================================================================
# --- TASK 4.1: DEFINE THE TWO BASELINES ---
# ==============================================================================
message("\nFitting 'Old' Baseline (Main Effects + Joint Phase Only)...")
mdl_base_old <- bam(Event_Coherence_Post ~ 
                      s(RipplePower_Z, k = 5) + 
                      s(SpindlePower_Match_Z, k = 5) + 
                      te(SOPhase_Match, SOPhase_NonMatch, bs = "cc", k = 8) + 
                      s(SessionID, bs = "re"), 
                    data = dat_clean, method = "fREML", discrete = TRUE, nthreads = 4)

message("Fitting 'New' Baseline (+ The Winning Local Ripple Gate)...")
mdl_base_new <- bam(Event_Coherence_Post ~ 
                      s(RipplePower_Z, k = 5) + 
                      s(SpindlePower_Match_Z, k = 5) + 
                      te(SOPhase_Match, SOPhase_NonMatch, bs = "cc", k = 8) + 
                      ti(RipplePower_Z, SOPhase_Match, bs = c("tp", "cc"), k = c(5, 8)) + 
                      s(SessionID, bs = "re"), 
                    data = dat_clean, method = "fREML", discrete = TRUE, nthreads = 4)

# ==============================================================================
# --- TASK 4.2: TEST SPINDLE GATING AGAINST OLD BASELINE ---
# ==============================================================================
message("\nTesting Spindle Gating (Against OLD Baseline)...")

mdl_spin_local_old <- bam(Event_Coherence_Post ~ 
                            s(RipplePower_Z, k = 5) + 
                            s(SpindlePower_Match_Z, k = 5) + 
                            te(SOPhase_Match, SOPhase_NonMatch, bs = "cc", k = 8) + 
                            ti(SpindlePower_Match_Z, SOPhase_Match, bs = c("tp", "cc"), k = c(5, 8)) + 
                            s(SessionID, bs = "re"), 
                          data = dat_clean, method = "fREML", discrete = TRUE, nthreads = 4)

mdl_spin_global_old <- bam(Event_Coherence_Post ~ 
                             s(RipplePower_Z, k = 5) + 
                             s(SpindlePower_Match_Z, k = 5) + 
                             te(SOPhase_Match, SOPhase_NonMatch, bs = "cc", k = 8) + 
                             ti(SpindlePower_Match_Z, SOPhase_Match, SOPhase_NonMatch, bs = c("tp", "cc", "cc"), k = c(5, 5, 5)) + 
                             s(SessionID, bs = "re"), 
                           data = dat_clean, method = "fREML", discrete = TRUE, nthreads = 4)

# ==============================================================================
# --- TASK 4.3: TEST SPINDLE GATING AGAINST NEW BASELINE ---
# ==============================================================================
message("\nTesting Spindle Gating (Against NEW Baseline)...")

mdl_spin_local_new <- bam(Event_Coherence_Post ~ 
                            s(RipplePower_Z, k = 5) + 
                            s(SpindlePower_Match_Z, k = 5) + 
                            te(SOPhase_Match, SOPhase_NonMatch, bs = "cc", k = 8) + 
                            ti(RipplePower_Z, SOPhase_Match, bs = c("tp", "cc"), k = c(5, 8)) + # The Ripple Gate
                            ti(SpindlePower_Match_Z, SOPhase_Match, bs = c("tp", "cc"), k = c(5, 8)) + # The Spindle Gate
                            s(SessionID, bs = "re"), 
                          data = dat_clean, method = "fREML", discrete = TRUE, nthreads = 4)

mdl_spin_global_new <- bam(Event_Coherence_Post ~ 
                             s(RipplePower_Z, k = 5) + 
                             s(SpindlePower_Match_Z, k = 5) + 
                             te(SOPhase_Match, SOPhase_NonMatch, bs = "cc", k = 8) + 
                             ti(RipplePower_Z, SOPhase_Match, bs = c("tp", "cc"), k = c(5, 8)) + # The Ripple Gate
                             ti(SpindlePower_Match_Z, SOPhase_Match, SOPhase_NonMatch, bs = c("tp", "cc", "cc"), k = c(5, 5, 5)) + 
                             s(SessionID, bs = "re"), 
                           data = dat_clean, method = "fREML", discrete = TRUE, nthreads = 4)

# ==============================================================================
# --- TASK 4.4: MODEL COMPARISON ---
# ==============================================================================
message("\n--- COMPARING AGAINST OLD BASELINE (Does Spindle Gating Exist?) ---")
print(anova(mdl_base_old, mdl_spin_local_old, mdl_spin_global_old, test = "Chisq"))

message("\n--- COMPARING AGAINST NEW BASELINE (Does Spindle Gating Survive?) ---")
print(anova(mdl_base_new, mdl_spin_local_new, mdl_spin_global_new, test = "Chisq"))

message("\n--- MASTER AIC TOURNAMENT ---")
# Let's pit all models against each other to see the ultimate winner
aic_results <- AIC(mdl_base_old, mdl_spin_local_old, mdl_spin_global_old, 
                   mdl_base_new, mdl_spin_local_new, mdl_spin_global_new)
aic_results$Delta_AIC <- aic_results$AIC - min(aic_results$AIC)
print(aic_results[order(aic_results$Delta_AIC), ]) # Sorts to put the winner at the top

# ==============================================================================
# --- STAGE 4 VISUALIZATIONS ---
# ==============================================================================
message("\nGenerating Interaction Plots for Spindle x Local Phase...")

dev.new(noRStudioGD = TRUE)
p_spindle_local <- draw(mdl_spin_local_old, select = "ti(SpindlePower_Match_Z,SOPhase_Match)", rug = FALSE) + 
  scale_fill_viridis_c(option = "plasma", name = "Interaction\nEffect") +
  scale_y_continuous(breaks = c(0, pi, 2*pi), labels = c("0", expression(pi), expression(2*pi))) +
  theme_bw() + 
  labs(title = "Stage 4: Local Spindle Gating",
       subtitle = "Does spindle power matter more at specific local SO phases?",
       x = "Match Spindle Power (Z)", y = "Match SO Phase (Rad)")

print(p_spindle_local)

message("Stage 4 Spindle test complete. Look at the top row of the MASTER AIC TOURNAMENT to see the overall winner.")