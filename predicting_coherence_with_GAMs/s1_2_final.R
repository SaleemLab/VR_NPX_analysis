# ==============================================================================
# V1-HC COHERENCE: STAGE 1 (FINAL) - REFINED BASELINE MODEL
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

# Ensure categorical variables are factors
dat$SessionID <- as.factor(dat$SessionID)
dat$AnimalID <- as.factor(dat$AnimalID)

# --- 3. Clean Data (Configurable Threshold) ---
z_thresh <- 2.56  

message(sprintf("Trimming extreme outliers beyond +/- %s Z-scores...", z_thresh))
dat_clean <- dat %>%
  filter(
    RipplePower_Z > -z_thresh & RipplePower_Z < z_thresh,
    SpindlePower_Match_Z > -z_thresh & SpindlePower_Match_Z < z_thresh,
    SpindlePower_NonMatch_Z > -z_thresh & SpindlePower_NonMatch_Z < z_thresh
  )

message(sprintf("Removed %d extreme data points.", nrow(dat) - nrow(dat_clean)))

# ==============================================================================
# --- TASK 1.4: FIT FINAL STAGE 1 MODEL ---
# ==============================================================================
message("\nFitting Final Stage 1 Model...")
message("-> Dropped: SpindlePower_NonMatch_Z (Redundant)")
message("-> Dropped: AnimalID (Non-informative)")

mdl <- bam(Event_Coherence_Post ~ 
                          # Surviving Power main effects
                          s(RipplePower_Z, k = 5) + 
                          s(SpindlePower_Match_Z, k = 5) + 
                          
                          # Independent 1D phase terms (to be tested against 2D in Stage 2)
                          s(SOPhase_Match, bs = "cc", k = 8) + 

                          # Surviving Random effect
                          s(SessionID, bs = "re"), 
                        data = dat_clean, 
                        method = "fREML", 
                        discrete = TRUE,
                        nthreads = 4)

message("\n--- Final Stage 1 Model Summary ---")
print(summary(mdl))

# ==============================================================================
# --- STAGE 1 (FINAL) VISUALIZATIONS ---
# ==============================================================================
message("Generating Partial Effect Plots for surviving 1D terms...")

# 1. Surviving Power Terms
dev.new(noRStudioGD = TRUE)
p1 <- draw(mdl, select = "s(RipplePower_Z)", residuals = FALSE, rug = FALSE) + 
  theme_bw() + labs(title = "Final Stage 1: Ripple Power")
print(p1)

dev.new(noRStudioGD = TRUE)
p2 <- draw(mdl, select = "s(SpindlePower_Match_Z)", residuals = FALSE, rug = FALSE) + 
  theme_bw() + labs(title = "Final Stage 1: Match Spindle Power")
print(p2)

# 2. Independent Phase Terms
dev.new(noRStudioGD = TRUE)
p3 <- draw(mdl, select = "s(SOPhase_Match)", residuals = FALSE, rug = FALSE) + 
  scale_x_continuous(breaks = c(0, pi, 2*pi), labels = c("0", expression(pi), expression(2*pi))) +
  theme_bw() + labs(title = "Final Stage 1: Match SO Phase (1D)")
print(p3)


message("Stage 1 Final execution complete. This model serves as the foundation for Stage 2.")