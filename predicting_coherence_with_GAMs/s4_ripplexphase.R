# ==============================================================================
# V1-HC COHERENCE: STAGE 4 - CROSS-MODALITY GATING (POWER x PHASE)
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
    SpindlePower_Match_Z > -z_thresh & SpindlePower_Match_Z < z_thresh
  )

# ==============================================================================
# --- TASK 4.1: FIT THE CROSS-MODALITY MODELS ---
# ==============================================================================
message("\nFitting Baseline Model (Main Effects + Joint Phase)...")
mdl_base <- bam(Event_Coherence_Post ~ 
                  s(RipplePower_Z, k = 5) + 
                  s(SpindlePower_Match_Z, k = 5) + 
                  te(SOPhase_Match, SOPhase_NonMatch, bs = "cc", k = 8) + 
                  s(SessionID, bs = "re"), 
                data = dat_clean, method = "fREML", discrete = TRUE, nthreads = 4)

message("Fitting Local Gating Model (+ Ripple x Match Phase)...")
mdl_local <- bam(Event_Coherence_Post ~ 
                   s(RipplePower_Z, k = 5) + 
                   s(SpindlePower_Match_Z, k = 5) + 
                   te(SOPhase_Match, SOPhase_NonMatch, bs = "cc", k = 8) + 
                   # Local interaction: bs=c("tp", "cc") handles continuous vs circular
                   ti(RipplePower_Z, SOPhase_Match, bs = c("tp", "cc"), k = c(5, 8)) + 
                   s(SessionID, bs = "re"), 
                 data = dat_clean, method = "fREML", discrete = TRUE, nthreads = 4)

message("Fitting Global Gating Model (+ Ripple x Joint Phase Landscape)...")
mdl_global <- bam(Event_Coherence_Post ~ 
                    s(RipplePower_Z, k = 5) + 
                    s(SpindlePower_Match_Z, k = 5) + 
                    te(SOPhase_Match, SOPhase_NonMatch, bs = "cc", k = 8) + 
                    # Global 3D interaction: k=c(5,5,5) used to prevent memory crashes
                    ti(RipplePower_Z, SOPhase_Match, SOPhase_NonMatch, bs = c("tp", "cc", "cc"), k = c(5, 5, 5)) + 
                    s(SessionID, bs = "re"), 
                  data = dat_clean, method = "fREML", discrete = TRUE, nthreads = 4)

# ==============================================================================
# --- TASK 4.2: MODEL COMPARISON ---
# ==============================================================================
message("\n--- MODEL COMPARISON: SUMMARY OF TI() TERMS ---")
message("Local Interaction p-value:")
print(summary(mdl_local)$s.table["ti(RipplePower_Z,SOPhase_Match)", "p-value"])

message("Global Interaction p-value:")
print(summary(mdl_global)$s.table["ti(RipplePower_Z,SOPhase_Match,SOPhase_NonMatch)", "p-value"])

message("\n--- MODEL COMPARISON: ANOVA ---")
print(anova(mdl_base, mdl_local, mdl_global, test = "Chisq"))

message("\n--- MODEL COMPARISON: AIC ---")
aic_results <- AIC(mdl_base, mdl_local, mdl_global)
aic_results$Delta_AIC <- aic_results$AIC - min(aic_results$AIC)
print(aic_results)

# ==============================================================================
# --- STAGE 4 VISUALIZATIONS ---
# ==============================================================================
message("\nGenerating Interaction Plots...")

# 1. Local Gating Plot (2D Surface)
dev.new(noRStudioGD = TRUE)
p_local <- draw(mdl_local, select = "ti(RipplePower_Z,SOPhase_Match)", rug = FALSE) + 
  scale_fill_viridis_c(option = "plasma", name = "Interaction\nEffect") +
  scale_y_continuous(breaks = c(0, pi, 2*pi), labels = c("0", expression(pi), expression(2*pi))) +
  theme_bw() + 
  labs(title = "Stage 4.1: Local Ripple Gating",
       subtitle = "Does ripple power matter more at specific local SO phases?",
       x = "Ripple Power (Z)", y = "Match SO Phase (Rad)")
print(p_local)

# 2. Global Gating Plot (Slices of 3D Surface)
# gratia handles 3D interactions by taking slices across the continuous predictor
dev.new(noRStudioGD = TRUE)
p_global <- draw(mdl_global, select = "ti(RipplePower_Z,SOPhase_Match,SOPhase_NonMatch)", rug = FALSE) +
  theme_bw() +
  labs(title = "Stage 4.2: Global Phase Gating",
       subtitle = "How Ripple Power interacts with the full 2D Phase Landscape")
print(p_global)

message("Stage 4 complete. Check the AIC hierarchy to determine if the Local or Global gate wins.")