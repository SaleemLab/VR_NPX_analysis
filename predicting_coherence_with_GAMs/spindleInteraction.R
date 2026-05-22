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
#dat <- read.csv("C:/Users/masah/Documents/GitHub/VR_NPX_analysis/predicting_coherence_with_GAMs/v1_hc_data_geo3.csv")
#dat <- read.csv("C:/Users/masah/Documents/GitHub/VR_NPX_analysis/predicting_coherence_with_GAMs/v1_hc_data_geo_100ms.csv")
#dat <- read.csv("C:/Users/masah/Documents/GitHub/VR_NPX_analysis/predicting_coherence_with_GAMs/v1_hc_data_geo_PRE.csv")
#dat <- read.csv("C:/Users/masah/Documents/GitHub/VR_NPX_analysis/predicting_coherence_with_GAMs/v1_hc_data_geo_100_200ms.csv")
#dat <- read.csv("C:/Users/masah/Documents/GitHub/VR_NPX_analysis/predicting_coherence_with_GAMs/v1_hc_data_geo_0_100ms_POST.csv")
#dat <- read.csv("C:/Users/masah/Documents/GitHub/VR_NPX_analysis/predicting_coherence_with_GAMs/v1_hc_data_geo_100_200ms_POST.csv")
#dat <- read.csv("C:/Users/masah/Documents/GitHub/VR_NPX_analysis/predicting_coherence_with_GAMs/v1_hc_data_geo_0_200ms_POST.csv")

dat <- read.csv("C:/Users/masah/Documents/GitHub/VR_NPX_analysis/predicting_coherence_with_GAMs/v1_hc_data_geo_-200_200ms.csv")
dat <- read.csv("C:/Users/masah/Documents/GitHub/VR_NPX_analysis/predicting_coherence_with_GAMs/v1_hc_data_geo_-100_100ms.csv")

#dat <- read.csv("C:/Users/masah/Documents/GitHub/VR_NPX_analysis/predicting_coherence_with_GAMs/v1_hc_data_geo_0_100ms_POST.csv")
#dat <- read.csv("C:/Users/masah/Documents/GitHub/VR_NPX_analysis/predicting_coherence_with_GAMs/v1_hc_data_geo_0_100ms_PRE2.csv")

dat$SessionID <- as.factor(dat$SessionID)

# --- 3. Clean Data ---
z_thresh <- 2.56  # include <99th centiles of data

message(sprintf("Trimming extreme outliers beyond +/- %s Z-scores...", z_thresh))
dat_clean <- dat %>%
  filter(
    RipplePower_Z > -z_thresh & RipplePower_Z < z_thresh,
    SpindlePower_Match_Z > -z_thresh & SpindlePower_Match_Z < z_thresh,
    SpindlePower_NonMatch_Z > -z_thresh & SpindlePower_NonMatch_Z < z_thresh,
  )
  # Standardize the existing Coherence column from the CSV

# Standardize the Dependent Variable (Coherence)
# This makes 'Partial Effect' units equal to Standard Deviations of Coherence
  #mutate(Event_Coherence_Post_Z = as.numeric(scale(Event_Coherence_Post)))

# ==============================================================================
# --- TASK 5.1: FIT THE MINIMAL ADEQUATE MODEL ---
# ==============================================================================
message("\nFitting the Final Minimal Adequate Model...")

mdl_base <- bam(Event_Coherence_Post_GeoMean ~ 
                       # 1. Surviving Main Power Effects
                       s(RipplePower_Z, k = 5) + 
                       s(SpindlePower_Match_Z, k = 5) + 
                       s(SpindlePower_NonMatch_Z, k = 5) + 
                       ti(SpindlePower_Match_Z, SpindlePower_NonMatch_Z, k = 5) +
                       #te(SpindlePower_Match_Z, SpindlePower_NonMatch_Z, k = c(5, 5))+
                       
                       # 2. The Baseline Phase Landscape (Synergy)
                       te(SOPhase_Match, SOPhase_NonMatch, bs = "cc", k = 8) + 
                       ti(SOPhase_Match, SOPhase_NonMatch, bs = "cc", k = 8) + 
                       s(SOPhase_Match, bs = "cc", k = 8) + 
                       #s(SOPhase_NonMatch, bs = "cc", k = 8) + 
                       
                       # 3. The Winning Interaction (Local Ripple Gate)
                       #ti(RipplePower_Z, SOPhase_Match, bs = c("tp", "cc"), k = c(5, 8)) + 
                       #ti(SpindlePower_Match_Z, SOPhase_Match, bs = c("tp", "cc"), k = c(5, 8)) + 
                       #ti(SpindlePower_Match_Z, RipplePower_Z, bs = c("tp", "tp"), k = c(5, 5)) + 
                       
                       
                       # 4. Control Term
                       s(AnimalID, bs = "re") +
                       s(SessionID, bs = "re"), 
                     
                     data = dat_clean, 
                     method = "fREML", 
                     discrete = TRUE, 
                     nthreads = 4)

message("\n--- FINAL MODEL SUMMARY ---")
print(summary(mdl_base))



# ==============================================================================
# --- Spindle SO interaction
# ==============================================================================

mdl_spindleSO <- bam(Event_Coherence_Post_GeoMean ~ 
                   # 1. Surviving Main Power Effects
                   s(RipplePower_Z, k = 5) + 
                   s(SpindlePower_Match_Z, k = 5) + 
                   s(SpindlePower_NonMatch_Z, k = 5) + 
                   ti(SpindlePower_Match_Z, SpindlePower_NonMatch_Z, k = 5) +
                   #te(SpindlePower_Match_Z, SpindlePower_NonMatch_Z, k = c(5, 5))+
                 
                   # 2. The Baseline Phase Landscape (Synergy)
                   te(SOPhase_Match, SOPhase_NonMatch, bs = "cc", k = 8) + 
                   ti(SOPhase_Match, SOPhase_NonMatch, bs = "cc", k = 8) + 
                   s(SOPhase_Match, bs = "cc", k = 8) + 
                   #s(SOPhase_NonMatch, bs = "cc", k = 8) + 
                   
                   # 3. The Winning Interaction (Local Ripple Gate)
                     #ti(RipplePower_Z, SOPhase_Match, bs = c("tp", "cc"), k = c(5, 8)) + 
                   ti(SpindlePower_Match_Z, SOPhase_Match, bs = c("tp", "cc"), k = c(5, 8)) + 
                   #ti(SpindlePower_Match_Z, RipplePower_Z, bs = c("tp", "tp"), k = c(5, 5)) + 
                   
                  
                   # 4. Control Term
                   s(AnimalID, bs = "re") +
                   s(SessionID, bs = "re"), 
                 
                 data = dat_clean, 
                 method = "fREML", 
                 discrete = TRUE, 
                 nthreads = 4)

message("\n--- FINAL MODEL SUMMARY ---")
print(summary(mdl_spindleSO))


message("\n--- MODEL COMPARISON: ANOVA ---")
print(anova(mdl_base, mdl_spindleSO, test = "Chisq"))


message("\n--- MODEL COMPARISON: AIC ---")
aic_results <- AIC(mdl_base, mdl_spindleSO)
aic_results$Delta_AIC <- aic_results$AIC - min(aic_results$AIC)
print(aic_results)


# ==============================================================================
# --- Spindle Ripple interaction
# ==============================================================================
mdl_spindleRipple <- bam(Event_Coherence_Post_GeoMean ~ 
                   # 1. Surviving Main Power Effects
                   s(RipplePower_Z, k = 5) + 
                   s(SpindlePower_Match_Z, k = 5) + 
                   s(SpindlePower_NonMatch_Z, k = 5) + 
                   ti(SpindlePower_Match_Z, SpindlePower_NonMatch_Z, k = 5) +
                   #te(SpindlePower_Match_Z, SpindlePower_NonMatch_Z, k = c(5, 5))+
                   
                   # 2. The Baseline Phase Landscape (Synergy)
                   te(SOPhase_Match, SOPhase_NonMatch, bs = "cc", k = 8) +
                   ti(SOPhase_Match, SOPhase_NonMatch, bs = "cc", k = 8) +
                   s(SOPhase_Match, bs = "cc", k = 8) +
                   #s(SOPhase_NonMatch, bs = "cc", k = 8) + 
                   
                   # 3. The Winning Interaction (Local Ripple Gate)
                     #ti(RipplePower_Z, SOPhase_Match, bs = c("tp", "cc"), k = c(5, 8)) + 
                   #ti(SpindlePower_Match_Z, SOPhase_Match, bs = c("tp", "cc"), k = c(5, 8)) + 
                   ti(SpindlePower_Match_Z, RipplePower_Z, bs = c("tp", "tp"), k = c(5, 5)) + 
                   
                   
                   # 4. Control Term
                   s(AnimalID, bs = "re") +
                   s(SessionID, bs = "re"), 
                 
                 data = dat_clean, 
                 method = "fREML", 
                 discrete = TRUE, 
                 nthreads = 4)

message("\n--- FINAL MODEL SUMMARY ---")
print(summary(mdl_spindleRipple))

message("\n--- MODEL COMPARISON: ANOVA ---")
print(anova(mdl_base, mdl_spindleRipple, test = "Chisq"))


message("\n--- MODEL COMPARISON: AIC ---")
aic_results <- AIC(mdl_base, mdl_spindleRipple)
aic_results$Delta_AIC <- aic_results$AIC - min(aic_results$AIC)
print(aic_results)

