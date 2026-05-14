# ==============================================================================
# V1-HC COHERENCE: STAGE 1 - REDUNDANCY & INDEPENDENCE AUDIT (V2)
# ==============================================================================

# --- 1. Load Required Libraries ---
packages <- c("mgcv", "gratia", "ggplot2", "dplyr", "corrplot")
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
dat$AnimalID  <- as.factor(dat$AnimalID)

# --- 3. Clean Data (Configurable Threshold) ---
z_thresh <- 2.56  # Configurable parameter for outlier exclusion

message(sprintf("Trimming extreme outliers beyond +/- %s Z-scores...", z_thresh))
dat_clean <- dat %>%
  filter(
    RipplePower_Z > -z_thresh & RipplePower_Z < z_thresh,
    SpindlePower_Match_Z > -z_thresh & SpindlePower_Match_Z < z_thresh,
    SpindlePower_NonMatch_Z > -z_thresh & SpindlePower_NonMatch_Z < z_thresh
  )

message(sprintf("Removed %d extreme data points.", nrow(dat) - nrow(dat_clean)))

# ==============================================================================
# --- TASK 1.1: COLLINEARITY (Spearman Correlation) ---
# ==============================================================================
message("\n--- Task 1.1: Spearman Correlation Matrix ---")
power_vars <- dat_clean %>% 
  select(RipplePower_Z, SpindlePower_Match_Z, SpindlePower_NonMatch_Z)

cor_matrix <- cor(power_vars, method = "spearman")
print(round(cor_matrix, 3))

# ==============================================================================
# --- TASK 1.3: PARTIAL EFFECTS (Fit Model) ---
# ==============================================================================
message("\nFitting Stage 1 Model with independent 1D terms and Random Effects...")

mdl <- bam(Event_Coherence_Post ~ 
             # Power main effects
             s(RipplePower_Z, k = 5) + 
             s(SpindlePower_Match_Z, k = 5) + 
             s(SpindlePower_NonMatch_Z, k = 5) + 
             # Independent 1D phase terms
             s(SOPhase_Match, bs = "cc", k = 8) + 
             s(SOPhase_NonMatch, bs = "cc", k = 8) +
             # Random effects
             s(SessionID, bs = "re"),
           data = dat_clean, 
           method = "fREML", 
           discrete = TRUE,
           nthreads = 4)

message("\n--- Task 1.3: Main Effects Summary ---")
print(summary(mdl))

# ==============================================================================
# --- TASK 1.2: CONCURVITY ---
# ==============================================================================
message("\n--- Task 1.2: Concurvity Assessment ---")
concurvity_results <- concurvity(mdl, full = FALSE)
print(round(concurvity_results$estimate, 3))

# ==============================================================================
# --- STAGE 1 VISUALIZATIONS: MAIN EFFECTS ---
# ==============================================================================
message("Generating Partial Effect Plots for all 1D terms...")

# Plotting the Power Terms
dev.new(noRStudioGD = TRUE)
p1 <- draw(mdl, select = "s(RipplePower_Z)", residuals = FALSE, rug = FALSE) + 
  theme_bw() + labs(title = "Stage 1: Ripple Power")
print(p1)

dev.new(noRStudioGD = TRUE)
p2 <- draw(mdl, select = "s(SpindlePower_Match_Z)", residuals = FALSE, rug = FALSE) + 
  theme_bw() + labs(title = "Stage 1: Match Spindle Power")
print(p2)

dev.new(noRStudioGD = TRUE)
p3 <- draw(mdl, select = "s(SpindlePower_NonMatch_Z)", residuals = FALSE, rug = FALSE) + 
  theme_bw() + labs(title = "Stage 1: Non-Match Spindle Power")
print(p3)

# Plotting the independent Phase Terms
dev.new(noRStudioGD = TRUE)
p4 <- draw(mdl, select = "s(SOPhase_Match)", residuals = FALSE, rug = FALSE) + 
  scale_x_continuous(breaks = c(0, pi, 2*pi), labels = c("0", expression(pi), expression(2*pi))) +
  theme_bw() + labs(title = "Stage 1: Match SO Phase (1D)")
print(p4)

dev.new(noRStudioGD = TRUE)
p5 <- draw(mdl, select = "s(SOPhase_NonMatch)", residuals = FALSE, rug = FALSE) + 
  scale_x_continuous(breaks = c(0, pi, 2*pi), labels = c("0", expression(pi), expression(2*pi))) +
  theme_bw() + labs(title = "Stage 1: Non-Match SO Phase (1D)")
print(p5)

message("Stage 1 execution complete.")