# ==============================================================================
# V1-HC COHERENCE: FINAL PRUNED GAM ANALYSIS SCRIPT 
# ==============================================================================

# --- 1. Load Required Libraries ---
packages <- c("mgcv", "gratia", "ggplot2", "viridis", "dplyr")
for (p in packages) {
  if (!require(p, character.only = TRUE)) install.packages(p)
}

library(mgcv)
library(gratia)
library(ggplot2)
library(viridis)
library(dplyr)

# --- 2. Load and Prepare Data ---
dat <- read.csv("D:/Code/MasaCoherenceDataAnalysis/v1_hc_data.csv")

# Ensure categorical variables are factors
dat$SessionID <- as.factor(dat$SessionID)
dat$AnimalID  <- as.factor(dat$AnimalID)

# --- 3. Clean Data ---
message("Trimming extreme outliers to stabilize multi-dimensional splines...")

dat_clean <- dat 
dat_clean <- dat %>%
  filter(
    RipplePower_Z > -2.6 & RipplePower_Z < 2.6,
    SpindlePower_Match_Z > -2.6 & SpindlePower_Match_Z < 2.6,
   SpindlePower_NonMatch_Z > -2.6 & SpindlePower_NonMatch_Z < 2.6
  )

message(sprintf("Removed %d extreme data points.", nrow(dat) - nrow(dat_clean)))

# --- 4. Fit the Pruned "Minimal Adequate" GAM ---
message("Fitting Pruned Interaction GAM...")

mdl <- bam(Event_Coherence_Post ~ 
                    # 1. Main Power Effects (Spindles kept as negative controls)
                    s(RipplePower_Z, k = 5) + 
                    s(SpindlePower_Match_Z, k = 5) + 
                    s(SpindlePower_NonMatch_Z, k = 5) + 
                    
                    # 2. The Base 2D Phase Landscape
                    te(SOPhase_Match, SOPhase_NonMatch, bs = "cc", k = 8) + 
                    
                    # 3. The ONLY significant 2D Interaction
                    ti(RipplePower_Z, SOPhase_Match, bs = c("tp", "cc"), k = c(5, 8)) +
                    
                    # 4. Random Effects
                    s(SessionID, bs = "re"),
                  data = dat_clean, 
                  method = "fREML", 
                  discrete = TRUE,
                  nthreads = 4)

print(summary(md))

# ==============================================================================
# --- 5. FIGURE GENERATION ---
# ==============================================================================
message("Generating Figures...")

# Exclusions for True Marginals 
# We explicitly turn off power and the 2D ripple/phase interaction to isolate pure phase
terms_to_exclude_1D <- c("s(SessionID)", "s(AnimalID)", "s(RipplePower_Z)", 
                         "s(SpindlePower_Match_Z)", "s(SpindlePower_NonMatch_Z)",
                         "ti(RipplePower_Z,SOPhase_Match)")

# ---------------------------------------------------------
# FIGURE 1: BASE 2D JOINT PHASE HEATMAP (Shifted to center pi)
# ---------------------------------------------------------
sm_2d <- smooth_estimates(mdl, smooth = "te(SOPhase_Match,SOPhase_NonMatch)")

sm_2d_shifted <- sm_2d %>%
  mutate(
    SOPhase_Match = ifelse(SOPhase_Match < 0, SOPhase_Match + 2*pi, SOPhase_Match),
    SOPhase_NonMatch = ifelse(SOPhase_NonMatch < 0, SOPhase_NonMatch + 2*pi, SOPhase_NonMatch)
  ) %>% arrange(SOPhase_Match, SOPhase_NonMatch)

p_heatmap <- ggplot(sm_2d_shifted, aes(x = SOPhase_Match, y = SOPhase_NonMatch, fill = .estimate)) +
  geom_tile(width = (2*pi)/99, height = (2*pi)/99) + 
  geom_contour(aes(z = .estimate), color = "black", alpha = 0.2) + 
  scale_fill_viridis_c(option = "D", name = "Effect") +
  scale_x_continuous(breaks = c(0, pi, 2*pi), labels = c("0", expression(pi), expression(2*pi))) +
  scale_y_continuous(breaks = c(0, pi, 2*pi), labels = c("0", expression(pi), expression(2*pi))) +
  theme_minimal() +
  labs(title = "Figure 1: Baseline Joint Phase Landscape",
       subtitle = "The standalone 2D Phase synergy (averaged across all power levels)",
       x = "Match Phase (Rad)", y = "Non-Match Phase (Rad)")

dev.new(noRStudioGD = TRUE)
print(p_heatmap)

# ---------------------------------------------------------
# FIGURE 2 & 3: TRUE MARGINAL EFFECTS (PHASE)
# ---------------------------------------------------------
grid_full <- data_slice(mdl, 
                        SOPhase_Match = evenly(SOPhase_Match, n = 100),
                        SOPhase_NonMatch = evenly(SOPhase_NonMatch, n = 100))

fit_full <- fitted_values(mdl, data = grid_full, exclude = terms_to_exclude_1D)

# Match Marginal
true_marg_match <- fit_full %>%
  group_by(SOPhase_Match) %>%
  summarise(.fitted = mean(.fitted), .lower_ci = mean(.lower_ci), .upper_ci = mean(.upper_ci))

p_marg_match <- ggplot(true_marg_match, aes(x = SOPhase_Match, y = .fitted)) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), alpha = 0.2, fill = "purple") +
  geom_line(color = "purple", linewidth = 1) + theme_bw() +
  labs(title = "Figure 2: True Marginal Effect (Match Phase)", x = "Phase (Rad)", y = "Effect")

dev.new(noRStudioGD = TRUE)
print(p_marg_match)

# Non-Match Marginal
true_marg_nonmatch <- fit_full %>%
  group_by(SOPhase_NonMatch) %>%
  summarise(.fitted = mean(.fitted), .lower_ci = mean(.lower_ci), .upper_ci = mean(.upper_ci))

p_marg_nonmatch <- ggplot(true_marg_nonmatch, aes(x = SOPhase_NonMatch, y = .fitted)) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), alpha = 0.2, fill = "forestgreen") +
  geom_line(color = "forestgreen", linewidth = 1) + theme_bw() +
  labs(title = "Figure 3: True Marginal Effect (Non-Match Phase)", x = "Phase (Rad)", y = "Effect")

dev.new(noRStudioGD = TRUE)
print(p_marg_nonmatch)

# ---------------------------------------------------------
# FIGURES 4 - 7: 1D MAIN EFFECTS & SIGNIFICANT 2D INTERACTION
# ---------------------------------------------------------
# rug = FALSE completely removes the dots and tick marks for lag-free rendering

dev.new(noRStudioGD = TRUE)
draw(mdl, select = "s(RipplePower_Z)", residuals = FALSE, rug = FALSE) + theme_bw()

dev.new(noRStudioGD = TRUE)
draw(mdl, select = "s(SpindlePower_Match_Z)", residuals = FALSE, rug = FALSE) + theme_bw()

#dev.new(noRStudioGD = TRUE)
#draw(mdl, select = "s(SpindlePower_NonMatch_Z)", residuals = FALSE, rug = FALSE) + theme_bw()

dev.new(noRStudioGD = TRUE)
draw(mdl, select = "ti(RipplePower_Z,SOPhase_Match)", rug = FALSE) + 
  scale_fill_viridis_c(option = "plasma") + theme_minimal() +
  labs(title = "Figure 7: Ripple Power vs Match SO Phase",
       subtitle = "The specific distortion of Match Phase by Ripple Power")


message("All analysis and plotting complete!")
