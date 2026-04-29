# ==============================================================================
# V1-HC COHERENCE: COMPLETE GAM ANALYSIS SCRIPT
# ==============================================================================

# --- 1. Load Required Libraries ---
# Automatically install missing packages
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

# --- 3. Clean Data (The "Trumpet" Fix) ---
message("Trimming extreme outliers to stabilize splines...")
dat_clean <- dat
# Filter out extreme power values (beyond +/- 3 Z-scores)
#dat_clean <- dat %>%
#  filter(
#    RipplePower_Z > -3.5 & RipplePower_Z < 3.5,
#    SpindlePower_Match_Z > -3.5 & SpindlePower_Match_Z < 3.5,
#    SpindlePower_NonMatch_Z > -3.5 & SpindlePower_NonMatch_Z < 3.5
#  )

message(sprintf("Removed %d extreme data points.", nrow(dat) - nrow(dat_clean)))

# --- 4. Fit the High-Speed Joint Model ---
message("Fitting Joint Phase GAM on cleaned data...")

# k = 5 for continuous power terms prevents "fewer unique covariate" errors
# nthreads = 4 speeds up the computation
mdl_fast <- bam(Event_Coherence_Post ~ 
                  s(RipplePower_Z, k = 2) + 
                  s(SpindlePower_Match_Z, k = 5) + 
                  s(SpindlePower_NonMatch_Z, k = 5) + 
                  te(SOPhase_Match, SOPhase_NonMatch, bs = "cc", k = 8) + 
                  s(SessionID, bs = "re") +
                  s(AnimalID, bs = "re"),
                data = dat_clean, 
                method = "fREML", 
                discrete = TRUE,
                nthreads = 4)

# Print stats to console
print(summary(mdl_fast))

# ==============================================================================
# --- 5. FIGURE GENERATION ---
# ==============================================================================
message("Generating Figures...")

# ---------------------------------------------------------
# FIGURE 1: SHIFTED 2D JOINT PHASE HEATMAP (Center = pi)
# ---------------------------------------------------------
message("Generating Shifted Heatmap...")

# 1. Extract the raw 2D grid math
sm_2d <- smooth_estimates(mdl_fast, smooth = "te(SOPhase_Match,SOPhase_NonMatch)")

# 2. Shift the phase coordinates and sort them
sm_2d_shifted <- sm_2d %>%
  mutate(
    SOPhase_Match = ifelse(SOPhase_Match < 0, SOPhase_Match + 2*pi, SOPhase_Match),
    SOPhase_NonMatch = ifelse(SOPhase_NonMatch < 0, SOPhase_NonMatch + 2*pi, SOPhase_NonMatch)
  ) %>%
  arrange(SOPhase_Match, SOPhase_NonMatch) # <-- CRITICAL: Sort the new grid

# 3. Plot the shifted heatmap
p_heatmap <- ggplot(sm_2d_shifted, aes(x = SOPhase_Match, y = SOPhase_NonMatch, fill = .estimate)) +
  # <-- CRITICAL: Force the tile width/height so it doesn't collapse to 0
  geom_tile(width = (2*pi)/99, height = (2*pi)/99) + 
  geom_contour(aes(z = .estimate), color = "black", alpha = 0.2) + 
  scale_fill_viridis_c(option = "D", name = "Effect") +
  scale_x_continuous(breaks = c(0, pi, 2*pi), labels = c("0", expression(pi), expression(2*pi))) +
  scale_y_continuous(breaks = c(0, pi, 2*pi), labels = c("0", expression(pi), expression(2*pi))) +
  theme_minimal() +
  labs(title = "Figure 1: Joint Phase Landscape",
       subtitle = "Continuous 2D Synergy (Shifted to center \u03c0)",
       x = "Match Phase (Rad)", 
       y = "Non-Match Phase (Rad)")

dev.new(noRStudioGD = TRUE)
print(p_heatmap)

# ---------------------------------------------------------
# FIGURE 2: TRUE MARGINAL EFFECT (MATCH PHASE)
# ---------------------------------------------------------
message("Calculating true marginal for Match Phase...")

# 1. Create the full 2D grid (100x100)
grid_full <- data_slice(mdl_fast, 
                        SOPhase_Match = evenly(SOPhase_Match, n = 100),
                        SOPhase_NonMatch = evenly(SOPhase_NonMatch, n = 100))

# 2. Predict the entire 2D surface
fit_full <- fitted_values(mdl_fast, 
                          data = grid_full, 
                          exclude = c("s(SessionID)", "s(AnimalID)", "s(RipplePower_Z)", 
                                      "s(SpindlePower_Match_Z)", "s(SpindlePower_NonMatch_Z)"))

# 3. Collapse (Average) across the Non-Match dimension using dplyr
true_marg_match <- fit_full %>%
  group_by(SOPhase_Match) %>%
  summarise(.fitted = mean(.fitted),
            .lower_ci = mean(.lower_ci),
            .upper_ci = mean(.upper_ci))

# 4. Plot
p_marg_match <- ggplot(true_marg_match, aes(x = SOPhase_Match, y = .fitted)) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), alpha = 0.2, fill = "purple") +
  geom_line(color = "purple", linewidth = 1) +
  theme_bw() +
  labs(title = "Figure 2: True Marginal Effect (Match Phase)", 
       subtitle = "Averaged across all Non-Match phases",
       x = "Match Phase (Rad)", y = "Overall Partial Effect")

dev.new(noRStudioGD = TRUE)
print(p_marg_match)

# ---------------------------------------------------------
# FIGURE 3: TRUE MARGINAL EFFECT (NON-MATCH PHASE)
# ---------------------------------------------------------
message("Calculating true marginal for Non-Match Phase...")

# 1. We already have the full 2D prediction (fit_full), so we just collapse the other way!
true_marg_nonmatch <- fit_full %>%
  group_by(SOPhase_NonMatch) %>%
  summarise(.fitted = mean(.fitted),
            .lower_ci = mean(.lower_ci),
            .upper_ci = mean(.upper_ci))

# 2. Plot
p_marg_nonmatch <- ggplot(true_marg_nonmatch, aes(x = SOPhase_NonMatch, y = .fitted)) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), alpha = 0.2, fill = "forestgreen") +
  geom_line(color = "forestgreen", linewidth = 1) +
  theme_bw() +
  labs(title = "Figure 3: True Marginal Effect (Non-Match Phase)", 
       subtitle = "Averaged across all Match phases",
       x = "Non-Match Phase (Rad)", y = "Overall Partial Effect")

dev.new(noRStudioGD = TRUE)
print(p_marg_nonmatch)

# ---------------------------------------------------------
# FIGURES 4, 5 & 6: POWER EFFECTS
# ---------------------------------------------------------
# residuals = FALSE removes the scatter dots
# rug = FALSE removes the tick marks on the x-axis

dev.new(noRStudioGD = TRUE)
draw(mdl_fast, select = "s(RipplePower_Z)", residuals = FALSE, rug = FALSE) + 
  theme_bw() + labs(title = "Figure 4: Ripple Power Effect")

dev.new(noRStudioGD = TRUE)
draw(mdl_fast, select = "s(SpindlePower_Match_Z)", residuals = FALSE, rug = FALSE) + 
  theme_bw() + labs(title = "Figure 5: Matched Spindle Power Effect")

dev.new(noRStudioGD = TRUE)
draw(mdl_fast, select = "s(SpindlePower_NonMatch_Z)", residuals = FALSE, rug = FALSE) + 
  theme_bw() + labs(title = "Figure 6: Non-Matched Spindle Power Effect")
message("All figures plotted successfully. Check your taskbar!")