# ==============================================================================
# DIAGNOSTIC PIPELINE: ANOVA DECOMPOSITION OF JOINT SO PHASE
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
dat <- read.csv("D:/Code/MasaCoherenceDataAnalysis/v1_hc_data_geo3.csv")
dat$SessionID <- as.factor(dat$SessionID)

# --- 3. Clean Data (Fixed Pipe!) ---
z_thresh <- 2.56  

message(sprintf("Trimming extreme outliers beyond +/- %s Z-scores...", z_thresh))
dat_clean <- dat %>%
  filter(
    RipplePower_Z > -z_thresh & RipplePower_Z < z_thresh,
    SpindlePower_Match_Z > -z_thresh & SpindlePower_Match_Z < z_thresh
  ) 
# Note: mutate() is safely commented out, and no hanging pipe above it.

# ==============================================================================
# --- 4. FIT THE DECOMPOSED MODEL ---
# ==============================================================================
message("\nFitting Decomposed Phase Model...")

mdl_decomp <- bam(Event_Coherence_Post ~ 
                    # 1. Main Power Effects
                    s(RipplePower_Z, k = 5) + 
                    s(SpindlePower_Match_Z, k = 5) + 
                    
                    # 2. DECOMPOSED JOINT PHASE
                    # Here we split the te() into two main effects and a ti()
                    s(SOPhase_Match, bs = "cc", k = 8) + 
                    s(SOPhase_NonMatch, bs = "cc", k = 8) + 
                    ti(SOPhase_Match, SOPhase_NonMatch, bs = c("cc", "cc"), k = c(8, 8)) + 
                    
                    # 3. Local Ripple Gate
                    ti(RipplePower_Z, SOPhase_Match, bs = c("tp", "cc"), k = c(5, 8)) + 
                    
                    # 4. Control Term
                    s(SessionID, bs = "re"), 
                  
                  data = dat_clean, 
                  method = "fREML", 
                  discrete = TRUE, 
                  nthreads = 4)

# ==============================================================================
# --- 5. CHECK THE WEIGHTS & P-VALUES ---
# ==============================================================================
message("\n==========================================================")
message("--- DECOMPOSED MODEL SUMMARY ---")
message("Look at the 'Approximate significance of smooth terms' table.")
message("Check the F-value and p-value for:")
message(" 1. s(SOPhase_Match)     <- The independent effect of local phase")
message(" 2. s(SOPhase_NonMatch)  <- The independent effect of distant phase")
message(" 3. ti(SOPhase_Match,..) <- The *pure* synergy between the two")
message("==========================================================\n")

print(summary(mdl_decomp))

# ==============================================================================
# --- 6. VISUALIZE ISOLATED EFFECTS (With Phase-Wrapping & Tight Borders) ---
# ==============================================================================
message("\nExtracting and plotting isolated phase components...")

# ---------------------------------------------------------
# Plot 1: Main Effect of Match Phase
# ---------------------------------------------------------
sm_match <- smooth_estimates(mdl_decomp, smooth = "s(SOPhase_Match)", n = 200) %>%
  mutate(
    # The wrap-around trick!
    SOPhase_Match = ifelse(SOPhase_Match < 0, SOPhase_Match + 2*pi, SOPhase_Match),
    .lower_ci = .estimate - (1.96 * .se),
    .upper_ci = .estimate + (1.96 * .se)
  ) %>% arrange(SOPhase_Match)

dev.new(noRStudioGD = TRUE)
p_match <- ggplot(sm_match, aes(x = SOPhase_Match, y = .estimate)) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), alpha = 0.2, fill = "dodgerblue") +
  geom_line(linewidth = 1, color = "dodgerblue") +
  theme_bw() + 
  scale_x_continuous(breaks = c(0, pi, 2*pi), labels = c("0", expression(pi), expression(2*pi)),
                     limits = c(0, 2*pi), expand = c(0, 0)) +
  labs(title = "Isolated Main Effect: Match SO Phase", 
       x = "Match SO Phase (Rad)", y = "Partial Effect")
print(p_match)

# ---------------------------------------------------------
# Plot 2: Main Effect of Non-Match Phase
# ---------------------------------------------------------
sm_nonmatch <- smooth_estimates(mdl_decomp, smooth = "s(SOPhase_NonMatch)", n = 200) %>%
  mutate(
    # The wrap-around trick!
    SOPhase_NonMatch = ifelse(SOPhase_NonMatch < 0, SOPhase_NonMatch + 2*pi, SOPhase_NonMatch),
    .lower_ci = .estimate - (1.96 * .se),
    .upper_ci = .estimate + (1.96 * .se)
  ) %>% arrange(SOPhase_NonMatch)

dev.new(noRStudioGD = TRUE)
p_nonmatch <- ggplot(sm_nonmatch, aes(x = SOPhase_NonMatch, y = .estimate)) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), alpha = 0.2, fill = "firebrick") +
  geom_line(linewidth = 1, color = "firebrick") +
  theme_bw() + 
  scale_x_continuous(breaks = c(0, pi, 2*pi), labels = c("0", expression(pi), expression(2*pi)),
                     limits = c(0, 2*pi), expand = c(0, 0)) +
  labs(title = "Isolated Main Effect: Non-Match SO Phase", 
       x = "Non-Match SO Phase (Rad)", y = "Partial Effect")
print(p_nonmatch)

# ---------------------------------------------------------
# Plot 3: The "Pure" Phase Interaction Synergy
# ---------------------------------------------------------
sm_syn <- smooth_estimates(mdl_decomp, smooth = "ti(SOPhase_Match,SOPhase_NonMatch)", n = 100) %>%
  mutate(
    # Wrap BOTH phases!
    SOPhase_Match = ifelse(SOPhase_Match < 0, SOPhase_Match + 2*pi, SOPhase_Match),
    SOPhase_NonMatch = ifelse(SOPhase_NonMatch < 0, SOPhase_NonMatch + 2*pi, SOPhase_NonMatch)
  ) %>% arrange(SOPhase_Match, SOPhase_NonMatch)

dev.new(noRStudioGD = TRUE)
p_synergy <- ggplot(sm_syn, aes(x = SOPhase_Match, y = SOPhase_NonMatch, fill = .estimate)) +
  geom_tile(width = (2*pi)/99, height = (2*pi)/99) + 
  geom_contour(aes(z = .estimate), color = "black", alpha = 0.2) + 
  scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "firebrick", midpoint = 0, name = "Interaction") +
  theme_bw() + 
  scale_x_continuous(breaks = c(0, pi, 2*pi), labels = c("0", expression(pi), expression(2*pi)),
                     limits = c(0, 2*pi), expand = c(0, 0)) +
  scale_y_continuous(breaks = c(0, pi, 2*pi), labels = c("0", expression(pi), expression(2*pi)),
                     limits = c(0, 2*pi), expand = c(0, 0)) +
  labs(title = "Pure Phase Synergy (Interaction Only)", 
       x = "Match SO Phase (Rad)", y = "Non-Match SO Phase (Rad)")
print(p_synergy)

# ==============================================================================
# --- 7. RECREATING THE te() LANDSCAPE FROM DECOMPOSED PARTS ---
# ==============================================================================
message("\nReconstructing Total Phase Synergy from Decomposed Terms...")

# 1. Create a clean 2D grid of all possible phase combinations
pred_grid_phase <- expand.grid(
  SOPhase_Match = seq(0, 2*pi, length.out = 100),
  SOPhase_NonMatch = seq(0, 2*pi, length.out = 100),
  # Provide baseline values for the other terms (though type="terms" ignores them anyway)
  RipplePower_Z = 0,
  SpindlePower_Match_Z = 0,
  SessionID = dat_clean$SessionID[1] 
)

# 2. Extract the individual mathematical terms
term_preds_phase <- predict(mdl_decomp, newdata = pred_grid_phase, type = "terms")

# 3. THE MAGIC MATH: Sum the two main effects and the pure interaction
pred_grid_phase$Reconstructed_TE <- 
  term_preds_phase[, "s(SOPhase_Match)"] + 
  term_preds_phase[, "s(SOPhase_NonMatch)"] + 
  term_preds_phase[, "ti(SOPhase_Match,SOPhase_NonMatch)"]

# 4. Plot the Reconstructed Landscape
dev.new(noRStudioGD = TRUE)
p_reconstructed <- ggplot(pred_grid_phase, aes(x = SOPhase_Match, y = SOPhase_NonMatch, fill = Reconstructed_TE)) +
  geom_tile(width = (2*pi)/99, height = (2*pi)/99) + 
  geom_contour(aes(z = Reconstructed_TE), color = "black", alpha = 0.2) + 
  scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "firebrick", midpoint = 0, name = "Total Effect") +
  theme_bw() + 
  scale_x_continuous(breaks = c(0, pi, 2*pi), labels = c("0", expression(pi), expression(2*pi)),
                     limits = c(0, 2*pi), expand = c(0, 0)) +
  scale_y_continuous(breaks = c(0, pi, 2*pi), labels = c("0", expression(pi), expression(2*pi)),
                     limits = c(0, 2*pi), expand = c(0, 0)) +
  labs(title = "Reconstructed Joint Phase Landscape", 
       subtitle = "Mathematically derived from: s(Match) + s(NonMatch) + ti(Match, NonMatch)",
       x = "Match SO Phase (Rad)", y = "Non-Match SO Phase (Rad)")
print(p_reconstructed)