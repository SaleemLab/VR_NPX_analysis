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
dat <- read.csv("C:/Users/masah/Documents/GitHub/VR_NPX_analysis/predicting_coherence_with_GAMs/v1_hc_data_geo_POST.csv")

dat <- read.csv("C:/Users/masah/Documents/GitHub/VR_NPX_analysis/predicting_coherence_with_GAMs/v1_hc_data_geo_0_100ms_POST.csv")

dat <- read.csv("C:/Users/masah/Documents/GitHub/VR_NPX_analysis/predicting_coherence_with_GAMs/v1_hc_data_geo_100_200ms_POST.csv")
dat <- read.csv("C:/Users/masah/Documents/GitHub/VR_NPX_analysis/predicting_coherence_with_GAMs/v1_hc_data_geo_0_200ms_POST.csv")

dat$SessionID <- as.factor(dat$SessionID)
dat$AnimalID <- as.factor(dat$AnimalID)

# --- 3. Clean Data ---
z_thresh <- 2.56  # include <99th centiles of data

message(sprintf("Trimming extreme outliers beyond +/- %s Z-scores...", z_thresh))
dat_clean <- dat %>%
  filter(
    RipplePower_Z > -z_thresh & RipplePower_Z < z_thresh,
    SpindlePower_Match_Z > -z_thresh & SpindlePower_Match_Z < z_thresh,
    SpindlePower_NonMatch_Z > -z_thresh & SpindlePower_NonMatch_Z < z_thresh,
    #SpindlePowerPRE_Match_Z > -z_thresh & SpindlePowerPRE_Match_Z < z_thresh,
    #SpindlePowerPRE_NonMatch_Z > -z_thresh & SpindlePowerPRE_NonMatch_Z < z_thresh,
    #SpindlePowerPOST_Match_Z > -z_thresh & SpindlePowerPOST_Match_Z < z_thresh,
    #SpindlePowerPOST_NonMatch_Z > -z_thresh & SpindlePowerPOST_NonMatch_Z < z_thresh
  )


library(ggplot2)

ggplot(dat_clean, aes(x = SpindlePowerPOST_Match_Z)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "white") +
  theme_minimal() +
  labs(title = "Distribution of Spindle Power POST (Match)",
       x = "Z-Score",
       y = "Count")

# Standardize the existing Coherence column from the CSV

# Standardize the Dependent Variable (Coherence)
# This makes 'Partial Effect' units equal to Standard Deviations of Coherence
#mutate(Event_Coherence_Post_Z = as.numeric(scale(Event_Coherence_Post)))

# ==============================================================================
# --- TASK 5.1: FIT THE MINIMAL ADEQUATE MODEL ---
# ==============================================================================
message("\nFitting the Final Minimal Adequate Model...")

mdl_final <- bam(Event_Coherence_Post_GeoMean ~ 
                   # 1. Surviving Main Power Effects
                   s(RipplePower_Z, k = 5) + 
                   s(SpindlePowerPOST_Match_Z, k = 5) + 
                   s(SpindlePowerPOST_NonMatch_Z, k = 5) + 
                   #ti(SpindlePowerPOST_Match_Z, SpindlePowerPOST_NonMatch_Z, k = 5) +
                   
                   # 2. The Baseline Phase Landscape (Synergy)
                   #te(SOPhasePOST_Match, SOPhasePOST_NonMatch, bs = "cc", k = 8) + 
                   ti(SOPhasePOST_Match, SOPhasePOST_NonMatch, bs = "cc", k = 8) + 
                   s(SOPhasePOST_Match, bs = "cc", k = 8) + 
                   s(SOPhasePOST_NonMatch, bs = "cc", k = 8) + 
                   # 3. The Winning Interaction (Local Ripple Gate)
                   ti(SpindlePowerPOST_Match_Z, SOPhasePOST_Match, bs = c("tp", "cc"), k = c(5, 8)) + 
                   ti(RipplePower_Z, SOPhasePOST_Match, bs = c("tp", "cc"), k = c(5, 8)) + 
                   
                   # 4. Control Term
                   s(AnimalID, bs = "re") +
                   s(SessionID, bs = "re"), 
                 
                 data = dat_clean, 
                 method = "fREML", 
                 discrete = TRUE, 
                 nthreads = 4)

message("\n--- FINAL MODEL SUMMARY ---")
print(summary(mdl_final))


# ==============================================================================
# --- TASK 5.2: FORMAL DIAGNOSTICS (GAM.CHECK) ---
# ==============================================================================
message("\n--- RUNNING DIAGNOSTICS (gam.check) ---")
# gam.check(mdl_final)

# ==============================================================================
# --- STAGE 5 VISUALIZATIONS: PUBLICATION FIGURES ---
# ==============================================================================
message("\nGenerating Final Publication Figures...")



# 1. Main Effect: Ripple Power
# Changed file to "post" and Title to "Post"
cairo_pdf("ripple_power_post_z.pdf", width = 4.3, height = 4.3)
p_ripple <- draw(mdl_final, select = "s(RipplePower_Z)", residuals = FALSE, rug = FALSE) + 
  theme_bw(base_family = "Arial") + 
  theme(aspect.ratio = 1) +
  labs(title = "Main Effect: Ripple Power (Post)", x = "Ripple Power (Z)", y = "Partial Effect")

print(p_ripple)
dev.off()

# Calculate scaling factors
raw_breaks <- c(5, 10, 15)
z_breaks <- sapply(raw_breaks, function(val) {
  dat_clean$RipplePower_Z[which.min(abs(dat_clean$RipplePower - val))]
})

cairo_pdf("ripple_power_post_z_raw.pdf", width = 4.3, height = 4.3)
p_ripple_raw <- draw(mdl_final, select = "s(RipplePower_Z)", residuals = FALSE, rug = FALSE) + 
  theme_bw(base_family = "Arial") + 
  theme(aspect.ratio = 1) +
  scale_x_continuous(breaks = z_breaks, labels = raw_breaks) +
  labs(
    title = "Ripple Power and Post coherence", 
    x = "Ripple Power (Raw Units)", 
    y = "Partial Effect"
  )
dev.new(noRStudioGD = TRUE)
print(p_ripple_raw)
dev.off()


# 2. Main Effect: Spindle Power
p_spindle <- draw(mdl_final, select = "s(SpindlePowerPOST_NonMatch_Z)", residuals = FALSE, rug = FALSE) + 
  theme_bw(base_family = "Arial") + 
  theme(aspect.ratio = 1) +
  labs(title = "Main Effect: Spindle Power NonMatched (POST)", x = "Spindle Power NonMatched (Z)", y = "Partial Effect")
dev.new(noRStudioGD = TRUE)
print(p_spindle)



p_spindle <- draw(mdl_final, select = "s(SpindlePowerPOST_Match_Z)", residuals = FALSE, rug = FALSE) + 
  theme_bw(base_family = "Arial") + 
  theme(aspect.ratio = 1) +
  labs(title = "Main Effect: Spindle Power Matched (Post)", x = "Spindle Power Matched (Z)", y = "Partial Effect")

dev.new(noRStudioGD = TRUE)
print(p_spindle)

cairo_pdf("spindle_power_matched_post_z.pdf", width = 4.3, height = 4.3)
print(p_spindle)
dev.off()

# Calculate scaling factors
raw_breaks_spindle <- c(-1, 1, 3)
z_breaks_spindle <- sapply(raw_breaks_spindle, function(val) {
  dat_clean$SpindlePower_Match_Z[which.min(abs(dat_clean$SpindlePower_Match - val))]
})

cairo_pdf("spindle_power_matched_post_raw_z.pdf", width = 4.3, height = 4.3)
p_spindle_raw <- draw(mdl_final, select = "s(SpindlePowerPOST_Match_Z)", residuals = FALSE, rug = FALSE) + 
  theme_bw(base_family = "Arial") + 
  theme(aspect.ratio = 1) +
  scale_x_continuous(breaks = z_breaks_spindle, labels = raw_breaks_spindle) +
  labs(
    title = "Spindle Power Matched and Post coherence", 
    x = "Spindle Power (Raw Units)", 
    y = "Partial Effect"
  )

print(p_spindle_raw)
dev.off()


# 3. The 2D Joint Phase Landscape
message("Rendering Joint Phase Landscape...")
sm_2d <- smooth_estimates(mdl_final, smooth = "te(SOPhasePOST_Match,SOPhasePOST_NonMatch)", n = 100) %>%
  mutate(
    .lower_ci = .estimate - (1.96 * .se),
    .upper_ci = .estimate + (1.96 * .se)
  )

p_joint <- sm_2d %>%
  mutate(
    SOPhase_Match = ifelse(SOPhasePOST_Match < 0, SOPhasePOST_Match + 2*pi, SOPhasePOST_Match),
    SOPhasePOST_NonMatch = ifelse(SOPhasePOST_NonMatch < 0, SOPhasePOST_NonMatch + 2*pi, SOPhasePOST_NonMatch)
  ) %>%
  arrange(SOPhasePOST_Match, SOPhasePOST_NonMatch) %>%
  ggplot(aes(x = SOPhasePOST_Match, y = SOPhasePOST_NonMatch, fill = .estimate)) +
  geom_tile(width = (2*pi)/99, height = (2*pi)/99) + 
  geom_contour(aes(z = .estimate), color = "black", alpha = 0.2) + 
  scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "firebrick", midpoint = 0, name = "Effect") +
  scale_x_continuous(breaks = c(0, pi, 2*pi), labels = c("0", expression(pi), expression(2*pi))) +
  scale_y_continuous(breaks = c(0, pi, 2*pi), labels = c("0", expression(pi), expression(2*pi))) +
  theme_minimal(base_family = "Arial") +
  theme(aspect.ratio = 1) +
  labs(title = "Inter-hemispheric Phase effect on Post coherence",
       x = "Match SO Phase (Rad)", y = "Non-Match SO Phase (Rad)")
dev.new(noRStudioGD = TRUE)
print(p_joint)

# Increased width slightly for the legend, but kept aspect ratio square
cairo_pdf("SO_phase_post_z.pdf", width = 5.2, height = 4.3)
print(p_joint)
dev.off()



# ==============================================================================
# VISUALIZING THE TOTAL SPINDLE EFFECT (Main Effect + ti() Gate) - POST
# ==============================================================================
# ==============================================================================
# TOTAL SPINDLE EFFECT (Main + ti) - POST
# ==============================================================================

# 1. Prediction Grid
pred_grid_spindle_post <- expand.grid(
  SpindlePowerPOST_Match_Z = seq(-2.5, 2.5, length.out = 100),
  SOPhasePOST_Match = seq(0, 2*pi, length.out = 100),
  
  # Hold others constant
  SpindlePowerPOST_NonMatch_Z = 0,
  RipplePower_Z = 0,
  SOPhasePOST_NonMatch = pi,
  SessionID = dat_clean$SessionID[1],
  AnimalID = dat_clean$AnimalID[1] 
)

# 2. Extract and Sum
term_preds_spindle <- predict(mdl_final, newdata = pred_grid_spindle_post, type = "terms")

pred_grid_spindle_post$Reconstructed_TE <- 
  term_preds_spindle[, "s(SpindlePowerPOST_Match_Z)"] + 
  term_preds_spindle[, "ti(SpindlePowerPOST_Match_Z,SOPhasePOST_Match)"]

# 3. Plot
p_spindle_total_post <- pred_grid_spindle_post %>%
  arrange(SOPhasePOST_Match, SpindlePowerPOST_Match_Z) %>%
  ggplot(aes(x = SOPhasePOST_Match, y = SpindlePowerPOST_Match_Z, fill = Reconstructed_TE)) +
  geom_tile(width = (2*pi)/99, height = 5/99) + 
  geom_contour(aes(z = Reconstructed_TE), color = "black", alpha = 0.2) + 
  scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "firebrick", midpoint = 0, name = "Effect") +
  scale_x_continuous(breaks = c(0, pi, 2*pi), labels = c("0", expression(pi), expression(2*pi)), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_minimal(base_family = "Arial") +
  theme(aspect.ratio = 1,
        axis.ticks = element_line(color = "grey80"),
        axis.ticks.length = unit(0.2, "cm"),
        axis.line = element_line(color = "grey80"),
        panel.grid.minor = element_blank()) +
  labs(title = "Total Effect of Spindle Power by SO Phase (POST)", 
       subtitle = "Reconstructed: s(Spindle) + ti(Spindle, Phase)",
       x = "Match SO Phase POST (Rad)", y = "Match Spindle Power POST (Z)")
dev.new(noRStudioGD = TRUE)
print(p_spindle_total_post)


# ==============================================================================
# INTERACTION ONLY (ti Gate) - POST
# ==============================================================================

# Use the same prediction grid as above
pred_grid_spindle_post$Interaction_Only <- 
  term_preds_spindle[, "ti(SpindlePowerPOST_Match_Z,SOPhasePOST_Match)"]

# Plotting
p_spindle_inter_post <- pred_grid_spindle_post %>%
  arrange(SOPhasePOST_Match, SpindlePowerPOST_Match_Z) %>%
  ggplot(aes(x = SOPhasePOST_Match, y = SpindlePowerPOST_Match_Z, fill = Interaction_Only)) +
  geom_tile(width = (2*pi)/99, height = 5/99) + 
  geom_contour(aes(z = Interaction_Only), color = "black", alpha = 0.2) + 
  scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "firebrick", midpoint = 0, name = "Effect") +
  scale_x_continuous(breaks = c(0, pi, 2*pi), labels = c("0", expression(pi), expression(2*pi)), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_minimal(base_family = "Arial") +
  theme(aspect.ratio = 1,
        axis.ticks = element_line(color = "grey80"),
        axis.ticks.length = unit(0.2, "cm"),
        axis.line = element_line(color = "grey80"),
        panel.grid.minor = element_blank()) +
  labs(title = "Interaction Effect of Spindle Power by SO Phase (POST)", 
       subtitle = "Term: ti(Spindle, Phase) Only",
       x = "Match SO Phase POST (Rad)", y = "Match Spindle Power POST (Z)")
dev.new(noRStudioGD = TRUE)
print(p_spindle_inter_post)


# ==============================================================================
# TOTAL SPINDLE EFFECT (Main + Interaction raw) - POST
# ==============================================================================

# 0. Get conversion parameters from your data
mean_spin <- mean(dat_clean$SpindlePowerPOST_Match, na.rm = TRUE)
sd_spin   <- sd(dat_clean$SpindlePowerPOST_Match, na.rm = TRUE)

# Define your desired RAW range
raw_min_spin <- -2
raw_max_spin <- 4

# Convert those raw limits to Z-scores for the prediction grid
z_min_spin <- (raw_min_spin - mean_spin) / sd_spin
z_max_spin <- (raw_max_spin - mean_spin) / sd_spin

# 1. Create the prediction grid using the converted Z-limits
pred_grid_spindle_post <- expand.grid(
  SpindlePowerPOST_Match_Z = seq(z_min_spin, z_max_spin, length.out = 100),
  SOPhasePOST_Match = seq(0, 2*pi, length.out = 100),
  
  # Hold other variables at baseline
  SpindlePowerPOST_NonMatch_Z = 0,
  RipplePower_Z = 0,
  SOPhasePOST_NonMatch = pi,
  SessionID = dat_clean$SessionID[1],
  AnimalID = dat_clean$AnimalID[1] 
)

# 2. Extract and Sum terms
term_preds_spindle <- predict(mdl_final, newdata = pred_grid_spindle_post, type = "terms")

pred_grid_spindle_post$Total_Effect <- 
  term_preds_spindle[, "s(SpindlePowerPOST_Match_Z)"] + 
  term_preds_spindle[, "ti(SpindlePowerPOST_Match_Z,SOPhasePOST_Match)"]

# 3. Prepare Axis Ticks (Mapping Z-score positions to Raw labels)
# Create a few labels between -1.5 and 3.5
raw_axis_labels <- seq(-2, 4, by = 1)
z_axis_breaks   <- (raw_axis_labels - mean_spin) / sd_spin

# 4. Plotting
p_spindle_total_raw <- pred_grid_spindle_post %>%
  arrange(SOPhasePOST_Match, SpindlePowerPOST_Match_Z) %>%
  ggplot(aes(x = SOPhasePOST_Match, y = SpindlePowerPOST_Match_Z, fill = Total_Effect)) +
  
  # Height is now based on the Z-range we calculated
  geom_tile(width = (2*pi)/99, height = (z_max_spin - z_min_spin)/99) + 
  geom_contour(aes(z = Total_Effect), color = "black", alpha = 0.2) + 
  
  scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "firebrick", 
                       midpoint = 0, name = "Effect") +
  
  scale_x_continuous(breaks = c(0, pi, 2*pi), 
                     labels = c("0", expression(pi), expression(2*pi)), 
                     expand = c(0,0)) +
  
  # Map the Y-axis back to Raw labels
  scale_y_continuous(breaks = z_axis_breaks, 
                     labels = raw_axis_labels, 
                     expand = c(0,0)) +
  
  theme_minimal(base_family = "Arial") +
  theme(
    aspect.ratio = 1,
    axis.ticks = element_line(color = "grey80"),
    axis.ticks.length = unit(0.2, "cm"),
    axis.line = element_line(color = "grey80"),
    panel.grid.minor = element_blank()
  ) +
  labs(title = "Total Effect of Spindle Power by SO Phase (POST)", 
       subtitle = "Reconstructed: s(Spindle) + ti(Spindle, Phase)",
       x = "Match SO Phase POST (Rad)", 
       y = "Match Spindle Power POST (Raw Units)")
dev.new(noRStudioGD = TRUE)
print(p_spindle_total_raw)


# ==============================================================================
# SPINDLE INTERACTION ONLY (ti Gate raw) - POST
# ==============================================================================
# 0. Conversion Parameters (using same mean/sd as Total Effect)
mean_spin <- mean(dat_clean$SpindlePowerPOST_Match, na.rm = TRUE)
sd_spin   <- sd(dat_clean$SpindlePowerPOST_Match, na.rm = TRUE)

# New Raw boundaries
raw_min_inter <- -2
raw_max_inter <- 4

# Convert to Z-scores for prediction grid
z_min_inter <- (raw_min_inter - mean_spin) / sd_spin
z_max_inter <- (raw_max_inter - mean_spin) / sd_spin

# 1. Create Grid
pred_grid_inter_post <- expand.grid(
  SpindlePowerPOST_Match_Z = seq(z_min_inter, z_max_inter, length.out = 100),
  SOPhasePOST_Match = seq(0, 2*pi, length.out = 100),
  
  # Hold others constant
  SpindlePowerPOST_NonMatch_Z = 0,
  RipplePower_Z = 0,
  SOPhasePOST_NonMatch = pi,
  SessionID = dat_clean$SessionID[1],
  AnimalID = dat_clean$AnimalID[1] 
)

# 2. Extract ti() term only
term_preds_inter <- predict(mdl_final, newdata = pred_grid_inter_post, type = "terms")
pred_grid_inter_post$Interaction_Only <- 
  term_preds_inter[, "ti(SpindlePowerPOST_Match_Z,SOPhasePOST_Match)"]

# 3. Axis Setup
raw_axis_inter <- seq(-2, 4, by = 1)
z_axis_inter   <- (raw_axis_inter - mean_spin) / sd_spin

# 4. Plotting
p_spindle_inter_raw <- pred_grid_inter_post %>%
  arrange(SOPhasePOST_Match, SpindlePowerPOST_Match_Z) %>%
  ggplot(aes(x = SOPhasePOST_Match, y = SpindlePowerPOST_Match_Z, fill = Interaction_Only)) +
  
  # Height is dynamic: (z_max - z_min) / 99
  geom_tile(width = (2*pi)/99, height = (z_max_inter - z_min_inter)/99) + 
  geom_contour(aes(z = Interaction_Only), color = "black", alpha = 0.2) + 
  
  scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "firebrick", 
                       midpoint = 0, name = "Effect") +
  
  scale_x_continuous(breaks = c(0, pi, 2*pi), 
                     labels = c("0", expression(pi), expression(2*pi)), 
                     expand = c(0,0)) +
  
  scale_y_continuous(breaks = z_axis_inter, 
                     labels = raw_axis_inter, 
                     expand = c(0,0)) +
  
  theme_minimal(base_family = "Arial") +
  theme(
    aspect.ratio = 1,
    axis.ticks = element_line(color = "grey80"),
    axis.ticks.length = unit(0.2, "cm"),
    axis.line = element_line(color = "grey80"),
    panel.grid.minor = element_blank()
  ) +
  labs(title = "Interaction Effect of Spindle Power by SO Phase (POST)", 
       subtitle = "Term: ti(Spindle, Phase) Only | Raw Range [-2, 4]",
       x = "Match SO Phase POST (Rad)", 
       y = "Match Spindle Power POST (Raw Units)")
dev.new(noRStudioGD = TRUE)
print(p_spindle_inter_raw)





# ==============================================================================
# TOTAL RIPPLE EFFECT (Main + Interaction) - POST
# ==============================================================================

# 0. Parameters for Raw Y-axis conversion
mean_rip <- mean(dat_clean$RipplePower, na.rm = TRUE)
sd_rip   <- sd(dat_clean$RipplePower, na.rm = TRUE)
raw_min_rip <- 5
raw_max_rip <- 16
z_min_rip <- (raw_min_rip - mean_rip) / sd_rip
z_max_rip <- (raw_max_rip - mean_rip) / sd_rip

# 1. Create the prediction grid
pred_grid_ripple_post <- expand.grid(
  RipplePower_Z = seq(z_min_rip, z_max_rip, length.out = 100),
  SOPhasePOST_Match = seq(0, 2*pi, length.out = 100),
  SpindlePowerPOST_Match_Z = 0,
  SpindlePowerPOST_NonMatch_Z = 0,
  SOPhasePOST_NonMatch = pi,
  SessionID = dat_clean$SessionID[1],
  AnimalID = dat_clean$AnimalID[1]
)

# 2. Extract terms
term_preds_ripple <- predict(mdl_final, newdata = pred_grid_ripple_post, type = "terms")

# 3. Sum: Main Effect s(Ripple) + Interaction ti(Ripple, Phase)
pred_grid_ripple_post$Total_Effect <- 
  term_preds_ripple[, "s(RipplePower_Z)"] + 
  term_preds_ripple[, "ti(RipplePower_Z,SOPhasePOST_Match)"]

# 4. Prepare axis breaks
raw_axis_rip <- seq(5, 15, by = 5)
z_axis_rip   <- (raw_axis_rip - mean_rip) / sd_rip

# 5. Plot Total Effect
p_ripple_total_post <- pred_grid_ripple_post %>%
  arrange(SOPhasePOST_Match, RipplePower_Z) %>%
  ggplot(aes(x = SOPhasePOST_Match, y = RipplePower_Z, fill = Total_Effect)) +
  geom_tile(width = (2*pi)/99, height = (z_max_rip - z_min_rip)/99) + 
  geom_contour(aes(z = Total_Effect), color = "black", alpha = 0.2) + 
  scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "firebrick", midpoint = 0, name = "Effect") +
  scale_x_continuous(breaks = c(0, pi, 2*pi), labels = c("0", expression(pi), expression(2*pi)), expand = c(0,0)) +
  scale_y_continuous(breaks = z_axis_rip, labels = raw_axis_rip, expand = c(0,0)) +
  theme_minimal(base_family = "Arial") +
  theme(aspect.ratio = 1,
        axis.ticks = element_line(color = "grey80"),
        axis.ticks.length = unit(0.2, "cm"),
        axis.line = element_line(color = "grey80"),
        panel.grid.minor = element_blank()) +
  labs(title = "Total Effect of Ripple Power by SO Phase (POST)", 
       subtitle = "Reconstructed: s(Ripple) + ti(Ripple, Phase)",
       x = "Match SO Phase POST (Rad)", y = "Ripple Power (Raw Units)")
dev.new(noRStudioGD = TRUE)
print(p_ripple_total_post)



# ==============================================================================
# RIPPLE INTERACTION ONLY (ti Gate) - POST
# ==============================================================================

# Use the extracted terms from above, but only take the ti() column
pred_grid_ripple_post$Interaction_Only <- 
  term_preds_ripple[, "ti(RipplePower_Z,SOPhasePOST_Match)"]

# Plot Interaction Only
p_ripple_inter_post <- pred_grid_ripple_post %>%
  arrange(SOPhasePOST_Match, RipplePower_Z) %>%
  ggplot(aes(x = SOPhasePOST_Match, y = RipplePower_Z, fill = Interaction_Only)) +
  geom_tile(width = (2*pi)/99, height = (z_max_rip - z_min_rip)/99) + 
  geom_contour(aes(z = Interaction_Only), color = "black", alpha = 0.2) + 
  scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "firebrick", midpoint = 0, name = "Effect") +
  scale_x_continuous(breaks = c(0, pi, 2*pi), labels = c("0", expression(pi), expression(2*pi)), expand = c(0,0)) +
  scale_y_continuous(breaks = z_axis_rip, labels = raw_axis_rip, expand = c(0,0)) +
  theme_minimal(base_family = "Arial") +
  theme(aspect.ratio = 1,
        axis.ticks = element_line(color = "grey80"),
        axis.ticks.length = unit(0.2, "cm"),
        axis.line = element_line(color = "grey80"),
        panel.grid.minor = element_blank()) +
  labs(title = "Interaction Effect of Ripple Power by SO Phase (POST)", 
       subtitle = "Term: ti(Ripple, Phase) Only",
       x = "Match SO Phase POST (Rad)", y = "Ripple Power (Raw Units)")
dev.new(noRStudioGD = TRUE)
print(p_ripple_inter_post)



# ==============================================================================
# --- 3. MULTI-METRIC EFFECT SIZE CALCULATIONS ---
# ==============================================================================
message("\nCalculating 4 Effect Size Metrics (This may take a minute)...")

# EXACT formula components to rebuild the models robustly
formula_terms <- c(
  "s(RipplePower_Z, k = 5)", 
  "s(SpindlePower_Match_Z, k = 5)", 
  "te(SOPhase_Match, SOPhase_NonMatch, bs = 'cc', k = 8)"
)

# EXACT labels output by summary() and smooth_estimates()
smooth_labels <- c(
  "s(RipplePower_Z)", 
  "s(SpindlePower_Match_Z)", 
  "te(SOPhase_Match,SOPhase_NonMatch)"
)

# Initialize results dataframe
effect_sizes <- data.frame(Term = smooth_labels, Partial_Deviance = NA, 
                           Eta_Sq_Partial = NA, Amplitude = NA, RMS = NA)

# 3A. Deviance Partitioning (Drop-One)
full_dev <- summary(mdl_final)$dev.expl * 100

message("  -> Running Drop-One Deviance Models:")
for(i in 1:length(formula_terms)) {
  message(sprintf("     Fitting reduced model without: %s", smooth_labels[i]))
  
  # Rebuild formula safely by taking all terms EXCEPT the current one, plus SessionID
  active_terms <- c(formula_terms[-i], "s(SessionID, bs = 're')")
  red_form <- as.formula(paste("Event_Coherence_Post_GeoMean ~", paste(active_terms, collapse = " + ")))
  
  # Fit reduced model
  mdl_red <- bam(red_form, data = dat_clean, method = "fREML", discrete = TRUE, nthreads = 4)
  
  # Calculate Dev Lost
  effect_sizes$Partial_Deviance[i] <- full_dev - (summary(mdl_red)$dev.expl * 100)
}

# 3B. Partial Eta-Squared from F-statistics
sum_tab <- as.data.frame(summary(mdl_final)$s.table)
sum_tab$Term <- rownames(sum_tab)
res_df <- mdl_final$df.residual

for(i in 1:length(smooth_labels)) {
  term_row <- sum_tab[sum_tab$Term == smooth_labels[i], ]
  if(nrow(term_row) == 1) {
    # Eta_p^2 = (F * edf) / (F * edf + df_residual)
    eta_sq <- (term_row$F * term_row$edf) / ((term_row$F * term_row$edf) + res_df)
    effect_sizes$Eta_Sq_Partial[i] <- eta_sq
  }
}

# 3C & 3D. Amplitude (Peak-to-Trough) and RMS (Average Impact)
for(i in 1:length(smooth_labels)) {
  sm <- smooth_estimates(mdl_final, smooth = smooth_labels[i], n = 100)
  effect_sizes$Amplitude[i] <- max(sm$.estimate) - min(sm$.estimate)
  effect_sizes$RMS[i] <- sqrt(mean(sm$.estimate^2))
}

# Print Results cleanly to console
message("\n--- COMPREHENSIVE EFFECT SIZES (RAW UNITS) ---")
print(effect_sizes %>% arrange(desc(Partial_Deviance)))



# ==============================================================================
# --- 4. EFFECT SIZE DASHBOARD VISUALIZATION ---
# ==============================================================================
# Reshape data for faceted plotting
eff_long <- effect_sizes %>%
  pivot_longer(cols = -Term, names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = factor(Metric, 
                         levels = c("Partial_Deviance", "Eta_Sq_Partial", "Amplitude", "RMS"),
                         labels = c("Deviance Explained (%)", "Partial Eta-Squared", 
                                    "Peak-to-Trough Amplitude", "RMS Effect")))

p_dashboard <- ggplot(eff_long, aes(x = reorder(Term, Value), y = Value, fill = Term)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  coord_flip() +
  facet_wrap(~Metric, scales = "free_x", ncol = 2) +
  scale_fill_viridis_d(option = "mako") +
  theme_bw() +
  labs(title = "Effect Size Dashboard",
       subtitle = "Comparing alternative metrics for biological importance",
       x = NULL, y = NULL) +
  theme(strip.text = element_text(face = "bold", size = 11))

print(p_dashboard)

cairo_pdf("Effect_Size_Post.pdf", width = 8, height = 6)
print(p_dashboard)
dev.off()

#dev.new(noRStudioGD = TRUE); print(p_dashboard)


# ==============================================================================
# --- 5. SAVE MODEL OUTPUT ---
# ==============================================================================
library(dplyr)

# 1. Clean up the summary table for merging
model_stats <- as.data.frame(summary(mdl_final)$s.table) %>%
  mutate(Term = rownames(.))

# 2. Join with your effect_sizes dataframe
# We use left_join so random effects (AnimalID/SessionID) stay in the table 
# even if you didn't calculate Amplitude/RMS for them.
combined_results <- model_stats %>%
  left_join(effect_sizes, by = "Term") %>%
  select(Term, everything()) # Put Term name in the first column

# 3. Save to a single CSV
write.csv(combined_results, "GAM_model_post.csv", row.names = FALSE)

message("Files saved: 'GAM_model_post.csv'")

message("\nAnalysis Pipeline Complete!")






# ---------------------------------------------------------
# Plot 1: Main Effect of POST Match Phase
# ---------------------------------------------------------
# Change smooth name to POST
sm_match_post <- smooth_estimates(mdl_final, smooth = "s(SOPhasePOST_Match)", n = 200) %>%
  mutate(
    SOPhasePOST_Match = ifelse(SOPhasePOST_Match < 0, SOPhasePOST_Match + 2*pi, SOPhasePOST_Match),
    .lower_ci = .estimate - (1.96 * .se),
    .upper_ci = .estimate + (1.96 * .se)
  ) %>% arrange(SOPhasePOST_Match)

p_match_post <- ggplot(sm_match_post, aes(x = SOPhasePOST_Match, y = .estimate)) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), alpha = 0.2, fill = "dodgerblue") +
  geom_line(linewidth = 1, color = "dodgerblue") +
  theme_bw() + 
  scale_x_continuous(breaks = c(0, pi, 2*pi), labels = c("0", expression(pi), expression(2*pi)),
                     limits = c(0, 2*pi), expand = c(0, 0)) +
  labs(title = "Isolated Main Effect: POST Match SO Phase", 
       x = "Match SO Phase (Rad)", y = "Partial Effect")

dev.new(noRStudioGD = TRUE)
print(p_match_post)
# ---------------------------------------------------------
# Plot 2: Main Effect of POST Non-Match Phase
# ---------------------------------------------------------
sm_nonmatch_post <- smooth_estimates(mdl_final, smooth = "s(SOPhasePOST_NonMatch)", n = 200) %>%
  mutate(
    SOPhasePOST_NonMatch = ifelse(SOPhasePOST_NonMatch < 0, SOPhasePOST_NonMatch + 2*pi, SOPhasePOST_NonMatch),
    .lower_ci = .estimate - (1.96 * .se),
    .upper_ci = .estimate + (1.96 * .se)
  ) %>% arrange(SOPhasePOST_NonMatch)

p_nonmatch_post <- ggplot(sm_nonmatch_post, aes(x = SOPhasePOST_NonMatch, y = .estimate)) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), alpha = 0.2, fill = "firebrick") +
  geom_line(linewidth = 1, color = "firebrick") +
  theme_bw() + 
  scale_x_continuous(breaks = c(0, pi, 2*pi), labels = c("0", expression(pi), expression(2*pi)),
                     limits = c(0, 2*pi), expand = c(0, 0)) +
  labs(title = "Isolated Main Effect: POST Non-Match SO Phase", 
       x = "Non-Match SO Phase (Rad)", y = "Partial Effect")

dev.new(noRStudioGD = TRUE)
print(p_nonmatch_post)

# ---------------------------------------------------------
# Plot 3: The POST Phase Interaction (ti term)
# ---------------------------------------------------------

sm_2d <- smooth_estimates(mdl_final, smooth = "ti(SOPhasePOST_Match,SOPhasePOST_NonMatch)", n = 100) %>%
  mutate(
    .lower_ci = .estimate - (1.96 * .se),
    .upper_ci = .estimate + (1.96 * .se)
  )

p_synergy <- sm_2d %>%
  mutate(
    SOPhasePOST_Match = ifelse(SOPhasePOST_Match < 0, SOPhasePOST_Match + 2*pi, SOPhasePOST_Match),
    SOPhasePOST_NonMatch = ifelse(SOPhasePOST_NonMatch < 0, SOPhasePOST_NonMatch + 2*pi, SOPhasePOST_NonMatch)
  ) %>%
  arrange(SOPhasePOST_Match, SOPhasePOST_NonMatch) %>%
  ggplot(aes(x = SOPhasePOST_Match, y = SOPhasePOST_NonMatch, fill = .estimate)) +
  geom_tile(width = (2*pi)/99, height = (2*pi)/99) + 
  geom_contour(aes(z = .estimate), color = "black", alpha = 0.2) + 
  scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "firebrick", midpoint = 0, name = "Effect") +
  scale_x_continuous(breaks = c(0, pi, 2*pi), labels = c("0", expression(pi), expression(2*pi))) +
  scale_y_continuous(breaks = c(0, pi, 2*pi), labels = c("0", expression(pi), expression(2*pi))) +
  theme_minimal(base_family = "Arial") +
  theme(aspect.ratio = 1) +
  labs(title = "Inter-hemispheric Phase interaction effect on Pre coherence",
       x = "Match Post SO Phase (Rad)", y = "Non-Match Post SO Phase (Rad)")

dev.new(noRStudioGD = TRUE)
print(p_synergy)



# ---------------------------------------------------------
# Total Landscape Reconstruction (POST)
# ---------------------------------------------------------
pred_grid_post <- expand.grid(
  SOPhasePOST_Match = seq(0, 2*pi, length.out = 100),
  SOPhasePOST_NonMatch = seq(0, 2*pi, length.out = 100),
  RipplePower_Z = 0,
  SpindlePowerPOST_Match_Z = 0, # Assuming you have POST versions of these too
  SpindlePowerPOST_NonMatch_Z = 0, # Assuming you have POST versions of these too
  
  SessionID = dat_clean$SessionID[1],
  AnimalID = dat_clean$AnimalID[1] 
)

term_preds_post <- predict(mdl_final, newdata = pred_grid_post, type = "terms")

# Reconstructing using POST terms
pred_grid_post$Reconstructed_TE <- 
  term_preds_post[, "s(SOPhasePOST_Match)"] + 
  term_preds_post[, "s(SOPhasePOST_NonMatch)"] + 
  term_preds_post[, "ti(SOPhasePOST_Match,SOPhasePOST_NonMatch)"]

p_reconstructed_post <- ggplot(pred_grid_post, aes(x = SOPhasePOST_Match, y = SOPhasePOST_NonMatch, fill = Reconstructed_TE)) +
  geom_tile() + 
  geom_contour(aes(z = Reconstructed_TE), color = "black", alpha = 0.2) + 
  scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "firebrick", midpoint = 0) +
  scale_x_continuous(breaks = c(0, pi, 2*pi), labels = c("0", expression(pi), expression(2*pi))) +
  scale_y_continuous(breaks = c(0, pi, 2*pi), labels = c("0", expression(pi), expression(2*pi))) +
  theme_minimal() +
  theme(aspect.ratio = 1) +
  labs(title = "Total Phase effect on coherence (POST)", 
       x = "Match SO Phase POST", y = "Non-Match SO Phase POST")

dev.new(noRStudioGD = TRUE)
print(p_reconstructed_post)