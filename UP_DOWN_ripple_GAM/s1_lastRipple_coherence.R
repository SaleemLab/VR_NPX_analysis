rm(list = ls())

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
dat <- read.csv("C:/Users/masah/Documents/GitHub/VR_NPX_analysis/UP_DOWN_ripple_GAM/UP_DOWN_info_GAM.csv")

dat$SessionID <- as.factor(dat$SessionID)
dat$AnimalID <- as.factor(dat$AnimalID)

dat <- dat %>%
  mutate(
    lastRippleNormalisedUP = case_when(
      lastRippleNormalisedUP < 0 ~ 0,
      lastRippleNormalisedUP > 1 ~ 1,
      TRUE ~ lastRippleNormalisedUP
      
    ),
    lastRippleNormalisedUP_z = as.numeric((lastRippleNormalisedUP))
  )

dat <- dat %>%
  mutate(
    firstRippleNormalisedUP = case_when(
      firstRippleNormalisedUP < 0 ~ 0,
      firstRippleNormalisedUP > 1 ~ 1,
      TRUE ~ firstRippleNormalisedUP
      
    ),
    firstRippleNormalisedUP_z = as.numeric((lastRippleNormalisedUP))
  )

# 1. Calculate Signed Geometric Mean (Updated names with HPC/V1 suffixes)
dat$geo_coherence <- sign(dat$lastRippleV1 * dat$lastRippleHPC) * 
  sqrt(abs(dat$lastRippleV1 * dat$lastRippleHPC))

dat$geo_coherenceLate <- sign(dat$lateUPV1 * dat$lateUPHPC) * 
  sqrt(abs(dat$lateUPV1 * dat$lateUPHPC))

dat$geo_coherenceRippleLate <- sign(dat$lateUPV1 * dat$lastRippleHPC) * 
  sqrt(abs(dat$lateUPV1 * dat$lastRippleHPC))

# 2. Calculate Signed Geometric Mean for PRE (Updated to use lastRippleHPCPRE)
dat$geo_coherencePRE <- sign(dat$lastRippleV1PRE * dat$lastRippleHPCPRE) * 
  sqrt(abs(dat$lastRippleV1PRE * dat$lastRippleHPCPRE))

# 1. Calculate Signed Geometric Mean (Updated names with HPC/V1 suffixes)
dat$geo_coherence1 <- sign(dat$firstRippleV1 * dat$firstRippleHPC) * 
  sqrt(abs(dat$firstRippleV1 * dat$firstRippleHPC))

# 2. Calculate Signed Geometric Mean for PRE (Updated to use lastRippleHPCPRE)
dat$geo_coherencePRE1 <- sign(dat$firstRippleV1PRE * dat$firstRippleHPC) * 
  sqrt(abs(dat$lastRippleV1PRE * dat$firstRippleHPC))

dat$geo_coherenceRippleNext <- sign(dat$lastRippleHPC * dat$nextUPV1) * 
  sqrt(abs(dat$lastRippleHPC * dat$nextUPV1))

dat$geo_coherenceNext <- sign(dat$lateUPHPC * dat$nextUPV1) * 
  sqrt(abs(dat$lateUPHPC * dat$nextUPV1))

dat$geo_coherenceNextV1 <- sign(dat$lateUPV1 * dat$nextUPV1) * 
  sqrt(abs(dat$lateUPV1 * dat$nextUPV1))


dat_Ripple <- dat %>%
  filter(
    # Keep if the value is <= 1s OR if the value is NaN/NA
    (RippleCounts>0) 
    #(UPDuration_Match >= 0.1    | is.na(UPDuration_Match)) &
    #(UPDuration_NonMatch >= 0.1 | is.na(UPDuration_NonMatch)) 
  )


dat_NoRipple <- dat %>%
  filter(
    # Keep if the value is <= 1s OR if the value is NaN/NA
    (is.na(RippleCounts)) 
    #(UPDuration_Match >= 0.1    | is.na(UPDuration_Match)) &
    #(UPDuration_NonMatch >= 0.1 | is.na(UPDuration_NonMatch)) 
  )

dat_lateUP_Ripple <- dat %>%
  filter(
    # Keep if the value is <= 1s OR if the value is NaN/NA
    #(UPDuration > 0.2    | is.na(UPDuration)) &
    (TimefromLastRipple < 0.1    | is.na(TimefromLastRipple)) 
    #(UPDuration_Match >= 0.1    | is.na(UPDuration_Match)) &
    #(UPDuration_NonMatch >= 0.1 | is.na(UPDuration_NonMatch)) 
  )

dat_lateUP_NoRipple <- dat %>%
  filter(
    # Keep if the value is <= 1s OR if the value is NaN/NA
    #(UPDuration > 0.2    | is.na(UPDuration)) &
    (TimefromLastRipple > 0.1    | is.na(TimefromLastRipple)) 
    #(UPDuration_Match >= 0.1    | is.na(UPDuration_Match)) &
    #(UPDuration_NonMatch >= 0.1 | is.na(UPDuration_NonMatch)) 
  )


dat_Time <- dat %>%
  filter(
    # Keep if the value is <= 1s OR if the value is NaN/NA
    (UPDuration > 0.2    | is.na(UPDuration)) 
    #(UPDuration_Match >= 0.1    | is.na(UPDuration_Match)) &
    #(UPDuration_NonMatch >= 0.1 | is.na(UPDuration_NonMatch)) 
  )


dat_Ripples <- dat %>%
  filter(
    # Keep if the value is <= 1s OR if the value is NaN/NA
    (RippleCounts>1) 
    #(UPDuration_Match >= 0.1    | is.na(UPDuration_Match)) &
    #(UPDuration_NonMatch >= 0.1 | is.na(UPDuration_NonMatch)) 
  )



# 3. Create Z-scores for all numeric columns
dat_clean <- dat %>%
  mutate(across(where(is.numeric), 
                ~ as.numeric(scale(.)), 
                .names = "{.col}_z"))

dat_Time <- dat_Time %>%
  mutate(across(where(is.numeric), 
                ~ as.numeric(scale(.)), 
                .names = "{.col}_z"))

dat_lateUP_Ripple <- dat_lateUP_Ripple %>%
  mutate(across(where(is.numeric), 
                ~ as.numeric(scale(.)), 
                .names = "{.col}_z"))
#
dat_lateUP_NoRipple <- dat_lateUP_NoRipple %>%
  mutate(across(where(is.numeric), 
                ~ as.numeric(scale(.)), 
                .names = "{.col}_z"))

dat_NoRipple <- dat_NoRipple %>%
  mutate(across(where(is.numeric), 
                ~ as.numeric(scale(.)), 
                .names = "{.col}_z"))

dat_Ripple <- dat_Ripple %>%
  mutate(across(where(is.numeric), 
                ~ as.numeric(scale(.)), 
                .names = "{.col}_z"))

dat_Ripples <- dat_Ripples %>%
  mutate(across(where(is.numeric), 
                ~ as.numeric(scale(.)), 
                .names = "{.col}_z"))

# --- 3. Clean Data ---
z_thresh <- 3.5  # include <99.9th centiles of data


########
######## Last ripple
########
z_cols <- c(
  "lastRippleV1_z", "lastRippleV1PRE_z",
  "lastRippleHPC_z", "lastRippleHPC_z",
  "TimefromLastRipple_z","lastRippleNormalisedUP_z"

)

dat_clean1 <- dat_clean %>%
  filter(
    if_all(all_of(z_cols), ~ abs(.) < z_thresh | is.na(.))
  )


########
######## Divide ripple events based on its recency to DOWN
########


# Time from last ripple boundary.
# q_boundaries <- c(0, quantile(dat_clean$TimefromLastRipple, probs = c(0.25, 0.5,0.75), na.rm = TRUE), Inf)
q_boundaries <- c(0, quantile(dat_clean$TimefromLastRipple, probs = c(0.2, 0.4,0.6,0.8), na.rm = TRUE), Inf)
# q_boundaries <- c(0, quantile(dat_clean$TimefromLastRipple, probs = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9), na.rm = TRUE), Inf)

# 2. Cut the continuous variable into 4 labeled factor levels
dat_clean$ripple_quantile <- cut(
  dat_clean$TimefromLastRipple,
  breaks = q_boundaries,
  # labels = c("Q1", "Q2", "Q3", "Q4"),
  labels = c("Q1", "Q2", "Q3", "Q4","Q5"),
  # labels = c("Q1", "Q2", "Q3", "Q4","Q5","Q6","Q7","Q8","Q9","Q10"),
  include.lowest = TRUE
)

# Get categorical ripple within and outside of 100ms
dat_clean$is_near_DOWN <- ifelse(dat_clean$TimefromLastRipple < 0.1, "yes", "no")
# dat_clean$is_near_DOWN <- ifelse(dat_clean$TimefromLastRipple < 0.1, "yes", "no")
dat_clean$is_near_DOWN <- as.factor(dat_clean$is_near_DOWN)

# Log transform with tiny offset to avoid inf
dat_clean$log_TimefromLastRipple <- log(dat_clean$TimefromLastRipple + 0.000001)
dat_clean <- dat_clean %>%
 mutate(log_TimefromLastRipple_z = as.numeric(scale(log_TimefromLastRipple)))



########
######## last ripple HC content (especially those close to DOWN state) predicts does not predict V1 ripple content
########

# List of the new Z-score columns to filter (Updated to include all relevant ripple predictors)
z_cols <- c(
  # "lastRippleV1_z", "lastRippleV1PRE_z",
  # "lastRippleHPC_z", "lastRippleHPCPRE_z",
  "geo_coherence_z","geo_coherenceNext_z",
  "log_TimefromLastRipple_z"
  # "lastRippleNormalisedUP_z"
  # "firstRippleV1_z", "firstRippleV1PRE_z",
  # "firstRippleHPC_z", "firstRippleHPCPRE_z"
)

# dat_clean1 <- dat_Ripples %>%
# dat_clean1 <- dat_lateUP_Ripple %>%

dat_clean1 <- dat_clean %>%
  
  filter(
    if_all(all_of(z_cols), ~ abs(.) < z_thresh | is.na(.))
  )

mdl_lastRipple <- bam(geo_coherence_z ~ 
                        s(log_TimefromLastRipple_z,k = 5) +
                        s(lastRippleNormalisedUP_z, bs = "tp",k = 5) +
                        # ti(lastRippleNormalisedUP_z,log_TimefromLastRipple_z,k=5)+

                        
                        # 4. Control Term
                        s(AnimalID, bs = "re") +
                        s(SessionID, bs = "re"), 
                      
                      data = dat_clean1, 
                      method = "fREML", 
                      discrete = TRUE, 
                      nthreads = 4)

message("\n--- FINAL MODEL SUMMARY ---")
print(summary(mdl_lastRipple))



### Time from last ripple to DOWN transition
# Calculate scaling factors
raw_breaks <- c(-6,-4,-2,0,2)
z_breaks <- sapply(raw_breaks, function(val) {
  dat_clean$log_TimefromLastRipple_z[which.min(abs(dat_clean$log_TimefromLastRipple - val))]
})

# cairo_pdf("lateUPHPC_nextUP", width = 4.3, height = 4.3)
p_TimefromLastRipple_raw <- draw(mdl_lastRipple, select = "s(log_TimefromLastRipple_z)", residuals = FALSE, rug = FALSE) + 
  theme_bw(base_family = "Arial") + 
  theme(aspect.ratio = 1) +
  scale_x_continuous(breaks = z_breaks, labels = raw_breaks) +
  labs(
    title = "Time from last ripple to DOWN and last ripple coherence", 
    x = "Time from last ripple to DOWN (log sec)", 
    y = "Partial Effect"
  )
dev.new(noRStudioGD = TRUE)
print(p_TimefromLastRipple_raw)

# 
# # ---------------------------------------------------------
# # Plot 3: Pure Interaction Tensor (ti) 
# # ---------------------------------------------------------
sm_2d <- smooth_estimates(mdl_pre_post, smooth = "ti(lastRippleHPC_z,TimefromLastRipple_z)", n = 100)

p_interaction <- ggplot(sm_2d, aes(x = V1_pre_z, y = HC_post_z, fill = .estimate)) +
  geom_tile() +
  geom_contour(aes(z = .estimate), color = "black", alpha = 0.2) +
  scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "firebrick", midpoint = 0, name = "Effect") +
  theme_minimal() +
  theme(aspect.ratio = 1) +
  labs(title = "Pure Tensor Interaction (ti Term)",
       subtitle = "Variance unique to the combination",
       x = "V1_pre_z", y = "HC_post_z")

print(p_interaction)

# 
# # ==============================================================================
# # --- RECREATING THE TOTAL LANDSCAPE ---
# # ==============================================================================
# 
# 1. Create a clean grid across the observed range of your data
# grid_range_v1 <- seq(min(dat_clean$V1_pre_z, na.rm=TRUE), max(dat_clean$V1_pre_z, na.rm=TRUE), length.out = 100)
# grid_range_hc <- seq(min(dat_clean$HC_post_z, na.rm=TRUE), max(dat_clean$HC_post_z, na.rm=TRUE), length.out = 100)
grid_range_v1 <- seq(-2.5,2.5, length.out = 40)
grid_range_hc <- seq(-2.5,2.5, length.out = 40)

pred_grid_total <- expand.grid(
  V1_pre_z  = grid_range_v1,
  HC_post_z = grid_range_hc,
  AnimalID  = dat_clean$AnimalID[1], # Held constant (ignored by type="terms" if excluded, but safe)
  SessionID  = dat_clean$SessionID[1]
)

# 2. Extract specific term components matrix
term_preds <- predict(mdl_pre_post, newdata = pred_grid_total, type = "terms")

# 3. Sum only your targets of interest
pred_grid_total$Reconstructed_Effect <-
  term_preds[, "s(HC_post_z)"] +
  term_preds[, "s(V1_pre_z)"]
# term_preds[, "ti(V1_pre_z,HC_post_z)"]

# 4. Plot Combined Surface
p_combined <- ggplot(pred_grid_total, aes(x = V1_pre_z, y = HC_post_z, fill = Reconstructed_Effect)) +
  geom_tile() +
  geom_contour(aes(z = Reconstructed_Effect), color = "black", alpha = 0.2) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "firebrick", midpoint = 0, name = "Total\nEffect") +
  theme_minimal() +
  # theme(aspect.ratio = 1) +
  labs(
    title = "Total Combined Effect Surface",
    subtitle = "Reconstructed: s(HC) + s(V1) + ti(V1, HC)",
    x = "V1_pre_z",
    y = "HC_post_z"
  )
dev.new(noRStudioGD = TRUE)
print(p_combined)