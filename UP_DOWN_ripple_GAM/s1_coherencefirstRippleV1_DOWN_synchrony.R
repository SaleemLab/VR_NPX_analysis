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
    lastRippleNormalisedUP_z = as.numeric(scale(lastRippleNormalisedUP))
  )

# 1. Calculate Signed Geometric Mean (Updated names with HPC/V1 suffixes)
dat$geo_coherenceFirst <- sign(dat$firstRippleV1 * dat$firstRippleHPC) * 
  sqrt(abs(dat$firstRippleV1 * dat$firstRippleHPC))

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

dat$geo_coherenceFirstRippleNext <- sign(dat$firstRippleHPC * dat$nextUPV1) * 
  sqrt(abs(dat$firstRippleHPC * dat$nextUPV1))

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
my_folder <- "C:/Users/masah/Documents/GitHub/VR_NPX_analysis/UP_DOWN_ripple_GAM/coherenceLastRippleV1"

# 2. Check if it exists; if not, create it
if (!dir.exists(my_folder)) {
  dir.create(my_folder, recursive = TRUE)
}

# 3. Change the working directory to that folder
setwd(my_folder)



# Time from last ripple boundary.
# q_boundaries <- c(0, quantile(dat_clean$TimefromLastRipple, probs = c(0.25, 0.5,0.75), na.rm = TRUE), Inf)
q_boundaries <- c(0, quantile(dat_clean$TimefromLastRipple, probs = c(0.2, 0.4,0.6,0.8), na.rm = TRUE), Inf)
# q_boundaries <- c(0, quantile(dat_clean$TimefromLastRipple, probs = c(0.2,0.4,0.6,0.8), na.rm = TRUE), Inf)

# 2. Cut the continuous variable into 4 labeled factor levels
dat_clean$ripple_quantile <- cut(
  dat_clean$TimefromLastRipple,
  breaks = q_boundaries,
  labels = c("Q1", "Q2", "Q3", "Q4","Q5"),
  include.lowest = TRUE
)

# Get categorical ripple within and outside of 100ms or 50ms (50%)
# dat_clean$is_near_DOWN <- ifelse(dat_clean$TimefromLastRipple < 0.05, "yes", "no")
dat_clean$is_near_DOWN <- ifelse(dat_clean$TimefromLastRipple < 0.1, "yes", "no")
dat_clean$is_near_DOWN <- as.factor(dat_clean$is_near_DOWN)

# Log transform with tiny offset to avoid inf
dat_clean$log_TimefromLastRipple <- log(dat_clean$TimefromLastRipple + 0.000001)
dat_clean <- dat_clean %>%
  mutate(log_TimefromLastRipple_z = as.numeric(scale(log_TimefromLastRipple)))



#########
######### 
# z_thresh = 3.5
# z_cols <- c(
#   "firstRipplePower_z",
#   # "nextDOWNlag_z","nextDOWNSOPower_z",
#   "lastRipplePower_z"
#   # "TimefromfirstRipple_z"
# )


z_thresh = 2.56
z_cols <- c(
  # "geo_coherenceFirst_z",
  # "geo_coherence_z","geo_coherenceRippleNext_z",
  "nextDOWNlag_z","nextDOWNSOPower_z",
  "firstRipplePower_z",
  "TimefromFirstRipple_z",
  "TimefromLastRipple_z",
  "lastRipplePower_z"
)

dat_clean1 <- dat_Ripples %>%
  filter(
    if_all(all_of(z_cols), ~ abs(.) < z_thresh | is.na(.))
  )

# dat_clean1 <- dat_clean %>%
#   filter(
#     # Keep if the value is <= 1s OR if the value is NaN/NA
#     (RippleCounts>1) 
#     #(UPDuration_Match >= 0.1    | is.na(UPDuration_Match)) &
#     #(UPDuration_NonMatch >= 0.1 | is.na(UPDuration_NonMatch)) 
#   )


# dat_clean1 <- dat_clean %>%
#   filter(
#     if_all(all_of(z_cols), ~ abs(.) < z_thresh | is.na(.))
#   )
# geo_coherenceFirst_z


# z_cols <- c(
#   "nextDOWNlag_z","nextDOWNSOPower_z","geo_coherence_z","geo_coherenceRippleNext_z",
#   "TimefromLastRipple_z","lastRipplePower_z","geo_coherenceFirst_z"
# )

# dat_clean1 <- dat_clean %>%
#   filter(
#     if_all(all_of(z_cols), ~ abs(.) < z_thresh | is.na(.))
#   )



mdl_final <- bam(geo_coherenceFirst_z ~ 
                   # 1. Surviving Main Power Effects
                   # s(nextDOWNDuration_z, k = 5) +
                   # s(nextDOWNSOPower_z, k = 5) +
                   # s(nextDOWNlag_z, k = 5) +
                   # s(nextDOWNDuration_z, k = 5) +
                   
                   # is_near_DOWN +
                   # s(lastRipplePower_z,by = is_near_DOWN,k=5)+
                   # is_near_DOWN +
                   # s(firstRipplePower_z,by = is_near_DOWN,k=5)+
                   s(firstRipplePower_z,k=5)+
                   
                   # s(lastRipplePower_z,k=5)+
                   # ti(lastRipplePower_z, nextDOWNSOPower_z, k = 5)+
                   # te(nextDOWNSOPower_z, nextDOWNlag_z, k = 5)+
                   # s(TimeToFirstRipple_z, k = 5) +
                   # ti(firstRipplePower_z, TimeToFirstRipple_z, k = 5)+
                   # s(TimefromFirstRipple_z, k = 5) +
                   # ti(firstRipplePower_z, TimefromFirstRipple_z, k = 5)+
                   
                   # s(TimefromLastRipple_z, k = 5) +
                   # ti(lastRipplePower_z, TimefromLastRipple_z, k = 5)+
                   # ti(firstRipplePower_z, TimefromFirstRipple_z, k = 5)+
                   # s(nextDOWNSOPower_z, k = 5) +
                   # s(nextDOWNlag_z, k = 5) +
                   # ti(nextDOWNSOPower_z, nextDOWNlag_z, k = 5)+
                   
                   # s(TimetoNextUP_z, k = 5) +
                   
                   # 4. Control Term
                   s(AnimalID, bs = "re") +
                   s(SessionID, bs = "re"), 
                 
                 data = dat_clean1 , 
                 method = "fREML", 
                 discrete = TRUE, 
                 nthreads = 4)

#message("\n--- FINAL MODEL SUMMARY ---")
print(summary(mdl_final))


### first ripple power -> first ripple coherence
# Calculate scaling factors
raw_breaks <- c(5,10,15)
z_breaks <- sapply(raw_breaks, function(val) {
  dat_clean1$firstRipplePower_z[which.min(abs(dat_clean1$firstRipplePower - val))]
})

# cairo_pdf("lateUPHPC_nextUP", width = 4.3, height = 4.3)
p_lag_raw <- draw(mdl_final, select = "s(firstRipplePower_z)", residuals = FALSE, rug = FALSE) + 
  theme_bw(base_family = "Arial") + 
  theme(aspect.ratio = 1) +
  scale_x_continuous(breaks = z_breaks, labels = raw_breaks) +
  coord_cartesian(xlim = c(-1.4,2.6),ylim=c(-0.1,0.1),expand = FALSE) + 
  labs(
    title = "Ripple power on first ripple coherence", 
    x = "Ripple power", 
    y = "Partial Effect"
  )
dev.new(noRStudioGD = TRUE)
print(p_lag_raw)


# p_hist_raw <- ggplot(dat_clean1, aes(x = firstRipplePower)) +
#   geom_histogram(bins = 40, fill = "grey70", color = "black", linewidth = 0.2) +
#   theme_bw(base_family = "Arial") +
#   theme(aspect.ratio = 1) +
#   coord_cartesian(xlim = c(min(raw_breaks) - 5, max(raw_breaks) + 5)) +
#   scale_x_continuous(breaks = raw_breaks) +
#   labs(
#     title = "Distribution of ripple power",
#     x = "Ripple power",
#     y = "Count"
#   )
# 
# dev.new(noRStudioGD = TRUE)
# print(p_hist_raw)




mdl_final <- bam(geo_coherence_z ~ 
                   # 1. Surviving Main Power Effects
                   # s(nextDOWNDuration_z, k = 5) +
                   # s(nextDOWNSOPower_z, k = 5) +
                   # s(nextDOWNlag_z, k = 5) +
                   # s(nextDOWNDuration_z, k = 5) +
                   
                   # is_near_DOWN +
                   # s(lastRipplePower_z,by = is_near_DOWN,k=5)+
                   # is_near_DOWN +
                   # s(firstRipplePower_z,by = is_near_DOWN,k=5)+
                   # s(firstRipplePower_z,k=5)+
                   
                   s(lastRipplePower_z,k=5)+
                   # ti(lastRipplePower_z, nextDOWNSOPower_z, k = 5)+
                   # te(nextDOWNSOPower_z, nextDOWNlag_z, k = 5)+
                   # s(TimefromFirstRipple_z, k = 5) +
                   # s(TimefromLastRipple_z, k = 5) +
                   # ti(lastRipplePower_z, TimefromLastRipple_z, k = 5)+
                   # ti(firstRipplePower_z, TimefromFirstRipple_z, k = 5)+
                   # s(nextDOWNSOPower_z, k = 5) +
                   # s(nextDOWNlag_z, k = 5) +
                   # ti(nextDOWNSOPower_z, nextDOWNlag_z, k = 5)+
                   
                   # s(TimetoNextUP_z, k = 5) +
                   
                   # 4. Control Term
                   s(AnimalID, bs = "re") +
                   s(SessionID, bs = "re"), 
                 
                 data = dat_clean1 , 
                 method = "fREML", 
                 discrete = TRUE, 
                 nthreads = 4)

#message("\n--- FINAL MODEL SUMMARY ---")
print(summary(mdl_final))


### last ripple power -> last ripple coherence
# Calculate scaling factors
raw_breaks <- c(5,10,15)
z_breaks <- sapply(raw_breaks, function(val) {
  dat_clean1$lastRipplePower_z[which.min(abs(dat_clean1$lastRipplePower - val))]
})

# cairo_pdf("lateUPHPC_nextUP", width = 4.3, height = 4.3)
p_lag_raw <- draw(mdl_final, select = "s(lastRipplePower_z)", residuals = FALSE, rug = FALSE) + 
  theme_bw(base_family = "Arial") + 
  theme(aspect.ratio = 1) +
  scale_x_continuous(breaks = z_breaks, labels = raw_breaks) +
  coord_cartesian(xlim = c(-1.3,2.3),ylim=c(-0.12,0.12),expand = FALSE) +
  labs(
    title = "Ripple power on last ripple coherence", 
    x = "Ripple power", 
    y = "Partial Effect"
  )
dev.new(noRStudioGD = TRUE)
print(p_lag_raw)


p_hist_raw <- ggplot(dat_clean1, aes(x = lastRipplePower)) +
  geom_histogram(bins = 40, fill = "grey70", color = "black", linewidth = 0.2) +
  theme_bw(base_family = "Arial") +
  theme(aspect.ratio = 1) +
  coord_cartesian(xlim = c(min(raw_breaks) - 5, max(raw_breaks) + 5)) +
  scale_x_continuous(breaks = raw_breaks) +
  labs(
    title = "Distribution of ripple power",
    x = "Ripple power",
    y = "Count"
  )

dev.new(noRStudioGD = TRUE)
print(p_hist_raw)




### HC bias
# Calculate scaling factors
raw_breaks <- c(0,0.1,0.2,0.3,0.4,0.5)
z_breaks <- sapply(raw_breaks, function(val) {
  dat_clean1$TimefromFirstRipple_z[which.min(abs(dat_clean1$TimefromFirstRipple - val))]
})

# cairo_pdf("lateUPHPC_nextUP", width = 4.3, height = 4.3)
p_lag_raw <- draw(mdl_final, select = "s(TimefromFirstRipple_z)", residuals = FALSE, rug = FALSE) + 
  theme_bw(base_family = "Arial") + 
  theme(aspect.ratio = 1) +
  scale_x_continuous(breaks = z_breaks, labels = raw_breaks) +
  coord_cartesian(xlim = c(-0.845,0.49)) + 
  labs(
    title = "TimefromFirstRipple on Ripple V1", 
    x = "TimefromLastRipple", 
    y = "Partial Effect"
  )
dev.new(noRStudioGD = TRUE)
print(p_lag_raw)






# --- TimefromLastRipple ---

# Define breaks based on data range (adjust binwidth as needed)
lastRipple_breaks <- pretty(dat_clean1$TimefromLastRipple)

p_hist_lastRipple <- ggplot(dat_clean1, aes(x = TimefromLastRipple)) +
  geom_histogram(bins = 40, fill = "grey70", color = "black", linewidth = 0.2) +
  theme_bw(base_family = "Arial") +
  theme(aspect.ratio = 1) +
  coord_cartesian(xlim = c(min(lastRipple_breaks), max(lastRipple_breaks))) +
  scale_x_continuous(breaks = lastRipple_breaks) +
  labs(
    title = "Distribution of time from last ripple",
    x = "Time from last ripple",
    y = "Count"
  )

dev.new(noRStudioGD = TRUE)
print(p_hist_lastRipple)


# --- TimefromFirstRipple ---

firstRipple_breaks <- pretty(dat_clean1$TimefromFirstRipple)

p_hist_firstRipple <- ggplot(dat_clean1, aes(x = TimefromFirstRipple)) +
  geom_histogram(bins = 40, fill = "grey70", color = "black", linewidth = 0.2) +
  theme_bw(base_family = "Arial") +
  theme(aspect.ratio = 1) +
  coord_cartesian(xlim = c(min(firstRipple_breaks), max(firstRipple_breaks))) +
  scale_x_continuous(breaks = firstRipple_breaks) +
  labs(
    title = "Distribution of time from first ripple",
    x = "Time from first ripple",
    y = "Count"
  )

dev.new(noRStudioGD = TRUE)
print(p_hist_firstRipple)









# 
# ### HC bias
# # Calculate scaling factors
# raw_breaks <- c(0,0.1,0.2,0.3,0.4,0.5)
# z_breaks <- sapply(raw_breaks, function(val) {
#   dat_clean1$TimefromLastRipple_z[which.min(abs(dat_clean1$TimefromLastRipple - val))]
# })
# 
# # cairo_pdf("lateUPHPC_nextUP", width = 4.3, height = 4.3)
# p_lag_raw <- draw(mdl_final, select = "s(TimefromLastRipple_z)", residuals = FALSE, rug = FALSE) + 
#   theme_bw(base_family = "Arial") + 
#   theme(aspect.ratio = 1) +
#   scale_x_continuous(breaks = z_breaks, labels = raw_breaks) +
#   coord_cartesian(xlim = c(-0.845,0.49)) + 
#   labs(
#     title = "TimefromLastRipple on nextUP V1", 
#     x = "TimefromLastRipple", 
#     y = "Partial Effect"
#   )
# dev.new(noRStudioGD = TRUE)
# print(p_lag_raw)
# 
# 


# ==============================================================================
# --- interaction between next DOWN SO power and next DOWN bilateral lag ---
# ==============================================================================
# 1. Define sequence ranges using your actual RAW data minimums and maximums
# (Adjust the min/max or length if you want rounded limits like seq(0, 100, by=1))
raw_range_so   <- seq(1, 
                      4, length.out = 100)

raw_range_lag  <- seq(0, 
                      0.2, length.out = 100)

# 2. Build the prediction grid using the RAW scales
pred_grid_ti <- expand.grid(
  lastRipplePower_z   = 0,
  TimefromLastRipple_z   = 0,
  nextDOWNSOPower   = raw_range_so,
  nextDOWNlag       = raw_range_lag,
  nextDOWNDuration_z = 0,   # Constant mean for other model variables
  AnimalID           = dat_clean1$AnimalID[1],
  SessionID          = dat_clean1$SessionID[1]
)

# 3. Create the Z-scored columns that the model expects 
# This dynamically matches how your data was scaled (Mean=0, SD=1)
pred_grid_ti$nextDOWNSOPower_z <- (pred_grid_ti$nextDOWNSOPower - mean(dat_clean1$nextDOWNSOPower, na.rm=TRUE)) / sd(dat_clean1$nextDOWNSOPower, na.rm=TRUE)
pred_grid_ti$nextDOWNlag_z     <- (pred_grid_ti$nextDOWNlag     - mean(dat_clean1$nextDOWNlag, na.rm=TRUE))     / sd(dat_clean1$nextDOWNlag, na.rm=TRUE)

# 4. Extract specific term components matrix using the z-scores
term_preds <- predict(mdl_final, newdata = pred_grid_ti, type = "terms")

# Isolate the pure interaction tensor term
pred_grid_ti$Interaction_Effect <- term_preds[, "ti(nextDOWNSOPower_z,nextDOWNlag_z)"]

# 5. Plot Pure Interaction Surface using RAW scales for X and Y axes
p_interaction_raw_scale <- ggplot(pred_grid_ti, aes(x = nextDOWNSOPower, y = nextDOWNlag, fill = Interaction_Effect)) +
  geom_tile() + 
  geom_contour(aes(z = Interaction_Effect), color = "black", alpha = 0.2) + 
  scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "firebrick", midpoint = 0, name = "Effect") +
  theme_minimal() +
  theme(aspect.ratio = 1) +
  labs(
    title = "Pure Tensor Interaction (ti Term)",
    subtitle = "Variance unique to the combination (Raw Scale Mapping)",
    x = "nextDOWNSOPower (Raw)", 
    y = "nextDOWNlag (Raw)"
  )

print(p_interaction_raw_scale)




# 
# # ==============================================================================
# # --- 2D TE of next DOWN SO power and next DOWN bilateral lag ---
# # ==============================================================================
# 
# # 1. Define sequence ranges using your actual RAW data minimums and maximums
# raw_range_so   <- seq(1, 4, length.out = 100)
# raw_range_lag  <- seq(0, 0.5, length.out = 100)
# 
# # 2. Build the prediction grid using the RAW scales
# pred_grid_te <- expand.grid(
#   nextDOWNSOPower    = raw_range_so,
#   nextDOWNlag        = raw_range_lag,
#   AnimalID           = dat_clean1$AnimalID[1],  # Constant random effect
#   SessionID          = dat_clean1$SessionID[1]   # Constant random effect
# )
# 
# # 3. Create the Z-scored columns that the model expects 
# pred_grid_te$nextDOWNSOPower_z <- (pred_grid_te$nextDOWNSOPower - mean(dat_clean1$nextDOWNSOPower, na.rm=TRUE)) / sd(dat_clean1$nextDOWNSOPower, na.rm=TRUE)
# pred_grid_te$nextDOWNlag_z     <- (pred_grid_te$nextDOWNlag     - mean(dat_clean1$nextDOWNlag, na.rm=TRUE))     / sd(dat_clean1$nextDOWNlag, na.rm=TRUE)
# 
# # 4. Extract specific term components matrix using the z-scores
# term_preds <- predict(mdl_final, newdata = pred_grid_te, type = "terms")
# 
# # --- CRITICAL FIX HERE ---
# # mgcv drops the spaces and the 'k=5' syntax in the term name matrix column!
# pred_grid_te$Interaction_Effect <- term_preds[, "te(nextDOWNSOPower_z,nextDOWNlag_z)"]
# 
# # 5. Plot Total Tensor Surface using RAW scales for X and Y axes
# p_interaction_raw_scale <- ggplot(pred_grid_te, aes(x = nextDOWNSOPower, y = nextDOWNlag, fill = Interaction_Effect)) +
#   geom_tile() + 
#   geom_contour(aes(z = Interaction_Effect), color = "black", alpha = 0.2) + 
#   scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "firebrick", midpoint = 0, name = "Effect") +
#   theme_minimal() +
#   theme(aspect.ratio = 1) +
#   labs(
#     title = "Full Tensor Product Surface (ti Term)",
#     subtitle = "Combined Main Effects + Interaction (Raw Scale Mapping)",
#     x = "nextDOWNSOPower (Raw)", 
#     y = "nextDOWNlag (Raw)"
#   )
# 
# print(p_interaction_raw_scale)




# 
# # ==============================================================================
# # --- 2D Ti of next DOWN SO power and last ripple  ---
# # ==============================================================================
# 
# # 1. Define sequence ranges using your actual RAW data minimums and maximums
# raw_range_so   <- seq(1, 4, length.out = 100)
# raw_range_ripplePower  <- seq(5, 18, length.out = 100)
# 
# # 2. Build the prediction grid using the RAW scales
# pred_grid_te <- expand.grid(
#   nextDOWNSOPower    = raw_range_so,
#   lastRipplePower        =raw_range_ripplePower,
#   TimefromLastRipple_z           = 0,  # other held at 0
#   nextDOWNlag_z           = 0,  # other held at 0
# 
#   AnimalID           = dat_clean1$AnimalID[1],  # Constant random effect
#   SessionID          = dat_clean1$SessionID[1]   # Constant random effect
# )
# 
# # 3. Create the Z-scored columns that the model expects 
# pred_grid_te$nextDOWNSOPower_z <- (pred_grid_te$nextDOWNSOPower - mean(dat_clean1$nextDOWNSOPower, na.rm=TRUE)) / sd(dat_clean1$nextDOWNSOPower, na.rm=TRUE)
# pred_grid_te$lastRipplePower_z     <- (pred_grid_te$lastRipplePower     - mean(dat_clean1$lastRipplePower, na.rm=TRUE))     / sd(dat_clean1$lastRipplePower, na.rm=TRUE)
# 
# # 4. Extract specific term components matrix using the z-scores
# term_preds <- predict(mdl_final, newdata = pred_grid_te, type = "terms")
# 
# # --- CRITICAL FIX HERE ---
# # mgcv drops the spaces and the 'k=5' syntax in the term name matrix column!
# pred_grid_te$Interaction_Effect <- term_preds[, "ti(lastRipplePower_z,nextDOWNSOPower_z)"]
# 
# # 5. Plot Total Tensor Surface using RAW scales for X and Y axes
# p_interaction_raw_scale <- ggplot(pred_grid_te, aes(x = nextDOWNSOPower, y = lastRipplePower, fill = Interaction_Effect)) +
#   geom_tile() + 
#   geom_contour(aes(z = Interaction_Effect), color = "black", alpha = 0.2) + 
#   scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "firebrick", midpoint = 0, name = "Effect") +
#   theme_minimal() +
#   theme(aspect.ratio = 1) +
#   labs(
#     title = "Full Tensor Product Surface (ti Term)",
#     subtitle = "Variance unique to the combination (lastRipplePower and nextDOWNSOPower)",
#     x = "nextDOWNSOPower (Raw)", 
#     y = "lastRipplePower (Raw)"
#   )
# 
# print(p_interaction_raw_scale)
# 



