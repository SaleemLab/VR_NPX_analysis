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
#dat <- read.csv("C:/Users/masah/Documents/GitHub/VR_NPX_analysis/UP_DOWN_ripple_GAM/UP_DOWN_info_GAM.csv")

dat <- read.csv("C:/Users/masah/Documents/GitHub/VR_NPX_analysis/UP_DOWN_ripple_GAM/UP_DOWN_info_GAM_peak.csv")

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

dat$geo_coherenceLateNextV1 <- sign(dat$lateUPV1 * dat$nextUPV1) * 
  sqrt(abs(dat$lateUPV1 * dat$nextUPV1))

dat$geo_coherenceRippleNextV1 <- sign(dat$lastRippleV1 * dat$nextUPV1) * 
  sqrt(abs(dat$lastRippleV1 * dat$nextUPV1))

dat$geo_coherenceRippleNextV1PRE <- sign(dat$lastRippleV1PRE * dat$nextUPV1) * 
  sqrt(abs(dat$lastRippleV1 * dat$nextUPV1))

dat$geo_coherencePRE <- sign(dat$lastRippleV1PRE * dat$lastRippleHPC) * 
  sqrt(abs(dat$lastRippleV1PRE * dat$lastRippleHPC))

dat$geo_coherenceEarlyUPLastRipple <- sign(dat$earlyUPV1 * dat$lastRippleHPC) * 
  sqrt(abs(dat$earlyUPV1 * dat$lastRippleHPC))

dat$geo_coherenceEarlyNextV1 <- sign(dat$earlyUPV1 * dat$nextUPV1) * 
  sqrt(abs(dat$earlyUPV1 * dat$nextUPV1))


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
my_folder <- "C:/Users/masah/Documents/GitHub/VR_NPX_analysis/UP_DOWN_ripple_GAM/coherenceLateUP_NextUP_V1"

# 2. Check if it exists; if not, create it
if (!dir.exists(my_folder)) {
  dir.create(my_folder, recursive = TRUE)
}

# 3. Change the working directory to that folder
setwd(my_folder)



########
######## Divide ripple events based on its recency to DOWN
########


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

# Get categorical ripple within and outside of 100ms or 178ms (50%)
dat_clean$is_near_DOWN <- ifelse(dat_clean$TimefromLastRipple < 0.1, "yes", "no")
# dat_clean$is_near_DOWN <- ifelse(dat_clean$TimefromLastRipple < 0.1, "yes", "no")
dat_clean$is_near_DOWN <- as.factor(dat_clean$is_near_DOWN)

# Log transform with tiny offset to avoid inf
dat_clean$log_TimefromLastRipple <- log10(dat_clean$TimefromLastRipple + 0.000001)
dat_clean <- dat_clean %>%
  mutate(log_TimefromLastRipple_z = as.numeric(scale(log_TimefromLastRipple)))

dat_clean$log_nextDOWNDuration <- log10(dat_clean$nextDOWNDuration + 0.000001)
dat_clean <- dat_clean %>%
  mutate(log_nextDOWNDuration_z = as.numeric(scale(log_nextDOWNDuration)))


#########
#########


z_thresh = 3

z_cols <- c(
  # "geo_coherence_z","geo_coherenceRippleNext_z","geo_coherenceNext_z",\
  # "geo_coherenceLate_z",
  "geo_coherencePRE_z",
  "geo_coherenceRippleNextV1PRE_z",
  # "geo_coherenceLateNextV1_z",
  "nextDOWNlag_z","nextDOWNSOPower_z",
  "log_TimefromLastRipple_z",
  "log_nextDOWNDuration_z",
  "lastRippleMUArate_z"
  # "lastRipplePower_z"
)

# z_cols <- c(
#   # "geo_coherence_z","geo_coherenceRippleNext_z","geo_coherenceNext_z",\
#   "geo_coherenceEarlyUPLastRipple_z",
#   "geo_coherenceEarlyNextV1_z",
#   # "geo_coherencePRE_z",
#   # "geo_coherenceRippleNextV1_z",
#   # "geo_coherenceLateNextV1_z",
#   "nextDOWNlag_z","nextDOWNSOPower_z",
#   "log_TimefromLastRipple_z",
#   "log_nextDOWNDuration_z",
#   "lastRippleMUArate_z"
#   # "lastRipplePower_z"
# )


dat_clean1 <- dat_clean %>%
  filter(
    if_all(all_of(z_cols), ~ abs(.) < z_thresh | is.na(.)),
    # (RippleCounts>0) 
  )


# # Simple histogram with density curve overlay
# hist(dat_clean1$log_TimefromLastRipple, 
#      main = "Histogram of log_TimefromLastRipple_z",
#      xlab = "log_TimefromLastRipple_z",
#      col = "lightblue",
#      border = "black",
#      prob = TRUE)

# geo_coherenceRippleNextV1PRE_z
# geo_coherenceLateNextV1_z
mdl_final <- bam(geo_coherenceRippleNextV1PRE_z ~ 
                   # is_near_DOWN +
                   # s(geo_coherenceLate_z,by = is_near_DOWN,k=5)+
                   # s(lastRipplePower_z,k=5)+
                   s(lastRippleMUArate_z,k=5)+
                   # s(log_TimefromLastRipple_z, k = 5) +
                   # s(log_nextDOWNDuration_z, k = 5) +
                   # s(geo_coherenceLate_z, k = 5) +
                   s(geo_coherencePRE_z, k = 5) +
                   ti(lastRippleMUArate_z, geo_coherencePRE_z, k = 5)+
                   # ti(log_TimefromLastRipple_z, geo_coherenceLate_z, k = 5)+
                   
                   # ti(lastRipplePower_z, log_nextDOWNDuration_z, k = 5)+
                   # ti(log_TimefromLastRipple_z, log_nextDOWNDuration_z, k = 5)+
                   # ti(log_TimefromLastRipple_z, lastRipplePower_z, k = 5)+
                   
                   # s(nextDOWNSOPower_z, k = 5) +
                   # s(nextDOWNlag_z, k = 5) +
                   # ti(nextDOWNSOPower_z, nextDOWNlag_z, k = 5)+
                   # ti(lastRipplePower_z, nextDOWNSOPower_z, k = 5)+
                   # ti(lastRipplePower_z, nextDOWNlag_z, k = 5)+
                   
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

library(ggplot2)
# 
# dat_clean1 <- dat_clean %>%
#   filter(
#     if_all(all_of(z_cols), ~ abs(.) < z_thresh | is.na(.)),
#     (RippleCounts>0)
#   )

ggplot(dat_clean1, aes(x = geo_coherenceLate, y = geo_coherenceLateNextV1)) +
  geom_hex(bins = 100) +
  scale_fill_viridis_c() +
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  theme_minimal() +
  labs(title = "Density of geo_coherenceLate vs geo_coherenceLateNextV1",
       x = "geo_coherenceLate", y = "geo_coherenceLateNextV1",
       fill = "Count")

### geo_coherenceLate_z -> geo_coherenceLateNextV1_z

# Calculate scaling factors
raw_breaks <- c(-2,-1,0,1,2)
z_breaks <- sapply(raw_breaks, function(val) {
  dat_clean1$geo_coherenceEarlyUPLastRipple_z[which.min(abs(dat_clean1$geo_coherenceEarlyUPLastRipple - val))]
})

# cairo_pdf("lateUPHPC_nextUP", width = 4.3, height = 4.3)
p_lag_raw <- draw(mdl_final, select = "s(geo_coherenceEarlyUPLastRipple_z)", residuals = FALSE, rug = FALSE) +
  # p_lag_raw <- draw(mdl_final, select = "s(lastRipplePower_z):is_near_DOWNno", residuals = FALSE, rug = FALSE) +
  theme_bw(base_family = "Arial") +
  theme(aspect.ratio = 1) +
  scale_x_continuous(breaks = z_breaks, labels = raw_breaks) +
  # coord_cartesian(xlim = c(-1.3,1.8),ylim=c(-0.15,0.12),expand = FALSE) +
  # coord_cartesian(xlim = c(-1.3,2),ylim=c(-0.1,0.06),expand = FALSE) +
  labs(
    title = "EarlyUP and last ripple HC coherence effect on earlyUP - NextUP V1 coherence",
    x = "Early UP HC-V1 coherence",
    y = "Partial Effect"
  )
# dev.new(noRStudioGD = TRUE)
print(p_lag_raw)


# 
# ### last ripple Mua -> NextUP V1
# 
# # Calculate scaling factors
# raw_breaks <- c(0.25,0.5,0.75,1)
# z_breaks <- sapply(raw_breaks, function(val) {
#   dat_clean1$lastRippleMUArate_z[which.min(abs(dat_clean1$lastRippleMUArate - val))]
# })
# 
# # cairo_pdf("lateUPHPC_nextUP", width = 4.3, height = 4.3)
# p_lag_raw <- draw(mdl_final, select = "s(lastRippleMUArate_z)", residuals = FALSE, rug = FALSE) +
#   # p_lag_raw <- draw(mdl_final, select = "s(lastRipplePower_z):is_near_DOWNno", residuals = FALSE, rug = FALSE) +
#   theme_bw(base_family = "Arial") +
#   theme(aspect.ratio = 1) +
#   scale_x_continuous(breaks = z_breaks, labels = raw_breaks) +
#   # coord_cartesian(xlim = c(-1.3,1.8),ylim=c(-0.15,0.12),expand = FALSE) +
#   # coord_cartesian(xlim = c(-1.3,2),ylim=c(-0.1,0.06),expand = FALSE) +
#   labs(
#     title = "Ripple HC MUA rate on lateUP to nextUP V1 coherence",
#     x = "Ripple HC MUA",
#     y = "Partial Effect"
#   )
# # dev.new(noRStudioGD = TRUE)
# print(p_lag_raw)



### last ripple power -> NextUP V1

# Calculate scaling factors
raw_breaks <- c(5,7,9,11,13,15)
z_breaks <- sapply(raw_breaks, function(val) {
  dat_clean1$lastRipplePower_z[which.min(abs(dat_clean1$lastRipplePower - val))]
})

# cairo_pdf("lateUPHPC_nextUP", width = 4.3, height = 4.3)
p_lag_raw <- draw(mdl_final, select = "s(lastRipplePower_z)", residuals = FALSE, rug = FALSE) +
  # p_lag_raw <- draw(mdl_final, select = "s(lastRipplePower_z):is_near_DOWNno", residuals = FALSE, rug = FALSE) +
  theme_bw(base_family = "Arial") +
  theme(aspect.ratio = 1) +
  scale_x_continuous(breaks = z_breaks, labels = raw_breaks) +
  # coord_cartesian(xlim = c(-1.3,1.8),ylim=c(-0.15,0.12),expand = FALSE) +
  # coord_cartesian(xlim = c(-1.3,2),ylim=c(-0.1,0.06),expand = FALSE) +
  labs(
    title = "Ripple power on last ripple coherence",
    x = "Ripple power",
    y = "Partial Effect"
  )
# dev.new(noRStudioGD = TRUE)
print(p_lag_raw)


### Log Time from last ripple effect
# Calculate scaling factors
raw_breaks <- c(-3,-2,-1,0)
raw_labels <- parse(text = paste0("10^", raw_breaks))
z_breaks <- sapply(raw_breaks, function(val) {
  dat_clean1$log_TimefromLastRipple_z[which.min(abs(dat_clean1$log_TimefromLastRipple - val))]
})

# 1. Define desired raw time values in seconds (powers of 10)
raw_seconds <- c(10^-3, 10^-2, 10^-1, 10^0)

# 2. Convert raw seconds to log scale
# Use log10() if log_TimefromLastRipple was log10; use log() if it was natural log (ln)
log_vals <- log10(raw_seconds)  # yields: -3, -2, -1, 0

# 3. Calculate exact theoretical z-scores based on your dataset's mean and SD
mean_log <- mean(dat_clean1$log_TimefromLastRipple, na.rm = TRUE)
sd_log   <- sd(dat_clean1$log_TimefromLastRipple, na.rm = TRUE)

z_breaks <- (log_vals - mean_log) / sd_log

# 4. Formatter for 10^-3, 10^-2 notation
raw_labels <- parse(text = paste0("10^", log10(raw_seconds)))

# cairo_pdf("lateUPHPC_nextUP", width = 4.3, height = 4.3)
p_lag_raw <- draw(mdl_final, select = "s(log_TimefromLastRipple_z)", residuals = FALSE, rug = FALSE) + 
  theme_bw(base_family = "Arial") + 
  theme(aspect.ratio = 1) +
  scale_x_continuous(breaks = z_breaks, labels = raw_labels) +
  coord_cartesian(xlim = c(-2.8, 1.41),ylim = c(-0.06, 0.11),expand = FALSE) +
  
  labs(
    title = "TimefromLastRipple on LateUP-NextUP V1 coherence", 
    x = "TimefromLastRipple", 
    y = "Partial Effect"
  )
print(p_lag_raw)

cairo_pdf("LogTimefromLastRipple_lastRippleNextUPCoherence.pdf", width = 4.3, height = 4.3)
# dev.new(noRStudioGD = TRUE)
print(p_lag_raw)
dev.off()




### Log Time DOWN Duration
# Calculate scaling factors
raw_breaks <- c(0.03,0.06,0.12,0.25,0.5,1)
raw_labels <- raw_breaks
# raw_labels <- parse(text = paste0("10^", raw_breaks))
# z_breaks <- sapply(raw_breaks, function(val) {
#   dat_clean1$log_nextDOWNDuration_z[which.min(abs(dat_clean1$log_nextDOWNDuration - val))]
# })

# # 1. Define desired raw time values in seconds (powers of 10)
# raw_seconds <- c(10^-3, 10^-2, 10^-1, 10^0)
raw_seconds <- c(0.03,0.06,0.12,0.25,0.5,1)

# 2. Convert raw seconds to log scale
# Use log10() if log_nextDOWNDuration was log10; use log() if it was natural log (ln)
log_vals <- log10(raw_seconds)  # yields: -3, -2, -1, 0

# 3. Calculate exact theoretical z-scores based on your dataset's mean and SD
mean_log <- mean(dat_clean1$log_nextDOWNDuration, na.rm = TRUE)
sd_log   <- sd(dat_clean1$log_nextDOWNDuration, na.rm = TRUE)

z_breaks <- (log_vals - mean_log) / sd_log

# 4. Formatter for 10^-3, 10^-2 notation
# raw_labels <- parse(text = paste0("10^", log10(raw_seconds)))

# cairo_pdf("lateUPHPC_nextUP", width = 4.3, height = 4.3)
p_lag_raw <- draw(mdl_final, select = "s(log_nextDOWNDuration_z)", residuals = FALSE, rug = FALSE) + 
  theme_bw(base_family = "Arial") + 
  theme(aspect.ratio = 1) +
  scale_x_continuous(breaks = z_breaks, labels = raw_labels) +
  coord_cartesian(xlim = c(-2.6, 1.7),ylim = c(-0.08, 0.11),expand = FALSE) +
  
  labs(
    title = "log_nextDOWNDuration on lateUP nextUP V1 coherence", 
    x = "Next DOWN duration", 
    y = "Partial Effect"
  )
print(p_lag_raw)


cairo_pdf("log_nextDOWNDuration_LateUP_NextUP_V1_coherence.pdf", width = 4.3, height = 4.3)
# dev.new(noRStudioGD = TRUE)
print(p_lag_raw)
dev.off()



# ==============================================================================
# --- ti interaction: log_TimefromLastRipple_z and log_nextDOWNDuration_z ---
# ==============================================================================

# 1. Define sequence ranges using your actual RAW data minimums and maximums
# (Adjust the min/max values to match your dataset's raw range)
# raw_range_ripple_time <- seq(
#   min(dat_clean1$TimefromLastRipple, na.rm = TRUE),
#   max(dat_clean1$TimefromLastRipple, na.rm = TRUE),
#   length.out = 100
# )
# 
# raw_range_down_dur <- seq(
#   min(dat_clean1$nextDOWNDuration, na.rm = TRUE),
#   max(dat_clean1$nextDOWNDuration, na.rm = TRUE),
#   length.out = 100
# )

raw_range_ripple_time   <- seq(0.02,
                      0.25, length.out = 100)

raw_range_down_dur  <- seq(0.02,
                      0.25, length.out = 100)

# 2. Build the prediction grid using the RAW scales
pred_grid_ti <- expand.grid(
  TimefromLastRipple   = raw_range_ripple_time,
  nextDOWNDuration     = raw_range_down_dur,
  nextDOWNSOPower_z    = 0, # Constant mean for other model variables
  nextDOWNlag_z        = 0, # Constant mean for other model variables
  lastRipplePower_z    = 0, # Constant mean for other model variables
  AnimalID             = dat_clean1$AnimalID[1],
  SessionID            = dat_clean1$SessionID[1]
)

# 3. Log-transform and Z-score columns dynamically to match model requirements
# (Replace log() with log10() if your original transformation used base 10)
pred_grid_ti$log_TimefromLastRipple <- log10(pred_grid_ti$TimefromLastRipple)
pred_grid_ti$log_nextDOWNDuration   <- log10(pred_grid_ti$nextDOWNDuration)

pred_grid_ti$log_TimefromLastRipple_z <- (
  pred_grid_ti$log_TimefromLastRipple - mean(dat_clean1$log_TimefromLastRipple, na.rm = TRUE)
) / sd(dat_clean1$log_TimefromLastRipple, na.rm = TRUE)

pred_grid_ti$log_nextDOWNDuration_z <- (
  pred_grid_ti$log_nextDOWNDuration - mean(dat_clean1$log_nextDOWNDuration, na.rm = TRUE)
) / sd(dat_clean1$log_nextDOWNDuration, na.rm = TRUE)

# 4. Extract specific term components matrix using the z-scores
term_preds <- predict(mdl_final, newdata = pred_grid_ti, type = "terms")

# Isolate the pure interaction tensor term
pred_grid_ti$Interaction_Effect <- term_preds[, "ti(log_TimefromLastRipple_z,log_nextDOWNDuration_z)"]

# 5. Plot Pure Interaction Surface using RAW scales for X and Y axes
p_interaction_raw_scale <- ggplot(
  pred_grid_ti, 
  aes(x = TimefromLastRipple, y = nextDOWNDuration, fill = Interaction_Effect)
) +
  geom_tile() + 
  geom_contour(aes(z = Interaction_Effect), color = "black", alpha = 0.2) + 
  scale_fill_gradient2(
    low = "dodgerblue", 
    mid = "white", 
    high = "firebrick", 
    midpoint = 0, 
    limits = c(0, 0.155), 
    oob = scales::squish,
    name = "Effect"
  ) +
  theme_minimal() +
  theme(aspect.ratio = 1) +
  labs(
    title = "Pure Tensor Interaction (ti Term)",
    subtitle = "Variance unique to the combination (Raw Scale Mapping)",
    x = "Time from Last Ripple (Raw)", 
    y = "Next DOWN Duration (Raw)"
  )

print(p_interaction_raw_scale)

# 6. Save figure
cairo_pdf("TimefromLastRipple_DOWNDuration_ti_lateUP_nextUP_V1_cohernece.pdf", width = 4.3, height = 4.3)
print(p_interaction_raw_scale)
dev.off()

# 
# cairo_pdf("DOWNPowerLagInteraction_lastRipplenextUPCoherence.pdf", width = 4.3, height = 4.3)
# print(p_interaction_raw_scale)
# dev.off()



# ==============================================================================
# --- Total : log_TimefromLastRipple_z and log_nextDOWNDuration_z ---
# ==============================================================================

# 1. Define sequence ranges using your actual RAW data minimums and maximums
# (Adjust the min/max values to match your dataset's raw range)
# raw_range_ripple_time <- seq(
#   min(dat_clean1$TimefromLastRipple, na.rm = TRUE),
#   max(dat_clean1$TimefromLastRipple, na.rm = TRUE),
#   length.out = 100
# )
# 
# raw_range_down_dur <- seq(
#   min(dat_clean1$nextDOWNDuration, na.rm = TRUE),
#   max(dat_clean1$nextDOWNDuration, na.rm = TRUE),
#   length.out = 100
# )

raw_range_ripple_time   <- seq(0.02,
                               0.25, length.out = 100)

raw_range_down_dur  <- seq(0.02,
                           0.25, length.out = 100)

# 2. Build the prediction grid using the RAW scales
pred_grid_ti <- expand.grid(
  TimefromLastRipple   = raw_range_ripple_time,
  nextDOWNDuration     = raw_range_down_dur,
  nextDOWNSOPower_z    = 0, # Constant mean for other model variables
  nextDOWNlag_z        = 0, # Constant mean for other model variables
  lastRipplePower_z    = 0, # Constant mean for other model variables
  AnimalID             = dat_clean1$AnimalID[1],
  SessionID            = dat_clean1$SessionID[1]
)

# 3. Log-transform and Z-score columns dynamically to match model requirements
# (Replace log() with log10() if your original transformation used base 10)
pred_grid_ti$log_TimefromLastRipple <- log10(pred_grid_ti$TimefromLastRipple)
pred_grid_ti$log_nextDOWNDuration   <- log10(pred_grid_ti$nextDOWNDuration)

pred_grid_ti$log_TimefromLastRipple_z <- (
  pred_grid_ti$log_TimefromLastRipple - mean(dat_clean1$log_TimefromLastRipple, na.rm = TRUE)
) / sd(dat_clean1$log_TimefromLastRipple, na.rm = TRUE)

pred_grid_ti$log_nextDOWNDuration_z <- (
  pred_grid_ti$log_nextDOWNDuration - mean(dat_clean1$log_nextDOWNDuration, na.rm = TRUE)
) / sd(dat_clean1$log_nextDOWNDuration, na.rm = TRUE)

# 4. Extract specific term components matrix using the z-scores
term_preds <- predict(mdl_final, newdata = pred_grid_ti, type = "terms")

# Isolate the pure interaction tensor term
# pred_grid_ti$Interaction_Effect <- term_preds[, "ti(log_TimefromLastRipple_z,log_nextDOWNDuration_z)"]
pred_grid_ti$Interaction_Effect <- 
  term_preds[, "s(log_TimefromLastRipple_z)"] + 
  term_preds[, "s(log_nextDOWNDuration_z)"] + 
  term_preds[, "ti(log_TimefromLastRipple_z,log_nextDOWNDuration_z)"]

# 5. Plot Pure Interaction Surface using RAW scales for X and Y axes
p_interaction_raw_scale <- ggplot(
  pred_grid_ti, 
  aes(x = TimefromLastRipple, y = nextDOWNDuration, fill = Interaction_Effect)
) +
  geom_tile() + 
  geom_contour(aes(z = Interaction_Effect), color = "black", alpha = 0.1) + 
  scale_fill_gradient2(
    low = "dodgerblue", 
    mid = "white", 
    high = "firebrick", 
    midpoint = 0, 
    limits = c(0, 0.155), 
    oob = scales::squish,
    name = "Effect"
  ) +
  theme_minimal() +
  theme(aspect.ratio = 1) +
  labs(
    title = "Pure Tensor Interaction (ti + s + s Term)",
    subtitle = "Total Variance (Raw Scale Mapping)",
    x = "Time from Last Ripple (Raw)", 
    y = "Next DOWN Duration (Raw)"
  )

print(p_interaction_raw_scale)

# 6. Save figure
cairo_pdf("TimefromLastRipple_DOWNDuration_total_lateUP_nextUP_V1_cohernece.pdf", width = 4.3, height = 4.3)
print(p_interaction_raw_scale)
dev.off()




# ==============================================================================
# --- Total : log_TimefromLastRipple_z and log_nextDOWNDuration_z ---
# ==============================================================================

# 1. Define sequence ranges using your actual RAW data minimums and maximums
# (Adjust the min/max values to match your dataset's raw range)
# raw_range_ripple_time <- seq(
#   min(dat_clean1$TimefromLastRipple, na.rm = TRUE),
#   max(dat_clean1$TimefromLastRipple, na.rm = TRUE),
#   length.out = 100
# )
# 
# raw_range_down_dur <- seq(
#   min(dat_clean1$nextDOWNDuration, na.rm = TRUE),
#   max(dat_clean1$nextDOWNDuration, na.rm = TRUE),
#   length.out = 100
# )

raw_range_ripple_time   <- seq(0.02,
                               0.3, length.out = 100)

raw_range_down_dur  <- seq(0.02,
                           0.3, length.out = 100)

# 2. Build the prediction grid using the RAW scales
pred_grid_ti <- expand.grid(
  TimefromLastRipple   = raw_range_ripple_time,
  nextDOWNDuration     = raw_range_down_dur,
  nextDOWNSOPower_z    = 0, # Constant mean for other model variables
  nextDOWNlag_z        = 0, # Constant mean for other model variables
  lastRipplePower_z    = 0, # Constant mean for other model variables
  AnimalID             = dat_clean1$AnimalID[1],
  SessionID            = dat_clean1$SessionID[1]
)

# 3. Log-transform and Z-score columns dynamically to match model requirements
# (Replace log() with log10() if your original transformation used base 10)
pred_grid_ti$log_TimefromLastRipple <- log10(pred_grid_ti$TimefromLastRipple)
pred_grid_ti$log_nextDOWNDuration   <- log10(pred_grid_ti$nextDOWNDuration)

pred_grid_ti$log_TimefromLastRipple_z <- (
  pred_grid_ti$log_TimefromLastRipple - mean(dat_clean1$log_TimefromLastRipple, na.rm = TRUE)
) / sd(dat_clean1$log_TimefromLastRipple, na.rm = TRUE)

pred_grid_ti$log_nextDOWNDuration_z <- (
  pred_grid_ti$log_nextDOWNDuration - mean(dat_clean1$log_nextDOWNDuration, na.rm = TRUE)
) / sd(dat_clean1$log_nextDOWNDuration, na.rm = TRUE)

# 4. Extract specific term components matrix using the z-scores
term_preds <- predict(mdl_final, newdata = pred_grid_ti, type = "terms")

# Isolate the pure interaction tensor term
# pred_grid_ti$Interaction_Effect <- term_preds[, "ti(log_TimefromLastRipple_z,log_nextDOWNDuration_z)"]
pred_grid_ti$Interaction_Effect <- 
  term_preds[, "s(log_TimefromLastRipple_z)"] + 
  term_preds[, "s(log_nextDOWNDuration_z)"]
  # term_preds[, "ti(log_TimefromLastRipple_z,log_nextDOWNDuration_z)"]

# 5. Plot Pure Interaction Surface using RAW scales for X and Y axes
p_interaction_raw_scale <- ggplot(
  pred_grid_ti, 
  aes(x = TimefromLastRipple, y = nextDOWNDuration, fill = Interaction_Effect)
) +
  geom_tile() + 
  geom_contour(aes(z = Interaction_Effect), color = "black", alpha = 0.2) + 
  scale_fill_gradient2(
    low = "dodgerblue", 
    mid = "white", 
    high = "firebrick", 
    midpoint = 0, 
    name = "Effect"
  ) +
  theme_minimal() +
  theme(aspect.ratio = 1) +
  labs(
    title = "Pure Tensor Interaction (s + s Term)",
    subtitle = "Total Variance (Raw Scale Mapping)",
    x = "Time from Last Ripple (Raw)", 
    y = "Next DOWN Duration (Raw)"
  )

print(p_interaction_raw_scale)

# 6. Save figure
cairo_pdf("TimefromLastRipple_DOWNDuration_total_without_ti_lateUP_nextUP_V1_cohernece.pdf", width = 4.3, height = 4.3)
print(p_interaction_raw_scale)
dev.off()

# 
# # ==============================================================================
# # --- Total effect for next DOWN SO power and next DOWN bilateral lag ---
# # ==============================================================================
# # 1. Define sequence ranges using your actual RAW data minimums and maximums
# # (Adjust the min/max or length if you want rounded limits like seq(0, 100, by=1))
# raw_range_so   <- seq(1, 
#                       3, length.out = 100)
# 
# raw_range_lag  <- seq(0, 
#                       0.2, length.out = 100)
# 
# # 2. Build the prediction grid using the RAW scales
# pred_grid_reconstruct <- expand.grid(
#   nextDOWNSOPower   = raw_range_so,
#   nextDOWNlag       = raw_range_lag,
#   lastRipplePower_z = 0,   # Constant mean for other model variables
#   log_TimefromLastRipple_z = 0,   # Constant mean for other model variables
#   AnimalID           = dat_clean1$AnimalID[1],
#   SessionID          = dat_clean1$SessionID[1]
# )
# 
# # 3. Create the Z-scored columns that the model expects 
# # This dynamically matches how your data was scaled (Mean=0, SD=1)
# pred_grid_reconstruct$nextDOWNSOPower_z <- (pred_grid_reconstruct$nextDOWNSOPower - mean(dat_clean1$nextDOWNSOPower, na.rm=TRUE)) / sd(dat_clean1$nextDOWNSOPower, na.rm=TRUE)
# pred_grid_reconstruct$nextDOWNlag_z     <- (pred_grid_reconstruct$nextDOWNlag     - mean(dat_clean1$nextDOWNlag, na.rm=TRUE))     / sd(dat_clean1$nextDOWNlag, na.rm=TRUE)
# 
# # 4. Extract specific term components matrix using the z-scores
# term_preds <- predict(mdl_final, newdata = pred_grid_reconstruct, type = "terms")
# 
# # Isolate the pure interaction tensor term
# # pred_grid_ti$Interaction_Effect <- term_preds[, "ti(nextDOWNSOPower_z,nextDOWNlag_z)"]
# pred_grid_reconstruct$Interaction_Effect <- 
#   term_preds[, "s(nextDOWNSOPower_z)"] + 
#   term_preds[, "s(nextDOWNlag_z)"] + 
#   term_preds[, "ti(nextDOWNSOPower_z,nextDOWNlag_z)"]
# 
# # 5. Plot Pure Interaction Surface using RAW scales for X and Y axes
# p_interaction_raw_scale <- ggplot(pred_grid_reconstruct, aes(x = nextDOWNSOPower, y = nextDOWNlag, fill = Interaction_Effect)) +
#   geom_tile() + 
#   geom_contour(aes(z = Interaction_Effect), color = "black", alpha = 0.2) + 
#   scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "firebrick", midpoint = 0, name = "Effect") +
#   theme_minimal() +
#   theme(aspect.ratio = 1) +
#   labs(
#     title = "Pure Tensor Interaction (ti Term)",
#     subtitle = "Variance unique to the combination (Raw Scale Mapping)",
#     x = "nextDOWNSOPower (Raw)", 
#     y = "nextDOWNlag (Raw)"
#   )
# print(p_interaction_raw_scale)
# 
# cairo_pdf("DOWNPowerLagTotal_lastRipplenextUPCoherence.pdf", width = 4.3, height = 4.3)
# print(p_interaction_raw_scale)
# dev.off()
# 




# ==============================================================================
# --- ti interaction: lastRippleMUArate_z and geo_coherencePRE_z ---
# ==============================================================================

# 1. Define sequence ranges using your actual RAW data minimums and maximums
# Adjust these to match dat_clean1$lastRippleMUArate and dat_clean1$geo_coherencePRE
raw_range_MUArate <- seq(
  min(dat_clean1$lastRippleMUArate, na.rm = TRUE),
  max(dat_clean1$lastRippleMUArate, na.rm = TRUE),
  length.out = 100
)

raw_range_coherence <- seq(
  min(dat_clean1$geo_coherencePRE, na.rm = TRUE),
  max(dat_clean1$geo_coherencePRE, na.rm = TRUE),
  length.out = 100
)

# 2. Build the prediction grid using the RAW scales
pred_grid_ti <- expand.grid(
  lastRippleMUArate    = raw_range_MUArate,
  geo_coherencePRE     = raw_range_coherence,
  # Add any OTHER predictors in mdl_final here, held at their mean/reference level
  # e.g. nextDOWNSOPower_z = 0, nextDOWNlag_z = 0, etc.
  AnimalID              = dat_clean1$AnimalID[1],
  SessionID             = dat_clean1$SessionID[1]
)

# 3. Z-score columns dynamically to match model requirements
# (No log transform applied here — remove/add log10() if your model actually used it)
pred_grid_ti$lastRippleMUArate_z <- (
  pred_grid_ti$lastRippleMUArate - mean(dat_clean1$lastRippleMUArate, na.rm = TRUE)
) / sd(dat_clean1$lastRippleMUArate, na.rm = TRUE)

pred_grid_ti$geo_coherencePRE_z <- (
  pred_grid_ti$geo_coherencePRE - mean(dat_clean1$geo_coherencePRE, na.rm = TRUE)
) / sd(dat_clean1$geo_coherencePRE, na.rm = TRUE)

# 4. Extract specific term components matrix using the z-scores
term_preds <- predict(mdl_final, newdata = pred_grid_ti, type = "terms")

# Isolate the pure interaction tensor term
pred_grid_ti$Interaction_Effect <- term_preds[, "ti(lastRippleMUArate_z,geo_coherencePRE_z)"] + term_preds[,"s(lastRippleMUArate_z)"] 
+term_preds[,"s(geo_coherencePRE_z)"]

# 5. Plot Pure Interaction Surface using RAW scales for X and Y axes
p_interaction_raw_scale <- ggplot(
  pred_grid_ti,
  aes(x = lastRippleMUArate, y = geo_coherencePRE, fill = Interaction_Effect)
) +
  geom_tile() +
  geom_contour(aes(z = Interaction_Effect), color = "black", alpha = 0.2) +
  scale_fill_gradient2(
    low = "dodgerblue",
    mid = "white",
    high = "firebrick",
    midpoint = 0,
    # limits = c(0, 0.155),  # <- recalibrate these to this model's effect range
    oob = scales::squish,
    name = "Effect"
  ) +
  theme_minimal() +
  theme(aspect.ratio = 1) +
  labs(
    title = "Ripple MUA and ripple preV1-postHC coherence predicting preV1-NextUPV1 coherence",
    subtitle = "Variance unique to the combination (Raw Scale Mapping)",
    x = "Last Ripple MUA Rate (Raw)",
    y = "Geodesic Coherence, Pre (Raw)"
  )

print(p_interaction_raw_scale)











# ==============================================================================
# --- 4. MULTI-METRIC EFFECT SIZE CALCULATIONS ---
# ==============================================================================
message("\nCalculating 4 Effect Size Metrics (This may take a minute)...")

# EXACT formula components to rebuild the models robustly
formula_terms <- c(
  "s(lastRipplePower_z, k = 5)",
  "s(log_TimefromLastRipple_z, k = 5)",
  "s(nextDOWNSOPower_z, k = 5)",
  "s(nextDOWNlag_z, k = 5)",
  "ti(nextDOWNSOPower_z, nextDOWNlag_z, k = 5)"
)

# EXACT labels output by summary() and smooth_estimates()
smooth_labels <- c(
  "s(lastRipplePower_z)",
  "s(log_TimefromLastRipple_z)",
  "s(nextDOWNSOPower_z)",
  "s(nextDOWNlag_z)",
  "ti(nextDOWNSOPower_z,nextDOWNlag_z)"
)

library(parallel)
library(pbapply)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(tidyverse)

B <- 1000  # Set to 1000 for final analysis run
RE_TERMS <- c("s(SessionID, bs = 're')", "s(AnimalID, bs = 're')")

message(sprintf("\nLaunching %d Case Bootstrap Replicates...", B))

run_one_bootstrap <- function(rep_id, original_data, formula_terms, smooth_labels) {
  
  # 1. Resample data with replacement
  boot_data <- original_data[sample(nrow(original_data), replace = TRUE), ]
  
  # 2. Fit Full Model
  full_form <- as.formula(paste(
    "geo_coherenceRippleNext_z ~",
    paste(c(formula_terms, RE_TERMS), collapse = " + ")
  ))
  
  mdl_f <- mgcv::bam(full_form, data = boot_data, method = "fREML", discrete = TRUE)
  
  f_sum <- summary(mdl_f)
  full_dev <- f_sum$dev.expl * 100
  res_df <- mdl_f$df.residual
  sum_tab <- as.data.frame(f_sum$s.table)
  sum_tab$Term <- rownames(sum_tab)
  
  # Pre-allocate container for this run's metrics
  run_res <- data.frame(Term = smooth_labels, Partial_Deviance = NA, Eta_Sq_Partial = NA, Amplitude = NA, RMS = NA, Rep = rep_id)
  
  # 3. Process every smooth term
  for(j in 1:length(smooth_labels)) {
    # A. Partial Deviance (Drop-One refit)
    act_terms <- c(formula_terms[-j], RE_TERMS)
    r_form <- as.formula(paste("geo_coherenceRippleNext_z ~", paste(act_terms, collapse = " + ")))
    mdl_r <- mgcv::bam(r_form, data = boot_data, method = "fREML", discrete = TRUE)
    run_res$Partial_Deviance[j] <- full_dev - (summary(mdl_r)$dev.expl * 100)
    
    # B. Partial Eta-Squared
    t_row <- sum_tab[sum_tab$Term == smooth_labels[j], ]
    if(nrow(t_row) == 1) {
      run_res$Eta_Sq_Partial[j] <- (t_row$F * t_row$edf) / ((t_row$F * t_row$edf) + res_df)
    }
    
    # C & D. Amplitude and RMS via completely self-contained custom grid
    vars_in_smooth <- all.vars(as.formula(paste("~", smooth_labels[j])))
    vars_in_smooth <- vars_in_smooth[!vars_in_smooth %in% c("s", "ti", "bs", "k")]
    
    # Build a clean evaluation sequence for variables in this smooth
    slice_args <- lapply(vars_in_smooth, function(v) {
      seq(min(boot_data[[v]], na.rm = TRUE), max(boot_data[[v]], na.rm = TRUE), length.out = 50)
    })
    names(slice_args) <- vars_in_smooth
    
    # Use base R expand.grid to keep data structure strictly local
    grid_clean <- do.call(expand.grid, slice_args)
    
    # Add dummy/mean columns for any other variables that predict() might structurally look for
    all_model_vars <- all.vars(full_form)[-1] 
    missing_vars <- setdiff(all_model_vars, colnames(grid_clean))
    
    for(mv in missing_vars) {
      if(mv == "SessionID") {
        grid_clean[[mv]] <- boot_data$SessionID[1] 
      } else if(mv == "AnimalID") {
        grid_clean[[mv]] <- boot_data$AnimalID[1]
      } else if(is.numeric(boot_data[[mv]])) {
        grid_clean[[mv]] <- mean(boot_data[[mv]], na.rm = TRUE)
      } else {
        grid_clean[[mv]] <- boot_data[[mv]][1]
      }
    }
    
    # Isolate target smooth prediction using the lpmatrix
    Xp <- predict(mdl_f, newdata = grid_clean, type = "lpmatrix")
    smooth_cols <- grep(smooth_labels[j], colnames(Xp), fixed = TRUE)
    
    Xp_isolated <- matrix(0, nrow = nrow(Xp), ncol = ncol(Xp))
    Xp_isolated[, smooth_cols] <- Xp[, smooth_cols]
    
    fit_isolated <- Xp_isolated %*% coef(mdl_f)
    run_res$Amplitude[j] <- max(fit_isolated) - min(fit_isolated)
    run_res$RMS[j]       <- sqrt(mean(fit_isolated^2))
  }
  
  return(run_res)
}

# Run loop with progress bar
set.seed(42)
boot_results_list <- pblapply(1:B, run_one_bootstrap, 
                              original_data = dat_clean1, 
                              formula_terms = formula_terms, 
                              smooth_labels = smooth_labels)

# Clean out any NULLs if any structural issues happen
boot_results_list <- boot_results_list[!sapply(boot_results_list, is.null)]

# Combine results safely
boot_df <- do.call(rbind, boot_results_list)

# Format column headers cleanly before exporting
raw_iterations_clean <- boot_df %>%
  select(Rep, Term, Partial_Deviance, Eta_Sq_Partial, Amplitude, RMS) %>%
  rename(
    Replicate        = Rep,
    Deviance_Value   = Partial_Deviance,
    Eta_Sq_Value     = Eta_Sq_Partial,
    Amplitude_Value  = Amplitude,
    RMS_Value        = RMS
  ) %>%
  arrange(Term, Replicate)

write.csv(raw_iterations_clean, "LastRippleNextUPCoherence_GAM_model_raw_bootstrap_iterations.csv", row.names = FALSE)
message("Raw bootstrap iteration file saved successfully: 'LastRippleNextUP_GAM_model_raw_bootstrap_iterations.csv'")

# Load raw bootstrap iteration results from CSV
raw_iterations_clean <- read.csv("LastRippleNextUPCoherence_GAM_model_raw_bootstrap_iterations.csv",
                                 stringsAsFactors = FALSE)

# ==============================================================================
# --- 5. GENERATE PLOTS AND EXPORT SUMMARIES ---
# ==============================================================================
plot_data <- raw_iterations_clean %>%
  group_by(Term) %>%
  summarise(
    Deviance_Value__Val = median(Deviance_Value, na.rm = TRUE),
    Deviance_Value__Lwr = quantile(Deviance_Value, probs = 0.025, na.rm = TRUE),
    Deviance_Value__Upr = quantile(Deviance_Value, probs = 0.975, na.rm = TRUE),
    
    Eta_Sq_Value__Val   = median(Eta_Sq_Value, na.rm = TRUE),
    Eta_Sq_Value__Lwr   = quantile(Eta_Sq_Value, probs = 0.025, na.rm = TRUE),
    Eta_Sq_Value__Upr   = quantile(Eta_Sq_Value, probs = 0.975, na.rm = TRUE),
    
    Amplitude_Value__Val = median(Amplitude_Value, na.rm = TRUE),
    Amplitude_Value__Lwr = quantile(Amplitude_Value, probs = 0.025, na.rm = TRUE),
    Amplitude_Value__Upr = quantile(Amplitude_Value, probs = 0.975, na.rm = TRUE),
    
    RMS_Value__Val       = median(RMS_Value, na.rm = TRUE),
    RMS_Value__Lwr       = quantile(RMS_Value, probs = 0.025, na.rm = TRUE),
    RMS_Value__Upr       = quantile(RMS_Value, probs = 0.975, na.rm = TRUE)
  ) %>%
  pivot_longer(cols = -Term, names_to = "Metric_Stat", values_to = "Val") %>%
  separate(Metric_Stat, into = c("Metric", "Stat"), sep = "__") %>%
  pivot_wider(names_from = Stat, values_from = Val) %>%
  mutate(
    Metric = factor(Metric,
                    levels = c("Deviance_Value", "Eta_Sq_Value", "Amplitude_Value", "RMS_Value"),
                    labels = c("Deviance Explained (%)", "Partial Eta-Squared",
                               "Peak-to-Trough Amplitude", "RMS Effect"))
  )

p_bars_with_ci <- ggplot(plot_data, aes(x = reorder(Term, Val), y = Val, fill = Term)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60", alpha = 0.5) +
  geom_bar(stat = "identity", show.legend = FALSE, alpha = 0.85, width = 0.75) +
  geom_errorbar(aes(ymin = Lwr, ymax = Upr), width = 0.25, color = "black", linewidth = 0.6) +
  coord_flip() +
  facet_wrap(~Metric, scales = "free_x", ncol = 2) +
  scale_fill_viridis_d(option = "mako", direction = -1) +
  theme_bw() +
  labs(
    title = "GAMM Effect Size & Variance Metrics",
    subtitle = "Bars denote bootstrap medians and 95% non-parametric CIs.",
    x = NULL, y = NULL
  ) +
  theme(
    strip.text = element_text(face = "bold", size = 10, colour = "black"),
    axis.text.y = element_text(size = 9, color = "black"),
    panel.spacing = unit(1.2, "lines")
  )

print(p_bars_with_ci)

cairo_pdf("LastRippleNextUPCoherence_Model_Effect_Sizes_With_CI.pdf", width = 10, height = 6.5)
print(p_bars_with_ci)
dev.off()

final_dashboard_data <- raw_iterations_clean %>%
  group_by(Term) %>%
  summarise(across(c(Deviance_Value, Eta_Sq_Value, Amplitude_Value, RMS_Value),
                   list(Val = ~median(.x, na.rm = TRUE),
                        Lwr = ~quantile(.x, probs = 0.025, na.rm = TRUE),
                        Upr = ~quantile(.x, probs = 0.975, na.rm = TRUE)),
                   .names = "{.col}__{.fn}"))

flat_bootstrap_results <- final_dashboard_data %>%
  pivot_longer(cols = -Term, names_to = "Combined", values_to = "Value") %>%
  separate(Combined, into = c("Metric", "Stat"), sep = "__") %>%
  mutate(New_Col_Name = paste0(Metric, "_", Stat)) %>%
  select(-Metric, -Stat) %>%
  pivot_wider(names_from = New_Col_Name, values_from = Value)

model_stats <- as.data.frame(summary(mdl_final)$s.table) %>%
  mutate(Term = rownames(.))

combined_results <- model_stats %>%
  left_join(flat_bootstrap_results, by = "Term")

write.csv(combined_results, "LastRippleNextUPCoherence_GAM_model_CI_output.csv", row.names = FALSE)
message("\nAnalysis Pipeline Complete! Saved 'LastRippleNextUPCoherence_GAM_model_CI_output.csv' and figures.")



