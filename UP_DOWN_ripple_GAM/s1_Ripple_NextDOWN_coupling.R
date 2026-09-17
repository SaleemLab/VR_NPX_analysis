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
#dat <- read.csv("C:/Users/masah/Documents/GitHub/VR_NPX_analysis/UP_DOWN_ripple_GAM/UP_DOWN_info_GAM_offset.csv")

#dat <- read.csv("C:/Users/masah/Documents/GitHub/VR_NPX_analysis/UP_DOWN_ripple_GAM/UP_DOWN_info_GAM_peak.csv")

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

# # List of the new Z-score columns to filter (Updated to include all relevant ripple predictors)
# z_cols <- c(
#   "lastRippleV1_z", "lastRippleV1PRE_z",
#   "lastRippleHPC_z", "lastRippleHPCPRE_z",
#   "firstRippleV1_z", "firstRippleV1PRE_z",
#   "firstRippleHPC_z", "firstRippleHPCPRE_z"
# )




my_folder <- "C:/Users/masah/Documents/GitHub/VR_NPX_analysis/UP_DOWN_ripple_GAM/lastRipple_DOWN_coupling"

# 2. Check if it exists; if not, create it
if (!dir.exists(my_folder)) {
  dir.create(my_folder, recursive = TRUE)
}

# 3. Change the working directory to that folder
setwd(my_folder)


z_cols <- c(
  "nextDOWNlag_z","nextDOWNSOPower_z","lastRipplePower_z",
  "TimefromLastRipple_z"
)

# dat_clean1 <- dat_clean %>%
#   filter(
#     if_all(all_of(z_cols), ~ abs(.) < z_thresh | is.na(.))
#   )

dat_clean1 <- dat_clean %>%
  filter(
    if_all(all_of(z_cols), ~ abs(.) < z_thresh | is.na(.)),
    #nextDOWNlag <0.3,
    # nextDOWNlagSigned <= 0
    # nextDOWNlagSigned >= 0
  )


mdl_final <- bam(nextDOWNSOPower_z ~ 
                   s(nextDOWNlag_z, k = 5) +
                   # s(nextDOWNDuration_z, k = 5) +
                   # ti(lastRipplePower_z, nextDOWNSOPower_z, k = 5)+
                   
                   # s(lastRipplePower_z, k = 5) +
                   
                   # 4. Control Term
                   s(AnimalID, bs = "re") +
                   s(SessionID, bs = "re"), 
                 
                 data = dat_clean1 , 
                 method = "fREML", 
                 discrete = TRUE, 
                 nthreads = 4)

#message("\n--- FINAL MODEL SUMMARY ---")
print(summary(mdl_final))


### nextDOWN lag -> nextDOWN power
# Calculate scaling factors
# raw_breaks <- c(0,0.1,0.2,0.3,0.4,0.5)
raw_breaks <- c(0,0.1,0.2,0.3)
z_breaks <- sapply(raw_breaks, function(val) {
  dat_clean1$nextDOWNlag_z[which.min(abs(dat_clean1$nextDOWNlag - val))]
})


p_lag_raw <- draw(mdl_final, select = "s(nextDOWNlag_z)", residuals = FALSE, rug = FALSE) +
  theme_bw(base_family = "Arial") +
  theme(aspect.ratio = 1) +
  scale_x_continuous(breaks = z_breaks, labels = raw_breaks) +
  # coord_cartesian() +
  coord_cartesian(xlim = c(-1.2,2.3),ylim = c(-0.17,0.27),expand = FALSE) +

  labs(
    title = "Next DOWN lag predicts next DOWN SO power",
    x = "DOWN lag",
    y = "Partial Effect"
  )
# dev.new(noRStudioGD = TRUE)
cairo_pdf("DownLag_and_SOPower.pdf", width = 4.3, height = 4.3)
print(p_lag_raw)
dev.off()

# 
# 
# ### nextDOWN duration -> nextDOWN power
# # Calculate scaling factors
# raw_breaks <- c(0,0.1,0.2,0.3)
# z_breaks <- sapply(raw_breaks, function(val) {
#   dat_clean1$nextDOWNDuration_z[which.min(abs(dat_clean1$nextDOWNDuration - val))]
# })
# 
# # cairo_pdf("lateUPHPC_nextUP", width = 4.3, height = 4.3)
# p_lag_raw <- draw(mdl_final, select = "s(nextDOWNDuration_z)", residuals = FALSE, rug = FALSE) +
#   theme_bw(base_family = "Arial") +
#   theme(aspect.ratio = 1) +
#   scale_x_continuous(breaks = z_breaks, labels = raw_breaks) +
#   # coord_cartesian(xlim = c(-1.2,2.3)) +
#   labs(
#     title = "Next DOWN duration predicts next DOWN SO power",
#     x = "DOWN lag",
#     y = "Partial Effect"
#   )
# dev.new(noRStudioGD = TRUE)
# print(p_lag_raw)

# 


######
mdl_final <- bam(nextDOWNlag_z ~
                   s(nextDOWNSOPower_z, k = 5) +
                   # ti(lastRipplePower_z, nextDOWNSOPower_z, k = 5)+

                   # s(lastRipplePower_z, k = 5) +

                   # 4. Control Term
                   s(AnimalID, bs = "re") +
                   s(SessionID, bs = "re"),

                 data = dat_clean1 ,
                 method = "fREML",
                 discrete = TRUE,
                 nthreads = 4)

#message("\n--- FINAL MODEL SUMMARY ---")
print(summary(mdl_final))

### nextDOWN power -> nextDOWN lag
# Calculate scaling factors
raw_breaks <- c(1,2,3,4)
z_breaks <- sapply(raw_breaks, function(val) {
  dat_clean1$nextDOWNSOPower_z[which.min(abs(dat_clean1$nextDOWNSOPower - val))]
})

# cairo_pdf("lateUPHPC_nextUP", width = 4.3, height = 4.3)
p_lag_raw <- draw(mdl_final, select = "s(nextDOWNSOPower_z)", residuals = FALSE, rug = FALSE) +
  theme_bw(base_family = "Arial") +
  theme(aspect.ratio = 1) +
  scale_x_continuous(breaks = z_breaks, labels = raw_breaks) +
  # coord_cartesian(xlim = c(-1.2,3)) +
  labs(
    title = "Next DOWN power predicts next DOWN lag",
    x = "DOWN SO power",
    y = "Partial Effect"
  )
cairo_pdf("SOPower_predicts_DownLag.pdf", width = 4.3, height = 4.3)
# dev.new(noRStudioGD = TRUE)
print(p_lag_raw)
dev.off()


# 
# ######
# ###### DOWN SO power + ripple power -> DOWN lag
# ######
# 
# mdl_final <- bam(nextDOWNlag_z ~ 
#                    s(nextDOWNSOPower_z, k = 5) +
#                    # s(nextDOWNDuration_z, k = 5) +
#                    ti(lastRipplePower_z, nextDOWNSOPower_z, k = 5)+
#                    
#                    s(lastRipplePower_z, k = 5) +
#                    
#                    # 4. Control Term
#                    s(AnimalID, bs = "re") +
#                    s(SessionID, bs = "re"), 
#                  
#                  data = dat_clean1 , 
#                  method = "fREML", 
#                  discrete = TRUE, 
#                  nthreads = 4)
# 
# #message("\n--- FINAL MODEL SUMMARY ---")
# print(summary(mdl_final))
# 
# 
# 
# # ==============================================================================
# # --- 2D Ti of next DOWN SO power and last ripple  ---
# # ==============================================================================
# 
# # 1. Define sequence ranges using your actual RAW data minimums and maximums
# raw_range_so   <- seq(1, 4, length.out = 100)
# raw_range_ripplePower  <- seq(5, 20, length.out = 100)
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
# pred_grid_te$Interaction_Effect <- 
# term_preds[, "s(lastRipplePower_z)"] + 
#   term_preds[, "s(nextDOWNSOPower_z)"]
#   term_preds[, "ti(lastRipplePower_z,nextDOWNSOPower_z)"]
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



####### 
####### Last ripple power and last ripple time from DOWN -> Next DOWN SO power 
#######
z_cols <- c(
  "nextDOWNlag_z","nextDOWNSOPower_z","lastRipplePower_z",
  "TimefromLastRipple_z"
)

dat_clean1 <- dat_clean %>%
  filter(
    if_all(all_of(z_cols), ~ abs(.) < z_thresh | is.na(.)),
    #nextDOWNlag <0.15,
  )


####
####
####  
mdl_final <- bam(nextDOWNSOPower_z ~ 
                     s(lastRipplePower_z,k=5)+
                     s(TimefromLastRipple_z, k = 5) +
                     # 
                     ti(lastRipplePower_z, TimefromLastRipple_z, k = 5)+
                   # te(lastRipplePower_z, TimefromLastRipple_z, k = 5)+
                     
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



### last ripple power -> SO power
# Calculate scaling factors
raw_breaks <- c(5,10,15,20)
z_breaks <- sapply(raw_breaks, function(val) {
  dat_clean1$lastRipplePower_z[which.min(abs(dat_clean1$lastRipplePower - val))]
})

# cairo_pdf("lateUPHPC_nextUP", width = 4.3, height = 4.3)
p_lag_raw <- draw(mdl_final, select = "s(lastRipplePower_z)", residuals = FALSE, rug = FALSE) + 
  theme_bw(base_family = "Arial") + 
  theme(aspect.ratio = 1) +
  scale_x_continuous(breaks = z_breaks, labels = raw_breaks) +
  coord_cartesian(xlim = c(-1.2,3)) + 
  labs(
    title = "Ripple power effect on nextUP V1 SO", 
    x = "Ripple power", 
    y = "Partial Effect"
  )
# cairo_pdf("lastRipplePower_predicts_SOPower_single", width = 4.3, height = 4.3)
dev.new(noRStudioGD = TRUE)
print(p_lag_raw)
# dev.off()

### Time from last ripple
# Calculate scaling factors
# Calculate scaling factors
raw_breaks <- c(0,0.1,0.2,0.3,0.4,0.5)
z_breaks <- sapply(raw_breaks, function(val) {
  dat_clean1$TimefromLastRipple_z[which.min(abs(dat_clean1$TimefromLastRipple - val))]
})

# cairo_pdf("lateUPHPC_nextUP", width = 4.3, height = 4.3)
p_lag_raw <- draw(mdl_final, select = "s(TimefromLastRipple_z)", residuals = FALSE, rug = FALSE) + 
  theme_bw(base_family = "Arial") + 
  theme(aspect.ratio = 1) +
  scale_x_continuous(breaks = z_breaks, labels = raw_breaks) +
  coord_cartesian(xlim = c(-0.845,0.49)) + 
  labs(
    title = "TimefromLastRipple on next DOWN power", 
    x = "TimefromLastRipple", 
    y = "Partial Effect"
  )

cairo_pdf("TimefromLastRipple_predicts_DownSOPower.pdf", width = 4.3, height = 4.3)
# dev.new(noRStudioGD = TRUE)
print(p_lag_raw)
dev.off()



# ==============================================================================
# --- interaction between ripple power and last ripple to DOWN time ---
# ==============================================================================
# 1. Define sequence ranges using your actual RAW data minimums and maximums
# (Adjust the min/max or length if you want rounded limits like seq(0, 100, by=1))
raw_range_RipplePower   <- seq(5, 
                      20, length.out = 100)

raw_range_Time  <- seq(0, 
                      0.5, length.out = 100)

# 2. Build the prediction grid using the RAW scales
pred_grid_ti <- expand.grid(
  lastRipplePower   = raw_range_RipplePower,
  TimefromLastRipple       = raw_range_Time,
  # nextDOWNDuration_z = 0,   # Constant mean for other model variables
  AnimalID           = dat_clean1$AnimalID[1],
  SessionID          = dat_clean1$SessionID[1]
)

# 3. Create the Z-scored columns that the model expects 
# This dynamically matches how your data was scaled (Mean=0, SD=1)
pred_grid_ti$lastRipplePower_z <- (pred_grid_ti$lastRipplePower - mean(dat_clean1$lastRipplePower, na.rm=TRUE)) / sd(dat_clean1$lastRipplePower, na.rm=TRUE)
pred_grid_ti$TimefromLastRipple_z     <- (pred_grid_ti$TimefromLastRipple     - mean(dat_clean1$TimefromLastRipple, na.rm=TRUE))     / sd(dat_clean1$TimefromLastRipple, na.rm=TRUE)

# 4. Extract specific term components matrix using the z-scores
term_preds <- predict(mdl_final, newdata = pred_grid_ti, type = "terms")

# Isolate the pure interaction tensor term
pred_grid_ti$Interaction_Effect <- term_preds[, "ti(lastRipplePower_z,TimefromLastRipple_z)"]

# 5. Plot Pure Interaction Surface using RAW scales for X and Y axes
p_interaction_raw_scale <- ggplot(pred_grid_ti, aes(x = TimefromLastRipple, y = lastRipplePower , fill = Interaction_Effect)) +
  geom_tile() + 
  geom_contour(aes(z = Interaction_Effect), color = "black", alpha = 0.2) + 
  scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "firebrick", midpoint = 0, name = "Effect") +
  theme_minimal() +
  theme(aspect.ratio = 1) +
  labs(
    title = "Pure Tensor Interaction (ti Term)",
    subtitle = "Variance unique to the combination (Raw Scale Mapping)",
    y = "lastRipplePower_z (Raw)", 
    x = "TimefromLastRipple_z (Raw)"
  )


cairo_pdf("TimefromLastRipple_and_lastRipplePower_ti_predicts_DownSOPower.pdf", width = 4.3, height = 4.3)
# dev.new(noRStudioGD = TRUE)
print(p_interaction_raw_scale)
dev.off()

# ==============================================================================
# --- ti + s Total effect ---
# ==============================================================================
# 1. Define sequence ranges using your actual RAW data minimums and maximums
# (Adjust the min/max or length if you want rounded limits like seq(0, 100, by=1))
raw_range_RipplePower   <- seq(5, 
                               20, length.out = 100)

raw_range_Time  <- seq(0, 
                       0.5, length.out = 100)

# 2. Build the prediction grid using the RAW scales
pred_grid_reconstruct <- expand.grid(
  lastRipplePower   = raw_range_RipplePower,
  TimefromLastRipple       = raw_range_Time,
  # nextDOWNDuration_z = 0,   # Constant mean for other model variables
  AnimalID           = dat_clean1$AnimalID[1],
  SessionID          = dat_clean1$SessionID[1]
)

# 3. Create the Z-scored columns that the model expects 
# This dynamically matches how your data was scaled (Mean=0, SD=1)
pred_grid_reconstruct$lastRipplePower_z <- (pred_grid_reconstruct$lastRipplePower - mean(dat_clean1$lastRipplePower, na.rm=TRUE)) / sd(dat_clean1$lastRipplePower, na.rm=TRUE)
pred_grid_reconstruct$TimefromLastRipple_z     <- (pred_grid_reconstruct$TimefromLastRipple     - mean(dat_clean1$TimefromLastRipple, na.rm=TRUE))     / sd(dat_clean1$TimefromLastRipple, na.rm=TRUE)

# 4. Extract specific term components matrix using the z-scores
term_preds <- predict(mdl_final, newdata = pred_grid_reconstruct, type = "terms")

# Isolate the pure interaction tensor term
# pred_grid_ti$Interaction_Effect <- term_preds[, "ti(lastRipplePower_z,TimefromLastRipple_z)"]
pred_grid_reconstruct$Reconstructed_effect <- 
  term_preds[, "s(lastRipplePower_z)"] + 
  term_preds[, "s(TimefromLastRipple_z)"] + 
  term_preds[, "ti(lastRipplePower_z,TimefromLastRipple_z)"]

# 5. Plot Pure Interaction Surface using RAW scales for X and Y axes
p_interaction_raw_scale <- ggplot(pred_grid_reconstruct, aes(x = TimefromLastRipple, y = lastRipplePower , fill = Reconstructed_effect)) +
  geom_tile() + 
  geom_contour(aes(z = Reconstructed_effect), color = "black", alpha = 0.2) + 
  scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "firebrick", midpoint = 0, name = "Effect") +
  theme_minimal() +
  theme(aspect.ratio = 1) +
  labs(
    title = "Total effect last ripple power and time from last ripple -> next DOWN power",
    subtitle = "Variance unique to the combination (Raw Scale Mapping)",
    y = "lastRipplePower_z (Raw)", 
    x = "TimefromLastRipple_z (Raw)"
  )
dev.new(noRStudioGD = TRUE)
cairo_pdf("TimefromLastRipple_and_lastRipplePower_total_effect_predicts_DownSOPower.pdf", width = 4.3, height = 4.3)

print(p_interaction_raw_scale)
dev.off()





z_cols <- c(
  "nextDOWNlag_z","nextDOWNSOPower_z","lastRipplePower_z",
  "TimefromLastRipple_z"
)

dat_clean1 <- dat_clean %>%
  filter(
    if_all(all_of(z_cols), ~ abs(.) < z_thresh | is.na(.)),
    #nextDOWNlag <0.3,
  )


###
### Next DOWN lag
mdl_final <- bam(nextDOWNlag_z ~ 
                   s(lastRipplePower_z,k=5)+
                   
                   s(TimefromLastRipple_z, k = 5) +
                   # s(nextDOWNlag_z, k = 5) +
                   ti(lastRipplePower_z, TimefromLastRipple_z, k = 5)+
                   
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



### last ripple power -> DOWN lag
# Calculate scaling factors
raw_breaks <- c(5,10,15,20)
z_breaks <- sapply(raw_breaks, function(val) {
  dat_clean1$lastRipplePower_z[which.min(abs(dat_clean1$lastRipplePower - val))]
})

# cairo_pdf("lateUPHPC_nextUP", width = 4.3, height = 4.3)
p_lag_raw <- draw(mdl_final, select = "s(lastRipplePower_z)", residuals = FALSE, rug = FALSE) + 
  theme_bw(base_family = "Arial") + 
  theme(aspect.ratio = 1) +
  scale_x_continuous(breaks = z_breaks, labels = raw_breaks) +
  coord_cartesian(xlim = c(-1.2,3)) + 
  labs(
    title = "Ripple power effect on nextUP V1 SO", 
    x = "Ripple power", 
    y = "Partial Effect"
  )
# cairo_pdf("LastRipplePower_predicts_DownLag_single.pdf", width = 4.3, height = 4.3)
dev.new(noRStudioGD = TRUE)
print(p_lag_raw)
# dev.off()



### Time from last ripple -> DOWN lag
# Calculate scaling factors
# Calculate scaling factors
raw_breaks <- c(0,0.1,0.2,0.3,0.4,0.5)
z_breaks <- sapply(raw_breaks, function(val) {
  dat_clean1$TimefromLastRipple_z[which.min(abs(dat_clean1$TimefromLastRipple - val))]
})

# cairo_pdf("lateUPHPC_nextUP", width = 4.3, height = 4.3)
p_lag_raw <- draw(mdl_final, select = "s(TimefromLastRipple_z)", residuals = FALSE, rug = FALSE) + 
  theme_bw(base_family = "Arial") + 
  theme(aspect.ratio = 1) +
  scale_x_continuous(breaks = z_breaks, labels = raw_breaks) +
  coord_cartesian(xlim = c(-0.845,0.49)) + 
  labs(
    title = "TimefromLastRipple on next DOWN lag", 
    x = "TimefromLastRipple", 
    y = "Partial Effect"
  )
# dev.new(noRStudioGD = TRUE)
cairo_pdf("TimefromLastRipple_predicts_DownSOLag.pdf", width = 4.3, height = 4.3)
print(p_lag_raw)
dev.off()




# ==============================================================================
# --- interaction between ripple power and ripple time to DOWN ---
# ==============================================================================
# 1. Define sequence ranges using your actual RAW data minimums and maximums
# (Adjust the min/max or length if you want rounded limits like seq(0, 100, by=1))
raw_range_RipplePower   <- seq(5, 
                               16, length.out = 100)

raw_range_Time  <- seq(0, 
                       0.5, length.out = 100)

# 2. Build the prediction grid using the RAW scales
pred_grid_ti <- expand.grid(
  lastRipplePower   = raw_range_RipplePower,
  TimefromLastRipple       = raw_range_Time,
  # nextDOWNDuration_z = 0,   # Constant mean for other model variables
  AnimalID           = dat_clean1$AnimalID[1],
  SessionID          = dat_clean1$SessionID[1]
)

# 3. Create the Z-scored columns that the model expects 
# This dynamically matches how your data was scaled (Mean=0, SD=1)
pred_grid_ti$lastRipplePower_z <- (pred_grid_ti$lastRipplePower - mean(dat_clean1$lastRipplePower, na.rm=TRUE)) / sd(dat_clean1$lastRipplePower, na.rm=TRUE)
pred_grid_ti$TimefromLastRipple_z     <- (pred_grid_ti$TimefromLastRipple     - mean(dat_clean1$TimefromLastRipple, na.rm=TRUE))     / sd(dat_clean1$TimefromLastRipple, na.rm=TRUE)

# 4. Extract specific term components matrix using the z-scores
term_preds <- predict(mdl_final, newdata = pred_grid_ti, type = "terms")

# Isolate the pure interaction tensor term
pred_grid_ti$Interaction_Effect <- term_preds[, "ti(lastRipplePower_z,TimefromLastRipple_z)"]

# 5. Plot Pure Interaction Surface using RAW scales for X and Y axes
p_interaction_raw_scale <- ggplot(pred_grid_ti, aes(x = TimefromLastRipple, y = lastRipplePower , fill = Interaction_Effect)) +
  geom_tile() + 
  geom_contour(aes(z = Interaction_Effect), color = "black", alpha = 0.2) + 
  scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "firebrick", midpoint = 0, name = "Effect") +
  scale_y_continuous(breaks = c(5,10,15)) +
  scale_x_continuous(breaks = c(0,0.1,0.2,0.3,0.4,0.5)) +
  theme_minimal() +
  theme(aspect.ratio = 1) +
  labs(
    title = "Pure Tensor Interaction (ti Term)",
    subtitle = "Variance unique to the combination (Raw Scale Mapping)",
    y = "lastRipplePower_z (Raw)", 
    x = "TimefromLastRipple_z (Raw)"
  )

cairo_pdf("TimefromLastRipple_and_lastRipplePower_ti_predicts_DownLag.pdf", width = 4.3, height = 4.3)
print(p_interaction_raw_scale)
dev.off()


# ==============================================================================
# --- ti + s Total effect ---
# ==============================================================================
# 1. Define sequence ranges using your actual RAW data minimums and maximums
# (Adjust the min/max or length if you want rounded limits like seq(0, 100, by=1))
raw_range_RipplePower   <- seq(5, 
                               16, length.out = 100)

raw_range_Time  <- seq(0, 
                       0.5, length.out = 100)

# 2. Build the prediction grid using the RAW scales
pred_grid_reconstruct <- expand.grid(
  lastRipplePower   = raw_range_RipplePower,
  TimefromLastRipple       = raw_range_Time,
  # nextDOWNDuration_z = 0,   # Constant mean for other model variables
  AnimalID           = dat_clean1$AnimalID[1],
  SessionID          = dat_clean1$SessionID[1]
)

# 3. Create the Z-scored columns that the model expects 
# This dynamically matches how your data was scaled (Mean=0, SD=1)
pred_grid_reconstruct$lastRipplePower_z <- (pred_grid_reconstruct$lastRipplePower - mean(dat_clean1$lastRipplePower, na.rm=TRUE)) / sd(dat_clean1$lastRipplePower, na.rm=TRUE)
pred_grid_reconstruct$TimefromLastRipple_z     <- (pred_grid_reconstruct$TimefromLastRipple     - mean(dat_clean1$TimefromLastRipple, na.rm=TRUE))     / sd(dat_clean1$TimefromLastRipple, na.rm=TRUE)

# 4. Extract specific term components matrix using the z-scores
term_preds <- predict(mdl_final, newdata = pred_grid_reconstruct, type = "terms")

# Isolate the pure interaction tensor term
# pred_grid_ti$Interaction_Effect <- term_preds[, "ti(lastRipplePower_z,TimefromLastRipple_z)"]
pred_grid_reconstruct$Reconstructed_effect <- 
  term_preds[, "s(lastRipplePower_z)"] + 
  term_preds[, "s(TimefromLastRipple_z)"] +
  term_preds[, "ti(lastRipplePower_z,TimefromLastRipple_z)"]

# 5. Plot Pure Interaction Surface using RAW scales for X and Y axes
p_interaction_raw_scale <- ggplot(pred_grid_reconstruct, aes(x = TimefromLastRipple, y = lastRipplePower , fill = Reconstructed_effect)) +
  geom_tile() + 
  geom_contour(aes(z = Reconstructed_effect), color = "black", alpha = 0.2) + 
  scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "firebrick", midpoint = 0, name = "Effect") +
  scale_y_continuous(breaks = c(5,10,15)) +
  scale_x_continuous(breaks = c(0,0.1,0.2,0.3,0.4,0.5)) +
  theme_minimal() +
  theme(aspect.ratio = 1) +
  labs(
    title = "Total effect last ripple power and time from last ripple -> next DOWN lag",
    # subtitle = "",
    y = "lastRipplePower_z (Raw)", 
    x = "TimefromLastRipple_z (Raw)"
  )

cairo_pdf("TimefromLastRipple_and_lastRipplePower_total_effect_predicts_DownLag.pdf", width = 4.3, height = 4.3)
print(p_interaction_raw_scale)
dev.off()



############################################################
############################################################

z_cols <- c(
  "nextDOWNlag_z","nextDOWNSOPower_z","lastRipplePower_z",
  "TimefromLastRipple_z"
)

dat_clean1 <- dat_clean %>%
  filter(
    if_all(all_of(z_cols), ~ abs(.) < z_thresh | is.na(.)),
  )



###
### Next DOWN lag and next SO 
mdl_final <- bam(lastRipplePower_z ~ 
                   s(nextDOWNSOPower_z, k = 5) +
                   s(nextDOWNlag_z, k = 5) +
                   ti(nextDOWNSOPower_z, nextDOWNlag_z, k = 5)+
                   
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



### nextDOWN power -> ripple power
# Calculate scaling factors
raw_breaks <- c(1,2,3,4)
z_breaks <- sapply(raw_breaks, function(val) {
  dat_clean1$nextDOWNSOPower_z[which.min(abs(dat_clean1$nextDOWNSOPower - val))]
})

# cairo_pdf("lateUPHPC_nextUP", width = 4.3, height = 4.3)
p_lag_raw <- draw(mdl_final, select = "s(nextDOWNSOPower_z)", residuals = FALSE, rug = FALSE) +
  theme_bw(base_family = "Arial") +
  theme(aspect.ratio = 1) +
  scale_x_continuous(breaks = z_breaks, labels = raw_breaks) +
  # coord_cartesian(xlim = c(-1.2,3)) +
  labs(
    title = "Next DOWN power predicts last ripple power",
    x = "DOWN SO power",
    y = "Partial Effect"
  )
dev.new(noRStudioGD = TRUE)
print(p_lag_raw)


# ==============================================================================
# --- interaction between next DOWN SO power and next DOWN bilateral lag ---
# ==============================================================================
# 1. Define sequence ranges using your actual RAW data minimums and maximums
# (Adjust the min/max or length if you want rounded limits like seq(0, 100, by=1))
raw_range_so   <- seq(1, 
                      3.5, length.out = 100)

raw_range_lag  <- seq(0, 
                      0.2, length.out = 100)

# 2. Build the prediction grid using the RAW scales
pred_grid_ti <- expand.grid(
  nextDOWNSOPower   = raw_range_so,
  nextDOWNlag       = raw_range_lag,

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






# ==============================================================================
# --- total effect of next DOWN SO power and next DOWN bilateral lag ---
# ==============================================================================
# 1. Define sequence ranges using your actual RAW data minimums and maximums
# (Adjust the min/max or length if you want rounded limits like seq(0, 100, by=1))
raw_range_so   <- seq(1, 
                      3.5, length.out = 100)

raw_range_lag  <- seq(0, 
                      0.2, length.out = 100)

# 2. Build the prediction grid using the RAW scales
pred_grid_ti <- expand.grid(
  nextDOWNSOPower   = raw_range_so,
  nextDOWNlag       = raw_range_lag,
  
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
pred_grid_ti$Interaction_Effect <- 
  term_preds[, "s(nextDOWNSOPower_z)"] + 
  #  term_preds[, "s(nextDOWNlag_z)"] +
  term_preds[, "ti(nextDOWNSOPower_z,nextDOWNlag_z)"]


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




