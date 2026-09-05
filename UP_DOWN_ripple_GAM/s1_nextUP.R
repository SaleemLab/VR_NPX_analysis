
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

# List of the new Z-score columns to filter (Updated to include all relevant ripple predictors)
z_cols <- c(
  "lastRippleV1_z", "lastRippleV1PRE_z",
  "lastRippleHPC_z", "lastRippleHPCPRE_z",
  "nextUPV1_z"
  # "firstRippleV1_z", "firstRippleV1PRE_z",
  # "firstRippleHPC_z", "firstRippleHPCPRE_z"
)



###### last ripple close to UP termination predicts nextUP V1
######

dat_lateUP_Ripple1 <- dat_lateUP_Ripple %>%
  filter(
    if_all(all_of(z_cols), ~ abs(.) < z_thresh | is.na(.))
  )


mdl_final <- bam(nextUPV1_z ~ 
                   # 1. Surviving Main Power Effects
                   s(lastRippleHPC_z, k = 5) +
                   s(lastRippleV1_z, k = 5) +
                   
                   # 4. Control Term
                   s(AnimalID, bs = "re") +
                   s(SessionID, bs = "re"), 
                 
                 data = dat_lateUP_Ripple1 , 
                 # data = dat_clean , 
                 method = "fREML", 
                 discrete = TRUE, 
                 nthreads = 4)

#message("\n--- FINAL MODEL SUMMARY ---")
print(summary(mdl_final))


### HC bias
# Calculate scaling factors
raw_breaks <- c(-2,-1,0,1,2)
z_breaks <- sapply(raw_breaks, function(val) {
  dat_clean$lastRippleHPC_z[which.min(abs(dat_clean$lastRippleHPC - val))]
})

# cairo_pdf("lateUPHPC_nextUP", width = 4.3, height = 4.3)
p_HC_raw <- draw(mdl_final, select = "s(lastRippleHPC_z)", residuals = FALSE, rug = FALSE) + 
  theme_bw(base_family = "Arial") + 
  theme(aspect.ratio = 1) +
  scale_x_continuous(breaks = z_breaks, labels = raw_breaks) +
  labs(
    title = "HC lastRippleHPC and nextUP V1", 
    x = "lateUP HC bias (ripple in 100ms)", 
    y = "Partial Effect"
  )
dev.new(noRStudioGD = TRUE)
print(p_HC_raw)



### HC bias
# Calculate scaling factors
raw_breaks <- c(-2,-1,0,1,2)
z_breaks <- sapply(raw_breaks, function(val) {
  dat_clean$lastRippleHPC_z[which.min(abs(dat_clean$lastRippleHPC - val))]
})

# cairo_pdf("lateUPHPC_nextUP", width = 4.3, height = 4.3)
p_HC_raw <- draw(mdl_final, select = "s(lastRippleHPC_z)", residuals = FALSE, rug = FALSE) + 
  theme_bw(base_family = "Arial") + 
  theme(aspect.ratio = 1) +
  scale_x_continuous(breaks = z_breaks, labels = raw_breaks) +
  labs(
    title = "HC lastRippleHPC and nextUP V1", 
    x = "lateUP HC bias (ripple in 100ms)", 
    y = "Partial Effect"
  )
dev.new(noRStudioGD = TRUE)
print(p_HC_raw)


###### last ripple that did not happen during late UP did not predict nextUP V1
######
mdl_final <- bam(nextUPV1_z ~ 
                   # 1. Surviving Main Power Effects
                   s(lastRippleHPC_z, k = 5) +
                   s(lastRippleV1_z, k = 5) +

                   # 4. Control Term
                   s(AnimalID, bs = "re") +
                   s(SessionID, bs = "re"), 
                 
                 data = dat_lateUP_NoRipple, 
                 method = "fREML", 
                 discrete = TRUE, 
                 nthreads = 4)

message("\n--- FINAL MODEL SUMMARY ---")
print(summary(mdl_final))

### HC bias
# Calculate scaling factors
raw_breaks <- c(-2,-1,0,1,2)
z_breaks <- sapply(raw_breaks, function(val) {
  dat_clean$lastRippleHPC_z[which.min(abs(dat_clean$lastRippleHPC - val))]
})

# cairo_pdf("lateUPHPC_nextUP", width = 4.3, height = 4.3)
p_HC_raw <- draw(mdl_final, select = "s(lastRippleHPC_z)", residuals = FALSE, rug = FALSE) + 
  theme_bw(base_family = "Arial") + 
  theme(aspect.ratio = 1) +
  scale_x_continuous(breaks = z_breaks, labels = raw_breaks) +
  labs(
    title = "HC lastRippleHPC and nextUP V1", 
    x = "lateUP HC bias (without ripple in 100ms)", 
    y = "Partial Effect"
  )
dev.new(noRStudioGD = TRUE)
print(p_HC_raw)



######
###### lateUP HC (100ms) when ripple happened predicts nextUP V1
######

# --- 3. Clean Data ---
z_thresh <- 3.5  # include <99.9th centiles of data

z_cols <- c(
  "lateUPHPC_z","lateUPV1_z"
)
dat_clean1 <- dat_clean %>%
  filter(
    if_all(all_of(z_cols), ~ abs(.) < z_thresh | is.na(.))
  )
# 
# mdl_final <- bam(lateUPV1_z ~
#                    # 1. Surviving Main Power Effects
#                    s(lateUPHPC_z, k = 5) +
#                    # s(lateUPV1_z, k = 5) +
# 
#                    # 4. Control Term
#                    s(AnimalID, bs = "re") +
#                    s(SessionID, bs = "re"),
# 
#                  data = dat_Ripple1 ,
#                  method = "fREML",
#                  discrete = TRUE,
#                  nthreads = 4)
# 
# #message("\n--- FINAL MODEL SUMMARY ---")
# print(summary(mdl_final))


dat_NoRipple1 <- dat_NoRipple %>%
  filter(
    if_all(all_of(z_cols), ~ abs(.) < z_thresh | is.na(.))
  )


dat_lateUP_NoRipple1 <- dat_lateUP_NoRipple %>%
  filter(
    if_all(all_of(z_cols), ~ abs(.) < z_thresh | is.na(.))
  )



mdl_final <- bam(nextUPV1_z ~ 
                   # 1. Surviving Main Power Effects
                   s(lateUPHPC_z, k = 5) +
                   s(lateUPV1_z, k = 5) +
                   
                   # 4. Control Term
                   s(AnimalID, bs = "re") +
                   s(SessionID, bs = "re"), 
                 
                 data = dat_lateUP_NoRipple1 , 
                 method = "fREML", 
                 discrete = TRUE, 
                 nthreads = 4)

#message("\n--- FINAL MODEL SUMMARY ---")
print(summary(mdl_final))




### HC bias
# Calculate scaling factors
raw_breaks <- c(-2,-1,0,1,2)
z_breaks <- sapply(raw_breaks, function(val) {
  dat_clean$lateUPHPC_z[which.min(abs(dat_clean$lateUPHPC - val))]
})

# cairo_pdf("lateUPHPC_nextUP", width = 4.3, height = 4.3)
p_HC_raw <- draw(mdl_final, select = "s(lateUPHPC_z)", residuals = FALSE, rug = FALSE) + 
  theme_bw(base_family = "Arial") + 
  theme(aspect.ratio = 1) +
  scale_x_continuous(breaks = z_breaks, labels = raw_breaks) +
  labs(
    title = "HC late UP bias and nextUP V1", 
    x = "lateUP HC bias (without ripple in 100ms)", 
    y = "Partial Effect"
  )
dev.new(noRStudioGD = TRUE)
print(p_HC_raw)


dat_lateUP_Ripple1 <- dat_lateUP_Ripple %>%
  filter(
    if_all(all_of(z_cols), ~ abs(.) < z_thresh | is.na(.))
  )

dat_Ripple1 <- dat_Ripple %>%
  filter(
    if_all(all_of(z_cols), ~ abs(.) < z_thresh | is.na(.))
  )

mdl_final <- bam(nextUPV1_z ~ 
                   # 1. Surviving Main Power Effects
                   s(lateUPHPC_z, k = 5) +
                   s(lateUPV1_z, k = 5) +
                   
                   # 4. Control Term
                   s(AnimalID, bs = "re") +
                   s(SessionID, bs = "re"), 
                 
                 data = dat_Ripple1  , 
                 method = "fREML", 
                 discrete = TRUE, 
                 nthreads = 4)

#message("\n--- FINAL MODEL SUMMARY ---")
print(summary(mdl_final))


### HC bias
# Calculate scaling factors
raw_breaks <- c(-2,-1,0,1,2)
z_breaks <- sapply(raw_breaks, function(val) {
  dat_clean$lateUPHPC_z[which.min(abs(dat_clean$lateUPHPC - val))]
})

# cairo_pdf("lateUPHPC_nextUP", width = 4.3, height = 4.3)
p_HC_raw <- draw(mdl_final, select = "s(lateUPHPC_z)", residuals = FALSE, rug = FALSE) + 
  theme_bw(base_family = "Arial") + 
  theme(aspect.ratio = 1) +
  scale_x_continuous(breaks = z_breaks, labels = raw_breaks) +
  labs(
    title = "HC late UP bias and nextUP V1", 
    x = "lateUP HC bias (with ripple)", 
    y = "Partial Effect"
  )
dev.new(noRStudioGD = TRUE)
print(p_HC_raw)



#########
######### 
z_cols <- c(
  "nextDOWNlag_z","nextDOWNSOPower_z","geo_coherenceRippleLate_z"
)

# z_cols <- c(
#   "nextDOWNlag_z","nextDOWNSOPower_z","geo_coherence_z"
# )
dat_clean1 <- dat_clean %>%
  filter(
    if_all(all_of(z_cols), ~ abs(.) < z_thresh | is.na(.))
  )


mdl_final <- bam(geo_coherenceRippleLate_z ~ 
                   # 1. Surviving Main Power Effects
                   # s(nextDOWNDuration_z, k = 5) +
                   # s(nextDOWNSOPower_z, k = 5) +
                   # s(nextDOWNlag_z, k = 5) +
                   # s(lastRipplePower_z,k=5)+
                   # te(nextDOWNSOPower_z, nextDOWNlag_z, k = 5)+
                   
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









z_cols <- c(
  "nextDOWNlag_z","nextDOWNSOPower_z","geo_coherenceRippleNext_z"
)
# z_cols <- c(
#   "lastRipplePower_z","nextDOWNlag_z","nextDOWNSOPower_z","geo_coherenceRippleNext_z"
# )
# 
# z_cols <- c(
#   "nextDOWNlag_z","nextDOWNSOPower_z","geo_coherenceNext_z"
# )
dat_clean1 <- dat_clean %>%
  filter(
    if_all(all_of(z_cols), ~ abs(.) < z_thresh | is.na(.))
  )


mdl_final <- bam(geo_coherenceRippleNext_z ~ 
                   # 1. Surviving Main Power Effects
                   # s(nextDOWNDuration_z, k = 5) +
                   # s(nextDOWNSOPower_z, k = 5) +
                   # s(nextDOWNlag_z, k = 5) +
                   # s(lastRipplePower_z,k=5)+
                   # te(nextDOWNSOPower_z, nextDOWNlag_z, k = 5)+
                   
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





### HC bias
# Calculate scaling factors
raw_breaks <- c(0,0.05,0.1,0.15,0.2)
z_breaks <- sapply(raw_breaks, function(val) {
  dat_clean$nextDOWNlag_z[which.min(abs(dat_clean1$nextDOWNlag - val))]
})

# cairo_pdf("lateUPHPC_nextUP", width = 4.3, height = 4.3)
p_lag_raw <- draw(mdl_final, select = "s(nextDOWNlag_z)", residuals = FALSE, rug = FALSE) + 
  theme_bw(base_family = "Arial") + 
  theme(aspect.ratio = 1) +
  scale_x_continuous(breaks = z_breaks, labels = raw_breaks) +
  labs(
    title = "DOWN synchrony (bilateral lag)", 
    x = "Bilateral DOWN transition lag", 
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
                      4, length.out = 100)

raw_range_lag  <- seq(0, 
                      0.2, length.out = 100)

# 2. Build the prediction grid using the RAW scales
pred_grid_ti <- expand.grid(
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





# ==============================================================================
# --- 2D TE of next DOWN SO power and next DOWN bilateral lag ---
# ==============================================================================

# 1. Define sequence ranges using your actual RAW data minimums and maximums
raw_range_so   <- seq(1, 4, length.out = 100)
raw_range_lag  <- seq(0, 0.5, length.out = 100)

# 2. Build the prediction grid using the RAW scales
pred_grid_te <- expand.grid(
  nextDOWNSOPower    = raw_range_so,
  nextDOWNlag        = raw_range_lag,
  AnimalID           = dat_clean1$AnimalID[1],  # Constant random effect
  SessionID          = dat_clean1$SessionID[1]   # Constant random effect
)

# 3. Create the Z-scored columns that the model expects 
pred_grid_te$nextDOWNSOPower_z <- (pred_grid_te$nextDOWNSOPower - mean(dat_clean1$nextDOWNSOPower, na.rm=TRUE)) / sd(dat_clean1$nextDOWNSOPower, na.rm=TRUE)
pred_grid_te$nextDOWNlag_z     <- (pred_grid_te$nextDOWNlag     - mean(dat_clean1$nextDOWNlag, na.rm=TRUE))     / sd(dat_clean1$nextDOWNlag, na.rm=TRUE)

# 4. Extract specific term components matrix using the z-scores
term_preds <- predict(mdl_final, newdata = pred_grid_te, type = "terms")

# --- CRITICAL FIX HERE ---
# mgcv drops the spaces and the 'k=5' syntax in the term name matrix column!
pred_grid_te$Interaction_Effect <- term_preds[, "te(nextDOWNSOPower_z,nextDOWNlag_z)"]

# 5. Plot Total Tensor Surface using RAW scales for X and Y axes
p_interaction_raw_scale <- ggplot(pred_grid_te, aes(x = nextDOWNSOPower, y = nextDOWNlag, fill = Interaction_Effect)) +
  geom_tile() + 
  geom_contour(aes(z = Interaction_Effect), color = "black", alpha = 0.2) + 
  scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "firebrick", midpoint = 0, name = "Effect") +
  theme_minimal() +
  theme(aspect.ratio = 1) +
  labs(
    title = "Full Tensor Product Surface (te Term)",
    subtitle = "Combined Main Effects + Interaction (Raw Scale Mapping)",
    x = "nextDOWNSOPower (Raw)", 
    y = "nextDOWNlag (Raw)"
  )

print(p_interaction_raw_scale)




###########
###########
###########

library(dplyr)

# Create an explicit index column first
dat_indexed <- dat_clean %>% 
  mutate(current_idx = row_number())

# Map by checking the immediate adjacent rows
final_sequence_map <- dat_indexed %>%
  mutate(
    # Look at the lateUPV1 from the row directly BEFORE
    actual_prev_late = lag(lateUPV1, n = 1),
    # Look at the earlyUPV1 from the row directly AFTER
    actual_next_early = lead(earlyUPV1, n = 1),
    
    # Store the row numbers of those neighbors
    prev_idx = lag(current_idx, n = 1),
    next_idx = lead(current_idx, n = 1)
  ) %>%
  # Filter where the fingerprints match exactly
  filter(
    previousUPV1 == actual_prev_late & 
      nextUPV1 == actual_next_early
  ) %>%
  # Keep only the triplet index map and the biases for your models
  select(prev_idx, current_idx, next_idx, 
         previousUPV1, lateUPV1, earlyUPV1, nextUPV1)

# Check your new map
head(final_sequence_map)

# 
# library(mgcv)
# library(mgcViz)
# # Fit the model where HC bias is predicted by the interaction of Pre-V1 and Time
# mdl_next <- bam(nextUPV1_z ~ 
#                       # s(HC_post_z,k = 5) +
#                       # s(V1_pre_z, k = 5) +
#                       te(lateUPV1_z, lateUPHPC_z, k = c(5,5)) +
#                       s(AnimalID, bs = "re")+
#                       s(SessionID, bs = "re"), 
#                     data = dat_Ripple1, method = "fREML", discrete = TRUE)
# 
# print(summary(mdl_next))
# 
# 
# 
# # Convert to an mgcViz object and plot
# viz <- getViz(mdl_next)
# plot(sm(viz, 1)) + 
#   l_fitRaster() + 
#   l_fitContour() + 
#   scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
#   labs(title = "Next V1 Bias as a function of late V1 Bias and late HC Bias",
#        x = "Late V1 Track Bias", y = "Late HC Track bias")

