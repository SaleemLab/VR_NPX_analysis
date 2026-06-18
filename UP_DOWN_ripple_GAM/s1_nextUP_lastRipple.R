
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
  #"lastRippleV1_z", "lastRippleV1PRE_z", 
  #"lastRippleHPC_z", "lastRippleHPCPRE_z",
  "nextUPV1_z", "earlyUPV1_z", "lateUPV1_z", "previousUPV1_z",
  #"firstRippleV1_z", "firstRippleV1PRE_z",
  "nextUPHPC_z", "earlyUPHPC_z", "lateUPHPC_z", "previousUPHPC_z",
  #"firstRippleHPC_z", "firstRippleHPCPRE_z"
)

message(sprintf("Trimming extreme outliers beyond +/- %s Z-scores for %d variables...", z_thresh, length(z_cols)))

dat_clean <- dat_clean %>%
  filter(
    if_all(all_of(z_cols), ~ abs(.) < z_thresh | is.na(.))
  )

dat_Time <- dat_Time %>%
  filter(
    if_all(all_of(z_cols), ~ abs(.) < z_thresh | is.na(.))
  )

dat_lateUP_NoRipple <- dat_lateUP_NoRipple %>%
  filter(
    if_all(all_of(z_cols), ~ abs(.) < z_thresh | is.na(.))
  )

dat_NoRipple <- dat_NoRipple %>%
  filter(
    if_all(all_of(z_cols), ~ abs(.) < z_thresh | is.na(.))
  )

dat_Ripple <- dat_Ripple %>%
  filter(
    if_all(all_of(z_cols), ~ abs(.) < z_thresh | is.na(.))
  )


dat_Ripples <- dat_Ripples %>%
  filter(
    if_all(all_of(z_cols), ~ abs(.) < z_thresh | is.na(.))
  )

###### last ripple predicts nextUP V1
######

mdl_final <- bam(nextUPV1_z ~ 
                   # 1. Surviving Main Power Effects
                   s(lastRippleHPC_z, k = 5) +
                   s(lastRippleV1_z, k = 5) +
                   
                   # 4. Control Term
                   s(AnimalID, bs = "re") +
                   s(SessionID, bs = "re"), 
                 
                 data = dat_clean , 
                 method = "fREML", 
                 discrete = TRUE, 
                 nthreads = 4)

#message("\n--- FINAL MODEL SUMMARY ---")
print(summary(mdl_final))


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


######
######

library(mgcv)
library(mgcViz)
# Fit the model where HC bias is predicted by the interaction of Pre-V1 and Time
mdl_next <- bam(nextUPV1_z ~ 
                      # s(HC_post_z,k = 5) +
                      # s(V1_pre_z, k = 5) +
                      te(lateUPV1_z, lateUPHPC_z, k = c(5,5)) +
                      s(AnimalID, bs = "re")+
                      s(SessionID, bs = "re"), 
                    data = dat_clean, method = "fREML", discrete = TRUE)
print(summary(mdl_next))



# Convert to an mgcViz object and plot
viz <- getViz(mdl_next)
plot(sm(viz, 1)) + 
  l_fitRaster() + 
  l_fitContour() + 
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  labs(title = "Post V1 Bias as a function of Pre-V1 Bias and HC",
       x = "Pre-V1 Track Bias", y = "HC Track bias")



######
###### lateUP HC predicts nextUP V1
######

# --- 3. Clean Data ---
z_thresh <- 3.5  # include <99.9th centiles of data

z_cols <- c(
  # "lastRippleV1_z", "lastRippleV1PRE_z", 
  # "lastRippleHPC_z", "lastRippleHPCPRE_z",
  "nextUPV1_z", "earlyUPV1_z", "lateUPV1_z", "previousUPV1_z",
  # "firstRippleV1_z", "firstRippleV1PRE_z",
  "nextUPHPC_z", "earlyUPHPC_z", "lateUPHPC_z", "previousUPHPC_z"
  # "firstRippleHPC_z", "firstRippleHPCPRE_z"
)
dat_clean1 <- dat_clean %>%
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
                 
                 data = dat_clean1 , 
                 method = "fREML", 
                 discrete = TRUE, 
                 nthreads = 4)

#message("\n--- FINAL MODEL SUMMARY ---")
print(summary(mdl_final))


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
                 
                 data = dat_lateUP_NoRipple1  , 
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



######
###### lateUP nextUP coherence
######

mdl_final <- bam(geo_coherenceNext_z ~ 
                   # 1. Surviving Main Power Effects
                   # s(UPDuration_z, k = 5) +
                   #s(lastRippleNormalisedUP_z, k = 5) +
                   s(TimefromLastRipple_z, k = 5) +
                   
                   # 4. Control Term
                   s(AnimalID, bs = "re") +
                   s(SessionID, bs = "re"), 
                 
                 data = dat_clean  , 
                 method = "fREML", 
                 discrete = TRUE, 
                 nthreads = 4)

#message("\n--- FINAL MODEL SUMMARY ---")
print(summary(mdl_final))



### Normalised duration UP

# Calculate scaling factors
raw_breaks <- c(0, 0.25, 0.5,0.75,1)
z_breaks <- sapply(raw_breaks, function(val) {
  dat_clean$NormalisedUP_NonMatch_z[which.min(abs(dat_clean$NormalisedUP_NonMatch - val))]
})

# cairo_pdf("ripple_power_post_z_raw.pdf", width = 4.3, height = 4.3)
p_norm_nonmatch <- draw(mdl_1, select = "s(NormalisedUP_NonMatch_z)", residuals = FALSE, rug = FALSE) +
  theme_bw(base_family = "Arial") +
  theme(aspect.ratio = 1) +
  scale_x_continuous(breaks = z_breaks, labels = raw_breaks) +
  labs(
    title = "NormalisedUP_NonMatch and Coherence",
    x = "lastRippleNormalisedUP",
    y = "Partial Effect"
  )

print(p_norm_nonmatch)
