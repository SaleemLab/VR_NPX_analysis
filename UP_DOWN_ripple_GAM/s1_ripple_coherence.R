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
#dat <- read.csv("C:/Users/masah/Documents/GitHub/VR_NPX_analysis/UP_DOWN_ripple_GAM/pre_post_normalised_UP_ripple_20_120ms.csv")
#dat <- read.csv("C:/Users/masah/Documents/GitHub/VR_NPX_analysis/UP_DOWN_ripple_GAM/pre_post_normalised_UP_ripple_-150_-50ms.csv")
dat <- read.csv("C:/Users/masah/Documents/GitHub/VR_NPX_analysis/UP_DOWN_ripple_GAM/pre_post_normalised_UP_ripple.csv")

dat$SessionID <- as.factor(dat$SessionID)
dat$AnimalID <- as.factor(dat$AnimalID)


# 2. Create the Match/NonMatch variables
dat$V1_mean <- rowMeans(dat[, c("V1_pre", "V1_post")], na.rm = TRUE)

dat <- dat %>%
  mutate(
    # Determine the side V1 is biased toward
    TargetSide = ifelse(V1_post > 0, "R", "L"),
    
    # --- Spindle Presence Variables ---
    Spindle_Match    = ifelse(TargetSide == "R", Spindle_R, Spindle_L),
    Spindle_NonMatch = ifelse(TargetSide == "R", Spindle_L, Spindle_R),
    
    # --- Spindle Power Variables ---
    SpindlePower_Match    = ifelse(TargetSide == "R", SpindlePower_R, SpindlePower_L),
    SpindlePower_NonMatch = ifelse(TargetSide == "R", SpindlePower_L, SpindlePower_R),
    
    # --- Normalised UP/DOWN Variables --- Old
    NormalisedUP_Match    = ifelse(TargetSide == "R", NormalisedUP_R, NormalisedUP_L),
    NormalisedUP_NonMatch = ifelse(TargetSide == "R", NormalisedUP_L, NormalisedUP_R),
    
    NormalisedDOWN_Match    = ifelse(TargetSide == "R", NormalisedDOWN_R, NormalisedDOWN_L),
    NormalisedDOWN_NonMatch = ifelse(TargetSide == "R", NormalisedDOWN_L, NormalisedDOWN_R),
    
    # --- Time From/To UP Variables ---
    TimeFromUP_Match    = ifelse(TargetSide == "R", TimeFromUP_R, TimeFromUP_L),
    TimeFromUP_NonMatch = ifelse(TargetSide == "R", TimeFromUP_L, TimeFromUP_R),
    
    TimeToDOWN_Match    = ifelse(TargetSide == "R", TimeToDOWN_R, TimeToDOWN_L),
    TimeToDOWN_NonMatch = ifelse(TargetSide == "R", TimeToDOWN_L, TimeToDOWN_R),
    
    # --- Time From/To DOWN Variables ---
    TimeFromDOWN_Match    = ifelse(TargetSide == "R", TimeFromDOWN_R, TimeFromDOWN_L),
    TimeFromDOWN_NonMatch = ifelse(TargetSide == "R", TimeFromDOWN_L, TimeFromDOWN_R),
    
    TimeToUP_Match    = ifelse(TargetSide == "R", TimeToUP_R, TimeToUP_L),
    TimeToUP_NonMatch = ifelse(TargetSide == "R", TimeToUP_L, TimeToUP_R),
    
    # --- Duration Variables ---
    UPDuration_Match    = ifelse(TargetSide == "R", UPDuration_R, UPDuration_L),
    UPDuration_NonMatch = ifelse(TargetSide == "R", UPDuration_L, UPDuration_R),
    
    DOWNDuration_Match    = ifelse(TargetSide == "R", DOWNDuration_R, DOWNDuration_L),
    DOWNDuration_NonMatch = ifelse(TargetSide == "R", DOWNDuration_L, DOWNDuration_R)
  ) %>%
  # Remove the helper column
  select(-TargetSide)


# Get UP DOWN normalised duration combined 
library(dplyr)

# Double check what your match/non-match column is actually called!
# colnames(dat_clean) 

dat <- dat %>%
  mutate(
    # 1. Create the separate MATCH lifecycle column (0 to 2)
    NormalisedUPDOWN_Match = case_when(
      !is.na(NormalisedUP_Match)   ~ NormalisedUP_Match,      # 0 to 1 if UP exists
      !is.na(NormalisedDOWN_Match) ~ 1 + NormalisedDOWN_Match, # 1 to 2 if DOWN exists
      TRUE                         ~ NA_real_
    ),
    
    # 2. Create the separate NON-MATCH lifecycle column (0 to 2)
    NormalisedUPDOWN_NonMatch = case_when(
      !is.na(NormalisedUP_NonMatch)   ~ NormalisedUP_NonMatch, 
      !is.na(NormalisedDOWN_NonMatch) ~ 1 + NormalisedDOWN_NonMatch,
      TRUE                            ~ NA_real_
    )
  )

# Verify the transformations worked cleanly
message("Match lifecycle range: ", paste(range(dat$NormalisedUPDOWN_Match, na.rm = TRUE), collapse = " to "))
message("Non-Match lifecycle range: ", paste(range(dat$NormalisedUPDOWN_NonMatch, na.rm = TRUE), collapse = " to "))


# 3. Scale all numeric columns and add "_z" suffix
# This will now include the new Match/NonMatch variables you just created
dat$Spindle_Match <- as.factor(dat$Spindle_Match)
dat$Spindle_NonMatch <- as.factor(dat$Spindle_NonMatch)

dat$geo_coherence <- sign(dat$V1_post * dat$HC_post) * 
  sqrt(abs(dat$V1_post * dat$HC_post))

# 2. Calculate Signed Geometric Mean for V1-PRE and HPC
dat$geo_coherencePRE <- sign(dat$V1_pre * dat$HC_post) * 
  sqrt(abs(dat$V1_pre * dat$HC_post))


# Create a cleaned data frame with no NAs in these two variables
# dat <- dat[!is.na(dat$NormalisedUP_L) & !is.na(dat$NormalisedUP_R), ]
# Keep rows where NOT (both are NA)
#dat <- dat[!(is.na(dat$NormalisedUP_L) & is.na(dat$NormalisedUP_R)), ]


dat_Time <- dat %>%
  filter(
    # Keep if the value is <= 1s OR if the value is NaN/NA
    (UPDuration_Match <= 2  ) |
      (UPDuration_NonMatch <= 2 ) 
      #(UPDuration_Match >= 0.1    | is.na(UPDuration_Match)) &
      #(UPDuration_NonMatch >= 0.1 | is.na(UPDuration_NonMatch)) 
  )
dat_Time <- dat_Time %>%
  mutate(across(where(is.numeric), 
                ~ as.numeric(scale(.)), 
                .names = "{.col}_z"))



dat_clean <- dat %>%
  mutate(across(where(is.numeric), 
                ~ as.numeric(scale(.)), 
                .names = "{.col}_z"))

# Check the results
head(dat_clean)

# Verify the new columns exist
colnames(dat_clean)
# 
# message(sprintf("Trimming extreme outliers beyond +/- %s Z-scores for %d variables...", z_thresh, length(z_cols)))
# 
# --- 3. Clean Data ---
z_thresh <- 3.5  # include <99th centiles of data

# List of the new Z-score columns to filter
z_cols <- c(
  "V1_pre_z", "V1_post_z", "HC_post_z" , "HC_pre_z" ,
  "geo_coherencePRE_z","geo_coherence_z" 
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





########
########
########


library(mgcv)
library(mgcViz)
# Fit the model where HC bias is predicted by the interaction of Pre-V1 and Time
mdl_pre_post <- bam(V1_post_z ~ 
                    # s(HC_post_z,k = 5) +
                    # s(V1_pre_z, k = 5) +
                    te(V1_pre_z, HC_post_z, k = c(5,5)) +
                    s(AnimalID, bs = "re")+
                    s(SessionID, bs = "re"), 
                  data = dat_clean, method = "fREML", discrete = TRUE)

print(summary(mdl_pre_post))



# Convert to an mgcViz object and plot
viz <- getViz(mdl_pre_post)
plot(sm(viz, 1)) + 
  l_fitRaster() + 
  l_fitContour() + 
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  labs(title = "Post V1 Bias as a function of Pre-V1 Bias and HC",
       x = "Pre-V1 Track Bias", y = "HC Track bias")




library(mgcv)
library(ggplot2)
library(dplyr)

# ---------------------------------------------------------
# Plot 1: Main Effect of HC_post_z
# ---------------------------------------------------------
sm_hc <- smooth_estimates(mdl_visual, smooth = "s(HC_post_z)", n = 100) %>%
  mutate(
    .lower_ci = .estimate - (1.96 * .se),
    .upper_ci = .estimate + (1.96 * .se)
  )

p_hc <- ggplot(sm_hc, aes(x = HC_post_z, y = .estimate)) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), alpha = 0.2, fill = "dodgerblue") +
  geom_line(linewidth = 1, color = "dodgerblue") +
  theme_bw() + 
  labs(title = "Isolated Main Effect: HC Post (Z)", 
       x = "HC_post_z", y = "Partial Effect")
dev.new(noRStudioGD = TRUE)
# cairo_pdf("SO_NonMatch_effect.pdf", width = 4.3, height = 4.3)
print(p_hc)
# ---------------------------------------------------------
# Plot 2: Main Effect of V1_pre_z
# ---------------------------------------------------------
sm_v1 <- smooth_estimates(mdl_visual, smooth = "s(V1_pre_z)", n = 200) %>%
  mutate(
    .lower_ci = .estimate - (1.96 * .se),
    .upper_ci = .estimate + (1.96 * .se)
  )

p_v1 <- ggplot(sm_v1, aes(x = V1_pre_z, y = .estimate)) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), alpha = 0.2, fill = "firebrick") +
  geom_line(linewidth = 1, color = "firebrick") +
  theme_bw() + 
  labs(title = "Isolated Main Effect: V1 Pre (Z)", 
       x = "V1_pre_z", y = "Partial Effect")

print(p_v1)


# ---------------------------------------------------------
# Plot 3: Pure Interaction Tensor (ti) 
# ---------------------------------------------------------
sm_2d <- smooth_estimates(mdl_pre_post, smooth = "ti(V1_pre_z,HC_post_z)", n = 100)

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


# ==============================================================================
# --- RECREATING THE TOTAL LANDSCAPE ---
# ==============================================================================

# 1. Create a clean grid across the observed range of your data
# grid_range_v1 <- seq(min(dat_clean$V1_pre_z, na.rm=TRUE), max(dat_clean$V1_pre_z, na.rm=TRUE), length.out = 100)
# grid_range_hc <- seq(min(dat_clean$HC_post_z, na.rm=TRUE), max(dat_clean$HC_post_z, na.rm=TRUE), length.out = 100)
grid_range_v1 <- seq(-3,3, length.out = 100)
grid_range_hc <- seq(-3,3, length.out = 100)

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
  scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "firebrick", midpoint = 0, name = "Total\nEffect") +
  theme_minimal() +
  theme(aspect.ratio = 1) +
  labs(
    title = "Total Combined Effect Surface", 
    subtitle = "Reconstructed: s(HC) + s(V1) + ti(V1, HC)",
    x = "V1_pre_z", 
    y = "HC_post_z"
  )
print(p_combined)


########
########
########

# Coherence does not depend on UP duration
mdl_1 <- bam(geo_coherence_z ~ 
               s(UPDuration_Match_z, k = 5) +
               s(UPDuration_NonMatch_z, k = 5) +
               #ti(UPDuration_Match_z, UPDuration_NonMatch_z, k = 5) +
               
               # 4. Control Term
               s(AnimalID, bs = "re") +
               s(SessionID, bs = "re"), 
             
             data = dat_clean, 
             method = "fREML", 
             discrete = TRUE, 
             nthreads = 4)

message("\n--- FINAL MODEL SUMMARY ---")
print(summary(mdl_1))

########
########


########
########
########

# Coherence depends on Normalised UP duration and Time To DOWN
mdl_1 <- bam(geo_coherence_z ~ 
               s(NormalisedUP_Match_z, k = 5) +
               s(NormalisedUP_NonMatch_z, k = 5) +
               #ti(NormalisedUP_Match_z, NormalisedUP_NonMatch_z, k = 5) +
    
                   # 4. Control Term
                   s(AnimalID, bs = "re") +
                   s(SessionID, bs = "re"), 
                 
                 data = dat_clean, 
                 method = "fREML", 
                 discrete = TRUE, 
                 nthreads = 4)

message("\n--- FINAL MODEL SUMMARY ---")
print(summary(mdl_1))


### Normalised duration

# Calculate scaling factors
raw_breaks <- c(0, 0.25, 0.5,0.75,1)
z_breaks <- sapply(raw_breaks, function(val) {
  dat_clean$NormalisedUP_Match_z[which.min(abs(dat_clean$NormalisedUP_Match - val))]
})

# cairo_pdf("ripple_power_post_z_raw.pdf", width = 4.3, height = 4.3)
p_norm_raw <- draw(mdl_1, select = "s(NormalisedUP_Match_z)", residuals = FALSE, rug = FALSE) + 
  theme_bw(base_family = "Arial") + 
  theme(aspect.ratio = 1) +
  scale_x_continuous(breaks = z_breaks, labels = raw_breaks) +
  labs(
    title = "NormalisedUP_Match and Coherence", 
    x = "All Ripples Normalised UP", 
    y = "Partial Effect"
  )

print(p_norm_raw)



### Normalised duration Non matched

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



# Coherence does not depend on exact Timing relative To DOWN transition
mdl_2 <- bam(geo_coherence_z ~ 

                   s(TimeToDOWN_Match_z, k = 5) +
                   s(TimeToDOWN_NonMatch_z, k = 5) + 
               #ti(TimeToDOWN_Match_z, TimeToDOWN_NonMatch_z, k = 5) +
                   
                   # 4. Control Term
                   s(AnimalID, bs = "re") +
                   s(SessionID, bs = "re"), 
                 
                 data = dat_clean, 
                 method = "fREML", 
                 discrete = TRUE, 
                 nthreads = 4)

message("\n--- FINAL MODEL SUMMARY ---")
print(summary(mdl_2))



# Coherence does not depend on Time from UP
mdl_FromUP <- bam(geo_coherence_z ~ 
               
               s(TimeFromUP_Match_z, k = 5) +
               s(TimeFromUP_NonMatch_z, k = 5) + 
                 #ti(TimeFromUP_Match_z, TimeFromUP_NonMatch_z, k = 5) +
               
               # 4. Control Term
               s(AnimalID, bs = "re") +
               s(SessionID, bs = "re"), 
             
             data = dat_clean, 
             method = "fREML", 
             discrete = TRUE, 
             nthreads = 4)

message("\n--- FINAL MODEL SUMMARY ---")
print(summary(mdl_FromUP))




### Time from UP Match

# 1. Define raw time breaks (e.g., 0 to 500ms in 125ms steps)
raw_breaks_time <- c(0, 0.25, 0.5, 1,2)

# 2. Map raw seconds to the nearest Z-score in the cleaned dataset
z_breaks_match <- sapply(raw_breaks_time, function(val) {
  dat_Time$TimeFromUP_Match_z[which.min(abs(dat_Time$TimeFromUP_Match - val))]
})

# 3. Plot
p_time_match <- draw(mdl_FromUP, select = "s(TimeFromUP_Match_z)", residuals = FALSE, rug = FALSE) + 
  theme_bw(base_family = "Arial") + 
  theme(aspect.ratio = 1) +
  scale_x_continuous(breaks = z_breaks_match, labels = raw_breaks_time) +
  labs(
    title = "TimeFromUP_Match (Match) and Coherence", 
    x = "Time from UP transition (s)", 
    y = "Partial Effect"
  )

print(p_time_match)





# Coherence depends on Normalised UP duration and Time To DOWN

mdl_3 <- bam(geo_coherence_z ~ 
                   # 1. Surviving Main Power Effects
                   # s(lastRippleV1PRE_z, k = 5) +
                   # s(lastRipple_z, k = 5) +
                   
                   s(NormalisedUP_Match_z, k = 5) +
               #s(NormalisedUP_NonMatch_z, k = 5) +
                   # ti(NormalisedUP_Match_z, NormalisedUP_NonMatch_z, k = 5) +
                   s(TimeToDOWN_Match_z, k = 5) +
               ti(NormalisedUP_Match_z, TimeToDOWN_Match_z, k = 5) +
               #s(TimeToDOWN_NonMatch_z, k = 5) +
                   # ti(TimeToDOWN_Match_z, TimeToDOWN_NonMatch_z, k = 5) +
                   
               #s(TimeFromUP_Match_z, k = 5) +
               #   ti(NormalisedUP_Match_z, TimeFromUP_Match_z, k = 5) +
 
                   # 4. Control Term
                   s(AnimalID, bs = "re") +
                   s(SessionID, bs = "re"), 
                 
                 data = dat_clean, 
                 method = "fREML", 
                 discrete = TRUE, 
                 nthreads = 4)

message("\n--- FINAL MODEL SUMMARY ---")
print(summary(mdl_3))



######
###### 
######
time_threshold <- 1


mdl_final <- bam(geo_coherence_z ~ 
                   # 1. Surviving Main Power Effects
                   # s(lastRippleV1PRE_z, k = 5) +
                   # s(lastRipple_z, k = 5) +
                   
                   s(NormalisedUP_Match_z, k = 5) +
                   s(TimeToDOWN_Match_z, k = 5) +
                   te(NormalisedUP_Match_z, TimeToDOWN_Match_z, k = 5) +
  
                   # 4. Control Term
                   s(AnimalID, bs = "re") +
                   s(SessionID, bs = "re"), 
                 
                 data = dat_clean, 
                 method = "fREML", 
                 discrete = TRUE, 
                 nthreads = 4)

message("\n--- FINAL MODEL SUMMARY ---")
print(summary(mdl_final))



########
######## Normalised Match VS TimeToDOWN Match (Interaction)
########
time_threshold = 2

# 1. Define Parameters for Scaling
mean_up_m   <- mean(dat_Time$NormalisedUP_Match, na.rm = TRUE)
sd_up_m     <- sd(dat_Time$NormalisedUP_Match, na.rm = TRUE)

mean_down_m <- mean(dat_Time$TimeToDOWN_Match, na.rm = TRUE)
sd_down_m   <- sd(dat_Time$TimeToDOWN_Match, na.rm = TRUE)

# 2. Define Raw Ranges and Convert to Z
# Note: TimeToDOWN is capped at 1s per your filter
raw_up_seq   <- seq(0, 1, length.out = 100)
raw_down_seq <- seq(0, time_threshold, length.out = 100) 

z_up_seq   <- (raw_up_seq - mean_up_m) / sd_up_m
z_down_seq <- (raw_down_seq - mean_down_m) / sd_down_m

# 3. Create the prediction grid
# Ensure variables match the exact names used in mdl_final
pred_grid_inter <- expand.grid(
  NormalisedUP_Match_z = z_up_seq,
  TimeToDOWN_Match_z   = z_down_seq,
  AnimalID  = dat_Time$AnimalID[1],
  SessionID = dat_Time$SessionID[1]
)

# 4. Extract Interaction Term
# Update the string "ti(...)" to match exactly how it appears in summary(mdl_final)
term_preds_inter <- predict(mdl_final, newdata = pred_grid_inter, type = "terms")
term_name <- "ti(NormalisedUP_Match_z,TimeToDOWN_Match_z)"

pred_grid_inter$Interaction_Effect <- term_preds_inter[, term_name]

# 5. Convert Z back to Raw for plotting
pred_grid_inter$UP_Raw   <- (pred_grid_inter$NormalisedUP_Match_z * sd_up_m) + mean_up_m
pred_grid_inter$DOWN_Raw <- (pred_grid_inter$TimeToDOWN_Match_z * sd_down_m) + mean_down_m

########
######## Plotting
########

p_inter_up_down <- pred_grid_inter %>%
  ggplot(aes(x = UP_Raw, y = DOWN_Raw, fill = Interaction_Effect)) +
  geom_tile() + 
  geom_contour(aes(z = Interaction_Effect), color = "black", alpha = 0.2) + 
  scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "firebrick", 
                       midpoint = 0, name = "Effect") +
  scale_x_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, time_threshold), expand = c(0,0)) +
  theme_minimal(base_family = "Arial") +
  theme(aspect.ratio = 1) +
  labs(title = "Interaction: Normalised Position vs Time to DOWN", 
       subtitle = paste("Term:", term_name),
       x = "Normalised UP Match (Raw)", 
       y = "Time to DOWN Match (s)")


print(p_inter_up_down)

dev.new(noRStudioGD = TRUE); print(p_inter_up_down)






ggplot(dat_clean, aes(x = NormalisedDOWN_Match, y = NormalisedDOWN_NonMatch)) +
  geom_bin2d(bins = 20) +
  scale_fill_viridis_c() +
  theme_minimal()



ggplot(dat_clean, aes(x = NormalisedUP_Match, y = TimeToDOWN_Match)) +
  geom_bin2d(bins = 20) +
  scale_fill_viridis_c() +
  theme_minimal()




library(ggplot2)
library(viridis)

p_density <- dat_clean %>%
  ggplot(aes(x = NormalisedUP_Match, y = TimeToDOWN_Match)) +
  geom_bin2d(bins = 20) + # Adjust bins for higher/lower resolution
  scale_fill_viridis_c(option = "plasma", name = "Event Count") +
  scale_x_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 2), expand = c(0,0)) +
  theme_minimal(base_family = "Arial") +
  theme(aspect.ratio = 1) +
  labs(
    title = "Data Density: Normalised Position vs Time to DOWN",
    subtitle = "Where do the actual ripple events sit?",
    x = "Normalised UP Match (Raw)",
    y = "Time to DOWN Match (s)"
  )
print(p_density)



library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

# 1. Define configuration parameters
time_threshold <- 2.0
x_bins         <- 25    
time_step_s    <- 0.05  # 50ms intervals

# 2. Break the continuous data into discrete, raw grid bins
dat_binned <- dat_clean %>%
  filter(NormalisedUP_Match >= 0 & NormalisedUP_Match <= 1,
         TimeToDOWN_Match  >= 0 & TimeToDOWN_Match  <= time_threshold) %>%
  mutate(
    X_bin = round(NormalisedUP_Match * x_bins) / x_bins,
    Y_bin = round(TimeToDOWN_Match / time_step_s) * time_step_s
  ) %>%
  count(X_bin, Y_bin, name = "Event_Count")

# 3. Create a complete rectangular grid skeleton
full_grid <- expand.grid(
  X_bin = seq(0, 1, by = 1 / x_bins),
  Y_bin = seq(0, time_threshold, by = time_step_s)
)

# 4. Merge data into the full grid and replace missing slots with a true 0
dat_plot <- full_grid %>%
  left_join(dat_binned, by = c("X_bin", "Y_bin")) %>%
  mutate(Event_Count = replace_na(Event_Count, 0))

# 5. Dynamically calculate the 20% and 80% limits for the color bar
max_nonzero_count <- max(dat_plot$Event_Count, na.rm = TRUE)
lower_limit       <- 0
upper_limit       <- max_nonzero_count

# 6. Plot the raw binned histogram with REVERSED colors (White background, Black bins)
p_raw_bw <- dat_plot %>%
  ggplot(aes(x = X_bin, y = Y_bin, fill = Event_Count)) +
  geom_tile() +
  
  # REVERSED Gradient: Low counts are White, High counts are Black
  scale_fill_gradient(
    low = "white", 
    high = "black", 
    name = "Event\nCount",
    limits = c(lower_limit, upper_limit),
    oob = scales::squish # Low counts squish to white, high counts squish to black
  ) +
  
  # Format axes to perfectly match your data boundaries
  scale_x_continuous(limits = c(0, 1), expand = c(0,0), breaks = seq(0, 1, 0.25)) +
  scale_y_continuous(limits = c(0, time_threshold), expand = c(0,0), breaks = seq(0, time_threshold, 0.25)) +
  
  # Clean, classic white background theme
  theme_classic(base_family = "Arial") +
  theme(
    aspect.ratio = 1,
    panel.grid.major = element_line(color = "gray90", size = 0.2), # Subtle light gridlines
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(size = 12),
    axis.line = element_line(color = "black")
  ) +
  
  labs(
    title = "Raw Data Density: Normalised Position vs Time to DOWN",
    subtitle = paste0("White background (0 counts) | High density cued by dark/black bins"),
    x = "Normalised UP Match (Raw)",
    y = "Time to DOWN Match (s)"
  )

print(p_raw_bw)




















########
######## Normalised Match VS Nonmatch (main + interaction)
########

# Parameters for Match
mean_match <- mean(dat_clean$NormalisedUP_Match, na.rm = TRUE)
sd_match   <- sd(dat_clean$NormalisedUP_Match, na.rm = TRUE)

# Parameters for Non-Match
mean_nonmatch <- mean(dat_clean$NormalisedUP_NonMatch, na.rm = TRUE)
sd_nonmatch   <- sd(dat_clean$NormalisedUP_NonMatch, na.rm = TRUE)

# Define the Raw range (0 to 1) and convert to Z for the grid
raw_seq <- seq(0, 1, length.out = 100)
z_match_seq <- (raw_seq - mean_match) / sd_match
z_nonmatch_seq <- (raw_seq - mean_nonmatch) / sd_nonmatch

# 1. Create the prediction grid
pred_grid_up <- expand.grid(
  NormalisedUP_Match_z    = z_match_seq,
  NormalisedUP_NonMatch_z = z_nonmatch_seq,
  AnimalID  = dat_clean$AnimalID[1],
  SessionID = dat_clean$SessionID[1]
)

# 2. Extract terms
term_preds_up <- predict(mdl_final, newdata = pred_grid_up, type = "terms")

# 3. Reconstruct Effects
# Total Effect: Match Smooth + Non-Match Smooth + Interaction
pred_grid_up$Total_Effect <- 
  term_preds_up[, "s(NormalisedUP_Match_z)"] + 
  term_preds_up[, "s(NormalisedUP_NonMatch_z)"] + 
  term_preds_up[, "ti(NormalisedUP_Match_z,NormalisedUP_NonMatch_z)"]

# Interaction Only
pred_grid_up$Interaction_Only <- 
  term_preds_up[, "ti(NormalisedUP_Match_z,NormalisedUP_NonMatch_z)"]

# Convert Z back to Raw for plotting axes
pred_grid_up$Match_Raw    <- (pred_grid_up$NormalisedUP_Match_z * sd_match) + mean_match
pred_grid_up$NonMatch_Raw <- (pred_grid_up$NormalisedUP_NonMatch_z * sd_nonmatch) + mean_nonmatch


########
######## Normalised Match VS Nonmatch (interaction only)
########
p_up_total <- pred_grid_up %>%
  ggplot(aes(x = Match_Raw, y = NonMatch_Raw, fill = Total_Effect)) +
  geom_tile() + 
  geom_contour(aes(z = Total_Effect), color = "black", alpha = 0.2) + 
  scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "firebrick", 
                       midpoint = 0, name = "Effect") +
  scale_x_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  theme_minimal(base_family = "Arial") +
  theme(aspect.ratio = 1) +
  labs(title = "Total Effect: Match vs Non-Match Bias", 
       subtitle = "Reconstructed: s(Match) + s(NonMatch) + ti(Match, NonMatch)",
       x = "Normalised UP Match (Raw)", y = "Normalised UP Non-Match (Raw)")

dev.new(noRStudioGD = TRUE); print(p_up_total)



p_up_inter_only <- pred_grid_up %>%
  ggplot(aes(x = Match_Raw, y = NonMatch_Raw, fill = Interaction_Only)) +
  geom_tile() + 
  geom_contour(aes(z = Interaction_Only), color = "black", alpha = 0.2) + 
  scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "firebrick", 
                       midpoint = 0, name = "Effect") +
  scale_x_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  theme_minimal(base_family = "Arial") +
  theme(aspect.ratio = 1) +
  labs(title = "Interaction Only: Match vs Non-Match Bias", 
       subtitle = "Term: ti(NormalisedUP_Match_z, NormalisedUP_NonMatch_z)",
       x = "Normalised UP Match (Raw)", y = "Normalised UP Non-Match (Raw)")

dev.new(noRStudioGD = TRUE); print(p_up_inter_only)






# --- 0. Get Parameters for Raw Scale (Match and Non-Match) ---
# Assuming dat_clean contains the original columns used to create the _z versions
mean_match <- mean(dat_clean$NormalisedUP_Match, na.rm = TRUE)
sd_match   <- sd(dat_clean$NormalisedUP_Match, na.rm = TRUE)

mean_nonmatch <- mean(dat_clean$NormalisedUP_NonMatch, na.rm = TRUE)
sd_nonmatch   <- sd(dat_clean$NormalisedUP_NonMatch, na.rm = TRUE)

# --- 1. Create 2D Prediction Grid (Still in Z-space for the model) ---
pred_grid_raw <- expand.grid(
  NormalisedUP_Match_z    = seq(-2.5, 2.5, length.out = 100),
  NormalisedUP_NonMatch_z = seq(-2.5, 2.5, length.out = 100),
  AnimalID  = dat_clean$AnimalID[1],
  SessionID = dat_clean$SessionID[1]
)

# --- 2. Extract specific terms ---
term_preds <- predict(mdl_final, newdata = pred_grid_raw, type = "terms")

# --- 3. Reconstruct Effect and Convert Grid to Raw Scale ---
pred_grid_raw$Total_Effect <- 
  term_preds[, "s(NormalisedUP_Match_z)"] + 
  term_preds[, "s(NormalisedUP_NonMatch_z)"] + 
  term_preds[, "ti(NormalisedUP_Match_z,NormalisedUP_NonMatch_z)"]

# Transform the grid coordinates back to 0-1 range
pred_grid_raw$Match_Raw <- (pred_grid_raw$NormalisedUP_Match_z * sd_match) + mean_match
pred_grid_raw$NonMatch_Raw <- (pred_grid_raw$NormalisedUP_NonMatch_z * sd_nonmatch) + mean_nonmatch

# --- 4. Plotting on Raw Scale ---
p_raw_interaction <- pred_grid_raw %>%
  ggplot(aes(x = Match_Raw, 
             y = NonMatch_Raw, 
             fill = Total_Effect)) +
  geom_tile() + 
  geom_contour(aes(z = Total_Effect), color = "black", alpha = 0.2) + 
  
  # Ensure the axes strictly respect the 0 to 1 range
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  
  scale_fill_gradient2(
    low = "dodgerblue", 
    mid = "white", 
    high = "firebrick", 
    midpoint = 0, 
    name = "Effect Size"
  ) +
  theme_minimal() +
  theme(aspect.ratio = 1) +
  labs(
    title = "Interaction on Raw Scale (0-1)",
    x = "Normalised UP Match (Raw)",
    y = "Normalised UP Non-Match (Raw)"
  )

print(p_raw_interaction)









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

print(p_ripple_raw)
dev.off()




# --- 0. Parameters for TimetoLastRipple ---
mean_time_to <- mean(dat_clean$TimeToLastRipple, na.rm = TRUE)
sd_time_to   <- sd(dat_clean$TimeToLastRipple, na.rm = TRUE)

# Define range based on your data distribution (example: -2 to 2 seconds)
raw_min_t <- -2
raw_max_t <- 2

z_min_t <- (raw_min_t - mean_time_to) / sd_time_to
z_max_t <- (raw_max_t - mean_time_to) / sd_time_to

# --- 1. Prediction Grid ---
pred_grid_time <- expand.grid(
  TimetoLastRipple_z = seq(z_min_t, z_max_t, length.out = 100),
  # Hold others at mean/reference
  lastRipple_z = 0,
  TimefromLastRipple_z = 0,
  SessionID = dat_clean$SessionID[1],
  AnimalID = dat_clean$AnimalID[1]
)

# --- 2. Extract and Plot ---
term_preds_time <- predict(mdl_final, newdata = pred_grid_time, type = "terms")
pred_grid_time$Effect <- term_preds_time[, "s(TimetoLastRipple_z)"]

ggplot(pred_grid_time, aes(x = TimetoLastRipple_z * sd_time_to + mean_time_to, y = Effect)) +
  geom_line(color = "firebrick", linewidth = 1) +
  theme_minimal() +
  labs(title = "Main Effect: Time to Last Ripple", x = "Time (s)", y = "Partial Effect")








# 1. Generate the smooth estimates for non-cyclic variables
message("Rendering 2D Joint Landscape (Linear)...")

# Note: Ensure 'smooth' matches the label in your mdl_final summary
sm_2d <- smooth_estimates(mdl_final, 
                          smooth = "te(lastRipple_z,lastRippleNormalisedUP_z)", 
                          n = 100) %>%
  mutate(
    .lower_ci = .estimate - (1.96 * .se),
    .upper_ci = .estimate + (1.96 * .se)
  )

# 2. Plotting without circular logic
p_joint <- sm_2d %>%
  ggplot(aes(x = lastRipple_z, y = lastRippleNormalisedUP_z, fill = .estimate)) +
  # Using geom_raster for better performance with 100x100 grids
  geom_raster(interpolate = TRUE) + 
  geom_contour(aes(z = .estimate), color = "black", alpha = 0.2, size = 0.3) + 
  # Gradient scale centered at 0
  scale_fill_gradient2(
    low = "dodgerblue", 
    mid = "white", 
    high = "firebrick", 
    midpoint = 0, 
    name = "Estimate"
  ) +
  # Standard linear axes
  labs(
    title = "Joint Effect of Ripple Metrics",
    subtitle = "Tensor Product Smooth Estimate",
    x = "Last Ripple (Z-score)", 
    y = "Normalised UP Ripple (Z-score)"
  ) +
  theme_minimal(base_family = "Arial") +
  theme(
    aspect.ratio = 1,
    legend.position = "right",
    panel.grid = element_blank() # Clean look for heatmaps
  )

# 3. Output
dev.new(noRStudioGD = TRUE)
print(p_joint)

# Increased width slightly for the legend, but kept aspect ratio square
cairo_pdf("SO_phase_post_z.pdf", width = 5.2, height = 4.3)
print(p_joint)
dev.off()


# ==============================================================================
# Total SPINDLE EFFECT (interaction + main)
# ==============================================================================

# 1. Create a dense prediction grid 
pred_grid_spindle <- expand.grid(
  SpindlePower_Match_Z = seq(-2.5, 2.5, length.out = 100),
  SOPhase_Match = seq(0, 2*pi, length.out = 100),
  
  # Hold other variables at baseline
  SpindlePower_NonMatch_Z = 0,
  RipplePower_Z = 0,
  SOPhase_NonMatch = pi,
  SessionID = dat_clean$SessionID[1],
  AnimalID = dat_clean$AnimalID[1] 
)

# 2. Extract and Sum terms
term_preds_spindle <- predict(mdl_final, newdata = pred_grid_spindle, type = "terms")
pred_grid_spindle$Reconstructed_TE <- 
  term_preds_spindle[, "s(SpindlePower_Match_Z)"] + 
  term_preds_spindle[, "ti(SpindlePower_Match_Z,SOPhase_Match)"]

# 3. Plot using the Phase Plot style
p_spindle_reconstructed <- pred_grid_spindle %>%
  arrange(SOPhase_Match, SpindlePower_Match_Z) %>%
  ggplot(aes(x = SOPhase_Match, y = SpindlePower_Match_Z, fill = Reconstructed_TE)) +
  
  # geom_tile logic from the phase plot
  geom_tile(width = (2*pi)/99, height = 5/99) + 
  geom_contour(aes(z = Reconstructed_TE), color = "black", alpha = 0.2) + 
  
  scale_fill_gradient2(
    low = "dodgerblue", 
    mid = "white", 
    high = "firebrick", 
    midpoint = 0, 
    name = "Effect"
  ) +
  
  # Visible grey ticks and lines
  scale_x_continuous(breaks = c(0, pi, 2*pi), labels = c("0", expression(pi), expression(2*pi)), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  
  theme_minimal(base_family = "Arial") +
  theme(
    aspect.ratio = 1,
    axis.ticks = element_line(color = "grey80"),
    axis.ticks.length = unit(0.2, "cm"),
    axis.line = element_line(color = "grey80"),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "Total Effect of Spindle Power by SO Phase", 
    subtitle = "Reconstructed: s(Spindle) + ti(Spindle, Phase)",
    x = "Match SO Phase (Rad)", 
    y = "Match Spindle Power (Z)"
  )

print(p_spindle_reconstructed)


# ==============================================================================
# Total RIPPLE EFFECT raw (main + interaction)
# ==============================================================================

# 0. Get conversion parameters for Raw Y-axis
mean_rip <- mean(dat_clean$RipplePower, na.rm = TRUE)
sd_rip   <- sd(dat_clean$RipplePower, na.rm = TRUE)

raw_min <- 5
raw_max <- 16

z_min <- (raw_min - mean_rip) / sd_rip
z_max <- (raw_max - mean_rip) / sd_rip

# 1. Create a clean 2D grid
pred_grid_ripple <- expand.grid(
  RipplePower_Z = seq(z_min, z_max, length.out = 100),
  SOPhase_Match = seq(0, 2*pi, length.out = 100),
  
  # Hold others constant
  SpindlePower_Match_Z = 0,
  SpindlePower_NonMatch_Z = 0,
  SOPhase_NonMatch = pi,
  SessionID = dat_clean$SessionID[1],
  AnimalID = dat_clean$AnimalID[1] 
)

# 2. Extract terms 
term_preds_ripple <- predict(mdl_final, newdata = pred_grid_ripple, type = "terms")

# 3. Sum the specific terms
# Matching the "Reconstructed" logic: Main Effect + Interaction
pred_grid_ripple$Reconstructed_TE <- 
  term_preds_ripple[, "s(RipplePower_Z)"] + 
  term_preds_ripple[, "ti(RipplePower_Z,SOPhase_Match)"]

# 4. Prepare Y-axis breaks
raw_axis_labels <- seq(5, 15, by = 5) 
z_axis_breaks   <- (raw_axis_labels - mean_rip) / sd_rip

# 5. Plotting using the exact Phase Plot style
p_ripple_reconstructed <- pred_grid_ripple %>%
  arrange(SOPhase_Match, RipplePower_Z) %>% # Maintaining the sort logic
  ggplot(aes(x = SOPhase_Match, y = RipplePower_Z, fill = Reconstructed_TE)) +
  
  # Using geom_tile with calculated widths for that specific "grid" look
  geom_tile(width = (2*pi)/99, height = (z_max - z_min)/99) + 
  geom_contour(aes(z = Reconstructed_TE), color = "black", alpha = 0.2) + 
  
  scale_fill_gradient2(
    low = "dodgerblue", 
    mid = "white", 
    high = "firebrick", 
    midpoint = 0, 
    name = "Effect"
  ) +
  
  # Match the axis ticks from the Phase plot
  scale_x_continuous(breaks = c(0, pi, 2*pi), 
                     labels = c("0", expression(pi), expression(2*pi)),
                     expand = c(0,0)) +
  
  # Raw units for Y-axis
  scale_y_continuous(breaks = z_axis_breaks, 
                     labels = raw_axis_labels, 
                     expand = c(0,0)) +
  
  theme_minimal(base_family = "Arial") +
  theme(
    aspect.ratio = 1,
    axis.ticks = element_line(color = "grey80"), # Visible ticks
    axis.ticks.length = unit(0.2, "cm"),
    axis.line = element_line(color = "grey80"),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "Total Effect of Ripple Power by SO Phase", 
    subtitle = "Reconstructed: s(Ripple) + ti(Ripple, Phase)",
    x = "Match SO Phase (Rad)", 
    y = "Ripple Power (Raw Units)"
  )

dev.new(noRStudioGD = TRUE)
print(p_ripple_reconstructed)


# ==============================================================================
# RIPPLE INTERACTION ONLY (ti Gate) raw
# ==============================================================================

# 0. Parameters for Raw Y-axis conversion (matches your previous code)
mean_rip <- mean(dat_clean$RipplePower, na.rm = TRUE)
sd_rip   <- sd(dat_clean$RipplePower, na.rm = TRUE)
raw_min <- 5
raw_max <- 16
z_min <- (raw_min - mean_rip) / sd_rip
z_max <- (raw_max - mean_rip) / sd_rip

# 1. Create the prediction grid
pred_grid_ripple_inter <- expand.grid(
  RipplePower_Z = seq(z_min, z_max, length.out = 100),
  SOPhase_Match = seq(0, 2*pi, length.out = 100),
  SpindlePower_Match_Z = 0,
  SpindlePower_NonMatch_Z = 0,
  SOPhase_NonMatch = pi,
  SessionID = dat_clean$SessionID[1],
  AnimalID = dat_clean$AnimalID[1]
)

# 2. Extract terms
term_preds_ripple <- predict(mdl_final, newdata = pred_grid_ripple_inter, type = "terms")

# 3. ISOLATE INTERACTION ONLY
# Note: We do NOT add the s(RipplePower_Z) term here
pred_grid_ripple_inter$Interaction_Only <- 
  term_preds_ripple[, "ti(RipplePower_Z,SOPhase_Match)"]

# 4. Prepare axis breaks
raw_axis_labels <- seq(5, 15, by = 5)
z_axis_breaks   <- (raw_axis_labels - mean_rip) / sd_rip

# 5. Plotting using the Phase Plot style
p_ripple_interaction <- pred_grid_ripple_inter %>%
  arrange(SOPhase_Match, RipplePower_Z) %>%
  ggplot(aes(x = SOPhase_Match, y = RipplePower_Z, fill = Interaction_Only)) +
  
  # geom_tile with the same logic as the phase Synergy plot
  geom_tile(width = (2*pi)/99, height = (z_max - z_min)/99) + 
  geom_contour(aes(z = Interaction_Only), color = "black", alpha = 0.2) + 
  
  scale_fill_gradient2(
    low = "dodgerblue", 
    mid = "white", 
    high = "firebrick", 
    midpoint = 0, 
    name = "Effect"
  ) +
  
  # Visible grey ticks and lines
  scale_x_continuous(breaks = c(0, pi, 2*pi), 
                     labels = c("0", expression(pi), expression(2*pi)),
                     expand = c(0,0)) +
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
  labs(
    title = "Interaction Effect of Ripple Power by SO Phase", 
    subtitle = "Term: ti(Ripple, Phase) Only (No Main Effect)",
    x = "Match SO Phase (Rad)", 
    y = "Ripple Power (Raw Units)"
  )

dev.new(noRStudioGD = TRUE)
print(p_ripple_interaction)

# ==============================================================================
# --- 3. MULTI-METRIC EFFECT SIZE CALCULATIONS ---
# ==============================================================================
message("\nCalculating 4 Effect Size Metrics (This may take a minute)...")

# EXACT formula components to rebuild the models robustly
formula_terms <- c(
  "s(RipplePower_Z, k = 5)", 
  "s(SpindlePower_Match_Z, k = 5)", 
  "s(SpindlePower_NonMatch_Z, k = 5)", 
  "ti(RipplePower_Z, SOPhase_Match, bs = c('tp', 'cc'), k = c(5, 8))", 
  "ti(SpindlePower_Match_Z, SOPhase_Match, bs = c('tp', 'cc'), k = c(5, 8))", 
  "te(SOPhase_Match, SOPhase_NonMatch, bs = 'cc', k = 8)"
)

# EXACT labels output by summary() and smooth_estimates()
smooth_labels <- c(
  "s(RipplePower_Z)", 
  "s(SpindlePower_Match_Z)", 
  "s(SpindlePower_NonMatch_Z)", 
  "ti(RipplePower_Z, SOPhase_Match", 
  "ti(SpindlePower_Match_Z, SOPhase_Match", 
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
# Plot 1: Main Effect of Match Phase
# ---------------------------------------------------------
sm_match <- smooth_estimates(mdl_final, smooth = "s(SOPhase_Match)", n = 200) %>%
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
sm_nonmatch <- smooth_estimates(mdl_final, smooth = "s(SOPhase_NonMatch)", n = 200) %>%
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
# Plot 3: The "Pure" Phase Interaction (ti term)
# ---------------------------------------------------------

# Note: Ensure your model (mdl_final) uses these renamed terms
sm_2d <- smooth_estimates(mdl_final, smooth = "ti(SOPhase_Match,SOPhase_NonMatch)", n = 100) %>%
  mutate(
    .lower_ci = .estimate - (1.96 * .se),
    .upper_ci = .estimate + (1.96 * .se)
  )

p_synergy <- sm_2d %>%
  mutate(
    SOPhase_Match = ifelse(SOPhase_Match < 0, SOPhase_Match + 2*pi, SOPhase_Match),
    SOPhase_NonMatch = ifelse(SOPhase_NonMatch < 0, SOPhase_NonMatch + 2*pi, SOPhase_NonMatch)
  ) %>%
  arrange(SOPhase_Match, SOPhase_NonMatch) %>%
  ggplot(aes(x = SOPhase_Match, y = SOPhase_NonMatch, fill = .estimate)) +
  geom_tile(width = (2*pi)/99, height = (2*pi)/99) + 
  geom_contour(aes(z = .estimate), color = "black", alpha = 0.2) + 
  scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "firebrick", midpoint = 0, name = "Effect") +
  scale_x_continuous(breaks = c(0, pi, 2*pi), labels = c("0", expression(pi), expression(2*pi))) +
  scale_y_continuous(breaks = c(0, pi, 2*pi), labels = c("0", expression(pi), expression(2*pi))) +
  theme_minimal(base_family = "Arial") +
  theme(aspect.ratio = 1) +
  labs(title = "Inter-hemispheric Phase interaction effect on coherence",
       x = "Match SO Phase (Rad)", y = "Non-Match SO Phase (Rad)")

dev.new(noRStudioGD = TRUE)
print(p_synergy)

# ==============================================================================
# --- RECREATING THE TOTAL LANDSCAPE ---
# ==============================================================================

# 1. Create a clean 2D grid
pred_grid_phase <- expand.grid(
  SOPhase_Match = seq(0, 2*pi, length.out = 100),
  SOPhase_NonMatch = seq(0, 2*pi, length.out = 100),
  RipplePower_Z = 0,
  SpindlePower_Match_Z = 0,
  SpindlePower_NonMatch_Z = 0,
  # Fixed the ID assignment
  SessionID = dat_clean$SessionID[1],
  AnimalID = dat_clean$AnimalID[1] 
)

# 2. Extract terms 
term_preds_phase <- predict(mdl_final, newdata = pred_grid_phase, type = "terms")

# 3. Sum the specific terms
# Note: Ensure these strings match your mdl_final formula exactly
pred_grid_phase$Reconstructed_TE <- 
  term_preds_phase[, "s(SOPhase_Match)"] + 
  term_preds_phase[, "s(SOPhase_NonMatch)"] + 
  term_preds_phase[, "ti(SOPhase_Match,SOPhase_NonMatch)"]

message("\nReconstructing Total Phase Synergy from Decomposed Terms...")

p_reconstructed <- pred_grid_phase %>%
  mutate(
    SOPhase_Match = ifelse(SOPhase_Match < 0, SOPhase_Match + 2*pi, SOPhase_Match),
    SOPhase_NonMatch = ifelse(SOPhase_NonMatch < 0, SOPhase_NonMatch + 2*pi, SOPhase_NonMatch)
  ) %>%
  arrange(SOPhase_Match, SOPhase_NonMatch) %>%
  ggplot(aes(x = SOPhase_Match, y = SOPhase_NonMatch, fill = Reconstructed_TE)) +
  geom_tile(width = (2*pi)/99, height = (2*pi)/99) + 
  geom_contour(aes(z = Reconstructed_TE), color = "black", alpha = 0.2) + 
  scale_fill_gradient2(
    low = "dodgerblue", 
    mid = "white", 
    high = "firebrick", 
    midpoint = 0, 
    name = "Effect"
  ) +
  scale_x_continuous(breaks = c(0, pi, 2*pi), labels = c("0", expression(pi), expression(2*pi))) +
  scale_y_continuous(breaks = c(0, pi, 2*pi), labels = c("0", expression(pi), expression(2*pi))) +
  theme_minimal(base_family = "Arial") +
  theme(aspect.ratio = 1) +
  labs(
    title = "Inter-hemispheric Phase effect on coherence", 
    subtitle = "Reconstructed: s(Match) + s(NonMatch) + ti(Match, NonMatch)",
    x = "Match SO Phase (Rad)", 
    y = "Non-Match SO Phase (Rad)"
  )

dev.new(noRStudioGD = TRUE)
print(p_reconstructed)