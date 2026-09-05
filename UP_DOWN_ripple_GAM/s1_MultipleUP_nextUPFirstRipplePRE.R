
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
  "firstRippleV1_z", "firstRippleV1PRE_z",
  "firstRippleHPC_z", "firstRippleHPCPRE_z"
)


my_folder <- "C:/Users/masah/Documents/GitHub/VR_NPX_analysis/UP_DOWN_ripple_GAM/MultipleUP_nextUP"

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


###########
########### Previous UP -> Current UP -> Next UP (Triplet)
###########
z_cols <- c(
  # "next_firstRippleV1PRE","prev_lastRippleHPC","curr_firstRippleV1PRE",
  "earlyUPV1","lastRippleHPC",'nextUPV1'
  # "next_firstRippleV1","curr_firstRippleV1"
)

dat_clean1 <- dat_clean %>%
  filter(
    if_all(all_of(z_cols), ~ abs(.) < z_thresh | is.na(.))
  )

library(dplyr)

# Create an explicit index column first
dat_indexed <- dat_clean1 %>% 
  mutate(current_idx = row_number())

# Map by checking the immediate adjacent rows
UP_chain_id <- dat_indexed %>%
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
  select(prev_idx, current_idx, next_idx)

# Check your new map
head(UP_chain_id)
# UP_chain_id <- UP_chain_id %>%
#   filter(
#     prev_idx %in% current_idx,
#     next_idx %in% current_idx
#   )


library(dplyr)

# Helper function keeping your EXACT column names with a prefix
pull_ripple_data <- function(indices, prefix) {
  dat_indexed %>%
    filter(current_idx %in% indices) %>%
    select(
      current_idx,
      firstRippleV1_z, firstRippleV1PRE_z, lastRippleV1_z, lastRippleV1PRE_z,
      nextUPV1_z, earlyUPV1_z, lateUPV1_z, previousUPV1_z,
      firstRippleHPC_z, firstRippleHPCPRE_z, lastRippleHPC_z, lastRippleHPCPRE_z,
      nextUPHPC_z, earlyUPHPC_z, lateUPHPC_z, previousUPHPC_z,
      firstRippleV1, firstRippleV1PRE, lastRippleV1, lastRippleV1PRE,
      nextUPV1, earlyUPV1, lateUPV1, previousUPV1,
      firstRippleHPC, firstRippleHPCPRE, lastRippleHPC, lastRippleHPCPRE,
      nextUPHPC, earlyUPHPC, lateUPHPC, previousUPHPC,is_near_DOWN,ripple_quantile,
      TimefromFirstRipple,TimefromLastRipple,TimeToFirstRipple,TimeToLastRipple,
      
      AnimalID,SessionID
    ) %>%
    rename_with(~ paste0(prefix, "_", .), -c(current_idx, AnimalID, SessionID))
  
}

# Extract data for each step in the chain
prev_ripples    <- pull_ripple_data(UP_chain_id$prev_idx, "prev")
curr_ripples    <- pull_ripple_data(UP_chain_id$current_idx, "curr")
next_ripples    <- pull_ripple_data(UP_chain_id$next_idx, "next")

# Join them into your final modeling dataframe
ripple_chain_dataset <- UP_chain_id %>%
  left_join(prev_ripples,    by = c("prev_idx" = "current_idx")) %>%
  left_join(curr_ripples,    by = c("current_idx" = "current_idx")) %>%
  left_join(next_ripples,    by = c("next_idx" = "current_idx"))

# Quick look at your new model-ready variables
colnames(ripple_chain_dataset)


# ###### Previous
# mdl_rippleEcho <- bam(next_firstRippleV1PRE ~
#                        s(prev_lastRippleHPC, k = 5) +
#                        s(AnimalID, bs = "re")+
#                        s(SessionID, bs = "re"),
#                      data = ripple_chain_dataset, method = "fREML", discrete = TRUE)
# print(summary(mdl_rippleEcho))


###### Previous
mdl_multipleUP <- bam(next_firstRippleV1_z ~
                        # prev_is_near_DOWN +
                        # s(prev_lastRippleHPC_z,by = prev_is_near_DOWN, k = 5) +
                        # prev_ripple_quantile +
                        # s(prev_lastRippleHPC_z,by = prev_ripple_quantile, k = 5) +
                        

                        # curr_is_near_DOWN+
                        # s(curr_lastRippleHPC_z, by = curr_is_near_DOWN,k = 5) +
                        
                        # curr_ripple_quantile +
                        # s(curr_lastRippleHPC_z, by = curr_ripple_quantile,k = 5) +
                        # s(curr_lastRippleV1PRE, k = 5) +
                        # s(prev_lastRippleV1PRE, k = 5) +

                        s(prev_lastRippleHPC_z, k = 5) +
                        s(curr_lastRippleHPC_z, k = 5) +
                        
                        s(AnimalID, bs = "re")+
                        s(SessionID, bs = "re"),
                      data = ripple_chain_dataset, method = "fREML", discrete = TRUE)
print(summary(mdl_multipleUP))



# 1. Subset your dataset into two groups
data_near_DOWN <- subset(ripple_chain_dataset, next_is_near_DOWN == 'yes')
data_far_DOWN  <- subset(ripple_chain_dataset, next_is_near_DOWN == 'no')

# 2. Run the BAM model on the "Near" group
mdl_near <- bam(next_firstRippleV1_z ~ 
                  # s(next_firstRippleV1PRE_z,k=5) +
                  prev_is_near_DOWN + s(prev_lastRippleV1_z, by = prev_is_near_DOWN, k = 5) +
                  curr_is_near_DOWN + s(curr_lastRippleV1_z, by = curr_is_near_DOWN, k = 5) +
                  s(AnimalID, bs = "re") + s(SessionID, bs = "re"),
                data = data_near_DOWN, method = "fREML", discrete = TRUE)

# 3. Run the BAM model on the "Far" group
mdl_far <- bam(next_firstRippleV1_z ~ 
                 # s(next_firstRippleV1PRE_z,k=5) +
                 prev_is_near_DOWN + s(prev_lastRippleV1_z, by = prev_is_near_DOWN, k = 5) +
                 curr_is_near_DOWN + s(curr_lastRippleV1_z, by = curr_is_near_DOWN, k = 5) +
                 s(AnimalID, bs = "re") + s(SessionID, bs = "re"),
               data = data_far_DOWN, method = "fREML", discrete = TRUE)

# 4. Compare the outputs
summary(mdl_near)

summary(mdl_far)






# Extract residuals from your existing model
ripple_chain_dataset$residuals <- residuals(mdl_multipleUP)

# Plot residuals against distance to DOWN
library(ggplot2)
ggplot(ripple_chain_dataset, aes(x = next_TimefromFirstRipple, y = residuals)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs")) +
  theme_minimal() +
  labs(title = "Model Residuals vs. Proximity to DOWN",
       x = "Distance to DOWN State", y = "Residuals")





# ---------------------------------------------------------
# Plot 1: Main Effect of last ripple current UP HC
# ---------------------------------------------------------
### HC post
# Calculate scaling factors
raw_breaks <- c(-2,-1,0,1,2)
z_breaks <- sapply(raw_breaks, function(val) {
  ripple_chain_dataset$curr_lastRippleHPC_z[which.min(abs(ripple_chain_dataset$curr_lastRippleHPC - val))]
})

cairo_pdf("currentUP_last_ripple_nearDOWN_nextUP", width = 4.3, height = 4.3)
p_HC_raw <- draw(mdl_multipleUP, select = "s(curr_lastRippleHPC_z):curr_is_near_DOWNyes", residuals = FALSE, rug = FALSE) + 
  theme_bw(base_family = "Arial") + 
  theme(aspect.ratio = 1) +
  scale_x_continuous(breaks = z_breaks, labels = raw_breaks) +
  labs(
    title = "last ripple current UP -> next Up V1", 
    x = "HC bias (current UP last ripple)", 
    y = "Partial Effect"
  )
# dev.new(noRStudioGD = TRUE)
print(p_HC_raw)
dev.off()


# ---------------------------------------------------------
# Plot 1: Main Effect of last ripple previous UP HC
# ---------------------------------------------------------
# Calculate scaling factors
raw_breaks <- c(-2,-1,0,1,2)
z_breaks <- sapply(raw_breaks, function(val) {
  ripple_chain_dataset$prev_lastRippleHPC_z[which.min(abs(ripple_chain_dataset$prev_lastRippleHPC - val))]
})

cairo_pdf("previousUP_last_ripple_nearDOWN_nextUP", width = 4.3, height = 4.3)
p_HC_raw <- draw(mdl_multipleUP, select = "s(prev_lastRippleHPC_z):prev_is_near_DOWNyes", residuals = FALSE, rug = FALSE) + 
  theme_bw(base_family = "Arial") + 
  theme(aspect.ratio = 1) +
  scale_x_continuous(breaks = z_breaks, labels = raw_breaks) +
  labs(
    title = "last ripple previous UP -> next Up V1", 
    x = "HC bias (previous UP last ripple)", 
    y = "Partial Effect"
  )
# dev.new(noRStudioGD = TRUE)
print(p_HC_raw)
dev.off()




# ==============================================================================
# --- MULTI-METRIC EFFECT SIZE CALCULATIONS FOR FINAL MODEL ---
# ==============================================================================
library(parallel)
library(pbapply)
library(dplyr)
library(stringr)
library(tidyr)
library(tidyverse)

message("\nCalculating 4 Effect Size Metrics for mdl_multipleUP (This may take a minute)...")

# --- 1. DEFINE MODEL ARCHITECTURE ---
BASE_TERMS <- c(
  "prev_is_near_DOWN",
  "curr_is_near_DOWN"
)

formula_terms <- c(
  "s(prev_lastRippleHPC_z, by = prev_is_near_DOWN, k = 5)",
  "s(curr_lastRippleHPC_z, by = curr_is_near_DOWN, k = 5)"
)

smooth_labels <- c(
  "s(prev_lastRippleHPC_z):prev_is_near_DOWNno",
  "s(prev_lastRippleHPC_z):prev_is_near_DOWNyes",
  "s(curr_lastRippleHPC_z):curr_is_near_DOWNno",
  "s(curr_lastRippleHPC_z):curr_is_near_DOWNyes"
)

RE_TERMS <- c("s(SessionID, bs = 're')", "s(AnimalID, bs = 're')")
B <- 1000  # Number of Bootstrap Replicates

message(sprintf("\nLaunching %d Case Bootstrap Replicates...", B))

# --- 2. BOOTSTRAP WORKER FUNCTION ---
run_one_bootstrap <- function(rep_id, original_data, base_terms, formula_terms, smooth_labels) {
  
  # 1. Resample data with replacement
  boot_data <- original_data[sample(nrow(original_data), replace = TRUE), ]
  
  # Ensure factor constraints match original configurations explicitly
  boot_data$prev_is_near_DOWN <- factor(boot_data$prev_is_near_DOWN)
  boot_data$curr_is_near_DOWN <- factor(boot_data$curr_is_near_DOWN)
  
  # 2. Fit Full Model
  full_form <- as.formula(paste(
    "next_earlyUPV1_z ~",
    paste(c(base_terms, formula_terms, RE_TERMS), collapse = " + ")
  ))
  
  mdl_f <- tryCatch({
    mgcv::bam(full_form, data = boot_data, method = "fREML", discrete = TRUE)
  }, error = function(e) return(NULL))
  
  if (is.null(mdl_f)) return(NULL)
  
  f_sum <- summary(mdl_f)
  full_dev <- f_sum$dev.expl * 100
  res_df <- mdl_f$df.residual
  sum_tab <- as.data.frame(f_sum$s.table)
  sum_tab$Term <- rownames(sum_tab)
  
  run_res <- data.frame(Term = smooth_labels, Partial_Deviance = NA, Eta_Sq_Partial = NA,
                        Amplitude = NA, RMS = NA, Rep = rep_id)
  
  # 3. Process every individual expanded smooth term
  for (j in seq_along(smooth_labels)) {
    target_term <- smooth_labels[j]
    
    # A. Partial Deviance (Drop-One Smooth Term block)
    if (grepl("prev_lastRippleHPC_z", target_term)) {
      act_terms <- c(base_terms, formula_terms[2], RE_TERMS) 
    } else {
      act_terms <- c(base_terms, formula_terms[1], RE_TERMS) 
    }
    
    r_form <- as.formula(paste("next_earlyUPV1_z ~", paste(act_terms, collapse = " + ")))
    mdl_r <- tryCatch({
      mgcv::bam(r_form, data = boot_data, method = "fREML", discrete = TRUE)
    }, error = function(e) return(NULL))
    
    if (!is.null(mdl_r)) {
      run_res$Partial_Deviance[j] <- full_dev - (summary(mdl_r)$dev.expl * 100)
    }
    
    # B. Partial Eta-Squared 
    t_row <- sum_tab[sum_tab$Term == target_term, ]
    if (nrow(t_row) == 1) {
      run_res$Eta_Sq_Partial[j] <- (t_row$F * t_row$edf) / ((t_row$F * t_row$edf) + res_df)
    }
    
    # C & D. Amplitude and RMS via prediction grid isolation
    f_level  <- if(grepl("yes", target_term)) "yes" else "no"
    
    # 1. Create a base grid for a single continuous sequence
    if (grepl("prev_lastRippleHPC_z", target_term)) {
      x_seq <- seq(min(boot_data$prev_lastRippleHPC_z, na.rm = TRUE), 
                   max(boot_data$prev_lastRippleHPC_z, na.rm = TRUE), 
                   length.out = 50)
      
      # We create rows for BOTH "no" and "yes" so mgcv doesn't see 0 rows for either smooth
      grid_clean <- expand.grid(
        prev_lastRippleHPC_z = x_seq,
        prev_is_near_DOWN    = factor(c("no", "yes"), levels = levels(boot_data$prev_is_near_DOWN)),
        curr_lastRippleHPC_z = mean(boot_data$curr_lastRippleHPC_z, na.rm = TRUE),
        curr_is_near_DOWN    = factor(c("no", "yes"), levels = levels(boot_data$curr_is_near_DOWN))[1], 
        SessionID            = boot_data$SessionID[1],
        AnimalID             = boot_data$AnimalID[1]
      )
      # Ensure the other factor also has both levels represented to prevent its own smooth from crashing
      grid_clean$curr_is_near_DOWN <- factor(rep(c("no", "yes"), length.out = nrow(grid_clean)), 
                                             levels = levels(boot_data$curr_is_near_DOWN))
      
      # Keep track of which rows match our target level for later slicing
      target_rows <- which(grid_clean$prev_is_near_DOWN == f_level)
      
    } else {
      x_seq <- seq(min(boot_data$curr_lastRippleHPC_z, na.rm = TRUE), 
                   max(boot_data$curr_lastRippleHPC_z, na.rm = TRUE), 
                   length.out = 50)
      
      grid_clean <- expand.grid(
        curr_lastRippleHPC_z = x_seq,
        curr_is_near_DOWN    = factor(c("no", "yes"), levels = levels(boot_data$curr_is_near_DOWN)),
        prev_lastRippleHPC_z = mean(boot_data$prev_lastRippleHPC_z, na.rm = TRUE),
        prev_is_near_DOWN    = factor(c("no", "yes"), levels = levels(boot_data$prev_is_near_DOWN))[1],
        SessionID            = boot_data$SessionID[1],
        AnimalID             = boot_data$AnimalID[1]
      )
      grid_clean$prev_is_near_DOWN <- factor(rep(c("no", "yes"), length.out = nrow(grid_clean)), 
                                             levels = levels(boot_data$prev_is_near_DOWN))
      
      target_rows <- which(grid_clean$curr_is_near_DOWN == f_level)
    }
    
    # 2. Extract predictions using the full non-empty design matrix
    Xp <- predict(mdl_f, newdata = grid_clean, type = "lpmatrix")
    smooth_cols <- grep(target_term, colnames(Xp), fixed = TRUE)
    
    if (length(smooth_cols) > 0) {
      Xp_isolated <- matrix(0, nrow = nrow(Xp), ncol = ncol(Xp))
      Xp_isolated[, smooth_cols] <- Xp[, smooth_cols]
      
      # 3. Calculate full isolated fit, then subset only the rows matching our target level
      fit_full <- Xp_isolated %*% coef(mdl_f)
      fit_isolated <- fit_full[target_rows]
      
      run_res$Amplitude[j] <- max(fit_isolated) - min(fit_isolated)
      run_res$RMS[j]       <- sqrt(mean(fit_isolated^2))
    }
  }
  
  return(run_res)
}


# --- 3. RUN SIMULATION LOOP ---
set.seed(42)
boot_results_list <- pblapply(1:B, run_one_bootstrap, 
                              original_data = ripple_chain_dataset, 
                              formula_terms = formula_terms, 
                              base_terms    = BASE_TERMS, 
                              smooth_labels = smooth_labels)

boot_results_list <- boot_results_list[!sapply(boot_results_list, is.null)]
boot_df <- do.call(rbind, boot_results_list)

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

write.csv(raw_iterations_clean, "MultipleUP_GAM_model_raw_bootstrap_iterations.csv", row.names = FALSE)
message("Raw bootstrap iteration file saved successfully.")

# ==============================================================================
# --- 4. GENERATE PLOTS AND EXPORT SUMMARIES ---
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
                               "Peak-to-Trough Amplitude", "RMS Effect")),
    Term = factor(Term)
  )

p_bars_with_ci <- ggplot(plot_data, aes(x = Term, y = Val, fill = Term)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60", alpha = 0.5) +
  geom_bar(stat = "identity", show.legend = FALSE, alpha = 0.85, width = 0.6) +
  geom_errorbar(aes(ymin = Lwr, ymax = Upr), width = 0.2, color = "black", linewidth = 0.6) +
  coord_flip() +
  facet_wrap(~Metric, scales = "free_x", ncol = 2) +
  scale_fill_viridis_d(option = "mako", direction = -1) +
  theme_bw() +
  labs(
    title = "mdl_multipleUP Effect Size Metrics",
    subtitle = "Denotes bootstrap medians and 95% non-parametric CIs.",
    x = NULL, y = NULL
  ) +
  theme(
    strip.text = element_text(face = "bold", size = 10, colour = "black"),
    axis.text.y = element_text(size = 9, color = "black"),
    panel.spacing = unit(1.2, "lines")
  )

print(p_bars_with_ci)

cairo_pdf("MultipleUP_Model_Effect_Sizes_With_CI.pdf", width = 10, height = 6.5)
print(p_bars_with_ci)
dev.off()

final_dashboard_data <- raw_iterations_clean %>%
  group_by(Term) %>%
  summarise(across(c(Deviance_Value, Eta_Sq_Value, Amplitude_Value, RMS_Value),
                   list(Val = ~median(.x, na.rm=TRUE),
                        Lwr = ~quantile(.x, probs=0.025, na.rm=TRUE),
                        Upr = ~quantile(.x, probs=0.975, na.rm=TRUE)),
                   .names = "{.col}__{.fn}"))

flat_bootstrap_results <- final_dashboard_data %>%
  pivot_longer(cols = -Term, names_to = "Combined", values_to = "Value") %>%
  separate(Combined, into = c("Metric", "Stat"), sep = "__") %>%
  mutate(New_Col_Name = paste0(Metric, "_", Stat)) %>%
  select(-Metric, -Stat) %>%
  pivot_wider(names_from = New_Col_Name, values_from = Value)

model_stats <- as.data.frame(summary(mdl_multipleUP)$s.table) %>%
  mutate(Term = rownames(.))

combined_results <- model_stats %>%
  left_join(flat_bootstrap_results, by = "Term")

write.csv(combined_results, "MultipleUP_GAM_model_CI_output.csv", row.names = FALSE)
message("\nAnalysis Pipeline Complete! Saved 'MultipleUP_GAM_model_CI_output.csv' and figures.")