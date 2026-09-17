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


my_folder <- "C:/Users/masah/Documents/GitHub/VR_NPX_analysis/UP_DOWN_ripple_GAM/V1_HC_lastRipple_NextUP"

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
dat_clean$log_TimefromLastRipple <- log(dat_clean$TimefromLastRipple + 0.000001)
dat_clean <- dat_clean %>%
  mutate(log_TimefromLastRipple_z = as.numeric(scale(log_TimefromLastRipple)))


########
######## last ripple HC content (especially those close to DOWN state) predicts predict next UP V1
########

# List of the new Z-score columns to filter (Updated to include all relevant ripple predictors)
z_cols <- c(
  "lastRippleV1_z", "lastRippleV1PRE_z",
  "lastRippleHPC_z", "lastRippleHPCPRE_z",
  "nextUPV1_z"
  # "firstRippleV1_z", "firstRippleV1PRE_z",
  # "firstRippleHPC_z", "firstRippleHPCPRE_z"
)

# dat_clean1 <- dat_Ripples %>%
# dat_clean1 <- dat_lateUP_Ripple %>%

dat_clean1 <- dat_clean %>%
  
  filter(
    if_all(all_of(z_cols), ~ abs(.) < z_thresh | is.na(.))
  )

mdl_lastRippleNext <- bam(nextUPV1_z ~ 
                        # s(lastRippleV1_z, bs = "tp",k = 5) +
                        # s(lastRippleV1_z, by = is_near_DOWN) +
                        s(lastRippleV1PRE_z, by = is_near_DOWN) +
                        is_near_DOWN +
                        s(lastRippleHPC_z, by = is_near_DOWN)+

                        # ripple_quantile +
                        # s(lastRippleHPC_z, by = ripple_quantile)+

                        # s(log_TimefromLastRipple_z, bs = "tp",k = 5) +
                        # s(lastRippleHPC_z, bs = "tp",k = 5) +
                        # te(lastRippleHPC_z,log_TimefromLastRipple_z,k=5,bs = c("tp", "tp"))+
                        
                        # 4. Control Term
                        s(AnimalID, bs = "re") +
                        s(SessionID, bs = "re"), 
                      
                      data = dat_clean1, 
                      method = "fREML", 
                      discrete = TRUE, 
                      nthreads = 4)

message("\n--- FINAL MODEL SUMMARY ---")
print(summary(mdl_lastRippleNext))




# ---------------------------------------------------------
# Plot 1: HC_post_z not near DOWN
# ---------------------------------------------------------
### HC post
# Calculate scaling factors
raw_breaks <- c(-2,-1,0,1,2)
z_breaks <- sapply(raw_breaks, function(val) {
  dat_clean$lastRippleHPC_z[which.min(abs(dat_clean$lastRippleHPC - val))]
})


p_HC_raw <- draw(mdl_lastRippleNext, select = "s(lastRippleHPC_z):is_near_DOWNno", residuals = FALSE, rug = FALSE) + 
  theme_bw(base_family = "Arial") + 
  theme(aspect.ratio = 1) +
  scale_x_continuous(breaks = z_breaks, labels = raw_breaks) +
  coord_cartesian(xlim = c(-3.8,3.8), ylim = c(-0.21, 0.21),expand = FALSE) +
  
  labs(
    title = "HC ripple bias predicting nextUP V1 track bias (last ripple >100ms away from DOWN)", 
    x = "HC bias", 
    y = "Partial Effect"
  )
# dev.new(noRStudioGD = TRUE)
print(p_HC_raw)
cairo_pdf("lastRippleHPC_nextUP_away_from_DOWN.pdf", width = 4.3, height = 4.3)
print(p_HC_raw)
dev.off()



# ---------------------------------------------------------
# Plot 1: HC_post_z near DOWN
# ---------------------------------------------------------

### HC post
# Calculate scaling factors
raw_breaks <- c(-2,-1,0,1,2)
z_breaks <- sapply(raw_breaks, function(val) {
  dat_clean$lastRippleHPC_z[which.min(abs(dat_clean$lastRippleHPC - val))]
})


p_HC_raw <- draw(mdl_lastRippleNext, select = "s(lastRippleHPC_z):is_near_DOWNyes", residuals = FALSE, rug = FALSE) + 
  theme_bw(base_family = "Arial") + 
  theme(aspect.ratio = 1) +
  scale_x_continuous(breaks = z_breaks, labels = raw_breaks) +
  coord_cartesian(xlim = c(-3.8,3.8), ylim = c(-0.21, 0.21),expand = FALSE) +
  labs(
    title = "HC ripple bias predicting nextUP V1 track bias (last ripple <100ms away from DOWN)", 
    x = "HC bias", 
    y = "Partial Effect"
  )
print(p_HC_raw)
# dev.new(noRStudioGD = TRUE)
print(p_HC_raw)
cairo_pdf("lastRippleHPC_nextUP_close_to_DOWN.pdf", width = 4.3, height = 4.3)
dev.off()





# ---------------------------------------------------------
# Plot 3: V1 PRE away from DOWN
# ---------------------------------------------------------

### HC post
# Calculate scaling factors
raw_breaks <- c(-2,-1,0,1,2)
z_breaks <- sapply(raw_breaks, function(val) {
  dat_clean$lastRippleV1PRE_z[which.min(abs(dat_clean$lastRippleV1PRE - val))]
})


p_HC_raw <- draw(mdl_lastRippleNext, select = "s(lastRippleV1PRE_z):is_near_DOWNno", residuals = FALSE, rug = FALSE) + 
  theme_bw(base_family = "Arial") + 
  theme(aspect.ratio = 1) +
  scale_x_continuous(breaks = z_breaks, labels = raw_breaks) +
  coord_cartesian(xlim = c(-3.8,3.8), ylim = c(-0.21, 0.21),expand = FALSE) +
  labs(
    title = "PRE V1 ripple bias predicting nextUP V1 track bias (last ripple >100ms away from DOWN)", 
    x = "HC bias", 
    y = "Partial Effect"
  )
print(p_HC_raw)

cairo_pdf("lastRippleV1PRE_nextUP_away_from_DOWN.pdf", width = 4.3, height = 4.3)
print(p_HC_raw)
dev.off()


# ---------------------------------------------------------
# Plot 4: V1 PRE near DOWN
# ---------------------------------------------------------

### HC post
# Calculate scaling factors
raw_breaks <- c(-2,-1,0,1,2)
z_breaks <- sapply(raw_breaks, function(val) {
  dat_clean$lastRippleV1PRE_z[which.min(abs(dat_clean$lastRippleV1PRE - val))]
})


p_HC_raw <- draw(mdl_lastRippleNext, select = "s(lastRippleV1PRE_z):is_near_DOWNyes", residuals = FALSE, rug = FALSE) + 
  theme_bw(base_family = "Arial") + 
  theme(aspect.ratio = 1) +
  scale_x_continuous(breaks = z_breaks, labels = raw_breaks) +
  coord_cartesian(xlim = c(-4,4), ylim = c(-0.4, 0.52),expand = FALSE) +
  labs(
    title = "PRE V1 ripple bias predicting nextUP V1 track bias (last ripple <100ms away from DOWN)", 
    x = "HC bias", 
    y = "Partial Effect"
  )
print(p_HC_raw)

cairo_pdf("lastRippleV1PRE_nextUP_close_to_DOWN.pdf", width = 4.3, height = 4.3)
print(p_HC_raw)
dev.off()



# ==============================================================================
# --- COMBINED EFFECT SURFACE: lastRippleV1PRE_z + lastRippleHPC_z | is_near_DOWN == "yes" ---
# ==============================================================================

# 1. Grid across observed (z-scored) range
grid_range_v1  <- seq(-2.5, 2.5, length.out = 50)
grid_range_hpc <- seq(-2.5, 2.5, length.out = 50)

pred_grid_total <- expand.grid(
  lastRippleV1PRE_z = grid_range_v1,
  lastRippleHPC_z   = grid_range_hpc,
  is_near_DOWN      = factor("yes", levels = levels(dat_clean1$is_near_DOWN)),
  AnimalID          = dat_clean1$AnimalID[1],
  SessionID         = dat_clean1$SessionID[1]
)

# 2. Extract term components matrix
term_preds <- predict(mdl_lastRippleNext, newdata = pred_grid_total, type = "terms")

# 3. Sum only the "yes" by-level smooths of interest
pred_grid_total$Reconstructed_Effect <-
  term_preds[, "s(lastRippleV1PRE_z):is_near_DOWNyes"] +
  term_preds[, "s(lastRippleHPC_z):is_near_DOWNyes"]

# 4. Plot combined surface
p_combined_yes <- ggplot(pred_grid_total, aes(x = lastRippleV1PRE_z, y = lastRippleHPC_z, fill = Reconstructed_Effect)) +
  geom_tile() +
  geom_contour(aes(z = Reconstructed_Effect), color = "black", alpha = 0.2) +
  scale_fill_gradient2(high = "blue", mid = "white", low = "firebrick", midpoint = 0, name = "Total\nEffect") +
  theme_minimal() +
  labs(
    title = "Combined Effect Surface (is_near_DOWN = yes)",
    subtitle = "Reconstructed: s(lastRippleV1PRE_z):yes + s(lastRippleHPC_z):yes",
    x = "lastRippleV1PRE_z (z-scored)",
    y = "lastRippleHPC_z (z-scored)"
  )
print(p_combined_yes)

cairo_pdf("lastRippleV1PRE_and_HPC_combined_nextUP_close_to_DOWN.pdf", width = 4.3, height = 4.3)
print(p_combined_yes)
dev.off()


# ==============================================================================
# --- COMBINED EFFECT SURFACE: lastRippleV1PRE_z + lastRippleHPC_z | is_near_DOWN == "no" ---
# ==============================================================================

# 1. Grid across observed (z-scored) range
grid_range_v1  <- seq(-2.5, 2.5, length.out = 50)
grid_range_hpc <- seq(-2.5, 2.5, length.out = 50)

pred_grid_total <- expand.grid(
  lastRippleV1PRE_z = grid_range_v1,
  lastRippleHPC_z   = grid_range_hpc,
  is_near_DOWN      = factor("no", levels = levels(dat_clean1$is_near_DOWN)),
  AnimalID          = dat_clean1$AnimalID[1],
  SessionID         = dat_clean1$SessionID[1]
)

# 2. Extract term components matrix
term_preds <- predict(mdl_lastRippleNext, newdata = pred_grid_total, type = "terms")

# 3. Sum only the "no" by-level smooths of interest
pred_grid_total$Reconstructed_Effect <-
  term_preds[, "s(lastRippleV1PRE_z):is_near_DOWNno"] +
  term_preds[, "s(lastRippleHPC_z):is_near_DOWNno"]

# 4. Plot combined surface
p_combined_no <- ggplot(pred_grid_total, aes(x = lastRippleV1PRE_z, y = lastRippleHPC_z, fill = Reconstructed_Effect)) +
  geom_tile() +
  geom_contour(aes(z = Reconstructed_Effect), color = "black", alpha = 0.2) +
  scale_fill_gradient2(high = "blue", mid = "white", low = "firebrick", midpoint = 0, name = "Total\nEffect") +
  theme_minimal() +
  labs(
    title = "Combined Effect Surface (is_near_DOWN = no)",
    subtitle = "Reconstructed: s(lastRippleV1PRE_z):no + s(lastRippleHPC_z):no",
    x = "lastRippleV1PRE_z (z-scored)",
    y = "lastRippleHPC_z (z-scored)"
  )

print(p_combined_no)



# ==============================================================================
# --- MULTI-METRIC EFFECT SIZE CALCULATIONS FOR FINAL MODEL ---
# ==============================================================================
library(parallel)
library(pbapply)
library(dplyr)
library(stringr)
library(tidyr)
library(tidyverse)

message("\nCalculating 4 Effect Size Metrics for mdl_lastRippleNext (This may take a minute)...")

# --- 0. PRE-PROCESS DATASET WITH ISOLATED NUMERIC DUMMIES ---
dat_clean1 <- dat_clean1 %>%
  mutate(
    DOWN_no  = as.numeric(is_near_DOWN == "no"),
    DOWN_yes = as.numeric(is_near_DOWN == "yes")
  )

# --- 1. DEFINE MODEL ARCHITECTURE ---
BASE_TERMS <- c(
  "is_near_DOWN"
)

formula_terms <- c(
  "s(lastRippleHPC_z, by = DOWN_no, k = 5)",
  "s(lastRippleHPC_z, by = DOWN_yes, k = 5)",
  "s(lastRippleV1PRE_z, by = DOWN_no, k = 5)",
  "s(lastRippleV1PRE_z, by = DOWN_yes, k = 5)"
)

smooth_labels <- c(
  "s(lastRippleHPC_z):DOWN_no",
  "s(lastRippleHPC_z):DOWN_yes",
  "s(lastRippleV1PRE_z):DOWN_no",
  "s(lastRippleV1PRE_z):DOWN_yes"
)

RE_TERMS <- c("s(SessionID, bs = 're')", "s(AnimalID, bs = 're')")
B <- 1000  # Number of Bootstrap Replicates

message(sprintf("\nLaunching %d Case Bootstrap Replicates...", B))

# --- 2. BOOTSTRAP WORKER FUNCTION ---
run_one_bootstrap <- function(rep_id, original_data, base_terms, formula_terms, smooth_labels) {
  
  boot_data <- original_data[sample(nrow(original_data), replace = TRUE), ]
  boot_data$is_near_DOWN <- factor(boot_data$is_near_DOWN)
  
  full_form <- as.formula(paste(
    "nextUPV1_z ~",
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
  
  for (j in seq_along(smooth_labels)) {
    target_term <- smooth_labels[j]
    
    act_terms <- c(base_terms, formula_terms[-j], RE_TERMS)
    
    r_form <- as.formula(paste("nextUPV1_z ~", paste(act_terms, collapse = " + ")))
    mdl_r <- tryCatch({
      mgcv::bam(r_form, data = boot_data, method = "fREML", discrete = TRUE)
    }, error = function(e) return(NULL))
    
    if (!is.null(mdl_r)) {
      run_res$Partial_Deviance[j] <- full_dev - (summary(mdl_r)$dev.expl * 100)
    }
    
    t_row <- sum_tab[sum_tab$Term == target_term, ]
    if (nrow(t_row) == 1) {
      run_res$Eta_Sq_Partial[j] <- (t_row$F * t_row$edf) / ((t_row$F * t_row$edf) + res_df)
    }
    
    v_smooth <- if (grepl("lastRippleV1PRE_z", target_term)) "lastRippleV1PRE_z" else "lastRippleHPC_z"
    v_dummy  <- if (grepl("yes", target_term)) "DOWN_yes" else "DOWN_no"
    
    x_seq <- seq(min(boot_data[[v_smooth]], na.rm = TRUE),
                 max(boot_data[[v_smooth]], na.rm = TRUE),
                 length.out = 50)
    
    grid_clean <- data.frame(x = x_seq)
    colnames(grid_clean) <- v_smooth
    
    all_model_vars <- all.vars(full_form)[-1]
    missing_vars <- setdiff(all_model_vars, colnames(grid_clean))
    
    dummy_cols <- c("DOWN_no", "DOWN_yes")
    
    for (mv in missing_vars) {
      if (mv == "SessionID") {
        grid_clean[[mv]] <- boot_data$SessionID[1]
      } else if (mv == "AnimalID") {
        grid_clean[[mv]] <- boot_data$AnimalID[1]
      } else if (mv == "is_near_DOWN") {
        f_val <- if (grepl("DOWN_yes", v_dummy)) "yes" else "no"
        grid_clean[[mv]] <- factor(f_val, levels = levels(boot_data$is_near_DOWN))
      } else if (mv %in% dummy_cols) {
        grid_clean[[mv]] <- as.numeric(mv == v_dummy)
      } else if (is.numeric(boot_data[[mv]])) {
        grid_clean[[mv]] <- mean(boot_data[[mv]], na.rm = TRUE)
      } else {
        grid_clean[[mv]] <- boot_data[[mv]][1]
      }
    }
    
    Xp <- predict(mdl_f, newdata = grid_clean, type = "lpmatrix")
    smooth_cols <- grep(target_term, colnames(Xp), fixed = TRUE)
    
    if (length(smooth_cols) > 0) {
      Xp_isolated <- matrix(0, nrow = nrow(Xp), ncol = ncol(Xp))
      Xp_isolated[, smooth_cols] <- Xp[, smooth_cols]
      
      fit_isolated <- Xp_isolated %*% coef(mdl_f)
      run_res$Amplitude[j] <- max(fit_isolated) - min(fit_isolated)
      run_res$RMS[j]       <- sqrt(mean(fit_isolated^2))
    }
  }
  
  return(run_res)
}

# --- 3. RUN SIMULATION LOOP ---
set.seed(42)
boot_results_list <- pblapply(1:B, run_one_bootstrap,
                              original_data = dat_clean1,
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

write.csv(raw_iterations_clean, "LastRippleNextUP_GAM_model_raw_bootstrap_iterations.csv", row.names = FALSE)
message("Raw bootstrap iteration file saved successfully.")


raw_iterations_clean <- read.csv("LastRippleNextUP_GAM_model_raw_bootstrap_iterations.csv",
                                 stringsAsFactors = FALSE)


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
    Term = case_when(
      Term == "s(lastRippleHPC_z):DOWN_no"    ~ "HPC | Near DOWN: No (>100ms away)",
      Term == "s(lastRippleHPC_z):DOWN_yes"   ~ "HPC | Near DOWN: Yes (<100ms away)",
      Term == "s(lastRippleV1PRE_z):DOWN_no"  ~ "V1 PRE | Near DOWN: No (>100ms away)",
      Term == "s(lastRippleV1PRE_z):DOWN_yes" ~ "V1 PRE | Near DOWN: Yes (<100ms away)",
      TRUE ~ as.character(Term)
    ),
    Term = factor(Term, levels = c(
      "HPC | Near DOWN: Yes (<100ms away)", "HPC | Near DOWN: No (>100ms away)",
      "V1 PRE | Near DOWN: Yes (<100ms away)", "V1 PRE | Near DOWN: No (>100ms away)"
    ))
  )


# Invisible anchor points to force each facet's value-axis range
limits_df <- data.frame(
  Metric = factor(
    rep(c("Deviance Explained (%)", "Partial Eta-Squared",
          "Peak-to-Trough Amplitude", "RMS Effect"), each = 2),
    levels = levels(plot_data$Metric)
  ),
  Val = c(0, 0.17,      # Deviance Explained (%)
          0, 0.0017,    # Partial Eta-Squared
          0, 0.61,      # Peak-to-Trough Amplitude
          0, 0.6),      # RMS Effect
  Term = plot_data$Term[1]  # dummy placeholder, any valid factor level works
)

p_bars_with_ci <- ggplot(plot_data, aes(x = Term, y = Val, fill = Term)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60", alpha = 0.5) +
  geom_bar(stat = "identity", show.legend = FALSE, alpha = 0.85, width = 0.6) +
  geom_errorbar(aes(ymin = Lwr, ymax = Upr), width = 0.2, color = "black", linewidth = 0.6) +
  geom_blank(data = limits_df, aes(x = Term, y = Val)) +   # forces axis range
  coord_flip() +
  facet_wrap(~Metric, scales = "free_x", ncol = 2) +
  scale_fill_viridis_d(option = "mako", direction = -1) +
  theme_bw() +
  labs(
    title = "mdl_lastRipple Effect Size Metrics (Isolated Dummy Matrix Fix)",
    subtitle = "Denotes bootstrap medians and true isolated 95% non-parametric CIs.",
    x = NULL, y = NULL
  ) +
  theme(
    strip.text = element_text(face = "bold", size = 10, colour = "black"),
    axis.text.y = element_text(size = 9, color = "black"),
    panel.spacing = unit(1.2, "lines")
  )

print(p_bars_with_ci)

cairo_pdf("V1_HC_combinedLastRippleNextUP_Model_Effect_Sizes_With_CI.pdf", width = 10, height = 6.5)
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

model_stats <- as.data.frame(summary(mdl_lastRippleNext)$s.table) %>%
  mutate(Term = rownames(.))

flat_bootstrap_results <- flat_bootstrap_results %>%
  mutate(Term = case_when(
    Term == "s(lastRippleHPC_z):DOWN_no"    ~ "s(lastRippleHPC_z):is_near_DOWNno",
    Term == "s(lastRippleHPC_z):DOWN_yes"   ~ "s(lastRippleHPC_z):is_near_DOWNyes",
    Term == "s(lastRippleV1PRE_z):DOWN_no"  ~ "s(lastRippleV1PRE_z):is_near_DOWNno",
    Term == "s(lastRippleV1PRE_z):DOWN_yes" ~ "s(lastRippleV1PRE_z):is_near_DOWNyes",
    TRUE ~ Term
  ))

combined_results <- model_stats %>%
  left_join(flat_bootstrap_results, by = "Term")

write.csv(combined_results, "V1_HC_combined_LastRippleNextUP_GAM_model_CI_output.csv", row.names = FALSE)
message("\nAnalysis Pipeline Complete! Saved 'V1_HC_combined_LastRippleNextUP_GAM_model_CI_output.csv' and figures.")
