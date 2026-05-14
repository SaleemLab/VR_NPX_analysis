%% GAM - matlab version uses trees, not great
load('V1_HPC_reactivation_coherence_lme1_raw.mat');


%% 1. Outlier Removal
% Define the continuous variables that are vulnerable to extreme noise
outlierVars = {'RipplePower', 'SpindlePower_Match', 'SpindlePower_NonMatch', ...
               'HPC_logodds', 'V1_logodds_PRE', 'V1_logodds'};

% Find outliers (trimming the extreme 0.1% tails)
outlierMask = isoutlier(tbl{:, outlierVars}, 'percentiles', [0.01, 99.99]);
rowsWithOutliers = any(outlierMask, 2);

% Create the clean table
cleanTbl = tbl(~rowsWithOutliers, :);
fprintf('Removed %d rows containing extreme outliers.\n', sum(rowsWithOutliers));

%% 2. Z-Score Continuous Predictors (CRITICAL for interactions)
% Now that outliers are gone, we can safely compute means and standard deviations
cleanTbl.RipplePower_Z = zscore(cleanTbl.RipplePower);
cleanTbl.SpindlePower_Match_Z = zscore(cleanTbl.SpindlePower_Match);
cleanTbl.SpindlePower_NonMatch_Z = zscore(cleanTbl.SpindlePower_NonMatch);
cleanTbl.HPC_logodds_Z = zscore(cleanTbl.HPC_logodds);
cleanTbl.V1_logodds_PRE_Z = zscore(cleanTbl.V1_logodds_PRE);
cleanTbl.V1_logodds_Z = zscore(cleanTbl.V1_logodds); % Post-ripple


%% Calculate coherence
% The manuscript specifically uses z-scored log-odds for its LME models
cleanTbl.HPC_logodds_Z = zscore(cleanTbl.HPC_logodds);
cleanTbl.V1_logodds_Z  = zscore(cleanTbl.V1_logodds);     % Post-ripple
cleanTbl.V1_logodds_PRE_Z = zscore(cleanTbl.V1_logodds_PRE); % Pre-ripple

% 2. Calculate Per-Ripple Coherence (Agreement Score)
% This score represents the "coupling strength" of each individual event.
% Positive = Coherent (Matching content)
% Negative = Incoherent (Mismatching content)
% Near Zero = One or both regions are neutral
cleanTbl.Event_Coherence_Post = cleanTbl.HPC_logodds_Z .* cleanTbl.V1_logodds_Z;
cleanTbl.Event_Coherence_Pre  = cleanTbl.HPC_logodds_Z .* cleanTbl.V1_logodds_PRE_Z;

% Verify the distribution
fprintf('Mean Event Coherence (Post): %.3f\n', mean(cleanTbl.Event_Coherence_Post));

%% GAM
% 1. PREPARE DATA AND VARIABLES
cleanTbl.JointState = categorical(cleanTbl.JointState);
cleanTbl.AnimalID   = categorical(cleanTbl.AnimalID);

predictorNames = {'RipplePower_Z', 'SpindlePower_Match_Z', ...
                  'SpindlePower_NonMatch_Z', 'JointState', 'AnimalID', 'SessionID'};

% 2. FIT THE GAM (Fixed Optimization Struct)
disp('Fitting Categorical GAM...');

% We use 'Verbose', 0 instead of 'Display', 'none' to fix the error
optOptions = struct('MaxObjectiveEvaluations', 5, 'Verbose', 0);

% We remove 'OptimizeHyperparameters' and set manual complexity
mdl_gam = fitrgam(cleanTbl, 'Event_Coherence_Post', ...
    'PredictorNames', predictorNames, ...
    'CategoricalPredictors', {'JointState', 'AnimalID', 'SessionID'}, ...
    'MaxNumSplitsPerPredictor', 1, ...  % Allows for a "staircase" instead of one step
    'NumTreesPerPredictor', 500, ...     % More trees = smoother, more detailed curves
    'InitialLearnRateForPredictors', 0.01); % Slow learning prevents "cliff" artifacts

% 3. CALCULATE FEATURE IMPORTANCE (Permutation/Shuffling)
disp('Calculating Permutation Importance (Shuffling Variables)...');
y_hat_base = predict(mdl_gam, cleanTbl);
base_mse = mean((cleanTbl.Event_Coherence_Post - y_hat_base).^2);
delta_mse = zeros(length(predictorNames), 1);

for i = 1:length(predictorNames)
    tempTbl = cleanTbl;
    % Shuffle to destroy relationship
    tempTbl.(predictorNames{i}) = tempTbl.(predictorNames{i})(randperm(height(tempTbl)));
    
    y_hat_shuffled = predict(mdl_gam, tempTbl);
    shuffled_mse = mean((cleanTbl.Event_Coherence_Post - y_hat_shuffled).^2);
    
    % Percentage increase in error
    delta_mse(i) = (shuffled_mse - base_mse) / base_mse * 100; 
end

% 4. VISUALIZATION
figure('Color', 'w', 'Position', [100, 100, 1100, 750]);

% --- PLOT 1: FEATURE IMPORTANCE ---
subplot(2,2,1);
[sortedDelta, idx] = sort(delta_mse, 'descend');
b = barh(sortedDelta, 'FaceColor', [0.2 0.5 0.7]);
set(gca, 'YTick', 1:length(predictorNames), 'YTickLabel', predictorNames(idx), 'YDir', 'reverse');
xlabel('% Increase in MSE (Model Degradation)', 'FontWeight', 'bold');
title('Feature Importance', 'FontSize', 13);
grid on; box off;

% --- PLOT 2: JOINT STATE EFFECT ---
subplot(2,2,2);
plotPartialDependence(mdl_gam, 'JointState');
title('Partial Effect: Joint State', 'FontSize', 13);
ylabel('Contribution to Coherence');
grid on;

% --- PLOT 3: RIPPLE POWER SMOOTH ---
subplot(2,2,3);
plotPartialDependence(mdl_gam, 'RipplePower_Z');
title('Effect: Ripple Power', 'FontSize', 13);
xlabel('Ripple Power (Z)');
grid on;

% --- PLOT 4: MATCHED SPINDLE POWER SMOOTH ---
subplot(2,2,4);
plotPartialDependence(mdl_gam, 'SpindlePower_Match_Z');
title('Effect: Matched Spindle', 'FontSize', 13);
xlabel('Spindle Power (Z)');
grid on;

sgtitle('GAM Analysis of V1-HC Coherence', 'FontSize', 16, 'FontWeight', 'bold');



%% GAM using phase

%% 1. Prepare Predictors
% We drop JointState and use the raw Phase values.
predictorNames = {'RipplePower_Z', 'SpindlePower_Match_Z', 'SpindlePower_NonMatch_Z', ...
                  'SOPhase_Match', 'SOPhase_NonMatch', 'AnimalID', 'SessionID'};

% Create a logical matrix: 7 predictors x 7 predictors
% By default, all are false (no interactions)
p = length(predictorNames);

% Set the specific interaction between Match Phase (4) and Non-Match Phase (5)
% It must be symmetric: (4,5) and (5,4)
interactionMatrix = logical([0 0 0 1 1 0 0]);

% 2. FIT THE PHASE-INTERACTION GAM
disp('Fitting Continuous Phase GAM with Interaction Matrix...');

mdl_gam_phase = fitrgam(cleanTbl, 'Event_Coherence_Post', ...
    'PredictorNames', predictorNames, ...
    'CategoricalPredictors', {'AnimalID', 'SessionID'}, ...
    'InitialLearnRateForPredictors', 0.01, ...
    'MaxNumSplitsPerPredictor', 1, ...
    'NumTreesPerPredictor', 500);%, ...
    % 'Interactions', interactionMatrix); % Now passing the logical matrix

% 3. Calculate Feature Importance
disp('Calculating Importance (Shuffling continuous phases)...');
y_hat_base = predict(mdl_gam_phase, cleanTbl);
base_mse = mean((cleanTbl.Event_Coherence_Post - y_hat_base).^2);

% Bio predictors for the plot
bioPredictors = {'RipplePower_Z', 'SpindlePower_Match_Z', 'SOPhase_Match', 'SOPhase_NonMatch'};
delta_mse = zeros(length(bioPredictors), 1);

for i = 1:length(bioPredictors)
    tempTbl = cleanTbl;
    tempTbl.(bioPredictors{i}) = tempTbl.(bioPredictors{i})(randperm(height(tempTbl)));
    y_hat_shuffled = predict(mdl_gam_phase, tempTbl);
    shuffled_mse = mean((cleanTbl.Event_Coherence_Post - y_hat_shuffled).^2);
    delta_mse(i) = (shuffled_mse - base_mse) / base_mse * 100;
end

disp('Plotting...')
% 1. Create a Subsample for Plotting (The Speed Booster)
% We only use a random subset of the data to calculate the marginal effects.
idxSample = randperm(height(cleanTbl), min(5000, height(cleanTbl)));
plotTbl = cleanTbl(idxSample, :);

% 2. Fast Visualization
figure('Color', 'w', 'Position', [100, 100, 1100, 750]);

% --- Plot 1: Ripple Power (Fast) ---
subplot(2,2,1);
% Notice plotTbl is passed directly as the 3rd argument
plotPartialDependence(mdl_gam_phase, 'RipplePower_Z', plotTbl); 
title('Ripple Power Effect'); grid on;

% --- Plot 2: Matched Spindle Power (Fast) ---
subplot(2,2,2);
plotPartialDependence(mdl_gam_phase, 'SpindlePower_Match_Z', plotTbl);
title('Matched Spindle Effect'); grid on;

% --- Plot 3: Match Phase (Fast & Error-Free) ---
subplot(2,2,3);
plotPartialDependence(mdl_gam_phase, 'SOPhase_Match', plotTbl);
title('V1-HC Coherence by Match Phase');
xlabel('Phase (Radians)'); grid on;

% --- Plot 4: Non-Match Phase (Fast & Error-Free) ---
subplot(2,2,4);
plotPartialDependence(mdl_gam_phase, 'SOPhase_NonMatch', plotTbl);
title('V1-HC Coherence by Non-Match Phase');
xlabel('Phase (Radians)'); grid on;

sgtitle('Continuous Phase Tuning (Optimized Plotting)', 'FontSize', 16, 'FontWeight', 'bold');