%% 0. Load
load('V1_HPC_reactivation_coherence_lme1_raw.mat');


%% 1. Outlier Removal
% Define the continuous variables that are vulnerable to extreme noise
outlierVars = {'RipplePower', 'SpindlePower_Match', 'SpindlePower_NonMatch', ...
               'HPC_logodds', 'V1_logodds_PRE', 'V1_logodds'};

% Find outliers (trimming the extreme 0.1% tails)
outlierMask = isoutlier(tbl{:, outlierVars}, 'percentiles', [0.1, 99.9]);
rowsWithOutliers = any(outlierMask, 2);

% Create the clean table
cleanTbl = tbl(~rowsWithOutliers, :);
fprintf('Removed %d rows containing extreme outliers.\n', sum(rowsWithOutliers));

%% 2. Z-Score Continuous Predictors (CRITICAL for interactions)
% Now that outliers are gone, we can safely compute means and standard deviations
cleanTbl.RipplePower_Z = zscore(cleanTbl.RipplePower);
cleanTbl.SpindlePower_Match_Z = zscore(cleanTbl.SpindlePower_Match);
cleanTbl.HPC_logodds_Z = zscore(cleanTbl.HPC_logodds);
cleanTbl.V1_logodds_PRE_Z = zscore(cleanTbl.V1_logodds_PRE);
cleanTbl.V1_logodds_Z = zscore(cleanTbl.V1_logodds); % Post-ripple




%% 3. Add Reactivation Classification (Top/Bottom 25%)
% The paper defines "Track L" and "Track R" events as the top 25% of positive 
% and negative hippocampal biases, respectively.
p75 = prctile(cleanTbl.HPC_logodds, 75);
p25 = prctile(cleanTbl.HPC_logodds, 25);

cleanTbl.Reactivation_Type = repmat("Neutral", height(cleanTbl), 1);
cleanTbl.Reactivation_Type(cleanTbl.HPC_logodds >= p75) = "Track L Bias";
cleanTbl.Reactivation_Type(cleanTbl.HPC_logodds <= p25) = "Track R Bias";
cleanTbl.Reactivation_Type = categorical(cleanTbl.Reactivation_Type);

%% 4. Add Ripple Power Quartiles
% Splitting data into Low vs High ripple power
p_power = prctile(cleanTbl.RipplePower, [25, 50, 75]);

cleanTbl.Ripple_Quartile = repmat("Low (0-25%)", height(cleanTbl), 1);
cleanTbl.Ripple_Quartile(cleanTbl.RipplePower > p_power(1) & cleanTbl.RipplePower <= p_power(2)) = "Mid-Low (25-50%)";
cleanTbl.Ripple_Quartile(cleanTbl.RipplePower > p_power(2) & cleanTbl.RipplePower <= p_power(3)) = "Mid-High (50-75%)";
cleanTbl.Ripple_Quartile(cleanTbl.RipplePower > p_power(3)) = "High (75-100%)";

% Convert to an ordered categorical variable for plotting
cleanTbl.Ripple_Quartile = categorical(cleanTbl.Ripple_Quartile, ...
    {'Low (0-25%)', 'Mid-Low (25-50%)', 'Mid-High (50-75%)', 'High (75-100%)'}, 'Ordinal', true);

%% 5. Add the "Signed" Trough Phase (-1 to 1)
% Linearized UP-state metric: Trough (+1) vs Peak (-1)
cleanTbl.SO_Match_sin = sin(cleanTbl.SOPhase_Match);
cleanTbl.SO_Match_cos = cos(cleanTbl.SOPhase_Match);
cleanTbl.SO_NonMatch_sin = sin(cleanTbl.SOPhase_NonMatch);
cleanTbl.SO_NonMatch_cos = cos(cleanTbl.SOPhase_NonMatch);

refTrough = pi; 
ang_M = atan2(cleanTbl.SO_Match_sin, cleanTbl.SO_Match_cos);
cleanTbl.SOPhase_Match_Linearized = cos(ang_M - refTrough);

ang_NM = atan2(cleanTbl.SO_NonMatch_sin, cleanTbl.SO_NonMatch_cos);
cleanTbl.SOPhase_NonMatch_Linearized = cos(ang_NM - refTrough);

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

%% 6. The "Master" Linear Mixed-Effects Model
% We want to see how HC Bias drives V1 Bias, and how that relationship is 
% modulated by the 3-way interaction of Ripple * Spindle * Trough.
% Note: Ensure all continuous inputs to the interaction are the Z-scored versions!

master_formula = ['V1_logodds_Z ~ HPC_logodds_Z * RipplePower_Z * ', ...
                  'SpindlePower_Match_Z * Match_trough + ', ...
                  '(1|AnimalID) + (1|SessionID)'];

lme_master = fitlme(cleanTbl, master_formula);

disp('--- MASTER MODEL: RIPPLE x SPINDLE x TROUGH ---');
disp(lme_master);

%% hierarchical model fitting

% 1. Ensure all continuous variables are Z-scored
cleanTbl.RipplePower_Z = zscore(cleanTbl.RipplePower);
cleanTbl.SpindlePower_Match_Z = zscore(cleanTbl.SpindlePower_Match);
cleanTbl.HPC_logodds_Z = zscore(cleanTbl.HPC_logodds);
cleanTbl.V1_logodds_Z = zscore(cleanTbl.V1_logodds);

% 2. Define the Hierarchical Models

% Model 0: Baseline communication (Does HC drive V1 at all?)
formula0 = 'V1_logodds_Z ~ HPC_logodds_Z + (1|AnimalID) + (1|SessionID)';

% Model 1: The "Ripple Power" effect (Does ripple loudness matter?)
formula1 = 'V1_logodds_Z ~ HPC_logodds_Z * RipplePower_Z + (1|AnimalID) + (1|SessionID)';

% Model 2: The "SO Trough" effect (Does the UP-state gate communication?)
% We add this on top of Ripple Power.
formula2 = 'V1_logodds_Z ~ HPC_logodds_Z * RipplePower_Z + HPC_logodds_Z * Match_trough + (1|AnimalID) + (1|SessionID)';

% Model 3: The "Spindle" effect (Do Spindles add value above SO and Ripples?)
formula3 = ['V1_logodds_Z ~ HPC_logodds_Z * RipplePower_Z + ', ...
            'HPC_logodds_Z * Match_trough + ', ...
            'HPC_logodds_Z * SpindlePower_Match_Z + (1|AnimalID) + (1|SessionID)'];

% Model 4: The "Triple Synergy" (The Global UP-State hypothesis)
% This tests if Ripples, Spindles, and Troughs interact synergistically.
formula4 = ['V1_logodds_Z ~ HPC_logodds_Z * RipplePower_Z * ', ...
            'SpindlePower_Match_Z * Match_trough + (1|AnimalID) + (1|SessionID)'];

% 3. Run Models and Compare Systematically
fprintf('Fitting models...\n');
m0 = fitlme(cleanTbl, formula0);
m1 = fitlme(cleanTbl, formula1);
m2 = fitlme(cleanTbl, formula2);
m3 = fitlme(cleanTbl, formula3);
m4 = fitlme(cleanTbl, formula4);

% Comparison Table
fprintf('\n--- Hierarchical Model Comparison ---\n');
comp1 = compare(m0, m1); % Does Ripple Power help?
comp2 = compare(m1, m2); % Does SO Trough help above Ripple?
comp3 = compare(m2, m3); % Does Spindle Power help above both?
comp4 = compare(m3, m4); % Is the complex interaction justified?

% 4. Print Results for Interpretation
disp('1. Effect of Ripple Power:'); disp(comp1);
disp('2. Effect of SO Trough:');   disp(comp2);
disp('3. Effect of Spindles:');    disp(comp3);
disp('4. Effect of Triple Interaction:'); disp(comp4);


%% predict coherence directly

% Now you can run a model predicting "Coherence" directly instead of V1_logodds
coherence_formula = ['Event_Coherence_Pre ~ RipplePower_Z * ', ...
                     'SpindlePower_Match_Z * Match_trough + ', ...
                     '(1|AnimalID) + (1|SessionID)'];

lme_coherence_direct = fitlme(cleanTbl, coherence_formula);
disp('--- MODEL PREDICTING COHERENCE STRENGTH ---');
disp(lme_coherence_direct);



%% 1. Define Hierarchical Models for Coherence Strength
% Model 0: Null Model (Random effects only - baseline variance)
f0 = 'Event_Coherence_Post ~ 1 + (1|AnimalID) + (1|SessionID)';

% Model 1: Add Ripple Power (Does loudness alone drive agreement? [cite: 109])
f1 = 'Event_Coherence_Post ~ RipplePower_Z + (1|AnimalID) + (1|SessionID)';

% Model 2: Add SO Trough (Does the cortical "UP-state" gate the agreement? [cite: 149, 150])
f2 = 'Event_Coherence_Post ~ RipplePower_Z + Match_trough + (1|AnimalID) + (1|SessionID)';

% Model 3: Add Spindle Power (Does local spindle activity boost agreement? [cite: 128, 133])
f3 = 'Event_Coherence_Post ~ RipplePower_Z + Match_trough + SpindlePower_Match_Z + (1|AnimalID) + (1|SessionID)';

% Model 4: Full Interaction (Do they interact synergistically? [cite: 8])
f4 = 'Event_Coherence_Post ~ RipplePower_Z * SpindlePower_Match_Z * Match_trough + (1|AnimalID) + (1|SessionID)';

% 2. Fit and Compare
fprintf('Fitting Coherence Models...\n');
mc0 = fitlme(cleanTbl, f0,'DummyVarCoding','effects','CovariancePattern', {'Diagonal', 'Diagonal'});
mc1 = fitlme(cleanTbl, f1,'DummyVarCoding','effects','CovariancePattern', {'Diagonal', 'Diagonal'});
mc2 = fitlme(cleanTbl, f2,'DummyVarCoding','effects','CovariancePattern', {'Diagonal', 'Diagonal'});
mc3 = fitlme(cleanTbl, f3,'DummyVarCoding','effects','CovariancePattern', {'Diagonal', 'Diagonal'});
mc4 = fitlme(cleanTbl, f4,'DummyVarCoding','effects','CovariancePattern', {'Diagonal', 'Diagonal'});

fprintf('\n--- Coherence Strength Model Comparison ---\n');
c1 = compare(mc0, mc1); % Effect of Ripple
c2 = compare(mc1, mc2); % Effect of SO Trough
c3 = compare(mc2, mc3); % Effect of Spindle
c4 = compare(mc3, mc4);% Effect of Synergy


%% Drop-a-Term Analysis with Likelihood Ratio Test & Significance Plotting
% Choose Target Window: 'Event_Coherence_Pre' or 'Event_Coherence_Post'
targetVar = 'Event_Coherence_Pre'; % <--- CHANGE THIS TO SWITCH WINDOWS

% 1. STRICT DATA ALIGNMENT
% Ensure we drop NaNs only for the variables we are currently using
varsToKeep = {targetVar, 'RipplePower_Z', 'SpindlePower_Match_Z', ...
              'Match_trough', 'AnimalID', 'SessionID'};
cleanTbl_Model = rmmissing(cleanTbl(:, varsToKeep));

% 2. FIT FULL MODEL WITH MAXIMUM LIKELIHOOD (ML)
explicit_formula = [targetVar, ' ~ ', ...
    'RipplePower_Z + SpindlePower_Match_Z + Match_trough + ', ...
    'RipplePower_Z:SpindlePower_Match_Z + RipplePower_Z:Match_trough + ', ...
    'SpindlePower_Match_Z:Match_trough + ', ...
    'RipplePower_Z:SpindlePower_Match_Z:Match_trough + ', ...
    '(1|AnimalID) + (1|SessionID)'];

mFull_Explicit = fitlme(cleanTbl_Model, explicit_formula, 'FitMethod', 'ML');
LL_full = mFull_Explicit.LogLikelihood;
df_full = mFull_Explicit.NumEstimatedCoefficients; % Degrees of freedom for full model

% 3. SYSTEMATIC DROP-TERM LOOP
components = {
    'RipplePower_Z', ...
    'SpindlePower_Match_Z', ...
    'Match_trough', ...
    'RipplePower_Z:SpindlePower_Match_Z', ...
    'RipplePower_Z:Match_trough', ...
    'SpindlePower_Match_Z:Match_trough', ...
    'RipplePower_Z:SpindlePower_Match_Z:Match_trough'
}';

deltaLL = zeros(length(components), 1);
pValues = zeros(length(components), 1); % Array to store p-values

for i = 1:length(components)
    keepIdx = true(size(components));
    keepIdx(i) = false;
    kept_terms = strjoin(components(keepIdx), ' + ');
    
    % Dynamically inject the target variable into the reduced formula
    reduced_formula = [targetVar, ' ~ ', kept_terms, ' + (1|AnimalID) + (1|SessionID)'];
    
    try
        mReduced = fitlme(cleanTbl_Model, reduced_formula, 'FitMethod', 'ML');
        deltaLL(i) = LL_full - mReduced.LogLikelihood;
        
        % --- Likelihood Ratio Test ---
        df_reduced = mReduced.NumEstimatedCoefficients;
        deltaDF = df_full - df_reduced; % How many parameters did we drop?
        LRStat = 2 * deltaLL(i);        % Chi-square statistic
        pValues(i) = chi2cdf(LRStat, deltaDF, 'upper'); % Calculate p-value
        
    catch
        warning('Failed to fit model when dropping: %s', components{i});
        deltaLL(i) = NaN;
        pValues(i) = NaN;
    end
end

% 4. PLOT RESULTS WITH SIGNIFICANCE STARS
figure('Color', 'w', 'Position', [100 100 750 500]);

% Flip the arrays so they plot top-to-bottom correctly
plot_deltaLL = flipud(deltaLL);
plot_pValues = flipud(pValues);
plot_components = flipud(components);

% Create the horizontal bar chart
b = barh(plot_deltaLL, 'FaceColor', [0.2 0.6 0.5]);
set(gca, 'YTick', 1:length(plot_components), 'YTickLabel', plot_components, ...
    'TickLabelInterpreter', 'none', 'FontSize', 10);
xlabel('\Delta Log-Likelihood (Information Loss)');
title(sprintf('Feature Importance: %s', targetVar), 'Interpreter', 'none');
grid on;

% Add the significance stars to the plot
max_x = max(plot_deltaLL);
offset = max_x * 0.02; % Small padding so text doesn't touch the bar

for i = 1:length(plot_deltaLL)
    p = plot_pValues(i);
    
    % Determine the star string
    if isnan(p)
        stars = ' N/A';
    elseif p < 0.001
        stars = ' ***';
    elseif p < 0.01
        stars = ' **';
    elseif p < 0.05
        stars = ' *';
    else
        stars = ' n.s.';
    end
    
    % Add the text to the plot
    text(plot_deltaLL(i) + offset, i, stars, ...
        'VerticalAlignment', 'middle', 'FontSize', 11, 'FontWeight', 'bold');
end

% Extend the X-axis limit slightly so the stars don't get cut off
xlim([0, max_x * 1.15]);