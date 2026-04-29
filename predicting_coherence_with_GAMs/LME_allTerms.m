%% Linear mixed effects model (grand model formulation)

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

%% 2. Z-Score Continuous Predictors 
% Now that outliers are gone, we can safely compute means and standard deviations
cleanTbl.RipplePower_Z = zscore(cleanTbl.RipplePower);
cleanTbl.SpindlePower_Match_Z = zscore(cleanTbl.SpindlePower_Match);
cleanTbl.SpindlePower_NonMatch_Z = zscore(cleanTbl.SpindlePower_NonMatch);
cleanTbl.HPC_logodds_Z = cleanTbl.HPC_logodds; % already z-scored
cleanTbl.V1_logodds_PRE_Z = cleanTbl.V1_logodds_PRE; % already z-scored
cleanTbl.V1_logodds_Z = cleanTbl.V1_logodds; % Post-ripple, already z-scored


%% Calculate coherence

% 2. Calculate Per-Ripple Coherence (Agreement Score)
% This score represents the "coupling strength" of each individual event.
% Positive = Coherent (Matching content)
% Negative = Incoherent (Mismatching content)
% Near Zero = One or both regions are neutral
cleanTbl.Event_Coherence_Post = cleanTbl.HPC_logodds_Z .* cleanTbl.V1_logodds_Z;
cleanTbl.Event_Coherence_Pre  = cleanTbl.HPC_logodds_Z .* cleanTbl.V1_logodds_PRE_Z;

% In MATLAB (when you create the CSV)
cleanTbl.Event_Coherence_Post_GeoMean = sign(cleanTbl.HPC_logodds_Z .* cleanTbl.V1_logodds_Z) .* sqrt(abs(cleanTbl.HPC_logodds_Z .* cleanTbl.V1_logodds_Z));
% Verify the distribution
fprintf('Mean Event Coherence (Post): %.3f\n', mean(cleanTbl.Event_Coherence_Post));

cleanTbl.Event_Coherence_Post = cleanTbl.Event_Coherence_Post_GeoMean;
% writetable(cleanTbl, 'v1_hc_data_geo2.csv'); % export table for R GAM analysis

%% calculate "probability of agreement" coherence

% 1. Convert log-odds back into probabilities for Track 1
p_V1  = 1 ./ (1 + exp(-cleanTbl.V1_logodds));
p_HPC = 1 ./ (1 + exp(-cleanTbl.HPC_logodds));


% 2. Create a grid of probabilities from 0 to 1
p = linspace(0, 1, 500);
[p_V1_grid, p_HPC_grid] = meshgrid(p, p);

% 3. Calculate the exact probability that both areas agree
prob_agreement = (p_V1_grid .* p_HPC_grid) + ((1 - p_V1_grid) .* (1 - p_HPC_grid));

% 4. Scale it to a -1 to +1 Coherence metric
coherence = (2 .* prob_agreement) - 1;

% 5. Create the figure
figure('Position', [100, 100, 700, 550]);

% 6. Plot the heatmap
imagesc([0 1], [0 1], coherence);
set(gca, 'YDir', 'normal'); % Ensure Y-axis goes from 0 at the bottom to 1 at the top

% Define a custom blue-white-red colormap (similar to 'coolwarm')
n = 256;
blue_to_white = [linspace(0,1,n/2)', linspace(0,1,n/2)', ones(n/2,1)];
white_to_red  = [ones(n/2,1), linspace(1,0,n/2)', linspace(1,0,n/2)'];
custom_cmap = [blue_to_white; white_to_red];
colormap(custom_cmap);

% Set color limits from -1 to 1 to perfectly center the white at 0
clim([-1, 1]); % Note: use caxis([-1 1]) if on an older version of MATLAB (pre-R2022a)

% Add colorbar
cb = colorbar;
ylabel(cb, 'Event Coherence Post');

% Add labels and title
xlabel('Probability Track 1 for V1 (p_{V1})');
ylabel('Probability Track 1 for HPC (p_{HPC})');
title('Heatmap of Event Coherence');

% Add contour lines on top for readability
hold on;
[C, h] = contour(p_V1_grid, p_HPC_grid, coherence, 10, 'LineColor', 'k');
clabel(C, h, 'FontSize', 8, 'Color', 'k');
% Make lines slightly transparent (requires MATLAB R2018b or later)
try h.LineAlpha = 0.3; catch; end 
hold off;

%% independent terms

cleanTbl.Event_Coherence_Post =  cleanTbl.Event_Coherence_Post_GeoMean;

% 1. Refit Model with Non-Match Spindle Power
f = ['Event_Coherence_Post ~ RipplePower_Z + SpindlePower_Match_Z + ', ...
     'SpindlePower_NonMatch_Z + JointState + (1|AnimalID) + (1|AnimalID:SessionID)'];

mdl = fitlme(cleanTbl, f, 'FitMethod', 'REML', 'DummyVarCoding', 'reference');
disp(mdl);

% 2. Prediction Table (Centered on Mean)
states = categorical({'0', '1', '2', '3'});
numStates = length(states);

% Fill predTbl with zeros for all continuous predictors
predTbl = table(states', ...
    zeros(numStates,1), zeros(numStates,1), zeros(numStates,1), ...
    repmat(cleanTbl.AnimalID(1), numStates, 1), ...
    repmat(cleanTbl.SessionID(1), numStates, 1), ...
    'VariableNames', {'JointState', 'RipplePower_Z', 'SpindlePower_Match_Z', ...
                      'SpindlePower_NonMatch_Z', 'AnimalID', 'SessionID'});

[predCoh, predCI] = predict(mdl, predTbl, 'Conditional', false);

% 3. Calculate Deviations
grandMean = mean(predCoh);
state_Est = predCoh - grandMean;
state_CI  = predCI - grandMean;

% Identify indices for Power variables
idx_Rip    = strcmp(mdl.CoefficientNames, 'RipplePower_Z');
idx_Match  = strcmp(mdl.CoefficientNames, 'SpindlePower_Match_Z');
idx_NonMat = strcmp(mdl.CoefficientNames, 'SpindlePower_NonMatch_Z');

try ci_all = coefCI(mdl); catch; ci_all = [mdl.Coefficients.Lower, mdl.Coefficients.Upper]; end

% Extract Estimates and CIs
rip_Est = mdl.Coefficients.Estimate(idx_Rip);
rip_CI  = ci_all(idx_Rip, :);

mat_Est = mdl.Coefficients.Estimate(idx_Match);
mat_CI  = ci_all(idx_Match, :);

non_Est = mdl.Coefficients.Estimate(idx_NonMat);
non_CI  = ci_all(idx_NonMat, :);

% 4. Combine for Plotting
plot_Est = [state_Est; rip_Est; mat_Est; non_Est];
plot_CI  = [state_CI; rip_CI; mat_CI; non_CI];
plot_Names = {'State 0 (Both Peak)', 'State 1 (Non-Mat Trough)', ...
              'State 2 (Mat Trough)', 'State 3 (Both Trough)', ...
              'Ripple Power (Z)', 'Matched Spindle Power (Z)', 'Non-Matched Spindle Power (Z)'};

numVars = length(plot_Est);
yPos = 1:numVars;

% 5. Generate Figure
figure('Color', 'w', 'Position', [100, 100, 900, 600]);
hold on; xline(0, 'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off');

for i = 1:numVars
    % Color Logic: States (Orange), Power (Blue)
    if i <= 4, mColor = '#D95319'; else, mColor = '#0072BD'; end
    
    % Significance Check: If CI crosses 0, make it Gray
    if plot_CI(i,1) <= 0 && plot_CI(i,2) >= 0, mColor = [0.6 0.6 0.6]; end
    
    % Plot CI and Point
    line([plot_CI(i,1), plot_CI(i,2)], [yPos(i), yPos(i)], 'Color', mColor, 'LineWidth', 2.5);
    scatter(plot_Est(i), yPos(i), 130, 'o', 'filled', 'MarkerFaceColor', mColor, ...
            'MarkerEdgeColor', 'w', 'LineWidth', 1.2);
end

% Formatting
set(gca, 'YDir', 'reverse', 'YTick', yPos, 'YTickLabel', plot_Names, ...
    'XGrid', 'on', 'YGrid', 'off', 'FontSize', 11, 'Box', 'off');
xlabel('\Delta Coherence (Effect Size)', 'FontWeight', 'bold');
title('Hierarchical Predictors of V1-HC Coherence', 'FontSize', 14);
hold off;


%% delta ll
% 1. Define the Full and Reduced Formulas
f_full = 'Event_Coherence_Post ~ RipplePower_Z + SpindlePower_Match_Z + SpindlePower_NonMatch_Z + JointState + (1|AnimalID) + (1|AnimalID:SessionID)';

% Define the "Reduced" versions (removing one term at a time)
terms_to_drop = {'RipplePower_Z', 'SpindlePower_Match_Z', 'SpindlePower_NonMatch_Z', 'JointState'};
labels        = {'Ripple Power', 'Matched Spindle', 'Non-Matched Spindle', 'Joint State (ALL)'};

% 2. Fit the Full Model first (using ML)
disp('Fitting Full Model (ML)...');
mdl_full = fitlme(cleanTbl, f_full, 'FitMethod', 'ML', 'DummyVarCoding', 'reference');
L_full = mdl_full.LogLikelihood;

% 3. Loop through and drop each term
lr_stats = zeros(length(terms_to_drop), 1);
p_vals   = zeros(length(terms_to_drop), 1);

for i = 1:length(terms_to_drop)
    fprintf('Testing importance of: %s...\n', labels{i});
    
    % Generate the reduced formula by removing the term
    % We use regexprep to remove the term and any stray '+' signs
    f_reduced = regexprep(f_full, [terms_to_drop{i}, ' \+ '], ''); 
    f_reduced = regexprep(f_reduced, [' \+ ', terms_to_drop{i}], '');
    
    % Fit reduced model
    mdl_red = fitlme(cleanTbl, f_reduced, 'FitMethod', 'ML', 'DummyVarCoding', 'reference');
    
    % Likelihood Ratio Test Math
    L_red = mdl_red.LogLikelihood;
    lr_stats(i) = 2 * (L_full - L_red);
    
    % DF is the difference in number of estimated coefficients
    df = mdl_full.NumCoefficients - mdl_red.NumCoefficients;
    p_vals(i) = 1 - chi2cdf(max(0, lr_stats(i)), df);
end

% 4. Generate the Importance Plot
figure('Color', 'w', 'Position', [150, 150, 700, 500]);
b = barh(lr_stats, 'FaceColor', 'flat');

% Color logic: highlight the most important feature
[~, maxIdx] = max(lr_stats);
b.CData = repmat([0.2 0.4 0.6], length(lr_stats), 1); % Default blue-gray
b.CData(maxIdx, :) = [0.8 0.3 0.1]; % Highlight top predictor in Rust

% Add p-values as text labels
for i = 1:length(lr_stats)
    if p_vals(i) < 0.001
        p_txt = 'p < 0.001';
    else
        p_txt = sprintf('p = %.3f', p_vals(i));
    end
    text(lr_stats(i) + (max(lr_stats)*0.02), i, p_txt, 'FontSize', 10, 'FontWeight', 'bold');
end

% Formatting
set(gca, 'YTick', 1:length(labels), 'YTickLabel', labels, 'YDir', 'reverse');
xlabel('Likelihood Ratio Statistic (\chi^2)', 'FontWeight', 'bold');
title('Feature Importance: Model Fit Degradation upon Removal', 'FontSize', 14);
grid on; set(gca, 'YGrid', 'off', 'XGrid', 'on', 'GridAlpha', 0.2);
box off;
xlim([0, max(lr_stats) * 1.25]); % Make room for text labels

hold off;

%% Test whether interactions can improve the model

%% The Full Interaction Model: Triple Coupling x State
% 1. Define Formulas
f_base = 'Event_Coherence_Post ~ RipplePower_Z + SpindlePower_Match_Z + JointState + (1|AnimalID) + (1|AnimalID:SessionID)';
f_full = 'Event_Coherence_Post ~ RipplePower_Z * SpindlePower_Match_Z * JointState + (1|AnimalID) + (1|AnimalID:SessionID)';

% 2. Fit Models (Using ML for valid comparison)
disp('Fitting Base Model (ML)...');
mdl_base = fitlme(cleanTbl, f_base, 'FitMethod', 'ML', 'DummyVarCoding', 'reference');

disp('Fitting Full Interaction Model (ML)... This may take a moment.');
mdl_full = fitlme(cleanTbl, f_full, 'FitMethod', 'ML', 'DummyVarCoding', 'reference');

% 3. Manual Likelihood Ratio Test
L0 = mdl_base.LogLikelihood;
L1 = mdl_full.LogLikelihood;
LR_stat = 2 * (L1 - L0);

% Corrected: Use NumCoefficients for LinearMixedModel objects
deltaDF = mdl_full.NumCoefficients - mdl_base.NumCoefficients;

p_val = 1 - chi2cdf(max(0, LR_stat), deltaDF); 

fprintf('\n======================================================\n');
fprintf('   FULL INTERACTION TEST: Ripple * Spindle * State\n');
fprintf('======================================================\n');
fprintf('Log-Likelihood (Base):     %.4f\n', L0);
fprintf('Log-Likelihood (Full):     %.4f\n', L1);
fprintf('LR Stat:                   %.4f\n', LR_stat);
fprintf('Degrees of Freedom Added:  %d\n', deltaDF);
fprintf('p-value:                   %.6e\n', p_val);
fprintf('BIC (Base vs Full):        %.1f vs %.1f\n', mdl_base.ModelCriterion.AIC, mdl_full.ModelCriterion.AIC);
fprintf('======================================================\n');

% 4. Identify the "Winning" interactions
disp('--- Significant Interaction Terms (p < 0.05) ---');
coeffs = mdl_full.Coefficients;
sigIdx = coeffs.pValue < 0.05 & contains(coeffs.Name, ':');
if any(sigIdx)
    disp(coeffs(sigIdx, :));
else
    disp('No interaction terms reached significance.');
end








%% Compare Random Effect Structures
% 1. STRICT DATA TYPE ENFORCEMENT
cleanTbl.AnimalID   = categorical(cleanTbl.AnimalID);
cleanTbl.SessionID  = categorical(cleanTbl.SessionID);
cleanTbl.JointState = categorical(cleanTbl.JointState);

% 2. DEFINE BASE FORMULA
f_base = 'Event_Coherence_Post ~ RipplePower_Z + SpindlePower_Match_Z + JointState';

% 3. FIT THE 4 COMPETING MODELS (Using ML for valid comparison)
disp('Fitting Model 1: Animal ID Only...');
mdl_anim = fitlme(cleanTbl, [f_base, ' + (1|AnimalID)'], ...
    'FitMethod', 'ML', 'DummyVarCoding', 'reference');

disp('Fitting Model 2: Session ID Only...');
mdl_sess = fitlme(cleanTbl, [f_base, ' + (1|SessionID)'], ...
    'FitMethod', 'ML', 'DummyVarCoding', 'reference');

disp('Fitting Model 3: Original / Crossed (Animal & Session Independent)...');
mdl_cross = fitlme(cleanTbl, [f_base, ' + (1|AnimalID) + (1|SessionID)'], ...
    'FitMethod', 'ML', 'DummyVarCoding', 'reference', ...
    'CovariancePattern', {'Diagonal', 'Diagonal'});

disp('Fitting Model 4: Nested (Session within Animal)...');
mdl_nest = fitlme(cleanTbl, [f_base, ' + (1|AnimalID) + (1|AnimalID:SessionID)'], ...
    'FitMethod', 'ML', 'DummyVarCoding', 'reference');

% 4. EXTRACT METRICS
LL_A = mdl_anim.LogLikelihood;   AIC_A = mdl_anim.ModelCriterion.AIC;   BIC_A = mdl_anim.ModelCriterion.BIC;
LL_S = mdl_sess.LogLikelihood;   AIC_S = mdl_sess.ModelCriterion.AIC;   BIC_S = mdl_sess.ModelCriterion.BIC;
LL_C = mdl_cross.LogLikelihood;  AIC_C = mdl_cross.ModelCriterion.AIC;  BIC_C = mdl_cross.ModelCriterion.BIC;
LL_N = mdl_nest.LogLikelihood;   AIC_N = mdl_nest.ModelCriterion.AIC;   BIC_N = mdl_nest.ModelCriterion.BIC;

% 5. DISPLAY MASTER TABLE
fprintf('\n=======================================================================\n');
fprintf('                   RANDOM EFFECTS 4-WAY SHOOTOUT\n');
fprintf('=======================================================================\n');
fprintf('Model                     | LogLikelihood |    AIC    |    BIC    \n');
fprintf('-----------------------------------------------------------------------\n');
fprintf('1. Animal Only            | %13.4f | %9.1f | %9.1f \n', LL_A, AIC_A, BIC_A);
fprintf('2. Session Only           | %13.4f | %9.1f | %9.1f \n', LL_S, AIC_S, BIC_S);
fprintf('3. Original (Crossed)     | %13.4f | %9.1f | %9.1f \n', LL_C, AIC_C, BIC_C);
fprintf('4. Nested (Anim+Sess)     | %13.4f | %9.1f | %9.1f \n', LL_N, AIC_N, BIC_N);
fprintf('=======================================================================\n');

% 6. MANUAL LIKELIHOOD RATIO TEST (Original vs Nested)
% We can test these directly against each other because they have the same 
% number of parameters (DF = 2 variance components each). If they have the
% same parameters, we just pick the one with the higher Log-Likelihood
if LL_N > LL_C
    fprintf('The NESTED model fits the data better than the ORIGINAL model.\n');
elseif LL_C > LL_N
    fprintf('The ORIGINAL model fits the data better than the NESTED model.\n');
else
    fprintf('Both models collapsed to the exact same fit.\n');
end
fprintf('=======================================================================\n');