%% export data table for R GAM analysis

% load('V1_HPC_reactivation_coherence_lme1_raw.mat');
cd('P:\corticohippocampal_replay\V1-HPC sleep reactivation')
load('V1_HPC_reactivation_coherence_lme_POST1_raw.mat');
load('V1_HPC_reactivation_coherence_lme_PRE1_raw.mat');

% load('V1_HPC_reactivation_coherence_lme1_raw_100ms.mat');
% load('V1_HPC_reactivation_coherence_lme1_raw_100_200ms.mat');
% load('V1_HPC_reactivation_coherence_lme1_PRE_raw_100_200ms.mat');

% 
% load('V1_HPC_reactivation_coherence_lme2_raw.mat');
% load('V1_HPC_reactivation_coherence_lme2_PRE_raw.mat');
% load('V1_HPC_reactivation_coherence_lme_POST2_raw.mat');
% load('V1_HPC_reactivation_coherence_lme_PRE2_raw.mat');


% load('V1_HPC_reactivation_coherence_lme_-100_100ms_raw.mat');
% load('V1_HPC_reactivation_coherence_lme_-200_200ms_raw.mat');

% load('V1_HPC_reactivation_coherence_lme_0_100ms_raw.mat');
% load('V1_HPC_reactivation_coherence_lme_0_100ms_PRE_raw.mat');
% load('V1_HPC_reactivation_coherence_lme_100_200ms_raw.mat');
% load('V1_HPC_reactivation_coherence_lme_100_200ms_PRE_raw.mat');
% load('V1_HPC_reactivation_coherence_lme_0_200ms_raw.mat');
% load('V1_HPC_reactivation_coherence_lme_0_200ms_PRE_raw.mat');

% load('V1_HPC_reactivation_coherence_lme_POST_0_100ms_raw.mat');
% load('V1_HPC_reactivation_coherence_lme_PRE_0_100ms_raw.mat');
% load('V1_HPC_reactivation_coherence_lme_POST_100_200ms_raw.mat');
% load('V1_HPC_reactivation_coherence_lme_PRE_100_200ms_raw.mat');
% load('V1_HPC_reactivation_coherence_lme_POST_0_200ms_raw.mat');
% load('V1_HPC_reactivation_coherence_lme_PRE_0_200ms_raw.mat');

%% 1. Outlier Removal
% Define the continuous variables that are vulnerable to extreme noise
if ismember('SpindlePowerPRE_Match', tbl.Properties.VariableNames)
    outlierVars = {'RipplePower', 'SpindlePower_Match', 'SpindlePower_NonMatch', ...
        'SpindlePowerPRE_Match', 'SpindlePowerPRE_NonMatch', ...
        'SpindlePowerPOST_Match', 'SpindlePowerPOST_NonMatch', ...
        'HPC_logodds', 'V1_logodds_PRE', 'V1_logodds'};
else
    outlierVars = {'RipplePower', 'SpindlePower_Match', 'SpindlePower_NonMatch', ...
        'HPC_logodds', 'V1_logodds_PRE', 'V1_logodds'};
end

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

if ismember('SpindlePowerPRE_Match', tbl.Properties.VariableNames)
    cleanTbl.SpindlePowerPRE_Match_Z = normalize(cleanTbl.SpindlePowerPRE_Match);
    cleanTbl.SpindlePowerPRE_NonMatch_Z = normalize(cleanTbl.SpindlePowerPRE_NonMatch);
    cleanTbl.SpindlePowerPOST_Match_Z = normalize(cleanTbl.SpindlePowerPOST_Match);
    cleanTbl.SpindlePowerPOST_NonMatch_Z = normalize(cleanTbl.SpindlePowerPOST_NonMatch);
end

cleanTbl.HPC_logodds_Z = cleanTbl.HPC_logodds; % already z-scored
cleanTbl.V1_logodds_PRE_Z = cleanTbl.V1_logodds_PRE; % already z-scored
cleanTbl.V1_logodds_Z = cleanTbl.V1_logodds; % Post-ripple, already z-scored
% cleanTbl.HPC_logodds_Z = zscore(cleanTbl.HPC_logodds); 
% cleanTbl.V1_logodds_PRE_Z = zscore(cleanTbl.V1_logodds_PRE); 
% cleanTbl.V1_logodds_Z = zscore(cleanTbl.V1_logodds); 

%% Calculate coherence

cleanTbl.Event_Coherence_Post = cleanTbl.HPC_logodds_Z .* cleanTbl.V1_logodds_Z;
cleanTbl.Event_Coherence_Pre  = cleanTbl.HPC_logodds_Z .* cleanTbl.V1_logodds_PRE_Z;

% In MATLAB (when you create the CSV)
cleanTbl.Event_Coherence_Post_GeoMean = sign(cleanTbl.HPC_logodds_Z .* cleanTbl.V1_logodds_Z) .* sqrt(abs(cleanTbl.HPC_logodds_Z .* cleanTbl.V1_logodds_Z));
cleanTbl.Event_Coherence_Pre_GeoMean = sign(cleanTbl.HPC_logodds_Z .* cleanTbl.V1_logodds_PRE_Z) .* sqrt(abs(cleanTbl.HPC_logodds_Z .* cleanTbl.V1_logodds_PRE_Z));

% Verify the distribution
fprintf('Mean Event Coherence (Post): %.3f\n', mean(cleanTbl.Event_Coherence_Post));

% cleanTbl.Event_Coherence_Post = cleanTbl.Event_Coherence_Post_GeoMean; % use geometric mean for now
% cleanTbl.Event_Coherence_Pre = cleanTbl.Event_Coherence_Pre_GeoMean; % use geometric mean for now

% writetable(cleanTbl, 'v1_hc_data_geo3.csv'); % export table for R GAM analysis
% writetable(cleanTbl, 'v1_hc_data_geo_POST.csv'); % export table for R GAM analysis
% writetable(cleanTbl, 'v1_hc_data_geo_PRE.csv'); % export table for R GAM analysis
% writetable(cleanTbl, 'v1_hc_data_geo_100ms.csv'); % export table for R GAM analysis
% writetable(cleanTbl, 'v1_hc_data_geo_100_200ms.csv'); % export table for R GAM analysis
% writetable(cleanTbl, 'v1_hc_data_geo_PRE_100_200ms.csv'); % export table for R GAM analysis
% writetable(cleanTbl, 'v1_hc_data_geo_POST.csv'); % export table for R GAM analysis
% writetable(cleanTbl, 'v1_hc_data_geo_PRE.csv'); % export table for R GAM analysis
% 
% writetable(cleanTbl, 'v1_hc_data_geo_POST2.csv'); % POST POST
% writetable(cleanTbl, 'v1_hc_data_geo_PRE2.csv'); % PRE PRE
% writetable(cleanTbl, 'v1_hc_data_geo_POST_SO.csv'); % POST ime3
% writetable(cleanTbl, 'v1_hc_data_geo_PRE_SO.csv'); % PRE ime3
% writetable(cleanTbl, 'v1_hc_data_geo_POST2_SO.csv'); % POST POST3
% writetable(cleanTbl, 'v1_hc_data_geo_PRE2_SO.csv'); % PRE PRE3

% writetable(cleanTbl, 'v1_hc_data_geo_-100_100ms.csv'); % -100_100ms_raw
% writetable(cleanTbl, 'v1_hc_data_geo_-200_200ms.csv'); % -200_200ms_raw

% writetable(cleanTbl, 'v1_hc_data_geo_0_100ms_POST.csv'); % 0 100
% writetable(cleanTbl, 'v1_hc_data_geo_0_100ms_PRE.csv'); % 0 100 PRE
% writetable(cleanTbl, 'v1_hc_data_geo_100_200ms_POST.csv'); % 100 200 
% writetable(cleanTbl, 'v1_hc_data_geo_100_200ms_PRE.csv'); % 100 200 pre
% writetable(cleanTbl, 'v1_hc_data_geo_0_200ms_POST.csv'); % 0 200 
% writetable(cleanTbl, 'v1_hc_data_geo_0_200ms_PRE.csv'); % 0 200 pre

% writetable(cleanTbl, 'v1_hc_data_geo_0_100ms_POST2.csv'); % 0 100
% writetable(cleanTbl, 'v1_hc_data_geo_0_100ms_PRE2.csv'); % 0 100 PRE
% writetable(cleanTbl, 'v1_hc_data_geo_100_200ms_POST2.csv'); % 100 200 
% writetable(cleanTbl, 'v1_hc_data_geo_100_200ms_PRE2.csv'); % 100 200 pre
% writetable(cleanTbl, 'v1_hc_data_geo_0_200ms_POST2.csv'); % 0 200 
% writetable(cleanTbl, 'v1_hc_data_geo_0_200ms_PRE2.csv'); % 0 200 pre

% fig = figure('Name','Geometric Mean reactivation coherence distribution (-100 to 100ms)','Position',[482 111 665 515]);
fig = figure('Name','Geometric Mean reactivation coherence distribution (-200 to 200ms)','Position',[482 111 665 515]);
sgtitle('Reactivation coherence distribution')

subplot(2,2,1)
figure
histogram(cleanTbl.Event_Coherence_Post_GeoMean,'Normalization','probability','EdgeColor','none')
xlim([-2 2])
xlabel('GeoMean Post coherence')
ylabel('Proportion')
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);


subplot(2,2,2)
% histogram(cleanTbl.Event_Coherence_Pre_GeoMean,'Normalization','probability','EdgeColor','none')
% xlim([-2 2])
% xlabel('GeoMean Pre coherence')
% ylabel('Proportion')
% set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);

subplot(2,2,3)
histogram(cleanTbl.Event_Coherence_Post,'Normalization','probability','EdgeColor','none')
xlim([-2 2])
xlabel('Post coherence')
ylabel('Proportion')
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);


subplot(2,2,4)
% histogram(cleanTbl.Event_Coherence_Pre,'Normalization','probability','EdgeColor','none')
% xlim([-2 2])
% xlabel('Pre coherence')
% ylabel('Proportion')
% set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);

% save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','UP-DOWN log odds'),[])

save_all_figures('C:\Users\masah\Documents\GitHub\VR_NPX_analysis\predicting_coherence_with_GAMs\full_model_200ms',[])



fig = figure('Name','Oscillation effect on geomean coherence per session (-100 to 100ms)','Position',[100 100 1000 850]);
fig.Position = [898 165.5000 528 681.5000];
tlo = tiledlayout(3, 2, 'TileSpacing', 'Compact', 'Padding', 'Compact');
title(tlo, 'Reactivation Coherence: Session-Averaged CDF vs Session Pairs', 'FontSize', 14, 'FontWeight', 'Bold');

% Define names and conditions for the loop to parse sequentially
featureNames = {
    'Ripple Power (Low vs High)';
    'Spindle Power Match (Low vs High)';
    'SO Phase Match (Trough vs Peak/Sides)'
    };
featureVars = {
    'RipplePower';
    'SpindlePower_Match';
    'SOPhase_Match'
    };

% Define a shared x-axis vector to evaluate all session CDF curves on
coherence_limit = 2;
x_query = linspace(-coherence_limit, coherence_limit, 200)'; 

for i = 1:3
    % 1. Filter events based on the current feature rule
    if i == 3
        % Special filtering logic for SO Phase Match (circular/bounds)
        low_events  = groupfilter(cleanTbl, "SessionID", @(x) x > -pi/2 & x < pi/2, featureVars{i});
        high_events = groupfilter(cleanTbl, "SessionID", @(x) x < -pi/2 | x > pi/2, featureVars{i});
    else
        % Percentile filtering logic for Ripple and Spindle power
        low_events  = groupfilter(cleanTbl, "SessionID", @(x) x < prctile(x, 30), featureVars{i});
        high_events = groupfilter(cleanTbl, "SessionID", @(x) x > prctile(x, 70), featureVars{i});
    end
    
    % 2. Aggregate to Session Mean (groupsummary automatically creates a 'GroupCount' column)
    low_sess  = groupsummary(low_events, "SessionID", "mean", "Event_Coherence_Post_GeoMean");
    high_sess = groupsummary(high_events, "SessionID", "mean", "Event_Coherence_Post_GeoMean");
    
    % Filter each session table independently using its OWN GroupCount
    low_sess  = low_sess(low_sess.GroupCount >= 100, :);
    high_sess = high_sess(high_sess.GroupCount >= 100, :);
    % 
    % 3. Find intersection to ensure paired sessions only
    [commonSessions, idxLow, idxHigh] = intersect(low_sess.SessionID, high_sess.SessionID);
    low_paired_means  = low_sess.mean_Event_Coherence_Post_GeoMean(idxLow);
    high_paired_means = high_sess.mean_Event_Coherence_Post_GeoMean(idxHigh);
    
    % 4. Perform Paired Wilcoxon Signed-Rank Test
    if ~isempty(commonSessions)
        [p, h, stats] = signrank(high_paired_means, low_paired_means, 'tail', 'right');
        p_str = sprintf('%.4f', p);
    else
        p = NaN;
        p_str = 'N/A';
    end
    
    % --- Plot 1: Session-Averaged CDF Line with Shaded Error Bands (±1 SE) ---
    nexttile
    
    nSessions = length(commonSessions);
    if nSessions > 0
        % Matrices to store the interpolated CDF values for each session
        cdfs_high = zeros(length(x_query), nSessions);
        cdfs_low  = zeros(length(x_query), nSessions);
        
        for s = 1:nSessions
            sessID = commonSessions(s);
            
            % Get event data for this single session
            h_vals = high_events.Event_Coherence_Post_GeoMean(high_events.SessionID == sessID);
            l_vals = low_events.Event_Coherence_Post_GeoMean(low_events.SessionID == sessID);
            
            % Bound data to [-1, 1]
            h_vals = h_vals(h_vals >= -coherence_limit & h_vals <= coherence_limit);
            l_vals = l_vals(l_vals >= -coherence_limit & l_vals <= coherence_limit);
            
            % Calculate individual empirical CDF for High condition
            [f_h, x_h] = ecdf(h_vals);
            [x_h_u, idx_u] = unique(x_h);
            cdfs_high(:, s) = interp1(x_h_u, f_h(idx_u), x_query, 'linear', 'extrap');
            
            % Calculate individual empirical CDF for Low condition
            [f_l, x_l] = ecdf(l_vals);
            [x_l_u, idx_u] = unique(x_l);
            cdfs_low(:, s) = interp1(x_l_u, f_l(idx_u), x_query, 'linear', 'extrap');
        end
        
        % Clean up any extrapolation artifacts forcing bounds between 0 and 1
        cdfs_high = max(0, min(1, cdfs_high));
        cdfs_low  = max(0, min(1, cdfs_low));
        
        % Calculate mean curves
        mean_cdf_high = mean(cdfs_high, 2);
        mean_cdf_low  = mean(cdfs_low, 2);
        
        % Calculate Standard Error (SE = STD / sqrt(N))
        se_cdf_high = std(cdfs_high, 0, 2) / sqrt(nSessions);
        se_cdf_low  = std(cdfs_low, 0, 2) / sqrt(nSessions);
        
        % Plot Shaded Error Band for High condition
        x_conf = [x_query; flipud(x_query)];
        y_conf_high = [mean_cdf_high + se_cdf_high; flipud(mean_cdf_high - se_cdf_high)];
        y_conf_high = max(0, min(1, y_conf_high)); % clamp error bands between 0 and 1
        fill(x_conf, y_conf_high, [0 0.4470 0.7410], 'FaceAlpha', 0.15, 'EdgeColor', 'none');
        hold on;
        
        % Plot Shaded Error Band for Low condition
        y_conf_low = [mean_cdf_low + se_cdf_low; flipud(mean_cdf_low - se_cdf_low)];
        y_conf_low = max(0, min(1, y_conf_low)); 
        fill(x_conf, y_conf_low, [0.8500 0.3250 0.0980], 'FaceAlpha', 0.15, 'EdgeColor', 'none');
        
        % Plot solid mean lines over the shaded regions
        p_high = plot(x_query, mean_cdf_high, 'Color', [0 0.4470 0.7410], 'LineWidth', 1);
        p_low  = plot(x_query, mean_cdf_low, 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1);
    end
    
    xlim([-1 1])
    ylim([0 1])
    xlabel('GeoMean Post coherence')
    ylabel('Mean Probability (CDF)')
    title([featureNames{i} ' - Events'])
    
    % Connect legend exclusively to the lines, skipping the shaded blocks
    if nSessions > 0
        legend([p_high, p_low], {'High','Low'}, 'Location', 'southeast', 'Box', 'off');
    end
    set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 10);
    hold off;
    
    % --- Plot 2: Paired Scatter Plot (Session Level) ---
    nexttile
    if nSessions > 0
        % Calculate a unique horizontal jitter value per session
        jitter = (rand(nSessions, 1) - 0.5) * 0.15;
        
        % Shift both the scatter points AND the line ends by this exact jitter amount
        x_low_jittered  = ones(nSessions, 1) + jitter;
        x_high_jittered = zeros(nSessions, 1) + 2 + jitter;
        
        % Draw connecting lines using the jittered X coordinates
        plot([x_low_jittered, x_high_jittered]', [low_paired_means, high_paired_means]', 'Color', [0.7 0.7 0.7 0.6], 'LineWidth', 1);
        hold on;
        
        % Overlay the scatter points precisely at those same jittered locations
        scatter(x_low_jittered,  low_paired_means,  50, [0.8500 0.3250 0.0980], 'filled', 'MarkerFaceAlpha', 0.8);
        scatter(x_high_jittered, high_paired_means, 50, [0 0.4470 0.7410], 'filled', 'MarkerFaceAlpha', 0.8);
    end
    xlim([0.5 2.5])
    xticks([1 2])
    xticklabels({'Low Condition', 'High Condition'})
    ylabel('Mean GeoMean Post coherence')
    set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 10);
    title(sprintf('Session Means (p=%s)', p_str));
    hold off;

    nSessions
end

save_all_figures('C:\Users\masah\Documents\GitHub\VR_NPX_analysis\predicting_coherence_with_GAMs\full_model',[])
