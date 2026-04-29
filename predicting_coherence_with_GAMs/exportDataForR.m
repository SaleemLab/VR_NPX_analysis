%% export data table for R GAM analysis

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

cleanTbl.Event_Coherence_Post = cleanTbl.HPC_logodds_Z .* cleanTbl.V1_logodds_Z;
cleanTbl.Event_Coherence_Pre  = cleanTbl.HPC_logodds_Z .* cleanTbl.V1_logodds_PRE_Z;

% In MATLAB (when you create the CSV)
cleanTbl.Event_Coherence_Post_GeoMean = sign(cleanTbl.HPC_logodds_Z .* cleanTbl.V1_logodds_Z) .* sqrt(abs(cleanTbl.HPC_logodds_Z .* cleanTbl.V1_logodds_Z));
cleanTbl.Event_Coherence_Pre_GeoMean = sign(cleanTbl.HPC_logodds_Z .* cleanTbl.V1_logodds_PRE_Z) .* sqrt(abs(cleanTbl.HPC_logodds_Z .* cleanTbl.V1_logodds_PRE_Z));

% Verify the distribution
fprintf('Mean Event Coherence (Post): %.3f\n', mean(cleanTbl.Event_Coherence_Post));

cleanTbl.Event_Coherence_Post = cleanTbl.Event_Coherence_Post_GeoMean; % use geometric mean for now
cleanTbl.Event_Coherence_Pre = cleanTbl.Event_Coherence_Post_GeoMean; % use geometric mean for now

writetable(cleanTbl, 'v1_hc_data_geo3.csv'); % export table for R GAM analysis