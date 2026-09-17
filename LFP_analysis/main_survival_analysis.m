% MAIN_SURVIVAL_ANALYSIS
% Counting-process Cox Proportional Hazards Survival Analysis for Cortical UP State Termination
% Evaluates individual ripple LFP/MUA, cumulative ripple activity, and baseline MUA.
clear all
addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\OneDrive\Documents\GitHub\VR_NPX_analysis'))

if exist('C:\Users\masah\OneDrive\Documents\corticohippocampal_replay', 'dir')
    analysis_folder = 'C:\Users\masah\OneDrive\Documents\corticohippocampal_replay';
elseif exist('D:\corticohippocampal_replay', 'dir')
    analysis_folder = 'D:\corticohippocampal_replay';
elseif exist('P:\corticohippocampal_replay', 'dir')
    analysis_folder = 'P:\corticohippocampal_replay';
else
    analysis_folder = pwd;
end

output_dir = fullfile(analysis_folder, 'V1-HPC bilateral interaction', 'survival analysis counting process');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% 1. Load Data
load(fullfile(analysis_folder, 'ripples_all_POST.mat'));
load(fullfile(analysis_folder, 'V1-HPC sleep interaction', 'SO_ripples_probability_whole_combined.mat'));
probability_psth_whole = probability;

load(fullfile(analysis_folder, 'slow_waves_all_POST.mat'));
load(fullfile(analysis_folder, 'V1-HPC sleep reactivation', 'UP_DOWN_info_100ms.mat'), 'UP_DOWN_info');
load(fullfile(analysis_folder, 'V1-HPC sleep reactivation', 'ripple_info.mat'), 'ripple_info');
load(fullfile(analysis_folder, 'V1-HPC sleep interaction', 'merged_UP_DOWN_ripples_event_info.mat'), 'merged_event_info');
load(fullfile(analysis_folder,'V1-HPC sleep interaction','UP_DOWN_ripples_event_info.mat'),'event_info');

%% 2. Process Session Timestamps and Merge Bilateral Ripples & Spike Times
UP_ints = [];
DOWN_ints = [];
ripple_peaktimes = [];
ripple_ints = [];
SO_ints = [];
V1_MUA_spiketimes = [];
HC_MUA_spiketimes = [];

sessions_to_process = 1:max(slow_waves_all(1).UP_session_count);

for nprobe = 1:2
    V1_MUA_spiketimes{nprobe} = [];
    HC_MUA_spiketimes{nprobe} = [];

    UP_ints{nprobe}    = slow_waves_all(nprobe).UP_ints;
    DOWN_ints{nprobe}  = slow_waves_all(nprobe).DOWN_ints;
    SO_ints{nprobe}    = slow_waves_all(nprobe).DOWN_intervals;
    ripple_peaktimes{nprobe} = ripples_all(nprobe).peaktimes(ripples_all(nprobe).SWS_index == 1);
    ripple_ints{nprobe}      = [ripples_all(nprobe).onset(ripples_all(nprobe).SWS_index == 1), ...
                                ripples_all(nprobe).offset(ripples_all(nprobe).SWS_index == 1)];

    for nsession = 1:length(sessions_to_process)
        sess_val = sessions_to_process(nsession);
        
        index = find(slow_waves_all(nprobe).DOWN_intervals_session == sess_val);
        SO_ints{nprobe}(index, :) = SO_ints{nprobe}(index, :) + nsession * 1000000;

        index = find(slow_waves_all(nprobe).UP_session_count == sess_val);
        UP_ints{nprobe}(index, :) = UP_ints{nprobe}(index, :) + nsession * 1000000;

        index = find(slow_waves_all(nprobe).DOWN_session_count == sess_val);
        DOWN_ints{nprobe}(index, :) = DOWN_ints{nprobe}(index, :) + nsession * 1000000;

        index = find(ripples_all(nprobe).session_count(ripples_all(nprobe).SWS_index == 1) == sess_val);
        ripple_ints{nprobe}(index, :) = ripple_ints{nprobe}(index, :) + nsession * 1000000;
        ripple_peaktimes{nprobe}(index, :) = ripple_peaktimes{nprobe}(index, :) + nsession * 1000000;

        if nsession <= length(slow_waves_all(nprobe).V1_MUA_spiketimes) && ~isempty(slow_waves_all(nprobe).V1_MUA_spiketimes{nsession})
            V1_MUA_spiketimes{nprobe} = [V1_MUA_spiketimes{nprobe}; slow_waves_all(nprobe).V1_MUA_spiketimes{nsession} + nsession * 1000000];
            HC_MUA_spiketimes{nprobe} = [HC_MUA_spiketimes{nprobe}; slow_waves_all(nprobe).HPC_MUA_spiketimes{nsession} + nsession * 1000000];
        end
    end
end

% Assemble merged events across probes
merged_event_info.UP_ints   = [UP_ints{1}(probability_psth_whole(1).UP_all_index, :); ...
                                UP_ints{2}(probability_psth_whole(2).UP_all_index, :)];
merged_event_info.DOWN_ints = [DOWN_ints{1}(probability_psth_whole(1).DOWN_all_index, :); ...
                                DOWN_ints{2}(probability_psth_whole(2).DOWN_all_index, :)];

merged_event_info.UP_hemisphere_id   = [ones(length(probability_psth_whole(1).UP_all_index), 1); ...
                                        2 * ones(length(probability_psth_whole(2).UP_all_index), 1)];
merged_event_info.DOWN_hemisphere_id = [ones(length(probability_psth_whole(1).DOWN_all_index), 1); ...
                                        2 * ones(length(probability_psth_whole(2).DOWN_all_index), 1)];

merged_event_info.ripples_peaktimes     = [ripple_peaktimes{1}; ripple_peaktimes{2}];
merged_event_info.ripples_ints          = [ripple_ints{1}; ripple_ints{2}];
merged_event_info.ripples_hemisphere_id = [ones(length(ripple_peaktimes{1}), 1); ...
                                            2 * ones(length(ripple_peaktimes{2}), 1)];

[event_ids_first, event_ids_second] = merge_bilateral_ripple_events(merged_event_info.ripples_hemisphere_id, ...
    merged_event_info.ripples_peaktimes, 0.05);

merged_event_info.ripples_hemisphere_id = merged_event_info.ripples_hemisphere_id(event_ids_first);
merged_event_info.ripples_peaktimes     = merged_event_info.ripples_peaktimes(event_ids_first, :);
merged_event_info.ripples_ints          = merged_event_info.ripples_ints(event_ids_first, :);

ripples_original_index = [find(ripples_all(1).SWS_index); find(ripples_all(2).SWS_index)];
merged_event_info.ripples_original_index = ripples_original_index(event_ids_first);

ripplePower = [ripples_all(1).peak_zscore(ripples_all(1).SWS_index); ripples_all(2).peak_zscore(ripples_all(2).SWS_index)];
ripplePower = mean([ripplePower(event_ids_first), ripplePower(event_ids_second)], 2);
merged_event_info.ripples_power = ripplePower;

% Session and subject metadata extraction
UP_session_count = [slow_waves_all(1).UP_session_count(probability_psth_whole(1).UP_all_index); ...
                    slow_waves_all(2).UP_session_count(probability_psth_whole(2).UP_all_index)];
subject_id = str2double(cellstr(slow_waves_all(1).subject(UP_session_count, end-1:end)));
[~, ~, subject_id] = unique(subject_id);

merged_event_info.session_id = UP_session_count;
merged_event_info.subject_id = subject_id;

%% 3. Build Counting Process Interval Table directly from raw spike times
fprintf('Building UP state counting process interval table directly from raw spiketimes...\n');
T = build_UP_counting_process_intervals(merged_event_info, V1_MUA_spiketimes, HC_MUA_spiketimes);
% T.upID


save(fullfile(output_dir, 'survival_analysis_interval_table.mat'), 'T');
fprintf('Interval table generated with %d sub-intervals across %d UP events.\n', size(T,1), max(T.upID));


%% Model 1 V1 baseline activity
T_ripples = T(T.numRipplesInUP>0,:);
model1_covariates = {'cumTotalV1_MUA_incl'};
model1_labels     = {'Cumulative V1 activity'};
% model1_covariates = {'start','cumTotalV1_MUA'};
% model1_labels     = {'Time elapsed since UP','Cumulative V1 activity'};


output_model1 = plot_UP_survival_counting_process(T_ripples, model1_covariates, ...
    'title_name', 'Baseline V1 Effect on UP Survival', ...
    'feature_labels', model1_labels, ...
    'strata_var', 'session_id',....
    'bootstrap',false);


%% 4. Model 2: Individual Ripple Effect on UP State Termination
T_ripples = T(T.numRipplesInUP>0,:);
% 
% unique(T_ripples.upID)
fprintf('Running Model 1: Individual Ripple Effect...\n');
% model1_covariates = {'inRipple', 'ripplePower', 'rippleHPC_MUA_sum', 'rippleHPC_MUA_mean'};
% model1_labels     = {'In Ripple', 'Ripple Power', 'Ripple HPC MUA Sum', 'Ripple HPC MUA Rate'};
% model1_covariates = {'inRipple', 'ripplePower', 'rippleHPC_MUA_sum'};
% model1_labels     = {'In Ripple', 'Ripple Power', 'Ripple HPC MUA Sum'};

% model1_covariates = {'rippleHPC_MUA_sum','nonRippleHPC_MUA_sum'};
% model1_labels     = {'Ripple HPC MUA','Non Ripple HPC MUA'};

% model1_covariates = {'rippleHPC_MUA_sum','cumRippleHPC_MUA','cumNonRippleHPC_MUA'};
% model1_labels     = {'Ripple HPC MUA','Past Ripple HPC MUA','Non Ripple HPC MUA'};
model1_covariates = {'cumTotalV1_MUA_incl','lastRipplePower','lastRippleHPC_MUA_mean','lastRippleHPC_MUA_sum'}; % Last here means ripple power during ripple as well as interval after ripple 
model1_labels     = {'Cumulative V1 activity','Ripple power','Ripple HPC MUA rate','Ripple HPC MUA sum'};
% model1_covariates = {'cumTotalV1_MUA_incl'};
% model1_labels     = {'Cumulative V1 activity'};

% model1_covariates = {'start','cumTotalV1_MUA','lastRipplePower','lastRippleHPC_MUA_mean','lastRippleHPC_MUA_sum','timeSinceLastRipple', 'cumRippleHPC_MUA', 'cumNonRippleHPC_MUA'};
% model1_labels     = {'Time elapsed since UP','Cumulative V1 activity','Ripple Power','Ripple HPC MUA rate','Ripple HPC MUA sum','Time since ripple', 'Cum Ripple HPC MUA', 'Cum Non-Ripple HPC MUA'};


output_model = plot_UP_survival_counting_process(T_ripples, model1_covariates, ...
    'title_name', 'Individual Ripple Effect on UP Survival', ...
    'feature_labels', model1_labels, ...
    'strata_var', 'session_id',....
    'is_multivariate', true, ...
    'bootstrap',false);

output_model_univariate = plot_UP_survival_counting_process(T_ripples, model1_covariates, ...
    'title_name', 'Individual Ripple Effect on UP Survival (univariate)', ...
    'feature_labels', model1_labels, ...
    'strata_var', 'session_id',....
    'is_multivariate', false, ...
    'bootstrap',false);

save(fullfile(output_dir, 'survival_model_counting_process_individual_ripples.mat'), 'output_model','output_model_univariate');
save_all_figures(output_dir, []);



model1_covariates = {'cumTotalV1_MUA_incl','cumRippleHPC_MUA', 'cumNonRippleHPC_MUA'};
model1_labels     = {'Cumulative V1 activity','Cum Ripple HPC MUA','Cum Non-Ripple HPC MUA'};
% model1_covariates = {'start','cumTotalV1_MUA','lastRipplePower','lastRippleHPC_MUA_mean','lastRippleHPC_MUA_sum','timeSinceLastRipple', 'cumRippleHPC_MUA', 'cumNonRippleHPC_MUA'};
% model1_labels     = {'Time elapsed since UP','Cumulative V1 activity','Ripple Power','Ripple HPC MUA rate','Ripple HPC MUA sum','Time since ripple', 'Cum Ripple HPC MUA', 'Cum Non-Ripple HPC MUA'};

% output_model1 = plot_UP_survival_counting_process(T_ripples, model1_covariates, ...
%     'title_name', 'past ripple and past non-ripple history on UP Survival (univariate)', ...
%     'feature_labels', model1_labels, ...
%     'strata_var', 'session_id',....
%     'is_multivariate', false, ...
%     'bootstrap',false);

output_model2 = plot_UP_survival_counting_process(T_ripples, model1_covariates, ...
    'title_name', 'past ripple and past non-ripple history on UP Survival', ...
    'feature_labels', model1_labels, ...
    'strata_var', 'session_id',....
    'bootstrap',false,...
    'plot_survival',true);




%% 4. Model 2: All effects
% T_ripples = T;
T_ripples = T(T.numRipplesInUP>0,:);
fprintf('Running Model 1: Individual Ripple Effect...\n');
% model1_covariates = {'inRipple', 'ripplePower', 'rippleHPC_MUA_sum', 'rippleHPC_MUA_mean'};
% model1_labels     = {'In Ripple', 'Ripple Power', 'Ripple HPC MUA Sum', 'Ripple HPC MUA Rate'};
% model1_covariates = {'inRipple', 'ripplePower', 'rippleHPC_MUA_sum'};
% model1_labels     = {'In Ripple', 'Ripple Power', 'Ripple HPC MUA Sum'};

% model1_covariates = {'rippleHPC_MUA_sum','nonRippleHPC_MUA_sum'};
% model1_labels     = {'Ripple HPC MUA','Non Ripple HPC MUA'};

% model1_covariates = {'rippleHPC_MUA_sum','cumRippleHPC_MUA','cumNonRippleHPC_MUA'};
% model1_labels     = {'Ripple HPC MUA','Past Ripple HPC MUA','Non Ripple HPC MUA'};
model1_covariates = {'cumTotalV1_MUA_incl','lastRipplePower','lastRippleHPC_MUA_mean','cumRippleHPC_MUA', 'cumNonRippleHPC_MUA_incl'};
model1_labels     = {'Cumulative V1 activity','Ripple Power','Ripple HPC MUA rate','Cum Ripple HPC MUA', 'Cum Non-Ripple HPC MUA'};
% model1_covariates = {'start','cumTotalV1_MUA','lastRipplePower','lastRippleHPC_MUA_mean','lastRippleHPC_MUA_sum','timeSinceLastRipple', 'cumRippleHPC_MUA', 'cumNonRippleHPC_MUA'};
% model1_labels     = {'Time elapsed since UP','Cumulative V1 activity','Ripple Power','Ripple HPC MUA rate','Ripple HPC MUA sum','Time since ripple', 'Cum Ripple HPC MUA', 'Cum Non-Ripple HPC MUA'};


output_model = plot_UP_survival_counting_process(T_ripples, model1_covariates, ...
    'title_name', 'Ripple Effect on UP Survival', ...
    'feature_labels', model1_labels, ...
    'strata_var', 'session_id',....
    'bootstrap',false);

output_model_univariate = plot_UP_survival_counting_process(T_ripples, model1_covariates, ...
    'title_name', 'Ripple Effect on UP Survival (univariate)', ...
    'feature_labels', model1_labels, ...
    'strata_var', 'session_id',....
    'bootstrap',false,...
    'is_multivariate', false, ...
    'plot_survival',false);

save(fullfile(output_dir, 'survival_model_counting_process_all.mat'), 'output_model','output_model_univariate');
save_all_figures(output_dir, []);


%% 4. Model 1: Individual Ripple Effect on UP State Termination
fprintf('Running Model 1: Individual Ripple Effect...\n');
% model1_covariates = {'inRipple', 'ripplePower', 'rippleHPC_MUA_sum', 'rippleHPC_MUA_mean'};
% model1_labels     = {'In Ripple', 'Ripple Power', 'Ripple HPC MUA Sum', 'Ripple HPC MUA Rate'};
% model1_covariates = {'inRipple', 'ripplePower', 'rippleHPC_MUA_sum'};
% model1_labels     = {'In Ripple', 'Ripple Power', 'Ripple HPC MUA Sum'};

% model1_covariates = {'rippleHPC_MUA_sum','nonRippleHPC_MUA_sum'};
% model1_labels     = {'Ripple HPC MUA','Non Ripple HPC MUA'};

% model1_covariates = {'rippleHPC_MUA_sum','cumRippleHPC_MUA','cumNonRippleHPC_MUA'};
% model1_labels     = {'Ripple HPC MUA','Past Ripple HPC MUA','Non Ripple HPC MUA'};
model1_covariates = {'cumTotalV1_MUA','lastRipplePower','lastRippleHPC_MUA_mean','lastRippleHPC_MUA_sum','timeSinceLastRipple', 'cumRippleHPC_MUA', 'cumNonRippleHPC_MUA'};
model1_labels     = {'Cumulative V1 activity','Ripple Power','Ripple HPC MUA rate','Ripple HPC MUA sum','Time since ripple', 'Cum Ripple HPC MUA', 'Cum Non-Ripple HPC MUA'};
% model1_covariates = {'start','cumTotalV1_MUA','lastRipplePower','lastRippleHPC_MUA_mean','lastRippleHPC_MUA_sum','timeSinceLastRipple', 'cumRippleHPC_MUA', 'cumNonRippleHPC_MUA'};
% model1_labels     = {'Time elapsed since UP','Cumulative V1 activity','Ripple Power','Ripple HPC MUA rate','Ripple HPC MUA sum','Time since ripple', 'Cum Ripple HPC MUA', 'Cum Non-Ripple HPC MUA'};


output_model1 = plot_UP_survival_counting_process(T, model1_covariates, ...
    'title_name', 'Individual Ripple Effect on UP Survival', ...
    'feature_labels', model1_labels, ...
    'strata_var', 'session_id',....
    'bootstrap',false);

% output_model1 = plot_UP_survival_counting_process(T, model1_covariates, ...
%     'title_name', 'Individual Ripple Effect on UP Survival', ...
%     'feature_labels', model1_labels, ...
%     'strata_var', 'session_id', ...
%     'stratify_feature', 'rippleHPC_MUA_sum');

save(fullfile(output_dir, 'model1_individual_ripple_survival.mat'), 'output_model1');

%% 5. Model 2: Cumulative Ripple Effect on UP State Termination
fprintf('Running Model 2: Cumulative Ripple Effect...\n');
model2_covariates = {'rippleHPC_MUA_sum','cumRippleHPC_MUA', 'cumRippleCount'};
model2_labels     = {'Cum Ripple HPC MUA', 'Cum Ripple Count'};

output_model2 = plot_UP_survival_counting_process(T, model2_covariates, ...
    'title_name', 'Cumulative Ripple Effect on UP Survival', ...
    'feature_labels', model2_labels, ...
    'strata_var', 'session_id', ...
    'stratify_feature', 'cumRippleHPC_MUA');

save(fullfile(output_dir, 'model2_cumulative_ripple_survival.mat'), 'output_model2');

%% 6. Model 3: Baseline & Non-Ripple Cumulative MUA Effect
fprintf('Running Model 3: Baseline & Non-Ripple MUA Effect...\n');
model3_covariates = {'cumNonRippleHPC_MUA', 'cumTotalHPC_MUA', 'cumTotalV1_MUA'};
model3_labels     = {'Cum Non-Ripple HPC MUA', 'Cum Total HPC MUA', 'Cum Ipsi V1 MUA'};

output_model3 = plot_UP_survival_counting_process(T, model3_covariates, ...
    'title_name', 'Baseline MUA Activity Effect on UP Survival', ...
    'feature_labels', model3_labels, ...
    'strata_var', 'session_id', ...
    'stratify_feature', 'cumTotalHPC_MUA');

save(fullfile(output_dir, 'model3_baseline_mua_survival.mat'), 'output_model3');

%% 7. Model 4: Full Multivariable Combined Model
fprintf('Running Model 4: Full Multivariable Combined Model...\n');
model4_covariates = {'inRipple', 'ripplePower', 'rippleHPC_MUA_sum', 'cumRippleHPC_MUA', 'cumNonRippleHPC_MUA', 'cumTotalV1_MUA'};
model4_labels     = {'In Ripple', 'Ripple Power', 'Ripple HPC MUA', 'Cum Ripple HPC MUA', 'Cum Non-Ripple HPC MUA', 'Cum Ipsi V1 MUA'};

output_model4 = plot_UP_survival_counting_process(T, model4_covariates, ...
    'title_name', 'Full Multivariable UP Survival Model', ...
    'feature_labels', model4_labels, ...
    'strata_var', 'session_id', ...
    'stratify_feature', 'cumRippleHPC_MUA');

save(fullfile(output_dir, 'model4_full_multivariable_survival.mat'), 'output_model4');

%% 8. Model 5: Carried-Forward Ripple & Post-Ripple Elapsing Time Model
fprintf('Running Model 5: Carried-Forward Ripple & Post-Ripple Elapsing Time Model...\n');
model5_covariates = {'lastRipplePower', 'lastRippleHPC_MUA_sum', 'timeSinceLastRipple', 'cumRippleHPC_MUA', 'cumNonRippleHPC_MUA'};
model5_labels     = {'Last Ripple Power', 'Last Ripple HPC MUA', 'Time Since Last Ripple', 'Cum Ripple HPC MUA', 'Cum Non-Ripple HPC MUA'};

output_model5 = plot_UP_survival_counting_process(T, model5_covariates, ...
    'title_name', 'Carried-Forward Ripple and Post-Ripple Delay Survival Model', ...
    'feature_labels', model5_labels, ...
    'strata_var', 'session_id', ...
    'stratify_feature', 'lastRipplePower');

save(fullfile(output_dir, 'model5_carried_forward_ripple_survival.mat'), 'output_model5');

%% 9. Save Figures
if exist('save_all_figures', 'file')
    save_all_figures(output_dir, []);
end

fprintf('Survival analysis counting process pipeline completed successfully.\n');
