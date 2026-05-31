%% Corticohippocampal reactivation during ripples

addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\OneDrive\Documents\GitHub\VR_NPX_analysis'))


if exist('C:\Users\masah\OneDrive\Documents\corticohippocampal_replay')
    analysis_folder = 'C:\Users\masah\OneDrive\Documents\corticohippocampal_replay';
elseif exist('D:\corticohippocampal_replay')>0
    analysis_folder = 'D:\corticohippocampal_replay';
elseif exist('P:\corticohippocampal_replay')>0
    analysis_folder = 'P:\corticohippocampal_replay';
end
load(fullfile(analysis_folder,'slow_waves_all_POST.mat'))
% load(fullfile(analysis_folder,'slow_waves_all_markov_POST.mat'))
load(fullfile(analysis_folder,'ripples_all_POST.mat'))
load(fullfile(analysis_folder,'spindles_all_POST.mat'))
load(fullfile(analysis_folder,'behavioural_state_merged_all_POST.mat'))
load(fullfile(analysis_folder,'V1-HPC sleep interaction','UP_DOWN_ripples_event_info.mat'),'event_info');
load(fullfile(analysis_folder,'V1-HPC sleep interaction','merged_UP_DOWN_ripples_event_info.mat'),'merged_event_info');

load(fullfile(analysis_folder,'V1-HPC sleep interaction','UP_DOWN_ripple_PSTH_MUA.mat'),'UP_DOWN_ripple_PSTH_MUA');
PSTH_MUA = UP_DOWN_ripple_PSTH_MUA;
load(fullfile(analysis_folder,'V1-HPC sleep interaction','UP_DOWN_ripple_PSTH_MUA_baseline.mat'));
PSTH_MUA_baseline = UP_DOWN_ripple_PSTH_MUA;

load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_whole.mat'));
probability_psth_whole = probability;

load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_whole.mat'));

load(fullfile(analysis_folder,'bayesian_reactivation_V1_all_POST.mat'))
load(fullfile(analysis_folder,'bayesian_reactivation_all_POST.mat'))

sessions_to_process = 1:max(slow_waves_all(1).UP_session_count);



load(fullfile(analysis_folder,'periripple_LFP_info_V1.mat'));

%% Extract key information for DOWN UP, ripples and spindles
 
%%%%%%% Find reference channel/shank
cortex_ref_shank = [];
HPC_ref_shank = [];

for nsession = 1:max(ripples_all(1).session_count)
    for probe_no = 1:2
        cortex_ref_shank(nsession,probe_no) = find(slow_waves_all(probe_no).shank_id{nsession} == slow_waves_all(probe_no).shank{nsession}(slow_waves_all(probe_no).channel{nsession} == slow_waves_all(probe_no).best_channel(nsession))...
            &slow_waves_all(probe_no).probe_hemisphere{nsession} == probe_no);
        % [~,idx] = min(abs(ripples_all(probe_no).SWR_peaktimes{nsession}' - ripples_all(probe_no).peaktimes(ripples_all(probe_no).session_count==nsession))');
        % ripple_counts = histcounts(idx,length(ripples_all(probe_no).shank_id{nsession}));
        % [~,HPC_ref_shank(nsession,probe_no)] = max(ripple_counts);

        shank_id = find(ripples_all(probe_no).probe_hemisphere{nsession} == probe_no);
        HPC_ref_shank(nsession,probe_no) = shank_id(ripples_all(probe_no).best_channel(nsession));

    end
end

%%%%%%%%%%%%%%%%%% Ripple info


% if merge events to bilateral events as one event
[event_ids_first,event_ids_second] = merge_bilateral_ripple_events(merged_event_info.ripples_hemisphere_id,merged_event_info.ripples_peaktimes,0.05);


% load(fullfile(analysis_folder,'V1-HPC sleep interaction','ripples_spindles_probability_whole.mat'));
session_count = [ripples_all(1).session_count(ripples_all(1).SWS_index==1); ripples_all(2).session_count(ripples_all(2).SWS_index==1)];
subject_id = str2double(cellstr(ripples_all(1).subject(session_count,end-1:end)));
[~, ~, subject_id] = unique(subject_id);

%%% ripple power
ripple_info.ripple_power = [ripples_all(1).peak_zscore(ripples_all(1).SWS_index==1); ripples_all(2).peak_zscore(ripples_all(2).SWS_index==1)];
ripple_info.ripple_power = mean([ripple_info.ripple_power(event_ids_first) ripple_info.ripple_power(event_ids_second)],2);


%%% spindle co-occurance
[~,spindle_index,~,index] =RestrictInts(merged_event_info.ripples_ints,merged_event_info.spindles_ints);
ripple_info.spindle_presence = spindle_index;
ripple_info.spindle_presence_hemi = zeros(size(spindle_index));
ripple_info.spindle_presence_hemi(find(spindle_index)) = merged_event_info.spindles_hemisphere_id(index);

ripple_info.spindle_presence = ripple_info.spindle_presence(event_ids_first);
c.spindle_presence_hemi = ripple_info.spindle_presence_hemi(event_ids_first);

%%% spindle phase
spindle_phase=[];
for probe_no = 1:2
    spindle_phase{probe_no}=[];
    for nsession = 1:length(sessions_to_process)
        spindle_phase{probe_no} = [spindle_phase{probe_no} ripples_all(probe_no).spindle_phase_ripple_peaktime{nsession}(cortex_ref_shank(nsession,:),:)];
    end
end

spindle_phase = [spindle_phase{1}(:,ripples_all(1).SWS_index==1) spindle_phase{2}(:,ripples_all(2).SWS_index==1)];

ripple_info.spindle_phase = spindle_phase';
ripple_info.spindle_phase = ripple_info.spindle_phase(event_ids_first,:);

%%% spindle power
spindle_amplitude=[];
for probe_no = 1:2
    spindle_amplitude{probe_no}=[];
    for nsession = 1:length(sessions_to_process)
        spindle_amplitude{probe_no} = [spindle_amplitude{probe_no} ripples_all(probe_no).spindle_amplitude_ripple_peaktime{nsession}(cortex_ref_shank(nsession,:),:)];
    end
end

spindle_amplitude = [spindle_amplitude{1}(:,ripples_all(1).SWS_index==1) spindle_amplitude{2}(:,ripples_all(2).SWS_index==1)];

ripple_info.spindle_amplitude = spindle_amplitude';
ripple_info.spindle_amplitude = ripple_info.spindle_amplitude(event_ids_first,:);

%%% SO phase
SO_phase=[];
for probe_no = 1:2
    SO_phase{probe_no}=[];
    for nsession = 1:length(sessions_to_process)
        SO_phase{probe_no} = [SO_phase{probe_no} ripples_all(probe_no).SO_phase_ripple_peaktime{nsession}(cortex_ref_shank(nsession,:),:)];
    end
end

SO_phase = [SO_phase{1}(:,ripples_all(1).SWS_index==1) SO_phase{2}(:,ripples_all(2).SWS_index==1)];

ripple_info.SO_phase = SO_phase';
ripple_info.SO_phase = ripple_info.SO_phase(event_ids_first,:);

%%% SO power
SO_amplitude=[];
for probe_no = 1:2
    SO_amplitude{probe_no}=[];
    for nsession = 1:length(sessions_to_process)
        SO_amplitude{probe_no} = [SO_amplitude{probe_no} ripples_all(probe_no).SO_amplitude_ripple_peaktime{nsession}(cortex_ref_shank(nsession,:),:)];
    end
end

SO_amplitude = [SO_amplitude{1}(:,ripples_all(1).SWS_index==1) SO_amplitude{2}(:,ripples_all(2).SWS_index==1)];

ripple_info.SO_amplitude = SO_amplitude';
ripple_info.SO_amplitude = ripple_info.SO_amplitude(event_ids_first,:);



%%% early UP transition co-occurance
[~,UP_index,~,index] = RestrictInts(merged_event_info.ripples_ints,[merged_event_info.UP_ints(:,1) merged_event_info.UP_ints(:,1)+0.2]);
ripple_info.early_UP_index = UP_index;
ripple_info.early_UP_index_hemi = zeros(size(UP_index));
% ripple_info.spindle_presence_hemi = zeros(size(spindle_index));
ripple_info.early_UP_index_hemi(find(UP_index)) = merged_event_info.UP_hemisphere_id(index);

ripple_info.early_UP_index = ripple_info.early_UP_index(event_ids_first);
ripple_info.early_UP_index_hemi = ripple_info.early_UP_index_hemi(event_ids_first);

%%% late UP transition co-occurance
[~,UP_index,~,index] = RestrictInts(merged_event_info.ripples_ints,[merged_event_info.DOWN_ints(:,1)-0.2 merged_event_info.DOWN_ints(:,1)]);
ripple_info.late_UP_index = UP_index;
ripple_info.late_UP_index_hemi = zeros(size(UP_index));
% ripple_info.spindle_presence_hemi = zeros(size(spindle_index));
ripple_info.late_UP_index_hemi(find(UP_index)) = merged_event_info.DOWN_hemisphere_id(index);

ripple_info.late_UP_index = ripple_info.late_UP_index(event_ids_first);
ripple_info.late_UP_index_hemi = ripple_info.late_UP_index_hemi(event_ids_first);

%%% UP transition co-occurance
[~,UP_index,~,index] = RestrictInts(merged_event_info.ripples_ints,[merged_event_info.UP_ints(:,1) merged_event_info.UP_ints(:,2)]);
ripple_info.UP_index = UP_index;
ripple_info.UP_index_hemi = zeros(size(UP_index));
% ripple_info.spindle_presence_hemi = zeros(size(spindle_index));
ripple_info.UP_index_hemi(find(UP_index)) = merged_event_info.UP_hemisphere_id(index);

ripple_info.UP_index = ripple_info.UP_index(event_ids_first);
ripple_info.UP_index_hemi = ripple_info.UP_index_hemi(event_ids_first);
%%% DOWN transition co-occurance
[~,DOWN_index,~,index] = RestrictInts(merged_event_info.ripples_ints,[merged_event_info.DOWN_ints(:,1) merged_event_info.DOWN_ints(:,2)]);
ripple_info.DOWN_index = DOWN_index;
ripple_info.DOWN_index_hemi = zeros(size(DOWN_index));
% ripple_info.spindle_presence_hemi = zeros(size(spindle_index));
ripple_info.DOWN_index_hemi(find(DOWN_index)) = merged_event_info.DOWN_hemisphere_id(index);

ripple_info.DOWN_index = ripple_info.DOWN_index(event_ids_first);
ripple_info.DOWN_index_hemi = ripple_info.DOWN_index_hemi(event_ids_first);


%%%%%%
session_count = [ripples_all(1).session_count(ripples_all(1).SWS_index==1); ripples_all(2).session_count(ripples_all(2).SWS_index==1)];
subject_id = str2double(cellstr(ripples_all(1).subject(session_count,end-1:end)));
[~, ~, subject_id] = unique(subject_id);


% singlet_index = logical(([1; diff(merged_event_info.ripples_peaktimes)>0.1]));
singlet_index = logical(ones(length(merged_event_info.ripples_peaktimes),1));



%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%
%%%%%%%%%%%%
%%%%%%%%%%%%
%%%%%%%%%%%%%%% KDE reactivation bias 
load(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_reactivation_ripples_PSTH.mat'))
% load(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_reactivation_content.mat'))

timebin = 0.01;
time_windows = [-1 1];
% Generate bin edges
bin_edges = time_windows(1):timebin:time_windows(2);
% Generate bin centers
bin_centers = bin_edges(1:end-1) + timebin/2;


% z_bias1 = KDE_reactivation_content.HPC_z_ripples';
% T1_events = find(nanmean(z_bias1(1:10,:))>0.5);
% T2_events = find(nanmean(z_bias1(1:10,:))<-0.5);

% z_bias = KDE_reactivation_ripples_PSTH.HPC_z_ripples';
% z_bias_V1 = KDE_reactivation_ripples_PSTH.V1_z_ripples';

z_bias = KDE_reactivation_ripples_PSTH.HPC_z_logodds_ripples' + KDE_reactivation_ripples_PSTH.nan_mask';
z_bias_V1 = KDE_reactivation_ripples_PSTH.V1_z_logodds_ripples' + KDE_reactivation_ripples_PSTH.nan_mask';

z_bias1 = z_bias(isfinite(z_bias));
z_bias(z_bias>=inf) = prctile(z_bias1,99.5);
z_bias(z_bias<=-inf) = prctile(z_bias1,0.5);

z_bias1 = z_bias(isfinite(z_bias_V1));
z_bias_V1(z_bias_V1>=inf) = prctile(z_bias1,99.5);
z_bias_V1(z_bias_V1<=-inf) = prctile(z_bias1,0.5);

% z_bias_V1 = [PSTH_MUA(1).L_V1_ripples-PSTH_MUA(1).R_V1_ripples; PSTH_MUA(2).R_V1_ripples-PSTH_MUA(2).L_V1_ripples]';

z_bias = z_bias(:,event_ids_first);
z_bias_V1 = z_bias_V1(:,event_ids_first);



%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

bins_to_use = bin_centers>0 & bin_centers<0.1;
session_count = [ripples_all(1).session_count(ripples_all(1).SWS_index==1); ripples_all(2).session_count(ripples_all(2).SWS_index==1)];
subject_id = str2double(cellstr(ripples_all(1).subject(session_count,end-1:end)));
[~, ~, subject_id] = unique(subject_id);


singlet_index = logical(([1; diff(merged_event_info.ripples_peaktimes)>0.1]));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Mean Bias Diff


% timebin = 0.01;
% time_windows = [0 0.1];
% % Generate bin edges
% bin_edges = time_windows(1):timebin:time_windows(2);
% % Generate bin centers
% bin_centers = bin_edges(1:end-1) + timebin/2;
% % 
% % 
% % % z_bias1 = KDE_reactivation_content.HPC_z_ripples';
% % % T1_events = find(nanmean(z_bias1(1:10,:))>0.5);
% % % T2_events = find(nanmean(z_bias1(1:10,:))<-0.5);
% % 
% % z_bias = KDE_reactivation_content.HPC_zshuff_ripples';
% % z_bias_V1 = KDE_reactivation_content.V1_zshuff_ripples';
% z_bias = KDE_reactivation_content.HPC_z_ripples';
% z_bias_V1 = KDE_reactivation_content.V1_z_ripples';


% event_index = (singlet_index+power_index)>1;

%% Temporal log odds AUC for all ripple events

% Time window parameters
win_size  = 0.1;   % 100 ms selection window for V1
step_size = 0.02;  % 50 ms step
time_bins = -1:step_size:1;
nTime = numel(time_bins);
nBoot = 1000;

% Fixed HPC selection window (always 0–0.1 s)
bins_to_use = bin_centers >= 0 & bin_centers < 0.1;

% Colour scheme (Using the darkest pink from your original scheme for all events)
colour_line = [231,  41, 138] / 256; 

% Storage (Simplified to 1D over time)
AUC.mean = nan(nTime, 1);
AUC.ci = nan(nTime, 2);
AUC.shifted_mean = nan(nTime, 1);
AUC.shifted_ci = nan(nTime, 2);

for t = 1:nTime
    t0 = time_bins(t) - win_size/2;
    t1 = time_bins(t) + win_size/2;
    
    % Sliding V1 window
    bins_to_select = bin_centers >= t0 & bin_centers < t1;
    fprintf('Processing V1 window %.3f–%.3f s (HPC fixed 0–0.1 s) [All Ripples]\n', t0, t1);
    
    % Select ALL ripple events
    event_index = true(size(ripple_info.ripple_power)); 
    total_events = sum(event_index);
    
    if total_events < 10, continue; end
    
    % Mean log-odds
    mean_bias = mean(z_bias(bins_to_use, event_index), 'omitnan');        % HPC fixed
    mean_bias_V1 = mean(z_bias_V1(bins_to_select, event_index), 'omitnan'); % V1 sliding
    
    % Quantile thresholds based on |HPC bias|
    thresholds = prctile(abs(mean_bias), 0:10:100);
    thresholds = thresholds(1:end-1);
    nThresh = numel(thresholds);
    
    % Bootstrap arrays
    bias_diff_boot = NaN(nBoot, nThresh);
    bias_diff_shift_boot = NaN(nBoot, nThresh);
    
    parfor iBoot = 1:nBoot
        s = RandStream('philox4x32_10', 'Seed', iBoot);
        idx = randi(s, total_events, total_events, 1); % resample event IDs
        
        boot_HPC = mean_bias(idx);
        boot_V1 = mean_bias_V1(idx);
        boot_bias_shifted = mean_bias; % event id randomised
        
        diff_tmp = NaN(1, nThresh);
        diff_tmp_shift = NaN(1, nThresh);
        
        for i = 1:nThresh
            th = thresholds(i);
            
            % Real HPC–V1 pairing
            t1 = boot_HPC >= th;
            t2 = boot_HPC <= -th;
            t1_V1 = boot_V1(t1);
            t2_V1 = boot_V1(t2);
            if ~isempty(t1_V1) && ~isempty(t2_V1)
                diff_tmp(i) = mean(t1_V1, 'omitnan') - mean(t2_V1, 'omitnan');
            end
            
            % Shifted pairing (HPC–V1 decoupled)
            t1s = boot_bias_shifted >= th;
            t2s = boot_bias_shifted <= -th;
            t1_V1s = boot_V1(t1s);
            t2_V1s = boot_V1(t2s);
            if ~isempty(t1_V1s) && ~isempty(t2_V1s)
                diff_tmp_shift(i) = mean(t1_V1s, 'omitnan') - mean(t2_V1s, 'omitnan');
            end
        end
        bias_diff_boot(iBoot, :) = diff_tmp;
        bias_diff_shift_boot(iBoot, :) = diff_tmp_shift;
    end
    
    % Quantile-based AUC (mean across thresholds)
    auc_boot = (trapz(thresholds, bias_diff_boot') / (max(thresholds)-min(thresholds)))';
    auc_shift_boot = (trapz(thresholds, bias_diff_shift_boot') / (max(thresholds)-min(thresholds)))';
    
    % Store summaries
    AUC.mean(t) = mean(auc_boot, 'omitnan');
    AUC.ci(t, :) = prctile(auc_boot, [2.5 97.5]);
    AUC.shifted_mean(t) = mean(auc_shift_boot, 'omitnan');
    AUC.shifted_ci(t, :) = prctile(auc_shift_boot, [2.5 97.5]);
end

%%% Save Data
save_dir = fullfile(analysis_folder, 'V1-HPC sleep reactivation');
if ~exist(save_dir, 'dir'), mkdir(save_dir); end

save(fullfile(save_dir, 'KDE_temporal_bias_all_ripples.mat'), 'AUC');

load(fullfile(save_dir, 'KDE_temporal_bias_all_ripples.mat'), 'AUC');

%%% Plot Temporal AUC Trace (All Events vs Shuffled)
fig = figure('Name', 'Temporal V1 log-odds AUC All Ripples', 'Position', [640 100 1100/3 900/4]);
% fig = figure('Name','Temporal V1 log-odds AUC low vs high ripple powers','Position',[640 100 1100/3 900/4]);
tiledlayout(nBins, 1, 'TileSpacing','compact');
hold on;

m  = AUC.mean;
ci = AUC.ci;
m_shift  = AUC.shifted_mean;
ci_shift = AUC.shifted_ci;
tvec = time_bins;

% Remove any NaN times where there weren't enough events
valid_idx = ~isnan(m);
tvec = tvec(valid_idx);
m = m(valid_idx);
ci = ci(valid_idx, :);
m_shift = m_shift(valid_idx);
ci_shift = ci_shift(valid_idx, :);

% Plot Shifted / Shuffled Baseline (Black)
fill([tvec, fliplr(tvec)], ...
     [ci_shift(:,1)', fliplr(ci_shift(:,2)')], ...
     [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.15);
plot(tvec, m_shift, 'k', 'LineWidth', 1.2, 'DisplayName', 'Shuffled Baseline');

% Plot Real Data (Coloured Pink)
fill([tvec, fliplr(tvec)], ...
     [ci(:,1)', fliplr(ci(:,2)')], ...
     colour_line, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
plot(tvec, m, 'Color', colour_line, 'LineWidth', 2, 'DisplayName', 'Real Data');

% Formatting
yline(0, '--r', 'LineWidth', 1);
xlabel('Time (s relative to ripple onset)');
ylabel('V1 bias AUC');
title('Temporal V1 Log-Odds AUC (All Ripple Events)');
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
xlim([-0.5 0.5]);
ylim([-0.1 0.2]);
legend('show', 'Location', 'best');

% Save Figures
save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation','temporal KDE bias difference'),[])



%% Temporal log odds AUC with different ripple power
% Ripple power binning
power_thresholds = prctile(ripple_info.ripple_power, 0:99.9/4:99.9);
nBins = length(power_thresholds) - 1;

% Time window parameters
win_size  = 0.1;   % 100 ms selection window for V1
step_size = 0.02;  % 50 ms step
time_bins = -1:step_size:1;
nTime = numel(time_bins);
nBoot = 1000;

% Fixed HPC selection window (always 0–0.1 s)
bins_to_use = bin_centers >= 0 & bin_centers < 0.1;

% Colour scheme
colour_lines = [ ...
    241, 182, 218;   % lightest
    226, 132, 187;   % mid-light
    212,  78, 156;   % mid-dark
    231,  41, 138    % darkest
] / 256;

% Storage
AUC.mean = nan(nTime, nBins);
AUC.ci = nan(nTime, nBins, 2);
AUC.shifted_mean = nan(nTime, nBins);
AUC.shifted_ci = nan(nTime, nBins, 2);

for t = 1:nTime
    t0 = time_bins(t)-win_size/2;
    t1 = time_bins(t) + win_size/2;

    % Sliding V1 window
    bins_to_select = bin_centers >= t0 & bin_centers < t1;

    fprintf('Processing V1 window %.3f–%.3f s (HPC fixed 0–0.1 s)\n', t0, t1);

    for npower = 1:nBins
        % Ripple power bin index
        power_index = ripple_info.ripple_power > power_thresholds(npower) & ...
                      ripple_info.ripple_power <= power_thresholds(npower+1);
        event_index = power_index > 0;
        total_events = sum(event_index);
        if total_events < 10, continue; end

        % Mean log-odds
        mean_bias = mean(z_bias(bins_to_use, event_index), 'omitnan');        % HPC fixed
        mean_bias_V1 = mean(z_bias_V1(bins_to_select, event_index), 'omitnan'); % V1 sliding

        % Quantile thresholds based on |HPC bias|
        thresholds = prctile(abs(mean_bias), 0:10:100);
        thresholds = thresholds(1:end-1);
        nThresh = numel(thresholds);

        % Bootstrap arrays
        bias_diff_boot = NaN(nBoot, nThresh);
        bias_diff_shift_boot = NaN(nBoot, nThresh);

        parfor iBoot = 1:nBoot
            s = RandStream('philox4x32_10', 'Seed', iBoot);
            idx = randi(s, total_events, total_events, 1); % resample event IDs
            boot_HPC = mean_bias(idx);
            boot_V1 = mean_bias_V1(idx);
            boot_bias_shifted = mean_bias; % event id randomised


            diff_tmp = NaN(1, nThresh);
            diff_tmp_shift = NaN(1, nThresh);

            for i = 1:nThresh
                th = thresholds(i);

                % Real HPC–V1 pairing
                t1 = boot_HPC >= th;
                t2 = boot_HPC <= -th;
                t1_V1 = boot_V1(t1);
                t2_V1 = boot_V1(t2);
                if ~isempty(t1_V1) && ~isempty(t2_V1)
                    diff_tmp(i) = mean(t1_V1, 'omitnan') - mean(t2_V1, 'omitnan');
                end

                % Shifted pairing (HPC–V1 decoupled)
                t1s = boot_bias_shifted >= th;
                t2s = boot_bias_shifted <= -th;
                t1_V1s = boot_V1(t1s);
                t2_V1s = boot_V1(t2s);
                if ~isempty(t1_V1s) && ~isempty(t2_V1s)
                    diff_tmp_shift(i) = mean(t1_V1s, 'omitnan') - mean(t2_V1s, 'omitnan');
                end
            end

            bias_diff_boot(iBoot, :) = diff_tmp;
            bias_diff_shift_boot(iBoot, :) = diff_tmp_shift;
        end

        % Quantile-based AUC (mean across thresholds)
        auc_boot = (trapz(thresholds, bias_diff_boot') / (max(thresholds)-min(thresholds)))';
        auc_shift_boot = (trapz(thresholds, bias_diff_shift_boot') / (max(thresholds)-min(thresholds)))';

        % Store summaries
        AUC.mean(t, npower) = mean(auc_boot, 'omitnan');
        AUC.ci(t, npower, :) = prctile(auc_boot, [2.5 97.5]);
        AUC.shifted_mean(t, npower) = mean(auc_shift_boot, 'omitnan');
        AUC.shifted_ci(t, npower, :) = prctile(auc_shift_boot, [2.5 97.5]);
    end
end

% Plot temporal AUC traces
fig = figure('Name','Temporal V1 log-odds AUC different ripple powers','Position',[640 100 1100/3 900]);
tiledlayout(nBins, 1, 'TileSpacing','compact');

for npower = 1:nBins
    nexttile; hold on;
    m  = AUC.mean(:,npower);
    ci = squeeze(AUC.ci(~isnan(m),npower,:));
    m_shift  = AUC.shifted_mean(~isnan(m),npower);
    ci_shift = squeeze(AUC.shifted_ci(~isnan(m),npower,:));
    tvec = time_bins(~isnan(m));
    m(isnan(m)) = [];


    % Real (coloured)
    fill([tvec fliplr(tvec)], ...
         [ci(:,1)' fliplr(ci(:,2)')], ...
         colour_lines(npower,:), 'EdgeColor','none','FaceAlpha',0.3);
    plot(tvec, m, 'Color', colour_lines(npower,:), 'LineWidth', 2);

    % Shifted (black)
    fill([tvec fliplr(tvec)], ...
         [ci_shift(:,1)' fliplr(ci_shift(:,2)')], ...
         [0 0 0], 'EdgeColor','none','FaceAlpha',0.15);
    plot(tvec, m_shift, 'k', 'LineWidth', 1.2);

    yline(0, '--r');
    xlabel('Time (s relative to ripple onset)');
    ylabel('V1 bias AUC');
    title(sprintf('Ripple power bin %d (%.2f–%.2f)', ...
          npower, power_thresholds(npower), power_thresholds(npower+1)));
    set(gca,'TickDir','out','Box','off','FontSize',12);
    xlim([-0.5 0.5]);
    ylim([-0.1 0.2])
end

save(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_temporal_bias_ripple_power.mat'),'AUC')
save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation','temporal KDE bias difference'),[])


load(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE_temporal_bias_ripple_power.mat'),'AUC')


% Plot temporal AUC traces
fig = figure('Name','Temporal V1 log-odds AUC low vs high ripple powers','Position',[640 100 1100/3 900/4]);
tiledlayout(nBins, 1, 'TileSpacing','compact');

for npower = [1 4]
    hold on;
    m  = AUC.mean(:,npower);
    ci = squeeze(AUC.ci(~isnan(m),npower,:));
    m_shift  = AUC.shifted_mean(~isnan(m),npower);
    ci_shift = squeeze(AUC.shifted_ci(~isnan(m),npower,:));
    tvec = time_bins(~isnan(m));
    m(isnan(m)) = [];


    % Real (coloured)
    fill([tvec fliplr(tvec)], ...
        [ci(:,1)' fliplr(ci(:,2)')], ...
        colour_lines(npower,:), 'EdgeColor','none','FaceAlpha',0.3);
    plot(tvec, m, 'Color', colour_lines(npower,:), 'LineWidth', 2);

    yline(0, '--r');
    xlabel('Time (s relative to ripple onset)');
    ylabel('V1 bias AUC');
    title(sprintf('Ripple power bin %d (%.2f–%.2f)', ...
        npower, power_thresholds(npower), power_thresholds(npower+1)));
    set(gca,'TickDir','out','Box','off','FontSize',12);
    xlim([-0.5 0.5]);
    ylim([-0.1 0.2])
end


% Shifted (black)
fill([tvec fliplr(tvec)], ...
    [ci_shift(:,1)' fliplr(ci_shift(:,2)')], ...
    [0 0 0], 'EdgeColor','none','FaceAlpha',0.15);
plot(tvec, m_shift, 'k', 'LineWidth', 1.2);
save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation','temporal KDE bias difference'),[])
