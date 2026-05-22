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

% load(fullfile(analysis_folder,'ripples_all_best_V1_SO_POST.mat'))

load(fullfile(analysis_folder,'spindles_all_POST.mat'))
load(fullfile(analysis_folder,'behavioural_state_merged_all_POST.mat'))
load(fullfile(analysis_folder,'V1-HPC sleep interaction','UP_DOWN_ripples_event_info.mat'),'event_info');
load(fullfile(analysis_folder,'V1-HPC sleep interaction','merged_UP_DOWN_ripples_event_info.mat'),'merged_event_info');

% load(fullfile(analysis_folder,'V1-HPC sleep interaction','UP_DOWN_ripple_PSTH_MUA.mat'),'UP_DOWN_ripple_PSTH_MUA');
% PSTH_MUA = UP_DOWN_ripple_PSTH_MUA;
% load(fullfile(analysis_folder,'V1-HPC sleep interaction','UP_DOWN_ripple_PSTH_MUA_baseline.mat'));
% PSTH_MUA_baseline = UP_DOWN_ripple_PSTH_MUA;

% load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_whole.mat'));
% probability_psth_whole = probability;

% load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_whole.mat'));

% load(fullfile(analysis_folder,'bayesian_reactivation_V1_all_POST.mat'))
% load(fullfile(analysis_folder,'bayesian_reactivation_all_POST.mat'))

sessions_to_process = 1:max(slow_waves_all(1).UP_session_count);
%% Extract key information for DOWN UP, ripples and spindles
 
% load((fullfile(analysis_folder,'cortex_SO_ref_shank_all.mat')),'cortex_SO_ref_shank_all');
%%%%%%% Find reference channel/shank
cortex_ref_shank = [];
HPC_ref_shank = [];

for nsession = 1:max(ripples_all(1).session_count)
    for probe_no = 1:2
        cortex_ref_shank(nsession,probe_no) = find(slow_waves_all(probe_no).shank_id{nsession} == slow_waves_all(probe_no).shank{nsession}(slow_waves_all(probe_no).channel{nsession} == slow_waves_all(probe_no).best_channel(nsession))...
            &slow_waves_all(probe_no).probe_hemisphere{nsession} == probe_no);

        % cortex_ref_shank(nsession,probe_no) = find(slow_waves_all(probe_no).shank_id{nsession}==cortex_SO_ref_shank_all(nsession,probe_no) ...
        %     & slow_waves_all(probe_no).probe_hemisphere{nsession}==probe_no);
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
ripple_info.spindle_presence_hemi = ripple_info.spindle_presence_hemi(event_ids_first);

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


%%% spindle power
spindle_amplitude=[];
for probe_no = 1:2
    spindle_amplitude{probe_no}=[];
    for nsession = 1:length(sessions_to_process)
        spindle_amplitude{probe_no} = [spindle_amplitude{probe_no} ripples_all(probe_no).spindle_amplitude_ripple_onset{nsession}(cortex_ref_shank(nsession,:),:)];
    end
end

spindle_amplitude = [spindle_amplitude{1}(:,ripples_all(1).SWS_index==1) spindle_amplitude{2}(:,ripples_all(2).SWS_index==1)];

ripple_info.spindle_amplitude_onset = spindle_amplitude';
ripple_info.spindle_amplitude_onset = ripple_info.spindle_amplitude_onset(event_ids_first,:);


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


%%%%%%
session_count = [ripples_all(1).session_count(ripples_all(1).SWS_index==1); ripples_all(2).session_count(ripples_all(2).SWS_index==1)];
subject_id = str2double(cellstr(ripples_all(1).subject(session_count,end-1:end)));
[~, ~, subject_id] = unique(subject_id);


% singlet_index = logical(([1; diff(merged_event_info.ripples_peaktimes)>0.1]));

% merged_event_info.ripples_hemisphere_id
% singlet_index = logical(ones(length(merged_event_info.ripples_peaktimes),1));



% singlet_index = logical(ones(length(merged_event_info.ripples_peaktimes),1));
cd(fullfile(analysis_folder,'V1-HPC sleep reactivation'))


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

z_bias = KDE_reactivation_ripples_PSTH.HPC_z_logodds_ripples';
z_bias_V1 = KDE_reactivation_ripples_PSTH.V1_z_logodds_ripples';

z_bias1 = z_bias(isfinite(z_bias));
z_bias(z_bias>=inf) = prctile(z_bias1,99.5);
z_bias(z_bias<=-inf) = prctile(z_bias1,0.5);

z_bias1 = z_bias(isfinite(z_bias_V1));
z_bias_V1(z_bias_V1>=inf) = prctile(z_bias1,99.5);
z_bias_V1(z_bias_V1<=-inf) = prctile(z_bias1,0.5);

%%%% Only grab unique ripple events
z_bias = z_bias(:,event_ids_first);
z_bias_V1 = z_bias_V1(:,event_ids_first);
% event_ids_first
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

bins_to_use = bin_centers>0 & bin_centers<0.2;
session_count = [ripples_all(1).session_count(ripples_all(1).SWS_index==1); ripples_all(2).session_count(ripples_all(2).SWS_index==1)];
subject_id = str2double(cellstr(ripples_all(1).subject(session_count,end-1:end)));
[~, ~, subject_id] = unique(subject_id);


singlet_index = logical(([1; diff(merged_event_info.ripples_peaktimes)>0.1]));



%%%%%%%%%%%%%%%%%%%%%%%% Spindle power binning

%%% spindle power
spindle_amplitude=[];
spindle_percentile = [];
for probe_no = 1:2
    spindle_amplitude{probe_no}=[];
    spindle_percentile{probe_no}=[];
    for nsession = 1:length(sessions_to_process)
        % spindle_amplitude{probe_no} = [spindle_amplitude{probe_no} ripples_all(probe_no).spindle_amplitude_ripple_onset{nsession}(cortex_ref_shank(nsession,:),:)];
        % spindle_amplitude{probe_no} = [spindle_amplitude{probe_no} ripples_all(probe_no).spindle_amplitude_ripple_peaktime{nsession}(cortex_ref_shank(nsession,:),:)];
        temp = ripples_all(probe_no).spindle_amplitude_ripple_peaktime{nsession}(cortex_ref_shank(nsession,:),:);
        % temp = ripples_all(probe_no).spindle_amplitude_ripple_peaktime{nsession}(cortex_ref_shank(nsession,:),:);
        percentiles = (tiedrank(temp')' - 0.5) / length(temp) * 100;
        spindle_percentile{probe_no} = [spindle_percentile{probe_no} percentiles];
    end
end
spindle_percentile = [spindle_percentile{1}(:,ripples_all(1).SWS_index==1) spindle_percentile{2}(:,ripples_all(2).SWS_index==1)];
ripple_info.spindle_percentile = spindle_percentile';
ripple_info.spindle_percentile = ripple_info.spindle_percentile(event_ids_first,:);
ripple_info.spindle_amplitude = ripple_info.spindle_percentile;


% Spindle power binning across both probes
% all_spindle_power = mean(ripple_info.spindle_amplitude, 1);  % avg of probe 1 and 2
spindle_thresholds = prctile(ripple_info.spindle_amplitude, 0:99.9/4:99.9);
nBins = length(spindle_thresholds) - 1;

bins_to_use = bin_centers>0 & bin_centers<0.1;
% bins_to_select = bin_centers>0 & bin_centers<0.2;

bins_to_select = bin_centers>-0.1 & bin_centers<0.1;
% bins_to_select = bin_centers>-0.1 & bin_centers<0.1;

bins_to_use_shifted = bin_centers > -1 & bin_centers < -0.9;
nBoot = 1000;
colour_lines = [ ...
    241, 182, 218;   % original end (lightest)
    226, 132, 187;   % interpolated 2/3
    212,  78, 156;   % interpolated 1/3
    231, 41, 138    % original start (darkest)
] / 256;

spindle_power_KDE_bias_difference = struct;

fig = figure;
fig.Position = [640 100 1100 650*2];
fig.Name = 'KDE bias difference in HPC with different spindle powers percentile (-0.2 to 0.2s)';
% fig.Name = 'KDE bias difference in HPC with different spindle powers percentile';
tiledlayout(nBins, 3, 'TileSpacing', 'compact');

for npower = 1:nBins
    tic
    % Select events in current spindle bin
    % power_index = ripple_info.spindle_amplitude(:,1) > spindle_thresholds(npower,1) & ...
    %               ripple_info.spindle_amplitude(:,1) <= spindle_thresholds(npower+1,1) |...
    %               ripple_info.spindle_amplitude(:,2) > spindle_thresholds(npower,2) & ...
    %               ripple_info.spindle_amplitude(:,2) <= spindle_thresholds(npower+1,2);

    event_index =1:length(z_bias);
    % event_index = find(power_index);

    mean_bias = mean(z_bias_V1(bins_to_select, event_index), 'omitnan');
    % mean_bias_shifted = mean(z_bias(bins_to_use_shifted, event_index), 'omitnan');
    mean_bias_V1 = mean(z_bias(bins_to_use, event_index), 'omitnan');
    selected_events = length(mean_bias);

    thresholds = prctile(abs(mean_bias), 0:10:100);
    thresholds = thresholds(1:end-1);
    nThresh = length(thresholds);

    bias_diff_boot = NaN(nBoot, nThresh);
    prop_events_boot = NaN(nBoot, nThresh);
    bias_diff_shifted_boot = NaN(nBoot, nThresh);
    prop_events_shifted_boot = NaN(nBoot, nThresh);


    parfor iBoot = 1:nBoot
        s = RandStream('philox4x32_10', 'Seed', iBoot);
        idx = randi(s, selected_events, selected_events, 1);

        true_idx = find(event_index);

        bb = mean_bias(idx);
        bb_shift = mean_bias;
        boot_V1 = mean_bias_V1(idx);

        diff_tmp = NaN(1, nThresh);
        prop_tmp = NaN(1, nThresh);
        diff_tmp_shifted = NaN(1, nThresh);
        prop_tmp_shifted = NaN(1, nThresh);

        for i = 1:nThresh
            th = thresholds(i);

            t1 = bb >= th;
            t2 = bb <= -th;



            t1 = t1' + (ripple_info.spindle_amplitude( true_idx(idx),2) > spindle_thresholds(npower,2) &...
                ripple_info.spindle_amplitude( true_idx(idx),2) <= spindle_thresholds(npower+1,2))>1';
            t2 = t2' + (ripple_info.spindle_amplitude( true_idx(idx),1) > spindle_thresholds(npower,1) &...
                ripple_info.spindle_amplitude( true_idx(idx),1) <= spindle_thresholds(npower+1,1))>1';

            t1_V1 = boot_V1(t1);
            t2_V1 = boot_V1(t2);
            if any(t1) && any(t2)
                diff_tmp(i) = mean(t1_V1) - mean(t2_V1);
            end

            total_events = mean([sum((ripple_info.spindle_amplitude( true_idx(idx),2) > spindle_thresholds(npower,2) &...
                ripple_info.spindle_amplitude( true_idx(idx),2) <= spindle_thresholds(npower+1,2))) ...
                sum((ripple_info.spindle_amplitude( true_idx(idx),1) > spindle_thresholds(npower,1) &...
                ripple_info.spindle_amplitude( true_idx(idx),1) <= spindle_thresholds(npower+1,1)))]);

            prop_tmp(i) = (sum(t1) + sum(t2)) / total_events;



            t1s = bb_shift >= th;
            t2s = bb_shift <= -th;
            t1s = t1s' + (ripple_info.spindle_amplitude( true_idx,2) > spindle_thresholds(npower,2) &...
                ripple_info.spindle_amplitude( true_idx,2) <= spindle_thresholds(npower+1,2))>1';
            t2s = t2s' + (ripple_info.spindle_amplitude( true_idx,1) > spindle_thresholds(npower,1) &...
                ripple_info.spindle_amplitude( true_idx,1) <= spindle_thresholds(npower+1,1))>1';
            t1_V1 = boot_V1(t1s);
            t2_V1 = boot_V1(t2s);

            if any(t1s) && any(t2s)
                diff_tmp_shifted(i) = mean(t1_V1) - mean(t2_V1);
            end
            prop_tmp_shifted(i) = (sum(t1s) + sum(t2s)) / total_events;
        end

        bias_diff_boot(iBoot, :) = diff_tmp;
        prop_events_boot(iBoot, :) = prop_tmp;
        bias_diff_shifted_boot(iBoot, :) = diff_tmp_shifted;
        prop_events_shifted_boot(iBoot, :) = prop_tmp_shifted;
    end

    % Stats
    bias_mean = mean(bias_diff_boot, 1, 'omitnan');
    bias_CI_lo = prctile(bias_diff_boot, 2.5, 1);
    bias_CI_hi = prctile(bias_diff_boot, 97.5, 1);
    prop_mean = mean(prop_events_boot, 1, 'omitnan');
    prop_CI_lo = prctile(prop_events_boot, 2.5, 1);
    prop_CI_hi = prctile(prop_events_boot, 97.5, 1);
    bias_shifted_mean = mean(bias_diff_shifted_boot, 1, 'omitnan');
    bias_shifted_CI_lo = prctile(bias_diff_shifted_boot, 2.5, 1);
    bias_shifted_CI_hi = prctile(bias_diff_shifted_boot, 97.5, 1);
    prop_shifted_mean = mean(prop_events_shifted_boot, 1, 'omitnan');
    prop_shifted_CI_lo = prctile(prop_events_shifted_boot, 2.5, 1);
    prop_shifted_CI_hi = prctile(prop_events_shifted_boot, 97.5, 1);

    % Store
    spindle_power_KDE_bias_difference(npower).power_range = ...
        [spindle_thresholds(npower), spindle_thresholds(npower+1)];
    spindle_power_KDE_bias_difference(npower).bias_diff_mean = bias_mean;
    spindle_power_KDE_bias_difference(npower).bias_diff_CI = [bias_CI_lo; bias_CI_hi];
    spindle_power_KDE_bias_difference(npower).prop_mean = prop_mean;
    spindle_power_KDE_bias_difference(npower).prop_CI = [prop_CI_lo; prop_CI_hi];
    spindle_power_KDE_bias_difference(npower).thresholds = thresholds;
    spindle_power_KDE_bias_difference(npower).bias_diff_shifted_mean = bias_shifted_mean;
    spindle_power_KDE_bias_difference(npower).bias_diff_shifted_CI = [bias_shifted_CI_lo; bias_shifted_CI_hi];
    spindle_power_KDE_bias_difference(npower).prop_shifted_mean = prop_shifted_mean;
    spindle_power_KDE_bias_difference(npower).prop_shifted_CI = [prop_shifted_CI_lo; prop_shifted_CI_hi];

    % store AUC
    auc_boot = (trapz(thresholds, bias_diff_boot') / (max(thresholds)-min(thresholds)))';
    auc_shift_boot = (trapz(thresholds, bias_diff_shifted_boot') / (max(thresholds)-min(thresholds)))';

    spindle_power_KDE_bias_difference(npower).AUC_mean = mean(auc_boot, 'omitnan');
    spindle_power_KDE_bias_difference(npower).AUC_CI = prctile(auc_boot, [2.5 97.5]);
    spindle_power_KDE_bias_difference(npower).AUC_mean_shuffled = mean(auc_shift_boot, 'omitnan');
    spindle_power_KDE_bias_difference(npower).AUC_CI_shuffled = prctile(auc_shift_boot, [2.5 97.5]);

    % ----------- PLOTS (A, B, C) ----------------

    % A: Bias vs. Threshold
    nexttile((npower-1)*3 + 1);
    hold on;
    fill([thresholds, fliplr(thresholds)], [bias_CI_lo, fliplr(bias_CI_hi)], ...
        colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    plot(thresholds, bias_mean, 'Color', colour_lines(npower,:), 'LineWidth', 2);
    fill([thresholds, fliplr(thresholds)], ...
         [bias_shifted_CI_lo, fliplr(bias_shifted_CI_hi)], ...
         [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    plot(thresholds, bias_shifted_mean, 'k-', 'LineWidth', 1.5);
    ylim([-0.15 0.35]); xlim([0 1]); yline(0, '--r');
    xlabel('V1 bias threshold'); ylabel('HPC bias diff (T1 - T2)');
    title(sprintf('Spindle bin %d: %.2f–%.2f', npower, ...
           spindle_thresholds(npower), spindle_thresholds(npower+1)));
    set(gca,"TickDir","out",'box','off','Color','none','FontSize',12)

    % B: Proportion vs Bias Diff
    nexttile((npower-1)*3 + 2);
    hold on;
    valid = isfinite(bias_mean) & isfinite(prop_mean);
    fill([bias_CI_lo(valid), fliplr(bias_CI_hi(valid))], ...
         [prop_mean(valid), fliplr(prop_mean(valid))], ...
         colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    plot(bias_mean(valid), prop_mean(valid), '-', 'Color', colour_lines(npower,:), 'LineWidth', 2);
    fill([bias_shifted_CI_lo(valid), fliplr(bias_shifted_CI_hi(valid))], ...
         [prop_shifted_mean(valid), fliplr(prop_shifted_mean(valid))], ...
         [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    plot(bias_shifted_mean(valid), prop_shifted_mean(valid), 'k-', 'LineWidth', 1.5);
    xlim([-0.1 0.35]); xline(0, '--r');
    xlabel('HPC bias diff (T1 - T2)'); ylabel('Proportion of events detected');
    title('Event Proportion vs. Bias Difference');
    set(gca,"TickDir","out",'box','off','Color','none','FontSize',12)

    % C: Proportion CI vs Bias Diff
    nexttile((npower-1)*3 + 3);
    hold on;
    fill([bias_mean(valid), fliplr(bias_mean(valid))], ...
         [prop_CI_lo(valid), fliplr(prop_CI_hi(valid))], ...
         colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    plot(bias_mean(valid), prop_mean(valid), '-', 'Color', colour_lines(npower,:), 'LineWidth', 2);
    fill([bias_shifted_mean(valid), fliplr(bias_shifted_mean(valid))], ...
         [prop_shifted_CI_lo(valid), fliplr(prop_shifted_CI_hi(valid))], ...
         [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    plot(bias_shifted_mean(valid), prop_shifted_mean(valid), 'k-', 'LineWidth', 1.5);
    xlim([-0.1 0.35]); xline(0, '--r');
    xlabel('HPC bias diff (T1 - T2)'); ylabel('Proportion of events detected');
    title('Proportion vs. HPC Bias Difference');
    set(gca,"TickDir","out",'box','off','Color','none','FontSize',12)
    toc
end

clear Fill

fig = figure;
fig.Position = [640 100 2*1100/3 650/2];
fig.Name = 'KDE bias difference in HPC low vs high spindle power percentile';
% fig.Name = 'KDE bias difference in HPC low vs high spindle power percentile (100ms example)';

% Select bins (1 = low, 4 = high)
nexttile;
for npower = [1 4]
    bias_mean = spindle_power_KDE_bias_difference(npower).bias_diff_mean;
    bias_CI_lo = spindle_power_KDE_bias_difference(npower).bias_diff_CI(1,:);
    bias_CI_hi = spindle_power_KDE_bias_difference(npower).bias_diff_CI(2,:);
    thresholds = spindle_power_KDE_bias_difference(npower).thresholds;

    hold on;
    x2 = [thresholds, fliplr(thresholds)];
    y2 = [bias_CI_lo, fliplr(bias_CI_hi)];
    Fill(npower) = fill(x2, y2, colour_lines(npower,:), ...
        'EdgeColor', 'none', 'FaceAlpha', 0.3);
    plot(thresholds, bias_mean, 'Color', colour_lines(npower,:), 'LineWidth', 2);

    xlabel('HPC Bias threshold');
    ylabel('V1 bias diff (T1 - T2)');
    set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12);
    ylim([-0.1 0.35]);
    % ylim([0 1])
end

bias_mean = spindle_power_KDE_bias_difference(npower).bias_diff_shifted_mean;
bias_CI_lo = spindle_power_KDE_bias_difference(npower).bias_diff_shifted_CI(1,:);
bias_CI_hi = spindle_power_KDE_bias_difference(npower).bias_diff_shifted_CI(2,:)

hold on;
x2 = [thresholds, fliplr(thresholds)];
y2 = [bias_CI_lo, fliplr(bias_CI_hi)];
Fill(end +1) = fill(x2, y2, 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
plot(thresholds, bias_mean, 'Color','k', 'LineWidth', 2);

yline(0, '--r');
legend(Fill([1 4 5]), {'Low spindle power', 'High spindle power','Shuffled'}, 'box', 'off');

nexttile;
for npower = [1 4]
    bias_mean = spindle_power_KDE_bias_difference(npower).bias_diff_mean;
    bias_CI_lo = spindle_power_KDE_bias_difference(npower).bias_diff_CI(1,:);
    bias_CI_hi = spindle_power_KDE_bias_difference(npower).bias_diff_CI(2,:);
    prop_mean = spindle_power_KDE_bias_difference(npower).prop_mean;

    hold on;
    y2 = [prop_mean, fliplr(prop_mean)];
    x2 = [bias_CI_lo, fliplr(bias_CI_hi)];
    Fill(npower) = fill(x2, y2, colour_lines(npower,:), ...
        'EdgeColor', 'none', 'FaceAlpha', 0.3);
    plot(bias_mean, prop_mean, 'Color', colour_lines(npower,:), 'LineWidth', 2);

    xlabel('V1 bias diff (T1 - T2)');
    ylabel('Proportion of events detected');
    set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12);
    xlim([-0.1 0.35]);
    ylim([0 1])
end


bias_mean = spindle_power_KDE_bias_difference(npower).bias_diff_shifted_mean;
bias_CI_lo = spindle_power_KDE_bias_difference(npower).bias_diff_shifted_CI(1,:);
bias_CI_hi = spindle_power_KDE_bias_difference(npower).bias_diff_shifted_CI(2,:)
prop_mean = spindle_power_KDE_bias_difference(npower).prop_shifted_mean;

hold on;
y2 = [prop_mean, fliplr(prop_mean)];
x2 = [bias_CI_lo, fliplr(bias_CI_hi)];
Fill(end + 1) = fill(x2, y2, 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
plot(bias_mean, prop_mean, 'Color','k', 'LineWidth', 2);


xline(0, '--r');
legend(Fill([1 4 5]), {'Low spindle power', 'High spindle power','Shuffled'}, 'box', 'off');




%%%%% AUC mean + CI bar plot
% Plot layout
fig = figure;
fig.Position = [640 100 281 325]
fig.Name = 'KDE bias V1 AUC spindle power percentile low vs high (-0.2 to 0.2s)';
% fig.Name = 'KDE bias V1 AUC spindle power percentile low vs high';
data = spindle_power_KDE_bias_difference;
n_bins = length(data);
bar_width = 0.3;      % Width of the bars
group_offset = 0.15;    % Distance from the center integer (half the gap between bars)
hold on;
clear BAR
for i = 1:4
    if i > length(data)
        x_shuf = i;
        y_shuf = data(1).AUC_mean_shuffled;
        bar(x_shuf, y_shuf, bar_width, ...
            'FaceColor', 'k', ...
            'FaceAlpha', 0.15, ...
            'EdgeColor', 'none');
    else


        % --- 1. Plot Shuffled Data (Left Bar) ---
        x_shuf = i - group_offset;
        y_shuf = data(i).AUC_mean_shuffled;

        % Calculate Error Deltas (Errorbar requires length relative to mean, not absolute values)
        % CI is [lower, upper]
        neg_err_shuf = y_shuf - data(i).AUC_CI_shuffled(1);
        pos_err_shuf = data(i).AUC_CI_shuffled(2) - y_shuf;

        % Plot Bar
        if i == 1
            BAR(1) = bar(x_shuf, y_shuf, bar_width, ...
                'FaceColor', 'k', ...
                'FaceAlpha', 0.15, ...
                'EdgeColor', 'none');
        else
            bar(x_shuf, y_shuf, bar_width, ...
                'FaceColor', 'k', ...
                'FaceAlpha', 0.15, ...
                'EdgeColor', 'none');
        end

        % Plot Error Bar
        E = errorbar(x_shuf, y_shuf, neg_err_shuf, pos_err_shuf, ...
            'Color', 'k', 'LineWidth', 1.5, 'CapSize', 8, 'LineStyle', 'none');


        % --- 2. Plot Real Data (Right Bar) ---
        x_real = i + group_offset;
        y_real = data(i).AUC_mean;

        neg_err_real = y_real - data(i).AUC_CI(1);
        pos_err_real = data(i).AUC_CI(2) - y_real;

        % Get specific color for this power bin
        % Assuming colour_lines is size [n_bins x 3]
        this_color = colour_lines(i, :);

        % Plot Bar
        BAR(i+1) = bar(x_real, y_real, bar_width, ...
            'FaceColor', this_color, ...
            'FaceAlpha', 0.3, ...
            'EdgeColor', 'none');

        errorbar(x_real, y_real, neg_err_real, pos_err_real, ...
            'Color', 'k', 'LineWidth', 1.5, 'CapSize', 8, 'LineStyle', 'none');
    end
end

hold off;
% Set X-ticks to be centered on the groups
set(gca, 'XTick', 1:nBins);
xlim([0.5, 4 + 0.5]);

% Labels
ylabel('V1 bias AUC');
xlabel('Power Bins');
legend([BAR(1:end)],{'Shuffled','Low','High'},'box','off')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)



% Save results
% save('spindle_power_percentile_KDE_bias_difference_based_on_V1_bias.mat', 'spindle_power_KDE_bias_difference');
save('spindle_power_percentile_KDE_bias_difference_based_on_V1_bias.mat', 'spindle_power_KDE_bias_difference_200ms');
save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE bias difference based on V1 bias'),[])


% load('spindle_power_percentile_KDE_bias_difference_based_on_V1_bias.mat', 'spindle_power_KDE_bias_difference');



%%%%%%%%%%%%%%%%%%%%%%%% Spindle power binning
% Spindle power binning across both probes
% all_spindle_power = mean(ripple_info.spindle_amplitude, 1);  % avg of probe 1 and 2
spindle_thresholds = prctile(ripple_info.spindle_amplitude, 0:99.9/4:99.9);
nBins = length(spindle_thresholds) - 1;

bins_to_use = bin_centers>0 & bin_centers<0.1;
bins_to_select = bin_centers>-0.2 & bin_centers<0;
bins_to_use_shifted = bin_centers > -1 & bin_centers < -0.9;
nBoot = 1000;
colour_lines = [ ...
    241, 182, 218;   % original end (lightest)
    226, 132, 187;   % interpolated 2/3
    212,  78, 156;   % interpolated 1/3
    231, 41, 138    % original start (darkest)
] / 256;

spindle_power_KDE_bias_difference = struct;

fig = figure;
fig.Position = [640 100 1100 650*2];
fig.Name = 'KDE bias difference in HPC with different spindle powers percentile (PRE ripple)';
tiledlayout(nBins, 3, 'TileSpacing', 'compact');

for npower = 1:nBins
    tic
    % Select events in current spindle bin
    % power_index = ripple_info.spindle_amplitude(:,1) > spindle_thresholds(npower,1) & ...
    %               ripple_info.spindle_amplitude(:,1) <= spindle_thresholds(npower+1,1) |...
    %               ripple_info.spindle_amplitude(:,2) > spindle_thresholds(npower,2) & ...
    %               ripple_info.spindle_amplitude(:,2) <= spindle_thresholds(npower+1,2);

    event_index =1:length(z_bias);
    % event_index = find(power_index);

    mean_bias = mean(z_bias_V1(bins_to_select, event_index), 'omitnan');
    % mean_bias_shifted = mean(z_bias(bins_to_use_shifted, event_index), 'omitnan');
    mean_bias_V1 = mean(z_bias(bins_to_use, event_index), 'omitnan');
    selected_events = length(mean_bias);

    thresholds = prctile(abs(mean_bias), 0:10:100);
    thresholds = thresholds(1:end-1);
    nThresh = length(thresholds);

    bias_diff_boot = NaN(nBoot, nThresh);
    prop_events_boot = NaN(nBoot, nThresh);
    bias_diff_shifted_boot = NaN(nBoot, nThresh);
    prop_events_shifted_boot = NaN(nBoot, nThresh);


    parfor iBoot = 1:nBoot
        s = RandStream('philox4x32_10', 'Seed', iBoot);
        idx = randi(s, selected_events, selected_events, 1);

        true_idx = find(event_index);

        bb = mean_bias(idx);
        bb_shift = mean_bias;
        boot_V1 = mean_bias_V1(idx);

        diff_tmp = NaN(1, nThresh);
        prop_tmp = NaN(1, nThresh);
        diff_tmp_shifted = NaN(1, nThresh);
        prop_tmp_shifted = NaN(1, nThresh);

        for i = 1:nThresh
            th = thresholds(i);

            t1 = bb >= th;
            t2 = bb <= -th;



            t1 = t1' + (ripple_info.spindle_amplitude( true_idx(idx),2) > spindle_thresholds(npower,2) &...
                ripple_info.spindle_amplitude( true_idx(idx),2) <= spindle_thresholds(npower+1,2))>1';
            t2 = t2' + (ripple_info.spindle_amplitude( true_idx(idx),1) > spindle_thresholds(npower,1) &...
                ripple_info.spindle_amplitude( true_idx(idx),1) <= spindle_thresholds(npower+1,1))>1';

            t1_V1 = boot_V1(t1);
            t2_V1 = boot_V1(t2);
            if any(t1) && any(t2)
                diff_tmp(i) = mean(t1_V1) - mean(t2_V1);
            end

            total_events = mean([sum((ripple_info.spindle_amplitude( true_idx(idx),2) > spindle_thresholds(npower,2) &...
                ripple_info.spindle_amplitude( true_idx(idx),2) <= spindle_thresholds(npower+1,2))) ...
                sum((ripple_info.spindle_amplitude( true_idx(idx),1) > spindle_thresholds(npower,1) &...
                ripple_info.spindle_amplitude( true_idx(idx),1) <= spindle_thresholds(npower+1,1)))]);

            prop_tmp(i) = (sum(t1) + sum(t2)) / total_events;



            t1s = bb_shift >= th;
            t2s = bb_shift <= -th;
            t1s = t1s' + (ripple_info.spindle_amplitude( true_idx,2) > spindle_thresholds(npower,2) &...
                ripple_info.spindle_amplitude( true_idx,2) <= spindle_thresholds(npower+1,2))>1';
            t2s = t2s' + (ripple_info.spindle_amplitude( true_idx,1) > spindle_thresholds(npower,1) &...
                ripple_info.spindle_amplitude( true_idx,1) <= spindle_thresholds(npower+1,1))>1';
            t1_V1 = boot_V1(t1s);
            t2_V1 = boot_V1(t2s);

            if any(t1s) && any(t2s)
                diff_tmp_shifted(i) = mean(t1_V1) - mean(t2_V1);
            end
            prop_tmp_shifted(i) = (sum(t1s) + sum(t2s)) / total_events;
        end

        bias_diff_boot(iBoot, :) = diff_tmp;
        prop_events_boot(iBoot, :) = prop_tmp;
        bias_diff_shifted_boot(iBoot, :) = diff_tmp_shifted;
        prop_events_shifted_boot(iBoot, :) = prop_tmp_shifted;
    end

    % Stats
    bias_mean = mean(bias_diff_boot, 1, 'omitnan');
    bias_CI_lo = prctile(bias_diff_boot, 2.5, 1);
    bias_CI_hi = prctile(bias_diff_boot, 97.5, 1);
    prop_mean = mean(prop_events_boot, 1, 'omitnan');
    prop_CI_lo = prctile(prop_events_boot, 2.5, 1);
    prop_CI_hi = prctile(prop_events_boot, 97.5, 1);
    bias_shifted_mean = mean(bias_diff_shifted_boot, 1, 'omitnan');
    bias_shifted_CI_lo = prctile(bias_diff_shifted_boot, 2.5, 1);
    bias_shifted_CI_hi = prctile(bias_diff_shifted_boot, 97.5, 1);
    prop_shifted_mean = mean(prop_events_shifted_boot, 1, 'omitnan');
    prop_shifted_CI_lo = prctile(prop_events_shifted_boot, 2.5, 1);
    prop_shifted_CI_hi = prctile(prop_events_shifted_boot, 97.5, 1);

    % Store
    spindle_power_KDE_bias_difference(npower).power_range = ...
        [spindle_thresholds(npower), spindle_thresholds(npower+1)];
    spindle_power_KDE_bias_difference(npower).bias_diff_mean = bias_mean;
    spindle_power_KDE_bias_difference(npower).bias_diff_CI = [bias_CI_lo; bias_CI_hi];
    spindle_power_KDE_bias_difference(npower).prop_mean = prop_mean;
    spindle_power_KDE_bias_difference(npower).prop_CI = [prop_CI_lo; prop_CI_hi];
    spindle_power_KDE_bias_difference(npower).thresholds = thresholds;
    spindle_power_KDE_bias_difference(npower).bias_diff_shifted_mean = bias_shifted_mean;
    spindle_power_KDE_bias_difference(npower).bias_diff_shifted_CI = [bias_shifted_CI_lo; bias_shifted_CI_hi];
    spindle_power_KDE_bias_difference(npower).prop_shifted_mean = prop_shifted_mean;
    spindle_power_KDE_bias_difference(npower).prop_shifted_CI = [prop_shifted_CI_lo; prop_shifted_CI_hi];

    % store AUC
    auc_boot = (trapz(thresholds, bias_diff_boot') / (max(thresholds)-min(thresholds)))';
    auc_shift_boot = (trapz(thresholds, bias_diff_shifted_boot') / (max(thresholds)-min(thresholds)))';

    spindle_power_KDE_bias_difference(npower).AUC_mean = mean(auc_boot, 'omitnan');
    spindle_power_KDE_bias_difference(npower).AUC_CI = prctile(auc_boot, [2.5 97.5]);
    spindle_power_KDE_bias_difference(npower).AUC_mean_shuffled = mean(auc_shift_boot, 'omitnan');
    spindle_power_KDE_bias_difference(npower).AUC_CI_shuffled = prctile(auc_shift_boot, [2.5 97.5]);

    % ----------- PLOTS (A, B, C) ----------------

    % A: Bias vs. Threshold
    nexttile((npower-1)*3 + 1);
    hold on;
    fill([thresholds, fliplr(thresholds)], [bias_CI_lo, fliplr(bias_CI_hi)], ...
        colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    plot(thresholds, bias_mean, 'Color', colour_lines(npower,:), 'LineWidth', 2);
    fill([thresholds, fliplr(thresholds)], ...
         [bias_shifted_CI_lo, fliplr(bias_shifted_CI_hi)], ...
         [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    plot(thresholds, bias_shifted_mean, 'k-', 'LineWidth', 1.5);
    ylim([-0.15 0.35]); xlim([0 1]); yline(0, '--r');
    xlabel('V1 bias threshold'); ylabel('HPC bias diff (T1 - T2)');
    title(sprintf('Spindle bin %d: %.2f–%.2f', npower, ...
           spindle_thresholds(npower), spindle_thresholds(npower+1)));
    set(gca,"TickDir","out",'box','off','Color','none','FontSize',12)

    % B: Proportion vs Bias Diff
    nexttile((npower-1)*3 + 2);
    hold on;
    valid = isfinite(bias_mean) & isfinite(prop_mean);
    fill([bias_CI_lo(valid), fliplr(bias_CI_hi(valid))], ...
         [prop_mean(valid), fliplr(prop_mean(valid))], ...
         colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    plot(bias_mean(valid), prop_mean(valid), '-', 'Color', colour_lines(npower,:), 'LineWidth', 2);
    fill([bias_shifted_CI_lo(valid), fliplr(bias_shifted_CI_hi(valid))], ...
         [prop_shifted_mean(valid), fliplr(prop_shifted_mean(valid))], ...
         [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    plot(bias_shifted_mean(valid), prop_shifted_mean(valid), 'k-', 'LineWidth', 1.5);
    xlim([-0.1 0.35]); xline(0, '--r');
    xlabel('HPC bias diff (T1 - T2)'); ylabel('Proportion of events detected');
    title('Event Proportion vs. Bias Difference');
    set(gca,"TickDir","out",'box','off','Color','none','FontSize',12)

    % C: Proportion CI vs Bias Diff
    nexttile((npower-1)*3 + 3);
    hold on;
    fill([bias_mean(valid), fliplr(bias_mean(valid))], ...
         [prop_CI_lo(valid), fliplr(prop_CI_hi(valid))], ...
         colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    plot(bias_mean(valid), prop_mean(valid), '-', 'Color', colour_lines(npower,:), 'LineWidth', 2);
    fill([bias_shifted_mean(valid), fliplr(bias_shifted_mean(valid))], ...
         [prop_shifted_CI_lo(valid), fliplr(prop_shifted_CI_hi(valid))], ...
         [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    plot(bias_shifted_mean(valid), prop_shifted_mean(valid), 'k-', 'LineWidth', 1.5);
    xlim([-0.1 0.35]); xline(0, '--r');
    xlabel('HPC bias diff (T1 - T2)'); ylabel('Proportion of events detected');
    title('Proportion vs. HPC Bias Difference');
    set(gca,"TickDir","out",'box','off','Color','none','FontSize',12)
    toc
end

clear Fill

fig = figure;
fig.Position = [640 100 2*1100/3 650/2];
fig.Name = 'KDE bias difference in HPC low vs high spindle power percentile (PRE ripple)';

% Select bins (1 = low, 4 = high)
nexttile;
for npower = [1 4]
    bias_mean = spindle_power_KDE_bias_difference(npower).bias_diff_mean;
    bias_CI_lo = spindle_power_KDE_bias_difference(npower).bias_diff_CI(1,:);
    bias_CI_hi = spindle_power_KDE_bias_difference(npower).bias_diff_CI(2,:);
    thresholds = spindle_power_KDE_bias_difference(npower).thresholds;

    hold on;
    x2 = [thresholds, fliplr(thresholds)];
    y2 = [bias_CI_lo, fliplr(bias_CI_hi)];
    Fill(npower) = fill(x2, y2, colour_lines(npower,:), ...
        'EdgeColor', 'none', 'FaceAlpha', 0.3);
    plot(thresholds, bias_mean, 'Color', colour_lines(npower,:), 'LineWidth', 2);

    xlabel('HPC Bias threshold');
    ylabel('V1 bias diff (T1 - T2)');
    set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12);
    ylim([-0.1 0.35]);
    % ylim([0 1])
end

bias_mean = spindle_power_KDE_bias_difference(npower).bias_diff_shifted_mean;
bias_CI_lo = spindle_power_KDE_bias_difference(npower).bias_diff_shifted_CI(1,:);
bias_CI_hi = spindle_power_KDE_bias_difference(npower).bias_diff_shifted_CI(2,:)

hold on;
x2 = [thresholds, fliplr(thresholds)];
y2 = [bias_CI_lo, fliplr(bias_CI_hi)];
Fill(end +1) = fill(x2, y2, 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
plot(thresholds, bias_mean, 'Color','k', 'LineWidth', 2);

yline(0, '--r');
legend(Fill([1 4 5]), {'Low spindle power', 'High spindle power','Shuffled'}, 'box', 'off');

nexttile;
for npower = [1 4]
    bias_mean = spindle_power_KDE_bias_difference(npower).bias_diff_mean;
    bias_CI_lo = spindle_power_KDE_bias_difference(npower).bias_diff_CI(1,:);
    bias_CI_hi = spindle_power_KDE_bias_difference(npower).bias_diff_CI(2,:);
    prop_mean = spindle_power_KDE_bias_difference(npower).prop_mean;

    hold on;
    y2 = [prop_mean, fliplr(prop_mean)];
    x2 = [bias_CI_lo, fliplr(bias_CI_hi)];
    Fill(npower) = fill(x2, y2, colour_lines(npower,:), ...
        'EdgeColor', 'none', 'FaceAlpha', 0.3);
    plot(bias_mean, prop_mean, 'Color', colour_lines(npower,:), 'LineWidth', 2);

    xlabel('V1 bias diff (T1 - T2)');
    ylabel('Proportion of events detected');
    set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12);
    xlim([-0.1 0.35]);
    ylim([0 1])
end


bias_mean = spindle_power_KDE_bias_difference(npower).bias_diff_shifted_mean;
bias_CI_lo = spindle_power_KDE_bias_difference(npower).bias_diff_shifted_CI(1,:);
bias_CI_hi = spindle_power_KDE_bias_difference(npower).bias_diff_shifted_CI(2,:)
prop_mean = spindle_power_KDE_bias_difference(npower).prop_shifted_mean;

hold on;
y2 = [prop_mean, fliplr(prop_mean)];
x2 = [bias_CI_lo, fliplr(bias_CI_hi)];
Fill(end + 1) = fill(x2, y2, 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
plot(bias_mean, prop_mean, 'Color','k', 'LineWidth', 2);


xline(0, '--r');
legend(Fill([1 4 5]), {'Low spindle power', 'High spindle power','Shuffled'}, 'box', 'off');




%%%%% AUC mean + CI bar plot
% Plot layout
fig = figure;
fig.Position = [640 100 281 325]
fig.Name = 'KDE bias V1 AUC spindle power percentile low vs high (PRE)';
data = spindle_power_KDE_bias_difference;
n_bins = length(data);
bar_width = 0.3;      % Width of the bars
group_offset = 0.15;    % Distance from the center integer (half the gap between bars)
hold on;
clear BAR
for i = 1:4
    if i > length(data)
        x_shuf = i;
        y_shuf = data(1).AUC_mean_shuffled;
        bar(x_shuf, y_shuf, bar_width, ...
            'FaceColor', 'k', ...
            'FaceAlpha', 0.15, ...
            'EdgeColor', 'none');
    else


        % --- 1. Plot Shuffled Data (Left Bar) ---
        x_shuf = i - group_offset;
        y_shuf = data(i).AUC_mean_shuffled;

        % Calculate Error Deltas (Errorbar requires length relative to mean, not absolute values)
        % CI is [lower, upper]
        neg_err_shuf = y_shuf - data(i).AUC_CI_shuffled(1);
        pos_err_shuf = data(i).AUC_CI_shuffled(2) - y_shuf;

        % Plot Bar
        if i == 1
            BAR(1) = bar(x_shuf, y_shuf, bar_width, ...
                'FaceColor', 'k', ...
                'FaceAlpha', 0.15, ...
                'EdgeColor', 'none');
        else
            bar(x_shuf, y_shuf, bar_width, ...
                'FaceColor', 'k', ...
                'FaceAlpha', 0.15, ...
                'EdgeColor', 'none');
        end

        % Plot Error Bar
        E = errorbar(x_shuf, y_shuf, neg_err_shuf, pos_err_shuf, ...
            'Color', 'k', 'LineWidth', 1.5, 'CapSize', 8, 'LineStyle', 'none');


        % --- 2. Plot Real Data (Right Bar) ---
        x_real = i + group_offset;
        y_real = data(i).AUC_mean;

        neg_err_real = y_real - data(i).AUC_CI(1);
        pos_err_real = data(i).AUC_CI(2) - y_real;

        % Get specific color for this power bin
        % Assuming colour_lines is size [n_bins x 3]
        this_color = colour_lines(i, :);

        % Plot Bar
        BAR(i+1) = bar(x_real, y_real, bar_width, ...
            'FaceColor', this_color, ...
            'FaceAlpha', 0.3, ...
            'EdgeColor', 'none');

        errorbar(x_real, y_real, neg_err_real, pos_err_real, ...
            'Color', 'k', 'LineWidth', 1.5, 'CapSize', 8, 'LineStyle', 'none');
    end
end

hold off;
% Set X-ticks to be centered on the groups
set(gca, 'XTick', 1:nBins);
xlim([0.5, 4 + 0.5]);

% Labels
ylabel('V1 bias AUC');
xlabel('Power Bins');
legend([BAR(1:end)],{'Shuffled','Low','High'},'box','off')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)



% Save results
save('spindle_power_percentile_KDE_bias_difference_based_on_PRE_V1_bias.mat', 'spindle_power_KDE_bias_difference');
save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE bias difference based on V1 bias'),[])




%%%%%%%%%%%%%%%%%%%%%%%% Spindle power binning unilateral vs bilateral
% Spindle power binning across both probes
% all_spindle_power = mean(ripple_info.spindle_amplitude, 1);  % avg of probe 1 and 2
spindle_thresholds = prctile(ripple_info.spindle_amplitude,[50]);
nBins = 2;

bins_to_use = bin_centers>0 & bin_centers<0.1;
bins_to_select = bin_centers>-0.2 & bin_centers<0;
bins_to_use_shifted = bin_centers > -1 & bin_centers < -0.9;
nBoot = 1000;
colour_lines = [ ...
    241, 182, 218;   % original end (lightest)
    % 226, 132, 187;   % interpolated 2/3
    % 212,  78, 156;   % interpolated 1/3
    231, 41, 138    % original start (darkest)
] / 256;

spindle_power_KDE_bias_difference = struct;

fig = figure;
fig.Position = [640 100 1100 650];
fig.Name = 'KDE bias difference in HPC with different spindle power synchrony percentile (PRE ripple)';
tiledlayout(nBins, 3, 'TileSpacing', 'compact');

for npower = 1:nBins
    tic
    % Select events in current spindle bin
    % power_index = ripple_info.spindle_amplitude(:,1) > spindle_thresholds(npower,1) & ...
    %               ripple_info.spindle_amplitude(:,1) <= spindle_thresholds(npower+1,1) |...
    %               ripple_info.spindle_amplitude(:,2) > spindle_thresholds(npower,2) & ...
    %               ripple_info.spindle_amplitude(:,2) <= spindle_thresholds(npower+1,2);

    event_index =1:length(z_bias);
    % event_index = find(power_index);

    mean_bias = mean(z_bias_V1(bins_to_select, event_index), 'omitnan');
    % mean_bias_shifted = mean(z_bias(bins_to_use_shifted, event_index), 'omitnan');
    mean_bias_V1 = mean(z_bias(bins_to_use, event_index), 'omitnan');
    selected_events = length(mean_bias);

    thresholds = prctile(abs(mean_bias), 0:10:100);
    thresholds = thresholds(1:end-1);
    nThresh = length(thresholds);

    bias_diff_boot = NaN(nBoot, nThresh);
    prop_events_boot = NaN(nBoot, nThresh);
    bias_diff_shifted_boot = NaN(nBoot, nThresh);
    prop_events_shifted_boot = NaN(nBoot, nThresh);


    parfor iBoot = 1:nBoot
        s = RandStream('philox4x32_10', 'Seed', iBoot);
        idx = randi(s, selected_events, selected_events, 1);

        true_idx = find(event_index);

        bb = mean_bias(idx);
        bb_shift = mean_bias;
        boot_V1 = mean_bias_V1(idx);

        diff_tmp = NaN(1, nThresh);
        prop_tmp = NaN(1, nThresh);
        diff_tmp_shifted = NaN(1, nThresh);
        prop_tmp_shifted = NaN(1, nThresh);
        total_events = [];
        for i = 1:nThresh
            th = thresholds(i);

            t1 = bb >= th;
            t2 = bb <= -th;


            if npower == 1 % only context matching high spindle
                t1 = t1' + (ripple_info.spindle_amplitude(true_idx(idx),2) > spindle_thresholds(1,2) & ...
                    ripple_info.spindle_amplitude(true_idx(idx),1) < spindle_thresholds(1,1) ) > 1';
                t2 = t2' + (ripple_info.spindle_amplitude(true_idx(idx),1) > spindle_thresholds(1,1) & ...
                    ripple_info.spindle_amplitude(true_idx(idx),2) < spindle_thresholds(1,2) ) > 1';

                total_events = mean([sum((ripple_info.spindle_amplitude( :,2) > spindle_thresholds(1,2) &...
                    ripple_info.spindle_amplitude( :,1) < spindle_thresholds(1,1))) ...
                    sum((ripple_info.spindle_amplitude( :,1) > spindle_thresholds(1,1) &...
                    ripple_info.spindle_amplitude(  :,2) < spindle_thresholds(1,2)))]);

            elseif npower == 2  % Both high spindle
                t1 = t1' + ( ripple_info.spindle_amplitude(true_idx(idx),2) > spindle_thresholds(1,2) & ...
                    ripple_info.spindle_amplitude(true_idx(idx),1) > spindle_thresholds(1,1) ) > 1';
                t2 = t2' + (ripple_info.spindle_amplitude(true_idx(idx),1) > spindle_thresholds(1,1) & ...
                    ripple_info.spindle_amplitude(true_idx(idx),2) > spindle_thresholds(1,2) ) > 1';

                total_events = mean([sum((ripple_info.spindle_amplitude( :,2) > spindle_thresholds(1,2) &...
                    ripple_info.spindle_amplitude( :,1) > spindle_thresholds(1,1)))]);
            end

            t1_V1 = boot_V1(t1);
            t2_V1 = boot_V1(t2);
            if any(t1) && any(t2)
                diff_tmp(i) = mean(t1_V1) - mean(t2_V1);
            end



            prop_tmp(i) = (sum(t1) + sum(t2)) / total_events;



            t1s = bb_shift >= th;
            t2s = bb_shift <= -th;
            if npower == 1
                t1s = t1s' + (ripple_info.spindle_amplitude(true_idx,2) > spindle_thresholds(1,2) & ...
                    ripple_info.spindle_amplitude(true_idx,1) < spindle_thresholds(1,1)) > 1';
                t2s = t2s' + (ripple_info.spindle_amplitude(true_idx,1) > spindle_thresholds(1,1) & ...
                    ripple_info.spindle_amplitude(true_idx,2) < spindle_thresholds(1,2)) > 1';
            else
                t1s = t1s' + (ripple_info.spindle_amplitude(true_idx,2) > spindle_thresholds(1,2) & ...
                    ripple_info.spindle_amplitude(true_idx,1) > spindle_thresholds(1,1)) > 1';
                t2s = t2s' + (ripple_info.spindle_amplitude(true_idx,1) > spindle_thresholds(1,1) & ...
                    ripple_info.spindle_amplitude(true_idx,2) > spindle_thresholds(1,2)) > 1';
            end
            t1_V1 = boot_V1(t1s);
            t2_V1 = boot_V1(t2s);

            if any(t1s) && any(t2s)
                diff_tmp_shifted(i) = mean(t1_V1) - mean(t2_V1);
            end
            prop_tmp_shifted(i) = (sum(t1s) + sum(t2s)) / total_events;
        end

        bias_diff_boot(iBoot, :) = diff_tmp;
        prop_events_boot(iBoot, :) = prop_tmp;
        bias_diff_shifted_boot(iBoot, :) = diff_tmp_shifted;
        prop_events_shifted_boot(iBoot, :) = prop_tmp_shifted;
    end

    % Stats
    bias_mean = mean(bias_diff_boot, 1, 'omitnan');
    bias_CI_lo = prctile(bias_diff_boot, 2.5, 1);
    bias_CI_hi = prctile(bias_diff_boot, 97.5, 1);
    prop_mean = mean(prop_events_boot, 1, 'omitnan');
    prop_CI_lo = prctile(prop_events_boot, 2.5, 1);
    prop_CI_hi = prctile(prop_events_boot, 97.5, 1);
    bias_shifted_mean = mean(bias_diff_shifted_boot, 1, 'omitnan');
    bias_shifted_CI_lo = prctile(bias_diff_shifted_boot, 2.5, 1);
    bias_shifted_CI_hi = prctile(bias_diff_shifted_boot, 97.5, 1);
    prop_shifted_mean = mean(prop_events_shifted_boot, 1, 'omitnan');
    prop_shifted_CI_lo = prctile(prop_events_shifted_boot, 2.5, 1);
    prop_shifted_CI_hi = prctile(prop_events_shifted_boot, 97.5, 1);

    % Store
    % spindle_power_KDE_bias_difference(npower).power_range = ...
    %     [spindle_thresholds(npower), spindle_thresholds(npower+1)];
    spindle_power_KDE_bias_difference(npower).bias_diff_mean = bias_mean;
    spindle_power_KDE_bias_difference(npower).bias_diff_CI = [bias_CI_lo; bias_CI_hi];
    spindle_power_KDE_bias_difference(npower).prop_mean = prop_mean;
    spindle_power_KDE_bias_difference(npower).prop_CI = [prop_CI_lo; prop_CI_hi];
    spindle_power_KDE_bias_difference(npower).thresholds = thresholds;
    spindle_power_KDE_bias_difference(npower).bias_diff_shifted_mean = bias_shifted_mean;
    spindle_power_KDE_bias_difference(npower).bias_diff_shifted_CI = [bias_shifted_CI_lo; bias_shifted_CI_hi];
    spindle_power_KDE_bias_difference(npower).prop_shifted_mean = prop_shifted_mean;
    spindle_power_KDE_bias_difference(npower).prop_shifted_CI = [prop_shifted_CI_lo; prop_shifted_CI_hi];

    % store AUC
    auc_boot = (trapz(thresholds, bias_diff_boot') / (max(thresholds)-min(thresholds)))';
    auc_shift_boot = (trapz(thresholds, bias_diff_shifted_boot') / (max(thresholds)-min(thresholds)))';

    spindle_power_KDE_bias_difference(npower).AUC_mean = mean(auc_boot, 'omitnan');
    spindle_power_KDE_bias_difference(npower).AUC_CI = prctile(auc_boot, [2.5 97.5]);
    spindle_power_KDE_bias_difference(npower).AUC_mean_shuffled = mean(auc_shift_boot, 'omitnan');
    spindle_power_KDE_bias_difference(npower).AUC_CI_shuffled = prctile(auc_shift_boot, [2.5 97.5]);

    % ----------- PLOTS (A, B, C) ----------------

    % A: Bias vs. Threshold
    nexttile((npower-1)*3 + 1);
    hold on;
    fill([thresholds, fliplr(thresholds)], [bias_CI_lo, fliplr(bias_CI_hi)], ...
        colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    plot(thresholds, bias_mean, 'Color', colour_lines(npower,:), 'LineWidth', 2);
    fill([thresholds, fliplr(thresholds)], ...
         [bias_shifted_CI_lo, fliplr(bias_shifted_CI_hi)], ...
         [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    plot(thresholds, bias_shifted_mean, 'k-', 'LineWidth', 1.5);
    ylim([-0.15 0.35]); xlim([0 1]); yline(0, '--r');
    xlabel('V1 bias threshold'); ylabel('HPC bias diff (T1 - T2)');
    % title(sprintf('Spindle bin %d: %.2f–%.2f', npower, ...
    %        spindle_thresholds(npower), spindle_thresholds(npower+1)));
    set(gca,"TickDir","out",'box','off','Color','none','FontSize',12)

    % B: Proportion vs Bias Diff
    nexttile((npower-1)*3 + 2);
    hold on;
    valid = isfinite(bias_mean) & isfinite(prop_mean);
    fill([bias_CI_lo(valid), fliplr(bias_CI_hi(valid))], ...
         [prop_mean(valid), fliplr(prop_mean(valid))], ...
         colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    plot(bias_mean(valid), prop_mean(valid), '-', 'Color', colour_lines(npower,:), 'LineWidth', 2);
    fill([bias_shifted_CI_lo(valid), fliplr(bias_shifted_CI_hi(valid))], ...
         [prop_shifted_mean(valid), fliplr(prop_shifted_mean(valid))], ...
         [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    plot(bias_shifted_mean(valid), prop_shifted_mean(valid), 'k-', 'LineWidth', 1.5);
    xlim([-0.1 0.35]); xline(0, '--r');
    xlabel('HPC bias diff (T1 - T2)'); ylabel('Proportion of events detected');
    title('Event Proportion vs. Bias Difference');
    set(gca,"TickDir","out",'box','off','Color','none','FontSize',12)

    % C: Proportion CI vs Bias Diff
    nexttile((npower-1)*3 + 3);
    hold on;
    fill([bias_mean(valid), fliplr(bias_mean(valid))], ...
         [prop_CI_lo(valid), fliplr(prop_CI_hi(valid))], ...
         colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    plot(bias_mean(valid), prop_mean(valid), '-', 'Color', colour_lines(npower,:), 'LineWidth', 2);
    fill([bias_shifted_mean(valid), fliplr(bias_shifted_mean(valid))], ...
         [prop_shifted_CI_lo(valid), fliplr(prop_shifted_CI_hi(valid))], ...
         [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    plot(bias_shifted_mean(valid), prop_shifted_mean(valid), 'k-', 'LineWidth', 1.5);
    xlim([-0.1 0.35]); xline(0, '--r');
    xlabel('HPC bias diff (T1 - T2)'); ylabel('Proportion of events detected');
    title('Proportion vs. HPC Bias Difference');
    set(gca,"TickDir","out",'box','off','Color','none','FontSize',12)
    toc
end

clear Fill

fig = figure;
fig.Position = [640 100 2*1100/3 650/2];
fig.Name = 'KDE bias difference in HPC unilateral vs bilateral spindle power percentile (PRE ripple)';

% Select bins (1 = low, 4 = high)
nexttile;
for npower = [1 2]
    bias_mean = spindle_power_KDE_bias_difference(npower).bias_diff_mean;
    bias_CI_lo = spindle_power_KDE_bias_difference(npower).bias_diff_CI(1,:);
    bias_CI_hi = spindle_power_KDE_bias_difference(npower).bias_diff_CI(2,:);
    thresholds = spindle_power_KDE_bias_difference(npower).thresholds;

    hold on;
    x2 = [thresholds, fliplr(thresholds)];
    y2 = [bias_CI_lo, fliplr(bias_CI_hi)];
    Fill(npower) = fill(x2, y2, colour_lines(npower,:), ...
        'EdgeColor', 'none', 'FaceAlpha', 0.3);
    plot(thresholds, bias_mean, 'Color', colour_lines(npower,:), 'LineWidth', 2);

    xlabel('HPC Bias threshold');
    ylabel('V1 bias diff (T1 - T2)');
    set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12);
    ylim([-0.1 0.35]);
    % ylim([0 1])
end

bias_mean = spindle_power_KDE_bias_difference(npower).bias_diff_shifted_mean;
bias_CI_lo = spindle_power_KDE_bias_difference(npower).bias_diff_shifted_CI(1,:);
bias_CI_hi = spindle_power_KDE_bias_difference(npower).bias_diff_shifted_CI(2,:)

hold on;
x2 = [thresholds, fliplr(thresholds)];
y2 = [bias_CI_lo, fliplr(bias_CI_hi)];
Fill(end +1) = fill(x2, y2, 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
plot(thresholds, bias_mean, 'Color','k', 'LineWidth', 2);

yline(0, '--r');
legend(Fill([1 2 3]), {'unilateral', 'bilateral','Shuffled'}, 'box', 'off');

nexttile;
for npower = [1 2]
    bias_mean = spindle_power_KDE_bias_difference(npower).bias_diff_mean;
    bias_CI_lo = spindle_power_KDE_bias_difference(npower).bias_diff_CI(1,:);
    bias_CI_hi = spindle_power_KDE_bias_difference(npower).bias_diff_CI(2,:);
    prop_mean = spindle_power_KDE_bias_difference(npower).prop_mean;

    hold on;
    y2 = [prop_mean, fliplr(prop_mean)];
    x2 = [bias_CI_lo, fliplr(bias_CI_hi)];
    Fill(npower) = fill(x2, y2, colour_lines(npower,:), ...
        'EdgeColor', 'none', 'FaceAlpha', 0.3);
    plot(bias_mean, prop_mean, 'Color', colour_lines(npower,:), 'LineWidth', 2);

    xlabel('V1 bias diff (T1 - T2)');
    ylabel('Proportion of events detected');
    set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12);
    xlim([-0.1 0.35]);
    ylim([0 1])
end


bias_mean = spindle_power_KDE_bias_difference(npower).bias_diff_shifted_mean;
bias_CI_lo = spindle_power_KDE_bias_difference(npower).bias_diff_shifted_CI(1,:);
bias_CI_hi = spindle_power_KDE_bias_difference(npower).bias_diff_shifted_CI(2,:)
prop_mean = spindle_power_KDE_bias_difference(npower).prop_shifted_mean;

hold on;
y2 = [prop_mean, fliplr(prop_mean)];
x2 = [bias_CI_lo, fliplr(bias_CI_hi)];
Fill(end + 1) = fill(x2, y2, 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
plot(bias_mean, prop_mean, 'Color','k', 'LineWidth', 2);


xline(0, '--r');
legend(Fill([1 2 3]), {'unilateral', 'bilateral','Shuffled'}, 'box', 'off');



%%%%% AUC mean + CI bar plot
% Plot layout
fig = figure;
fig.Position = [640 100 281 325]
fig.Name = 'KDE bias V1 AUC unilateral vs bilateral spindle power percentile low vs high (PRE)';
data = spindle_power_KDE_bias_difference;
n_bins = length(data);
bar_width = 0.3;      % Width of the bars
group_offset = 0.15;    % Distance from the center integer (half the gap between bars)
hold on;
clear BAR
for i = 1:4
    if i > length(data)
        x_shuf = i;
        y_shuf = data(1).AUC_mean_shuffled;
        bar(x_shuf, y_shuf, bar_width, ...
            'FaceColor', 'k', ...
            'FaceAlpha', 0.15, ...
            'EdgeColor', 'none');
    else


        % --- 1. Plot Shuffled Data (Left Bar) ---
        x_shuf = i - group_offset;
        y_shuf = data(i).AUC_mean_shuffled;

        % Calculate Error Deltas (Errorbar requires length relative to mean, not absolute values)
        % CI is [lower, upper]
        neg_err_shuf = y_shuf - data(i).AUC_CI_shuffled(1);
        pos_err_shuf = data(i).AUC_CI_shuffled(2) - y_shuf;

        % Plot Bar
        if i == 1
            BAR(1) = bar(x_shuf, y_shuf, bar_width, ...
                'FaceColor', 'k', ...
                'FaceAlpha', 0.15, ...
                'EdgeColor', 'none');
        else
            bar(x_shuf, y_shuf, bar_width, ...
                'FaceColor', 'k', ...
                'FaceAlpha', 0.15, ...
                'EdgeColor', 'none');
        end

        % Plot Error Bar
        E = errorbar(x_shuf, y_shuf, neg_err_shuf, pos_err_shuf, ...
            'Color', 'k', 'LineWidth', 1.5, 'CapSize', 8, 'LineStyle', 'none');


        % --- 2. Plot Real Data (Right Bar) ---
        x_real = i + group_offset;
        y_real = data(i).AUC_mean;

        neg_err_real = y_real - data(i).AUC_CI(1);
        pos_err_real = data(i).AUC_CI(2) - y_real;

        % Get specific color for this power bin
        % Assuming colour_lines is size [n_bins x 3]
        this_color = colour_lines(i, :);

        % Plot Bar
        BAR(i+1) = bar(x_real, y_real, bar_width, ...
            'FaceColor', this_color, ...
            'FaceAlpha', 0.3, ...
            'EdgeColor', 'none');

        errorbar(x_real, y_real, neg_err_real, pos_err_real, ...
            'Color', 'k', 'LineWidth', 1.5, 'CapSize', 8, 'LineStyle', 'none');
    end
end

hold off;
% Set X-ticks to be centered on the groups
set(gca, 'XTick', 1:nBins);
xlim([0.5, 4 + 0.5]);

% Labels
ylabel('V1 bias AUC');
xlabel('Power Bins');
legend([BAR(1:end)],{'Shuffled','unilateral','bilateral'},'box','off')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)



% Save results
save('spindle_power_synchrony_percentile_KDE_bias_difference_based_on_PRE_V1_bias.mat', 'spindle_power_KDE_bias_difference');
save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE bias difference based on V1 bias'),[])










%%%%%%%%%%%%%%%%%%%%%%%% delta power binning
% Spindle power binning across both probes
% all_spindle_power = mean(ripple_info.spindle_amplitude, 1);  % avg of probe 1 and 2
SO_thresholds = prctile(ripple_info.SO_amplitude, 0:99.9/4:99.9);
nBins = length(SO_thresholds) - 1;

bins_to_use = bin_centers>0 & bin_centers<0.1;
bins_to_select = bin_centers>0 & bin_centers<0.2;
bins_to_use_shifted = bin_centers > -1 & bin_centers < -0.9;
nBoot = 1000;

SO_power_KDE_bias_difference = struct;

colour_lines = [ ...
    241, 182, 218;   % original end (lightest)
    226, 132, 187;   % interpolated 2/3
    212,  78, 156;   % interpolated 1/3
    231, 41, 138    % original start (darkest)
] / 256;


fig = figure;
fig.Position = [640 100 1100 650*2];
fig.Name = 'KDE bias difference in HPC with different SO powers';
tiledlayout(nBins, 3, 'TileSpacing', 'compact');

for npower = 1:nBins
    tic
    % Select events in current spindle bin
    % power_index = ripple_info.spindle_amplitude(:,1) > spindle_thresholds(npower,1) & ...
    %               ripple_info.spindle_amplitude(:,1) <= spindle_thresholds(npower+1,1) |...
    %               ripple_info.spindle_amplitude(:,2) > spindle_thresholds(npower,2) & ...
    %               ripple_info.spindle_amplitude(:,2) <= spindle_thresholds(npower+1,2);

    event_index =1:length(z_bias);
    % event_index = find(power_index);

    mean_bias = mean(z_bias_V1(bins_to_select, event_index), 'omitnan');
    % mean_bias_shifted = mean(z_bias(bins_to_use_shifted, event_index), 'omitnan');
    mean_bias_V1 = mean(z_bias(bins_to_use, event_index), 'omitnan');
    selected_events = length(mean_bias);

    thresholds = prctile(abs(mean_bias), 0:10:100);
    thresholds = thresholds(1:end-1);
    nThresh = length(thresholds);

    bias_diff_boot = NaN(nBoot, nThresh);
    prop_events_boot = NaN(nBoot, nThresh);
    bias_diff_shifted_boot = NaN(nBoot, nThresh);
    prop_events_shifted_boot = NaN(nBoot, nThresh);


    parfor iBoot = 1:nBoot
        s = RandStream('philox4x32_10', 'Seed', iBoot);
        idx = randi(s, selected_events, selected_events, 1);

        true_idx = find(event_index);

        bb = mean_bias(idx);
        bb_shift = mean_bias;
        boot_V1 = mean_bias_V1(idx);

        diff_tmp = NaN(1, nThresh);
        prop_tmp = NaN(1, nThresh);
        diff_tmp_shifted = NaN(1, nThresh);
        prop_tmp_shifted = NaN(1, nThresh);

        for i = 1:nThresh
            th = thresholds(i);

            t1 = bb >= th;
            t2 = bb <= -th;



            t1 = t1' + (ripple_info.SO_amplitude( true_idx(idx),2) > SO_thresholds(npower,2) &...
                ripple_info.SO_amplitude( true_idx(idx),2) <= SO_thresholds(npower+1,2))>1';
            t2 = t2' + (ripple_info.SO_amplitude( true_idx(idx),1) > SO_thresholds(npower,1) &...
                ripple_info.SO_amplitude( true_idx(idx),1) <= SO_thresholds(npower+1,1))>1';

            t1_V1 = boot_V1(t1);
            t2_V1 = boot_V1(t2);
            if any(t1) && any(t2)
                diff_tmp(i) = mean(t1_V1) - mean(t2_V1);
            end

            total_events = mean([sum((ripple_info.SO_amplitude( true_idx(idx),2) > SO_thresholds(npower,2) &...
                ripple_info.SO_amplitude( true_idx(idx),2) <= SO_thresholds(npower+1,2))) ...
                sum((ripple_info.SO_amplitude( true_idx(idx),1) > SO_thresholds(npower,1) &...
                ripple_info.SO_amplitude( true_idx(idx),1) <= SO_thresholds(npower+1,1)))]);

            prop_tmp(i) = (sum(t1) + sum(t2)) / total_events;



            t1s = bb_shift >= th;
            t2s = bb_shift <= -th;
            t1s = t1s' + (ripple_info.SO_amplitude( true_idx,2) > SO_thresholds(npower,2) &...
                ripple_info.SO_amplitude( true_idx,2) <= SO_thresholds(npower+1,2))>1';
            t2s = t2s' + (ripple_info.SO_amplitude( true_idx,1) > SO_thresholds(npower,1) &...
                ripple_info.SO_amplitude( true_idx,1) <= SO_thresholds(npower+1,1))>1';
            t1_V1 = boot_V1(t1s);
            t2_V1 = boot_V1(t2s);

            if any(t1s) && any(t2s)
                diff_tmp_shifted(i) = mean(t1_V1) - mean(t2_V1);
            end
            prop_tmp_shifted(i) = (sum(t1s) + sum(t2s)) / total_events;
        end

        bias_diff_boot(iBoot, :) = diff_tmp;
        prop_events_boot(iBoot, :) = prop_tmp;
        bias_diff_shifted_boot(iBoot, :) = diff_tmp_shifted;
        prop_events_shifted_boot(iBoot, :) = prop_tmp_shifted;
    end

    % Stats
    bias_mean = mean(bias_diff_boot, 1, 'omitnan');
    bias_CI_lo = prctile(bias_diff_boot, 2.5, 1);
    bias_CI_hi = prctile(bias_diff_boot, 97.5, 1);
    prop_mean = mean(prop_events_boot, 1, 'omitnan');
    prop_CI_lo = prctile(prop_events_boot, 2.5, 1);
    prop_CI_hi = prctile(prop_events_boot, 97.5, 1);
    bias_shifted_mean = mean(bias_diff_shifted_boot, 1, 'omitnan');
    bias_shifted_CI_lo = prctile(bias_diff_shifted_boot, 2.5, 1);
    bias_shifted_CI_hi = prctile(bias_diff_shifted_boot, 97.5, 1);
    prop_shifted_mean = mean(prop_events_shifted_boot, 1, 'omitnan');
    prop_shifted_CI_lo = prctile(prop_events_shifted_boot, 2.5, 1);
    prop_shifted_CI_hi = prctile(prop_events_shifted_boot, 97.5, 1);

    % Store
    SO_power_KDE_bias_difference(npower).power_range = ...
        [SO_thresholds(npower), SO_thresholds(npower+1)];
    SO_power_KDE_bias_difference(npower).bias_diff_mean = bias_mean;
    SO_power_KDE_bias_difference(npower).bias_diff_CI = [bias_CI_lo; bias_CI_hi];
    SO_power_KDE_bias_difference(npower).prop_mean = prop_mean;
    SO_power_KDE_bias_difference(npower).prop_CI = [prop_CI_lo; prop_CI_hi];
    SO_power_KDE_bias_difference(npower).thresholds = thresholds;
    SO_power_KDE_bias_difference(npower).bias_diff_shifted_mean = bias_shifted_mean;
    SO_power_KDE_bias_difference(npower).bias_diff_shifted_CI = [bias_shifted_CI_lo; bias_shifted_CI_hi];
    SO_power_KDE_bias_difference(npower).prop_shifted_mean = prop_shifted_mean;
    SO_power_KDE_bias_difference(npower).prop_shifted_CI = [prop_shifted_CI_lo; prop_shifted_CI_hi];


    % store AUC
    auc_boot = (trapz(thresholds, bias_diff_boot') / (max(thresholds)-min(thresholds)))';
    auc_shift_boot = (trapz(thresholds, bias_diff_shifted_boot') / (max(thresholds)-min(thresholds)))';

    SO_power_KDE_bias_difference(npower).AUC_mean = mean(auc_boot, 'omitnan');
    SO_power_KDE_bias_difference(npower).AUC_CI = prctile(auc_boot, [2.5 97.5]);
    SO_power_KDE_bias_difference(npower).AUC_mean_shuffled = mean(auc_shift_boot, 'omitnan');
    SO_power_KDE_bias_difference(npower).AUC_CI_shuffled = prctile(auc_shift_boot, [2.5 97.5]);

    % ----------- PLOTS (A, B, C) ----------------

    % A: Bias vs. Threshold
    nexttile((npower-1)*3 + 1);
    hold on;
    fill([thresholds, fliplr(thresholds)], [bias_CI_lo, fliplr(bias_CI_hi)], ...
        colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    plot(thresholds, bias_mean, 'Color', colour_lines(npower,:), 'LineWidth', 2);
    fill([thresholds, fliplr(thresholds)], ...
         [bias_shifted_CI_lo, fliplr(bias_shifted_CI_hi)], ...
         [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    plot(thresholds, bias_shifted_mean, 'k-', 'LineWidth', 1.5);
    ylim([-0.15 0.35]); xlim([0 1]); yline(0, '--r');
    xlabel('V1 bias threshold'); ylabel('HPC bias diff (T1 - T2)');
    title(sprintf('Delta bin %d: %.2f–%.2f', npower, ...
           SO_thresholds(npower), SO_thresholds(npower+1)));
    set(gca,"TickDir","out",'box','off','Color','none','FontSize',12)

    % B: Proportion vs Bias Diff
    nexttile((npower-1)*3 + 2);
    hold on;
    valid = isfinite(bias_mean) & isfinite(prop_mean);
    fill([bias_CI_lo(valid), fliplr(bias_CI_hi(valid))], ...
         [prop_mean(valid), fliplr(prop_mean(valid))], ...
         colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    plot(bias_mean(valid), prop_mean(valid), '-', 'Color', colour_lines(npower,:), 'LineWidth', 2);
    fill([bias_shifted_CI_lo(valid), fliplr(bias_shifted_CI_hi(valid))], ...
         [prop_shifted_mean(valid), fliplr(prop_shifted_mean(valid))], ...
         [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    plot(bias_shifted_mean(valid), prop_shifted_mean(valid), 'k-', 'LineWidth', 1.5);
    xlim([-0.1 0.35]); xline(0, '--r');
    xlabel('HPC bias diff (T1 - T2)'); ylabel('Proportion of events detected');
    title('Event Proportion vs. Bias Difference');
    set(gca,"TickDir","out",'box','off','Color','none','FontSize',12)

    % C: Proportion CI vs Bias Diff
    nexttile((npower-1)*3 + 3);
    hold on;
    fill([bias_mean(valid), fliplr(bias_mean(valid))], ...
         [prop_CI_lo(valid), fliplr(prop_CI_hi(valid))], ...
         colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    plot(bias_mean(valid), prop_mean(valid), '-', 'Color', colour_lines(npower,:), 'LineWidth', 2);
    fill([bias_shifted_mean(valid), fliplr(bias_shifted_mean(valid))], ...
         [prop_shifted_CI_lo(valid), fliplr(prop_shifted_CI_hi(valid))], ...
         [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    plot(bias_shifted_mean(valid), prop_shifted_mean(valid), 'k-', 'LineWidth', 1.5);
    xlim([-0.1 0.35]); xline(0, '--r');
    xlabel('HPC bias diff (T1 - T2)'); ylabel('Proportion of events detected');
    title('Proportion vs. HPC Bias Difference');
    set(gca,"TickDir","out",'box','off','Color','none','FontSize',12)
    toc
end

clear Fill

fig = figure;
fig.Position = [640 100 2*1100/3 650/2];
fig.Name = 'KDE bias difference in HPC low vs high SO power';

% Select bins (1 = low, 4 = high)
nexttile;
for npower = [1 4]
    bias_mean = SO_power_KDE_bias_difference(npower).bias_diff_mean;
    bias_CI_lo = SO_power_KDE_bias_difference(npower).bias_diff_CI(1,:);
    bias_CI_hi = SO_power_KDE_bias_difference(npower).bias_diff_CI(2,:);
    thresholds = SO_power_KDE_bias_difference(npower).thresholds;

    hold on;
    x2 = [thresholds, fliplr(thresholds)];
    y2 = [bias_CI_lo, fliplr(bias_CI_hi)];
    Fill(npower) = fill(x2, y2, colour_lines(npower,:), ...
        'EdgeColor', 'none', 'FaceAlpha', 0.3);
    plot(thresholds, bias_mean, 'Color', colour_lines(npower,:), 'LineWidth', 2);

    xlabel('HPC Bias threshold');
    ylabel('V1 bias diff (T1 - T2)');
    set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12);
    ylim([-0.1 0.35]);
    % ylim([0 1])
end

bias_mean = SO_power_KDE_bias_difference(npower).bias_diff_shifted_mean;
bias_CI_lo = SO_power_KDE_bias_difference(npower).bias_diff_shifted_CI(1,:);
bias_CI_hi = SO_power_KDE_bias_difference(npower).bias_diff_shifted_CI(2,:)

hold on;
x2 = [thresholds, fliplr(thresholds)];
y2 = [bias_CI_lo, fliplr(bias_CI_hi)];
Fill(end +1) = fill(x2, y2, 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
plot(thresholds, bias_mean, 'Color','k', 'LineWidth', 2);

yline(0, '--r');
legend(Fill([1 4 5]), {'Low SO power', 'High SO power','Shuffled'}, 'box', 'off');

nexttile;
for npower = [1 4]
    bias_mean = SO_power_KDE_bias_difference(npower).bias_diff_mean;
    bias_CI_lo = SO_power_KDE_bias_difference(npower).bias_diff_CI(1,:);
    bias_CI_hi = SO_power_KDE_bias_difference(npower).bias_diff_CI(2,:);
    prop_mean = SO_power_KDE_bias_difference(npower).prop_mean;

    hold on;
    y2 = [prop_mean, fliplr(prop_mean)];
    x2 = [bias_CI_lo, fliplr(bias_CI_hi)];
    Fill(npower) = fill(x2, y2, colour_lines(npower,:), ...
        'EdgeColor', 'none', 'FaceAlpha', 0.3);
    plot(bias_mean, prop_mean, 'Color', colour_lines(npower,:), 'LineWidth', 2);

    xlabel('V1 bias diff (T1 - T2)');
    ylabel('Proportion of events detected');
    set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12);
    xlim([-0.1 0.35]);
    ylim([0 1])
end


bias_mean = SO_power_KDE_bias_difference(npower).bias_diff_shifted_mean;
bias_CI_lo = SO_power_KDE_bias_difference(npower).bias_diff_shifted_CI(1,:);
bias_CI_hi = SO_power_KDE_bias_difference(npower).bias_diff_shifted_CI(2,:)
prop_mean = SO_power_KDE_bias_difference(npower).prop_shifted_mean;

hold on;
y2 = [prop_mean, fliplr(prop_mean)];
x2 = [bias_CI_lo, fliplr(bias_CI_hi)];
Fill(end + 1) = fill(x2, y2, 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
plot(bias_mean, prop_mean, 'Color','k', 'LineWidth', 2);


xline(0, '--r');
legend(Fill([1 4 5]), {'Low SO power', 'High SO power','Shuffled'}, 'box', 'off');

% Save results
save('SO_power_KDE_bias_difference_based_on_V1_bias.mat', 'SO_power_KDE_bias_difference');





%%%%% AUC mean + CI bar plot
% Plot layout
fig = figure;
fig.Position = [640 100 281 325]
fig.Name = 'KDE bias V1 AUC low vs high SO power';
data = SO_power_KDE_bias_difference;
n_bins = length(data);
bar_width = 0.3;      % Width of the bars
group_offset = 0.15;    % Distance from the center integer (half the gap between bars)
hold on;
clear BAR
for i = 1:nBins
    % --- 1. Plot Shuffled Data (Left Bar) ---
    x_shuf = i - group_offset;
    y_shuf = data(i).AUC_mean_shuffled;

    % Calculate Error Deltas (Errorbar requires length relative to mean, not absolute values)
    % CI is [lower, upper]
    neg_err_shuf = y_shuf - data(i).AUC_CI_shuffled(1);
    pos_err_shuf = data(i).AUC_CI_shuffled(2) - y_shuf;

    % Plot Bar
    if i == 1
        BAR(1) = bar(x_shuf, y_shuf, bar_width, ...
            'FaceColor', 'k', ...
            'FaceAlpha', 0.15, ...
            'EdgeColor', 'none');
    else
        bar(x_shuf, y_shuf, bar_width, ...
            'FaceColor', 'k', ...
            'FaceAlpha', 0.15, ...
            'EdgeColor', 'none');
    end

    % Plot Error Bar
    E = errorbar(x_shuf, y_shuf, neg_err_shuf, pos_err_shuf, ...
        'Color', 'k', 'LineWidth', 1.5, 'CapSize', 8, 'LineStyle', 'none');


    % --- 2. Plot Real Data (Right Bar) ---
    x_real = i + group_offset;
    y_real = data(i).AUC_mean;

    neg_err_real = y_real - data(i).AUC_CI(1);
    pos_err_real = data(i).AUC_CI(2) - y_real;

    % Get specific color for this power bin
    % Assuming colour_lines is size [n_bins x 3]
    this_color = colour_lines(i, :);

    % Plot Bar
    BAR(i+1) = bar(x_real, y_real, bar_width, ...
        'FaceColor', this_color, ...
        'FaceAlpha', 0.3, ...
        'EdgeColor', 'none');

    errorbar(x_real, y_real, neg_err_real, pos_err_real, ...
        'Color', 'k', 'LineWidth', 1.5, 'CapSize', 8, 'LineStyle', 'none');

end

hold off;
% Set X-ticks to be centered on the groups
set(gca, 'XTick', 1:nBins);
xlim([0.5, nBins + 0.5]);

% Labels
ylabel('V1 bias AUC');
xlabel('Power Bins');
legend([BAR(1:end)],{'Shuffled','0-25','25-50','50-75','75-100'},'box','off')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE bias difference based on V1 bias'),[])



%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%% spindle phase binning
% Spindle power binning across both probes
% all_spindle_power = mean(ripple_info.spindle_amplitude, 1);  % avg of probe 1 and 2
% SO_thresholds = prctile(ripple_info.SO_phase, 0:99.9/4:99.9);
SO_thresholds = {'peak','trough'};
nBins = 2;

bins_to_use = bin_centers>0 & bin_centers<0.1;
bins_to_select = bin_centers>0 & bin_centers<0.2;
bins_to_use_shifted = bin_centers > -1 & bin_centers < -0.9;
nBoot = 1000;

SO_power_KDE_bias_difference = struct;

colour_lines = [ ...
    241, 182, 218;   % original end (lightest)
    231, 41, 138    % original start (darkest)
] / 256;

fig = figure;
fig.Position = [640 100 1100 650];
fig.Name = 'KDE bias difference in HPC with different spindle phase';
tiledlayout(nBins, 3, 'TileSpacing', 'compact');

for npower = 1:nBins
    tic
    % Select events in current spindle bin
    % power_index = ripple_info.spindle_amplitude(:,1) > spindle_thresholds(npower,1) & ...
    %               ripple_info.spindle_amplitude(:,1) <= spindle_thresholds(npower+1,1) |...
    %               ripple_info.spindle_amplitude(:,2) > spindle_thresholds(npower,2) & ...
    %               ripple_info.spindle_amplitude(:,2) <= spindle_thresholds(npower+1,2);

    event_index =1:length(z_bias);
    % event_index = find(power_index);

    mean_bias = mean(z_bias_V1(bins_to_select, event_index), 'omitnan');
    % mean_bias_shifted = mean(z_bias(bins_to_use_shifted, event_index), 'omitnan');
    mean_bias_V1 = mean(z_bias(bins_to_use, event_index), 'omitnan');
    selected_events = length(mean_bias);

    thresholds = prctile(abs(mean_bias), 0:10:100);
    thresholds = thresholds(1:end-1);
    nThresh = length(thresholds);

    bias_diff_boot = NaN(nBoot, nThresh);
    prop_events_boot = NaN(nBoot, nThresh);
    bias_diff_shifted_boot = NaN(nBoot, nThresh);
    prop_events_shifted_boot = NaN(nBoot, nThresh);


    % event_phase = ripple_info.SO_phase';


    parfor iBoot = 1:nBoot
        s = RandStream('philox4x32_10', 'Seed', iBoot);
        idx = randi(s, selected_events, selected_events, 1);

        true_idx = find(event_index);

        bb = mean_bias(idx);
        bb_shift = mean_bias;
        boot_V1 = mean_bias_V1(idx);

        diff_tmp = NaN(1, nThresh);
        prop_tmp = NaN(1, nThresh);
        diff_tmp_shifted = NaN(1, nThresh);
        prop_tmp_shifted = NaN(1, nThresh);
        event_phase = ripple_info.spindle_phase(true_idx(idx),:)';
        event_phase_shifted = ripple_info.spindle_phase(true_idx,:)';

        for i = 1:nThresh
            th = thresholds(i);

            t1 = bb >= th;
            t2 = bb <= -th;

            if npower == 1 % if phase peak


                t1 = t1 & (event_phase(2,:) >= -pi/2 & event_phase(2,:) <= pi/2);
                t2 = t2 & (event_phase(1,:) >= -pi/2 & event_phase(1,:) <= pi/2);

                t1_V1 = boot_V1(t1);
                t2_V1 = boot_V1(t2);
                if any(t1) && any(t2)
                    diff_tmp(i) = mean(t1_V1) - mean(t2_V1);
                end

                total_events = mean([sum(event_phase(2,:) >= -pi/2 & event_phase(2,:) <= pi/2) ...
                    sum(event_phase(1,:) >= -pi/2 & event_phase(1,:) <= pi/2)]);

                prop_tmp(i) = (sum(t1) + sum(t2)) / total_events;



                t1s = bb_shift >= th;
                t2s = bb_shift <= -th;
                t1s = t1s & (event_phase_shifted(2,:) >= -pi/2 & event_phase_shifted(2,:) <= pi/2);
                t2s = t2s & (event_phase_shifted(1,:) >= -pi/2 & event_phase_shifted(1,:) <= pi/2);

                t1_V1 = boot_V1(t1s);
                t2_V1 = boot_V1(t2s);

                if any(t1s) && any(t2s)
                    diff_tmp_shifted(i) = mean(t1_V1) - mean(t2_V1);
                end
                prop_tmp_shifted(i) = (sum(t1s) + sum(t2s)) / total_events;

            elseif npower == 2 % if phase trough

                t1 = t1 & (event_phase(2,:) >= -pi & event_phase(2,:) <= -pi/2 | event_phase(2,:) >= pi/2 & event_phase(2,:) <= pi);
                t2 = t2 & (event_phase(1,:) >= -pi & event_phase(1,:) <= -pi/2 | event_phase(1,:) >= pi/2 & event_phase(1,:) <= pi);

                t1_V1 = boot_V1(t1);
                t2_V1 = boot_V1(t2);
                if any(t1) && any(t2)
                    diff_tmp(i) = mean(t1_V1) - mean(t2_V1);
                end

                total_events = mean([sum(event_phase(2,:) >= -pi & event_phase(2,:) <= -pi/2 | event_phase(2,:) >= pi/2 & event_phase(2,:) <= pi) ...
                    sum(event_phase(1,:) >= -pi & event_phase(1,:) <= -pi/2 | event_phase(1,:) >= pi/2 & event_phase(1,:) <= pi)]);

                prop_tmp(i) = (sum(t1) + sum(t2)) / total_events;



                t1s = bb_shift >= th;
                t2s = bb_shift <= -th;
                t1s = t1s & (event_phase_shifted(2,:) >= -pi & event_phase_shifted(2,:) <= -pi/2 | event_phase_shifted(2,:) >= pi/2 & event_phase_shifted(2,:) <= pi);
                t2s = t2s & (event_phase_shifted(1,:) >= -pi & event_phase_shifted(1,:) <= -pi/2 | event_phase_shifted(1,:) >= pi/2 & event_phase_shifted(1,:) <= pi);

                t1_V1 = boot_V1(t1s);
                t2_V1 = boot_V1(t2s);

                if any(t1s) && any(t2s)
                    diff_tmp_shifted(i) = mean(t1_V1) - mean(t2_V1);
                end
                prop_tmp_shifted(i) = (sum(t1s) + sum(t2s)) / total_events;

            end



        end

        bias_diff_boot(iBoot, :) = diff_tmp;
        prop_events_boot(iBoot, :) = prop_tmp;
        bias_diff_shifted_boot(iBoot, :) = diff_tmp_shifted;
        prop_events_shifted_boot(iBoot, :) = prop_tmp_shifted;
    end

    % Stats
    bias_mean = mean(bias_diff_boot, 1, 'omitnan');
    bias_CI_lo = prctile(bias_diff_boot, 2.5, 1);
    bias_CI_hi = prctile(bias_diff_boot, 97.5, 1);
    prop_mean = mean(prop_events_boot, 1, 'omitnan');
    prop_CI_lo = prctile(prop_events_boot, 2.5, 1);
    prop_CI_hi = prctile(prop_events_boot, 97.5, 1);
    bias_shifted_mean = mean(bias_diff_shifted_boot, 1, 'omitnan');
    bias_shifted_CI_lo = prctile(bias_diff_shifted_boot, 2.5, 1);
    bias_shifted_CI_hi = prctile(bias_diff_shifted_boot, 97.5, 1);
    prop_shifted_mean = mean(prop_events_shifted_boot, 1, 'omitnan');
    prop_shifted_CI_lo = prctile(prop_events_shifted_boot, 2.5, 1);
    prop_shifted_CI_hi = prctile(prop_events_shifted_boot, 97.5, 1);

    % Store
    spindle_phase_KDE_bias_difference(npower).phase_range{npower} = ...
        SO_thresholds{npower};
    spindle_phase_KDE_bias_difference(npower).bias_diff_mean = bias_mean;
    spindle_phase_KDE_bias_difference(npower).bias_diff_CI = [bias_CI_lo; bias_CI_hi];
    spindle_phase_KDE_bias_difference(npower).prop_mean = prop_mean;
    spindle_phase_KDE_bias_difference(npower).prop_CI = [prop_CI_lo; prop_CI_hi];
    spindle_phase_KDE_bias_difference(npower).thresholds = thresholds;
    spindle_phase_KDE_bias_difference(npower).bias_diff_shifted_mean = bias_shifted_mean;
    spindle_phase_KDE_bias_difference(npower).bias_diff_shifted_CI = [bias_shifted_CI_lo; bias_shifted_CI_hi];
    spindle_phase_KDE_bias_difference(npower).prop_shifted_mean = prop_shifted_mean;
    spindle_phase_KDE_bias_difference(npower).prop_shifted_CI = [prop_shifted_CI_lo; prop_shifted_CI_hi];

    % ----------- PLOTS (A, B, C) ----------------

    % A: Bias vs. Threshold
    nexttile((npower-1)*3 + 1);
    hold on;
    fill([thresholds, fliplr(thresholds)], [bias_CI_lo, fliplr(bias_CI_hi)], ...
        colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    plot(thresholds, bias_mean, 'Color', colour_lines(npower,:), 'LineWidth', 2);
    fill([thresholds, fliplr(thresholds)], ...
         [bias_shifted_CI_lo, fliplr(bias_shifted_CI_hi)], ...
         [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    plot(thresholds, bias_shifted_mean, 'k-', 'LineWidth', 1.5);
    ylim([-0.15 0.35]); xlim([0 1]); yline(0, '--r');
    xlabel('V1 bias threshold'); ylabel('HPC bias diff (T1 - T2)');
    title(sprintf('spindle phase bin %s', SO_thresholds{npower}));
    set(gca,"TickDir","out",'box','off','Color','none','FontSize',12)

    % B: Proportion vs Bias Diff
    nexttile((npower-1)*3 + 2);
    hold on;
    valid = isfinite(bias_mean) & isfinite(prop_mean);
    fill([bias_CI_lo(valid), fliplr(bias_CI_hi(valid))], ...
         [prop_mean(valid), fliplr(prop_mean(valid))], ...
         colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    plot(bias_mean(valid), prop_mean(valid), '-', 'Color', colour_lines(npower,:), 'LineWidth', 2);
    fill([bias_shifted_CI_lo(valid), fliplr(bias_shifted_CI_hi(valid))], ...
         [prop_shifted_mean(valid), fliplr(prop_shifted_mean(valid))], ...
         [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    plot(bias_shifted_mean(valid), prop_shifted_mean(valid), 'k-', 'LineWidth', 1.5);
    xlim([-0.1 0.35]); xline(0, '--r');
    xlabel('HPC bias diff (T1 - T2)'); ylabel('Proportion of events detected');
    title('Event Proportion vs. Bias Difference');
    set(gca,"TickDir","out",'box','off','Color','none','FontSize',12)

    % C: Proportion CI vs Bias Diff
    nexttile((npower-1)*3 + 3);
    hold on;
    fill([bias_mean(valid), fliplr(bias_mean(valid))], ...
         [prop_CI_lo(valid), fliplr(prop_CI_hi(valid))], ...
         colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    plot(bias_mean(valid), prop_mean(valid), '-', 'Color', colour_lines(npower,:), 'LineWidth', 2);
    fill([bias_shifted_mean(valid), fliplr(bias_shifted_mean(valid))], ...
         [prop_shifted_CI_lo(valid), fliplr(prop_shifted_CI_hi(valid))], ...
         [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    plot(bias_shifted_mean(valid), prop_shifted_mean(valid), 'k-', 'LineWidth', 1.5);
    xlim([-0.1 0.35]); xline(0, '--r');
    xlabel('HPC bias diff (T1 - T2)'); ylabel('Proportion of events detected');
    title('Proportion vs. HPC Bias Difference');
    set(gca,"TickDir","out",'box','off','Color','none','FontSize',12)
    toc
end

clear Fill

fig = figure;
fig.Position = [640 100 2*1100/3 650/2];
% fig.Name = 'KDE bias difference in HPC low vs high SO power (PRE ripple)';
fig.Name = 'KDE bias difference in HPC spindle phase';

% Select bins (1 = low, 4 = high)
nexttile;
for npower = [1 2]
    bias_mean = spindle_phase_KDE_bias_difference(npower).bias_diff_mean;
    bias_CI_lo = spindle_phase_KDE_bias_difference(npower).bias_diff_CI(1,:);
    bias_CI_hi = spindle_phase_KDE_bias_difference(npower).bias_diff_CI(2,:);
    thresholds = spindle_phase_KDE_bias_difference(npower).thresholds;

    hold on;
    x2 = [thresholds, fliplr(thresholds)];
    y2 = [bias_CI_lo, fliplr(bias_CI_hi)];
    Fill(npower) = fill(x2, y2, colour_lines(npower,:), ...
        'EdgeColor', 'none', 'FaceAlpha', 0.3);
    plot(thresholds, bias_mean, 'Color', colour_lines(npower,:), 'LineWidth', 2);

    xlabel('HPC Bias threshold');
    ylabel('V1 bias diff (T1 - T2)');
    set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12);
    ylim([-0.1 0.35]);
    % ylim([0 1])
end

bias_mean = spindle_phase_KDE_bias_difference(npower).bias_diff_shifted_mean;
bias_CI_lo = spindle_phase_KDE_bias_difference(npower).bias_diff_shifted_CI(1,:);
bias_CI_hi = spindle_phase_KDE_bias_difference(npower).bias_diff_shifted_CI(2,:)

hold on;
x2 = [thresholds, fliplr(thresholds)];
y2 = [bias_CI_lo, fliplr(bias_CI_hi)];
Fill(end +1) = fill(x2, y2, 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
plot(thresholds, bias_mean, 'Color','k', 'LineWidth', 2);

yline(0, '--r');
legend(Fill([1 2 3]), {'spindle peak', 'spindle trough','Shuffled'}, 'box', 'off');

nexttile;
for npower = [1 2]
    bias_mean = spindle_phase_KDE_bias_difference(npower).bias_diff_mean;
    bias_CI_lo = spindle_phase_KDE_bias_difference(npower).bias_diff_CI(1,:);
    bias_CI_hi = spindle_phase_KDE_bias_difference(npower).bias_diff_CI(2,:);
    prop_mean = spindle_phase_KDE_bias_difference(npower).prop_mean;

    hold on;
    y2 = [prop_mean, fliplr(prop_mean)];
    x2 = [bias_CI_lo, fliplr(bias_CI_hi)];
    Fill(npower) = fill(x2, y2, colour_lines(npower,:), ...
        'EdgeColor', 'none', 'FaceAlpha', 0.3);
    plot(bias_mean, prop_mean, 'Color', colour_lines(npower,:), 'LineWidth', 2);

    xlabel('V1 bias diff (T1 - T2)');
    ylabel('Proportion of events detected');
    set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12);
    xlim([-0.1 0.35]);
    ylim([0 1])
end


bias_mean = spindle_phase_KDE_bias_difference(npower).bias_diff_shifted_mean;
bias_CI_lo = spindle_phase_KDE_bias_difference(npower).bias_diff_shifted_CI(1,:);
bias_CI_hi = spindle_phase_KDE_bias_difference(npower).bias_diff_shifted_CI(2,:)
prop_mean = spindle_phase_KDE_bias_difference(npower).prop_shifted_mean;

hold on;
y2 = [prop_mean, fliplr(prop_mean)];
x2 = [bias_CI_lo, fliplr(bias_CI_hi)];
Fill(end + 1) = fill(x2, y2, 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
plot(bias_mean, prop_mean, 'Color','k', 'LineWidth', 2);


xline(0, '--r');
legend(Fill([1 2 3]), {'spindle peak', 'spindle trough','Shuffled'}, 'box', 'off');

% Save results
save('spindle_phase_KDE_bias_difference_based_on_V1_bias.mat', 'spindle_phase_KDE_bias_difference');
save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE bias difference based on V1 bias'),[])








%%%%%%%%%%%%%%%%%%%%%%%% spindle phase binning
% Spindle power binning across both probes
% all_spindle_power = mean(ripple_info.spindle_amplitude, 1);  % avg of probe 1 and 2
% SO_thresholds = prctile(ripple_info.SO_phase, 0:99.9/4:99.9);
SO_thresholds = {'peak','trough'};
nBins = 2;

bins_to_use = bin_centers>0 & bin_centers<0.1;
bins_to_select = bin_centers>-0.2 & bin_centers<0;
bins_to_use_shifted = bin_centers > -1 & bin_centers < -0.9;
nBoot = 1000;

SO_power_KDE_bias_difference = struct;

colour_lines = [ ...
    241, 182, 218;   % original end (lightest)
    226, 132, 187;   % interpolated 2/3
    212,  78, 156;   % interpolated 1/3
    231, 41, 138    % original start (darkest)
] / 256;

fig = figure;
fig.Position = [640 100 1100 650];
fig.Name = 'KDE bias difference in HPC with different spindle phase (PRE ripple)';
tiledlayout(nBins, 3, 'TileSpacing', 'compact');

for npower = 1:nBins
    tic
    % Select events in current spindle bin
    % power_index = ripple_info.spindle_amplitude(:,1) > spindle_thresholds(npower,1) & ...
    %               ripple_info.spindle_amplitude(:,1) <= spindle_thresholds(npower+1,1) |...
    %               ripple_info.spindle_amplitude(:,2) > spindle_thresholds(npower,2) & ...
    %               ripple_info.spindle_amplitude(:,2) <= spindle_thresholds(npower+1,2);

    event_index =1:length(z_bias);
    % event_index = find(power_index);

    mean_bias = mean(z_bias_V1(bins_to_select, event_index), 'omitnan');
    % mean_bias_shifted = mean(z_bias(bins_to_use_shifted, event_index), 'omitnan');
    mean_bias_V1 = mean(z_bias(bins_to_use, event_index), 'omitnan');
    selected_events = length(mean_bias);

    thresholds = prctile(abs(mean_bias), 0:10:100);
    thresholds = thresholds(1:end-1);
    nThresh = length(thresholds);

    bias_diff_boot = NaN(nBoot, nThresh);
    prop_events_boot = NaN(nBoot, nThresh);
    bias_diff_shifted_boot = NaN(nBoot, nThresh);
    prop_events_shifted_boot = NaN(nBoot, nThresh);


    % event_phase = ripple_info.SO_phase';


    parfor iBoot = 1:nBoot
        s = RandStream('philox4x32_10', 'Seed', iBoot);
        idx = randi(s, selected_events, selected_events, 1);

        true_idx = find(event_index);

        bb = mean_bias(idx);
        bb_shift = mean_bias;
        boot_V1 = mean_bias_V1(idx);

        diff_tmp = NaN(1, nThresh);
        prop_tmp = NaN(1, nThresh);
        diff_tmp_shifted = NaN(1, nThresh);
        prop_tmp_shifted = NaN(1, nThresh);
        event_phase = ripple_info.spindle_phase(true_idx(idx),:)';
        event_phase_shifted = ripple_info.spindle_phase(true_idx,:)';

        for i = 1:nThresh
            th = thresholds(i);

            t1 = bb >= th;
            t2 = bb <= -th;

            if npower == 1 % if phase peak


                t1 = t1 & (event_phase(2,:) >= -pi/2 & event_phase(2,:) <= pi/2);
                t2 = t2 & (event_phase(1,:) >= -pi/2 & event_phase(1,:) <= pi/2);

                t1_V1 = boot_V1(t1);
                t2_V1 = boot_V1(t2);
                if any(t1) && any(t2)
                    diff_tmp(i) = mean(t1_V1) - mean(t2_V1);
                end

                total_events = mean([sum(event_phase(2,:) >= -pi/2 & event_phase(2,:) <= pi/2) ...
                    sum(event_phase(1,:) >= -pi/2 & event_phase(1,:) <= pi/2)]);

                prop_tmp(i) = (sum(t1) + sum(t2)) / total_events;



                t1s = bb_shift >= th;
                t2s = bb_shift <= -th;
                t1s = t1s & (event_phase_shifted(2,:) >= -pi/2 & event_phase_shifted(2,:) <= pi/2);
                t2s = t2s & (event_phase_shifted(1,:) >= -pi/2 & event_phase_shifted(1,:) <= pi/2);

                t1_V1 = boot_V1(t1s);
                t2_V1 = boot_V1(t2s);

                if any(t1s) && any(t2s)
                    diff_tmp_shifted(i) = mean(t1_V1) - mean(t2_V1);
                end
                prop_tmp_shifted(i) = (sum(t1s) + sum(t2s)) / total_events;

            elseif npower == 2 % if phase trough

                t1 = t1 & (event_phase(2,:) >= -pi & event_phase(2,:) <= -pi/2 | event_phase(2,:) >= pi/2 & event_phase(2,:) <= pi);
                t2 = t2 & (event_phase(1,:) >= -pi & event_phase(1,:) <= -pi/2 | event_phase(1,:) >= pi/2 & event_phase(1,:) <= pi);

                t1_V1 = boot_V1(t1);
                t2_V1 = boot_V1(t2);
                if any(t1) && any(t2)
                    diff_tmp(i) = mean(t1_V1) - mean(t2_V1);
                end

                total_events = mean([sum(event_phase(2,:) >= -pi & event_phase(2,:) <= -pi/2 | event_phase(2,:) >= pi/2 & event_phase(2,:) <= pi) ...
                    sum(event_phase(1,:) >= -pi & event_phase(1,:) <= -pi/2 | event_phase(1,:) >= pi/2 & event_phase(1,:) <= pi)]);

                prop_tmp(i) = (sum(t1) + sum(t2)) / total_events;



                t1s = bb_shift >= th;
                t2s = bb_shift <= -th;
                t1s = t1s & (event_phase_shifted(2,:) >= -pi & event_phase_shifted(2,:) <= -pi/2 | event_phase_shifted(2,:) >= pi/2 & event_phase_shifted(2,:) <= pi);
                t2s = t2s & (event_phase_shifted(1,:) >= -pi & event_phase_shifted(1,:) <= -pi/2 | event_phase_shifted(1,:) >= pi/2 & event_phase_shifted(1,:) <= pi);

                t1_V1 = boot_V1(t1s);
                t2_V1 = boot_V1(t2s);

                if any(t1s) && any(t2s)
                    diff_tmp_shifted(i) = mean(t1_V1) - mean(t2_V1);
                end
                prop_tmp_shifted(i) = (sum(t1s) + sum(t2s)) / total_events;

            end



        end

        bias_diff_boot(iBoot, :) = diff_tmp;
        prop_events_boot(iBoot, :) = prop_tmp;
        bias_diff_shifted_boot(iBoot, :) = diff_tmp_shifted;
        prop_events_shifted_boot(iBoot, :) = prop_tmp_shifted;
    end

    % Stats
    bias_mean = mean(bias_diff_boot, 1, 'omitnan');
    bias_CI_lo = prctile(bias_diff_boot, 2.5, 1);
    bias_CI_hi = prctile(bias_diff_boot, 97.5, 1);
    prop_mean = mean(prop_events_boot, 1, 'omitnan');
    prop_CI_lo = prctile(prop_events_boot, 2.5, 1);
    prop_CI_hi = prctile(prop_events_boot, 97.5, 1);
    bias_shifted_mean = mean(bias_diff_shifted_boot, 1, 'omitnan');
    bias_shifted_CI_lo = prctile(bias_diff_shifted_boot, 2.5, 1);
    bias_shifted_CI_hi = prctile(bias_diff_shifted_boot, 97.5, 1);
    prop_shifted_mean = mean(prop_events_shifted_boot, 1, 'omitnan');
    prop_shifted_CI_lo = prctile(prop_events_shifted_boot, 2.5, 1);
    prop_shifted_CI_hi = prctile(prop_events_shifted_boot, 97.5, 1);

    % Store
    spindle_phase_KDE_bias_difference(npower).phase_range{npower} = ...
        SO_thresholds{npower};
    spindle_phase_KDE_bias_difference(npower).bias_diff_mean = bias_mean;
    spindle_phase_KDE_bias_difference(npower).bias_diff_CI = [bias_CI_lo; bias_CI_hi];
    spindle_phase_KDE_bias_difference(npower).prop_mean = prop_mean;
    spindle_phase_KDE_bias_difference(npower).prop_CI = [prop_CI_lo; prop_CI_hi];
    spindle_phase_KDE_bias_difference(npower).thresholds = thresholds;
    spindle_phase_KDE_bias_difference(npower).bias_diff_shifted_mean = bias_shifted_mean;
    spindle_phase_KDE_bias_difference(npower).bias_diff_shifted_CI = [bias_shifted_CI_lo; bias_shifted_CI_hi];
    spindle_phase_KDE_bias_difference(npower).prop_shifted_mean = prop_shifted_mean;
    spindle_phase_KDE_bias_difference(npower).prop_shifted_CI = [prop_shifted_CI_lo; prop_shifted_CI_hi];

    % ----------- PLOTS (A, B, C) ----------------

    % A: Bias vs. Threshold
    nexttile((npower-1)*3 + 1);
    hold on;
    fill([thresholds, fliplr(thresholds)], [bias_CI_lo, fliplr(bias_CI_hi)], ...
        colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    plot(thresholds, bias_mean, 'Color', colour_lines(npower,:), 'LineWidth', 2);
    fill([thresholds, fliplr(thresholds)], ...
         [bias_shifted_CI_lo, fliplr(bias_shifted_CI_hi)], ...
         [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    plot(thresholds, bias_shifted_mean, 'k-', 'LineWidth', 1.5);
    ylim([-0.15 0.35]); xlim([0 1]); yline(0, '--r');
    xlabel('V1 bias threshold'); ylabel('HPC bias diff (T1 - T2)');
    title(sprintf('spindle phase bin %s', SO_thresholds{npower}));
    set(gca,"TickDir","out",'box','off','Color','none','FontSize',12)

    % B: Proportion vs Bias Diff
    nexttile((npower-1)*3 + 2);
    hold on;
    valid = isfinite(bias_mean) & isfinite(prop_mean);
    fill([bias_CI_lo(valid), fliplr(bias_CI_hi(valid))], ...
         [prop_mean(valid), fliplr(prop_mean(valid))], ...
         colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    plot(bias_mean(valid), prop_mean(valid), '-', 'Color', colour_lines(npower,:), 'LineWidth', 2);
    fill([bias_shifted_CI_lo(valid), fliplr(bias_shifted_CI_hi(valid))], ...
         [prop_shifted_mean(valid), fliplr(prop_shifted_mean(valid))], ...
         [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    plot(bias_shifted_mean(valid), prop_shifted_mean(valid), 'k-', 'LineWidth', 1.5);
    xlim([-0.1 0.35]); xline(0, '--r');
    xlabel('HPC bias diff (T1 - T2)'); ylabel('Proportion of events detected');
    title('Event Proportion vs. Bias Difference');
    set(gca,"TickDir","out",'box','off','Color','none','FontSize',12)

    % C: Proportion CI vs Bias Diff
    nexttile((npower-1)*3 + 3);
    hold on;
    fill([bias_mean(valid), fliplr(bias_mean(valid))], ...
         [prop_CI_lo(valid), fliplr(prop_CI_hi(valid))], ...
         colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    plot(bias_mean(valid), prop_mean(valid), '-', 'Color', colour_lines(npower,:), 'LineWidth', 2);
    fill([bias_shifted_mean(valid), fliplr(bias_shifted_mean(valid))], ...
         [prop_shifted_CI_lo(valid), fliplr(prop_shifted_CI_hi(valid))], ...
         [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    plot(bias_shifted_mean(valid), prop_shifted_mean(valid), 'k-', 'LineWidth', 1.5);
    xlim([-0.1 0.35]); xline(0, '--r');
    xlabel('HPC bias diff (T1 - T2)'); ylabel('Proportion of events detected');
    title('Proportion vs. HPC Bias Difference');
    set(gca,"TickDir","out",'box','off','Color','none','FontSize',12)
    toc
end

clear Fill

fig = figure;
fig.Position = [640 100 2*1100/3 650/2];
% fig.Name = 'KDE bias difference in HPC low vs high SO power (PRE ripple)';
fig.Name = 'KDE bias difference in HPC spindle phase (PRE ripple)';

% Select bins (1 = low, 4 = high)
nexttile;
for npower = [1 2]
    bias_mean = spindle_phase_KDE_bias_difference(npower).bias_diff_mean;
    bias_CI_lo = spindle_phase_KDE_bias_difference(npower).bias_diff_CI(1,:);
    bias_CI_hi = spindle_phase_KDE_bias_difference(npower).bias_diff_CI(2,:);
    thresholds = spindle_phase_KDE_bias_difference(npower).thresholds;

    hold on;
    x2 = [thresholds, fliplr(thresholds)];
    y2 = [bias_CI_lo, fliplr(bias_CI_hi)];
    Fill(npower) = fill(x2, y2, colour_lines(npower,:), ...
        'EdgeColor', 'none', 'FaceAlpha', 0.3);
    plot(thresholds, bias_mean, 'Color', colour_lines(npower,:), 'LineWidth', 2);

    xlabel('HPC Bias threshold');
    ylabel('V1 bias diff (T1 - T2)');
    set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12);
    ylim([-0.1 0.35]);
    % ylim([0 1])
end

bias_mean = spindle_phase_KDE_bias_difference(npower).bias_diff_shifted_mean;
bias_CI_lo = spindle_phase_KDE_bias_difference(npower).bias_diff_shifted_CI(1,:);
bias_CI_hi = spindle_phase_KDE_bias_difference(npower).bias_diff_shifted_CI(2,:)

hold on;
x2 = [thresholds, fliplr(thresholds)];
y2 = [bias_CI_lo, fliplr(bias_CI_hi)];
Fill(end +1) = fill(x2, y2, 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
plot(thresholds, bias_mean, 'Color','k', 'LineWidth', 2);

yline(0, '--r');
legend(Fill([1 2 3]), {'spindle peak', 'spindle trough','Shuffled'}, 'box', 'off');

nexttile;
for npower = [1 2]
    bias_mean = spindle_phase_KDE_bias_difference(npower).bias_diff_mean;
    bias_CI_lo = spindle_phase_KDE_bias_difference(npower).bias_diff_CI(1,:);
    bias_CI_hi = spindle_phase_KDE_bias_difference(npower).bias_diff_CI(2,:);
    prop_mean = spindle_phase_KDE_bias_difference(npower).prop_mean;

    hold on;
    y2 = [prop_mean, fliplr(prop_mean)];
    x2 = [bias_CI_lo, fliplr(bias_CI_hi)];
    Fill(npower) = fill(x2, y2, colour_lines(npower,:), ...
        'EdgeColor', 'none', 'FaceAlpha', 0.3);
    plot(bias_mean, prop_mean, 'Color', colour_lines(npower,:), 'LineWidth', 2);

    xlabel('V1 bias diff (T1 - T2)');
    ylabel('Proportion of events detected');
    set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12);
    xlim([-0.1 0.35]);
    ylim([0 1])
end


bias_mean = spindle_phase_KDE_bias_difference(npower).bias_diff_shifted_mean;
bias_CI_lo = spindle_phase_KDE_bias_difference(npower).bias_diff_shifted_CI(1,:);
bias_CI_hi = spindle_phase_KDE_bias_difference(npower).bias_diff_shifted_CI(2,:)
prop_mean = spindle_phase_KDE_bias_difference(npower).prop_shifted_mean;

hold on;
y2 = [prop_mean, fliplr(prop_mean)];
x2 = [bias_CI_lo, fliplr(bias_CI_hi)];
Fill(end + 1) = fill(x2, y2, 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
plot(bias_mean, prop_mean, 'Color','k', 'LineWidth', 2);


xline(0, '--r');
legend(Fill([1 2 3]), {'spindle peak', 'spindle trough','Shuffled'}, 'box', 'off');

% Save results
save('spindle_phase_KDE_bias_difference_based_on_PRE_V1_bias.mat', 'spindle_phase_KDE_bias_difference');
save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE bias difference based on V1 bias'),[])




%% SO
%%%%%%%%%%%%%%%%%%%%%%%% delta power binning
% Spindle power binning across both probes
% all_spindle_power = mean(ripple_info.spindle_amplitude, 1);  % avg of probe 1 and 2
SO_thresholds = prctile(ripple_info.SO_amplitude, 0:99.9/4:99.9);
nBins = length(SO_thresholds) - 1;

bins_to_use = bin_centers>0 & bin_centers<0.1;
bins_to_select = bin_centers>-0.2 & bin_centers<0;
bins_to_use_shifted = bin_centers > -1 & bin_centers < -0.9;
nBoot = 1000;

SO_power_KDE_bias_difference = struct;

colour_lines = [ ...
    241, 182, 218;   % original end (lightest)
    226, 132, 187;   % interpolated 2/3
    212,  78, 156;   % interpolated 1/3
    231, 41, 138    % original start (darkest)
] / 256;



fig = figure;
fig.Position = [640 100 1100 650*2];
fig.Name = 'KDE bias difference in HPC with different SO powers (PRE ripple)';
tiledlayout(nBins, 3, 'TileSpacing', 'compact');

for npower = 1:nBins
    tic
    % Select events in current spindle bin
    % power_index = ripple_info.spindle_amplitude(:,1) > spindle_thresholds(npower,1) & ...
    %               ripple_info.spindle_amplitude(:,1) <= spindle_thresholds(npower+1,1) |...
    %               ripple_info.spindle_amplitude(:,2) > spindle_thresholds(npower,2) & ...
    %               ripple_info.spindle_amplitude(:,2) <= spindle_thresholds(npower+1,2);

    event_index =1:length(z_bias);
    % event_index = find(power_index);

    mean_bias = mean(z_bias_V1(bins_to_select, event_index), 'omitnan');
    % mean_bias_shifted = mean(z_bias(bins_to_use_shifted, event_index), 'omitnan');
    mean_bias_V1 = mean(z_bias(bins_to_use, event_index), 'omitnan');
    selected_events = length(mean_bias);

    thresholds = prctile(abs(mean_bias), 0:10:100);
    thresholds = thresholds(1:end-1);
    nThresh = length(thresholds);

    bias_diff_boot = NaN(nBoot, nThresh);
    prop_events_boot = NaN(nBoot, nThresh);
    bias_diff_shifted_boot = NaN(nBoot, nThresh);
    prop_events_shifted_boot = NaN(nBoot, nThresh);


    parfor iBoot = 1:nBoot
        s = RandStream('philox4x32_10', 'Seed', iBoot);
        idx = randi(s, selected_events, selected_events, 1);

        true_idx = find(event_index);

        bb = mean_bias(idx);
        bb_shift = mean_bias;
        boot_V1 = mean_bias_V1(idx);

        diff_tmp = NaN(1, nThresh);
        prop_tmp = NaN(1, nThresh);
        diff_tmp_shifted = NaN(1, nThresh);
        prop_tmp_shifted = NaN(1, nThresh);

        for i = 1:nThresh
            th = thresholds(i);

            t1 = bb >= th;
            t2 = bb <= -th;



            t1 = t1' + (ripple_info.SO_amplitude( true_idx(idx),2) > SO_thresholds(npower,2) &...
                ripple_info.SO_amplitude( true_idx(idx),2) <= SO_thresholds(npower+1,2))>1';
            t2 = t2' + (ripple_info.SO_amplitude( true_idx(idx),1) > SO_thresholds(npower,1) &...
                ripple_info.SO_amplitude( true_idx(idx),1) <= SO_thresholds(npower+1,1))>1';

            t1_V1 = boot_V1(t1);
            t2_V1 = boot_V1(t2);
            if any(t1) && any(t2)
                diff_tmp(i) = mean(t1_V1) - mean(t2_V1);
            end

            total_events = mean([sum((ripple_info.SO_amplitude( true_idx(idx),2) > SO_thresholds(npower,2) &...
                ripple_info.SO_amplitude( true_idx(idx),2) <= SO_thresholds(npower+1,2))) ...
                sum((ripple_info.SO_amplitude( true_idx(idx),1) > SO_thresholds(npower,1) &...
                ripple_info.SO_amplitude( true_idx(idx),1) <= SO_thresholds(npower+1,1)))]);

            prop_tmp(i) = (sum(t1) + sum(t2)) / total_events;



            t1s = bb_shift >= th;
            t2s = bb_shift <= -th;
            t1s = t1s' + (ripple_info.SO_amplitude( true_idx,2) > SO_thresholds(npower,2) &...
                ripple_info.SO_amplitude( true_idx,2) <= SO_thresholds(npower+1,2))>1';
            t2s = t2s' + (ripple_info.SO_amplitude( true_idx,1) > SO_thresholds(npower,1) &...
                ripple_info.SO_amplitude( true_idx,1) <= SO_thresholds(npower+1,1))>1';
            t1_V1 = boot_V1(t1s);
            t2_V1 = boot_V1(t2s);

            if any(t1s) && any(t2s)
                diff_tmp_shifted(i) = mean(t1_V1) - mean(t2_V1);
            end
            prop_tmp_shifted(i) = (sum(t1s) + sum(t2s)) / total_events;
        end

        bias_diff_boot(iBoot, :) = diff_tmp;
        prop_events_boot(iBoot, :) = prop_tmp;
        bias_diff_shifted_boot(iBoot, :) = diff_tmp_shifted;
        prop_events_shifted_boot(iBoot, :) = prop_tmp_shifted;
    end

    % Stats
    bias_mean = mean(bias_diff_boot, 1, 'omitnan');
    bias_CI_lo = prctile(bias_diff_boot, 2.5, 1);
    bias_CI_hi = prctile(bias_diff_boot, 97.5, 1);
    prop_mean = mean(prop_events_boot, 1, 'omitnan');
    prop_CI_lo = prctile(prop_events_boot, 2.5, 1);
    prop_CI_hi = prctile(prop_events_boot, 97.5, 1);
    bias_shifted_mean = mean(bias_diff_shifted_boot, 1, 'omitnan');
    bias_shifted_CI_lo = prctile(bias_diff_shifted_boot, 2.5, 1);
    bias_shifted_CI_hi = prctile(bias_diff_shifted_boot, 97.5, 1);
    prop_shifted_mean = mean(prop_events_shifted_boot, 1, 'omitnan');
    prop_shifted_CI_lo = prctile(prop_events_shifted_boot, 2.5, 1);
    prop_shifted_CI_hi = prctile(prop_events_shifted_boot, 97.5, 1);

    % Store
    SO_power_KDE_bias_difference(npower).power_range = ...
        [SO_thresholds(npower), SO_thresholds(npower+1)];
    SO_power_KDE_bias_difference(npower).bias_diff_mean = bias_mean;
    SO_power_KDE_bias_difference(npower).bias_diff_CI = [bias_CI_lo; bias_CI_hi];
    SO_power_KDE_bias_difference(npower).prop_mean = prop_mean;
    SO_power_KDE_bias_difference(npower).prop_CI = [prop_CI_lo; prop_CI_hi];
    SO_power_KDE_bias_difference(npower).thresholds = thresholds;
    SO_power_KDE_bias_difference(npower).bias_diff_shifted_mean = bias_shifted_mean;
    SO_power_KDE_bias_difference(npower).bias_diff_shifted_CI = [bias_shifted_CI_lo; bias_shifted_CI_hi];
    SO_power_KDE_bias_difference(npower).prop_shifted_mean = prop_shifted_mean;
    SO_power_KDE_bias_difference(npower).prop_shifted_CI = [prop_shifted_CI_lo; prop_shifted_CI_hi];

    % ----------- PLOTS (A, B, C) ----------------

    % A: Bias vs. Threshold
    nexttile((npower-1)*3 + 1);
    hold on;
    fill([thresholds, fliplr(thresholds)], [bias_CI_lo, fliplr(bias_CI_hi)], ...
        colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    plot(thresholds, bias_mean, 'Color', colour_lines(npower,:), 'LineWidth', 2);
    fill([thresholds, fliplr(thresholds)], ...
         [bias_shifted_CI_lo, fliplr(bias_shifted_CI_hi)], ...
         [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    plot(thresholds, bias_shifted_mean, 'k-', 'LineWidth', 1.5);
    ylim([-0.15 0.35]); xlim([0 1]); yline(0, '--r');
    xlabel('V1 bias threshold'); ylabel('HPC bias diff (T1 - T2)');
    title(sprintf('Delta bin %d: %.2f–%.2f', npower, ...
           SO_thresholds(npower), SO_thresholds(npower+1)));
    set(gca,"TickDir","out",'box','off','Color','none','FontSize',12)

    % B: Proportion vs Bias Diff
    nexttile((npower-1)*3 + 2);
    hold on;
    valid = isfinite(bias_mean) & isfinite(prop_mean);
    fill([bias_CI_lo(valid), fliplr(bias_CI_hi(valid))], ...
         [prop_mean(valid), fliplr(prop_mean(valid))], ...
         colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    plot(bias_mean(valid), prop_mean(valid), '-', 'Color', colour_lines(npower,:), 'LineWidth', 2);
    fill([bias_shifted_CI_lo(valid), fliplr(bias_shifted_CI_hi(valid))], ...
         [prop_shifted_mean(valid), fliplr(prop_shifted_mean(valid))], ...
         [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    plot(bias_shifted_mean(valid), prop_shifted_mean(valid), 'k-', 'LineWidth', 1.5);
    xlim([-0.1 0.35]); xline(0, '--r');
    xlabel('HPC bias diff (T1 - T2)'); ylabel('Proportion of events detected');
    title('Event Proportion vs. Bias Difference');
    set(gca,"TickDir","out",'box','off','Color','none','FontSize',12)

    % C: Proportion CI vs Bias Diff
    nexttile((npower-1)*3 + 3);
    hold on;
    fill([bias_mean(valid), fliplr(bias_mean(valid))], ...
         [prop_CI_lo(valid), fliplr(prop_CI_hi(valid))], ...
         colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    plot(bias_mean(valid), prop_mean(valid), '-', 'Color', colour_lines(npower,:), 'LineWidth', 2);
    fill([bias_shifted_mean(valid), fliplr(bias_shifted_mean(valid))], ...
         [prop_shifted_CI_lo(valid), fliplr(prop_shifted_CI_hi(valid))], ...
         [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    plot(bias_shifted_mean(valid), prop_shifted_mean(valid), 'k-', 'LineWidth', 1.5);
    xlim([-0.1 0.35]); xline(0, '--r');
    xlabel('HPC bias diff (T1 - T2)'); ylabel('Proportion of events detected');
    title('Proportion vs. HPC Bias Difference');
    set(gca,"TickDir","out",'box','off','Color','none','FontSize',12)
    toc
end

clear Fill

fig = figure;
fig.Position = [640 100 2*1100/3 650/2];
fig.Name = 'KDE bias difference in HPC low vs high SO power (PRE ripple)';

% Select bins (1 = low, 4 = high)
nexttile;
for npower = [1 4]
    bias_mean = SO_power_KDE_bias_difference(npower).bias_diff_mean;
    bias_CI_lo = SO_power_KDE_bias_difference(npower).bias_diff_CI(1,:);
    bias_CI_hi = SO_power_KDE_bias_difference(npower).bias_diff_CI(2,:);
    thresholds = SO_power_KDE_bias_difference(npower).thresholds;

    hold on;
    x2 = [thresholds, fliplr(thresholds)];
    y2 = [bias_CI_lo, fliplr(bias_CI_hi)];
    Fill(npower) = fill(x2, y2, colour_lines(npower,:), ...
        'EdgeColor', 'none', 'FaceAlpha', 0.3);
    plot(thresholds, bias_mean, 'Color', colour_lines(npower,:), 'LineWidth', 2);

    xlabel('HPC Bias threshold');
    ylabel('V1 bias diff (T1 - T2)');
    set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12);
    ylim([-0.1 0.35]);
    % ylim([0 1])
end

bias_mean = SO_power_KDE_bias_difference(npower).bias_diff_shifted_mean;
bias_CI_lo = SO_power_KDE_bias_difference(npower).bias_diff_shifted_CI(1,:);
bias_CI_hi = SO_power_KDE_bias_difference(npower).bias_diff_shifted_CI(2,:)

hold on;
x2 = [thresholds, fliplr(thresholds)];
y2 = [bias_CI_lo, fliplr(bias_CI_hi)];
Fill(end +1) = fill(x2, y2, 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
plot(thresholds, bias_mean, 'Color','k', 'LineWidth', 2);

yline(0, '--r');
legend(Fill([1 4 5]), {'Low SO power', 'High SO power','Shuffled'}, 'box', 'off');

nexttile;
for npower = [1 4]
    bias_mean = SO_power_KDE_bias_difference(npower).bias_diff_mean;
    bias_CI_lo = SO_power_KDE_bias_difference(npower).bias_diff_CI(1,:);
    bias_CI_hi = SO_power_KDE_bias_difference(npower).bias_diff_CI(2,:);
    prop_mean = SO_power_KDE_bias_difference(npower).prop_mean;

    hold on;
    y2 = [prop_mean, fliplr(prop_mean)];
    x2 = [bias_CI_lo, fliplr(bias_CI_hi)];
    Fill(npower) = fill(x2, y2, colour_lines(npower,:), ...
        'EdgeColor', 'none', 'FaceAlpha', 0.3);
    plot(bias_mean, prop_mean, 'Color', colour_lines(npower,:), 'LineWidth', 2);

    xlabel('V1 bias diff (T1 - T2)');
    ylabel('Proportion of events detected');
    set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12);
    xlim([-0.1 0.35]);
    ylim([0 1])
end


bias_mean = SO_power_KDE_bias_difference(npower).bias_diff_shifted_mean;
bias_CI_lo = SO_power_KDE_bias_difference(npower).bias_diff_shifted_CI(1,:);
bias_CI_hi = SO_power_KDE_bias_difference(npower).bias_diff_shifted_CI(2,:)
prop_mean = SO_power_KDE_bias_difference(npower).prop_shifted_mean;

hold on;
y2 = [prop_mean, fliplr(prop_mean)];
x2 = [bias_CI_lo, fliplr(bias_CI_hi)];
Fill(end + 1) = fill(x2, y2, 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
plot(bias_mean, prop_mean, 'Color','k', 'LineWidth', 2);


xline(0, '--r');
legend(Fill([1 4 5]), {'Low SO power', 'High SO power','Shuffled'}, 'box', 'off');

% Save results
save('SO_power_KDE_bias_difference_based_on_PRE_V1_bias.mat', 'SO_power_KDE_bias_difference');
save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE bias difference based on V1 bias'),[])




%% SO phases
%%%%%%%%%%%%%%%%
% SO phase binning across both probes
% all_spindle_power = mean(ripple_info.spindle_amplitude, 1);  % avg of probe 1 and 2
% SO_thresholds = prctile(ripple_info.SO_phase, 0:99.9/4:99.9);
SO_thresholds = {'peak','trough'};
nBins = 2;

bins_to_use = bin_centers>0 & bin_centers<0.1;
bins_to_select = bin_centers>0 & bin_centers<0.2;
bins_to_use_shifted = bin_centers > -1 & bin_centers < -0.9;
nBoot = 1000;

SO_power_KDE_bias_difference = struct;

colour_lines = [ ...
    241, 182, 218;   % original end (lightest)
    231, 41, 138    % original start (darkest)
] / 256;

fig = figure;
fig.Position = [640 100 1100 650];
fig.Name = 'KDE bias difference in HPC with different SO phase';
tiledlayout(nBins, 3, 'TileSpacing', 'compact');

for npower = 1:nBins
    tic
    % Select events in current spindle bin
    % power_index = ripple_info.spindle_amplitude(:,1) > spindle_thresholds(npower,1) & ...
    %               ripple_info.spindle_amplitude(:,1) <= spindle_thresholds(npower+1,1) |...
    %               ripple_info.spindle_amplitude(:,2) > spindle_thresholds(npower,2) & ...
    %               ripple_info.spindle_amplitude(:,2) <= spindle_thresholds(npower+1,2);

    event_index =1:length(z_bias);
    % event_index = find(power_index);

    mean_bias = mean(z_bias_V1(bins_to_select, event_index), 'omitnan');
    % mean_bias_shifted = mean(z_bias(bins_to_use_shifted, event_index), 'omitnan');
    mean_bias_V1 = mean(z_bias(bins_to_use, event_index), 'omitnan');
    selected_events = length(mean_bias);

    thresholds = prctile(abs(mean_bias), 0:10:100);
    thresholds = thresholds(1:end-1);
    nThresh = length(thresholds);

    bias_diff_boot = NaN(nBoot, nThresh);
    prop_events_boot = NaN(nBoot, nThresh);
    bias_diff_shifted_boot = NaN(nBoot, nThresh);
    prop_events_shifted_boot = NaN(nBoot, nThresh);


    parfor iBoot = 1:nBoot
        s = RandStream('philox4x32_10', 'Seed', iBoot);
        idx = randi(s, selected_events, selected_events, 1);

        true_idx = find(event_index);

        bb = mean_bias(idx);
        bb_shift = mean_bias;
        boot_V1 = mean_bias_V1(idx);

        diff_tmp = NaN(1, nThresh);
        prop_tmp = NaN(1, nThresh);
        diff_tmp_shifted = NaN(1, nThresh);
        prop_tmp_shifted = NaN(1, nThresh);
        event_phase = ripple_info.SO_phase(true_idx(idx),:)';
        event_phase_shifted = ripple_info.SO_phase(true_idx,:)';

        for i = 1:nThresh
            th = thresholds(i);

            t1 = bb >= th;
            t2 = bb <= -th;

            if npower == 1 % if phase peak


                t1 = t1 & (event_phase(2,:) >= -pi/2 & event_phase(2,:) <= pi/2);
                t2 = t2 & (event_phase(1,:) >= -pi/2 & event_phase(1,:) <= pi/2);

                t1_V1 = boot_V1(t1);
                t2_V1 = boot_V1(t2);
                if any(t1) && any(t2)
                    diff_tmp(i) = mean(t1_V1) - mean(t2_V1);
                end

                total_events = mean([sum(event_phase(2,:) >= -pi/2 & event_phase(2,:) <= pi/2) ...
                    sum(event_phase(1,:) >= -pi/2 & event_phase(1,:) <= pi/2)]);

                prop_tmp(i) = (sum(t1) + sum(t2)) / total_events;



                t1s = bb_shift >= th;
                t2s = bb_shift <= -th;
                t1s = t1s & (event_phase_shifted(2,:) >= -pi/2 & event_phase_shifted(2,:) <= pi/2);
                t2s = t2s & (event_phase_shifted(1,:) >= -pi/2 & event_phase_shifted(1,:) <= pi/2);

                t1_V1 = boot_V1(t1s);
                t2_V1 = boot_V1(t2s);

                if any(t1s) && any(t2s)
                    diff_tmp_shifted(i) = mean(t1_V1) - mean(t2_V1);
                end
                prop_tmp_shifted(i) = (sum(t1s) + sum(t2s)) / total_events;

            elseif npower == 2 % if phase trough

                t1 = t1 & (event_phase(2,:) >= -pi & event_phase(2,:) <= -pi/2 | event_phase(2,:) >= pi/2 & event_phase(2,:) <= pi);
                t2 = t2 & (event_phase(1,:) >= -pi & event_phase(1,:) <= -pi/2 | event_phase(1,:) >= pi/2 & event_phase(1,:) <= pi);

                t1_V1 = boot_V1(t1);
                t2_V1 = boot_V1(t2);
                if any(t1) && any(t2)
                    diff_tmp(i) = mean(t1_V1) - mean(t2_V1);
                end

                total_events = mean([sum(event_phase(2,:) >= -pi & event_phase(2,:) <= -pi/2 | event_phase(2,:) >= pi/2 & event_phase(2,:) <= pi) ...
                    sum(event_phase(1,:) >= -pi & event_phase(1,:) <= -pi/2 | event_phase(1,:) >= pi/2 & event_phase(1,:) <= pi)]);

                prop_tmp(i) = (sum(t1) + sum(t2)) / total_events;



                t1s = bb_shift >= th;
                t2s = bb_shift <= -th;
                t1s = t1s & (event_phase_shifted(2,:) >= -pi & event_phase_shifted(2,:) <= -pi/2 | event_phase_shifted(2,:) >= pi/2 & event_phase_shifted(2,:) <= pi);
                t2s = t2s & (event_phase_shifted(1,:) >= -pi & event_phase_shifted(1,:) <= -pi/2 | event_phase_shifted(1,:) >= pi/2 & event_phase_shifted(1,:) <= pi);

                t1_V1 = boot_V1(t1s);
                t2_V1 = boot_V1(t2s);

                if any(t1s) && any(t2s)
                    diff_tmp_shifted(i) = mean(t1_V1) - mean(t2_V1);
                end
                prop_tmp_shifted(i) = (sum(t1s) + sum(t2s)) / total_events;

            end



        end

        bias_diff_boot(iBoot, :) = diff_tmp;
        prop_events_boot(iBoot, :) = prop_tmp;
        bias_diff_shifted_boot(iBoot, :) = diff_tmp_shifted;
        prop_events_shifted_boot(iBoot, :) = prop_tmp_shifted;
    end

    % Stats
    bias_mean = mean(bias_diff_boot, 1, 'omitnan');
    bias_CI_lo = prctile(bias_diff_boot, 2.5, 1);
    bias_CI_hi = prctile(bias_diff_boot, 97.5, 1);
    prop_mean = mean(prop_events_boot, 1, 'omitnan');
    prop_CI_lo = prctile(prop_events_boot, 2.5, 1);
    prop_CI_hi = prctile(prop_events_boot, 97.5, 1);
    bias_shifted_mean = mean(bias_diff_shifted_boot, 1, 'omitnan');
    bias_shifted_CI_lo = prctile(bias_diff_shifted_boot, 2.5, 1);
    bias_shifted_CI_hi = prctile(bias_diff_shifted_boot, 97.5, 1);
    prop_shifted_mean = mean(prop_events_shifted_boot, 1, 'omitnan');
    prop_shifted_CI_lo = prctile(prop_events_shifted_boot, 2.5, 1);
    prop_shifted_CI_hi = prctile(prop_events_shifted_boot, 97.5, 1);

    % Store
    SO_phase_KDE_bias_difference(npower).phase_range{npower} = ...
        SO_thresholds{npower};
    SO_phase_KDE_bias_difference(npower).bias_diff_mean = bias_mean;
    SO_phase_KDE_bias_difference(npower).bias_diff_CI = [bias_CI_lo; bias_CI_hi];
    SO_phase_KDE_bias_difference(npower).prop_mean = prop_mean;
    SO_phase_KDE_bias_difference(npower).prop_CI = [prop_CI_lo; prop_CI_hi];
    SO_phase_KDE_bias_difference(npower).thresholds = thresholds;
    SO_phase_KDE_bias_difference(npower).bias_diff_shifted_mean = bias_shifted_mean;
    SO_phase_KDE_bias_difference(npower).bias_diff_shifted_CI = [bias_shifted_CI_lo; bias_shifted_CI_hi];
    SO_phase_KDE_bias_difference(npower).prop_shifted_mean = prop_shifted_mean;
    SO_phase_KDE_bias_difference(npower).prop_shifted_CI = [prop_shifted_CI_lo; prop_shifted_CI_hi];

    % store AUC
    auc_boot = (trapz(thresholds, bias_diff_boot') / (max(thresholds)-min(thresholds)))';
    auc_shift_boot = (trapz(thresholds, bias_diff_shifted_boot') / (max(thresholds)-min(thresholds)))';

    SO_phase_KDE_bias_difference(npower).AUC_mean = mean(auc_boot, 'omitnan');
    SO_phase_KDE_bias_difference(npower).AUC_CI = prctile(auc_boot, [2.5 97.5]);
    SO_phase_KDE_bias_difference(npower).AUC_mean_shuffled = mean(auc_shift_boot, 'omitnan');
    SO_phase_KDE_bias_difference(npower).AUC_CI_shuffled = prctile(auc_shift_boot, [2.5 97.5]);


    % ----------- PLOTS (A, B, C) ----------------

    % A: Bias vs. Threshold
    nexttile((npower-1)*3 + 1);
    hold on;
    fill([thresholds, fliplr(thresholds)], [bias_CI_lo, fliplr(bias_CI_hi)], ...
        colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    plot(thresholds, bias_mean, 'Color', colour_lines(npower,:), 'LineWidth', 2);
    fill([thresholds, fliplr(thresholds)], ...
         [bias_shifted_CI_lo, fliplr(bias_shifted_CI_hi)], ...
         [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    plot(thresholds, bias_shifted_mean, 'k-', 'LineWidth', 1.5);
    ylim([-0.15 0.35]); xlim([0 1]); yline(0, '--r');
    xlabel('V1 bias threshold'); ylabel('HPC bias diff (T1 - T2)');
    title(sprintf('SO phase bin %s', SO_thresholds{npower}));
    set(gca,"TickDir","out",'box','off','Color','none','FontSize',12)

    % B: Proportion vs Bias Diff
    nexttile((npower-1)*3 + 2);
    hold on;
    valid = isfinite(bias_mean) & isfinite(prop_mean);
    fill([bias_CI_lo(valid), fliplr(bias_CI_hi(valid))], ...
         [prop_mean(valid), fliplr(prop_mean(valid))], ...
         colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    plot(bias_mean(valid), prop_mean(valid), '-', 'Color', colour_lines(npower,:), 'LineWidth', 2);
    fill([bias_shifted_CI_lo(valid), fliplr(bias_shifted_CI_hi(valid))], ...
         [prop_shifted_mean(valid), fliplr(prop_shifted_mean(valid))], ...
         [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    plot(bias_shifted_mean(valid), prop_shifted_mean(valid), 'k-', 'LineWidth', 1.5);
    xlim([-0.1 0.35]); xline(0, '--r');
    xlabel('HPC bias diff (T1 - T2)'); ylabel('Proportion of events detected');
    title('Event Proportion vs. Bias Difference');
    set(gca,"TickDir","out",'box','off','Color','none','FontSize',12)

    % C: Proportion CI vs Bias Diff
    nexttile((npower-1)*3 + 3);
    hold on;
    fill([bias_mean(valid), fliplr(bias_mean(valid))], ...
         [prop_CI_lo(valid), fliplr(prop_CI_hi(valid))], ...
         colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    plot(bias_mean(valid), prop_mean(valid), '-', 'Color', colour_lines(npower,:), 'LineWidth', 2);
    fill([bias_shifted_mean(valid), fliplr(bias_shifted_mean(valid))], ...
         [prop_shifted_CI_lo(valid), fliplr(prop_shifted_CI_hi(valid))], ...
         [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    plot(bias_shifted_mean(valid), prop_shifted_mean(valid), 'k-', 'LineWidth', 1.5);
    xlim([-0.1 0.35]); xline(0, '--r');
    xlabel('HPC bias diff (T1 - T2)'); ylabel('Proportion of events detected');
    title('Proportion vs. HPC Bias Difference');
    set(gca,"TickDir","out",'box','off','Color','none','FontSize',12)
    toc
end

clear Fill

fig = figure;
fig.Position = [640 100 2*1100/3 650/2];
fig.Name = 'KDE bias difference in HPC SO phase';

% Select bins (1 = low, 4 = high)
nexttile;
for npower = [1 2]
    bias_mean = SO_phase_KDE_bias_difference(npower).bias_diff_mean;
    bias_CI_lo = SO_phase_KDE_bias_difference(npower).bias_diff_CI(1,:);
    bias_CI_hi = SO_phase_KDE_bias_difference(npower).bias_diff_CI(2,:);
    thresholds = SO_phase_KDE_bias_difference(npower).thresholds;

    hold on;
    x2 = [thresholds, fliplr(thresholds)];
    y2 = [bias_CI_lo, fliplr(bias_CI_hi)];
    Fill(npower) = fill(x2, y2, colour_lines(npower,:), ...
        'EdgeColor', 'none', 'FaceAlpha', 0.3);
    plot(thresholds, bias_mean, 'Color', colour_lines(npower,:), 'LineWidth', 2);

    xlabel('HPC Bias threshold');
    ylabel('V1 bias diff (T1 - T2)');
    set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12);
    ylim([-0.1 0.35]);
    % ylim([0 1])
end

bias_mean = SO_phase_KDE_bias_difference(npower).bias_diff_shifted_mean;
bias_CI_lo = SO_phase_KDE_bias_difference(npower).bias_diff_shifted_CI(1,:);
bias_CI_hi = SO_phase_KDE_bias_difference(npower).bias_diff_shifted_CI(2,:)

hold on;
x2 = [thresholds, fliplr(thresholds)];
y2 = [bias_CI_lo, fliplr(bias_CI_hi)];
Fill(end +1) = fill(x2, y2, 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
plot(thresholds, bias_mean, 'Color','k', 'LineWidth', 2);

yline(0, '--r');
legend(Fill([1 2 3]), {'SO peak', 'SO trough','Shuffled'}, 'box', 'off');

nexttile;
for npower = [1 2]
    bias_mean = SO_phase_KDE_bias_difference(npower).bias_diff_mean;
    bias_CI_lo = SO_phase_KDE_bias_difference(npower).bias_diff_CI(1,:);
    bias_CI_hi = SO_phase_KDE_bias_difference(npower).bias_diff_CI(2,:);
    prop_mean = SO_phase_KDE_bias_difference(npower).prop_mean;

    hold on;
    y2 = [prop_mean, fliplr(prop_mean)];
    x2 = [bias_CI_lo, fliplr(bias_CI_hi)];
    Fill(npower) = fill(x2, y2, colour_lines(npower,:), ...
        'EdgeColor', 'none', 'FaceAlpha', 0.3);
    plot(bias_mean, prop_mean, 'Color', colour_lines(npower,:), 'LineWidth', 2);

    xlabel('V1 bias diff (T1 - T2)');
    ylabel('Proportion of events detected');
    set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12);
    xlim([-0.1 0.35]);
    ylim([0 1])
end


bias_mean = SO_phase_KDE_bias_difference(npower).bias_diff_shifted_mean;
bias_CI_lo = SO_phase_KDE_bias_difference(npower).bias_diff_shifted_CI(1,:);
bias_CI_hi = SO_phase_KDE_bias_difference(npower).bias_diff_shifted_CI(2,:)
prop_mean = SO_phase_KDE_bias_difference(npower).prop_shifted_mean;

hold on;
y2 = [prop_mean, fliplr(prop_mean)];
x2 = [bias_CI_lo, fliplr(bias_CI_hi)];
Fill(end + 1) = fill(x2, y2, 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
plot(bias_mean, prop_mean, 'Color','k', 'LineWidth', 2);


xline(0, '--r');
legend(Fill([1 2 3]), {'SO peak', 'SO trough','Shuffled'}, 'box', 'off');

% Save results
save('SO_phase_KDE_bias_difference_based_on_V1_bias.mat', 'SO_phase_KDE_bias_difference');

load('SO_phase_KDE_bias_difference_based_on_V1_bias.mat', 'SO_phase_KDE_bias_difference');



%%%%% AUC mean + CI bar plot
% Plot layout
fig = figure;
fig.Position = [640 100 281 325]
fig.Name = 'KDE bias V1 AUC SO peak vs trough';
data = SO_phase_KDE_bias_difference;
n_bins = length(data);
bar_width = 0.3;      % Width of the bars
group_offset = 0.15;    % Distance from the center integer (half the gap between bars)
hold on;
clear BAR
for i = 1:4
    if i > length(data)
        x_shuf = i;
        y_shuf = data(1).AUC_mean_shuffled;
        bar(x_shuf, y_shuf, bar_width, ...
            'FaceColor', 'k', ...
            'FaceAlpha', 0.15, ...
            'EdgeColor', 'none');
    else


        % --- 1. Plot Shuffled Data (Left Bar) ---
        x_shuf = i - group_offset;
        y_shuf = data(i).AUC_mean_shuffled;

        % Calculate Error Deltas (Errorbar requires length relative to mean, not absolute values)
        % CI is [lower, upper]
        neg_err_shuf = y_shuf - data(i).AUC_CI_shuffled(1);
        pos_err_shuf = data(i).AUC_CI_shuffled(2) - y_shuf;

        % Plot Bar
        if i == 1
            BAR(1) = bar(x_shuf, y_shuf, bar_width, ...
                'FaceColor', 'k', ...
                'FaceAlpha', 0.15, ...
                'EdgeColor', 'none');
        else
            bar(x_shuf, y_shuf, bar_width, ...
                'FaceColor', 'k', ...
                'FaceAlpha', 0.15, ...
                'EdgeColor', 'none');
        end

        % Plot Error Bar
        E = errorbar(x_shuf, y_shuf, neg_err_shuf, pos_err_shuf, ...
            'Color', 'k', 'LineWidth', 1.5, 'CapSize', 8, 'LineStyle', 'none');


        % --- 2. Plot Real Data (Right Bar) ---
        x_real = i + group_offset;
        y_real = data(i).AUC_mean;

        neg_err_real = y_real - data(i).AUC_CI(1);
        pos_err_real = data(i).AUC_CI(2) - y_real;

        % Get specific color for this power bin
        % Assuming colour_lines is size [n_bins x 3]
        this_color = colour_lines(i, :);

        % Plot Bar
        BAR(i+1) = bar(x_real, y_real, bar_width, ...
            'FaceColor', this_color, ...
            'FaceAlpha', 0.3, ...
            'EdgeColor', 'none');

        errorbar(x_real, y_real, neg_err_real, pos_err_real, ...
            'Color', 'k', 'LineWidth', 1.5, 'CapSize', 8, 'LineStyle', 'none');
    end
end

hold off;
% Set X-ticks to be centered on the groups
set(gca, 'XTick', 1:nBins);
xlim([0.5, 4 + 0.5]);

% Labels
ylabel('V1 bias AUC');
xlabel('Power Bins');
legend([BAR(1:end)],{'Shuffled','SO peak','SO trough'},'box','off')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE bias difference based on V1 bias'),[])





%%%%%%%%%%%%%%%%%%%%%%%% SO phase binning PRE
% all_spindle_power = mean(ripple_info.spindle_amplitude, 1);  % avg of probe 1 and 2
% SO_thresholds = prctile(ripple_info.SO_phase, 0:99.9/4:99.9);
SO_thresholds = {'peak','trough'};
nBins = 2;

bins_to_use = bin_centers>0 & bin_centers<0.1;
bins_to_select = bin_centers>-0.2 & bin_centers<0;
bins_to_use_shifted = bin_centers > -1 & bin_centers < -0.9;
nBoot = 1000;

SO_phase_KDE_bias_difference = struct;

colour_lines = [ ...
    241, 182, 218;   % original end (lightest)
    231, 41, 138    % original start (darkest)
] / 256;

fig = figure;
fig.Position = [640 100 1100 650];
fig.Name = 'KDE bias difference in HPC with different SO phase (PRE ripple)';
tiledlayout(nBins, 3, 'TileSpacing', 'compact');

for npower = 1:nBins
    tic
    % Select events in current spindle bin
    % power_index = ripple_info.spindle_amplitude(:,1) > spindle_thresholds(npower,1) & ...
    %               ripple_info.spindle_amplitude(:,1) <= spindle_thresholds(npower+1,1) |...
    %               ripple_info.spindle_amplitude(:,2) > spindle_thresholds(npower,2) & ...
    %               ripple_info.spindle_amplitude(:,2) <= spindle_thresholds(npower+1,2);

    event_index =1:length(z_bias);
    % event_index = find(power_index);

    mean_bias = mean(z_bias_V1(bins_to_select, event_index), 'omitnan');
    % mean_bias_shifted = mean(z_bias(bins_to_use_shifted, event_index), 'omitnan');
    mean_bias_V1 = mean(z_bias(bins_to_use, event_index), 'omitnan');
    selected_events = length(mean_bias);

    thresholds = prctile(abs(mean_bias), 0:10:100);
    thresholds = thresholds(1:end-1);
    nThresh = length(thresholds);

    bias_diff_boot = NaN(nBoot, nThresh);
    prop_events_boot = NaN(nBoot, nThresh);
    bias_diff_shifted_boot = NaN(nBoot, nThresh);
    prop_events_shifted_boot = NaN(nBoot, nThresh);


    % event_phase = ripple_info.SO_phase';


    parfor iBoot = 1:nBoot
        s = RandStream('philox4x32_10', 'Seed', iBoot);
        idx = randi(s, selected_events, selected_events, 1);

        true_idx = find(event_index);

        bb = mean_bias(idx);
        bb_shift = mean_bias;
        boot_V1 = mean_bias_V1(idx);

        diff_tmp = NaN(1, nThresh);
        prop_tmp = NaN(1, nThresh);
        diff_tmp_shifted = NaN(1, nThresh);
        prop_tmp_shifted = NaN(1, nThresh);
        event_phase = ripple_info.SO_phase(true_idx(idx),:)';
        event_phase_shifted = ripple_info.SO_phase(true_idx,:)';

        for i = 1:nThresh
            th = thresholds(i);

            t1 = bb >= th;
            t2 = bb <= -th;

            if npower == 1 % if phase peak


                t1 = t1 & (event_phase(2,:) >= -pi/2 & event_phase(2,:) <= pi/2);
                t2 = t2 & (event_phase(1,:) >= -pi/2 & event_phase(1,:) <= pi/2);

                t1_V1 = boot_V1(t1);
                t2_V1 = boot_V1(t2);
                if any(t1) && any(t2)
                    diff_tmp(i) = mean(t1_V1) - mean(t2_V1);
                end

                total_events = mean([sum(event_phase(2,:) >= -pi/2 & event_phase(2,:) <= pi/2) ...
                    sum(event_phase(1,:) >= -pi/2 & event_phase(1,:) <= pi/2)]);

                prop_tmp(i) = (sum(t1) + sum(t2)) / total_events;



                t1s = bb_shift >= th;
                t2s = bb_shift <= -th;
                t1s = t1s & (event_phase_shifted(2,:) >= -pi/2 & event_phase_shifted(2,:) <= pi/2);
                t2s = t2s & (event_phase_shifted(1,:) >= -pi/2 & event_phase_shifted(1,:) <= pi/2);

                t1_V1 = boot_V1(t1s);
                t2_V1 = boot_V1(t2s);

                if any(t1s) && any(t2s)
                    diff_tmp_shifted(i) = mean(t1_V1) - mean(t2_V1);
                end
                prop_tmp_shifted(i) = (sum(t1s) + sum(t2s)) / total_events;

            elseif npower == 2 % if phase trough

                t1 = t1 & (event_phase(2,:) >= -pi & event_phase(2,:) <= -pi/2 | event_phase(2,:) >= pi/2 & event_phase(2,:) <= pi);
                t2 = t2 & (event_phase(1,:) >= -pi & event_phase(1,:) <= -pi/2 | event_phase(1,:) >= pi/2 & event_phase(1,:) <= pi);

                t1_V1 = boot_V1(t1);
                t2_V1 = boot_V1(t2);
                if any(t1) && any(t2)
                    diff_tmp(i) = mean(t1_V1) - mean(t2_V1);
                end

                total_events = mean([sum(event_phase(2,:) >= -pi & event_phase(2,:) <= -pi/2 | event_phase(2,:) >= pi/2 & event_phase(2,:) <= pi) ...
                    sum(event_phase(1,:) >= -pi & event_phase(1,:) <= -pi/2 | event_phase(1,:) >= pi/2 & event_phase(1,:) <= pi)]);

                prop_tmp(i) = (sum(t1) + sum(t2)) / total_events;



                t1s = bb_shift >= th;
                t2s = bb_shift <= -th;
                t1s = t1s & (event_phase_shifted(2,:) >= -pi & event_phase_shifted(2,:) <= -pi/2 | event_phase_shifted(2,:) >= pi/2 & event_phase_shifted(2,:) <= pi);
                t2s = t2s & (event_phase_shifted(1,:) >= -pi & event_phase_shifted(1,:) <= -pi/2 | event_phase_shifted(1,:) >= pi/2 & event_phase_shifted(1,:) <= pi);

                t1_V1 = boot_V1(t1s);
                t2_V1 = boot_V1(t2s);

                if any(t1s) && any(t2s)
                    diff_tmp_shifted(i) = mean(t1_V1) - mean(t2_V1);
                end
                prop_tmp_shifted(i) = (sum(t1s) + sum(t2s)) / total_events;

            end



        end

        bias_diff_boot(iBoot, :) = diff_tmp;
        prop_events_boot(iBoot, :) = prop_tmp;
        bias_diff_shifted_boot(iBoot, :) = diff_tmp_shifted;
        prop_events_shifted_boot(iBoot, :) = prop_tmp_shifted;
    end

    % Stats
    bias_mean = mean(bias_diff_boot, 1, 'omitnan');
    bias_CI_lo = prctile(bias_diff_boot, 2.5, 1);
    bias_CI_hi = prctile(bias_diff_boot, 97.5, 1);
    prop_mean = mean(prop_events_boot, 1, 'omitnan');
    prop_CI_lo = prctile(prop_events_boot, 2.5, 1);
    prop_CI_hi = prctile(prop_events_boot, 97.5, 1);
    bias_shifted_mean = mean(bias_diff_shifted_boot, 1, 'omitnan');
    bias_shifted_CI_lo = prctile(bias_diff_shifted_boot, 2.5, 1);
    bias_shifted_CI_hi = prctile(bias_diff_shifted_boot, 97.5, 1);
    prop_shifted_mean = mean(prop_events_shifted_boot, 1, 'omitnan');
    prop_shifted_CI_lo = prctile(prop_events_shifted_boot, 2.5, 1);
    prop_shifted_CI_hi = prctile(prop_events_shifted_boot, 97.5, 1);

    % Store
    SO_phase_KDE_bias_difference(npower).phase_range{npower} = ...
        SO_thresholds{npower};
    SO_phase_KDE_bias_difference(npower).bias_diff_mean = bias_mean;
    SO_phase_KDE_bias_difference(npower).bias_diff_CI = [bias_CI_lo; bias_CI_hi];
    SO_phase_KDE_bias_difference(npower).prop_mean = prop_mean;
    SO_phase_KDE_bias_difference(npower).prop_CI = [prop_CI_lo; prop_CI_hi];
    SO_phase_KDE_bias_difference(npower).thresholds = thresholds;
    SO_phase_KDE_bias_difference(npower).bias_diff_shifted_mean = bias_shifted_mean;
    SO_phase_KDE_bias_difference(npower).bias_diff_shifted_CI = [bias_shifted_CI_lo; bias_shifted_CI_hi];
    SO_phase_KDE_bias_difference(npower).prop_shifted_mean = prop_shifted_mean;
    SO_phase_KDE_bias_difference(npower).prop_shifted_CI = [prop_shifted_CI_lo; prop_shifted_CI_hi];


    % store AUC
    auc_boot = (trapz(thresholds, bias_diff_boot') / (max(thresholds)-min(thresholds)))';
    auc_shift_boot = (trapz(thresholds, bias_diff_shifted_boot') / (max(thresholds)-min(thresholds)))';

    SO_phase_KDE_bias_difference(npower).AUC_mean = mean(auc_boot, 'omitnan');
    SO_phase_KDE_bias_difference(npower).AUC_CI = prctile(auc_boot, [2.5 97.5]);
    SO_phase_KDE_bias_difference(npower).AUC_mean_shuffled = mean(auc_shift_boot, 'omitnan');
    SO_phase_KDE_bias_difference(npower).AUC_CI_shuffled = prctile(auc_shift_boot, [2.5 97.5]);


    % ----------- PLOTS (A, B, C) ----------------

    % A: Bias vs. Threshold
    nexttile((npower-1)*3 + 1);
    hold on;
    fill([thresholds, fliplr(thresholds)], [bias_CI_lo, fliplr(bias_CI_hi)], ...
        colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    plot(thresholds, bias_mean, 'Color', colour_lines(npower,:), 'LineWidth', 2);
    fill([thresholds, fliplr(thresholds)], ...
         [bias_shifted_CI_lo, fliplr(bias_shifted_CI_hi)], ...
         [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    plot(thresholds, bias_shifted_mean, 'k-', 'LineWidth', 1.5);
    ylim([-0.15 0.35]); xlim([0 1]); yline(0, '--r');
    xlabel('V1 bias threshold'); ylabel('HPC bias diff (T1 - T2)');
    title(sprintf('SO phase bin %s', SO_thresholds{npower}));
    set(gca,"TickDir","out",'box','off','Color','none','FontSize',12)

    % B: Proportion vs Bias Diff
    nexttile((npower-1)*3 + 2);
    hold on;
    valid = isfinite(bias_mean) & isfinite(prop_mean);
    fill([bias_CI_lo(valid), fliplr(bias_CI_hi(valid))], ...
         [prop_mean(valid), fliplr(prop_mean(valid))], ...
         colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    plot(bias_mean(valid), prop_mean(valid), '-', 'Color', colour_lines(npower,:), 'LineWidth', 2);
    fill([bias_shifted_CI_lo(valid), fliplr(bias_shifted_CI_hi(valid))], ...
         [prop_shifted_mean(valid), fliplr(prop_shifted_mean(valid))], ...
         [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    plot(bias_shifted_mean(valid), prop_shifted_mean(valid), 'k-', 'LineWidth', 1.5);
    xlim([-0.1 0.35]); xline(0, '--r');
    xlabel('HPC bias diff (T1 - T2)'); ylabel('Proportion of events detected');
    title('Event Proportion vs. Bias Difference');
    set(gca,"TickDir","out",'box','off','Color','none','FontSize',12)

    % C: Proportion CI vs Bias Diff
    nexttile((npower-1)*3 + 3);
    hold on;
    fill([bias_mean(valid), fliplr(bias_mean(valid))], ...
         [prop_CI_lo(valid), fliplr(prop_CI_hi(valid))], ...
         colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    plot(bias_mean(valid), prop_mean(valid), '-', 'Color', colour_lines(npower,:), 'LineWidth', 2);
    fill([bias_shifted_mean(valid), fliplr(bias_shifted_mean(valid))], ...
         [prop_shifted_CI_lo(valid), fliplr(prop_shifted_CI_hi(valid))], ...
         [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    plot(bias_shifted_mean(valid), prop_shifted_mean(valid), 'k-', 'LineWidth', 1.5);
    xlim([-0.1 0.35]); xline(0, '--r');
    xlabel('HPC bias diff (T1 - T2)'); ylabel('Proportion of events detected');
    title('Proportion vs. HPC Bias Difference');
    set(gca,"TickDir","out",'box','off','Color','none','FontSize',12)
    toc
end

clear Fill

fig = figure;
fig.Position = [640 100 2*1100/3 650/2];
% fig.Name = 'KDE bias difference in HPC low vs high SO power (PRE ripple)';
fig.Name = 'KDE bias difference in HPC SO phase (PRE ripple)';

% Select bins (1 = low, 4 = high)
nexttile;
for npower = [1 2]
    bias_mean = SO_phase_KDE_bias_difference(npower).bias_diff_mean;
    bias_CI_lo = SO_phase_KDE_bias_difference(npower).bias_diff_CI(1,:);
    bias_CI_hi = SO_phase_KDE_bias_difference(npower).bias_diff_CI(2,:);
    thresholds = SO_phase_KDE_bias_difference(npower).thresholds;

    hold on;
    x2 = [thresholds, fliplr(thresholds)];
    y2 = [bias_CI_lo, fliplr(bias_CI_hi)];
    Fill(npower) = fill(x2, y2, colour_lines(npower,:), ...
        'EdgeColor', 'none', 'FaceAlpha', 0.3);
    plot(thresholds, bias_mean, 'Color', colour_lines(npower,:), 'LineWidth', 2);

    xlabel('HPC Bias threshold');
    ylabel('V1 bias diff (T1 - T2)');
    set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12);
    ylim([-0.1 0.35]);
    % ylim([0 1])
end

bias_mean = SO_phase_KDE_bias_difference(npower).bias_diff_shifted_mean;
bias_CI_lo = SO_phase_KDE_bias_difference(npower).bias_diff_shifted_CI(1,:);
bias_CI_hi = SO_phase_KDE_bias_difference(npower).bias_diff_shifted_CI(2,:)

hold on;
x2 = [thresholds, fliplr(thresholds)];
y2 = [bias_CI_lo, fliplr(bias_CI_hi)];
Fill(end +1) = fill(x2, y2, 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
plot(thresholds, bias_mean, 'Color','k', 'LineWidth', 2);

yline(0, '--r');
legend(Fill([1 2 3]), {'SO peak', 'SO trough','Shuffled'}, 'box', 'off');

nexttile;
for npower = [1 2]
    bias_mean = SO_phase_KDE_bias_difference(npower).bias_diff_mean;
    bias_CI_lo = SO_phase_KDE_bias_difference(npower).bias_diff_CI(1,:);
    bias_CI_hi = SO_phase_KDE_bias_difference(npower).bias_diff_CI(2,:);
    prop_mean = SO_phase_KDE_bias_difference(npower).prop_mean;

    hold on;
    y2 = [prop_mean, fliplr(prop_mean)];
    x2 = [bias_CI_lo, fliplr(bias_CI_hi)];
    Fill(npower) = fill(x2, y2, colour_lines(npower,:), ...
        'EdgeColor', 'none', 'FaceAlpha', 0.3);
    plot(bias_mean, prop_mean, 'Color', colour_lines(npower,:), 'LineWidth', 2);

    xlabel('V1 bias diff (T1 - T2)');
    ylabel('Proportion of events detected');
    set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12);
    xlim([-0.1 0.35]);
    ylim([0 1])
end


bias_mean = SO_phase_KDE_bias_difference(npower).bias_diff_shifted_mean;
bias_CI_lo = SO_phase_KDE_bias_difference(npower).bias_diff_shifted_CI(1,:);
bias_CI_hi = SO_phase_KDE_bias_difference(npower).bias_diff_shifted_CI(2,:)
prop_mean = SO_phase_KDE_bias_difference(npower).prop_shifted_mean;

hold on;
y2 = [prop_mean, fliplr(prop_mean)];
x2 = [bias_CI_lo, fliplr(bias_CI_hi)];
Fill(end + 1) = fill(x2, y2, 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
plot(bias_mean, prop_mean, 'Color','k', 'LineWidth', 2);


xline(0, '--r');
legend(Fill([1 2 3]), {'SO peak', 'SO trough','Shuffled'}, 'box', 'off');

% Save results
save('SO_phase_KDE_bias_difference_based_on_PRE_V1_bias.mat', 'SO_phase_KDE_bias_difference');




%%%%% AUC mean + CI bar plot
% Plot layout
fig = figure;
fig.Position = [640 100 281 325]
fig.Name = 'KDE bias V1 AUC SO peak vs trough (PRE)';
data = SO_phase_KDE_bias_difference;
n_bins = length(data);
bar_width = 0.3;      % Width of the bars
group_offset = 0.15;    % Distance from the center integer (half the gap between bars)
hold on;
clear BAR
for i = 1:4
    if i > length(data)
        x_shuf = i;
        y_shuf = data(1).AUC_mean_shuffled;
        bar(x_shuf, y_shuf, bar_width, ...
            'FaceColor', 'k', ...
            'FaceAlpha', 0.15, ...
            'EdgeColor', 'none');
    else


        % --- 1. Plot Shuffled Data (Left Bar) ---
        x_shuf = i - group_offset;
        y_shuf = data(i).AUC_mean_shuffled;

        % Calculate Error Deltas (Errorbar requires length relative to mean, not absolute values)
        % CI is [lower, upper]
        neg_err_shuf = y_shuf - data(i).AUC_CI_shuffled(1);
        pos_err_shuf = data(i).AUC_CI_shuffled(2) - y_shuf;

        % Plot Bar
        if i == 1
            BAR(1) = bar(x_shuf, y_shuf, bar_width, ...
                'FaceColor', 'k', ...
                'FaceAlpha', 0.15, ...
                'EdgeColor', 'none');
        else
            bar(x_shuf, y_shuf, bar_width, ...
                'FaceColor', 'k', ...
                'FaceAlpha', 0.15, ...
                'EdgeColor', 'none');
        end

        % Plot Error Bar
        E = errorbar(x_shuf, y_shuf, neg_err_shuf, pos_err_shuf, ...
            'Color', 'k', 'LineWidth', 1.5, 'CapSize', 8, 'LineStyle', 'none');


        % --- 2. Plot Real Data (Right Bar) ---
        x_real = i + group_offset;
        y_real = data(i).AUC_mean;

        neg_err_real = y_real - data(i).AUC_CI(1);
        pos_err_real = data(i).AUC_CI(2) - y_real;

        % Get specific color for this power bin
        % Assuming colour_lines is size [n_bins x 3]
        this_color = colour_lines(i, :);

        % Plot Bar
        BAR(i+1) = bar(x_real, y_real, bar_width, ...
            'FaceColor', this_color, ...
            'FaceAlpha', 0.3, ...
            'EdgeColor', 'none');

        errorbar(x_real, y_real, neg_err_real, pos_err_real, ...
            'Color', 'k', 'LineWidth', 1.5, 'CapSize', 8, 'LineStyle', 'none');
    end
end

hold off;
% Set X-ticks to be centered on the groups
set(gca, 'XTick', 1:nBins);
xlim([0.5, 4 + 0.5]);

% Labels
ylabel('V1 bias AUC');
xlabel('Power Bins');
legend([BAR(1:end)],{'Shuffled','SO peak','SO trough'},'box','off')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE bias difference based on V1 bias'),[])



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% SO phases
%%%%%%%%%%%%%%%%
% SO phase binning across both probes
% all_spindle_power = mean(ripple_info.spindle_amplitude, 1);  % avg of probe 1 and 2
% SO_thresholds = prctile(ripple_info.SO_phase, 0:99.9/4:99.9);
SO_thresholds = {'peak','trough'};
nBins = 2;

bins_to_use = bin_centers>0 & bin_centers<0.1;
bins_to_select = bin_centers>-0.05 & bin_centers<0.05;
bins_to_use_shifted = bin_centers > -1 & bin_centers < -0.9;
nBoot = 1000;

SO_power_KDE_bias_difference = struct;

colour_lines = [ ...
    241, 182, 218;   % original end (lightest)
    231, 41, 138    % original start (darkest)
] / 256;

fig = figure;
fig.Position = [640 100 1100 650];
fig.Name = 'KDE bias difference in HPC with different SO phase';
tiledlayout(nBins, 3, 'TileSpacing', 'compact');

for npower = 1:nBins
    tic
    % Select events in current spindle bin
    % power_index = ripple_info.spindle_amplitude(:,1) > spindle_thresholds(npower,1) & ...
    %               ripple_info.spindle_amplitude(:,1) <= spindle_thresholds(npower+1,1) |...
    %               ripple_info.spindle_amplitude(:,2) > spindle_thresholds(npower,2) & ...
    %               ripple_info.spindle_amplitude(:,2) <= spindle_thresholds(npower+1,2);

    event_index =1:length(z_bias);
    % event_index = find(power_index);

    mean_bias = mean(z_bias_V1(bins_to_select, event_index), 'omitnan');
    % mean_bias_shifted = mean(z_bias(bins_to_use_shifted, event_index), 'omitnan');
    mean_bias_V1 = mean(z_bias(bins_to_use, event_index), 'omitnan');
    selected_events = length(mean_bias);

    thresholds = prctile(abs(mean_bias), 0:10:100);
    thresholds = thresholds(1:end-1);
    nThresh = length(thresholds);

    bias_diff_boot = NaN(nBoot, nThresh);
    prop_events_boot = NaN(nBoot, nThresh);
    bias_diff_shifted_boot = NaN(nBoot, nThresh);
    prop_events_shifted_boot = NaN(nBoot, nThresh);


    parfor iBoot = 1:nBoot
        s = RandStream('philox4x32_10', 'Seed', iBoot);
        idx = randi(s, selected_events, selected_events, 1);

        true_idx = find(event_index);

        bb = mean_bias(idx);
        bb_shift = mean_bias;
        boot_V1 = mean_bias_V1(idx);

        diff_tmp = NaN(1, nThresh);
        prop_tmp = NaN(1, nThresh);
        diff_tmp_shifted = NaN(1, nThresh);
        prop_tmp_shifted = NaN(1, nThresh);
        event_phase = ripple_info.SO_phase(true_idx(idx),:)';
        event_phase_shifted = ripple_info.SO_phase(true_idx,:)';

        for i = 1:nThresh
            th = thresholds(i);

            t1 = bb >= th;
            t2 = bb <= -th;

            if npower == 1 % if phase peak


                t1 = t1 & (event_phase(2,:) >= -pi/2 & event_phase(2,:) <= pi/2);
                t2 = t2 & (event_phase(1,:) >= -pi/2 & event_phase(1,:) <= pi/2);

                t1_V1 = boot_V1(t1);
                t2_V1 = boot_V1(t2);
                if any(t1) && any(t2)
                    diff_tmp(i) = mean(t1_V1) - mean(t2_V1);
                end

                total_events = mean([sum(event_phase(2,:) >= -pi/2 & event_phase(2,:) <= pi/2) ...
                    sum(event_phase(1,:) >= -pi/2 & event_phase(1,:) <= pi/2)]);

                prop_tmp(i) = (sum(t1) + sum(t2)) / total_events;



                t1s = bb_shift >= th;
                t2s = bb_shift <= -th;
                t1s = t1s & (event_phase_shifted(2,:) >= -pi/2 & event_phase_shifted(2,:) <= pi/2);
                t2s = t2s & (event_phase_shifted(1,:) >= -pi/2 & event_phase_shifted(1,:) <= pi/2);

                t1_V1 = boot_V1(t1s);
                t2_V1 = boot_V1(t2s);

                if any(t1s) && any(t2s)
                    diff_tmp_shifted(i) = mean(t1_V1) - mean(t2_V1);
                end
                prop_tmp_shifted(i) = (sum(t1s) + sum(t2s)) / total_events;

            elseif npower == 2 % if phase trough

                t1 = t1 & (event_phase(2,:) >= -pi & event_phase(2,:) <= -pi/2 | event_phase(2,:) >= pi/2 & event_phase(2,:) <= pi);
                t2 = t2 & (event_phase(1,:) >= -pi & event_phase(1,:) <= -pi/2 | event_phase(1,:) >= pi/2 & event_phase(1,:) <= pi);

                t1_V1 = boot_V1(t1);
                t2_V1 = boot_V1(t2);
                if any(t1) && any(t2)
                    diff_tmp(i) = mean(t1_V1) - mean(t2_V1);
                end

                total_events = mean([sum(event_phase(2,:) >= -pi & event_phase(2,:) <= -pi/2 | event_phase(2,:) >= pi/2 & event_phase(2,:) <= pi) ...
                    sum(event_phase(1,:) >= -pi & event_phase(1,:) <= -pi/2 | event_phase(1,:) >= pi/2 & event_phase(1,:) <= pi)]);

                prop_tmp(i) = (sum(t1) + sum(t2)) / total_events;



                t1s = bb_shift >= th;
                t2s = bb_shift <= -th;
                t1s = t1s & (event_phase_shifted(2,:) >= -pi & event_phase_shifted(2,:) <= -pi/2 | event_phase_shifted(2,:) >= pi/2 & event_phase_shifted(2,:) <= pi);
                t2s = t2s & (event_phase_shifted(1,:) >= -pi & event_phase_shifted(1,:) <= -pi/2 | event_phase_shifted(1,:) >= pi/2 & event_phase_shifted(1,:) <= pi);

                t1_V1 = boot_V1(t1s);
                t2_V1 = boot_V1(t2s);

                if any(t1s) && any(t2s)
                    diff_tmp_shifted(i) = mean(t1_V1) - mean(t2_V1);
                end
                prop_tmp_shifted(i) = (sum(t1s) + sum(t2s)) / total_events;

            end



        end

        bias_diff_boot(iBoot, :) = diff_tmp;
        prop_events_boot(iBoot, :) = prop_tmp;
        bias_diff_shifted_boot(iBoot, :) = diff_tmp_shifted;
        prop_events_shifted_boot(iBoot, :) = prop_tmp_shifted;
    end

    % Stats
    bias_mean = mean(bias_diff_boot, 1, 'omitnan');
    bias_CI_lo = prctile(bias_diff_boot, 2.5, 1);
    bias_CI_hi = prctile(bias_diff_boot, 97.5, 1);
    prop_mean = mean(prop_events_boot, 1, 'omitnan');
    prop_CI_lo = prctile(prop_events_boot, 2.5, 1);
    prop_CI_hi = prctile(prop_events_boot, 97.5, 1);
    bias_shifted_mean = mean(bias_diff_shifted_boot, 1, 'omitnan');
    bias_shifted_CI_lo = prctile(bias_diff_shifted_boot, 2.5, 1);
    bias_shifted_CI_hi = prctile(bias_diff_shifted_boot, 97.5, 1);
    prop_shifted_mean = mean(prop_events_shifted_boot, 1, 'omitnan');
    prop_shifted_CI_lo = prctile(prop_events_shifted_boot, 2.5, 1);
    prop_shifted_CI_hi = prctile(prop_events_shifted_boot, 97.5, 1);

    % Store
    SO_phase_KDE_bias_difference(npower).phase_range{npower} = ...
        SO_thresholds{npower};
    SO_phase_KDE_bias_difference(npower).bias_diff_mean = bias_mean;
    SO_phase_KDE_bias_difference(npower).bias_diff_CI = [bias_CI_lo; bias_CI_hi];
    SO_phase_KDE_bias_difference(npower).prop_mean = prop_mean;
    SO_phase_KDE_bias_difference(npower).prop_CI = [prop_CI_lo; prop_CI_hi];
    SO_phase_KDE_bias_difference(npower).thresholds = thresholds;
    SO_phase_KDE_bias_difference(npower).bias_diff_shifted_mean = bias_shifted_mean;
    SO_phase_KDE_bias_difference(npower).bias_diff_shifted_CI = [bias_shifted_CI_lo; bias_shifted_CI_hi];
    SO_phase_KDE_bias_difference(npower).prop_shifted_mean = prop_shifted_mean;
    SO_phase_KDE_bias_difference(npower).prop_shifted_CI = [prop_shifted_CI_lo; prop_shifted_CI_hi];

    % store AUC
    auc_boot = (trapz(thresholds, bias_diff_boot') / (max(thresholds)-min(thresholds)))';
    auc_shift_boot = (trapz(thresholds, bias_diff_shifted_boot') / (max(thresholds)-min(thresholds)))';

    SO_phase_KDE_bias_difference(npower).AUC_mean = mean(auc_boot, 'omitnan');
    SO_phase_KDE_bias_difference(npower).AUC_CI = prctile(auc_boot, [2.5 97.5]);
    SO_phase_KDE_bias_difference(npower).AUC_mean_shuffled = mean(auc_shift_boot, 'omitnan');
    SO_phase_KDE_bias_difference(npower).AUC_CI_shuffled = prctile(auc_shift_boot, [2.5 97.5]);


    % ----------- PLOTS (A, B, C) ----------------

    % A: Bias vs. Threshold
    nexttile((npower-1)*3 + 1);
    hold on;
    fill([thresholds, fliplr(thresholds)], [bias_CI_lo, fliplr(bias_CI_hi)], ...
        colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    plot(thresholds, bias_mean, 'Color', colour_lines(npower,:), 'LineWidth', 2);
    fill([thresholds, fliplr(thresholds)], ...
         [bias_shifted_CI_lo, fliplr(bias_shifted_CI_hi)], ...
         [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    plot(thresholds, bias_shifted_mean, 'k-', 'LineWidth', 1.5);
    ylim([-0.15 0.35]); xlim([0 1]); yline(0, '--r');
    xlabel('V1 bias threshold'); ylabel('HPC bias diff (T1 - T2)');
    title(sprintf('SO phase bin %s', SO_thresholds{npower}));
    set(gca,"TickDir","out",'box','off','Color','none','FontSize',12)

    % B: Proportion vs Bias Diff
    nexttile((npower-1)*3 + 2);
    hold on;
    valid = isfinite(bias_mean) & isfinite(prop_mean);
    fill([bias_CI_lo(valid), fliplr(bias_CI_hi(valid))], ...
         [prop_mean(valid), fliplr(prop_mean(valid))], ...
         colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    plot(bias_mean(valid), prop_mean(valid), '-', 'Color', colour_lines(npower,:), 'LineWidth', 2);
    fill([bias_shifted_CI_lo(valid), fliplr(bias_shifted_CI_hi(valid))], ...
         [prop_shifted_mean(valid), fliplr(prop_shifted_mean(valid))], ...
         [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    plot(bias_shifted_mean(valid), prop_shifted_mean(valid), 'k-', 'LineWidth', 1.5);
    xlim([-0.1 0.35]); xline(0, '--r');
    xlabel('HPC bias diff (T1 - T2)'); ylabel('Proportion of events detected');
    title('Event Proportion vs. Bias Difference');
    set(gca,"TickDir","out",'box','off','Color','none','FontSize',12)

    % C: Proportion CI vs Bias Diff
    nexttile((npower-1)*3 + 3);
    hold on;
    fill([bias_mean(valid), fliplr(bias_mean(valid))], ...
         [prop_CI_lo(valid), fliplr(prop_CI_hi(valid))], ...
         colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    plot(bias_mean(valid), prop_mean(valid), '-', 'Color', colour_lines(npower,:), 'LineWidth', 2);
    fill([bias_shifted_mean(valid), fliplr(bias_shifted_mean(valid))], ...
         [prop_shifted_CI_lo(valid), fliplr(prop_shifted_CI_hi(valid))], ...
         [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    plot(bias_shifted_mean(valid), prop_shifted_mean(valid), 'k-', 'LineWidth', 1.5);
    xlim([-0.1 0.35]); xline(0, '--r');
    xlabel('HPC bias diff (T1 - T2)'); ylabel('Proportion of events detected');
    title('Proportion vs. HPC Bias Difference');
    set(gca,"TickDir","out",'box','off','Color','none','FontSize',12)
    toc
end

clear Fill

fig = figure;
fig.Position = [640 100 2*1100/3 650/2];
fig.Name = 'KDE bias difference in HPC SO phase (100ms example)';

% Select bins (1 = low, 4 = high)
nexttile;
for npower = [1 2]
    bias_mean = SO_phase_KDE_bias_difference(npower).bias_diff_mean;
    bias_CI_lo = SO_phase_KDE_bias_difference(npower).bias_diff_CI(1,:);
    bias_CI_hi = SO_phase_KDE_bias_difference(npower).bias_diff_CI(2,:);
    thresholds = SO_phase_KDE_bias_difference(npower).thresholds;

    hold on;
    x2 = [thresholds, fliplr(thresholds)];
    y2 = [bias_CI_lo, fliplr(bias_CI_hi)];
    Fill(npower) = fill(x2, y2, colour_lines(npower,:), ...
        'EdgeColor', 'none', 'FaceAlpha', 0.3);
    plot(thresholds, bias_mean, 'Color', colour_lines(npower,:), 'LineWidth', 2);

    xlabel('HPC Bias threshold');
    ylabel('V1 bias diff (T1 - T2)');
    set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12);
    ylim([-0.1 0.35]);
    % ylim([0 1])
end

bias_mean = SO_phase_KDE_bias_difference(npower).bias_diff_shifted_mean;
bias_CI_lo = SO_phase_KDE_bias_difference(npower).bias_diff_shifted_CI(1,:);
bias_CI_hi = SO_phase_KDE_bias_difference(npower).bias_diff_shifted_CI(2,:)

hold on;
x2 = [thresholds, fliplr(thresholds)];
y2 = [bias_CI_lo, fliplr(bias_CI_hi)];
Fill(end +1) = fill(x2, y2, 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
plot(thresholds, bias_mean, 'Color','k', 'LineWidth', 2);

yline(0, '--r');
legend(Fill([1 2 3]), {'SO peak', 'SO trough','Shuffled'}, 'box', 'off');

nexttile;
for npower = [1 2]
    bias_mean = SO_phase_KDE_bias_difference(npower).bias_diff_mean;
    bias_CI_lo = SO_phase_KDE_bias_difference(npower).bias_diff_CI(1,:);
    bias_CI_hi = SO_phase_KDE_bias_difference(npower).bias_diff_CI(2,:);
    prop_mean = SO_phase_KDE_bias_difference(npower).prop_mean;

    hold on;
    y2 = [prop_mean, fliplr(prop_mean)];
    x2 = [bias_CI_lo, fliplr(bias_CI_hi)];
    Fill(npower) = fill(x2, y2, colour_lines(npower,:), ...
        'EdgeColor', 'none', 'FaceAlpha', 0.3);
    plot(bias_mean, prop_mean, 'Color', colour_lines(npower,:), 'LineWidth', 2);

    xlabel('V1 bias diff (T1 - T2)');
    ylabel('Proportion of events detected');
    set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12);
    xlim([-0.1 0.35]);
    ylim([0 1])
end


bias_mean = SO_phase_KDE_bias_difference(npower).bias_diff_shifted_mean;
bias_CI_lo = SO_phase_KDE_bias_difference(npower).bias_diff_shifted_CI(1,:);
bias_CI_hi = SO_phase_KDE_bias_difference(npower).bias_diff_shifted_CI(2,:)
prop_mean = SO_phase_KDE_bias_difference(npower).prop_shifted_mean;

hold on;
y2 = [prop_mean, fliplr(prop_mean)];
x2 = [bias_CI_lo, fliplr(bias_CI_hi)];
Fill(end + 1) = fill(x2, y2, 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
plot(bias_mean, prop_mean, 'Color','k', 'LineWidth', 2);


xline(0, '--r');
legend(Fill([1 2 3]), {'SO peak', 'SO trough','Shuffled'}, 'box', 'off');


save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE bias difference based on V1 bias'),[])






%% SO peak synchrony
% Spindle power binning across both probes
% all_spindle_power = mean(ripple_info.spindle_amplitude, 1);  % avg of probe 1 and 2
% SO_thresholds = prctile(ripple_info.SO_phase, 0:99.9/4:99.9);
SO_thresholds = {'peak not sync','peak sync'};
nBins = 2;

bins_to_use = bin_centers>0 & bin_centers<0.1;
bins_to_select = bin_centers>=0 & bin_centers<0.2;
bins_to_use_shifted = bin_centers > -1 & bin_centers < -0.9;
nBoot = 1000;

SO_power_KDE_bias_difference = struct;

colour_lines = [ ...
    241, 182, 218;   % original end (lightest)
    231, 41, 138    % original start (darkest)
] / 256;

fig = figure;
fig.Position = [640 100 1100 650];
fig.Name = 'KDE bias difference in HPC with different SO peak synchrony';
tiledlayout(nBins, 3, 'TileSpacing', 'compact');

for npower = 1:nBins
    tic
    % Select events in current spindle bin
    % power_index = ripple_info.spindle_amplitude(:,1) > spindle_thresholds(npower,1) & ...
    %               ripple_info.spindle_amplitude(:,1) <= spindle_thresholds(npower+1,1) |...
    %               ripple_info.spindle_amplitude(:,2) > spindle_thresholds(npower,2) & ...
    %               ripple_info.spindle_amplitude(:,2) <= spindle_thresholds(npower+1,2);

    event_index =1:length(z_bias);
    % event_index = find(power_index);

    mean_bias = mean(z_bias_V1(bins_to_select, event_index), 'omitnan');
    % mean_bias_shifted = mean(z_bias(bins_to_use_shifted, event_index), 'omitnan');
    mean_bias_V1 = mean(z_bias(bins_to_use, event_index), 'omitnan');
    selected_events = length(mean_bias);

    thresholds = prctile(abs(mean_bias), 0:10:100);
    thresholds = thresholds(1:end-1);
    nThresh = length(thresholds);

    bias_diff_boot = NaN(nBoot, nThresh);
    prop_events_boot = NaN(nBoot, nThresh);
    bias_diff_shifted_boot = NaN(nBoot, nThresh);
    prop_events_shifted_boot = NaN(nBoot, nThresh);

    % event_phase = ripple_info.SO_phase';

    parfor iBoot = 1:nBoot
        s = RandStream('philox4x32_10', 'Seed', iBoot);
        idx = randi(s, selected_events, selected_events, 1);

        true_idx = find(event_index);

        bb = mean_bias(idx);
        bb_shift = mean_bias;
        boot_V1 = mean_bias_V1(idx);

        diff_tmp = NaN(1, nThresh);
        prop_tmp = NaN(1, nThresh);
        diff_tmp_shifted = NaN(1, nThresh);
        prop_tmp_shifted = NaN(1, nThresh);
        event_phase = ripple_info.SO_phase(true_idx(idx),:)';
        event_phase_shifted = ripple_info.SO_phase(true_idx,:)';

        for i = 1:nThresh
            th = thresholds(i);

            t1 = bb >= th;
            t2 = bb <= -th;

            if npower == 1 % if peak not sync


                t1 = t1 & (event_phase(2,:) >= -pi/2 & event_phase(2,:) <= pi/2) & (event_phase(1,:) >= -pi & event_phase(1,:) < -pi/2 | event_phase(1,:) > pi/2 & event_phase(1,:) <= pi);
                t2 = t2 & (event_phase(1,:) >= -pi/2 & event_phase(1,:) <= pi/2) & (event_phase(2,:) >= -pi & event_phase(2,:) < -pi/2 | event_phase(2,:) > pi/2 & event_phase(2,:) <= pi);

                t1_V1 = boot_V1(t1);
                t2_V1 = boot_V1(t2);
                if any(t1) && any(t2)
                    diff_tmp(i) = mean(t1_V1) - mean(t2_V1);
                end

                total_events = mean([sum((event_phase(2,:) >= -pi/2 & event_phase(2,:) <= pi/2) & (event_phase(1,:) >= -pi & event_phase(1,:) <= -pi/2 | event_phase(1,:) >= pi/2 & event_phase(1,:) <= pi)) ...
                    sum((event_phase(1,:) >= -pi/2 & event_phase(1,:) <= pi/2) & (event_phase(2,:) >= -pi & event_phase(2,:) <= -pi/2 | event_phase(2,:) >= pi/2 & event_phase(2,:) <= pi))]);

                prop_tmp(i) = (sum(t1) + sum(t2)) / total_events;



                t1s = bb_shift >= th;
                t2s = bb_shift <= -th;
                t1s = t1s & (event_phase_shifted(2,:) >= -pi/2 & event_phase_shifted(2,:) <= pi/2) & (event_phase_shifted(1,:) >= -pi & event_phase_shifted(1,:) < -pi/2 | event_phase_shifted(1,:) > pi/2 & event_phase_shifted(1,:) <= pi);
                t2s = t2s & (event_phase_shifted(1,:) >= -pi/2 & event_phase_shifted(1,:) <= pi/2) & (event_phase_shifted(2,:) >= -pi & event_phase_shifted(2,:) < -pi/2 | event_phase_shifted(2,:) > pi/2 & event_phase_shifted(2,:) <= pi);

                t1_V1 = boot_V1(t1s);
                t2_V1 = boot_V1(t2s);

                if any(t1s) && any(t2s)
                    diff_tmp_shifted(i) = mean(t1_V1) - mean(t2_V1);
                end
                prop_tmp_shifted(i) = (sum(t1s) + sum(t2s)) / total_events;

            elseif npower == 2 % if peak sync

                t1 = t1 & (event_phase(2,:) >= -pi/2 & event_phase(2,:) <= pi/2) & (event_phase(1,:) >= -pi/2 & event_phase(1,:) <= pi/2);
                t2 = t2 & (event_phase(2,:) >= -pi/2 & event_phase(2,:) <= pi/2) & (event_phase(1,:) >= -pi/2 & event_phase(1,:) <= pi/2);

                t1_V1 = boot_V1(t1);
                t2_V1 = boot_V1(t2);
                if any(t1) && any(t2)
                    diff_tmp(i) = mean(t1_V1) - mean(t2_V1);
                end

                total_events = sum((event_phase(2,:) >= -pi/2 & event_phase(2,:) <= pi/2) & (event_phase(1,:) >= -pi/2 & event_phase(1,:) <= pi/2));

                prop_tmp(i) = (sum(t1) + sum(t2)) / total_events;



                t1s = bb_shift >= th;
                t2s = bb_shift <= -th;
                t1s = t1s & (event_phase_shifted(2,:) >= -pi/2 & event_phase_shifted(2,:) <= pi/2) & (event_phase_shifted(1,:) >= -pi/2 & event_phase_shifted(1,:) <= pi/2);
                t2s = t2s & (event_phase_shifted(2,:) >= -pi/2 & event_phase_shifted(2,:) <= pi/2) & (event_phase_shifted(1,:) >= -pi/2 & event_phase_shifted(1,:) <= pi/2);


                t1_V1 = boot_V1(t1s);
                t2_V1 = boot_V1(t2s);

                if any(t1s) && any(t2s)
                    diff_tmp_shifted(i) = mean(t1_V1) - mean(t2_V1);
                end
                prop_tmp_shifted(i) = (sum(t1s) + sum(t2s)) / total_events;

            end



        end

        bias_diff_boot(iBoot, :) = diff_tmp;
        prop_events_boot(iBoot, :) = prop_tmp;
        bias_diff_shifted_boot(iBoot, :) = diff_tmp_shifted;
        prop_events_shifted_boot(iBoot, :) = prop_tmp_shifted;
    end

    % Stats
    bias_mean = mean(bias_diff_boot, 1, 'omitnan');
    bias_CI_lo = prctile(bias_diff_boot, 2.5, 1);
    bias_CI_hi = prctile(bias_diff_boot, 97.5, 1);
    prop_mean = mean(prop_events_boot, 1, 'omitnan');
    prop_CI_lo = prctile(prop_events_boot, 2.5, 1);
    prop_CI_hi = prctile(prop_events_boot, 97.5, 1);
    bias_shifted_mean = mean(bias_diff_shifted_boot, 1, 'omitnan');
    bias_shifted_CI_lo = prctile(bias_diff_shifted_boot, 2.5, 1);
    bias_shifted_CI_hi = prctile(bias_diff_shifted_boot, 97.5, 1);
    prop_shifted_mean = mean(prop_events_shifted_boot, 1, 'omitnan');
    prop_shifted_CI_lo = prctile(prop_events_shifted_boot, 2.5, 1);
    prop_shifted_CI_hi = prctile(prop_events_shifted_boot, 97.5, 1);

    % Store
    SO_phase_KDE_bias_difference(npower).phase_range{npower} = ...
        SO_thresholds{npower};
    SO_phase_KDE_bias_difference(npower).bias_diff_mean = bias_mean;
    SO_phase_KDE_bias_difference(npower).bias_diff_CI = [bias_CI_lo; bias_CI_hi];
    SO_phase_KDE_bias_difference(npower).prop_mean = prop_mean;
    SO_phase_KDE_bias_difference(npower).prop_CI = [prop_CI_lo; prop_CI_hi];
    SO_phase_KDE_bias_difference(npower).thresholds = thresholds;
    SO_phase_KDE_bias_difference(npower).bias_diff_shifted_mean = bias_shifted_mean;
    SO_phase_KDE_bias_difference(npower).bias_diff_shifted_CI = [bias_shifted_CI_lo; bias_shifted_CI_hi];
    SO_phase_KDE_bias_difference(npower).prop_shifted_mean = prop_shifted_mean;
    SO_phase_KDE_bias_difference(npower).prop_shifted_CI = [prop_shifted_CI_lo; prop_shifted_CI_hi];
    SO_phase_KDE_bias_difference(npower).total_events = total_events;

    % store AUC
    auc_boot = (trapz(thresholds, bias_diff_boot') / (max(thresholds)-min(thresholds)))';
    auc_shift_boot = (trapz(thresholds, bias_diff_shifted_boot') / (max(thresholds)-min(thresholds)))';

    SO_phase_KDE_bias_difference(npower).AUC_mean = mean(auc_boot, 'omitnan');
    SO_phase_KDE_bias_difference(npower).AUC_CI = prctile(auc_boot, [2.5 97.5]);
    SO_phase_KDE_bias_difference(npower).AUC_mean_shuffled = mean(auc_shift_boot, 'omitnan');
    SO_phase_KDE_bias_difference(npower).AUC_CI_shuffled = prctile(auc_shift_boot, [2.5 97.5]);

    % ----------- PLOTS (A, B, C) ----------------

    % A: Bias vs. Threshold
    nexttile((npower-1)*3 + 1);
    hold on;
    fill([thresholds, fliplr(thresholds)], [bias_CI_lo, fliplr(bias_CI_hi)], ...
        colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    plot(thresholds, bias_mean, 'Color', colour_lines(npower,:), 'LineWidth', 2);
    fill([thresholds, fliplr(thresholds)], ...
         [bias_shifted_CI_lo, fliplr(bias_shifted_CI_hi)], ...
         [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    plot(thresholds, bias_shifted_mean, 'k-', 'LineWidth', 1.5);
    ylim([-0.15 0.35]); xlim([0 1]); yline(0, '--r');
    xlabel('V1 bias threshold'); ylabel('HPC bias diff (T1 - T2)');
    title(sprintf('SO phase bin %s', SO_thresholds{npower}));
    set(gca,"TickDir","out",'box','off','Color','none','FontSize',12)

    % B: Proportion vs Bias Diff
    nexttile((npower-1)*3 + 2);
    hold on;
    valid = isfinite(bias_mean) & isfinite(prop_mean);
    fill([bias_CI_lo(valid), fliplr(bias_CI_hi(valid))], ...
         [prop_mean(valid), fliplr(prop_mean(valid))], ...
         colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    plot(bias_mean(valid), prop_mean(valid), '-', 'Color', colour_lines(npower,:), 'LineWidth', 2);
    fill([bias_shifted_CI_lo(valid), fliplr(bias_shifted_CI_hi(valid))], ...
         [prop_shifted_mean(valid), fliplr(prop_shifted_mean(valid))], ...
         [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    plot(bias_shifted_mean(valid), prop_shifted_mean(valid), 'k-', 'LineWidth', 1.5);
    xlim([-0.1 0.35]); xline(0, '--r');
    xlabel('HPC bias diff (T1 - T2)'); ylabel('Proportion of events detected');
    title('Event Proportion vs. Bias Difference');
    set(gca,"TickDir","out",'box','off','Color','none','FontSize',12)

    % C: Proportion CI vs Bias Diff
    nexttile((npower-1)*3 + 3);
    hold on;
    fill([bias_mean(valid), fliplr(bias_mean(valid))], ...
         [prop_CI_lo(valid), fliplr(prop_CI_hi(valid))], ...
         colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    plot(bias_mean(valid), prop_mean(valid), '-', 'Color', colour_lines(npower,:), 'LineWidth', 2);
    fill([bias_shifted_mean(valid), fliplr(bias_shifted_mean(valid))], ...
         [prop_shifted_CI_lo(valid), fliplr(prop_shifted_CI_hi(valid))], ...
         [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    plot(bias_shifted_mean(valid), prop_shifted_mean(valid), 'k-', 'LineWidth', 1.5);
    xlim([-0.1 0.35]); xline(0, '--r');
    xlabel('HPC bias diff (T1 - T2)'); ylabel('Proportion of events detected');
    title('Proportion vs. HPC Bias Difference');
    set(gca,"TickDir","out",'box','off','Color','none','FontSize',12)
    toc
end

clear Fill

fig = figure;
fig.Position = [640 100 2*1100/3 650/2];
% fig.Name = 'KDE bias difference in HPC low vs high SO power (PRE ripple)';
fig.Name = 'KDE bias difference in HPC SO peak synchrony';

% Select bins (1 = low, 4 = high)
nexttile;
for npower = [1 2]
    bias_mean = SO_phase_KDE_bias_difference(npower).bias_diff_mean;
    bias_CI_lo = SO_phase_KDE_bias_difference(npower).bias_diff_CI(1,:);
    bias_CI_hi = SO_phase_KDE_bias_difference(npower).bias_diff_CI(2,:);
    thresholds = SO_phase_KDE_bias_difference(npower).thresholds;

    hold on;
    x2 = [thresholds, fliplr(thresholds)];
    y2 = [bias_CI_lo, fliplr(bias_CI_hi)];
    Fill(npower) = fill(x2, y2, colour_lines(npower,:), ...
        'EdgeColor', 'none', 'FaceAlpha', 0.3);
    plot(thresholds, bias_mean, 'Color', colour_lines(npower,:), 'LineWidth', 2);

    xlabel('V1 Bias threshold');
    ylabel('HPC bias diff (T1 - T2)');
    set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12);
    ylim([-0.1 0.35]);
    % ylim([0 1])
end

bias_mean = SO_phase_KDE_bias_difference(2).bias_diff_shifted_mean;
bias_CI_lo = SO_phase_KDE_bias_difference(2).bias_diff_shifted_CI(1,:);
bias_CI_hi = SO_phase_KDE_bias_difference(2).bias_diff_shifted_CI(2,:)

hold on;
x2 = [thresholds, fliplr(thresholds)];
y2 = [bias_CI_lo, fliplr(bias_CI_hi)];
Fill(end +1) = fill(x2, y2, 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
plot(thresholds, bias_mean, 'Color','k', 'LineWidth', 2);

yline(0, '--r');
legend(Fill([1 2 3]), {'SO peak not sync', 'SO peak sync','Shuffled'}, 'box', 'off');

nexttile;
for npower = [1 2]
    bias_mean = SO_phase_KDE_bias_difference(npower).bias_diff_mean;
    bias_CI_lo = SO_phase_KDE_bias_difference(npower).bias_diff_CI(1,:);
    bias_CI_hi = SO_phase_KDE_bias_difference(npower).bias_diff_CI(2,:);
    prop_mean = SO_phase_KDE_bias_difference(npower).prop_mean;

    hold on;
    y2 = [prop_mean, fliplr(prop_mean)];
    x2 = [bias_CI_lo, fliplr(bias_CI_hi)];
    Fill(npower) = fill(x2, y2, colour_lines(npower,:), ...
        'EdgeColor', 'none', 'FaceAlpha', 0.3);
    plot(bias_mean, prop_mean, 'Color', colour_lines(npower,:), 'LineWidth', 2);

    xlabel('HPC bias diff (T1 - T2)');
    ylabel('Proportion of events detected');
    set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12);
    xlim([-0.1 0.35]);
    ylim([0 1])
end


bias_mean = SO_phase_KDE_bias_difference(2).bias_diff_shifted_mean;
bias_CI_lo = SO_phase_KDE_bias_difference(2).bias_diff_shifted_CI(1,:);
bias_CI_hi = SO_phase_KDE_bias_difference(2).bias_diff_shifted_CI(2,:)
prop_mean = SO_phase_KDE_bias_difference(2).prop_shifted_mean;

hold on;
y2 = [prop_mean, fliplr(prop_mean)];
x2 = [bias_CI_lo, fliplr(bias_CI_hi)];
Fill(end + 1) = fill(x2, y2, 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
plot(bias_mean, prop_mean, 'Color','k', 'LineWidth', 2);


xline(0, '--r');
legend(Fill([1 2 3]), {'SO peak not sync', 'SO peak sync','Shuffled'}, 'box', 'off');

% Save results
save('SO_peak_KDE_bias_difference_based_on_V1_bias.mat', 'SO_phase_KDE_bias_difference');




%%%%% AUC mean + CI bar plot
% Plot layout
fig = figure;
fig.Position = [640 100 281 325]
fig.Name = 'KDE bias V1 AUC SO peak synchrony';
data = SO_phase_KDE_bias_difference;
n_bins = length(data);
bar_width = 0.3;      % Width of the bars
group_offset = 0.15;    % Distance from the center integer (half the gap between bars)
hold on;
clear BAR
for i = 1:4
    if i > length(data)
        x_shuf = i;
        y_shuf = data(1).AUC_mean_shuffled;
        bar(x_shuf, y_shuf, bar_width, ...
            'FaceColor', 'k', ...
            'FaceAlpha', 0.15, ...
            'EdgeColor', 'none');
    else


        % --- 1. Plot Shuffled Data (Left Bar) ---
        x_shuf = i - group_offset;
        y_shuf = data(i).AUC_mean_shuffled;

        % Calculate Error Deltas (Errorbar requires length relative to mean, not absolute values)
        % CI is [lower, upper]
        neg_err_shuf = y_shuf - data(i).AUC_CI_shuffled(1);
        pos_err_shuf = data(i).AUC_CI_shuffled(2) - y_shuf;

        % Plot Bar
        if i == 1
            BAR(1) = bar(x_shuf, y_shuf, bar_width, ...
                'FaceColor', 'k', ...
                'FaceAlpha', 0.15, ...
                'EdgeColor', 'none');
        else
            bar(x_shuf, y_shuf, bar_width, ...
                'FaceColor', 'k', ...
                'FaceAlpha', 0.15, ...
                'EdgeColor', 'none');
        end

        % Plot Error Bar
        E = errorbar(x_shuf, y_shuf, neg_err_shuf, pos_err_shuf, ...
            'Color', 'k', 'LineWidth', 1.5, 'CapSize', 8, 'LineStyle', 'none');


        % --- 2. Plot Real Data (Right Bar) ---
        x_real = i + group_offset;
        y_real = data(i).AUC_mean;

        neg_err_real = y_real - data(i).AUC_CI(1);
        pos_err_real = data(i).AUC_CI(2) - y_real;

        % Get specific color for this power bin
        % Assuming colour_lines is size [n_bins x 3]
        this_color = colour_lines(i, :);

        % Plot Bar
        BAR(i+1) = bar(x_real, y_real, bar_width, ...
            'FaceColor', this_color, ...
            'FaceAlpha', 0.3, ...
            'EdgeColor', 'none');

        errorbar(x_real, y_real, neg_err_real, pos_err_real, ...
            'Color', 'k', 'LineWidth', 1.5, 'CapSize', 8, 'LineStyle', 'none');
    end
end

hold off;
% Set X-ticks to be centered on the groups
set(gca, 'XTick', 1:nBins);
xlim([0.5, 4 + 0.5]);

% Labels
ylabel('V1 bias AUC');
xlabel('Power Bins');
legend([BAR(1:end)],{'Shuffled','SO peak unilateral','SO peak bilateral'},'box','off')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE bias difference based on V1 bias'),[])






%%%%%%%%%%%%%%%%%%%%%%%% SO peak synchrony PRE
% Spindle power binning across both probes
% all_spindle_power = mean(ripple_info.spindle_amplitude, 1);  % avg of probe 1 and 2
% SO_thresholds = prctile(ripple_info.SO_phase, 0:99.9/4:99.9);
SO_thresholds = {'peak not sync','peak sync'};
nBins = 2;

bins_to_use = bin_centers>0 & bin_centers<0.1;
bins_to_select = bin_centers>-0.2 & bin_centers<0;
bins_to_use_shifted = bin_centers > -1 & bin_centers < -0.9;
nBoot = 1000;

SO_power_KDE_bias_difference = struct;

colour_lines = [ ...
    241, 182, 218;   % original end (lightest)
    231, 41, 138    % original start (darkest)
] / 256;
fig = figure;
fig.Position = [640 100 1100 650];
fig.Name = 'KDE bias difference in HPC with different SO peak synchrony (PRE ripple)';
tiledlayout(nBins, 3, 'TileSpacing', 'compact');

for npower = 1:nBins
    tic
    % Select events in current spindle bin
    % power_index = ripple_info.spindle_amplitude(:,1) > spindle_thresholds(npower,1) & ...
    %               ripple_info.spindle_amplitude(:,1) <= spindle_thresholds(npower+1,1) |...
    %               ripple_info.spindle_amplitude(:,2) > spindle_thresholds(npower,2) & ...
    %               ripple_info.spindle_amplitude(:,2) <= spindle_thresholds(npower+1,2);

    event_index =1:length(z_bias);
    % event_index = find(power_index);

    mean_bias = mean(z_bias_V1(bins_to_select, event_index), 'omitnan');
    % mean_bias_shifted = mean(z_bias(bins_to_use_shifted, event_index), 'omitnan');
    mean_bias_V1 = mean(z_bias(bins_to_use, event_index), 'omitnan');
    selected_events = length(mean_bias);

    thresholds = prctile(abs(mean_bias), 0:10:100);
    thresholds = thresholds(1:end-1);
    nThresh = length(thresholds);

    bias_diff_boot = NaN(nBoot, nThresh);
    prop_events_boot = NaN(nBoot, nThresh);
    bias_diff_shifted_boot = NaN(nBoot, nThresh);
    prop_events_shifted_boot = NaN(nBoot, nThresh);


    % event_phase = ripple_info.SO_phase';


    parfor iBoot = 1:nBoot
        s = RandStream('philox4x32_10', 'Seed', iBoot);
        idx = randi(s, selected_events, selected_events, 1);

        true_idx = find(event_index);

        bb = mean_bias(idx);
        bb_shift = mean_bias;
        boot_V1 = mean_bias_V1(idx);

        diff_tmp = NaN(1, nThresh);
        prop_tmp = NaN(1, nThresh);
        diff_tmp_shifted = NaN(1, nThresh);
        prop_tmp_shifted = NaN(1, nThresh);
        event_phase = ripple_info.SO_phase(true_idx(idx),:)';
        event_phase_shifted = ripple_info.SO_phase(true_idx,:)';

        for i = 1:nThresh
            th = thresholds(i);

            t1 = bb >= th;
            t2 = bb <= -th;

            if npower == 1 % if peak not sync


                t1 = t1 & (event_phase(2,:) >= -pi/2 & event_phase(2,:) <= pi/2) & (event_phase(1,:) >= -pi & event_phase(1,:) <= -pi/2 | event_phase(1,:) >= pi/2 & event_phase(1,:) <= pi);
                t2 = t2 & (event_phase(1,:) >= -pi/2 & event_phase(1,:) <= pi/2) & (event_phase(2,:) >= -pi & event_phase(2,:) <= -pi/2 | event_phase(2,:) >= pi/2 & event_phase(2,:) <= pi);

                t1_V1 = boot_V1(t1);
                t2_V1 = boot_V1(t2);
                if any(t1) && any(t2)
                    diff_tmp(i) = mean(t1_V1) - mean(t2_V1);
                end

                total_events = mean([sum((event_phase(2,:) >= -pi/2 & event_phase(2,:) <= pi/2) & (event_phase(1,:) >= -pi & event_phase(1,:) <= -pi/2 | event_phase(1,:) >= pi/2 & event_phase(1,:) <= pi)) ...
                    sum((event_phase(1,:) >= -pi/2 & event_phase(1,:) <= pi/2) & (event_phase(2,:) >= -pi & event_phase(2,:) <= -pi/2 | event_phase(2,:) >= pi/2 & event_phase(2,:) <= pi))]);

                prop_tmp(i) = (sum(t1) + sum(t2)) / total_events;



                t1s = bb_shift >= th;
                t2s = bb_shift <= -th;
                t1s = t1s & (event_phase_shifted(2,:) >= -pi/2 & event_phase_shifted(2,:) <= pi/2) & (event_phase_shifted(1,:) >= -pi & event_phase_shifted(1,:) <= -pi/2 | event_phase_shifted(1,:) >= pi/2 & event_phase_shifted(1,:) <= pi);
                t2s = t2s & (event_phase_shifted(1,:) >= -pi/2 & event_phase_shifted(1,:) <= pi/2) & (event_phase_shifted(2,:) >= -pi & event_phase_shifted(2,:) <= -pi/2 | event_phase_shifted(2,:) >= pi/2 & event_phase_shifted(2,:) <= pi);

                t1_V1 = boot_V1(t1s);
                t2_V1 = boot_V1(t2s);

                if any(t1s) && any(t2s)
                    diff_tmp_shifted(i) = mean(t1_V1) - mean(t2_V1);
                end
                prop_tmp_shifted(i) = (sum(t1s) + sum(t2s)) / total_events;

            elseif npower == 2 % if phase trough

                t1 = t1 & (event_phase(2,:) >= -pi/2 & event_phase(2,:) <= pi/2) & (event_phase(1,:) >= -pi/2 & event_phase(1,:) <= pi/2);
                t2 = t2 & (event_phase(2,:) >= -pi/2 & event_phase(2,:) <= pi/2) & (event_phase(1,:) >= -pi/2 & event_phase(1,:) <= pi/2);

                t1_V1 = boot_V1(t1);
                t2_V1 = boot_V1(t2);
                if any(t1) && any(t2)
                    diff_tmp(i) = mean(t1_V1) - mean(t2_V1);
                end

                total_events = sum((event_phase(2,:) >= -pi/2 & event_phase(2,:) <= pi/2) & (event_phase(1,:) >= -pi/2 & event_phase(1,:) <= pi/2));

                prop_tmp(i) = (sum(t1) + sum(t2)) / total_events;



                t1s = bb_shift >= th;
                t2s = bb_shift <= -th;
                t1s = t1s & (event_phase_shifted(2,:) >= -pi/2 & event_phase_shifted(2,:) <= pi/2) & (event_phase_shifted(1,:) >= -pi/2 & event_phase_shifted(1,:) <= pi/2);
                t2s = t2s & (event_phase_shifted(2,:) >= -pi/2 & event_phase_shifted(2,:) <= pi/2) & (event_phase_shifted(1,:) >= -pi/2 & event_phase_shifted(1,:) <= pi/2);


                t1_V1 = boot_V1(t1s);
                t2_V1 = boot_V1(t2s);

                if any(t1s) && any(t2s)
                    diff_tmp_shifted(i) = mean(t1_V1) - mean(t2_V1);
                end
                prop_tmp_shifted(i) = (sum(t1s) + sum(t2s)) / total_events;

            end



        end

        bias_diff_boot(iBoot, :) = diff_tmp;
        prop_events_boot(iBoot, :) = prop_tmp;
        bias_diff_shifted_boot(iBoot, :) = diff_tmp_shifted;
        prop_events_shifted_boot(iBoot, :) = prop_tmp_shifted;
    end

    % Stats
    bias_mean = mean(bias_diff_boot, 1, 'omitnan');
    bias_CI_lo = prctile(bias_diff_boot, 2.5, 1);
    bias_CI_hi = prctile(bias_diff_boot, 97.5, 1);
    prop_mean = mean(prop_events_boot, 1, 'omitnan');
    prop_CI_lo = prctile(prop_events_boot, 2.5, 1);
    prop_CI_hi = prctile(prop_events_boot, 97.5, 1);
    bias_shifted_mean = mean(bias_diff_shifted_boot, 1, 'omitnan');
    bias_shifted_CI_lo = prctile(bias_diff_shifted_boot, 2.5, 1);
    bias_shifted_CI_hi = prctile(bias_diff_shifted_boot, 97.5, 1);
    prop_shifted_mean = mean(prop_events_shifted_boot, 1, 'omitnan');
    prop_shifted_CI_lo = prctile(prop_events_shifted_boot, 2.5, 1);
    prop_shifted_CI_hi = prctile(prop_events_shifted_boot, 97.5, 1);

    % Store
    SO_phase_KDE_bias_difference(npower).phase_range{npower} = ...
        SO_thresholds{npower};
    SO_phase_KDE_bias_difference(npower).bias_diff_mean = bias_mean;
    SO_phase_KDE_bias_difference(npower).bias_diff_CI = [bias_CI_lo; bias_CI_hi];
    SO_phase_KDE_bias_difference(npower).prop_mean = prop_mean;
    SO_phase_KDE_bias_difference(npower).prop_CI = [prop_CI_lo; prop_CI_hi];
    SO_phase_KDE_bias_difference(npower).thresholds = thresholds;
    SO_phase_KDE_bias_difference(npower).bias_diff_shifted_mean = bias_shifted_mean;
    SO_phase_KDE_bias_difference(npower).bias_diff_shifted_CI = [bias_shifted_CI_lo; bias_shifted_CI_hi];
    SO_phase_KDE_bias_difference(npower).prop_shifted_mean = prop_shifted_mean;
    SO_phase_KDE_bias_difference(npower).prop_shifted_CI = [prop_shifted_CI_lo; prop_shifted_CI_hi];
    SO_phase_KDE_bias_difference(npower).total_events = total_events;


        % store AUC
    auc_boot = (trapz(thresholds, bias_diff_boot') / (max(thresholds)-min(thresholds)))';
    auc_shift_boot = (trapz(thresholds, bias_diff_shifted_boot') / (max(thresholds)-min(thresholds)))';

    SO_phase_KDE_bias_difference(npower).AUC_mean = mean(auc_boot, 'omitnan');
    SO_phase_KDE_bias_difference(npower).AUC_CI = prctile(auc_boot, [2.5 97.5]);
    SO_phase_KDE_bias_difference(npower).AUC_mean_shuffled = mean(auc_shift_boot, 'omitnan');
    SO_phase_KDE_bias_difference(npower).AUC_CI_shuffled = prctile(auc_shift_boot, [2.5 97.5]);

    % ----------- PLOTS (A, B, C) ----------------

    % A: Bias vs. Threshold
    nexttile((npower-1)*3 + 1);
    hold on;
    fill([thresholds, fliplr(thresholds)], [bias_CI_lo, fliplr(bias_CI_hi)], ...
        colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    plot(thresholds, bias_mean, 'Color', colour_lines(npower,:), 'LineWidth', 2);
    fill([thresholds, fliplr(thresholds)], ...
         [bias_shifted_CI_lo, fliplr(bias_shifted_CI_hi)], ...
         [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    plot(thresholds, bias_shifted_mean, 'k-', 'LineWidth', 1.5);
    ylim([-0.15 0.35]); xlim([0 1]); yline(0, '--r');
    xlabel('V1 bias threshold'); ylabel('HPC bias diff (T1 - T2)');
    title(sprintf('SO phase bin %s', SO_thresholds{npower}));
    set(gca,"TickDir","out",'box','off','Color','none','FontSize',12)

    % B: Proportion vs Bias Diff
    nexttile((npower-1)*3 + 2);
    hold on;
    valid = isfinite(bias_mean) & isfinite(prop_mean);
    fill([bias_CI_lo(valid), fliplr(bias_CI_hi(valid))], ...
         [prop_mean(valid), fliplr(prop_mean(valid))], ...
         colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    plot(bias_mean(valid), prop_mean(valid), '-', 'Color', colour_lines(npower,:), 'LineWidth', 2);
    fill([bias_shifted_CI_lo(valid), fliplr(bias_shifted_CI_hi(valid))], ...
         [prop_shifted_mean(valid), fliplr(prop_shifted_mean(valid))], ...
         [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    plot(bias_shifted_mean(valid), prop_shifted_mean(valid), 'k-', 'LineWidth', 1.5);
    xlim([-0.1 0.35]); xline(0, '--r');
    xlabel('HPC bias diff (T1 - T2)'); ylabel('Proportion of events detected');
    title('Event Proportion vs. Bias Difference');
    set(gca,"TickDir","out",'box','off','Color','none','FontSize',12)

    % C: Proportion CI vs Bias Diff
    nexttile((npower-1)*3 + 3);
    hold on;
    fill([bias_mean(valid), fliplr(bias_mean(valid))], ...
         [prop_CI_lo(valid), fliplr(prop_CI_hi(valid))], ...
         colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    plot(bias_mean(valid), prop_mean(valid), '-', 'Color', colour_lines(npower,:), 'LineWidth', 2);
    fill([bias_shifted_mean(valid), fliplr(bias_shifted_mean(valid))], ...
         [prop_shifted_CI_lo(valid), fliplr(prop_shifted_CI_hi(valid))], ...
         [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    plot(bias_shifted_mean(valid), prop_shifted_mean(valid), 'k-', 'LineWidth', 1.5);
    xlim([-0.1 0.35]); xline(0, '--r');
    xlabel('HPC bias diff (T1 - T2)'); ylabel('Proportion of events detected');
    title('Proportion vs. HPC Bias Difference');
    set(gca,"TickDir","out",'box','off','Color','none','FontSize',12)
    toc
end

clear Fill

fig = figure;
fig.Position = [640 100 2*1100/3 650/2];
% fig.Name = 'KDE bias difference in HPC low vs high SO power (PRE ripple)';
fig.Name = 'KDE bias difference in HPC SO peak synchrony (PRE ripple)';

% Select bins (1 = low, 4 = high)
nexttile;
for npower = [1 2]
    bias_mean = SO_phase_KDE_bias_difference(npower).bias_diff_mean;
    bias_CI_lo = SO_phase_KDE_bias_difference(npower).bias_diff_CI(1,:);
    bias_CI_hi = SO_phase_KDE_bias_difference(npower).bias_diff_CI(2,:);
    thresholds = SO_phase_KDE_bias_difference(npower).thresholds;

    hold on;
    x2 = [thresholds, fliplr(thresholds)];
    y2 = [bias_CI_lo, fliplr(bias_CI_hi)];
    Fill(npower) = fill(x2, y2, colour_lines(npower,:), ...
        'EdgeColor', 'none', 'FaceAlpha', 0.3);
    plot(thresholds, bias_mean, 'Color', colour_lines(npower,:), 'LineWidth', 2);

    xlabel('V1 Bias threshold');
    ylabel('HPC bias diff (T1 - T2)');
    set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12);
    ylim([-0.1 0.35]);
    % ylim([0 1])
end

bias_mean = SO_phase_KDE_bias_difference(2).bias_diff_shifted_mean;
bias_CI_lo = SO_phase_KDE_bias_difference(2).bias_diff_shifted_CI(1,:);
bias_CI_hi = SO_phase_KDE_bias_difference(2).bias_diff_shifted_CI(2,:)

hold on;
x2 = [thresholds, fliplr(thresholds)];
y2 = [bias_CI_lo, fliplr(bias_CI_hi)];
Fill(end +1) = fill(x2, y2, 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
plot(thresholds, bias_mean, 'Color','k', 'LineWidth', 2);

yline(0, '--r');
legend(Fill([1 2 3]), {'SO peak not sync', 'SO peak sync','Shuffled'}, 'box', 'off');

nexttile;
for npower = [1 2]
    bias_mean = SO_phase_KDE_bias_difference(npower).bias_diff_mean;
    bias_CI_lo = SO_phase_KDE_bias_difference(npower).bias_diff_CI(1,:);
    bias_CI_hi = SO_phase_KDE_bias_difference(npower).bias_diff_CI(2,:);
    prop_mean = SO_phase_KDE_bias_difference(npower).prop_mean;

    hold on;
    y2 = [prop_mean, fliplr(prop_mean)];
    x2 = [bias_CI_lo, fliplr(bias_CI_hi)];
    Fill(npower) = fill(x2, y2, colour_lines(npower,:), ...
        'EdgeColor', 'none', 'FaceAlpha', 0.3);
    plot(bias_mean, prop_mean, 'Color', colour_lines(npower,:), 'LineWidth', 2);

    xlabel('HPC bias diff (T1 - T2)');
    ylabel('Proportion of events detected');
    set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12);
    xlim([-0.1 0.35]);
    ylim([0 1])
end


bias_mean = SO_phase_KDE_bias_difference(2).bias_diff_shifted_mean;
bias_CI_lo = SO_phase_KDE_bias_difference(2).bias_diff_shifted_CI(1,:);
bias_CI_hi = SO_phase_KDE_bias_difference(2).bias_diff_shifted_CI(2,:)
prop_mean = SO_phase_KDE_bias_difference(2).prop_shifted_mean;

hold on;
y2 = [prop_mean, fliplr(prop_mean)];
x2 = [bias_CI_lo, fliplr(bias_CI_hi)];
Fill(end + 1) = fill(x2, y2, 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
plot(bias_mean, prop_mean, 'Color','k', 'LineWidth', 2);


xline(0, '--r');
legend(Fill([1 2 3]), {'SO peak not sync', 'SO peak sync','Shuffled'}, 'box', 'off');

% Save results
save('SO_peak_KDE_bias_difference_based_on_PRE_V1_bias.mat', 'SO_phase_KDE_bias_difference');




%%%%% AUC mean + CI bar plot
% Plot layout
fig = figure;
fig.Position = [640 100 281 325]
fig.Name = 'KDE bias V1 AUC SO peak synchrony (PRE)';
data = SO_phase_KDE_bias_difference;
n_bins = length(data);
bar_width = 0.3;      % Width of the bars
group_offset = 0.15;    % Distance from the center integer (half the gap between bars)
hold on;
clear BAR
for i = 1:4
    if i > length(data)
        x_shuf = i;
        y_shuf = data(1).AUC_mean_shuffled;
        bar(x_shuf, y_shuf, bar_width, ...
            'FaceColor', 'k', ...
            'FaceAlpha', 0.15, ...
            'EdgeColor', 'none');
    else


        % --- 1. Plot Shuffled Data (Left Bar) ---
        x_shuf = i - group_offset;
        y_shuf = data(i).AUC_mean_shuffled;

        % Calculate Error Deltas (Errorbar requires length relative to mean, not absolute values)
        % CI is [lower, upper]
        neg_err_shuf = y_shuf - data(i).AUC_CI_shuffled(1);
        pos_err_shuf = data(i).AUC_CI_shuffled(2) - y_shuf;

        % Plot Bar
        if i == 1
            BAR(1) = bar(x_shuf, y_shuf, bar_width, ...
                'FaceColor', 'k', ...
                'FaceAlpha', 0.15, ...
                'EdgeColor', 'none');
        else
            bar(x_shuf, y_shuf, bar_width, ...
                'FaceColor', 'k', ...
                'FaceAlpha', 0.15, ...
                'EdgeColor', 'none');
        end

        % Plot Error Bar
        E = errorbar(x_shuf, y_shuf, neg_err_shuf, pos_err_shuf, ...
            'Color', 'k', 'LineWidth', 1.5, 'CapSize', 8, 'LineStyle', 'none');


        % --- 2. Plot Real Data (Right Bar) ---
        x_real = i + group_offset;
        y_real = data(i).AUC_mean;

        neg_err_real = y_real - data(i).AUC_CI(1);
        pos_err_real = data(i).AUC_CI(2) - y_real;

        % Get specific color for this power bin
        % Assuming colour_lines is size [n_bins x 3]
        this_color = colour_lines(i, :);

        % Plot Bar
        BAR(i+1) = bar(x_real, y_real, bar_width, ...
            'FaceColor', this_color, ...
            'FaceAlpha', 0.3, ...
            'EdgeColor', 'none');

        errorbar(x_real, y_real, neg_err_real, pos_err_real, ...
            'Color', 'k', 'LineWidth', 1.5, 'CapSize', 8, 'LineStyle', 'none');
    end
end

hold off;
% Set X-ticks to be centered on the groups
set(gca, 'XTick', 1:nBins);
xlim([0.5, 4 + 0.5]);

% Labels
ylabel('V1 bias AUC');
xlabel('Power Bins');
legend([BAR(1:end)],{'Shuffled','SO peak unilateral','SO peak bilateral'},'box','off')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE bias difference based on V1 bias'),[])


%% SO trough synchrony
% Spindle power binning across both probes
% all_spindle_power = mean(ripple_info.spindle_amplitude, 1);  % avg of probe 1 and 2
% SO_thresholds = prctile(ripple_info.SO_phase, 0:99.9/4:99.9);
% Spindle power binning across both probes
% all_spindle_power = mean(ripple_info.spindle_amplitude, 1);  % avg of probe 1 and 2
% SO_thresholds = prctile(ripple_info.SO_phase, 0:99.9/4:99.9);
SO_thresholds = {'trough not sync','trough sync'};
nBins = 2;

bins_to_use = bin_centers>0 & bin_centers<0.1;
bins_to_select = bin_centers>0 & bin_centers<0.2;
bins_to_use_shifted = bin_centers > -1 & bin_centers < -0.9;
nBoot = 1000;

SO_power_KDE_bias_difference = struct;

colour_lines = [ ...
    241, 182, 218;   % original end (lightest)
    231, 41, 138    % original start (darkest)
] / 256;


fig = figure;
fig.Position = [640 100 1100 650];
fig.Name = 'KDE bias difference in HPC with different SO trough synchrony';
tiledlayout(nBins, 3, 'TileSpacing', 'compact');

for npower = 1:nBins
    tic
    % Select events in current spindle bin
    % power_index = ripple_info.spindle_amplitude(:,1) > spindle_thresholds(npower,1) & ...
    %               ripple_info.spindle_amplitude(:,1) <= spindle_thresholds(npower+1,1) |...
    %               ripple_info.spindle_amplitude(:,2) > spindle_thresholds(npower,2) & ...
    %               ripple_info.spindle_amplitude(:,2) <= spindle_thresholds(npower+1,2);

    event_index =1:length(z_bias);
    % event_index = find(power_index);

    mean_bias = mean(z_bias_V1(bins_to_select, event_index), 'omitnan');
    % mean_bias_shifted = mean(z_bias(bins_to_use_shifted, event_index), 'omitnan');
    mean_bias_V1 = mean(z_bias(bins_to_use, event_index), 'omitnan');
    selected_events = length(mean_bias);

    thresholds = prctile(abs(mean_bias), 0:10:100);
    thresholds = thresholds(1:end-1);
    nThresh = length(thresholds);

    bias_diff_boot = NaN(nBoot, nThresh);
    prop_events_boot = NaN(nBoot, nThresh);
    bias_diff_shifted_boot = NaN(nBoot, nThresh);
    prop_events_shifted_boot = NaN(nBoot, nThresh);


    % event_phase = ripple_info.SO_phase';


    parfor iBoot = 1:nBoot
        s = RandStream('philox4x32_10', 'Seed', iBoot);
        idx = randi(s, selected_events, selected_events, 1);

        true_idx = find(event_index);

        bb = mean_bias(idx);
        bb_shift = mean_bias;
        boot_V1 = mean_bias_V1(idx);

        diff_tmp = NaN(1, nThresh);
        prop_tmp = NaN(1, nThresh);
        diff_tmp_shifted = NaN(1, nThresh);
        prop_tmp_shifted = NaN(1, nThresh);
        event_phase = ripple_info.SO_phase(true_idx(idx),:)';
        event_phase_shifted = ripple_info.SO_phase(true_idx,:)';

        for i = 1:nThresh
            th = thresholds(i);

            t1 = bb >= th;
            t2 = bb <= -th;

            if npower == 1 % if trough not sync


                t1 = t1 & (event_phase(1,:) > -pi/2 & event_phase(1,:) < pi/2) & (event_phase(2,:) >= -pi & event_phase(2,:) < -pi/2 | event_phase(2,:) > pi/2 & event_phase(2,:) <= pi);
                t2 = t2 & (event_phase(2,:) > -pi/2 & event_phase(2,:) < pi/2) & (event_phase(1,:) >= -pi & event_phase(1,:) < -pi/2 | event_phase(1,:) > pi/2 & event_phase(1,:) <= pi);

                t1_V1 = boot_V1(t1);
                t2_V1 = boot_V1(t2);
                if any(t1) && any(t2)
                    diff_tmp(i) = mean(t1_V1) - mean(t2_V1);
                end

                total_events = mean([sum((event_phase(2,:) >= -pi/2 & event_phase(2,:) <= pi/2) & (event_phase(1,:) >= -pi & event_phase(1,:) <= -pi/2 | event_phase(1,:) >= pi/2 & event_phase(1,:) <= pi)) ...
                    sum((event_phase(1,:) >= -pi/2 & event_phase(1,:) <= pi/2) & (event_phase(2,:) >= -pi & event_phase(2,:) <= -pi/2 | event_phase(2,:) >= pi/2 & event_phase(2,:) <= pi))]);

                prop_tmp(i) = (sum(t1) + sum(t2)) / total_events;



                t1s = bb_shift >= th;
                t2s = bb_shift <= -th;
                t1s = t1s & (event_phase_shifted(1,:) >= -pi/2 & event_phase_shifted(1,:) <= pi/2) & (event_phase_shifted(2,:) >= -pi & event_phase_shifted(2,:) <= -pi/2 | event_phase_shifted(2,:) >= pi/2 & event_phase_shifted(2,:) <= pi);
                t2s = t2s & (event_phase_shifted(2,:) >= -pi/2 & event_phase_shifted(2,:) <= pi/2) & (event_phase_shifted(1,:) >= -pi & event_phase_shifted(1,:) <= -pi/2 | event_phase_shifted(1,:) >= pi/2 & event_phase_shifted(1,:) <= pi);

                t1_V1 = boot_V1(t1s);
                t2_V1 = boot_V1(t2s);

                if any(t1s) && any(t2s)
                    diff_tmp_shifted(i) = mean(t1_V1) - mean(t2_V1);
                end
                prop_tmp_shifted(i) = (sum(t1s) + sum(t2s)) / total_events;

            elseif npower == 2 % if phase trough

                t1 = t1 & (event_phase(2,:) >= -pi & event_phase(2,:) <= -pi/2 | event_phase(2,:) >= pi/2 & event_phase(2,:) <= pi) & (event_phase(1,:) >= -pi & event_phase(1,:) <= -pi/2 | event_phase(1,:) >= pi/2 & event_phase(1,:) <= pi);
                t2 = t2 & (event_phase(2,:) >= -pi & event_phase(2,:) <= -pi/2 | event_phase(2,:) >= pi/2 & event_phase(2,:) <= pi) & (event_phase(1,:) >= -pi & event_phase(1,:) <= -pi/2 | event_phase(1,:) >= pi/2 & event_phase(1,:) <= pi);

                t1_V1 = boot_V1(t1);
                t2_V1 = boot_V1(t2);
                if any(t1) && any(t2)
                    diff_tmp(i) = mean(t1_V1) - mean(t2_V1);
                end

                total_events = sum((event_phase(2,:) >= -pi & event_phase(2,:) <= -pi/2 | event_phase(2,:) >= pi/2 & event_phase(2,:) <= pi) & (event_phase(1,:) >= -pi & event_phase(1,:) <= -pi/2 | event_phase(1,:) >= pi/2 & event_phase(1,:) <= pi));

                prop_tmp(i) = (sum(t1) + sum(t2)) / total_events;



                t1s = bb_shift >= th;
                t2s = bb_shift <= -th;
                t1s = t1s & (event_phase_shifted(2,:) >= -pi & event_phase_shifted(2,:) <= -pi/2 | event_phase_shifted(2,:) >= pi/2 & event_phase_shifted(2,:) <= pi) & (event_phase_shifted(1,:) >= -pi & event_phase_shifted(1,:) <= -pi/2 | event_phase_shifted(1,:) >= pi/2 & event_phase_shifted(1,:) <= pi);
                t2s = t2s & (event_phase_shifted(2,:) >= -pi & event_phase_shifted(2,:) <= -pi/2 | event_phase_shifted(2,:) >= pi/2 & event_phase_shifted(2,:) <= pi) & (event_phase_shifted(1,:) >= -pi & event_phase_shifted(1,:) <= -pi/2 | event_phase_shifted(1,:) >= pi/2 & event_phase_shifted(1,:) <= pi);


                t1_V1 = boot_V1(t1s);
                t2_V1 = boot_V1(t2s);

                if any(t1s) && any(t2s)
                    diff_tmp_shifted(i) = mean(t1_V1) - mean(t2_V1);
                end
                prop_tmp_shifted(i) = (sum(t1s) + sum(t2s)) / total_events;

            end



        end

        bias_diff_boot(iBoot, :) = diff_tmp;
        prop_events_boot(iBoot, :) = prop_tmp;
        bias_diff_shifted_boot(iBoot, :) = diff_tmp_shifted;
        prop_events_shifted_boot(iBoot, :) = prop_tmp_shifted;
    end

    % Stats
    bias_mean = mean(bias_diff_boot, 1, 'omitnan');
    bias_CI_lo = prctile(bias_diff_boot, 2.5, 1);
    bias_CI_hi = prctile(bias_diff_boot, 97.5, 1);
    prop_mean = mean(prop_events_boot, 1, 'omitnan');
    prop_CI_lo = prctile(prop_events_boot, 2.5, 1);
    prop_CI_hi = prctile(prop_events_boot, 97.5, 1);
    bias_shifted_mean = mean(bias_diff_shifted_boot, 1, 'omitnan');
    bias_shifted_CI_lo = prctile(bias_diff_shifted_boot, 2.5, 1);
    bias_shifted_CI_hi = prctile(bias_diff_shifted_boot, 97.5, 1);
    prop_shifted_mean = mean(prop_events_shifted_boot, 1, 'omitnan');
    prop_shifted_CI_lo = prctile(prop_events_shifted_boot, 2.5, 1);
    prop_shifted_CI_hi = prctile(prop_events_shifted_boot, 97.5, 1);

    % Store
    SO_phase_KDE_bias_difference(npower).phase_range{npower} = ...
        SO_thresholds{npower};
    SO_phase_KDE_bias_difference(npower).bias_diff_mean = bias_mean;
    SO_phase_KDE_bias_difference(npower).bias_diff_CI = [bias_CI_lo; bias_CI_hi];
    SO_phase_KDE_bias_difference(npower).prop_mean = prop_mean;
    SO_phase_KDE_bias_difference(npower).prop_CI = [prop_CI_lo; prop_CI_hi];
    SO_phase_KDE_bias_difference(npower).thresholds = thresholds;
    SO_phase_KDE_bias_difference(npower).bias_diff_shifted_mean = bias_shifted_mean;
    SO_phase_KDE_bias_difference(npower).bias_diff_shifted_CI = [bias_shifted_CI_lo; bias_shifted_CI_hi];
    SO_phase_KDE_bias_difference(npower).prop_shifted_mean = prop_shifted_mean;
    SO_phase_KDE_bias_difference(npower).prop_shifted_CI = [prop_shifted_CI_lo; prop_shifted_CI_hi];
    SO_phase_KDE_bias_difference(npower).total_events = total_events;

    % store AUC
    auc_boot = (trapz(thresholds, bias_diff_boot') / (max(thresholds)-min(thresholds)))';
    auc_shift_boot = (trapz(thresholds, bias_diff_shifted_boot') / (max(thresholds)-min(thresholds)))';

    SO_phase_KDE_bias_difference(npower).AUC_mean = mean(auc_boot, 'omitnan');
    SO_phase_KDE_bias_difference(npower).AUC_CI = prctile(auc_boot, [2.5 97.5]);
    SO_phase_KDE_bias_difference(npower).AUC_mean_shuffled = mean(auc_shift_boot, 'omitnan');
    SO_phase_KDE_bias_difference(npower).AUC_CI_shuffled = prctile(auc_shift_boot, [2.5 97.5]);


    % ----------- PLOTS (A, B, C) ----------------

    % A: Bias vs. Threshold
    nexttile((npower-1)*3 + 1);
    hold on;
    fill([thresholds, fliplr(thresholds)], [bias_CI_lo, fliplr(bias_CI_hi)], ...
        colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    plot(thresholds, bias_mean, 'Color', colour_lines(npower,:), 'LineWidth', 2);
    fill([thresholds, fliplr(thresholds)], ...
         [bias_shifted_CI_lo, fliplr(bias_shifted_CI_hi)], ...
         [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    plot(thresholds, bias_shifted_mean, 'k-', 'LineWidth', 1.5);
    ylim([-0.15 0.35]); xlim([0 1]); yline(0, '--r');
    xlabel('V1 bias threshold'); ylabel('HPC bias diff (T1 - T2)');
    title(sprintf('SO phase bin %s', SO_thresholds{npower}));
    set(gca,"TickDir","out",'box','off','Color','none','FontSize',12)

    % B: Proportion vs Bias Diff
    nexttile((npower-1)*3 + 2);
    hold on;
    valid = isfinite(bias_mean) & isfinite(prop_mean);
    fill([bias_CI_lo(valid), fliplr(bias_CI_hi(valid))], ...
         [prop_mean(valid), fliplr(prop_mean(valid))], ...
         colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    plot(bias_mean(valid), prop_mean(valid), '-', 'Color', colour_lines(npower,:), 'LineWidth', 2);
    fill([bias_shifted_CI_lo(valid), fliplr(bias_shifted_CI_hi(valid))], ...
         [prop_shifted_mean(valid), fliplr(prop_shifted_mean(valid))], ...
         [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    plot(bias_shifted_mean(valid), prop_shifted_mean(valid), 'k-', 'LineWidth', 1.5);
    xlim([-0.1 0.35]); xline(0, '--r');
    xlabel('HPC bias diff (T1 - T2)'); ylabel('Proportion of events detected');
    title('Event Proportion vs. Bias Difference');
    set(gca,"TickDir","out",'box','off','Color','none','FontSize',12)

    % C: Proportion CI vs Bias Diff
    nexttile((npower-1)*3 + 3);
    hold on;
    fill([bias_mean(valid), fliplr(bias_mean(valid))], ...
         [prop_CI_lo(valid), fliplr(prop_CI_hi(valid))], ...
         colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    plot(bias_mean(valid), prop_mean(valid), '-', 'Color', colour_lines(npower,:), 'LineWidth', 2);
    fill([bias_shifted_mean(valid), fliplr(bias_shifted_mean(valid))], ...
         [prop_shifted_CI_lo(valid), fliplr(prop_shifted_CI_hi(valid))], ...
         [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    plot(bias_shifted_mean(valid), prop_shifted_mean(valid), 'k-', 'LineWidth', 1.5);
    xlim([-0.1 0.35]); xline(0, '--r');
    xlabel('HPC bias diff (T1 - T2)'); ylabel('Proportion of events detected');
    title('Proportion vs. HPC Bias Difference');
    set(gca,"TickDir","out",'box','off','Color','none','FontSize',12)
    toc
end

clear Fill

fig = figure;
fig.Position = [640 100 2*1100/3 650/2];
% fig.Name = 'KDE bias difference in HPC low vs high SO power (PRE ripple)';
fig.Name = 'KDE bias difference in HPC SO trough synchrony';

% Select bins (1 = low, 4 = high)
nexttile;
for npower = [1 2]
    bias_mean = SO_phase_KDE_bias_difference(npower).bias_diff_mean;
    bias_CI_lo = SO_phase_KDE_bias_difference(npower).bias_diff_CI(1,:);
    bias_CI_hi = SO_phase_KDE_bias_difference(npower).bias_diff_CI(2,:);
    thresholds = SO_phase_KDE_bias_difference(npower).thresholds;

    hold on;
    x2 = [thresholds, fliplr(thresholds)];
    y2 = [bias_CI_lo, fliplr(bias_CI_hi)];
    Fill(npower) = fill(x2, y2, colour_lines(npower,:), ...
        'EdgeColor', 'none', 'FaceAlpha', 0.3);
    plot(thresholds, bias_mean, 'Color', colour_lines(npower,:), 'LineWidth', 2);

    xlabel('V1 Bias threshold');
    ylabel('HPC bias diff (T1 - T2)');
    set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12);
    ylim([-0.1 0.35]);
    % ylim([0 1])
end

bias_mean = SO_phase_KDE_bias_difference(2).bias_diff_shifted_mean;
bias_CI_lo = SO_phase_KDE_bias_difference(2).bias_diff_shifted_CI(1,:);
bias_CI_hi = SO_phase_KDE_bias_difference(2).bias_diff_shifted_CI(2,:)

hold on;
x2 = [thresholds, fliplr(thresholds)];
y2 = [bias_CI_lo, fliplr(bias_CI_hi)];
Fill(end +1) = fill(x2, y2, 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
plot(thresholds, bias_mean, 'Color','k', 'LineWidth', 2);

yline(0, '--r');
legend(Fill([1 2 3]), {'SO trough not sync', 'SO trough sync','Shuffled'}, 'box', 'off');

nexttile;
for npower = [1 2]
    bias_mean = SO_phase_KDE_bias_difference(npower).bias_diff_mean;
    bias_CI_lo = SO_phase_KDE_bias_difference(npower).bias_diff_CI(1,:);
    bias_CI_hi = SO_phase_KDE_bias_difference(npower).bias_diff_CI(2,:);
    prop_mean = SO_phase_KDE_bias_difference(npower).prop_mean;

    hold on;
    y2 = [prop_mean, fliplr(prop_mean)];
    x2 = [bias_CI_lo, fliplr(bias_CI_hi)];
    Fill(npower) = fill(x2, y2, colour_lines(npower,:), ...
        'EdgeColor', 'none', 'FaceAlpha', 0.3);
    plot(bias_mean, prop_mean, 'Color', colour_lines(npower,:), 'LineWidth', 2);

    xlabel('HPC bias diff (T1 - T2)');
    ylabel('Proportion of events detected');
    set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12);
    xlim([-0.1 0.35]);
    ylim([0 1])
end


bias_mean = SO_phase_KDE_bias_difference(2).bias_diff_shifted_mean;
bias_CI_lo = SO_phase_KDE_bias_difference(2).bias_diff_shifted_CI(1,:);
bias_CI_hi = SO_phase_KDE_bias_difference(2).bias_diff_shifted_CI(2,:)
prop_mean = SO_phase_KDE_bias_difference(2).prop_shifted_mean;

hold on;
y2 = [prop_mean, fliplr(prop_mean)];
x2 = [bias_CI_lo, fliplr(bias_CI_hi)];
Fill(end + 1) = fill(x2, y2, 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
plot(bias_mean, prop_mean, 'Color','k', 'LineWidth', 2);


xline(0, '--r');
legend(Fill([1 2 3]), {'SO trough not sync', 'SO trough sync','Shuffled'}, 'box', 'off');

% Save results
save('SO_trough_KDE_bias_difference_based_on_V1_bias.mat', 'SO_phase_KDE_bias_difference');



%%%%% AUC mean + CI bar plot
% Plot layout
fig = figure;
fig.Position = [640 100 281 325]
fig.Name = 'KDE bias V1 AUC SO trough synchrony';
data = SO_phase_KDE_bias_difference;
n_bins = length(data);
bar_width = 0.3;      % Width of the bars
group_offset = 0.15;    % Distance from the center integer (half the gap between bars)
hold on;
clear BAR
for i = 1:4
    if i > length(data)
        x_shuf = i;
        y_shuf = data(1).AUC_mean_shuffled;
        bar(x_shuf, y_shuf, bar_width, ...
            'FaceColor', 'k', ...
            'FaceAlpha', 0.15, ...
            'EdgeColor', 'none');
    else


        % --- 1. Plot Shuffled Data (Left Bar) ---
        x_shuf = i - group_offset;
        y_shuf = data(i).AUC_mean_shuffled;

        % Calculate Error Deltas (Errorbar requires length relative to mean, not absolute values)
        % CI is [lower, upper]
        neg_err_shuf = y_shuf - data(i).AUC_CI_shuffled(1);
        pos_err_shuf = data(i).AUC_CI_shuffled(2) - y_shuf;

        % Plot Bar
        if i == 1
            BAR(1) = bar(x_shuf, y_shuf, bar_width, ...
                'FaceColor', 'k', ...
                'FaceAlpha', 0.15, ...
                'EdgeColor', 'none');
        else
            bar(x_shuf, y_shuf, bar_width, ...
                'FaceColor', 'k', ...
                'FaceAlpha', 0.15, ...
                'EdgeColor', 'none');
        end

        % Plot Error Bar
        E = errorbar(x_shuf, y_shuf, neg_err_shuf, pos_err_shuf, ...
            'Color', 'k', 'LineWidth', 1.5, 'CapSize', 8, 'LineStyle', 'none');


        % --- 2. Plot Real Data (Right Bar) ---
        x_real = i + group_offset;
        y_real = data(i).AUC_mean;

        neg_err_real = y_real - data(i).AUC_CI(1);
        pos_err_real = data(i).AUC_CI(2) - y_real;

        % Get specific color for this power bin
        % Assuming colour_lines is size [n_bins x 3]
        this_color = colour_lines(i, :);

        % Plot Bar
        BAR(i+1) = bar(x_real, y_real, bar_width, ...
            'FaceColor', this_color, ...
            'FaceAlpha', 0.3, ...
            'EdgeColor', 'none');

        errorbar(x_real, y_real, neg_err_real, pos_err_real, ...
            'Color', 'k', 'LineWidth', 1.5, 'CapSize', 8, 'LineStyle', 'none');
    end
end

hold off;
% Set X-ticks to be centered on the groups
set(gca, 'XTick', 1:nBins);
xlim([0.5, 4 + 0.5]);

% Labels
ylabel('V1 bias AUC');
xlabel('Power Bins');
legend([BAR(1:end)],{'Shuffled','SO trough unilateral','SO trough bilateral'},'box','off')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE bias difference based on V1 bias'),[])






%%%%%%%%%%%%%%%%%%%%%%%% SO phase synchrony
% Spindle power binning across both probes
% all_spindle_power = mean(ripple_info.spindle_amplitude, 1);  % avg of probe 1 and 2
% SO_thresholds = prctile(ripple_info.SO_phase, 0:99.9/4:99.9);
SO_thresholds = {'trough not sync','trough sync'};
nBins = 2;

bins_to_use = bin_centers>0 & bin_centers<0.1;
bins_to_select = bin_centers>-0.2 & bin_centers<0;
bins_to_use_shifted = bin_centers > -1 & bin_centers < -0.9;
nBoot = 1000;

SO_power_KDE_bias_difference = struct;

colour_lines = [ ...
    241, 182, 218;   % original end (lightest)
    231, 41, 138    % original start (darkest)
] / 256;



fig = figure;
fig.Position = [640 100 1100 650];
fig.Name = 'KDE bias difference in HPC with different SO trough synchrony (PRE ripple)';
tiledlayout(nBins, 3, 'TileSpacing', 'compact');

for npower = 1:nBins
    tic
    % Select events in current spindle bin
    % power_index = ripple_info.spindle_amplitude(:,1) > spindle_thresholds(npower,1) & ...
    %               ripple_info.spindle_amplitude(:,1) <= spindle_thresholds(npower+1,1) |...
    %               ripple_info.spindle_amplitude(:,2) > spindle_thresholds(npower,2) & ...
    %               ripple_info.spindle_amplitude(:,2) <= spindle_thresholds(npower+1,2);

    event_index =1:length(z_bias);
    % event_index = find(power_index);

    mean_bias = mean(z_bias_V1(bins_to_select, event_index), 'omitnan');
    % mean_bias_shifted = mean(z_bias(bins_to_use_shifted, event_index), 'omitnan');
    mean_bias_V1 = mean(z_bias(bins_to_use, event_index), 'omitnan');
    selected_events = length(mean_bias);

    thresholds = prctile(abs(mean_bias), 0:10:100);
    thresholds = thresholds(1:end-1);
    nThresh = length(thresholds);

    bias_diff_boot = NaN(nBoot, nThresh);
    prop_events_boot = NaN(nBoot, nThresh);
    bias_diff_shifted_boot = NaN(nBoot, nThresh);
    prop_events_shifted_boot = NaN(nBoot, nThresh);


    % event_phase = ripple_info.SO_phase';


    parfor iBoot = 1:nBoot
        s = RandStream('philox4x32_10', 'Seed', iBoot);
        idx = randi(s, selected_events, selected_events, 1);

        true_idx = find(event_index);

        bb = mean_bias(idx);
        bb_shift = mean_bias;
        boot_V1 = mean_bias_V1(idx);

        diff_tmp = NaN(1, nThresh);
        prop_tmp = NaN(1, nThresh);
        diff_tmp_shifted = NaN(1, nThresh);
        prop_tmp_shifted = NaN(1, nThresh);
        event_phase = ripple_info.SO_phase(true_idx(idx),:)';
        event_phase_shifted = ripple_info.SO_phase(true_idx,:)';

        for i = 1:nThresh
            th = thresholds(i);

            t1 = bb >= th;
            t2 = bb <= -th;

            if npower == 1 % if trough not sync


                t1 = t1 & (event_phase(1,:) >= -pi/2 & event_phase(1,:) <= pi/2) & (event_phase(2,:) >= -pi & event_phase(2,:) <= -pi/2 | event_phase(2,:) >= pi/2 & event_phase(2,:) <= pi);
                t2 = t2 & (event_phase(2,:) >= -pi/2 & event_phase(2,:) <= pi/2) & (event_phase(1,:) >= -pi & event_phase(1,:) <= -pi/2 | event_phase(1,:) >= pi/2 & event_phase(1,:) <= pi);

                t1_V1 = boot_V1(t1);
                t2_V1 = boot_V1(t2);
                if any(t1) && any(t2)
                    diff_tmp(i) = mean(t1_V1) - mean(t2_V1);
                end

                total_events = mean([sum((event_phase(2,:) >= -pi/2 & event_phase(2,:) <= pi/2) & (event_phase(1,:) >= -pi & event_phase(1,:) <= -pi/2 | event_phase(1,:) >= pi/2 & event_phase(1,:) <= pi)) ...
                    sum((event_phase(1,:) >= -pi/2 & event_phase(1,:) <= pi/2) & (event_phase(2,:) >= -pi & event_phase(2,:) <= -pi/2 | event_phase(2,:) >= pi/2 & event_phase(2,:) <= pi))]);

                prop_tmp(i) = (sum(t1) + sum(t2)) / total_events;



                t1s = bb_shift >= th;
                t2s = bb_shift <= -th;
                t1s = t1s & (event_phase_shifted(1,:) >= -pi/2 & event_phase_shifted(1,:) <= pi/2) & (event_phase_shifted(2,:) >= -pi & event_phase_shifted(2,:) <= -pi/2 | event_phase_shifted(2,:) >= pi/2 & event_phase_shifted(2,:) <= pi);
                t2s = t2s & (event_phase_shifted(2,:) >= -pi/2 & event_phase_shifted(2,:) <= pi/2) & (event_phase_shifted(1,:) >= -pi & event_phase_shifted(1,:) <= -pi/2 | event_phase_shifted(1,:) >= pi/2 & event_phase_shifted(1,:) <= pi);

                t1_V1 = boot_V1(t1s);
                t2_V1 = boot_V1(t2s);

                if any(t1s) && any(t2s)
                    diff_tmp_shifted(i) = mean(t1_V1) - mean(t2_V1);
                end
                prop_tmp_shifted(i) = (sum(t1s) + sum(t2s)) / total_events;

            elseif npower == 2 % if phase trough

                t1 = t1 & (event_phase(2,:) >= -pi & event_phase(2,:) <= -pi/2 | event_phase(2,:) >= pi/2 & event_phase(2,:) <= pi) & (event_phase(1,:) >= -pi & event_phase(1,:) <= -pi/2 | event_phase(1,:) >= pi/2 & event_phase(1,:) <= pi);
                t2 = t2 & (event_phase(2,:) >= -pi & event_phase(2,:) <= -pi/2 | event_phase(2,:) >= pi/2 & event_phase(2,:) <= pi) & (event_phase(1,:) >= -pi & event_phase(1,:) <= -pi/2 | event_phase(1,:) >= pi/2 & event_phase(1,:) <= pi);

                t1_V1 = boot_V1(t1);
                t2_V1 = boot_V1(t2);
                if any(t1) && any(t2)
                    diff_tmp(i) = mean(t1_V1) - mean(t2_V1);
                end

                total_events = sum((event_phase(2,:) >= -pi & event_phase(2,:) <= -pi/2 | event_phase(2,:) >= pi/2 & event_phase(2,:) <= pi) & (event_phase(1,:) >= -pi & event_phase(1,:) <= -pi/2 | event_phase(1,:) >= pi/2 & event_phase(1,:) <= pi));

                prop_tmp(i) = (sum(t1) + sum(t2)) / total_events;



                t1s = bb_shift >= th;
                t2s = bb_shift <= -th;
                t1s = t1s & (event_phase_shifted(2,:) >= -pi & event_phase_shifted(2,:) <= -pi/2 | event_phase_shifted(2,:) >= pi/2 & event_phase_shifted(2,:) <= pi) & (event_phase_shifted(1,:) >= -pi & event_phase_shifted(1,:) <= -pi/2 | event_phase_shifted(1,:) >= pi/2 & event_phase_shifted(1,:) <= pi);
                t2s = t2s & (event_phase_shifted(2,:) >= -pi & event_phase_shifted(2,:) <= -pi/2 | event_phase_shifted(2,:) >= pi/2 & event_phase_shifted(2,:) <= pi) & (event_phase_shifted(1,:) >= -pi & event_phase_shifted(1,:) <= -pi/2 | event_phase_shifted(1,:) >= pi/2 & event_phase_shifted(1,:) <= pi);


                t1_V1 = boot_V1(t1s);
                t2_V1 = boot_V1(t2s);

                if any(t1s) && any(t2s)
                    diff_tmp_shifted(i) = mean(t1_V1) - mean(t2_V1);
                end
                prop_tmp_shifted(i) = (sum(t1s) + sum(t2s)) / total_events;

            end



        end

        bias_diff_boot(iBoot, :) = diff_tmp;
        prop_events_boot(iBoot, :) = prop_tmp;
        bias_diff_shifted_boot(iBoot, :) = diff_tmp_shifted;
        prop_events_shifted_boot(iBoot, :) = prop_tmp_shifted;
    end

    % Stats
    bias_mean = mean(bias_diff_boot, 1, 'omitnan');
    bias_CI_lo = prctile(bias_diff_boot, 2.5, 1);
    bias_CI_hi = prctile(bias_diff_boot, 97.5, 1);
    prop_mean = mean(prop_events_boot, 1, 'omitnan');
    prop_CI_lo = prctile(prop_events_boot, 2.5, 1);
    prop_CI_hi = prctile(prop_events_boot, 97.5, 1);
    bias_shifted_mean = mean(bias_diff_shifted_boot, 1, 'omitnan');
    bias_shifted_CI_lo = prctile(bias_diff_shifted_boot, 2.5, 1);
    bias_shifted_CI_hi = prctile(bias_diff_shifted_boot, 97.5, 1);
    prop_shifted_mean = mean(prop_events_shifted_boot, 1, 'omitnan');
    prop_shifted_CI_lo = prctile(prop_events_shifted_boot, 2.5, 1);
    prop_shifted_CI_hi = prctile(prop_events_shifted_boot, 97.5, 1);

    % Store
    SO_phase_KDE_bias_difference(npower).phase_range{npower} = ...
        SO_thresholds{npower};
    SO_phase_KDE_bias_difference(npower).bias_diff_mean = bias_mean;
    SO_phase_KDE_bias_difference(npower).bias_diff_CI = [bias_CI_lo; bias_CI_hi];
    SO_phase_KDE_bias_difference(npower).prop_mean = prop_mean;
    SO_phase_KDE_bias_difference(npower).prop_CI = [prop_CI_lo; prop_CI_hi];
    SO_phase_KDE_bias_difference(npower).thresholds = thresholds;
    SO_phase_KDE_bias_difference(npower).bias_diff_shifted_mean = bias_shifted_mean;
    SO_phase_KDE_bias_difference(npower).bias_diff_shifted_CI = [bias_shifted_CI_lo; bias_shifted_CI_hi];
    SO_phase_KDE_bias_difference(npower).prop_shifted_mean = prop_shifted_mean;
    SO_phase_KDE_bias_difference(npower).prop_shifted_CI = [prop_shifted_CI_lo; prop_shifted_CI_hi];
    SO_phase_KDE_bias_difference(npower).total_events = total_events;

    % store AUC
    auc_boot = (trapz(thresholds, bias_diff_boot') / (max(thresholds)-min(thresholds)))';
    auc_shift_boot = (trapz(thresholds, bias_diff_shifted_boot') / (max(thresholds)-min(thresholds)))';

    SO_phase_KDE_bias_difference(npower).AUC_mean = mean(auc_boot, 'omitnan');
    SO_phase_KDE_bias_difference(npower).AUC_CI = prctile(auc_boot, [2.5 97.5]);
    SO_phase_KDE_bias_difference(npower).AUC_mean_shuffled = mean(auc_shift_boot, 'omitnan');
    SO_phase_KDE_bias_difference(npower).AUC_CI_shuffled = prctile(auc_shift_boot, [2.5 97.5]);

    % ----------- PLOTS (A, B, C) ----------------

    % A: Bias vs. Threshold
    nexttile((npower-1)*3 + 1);
    hold on;
    fill([thresholds, fliplr(thresholds)], [bias_CI_lo, fliplr(bias_CI_hi)], ...
        colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    plot(thresholds, bias_mean, 'Color', colour_lines(npower,:), 'LineWidth', 2);
    fill([thresholds, fliplr(thresholds)], ...
         [bias_shifted_CI_lo, fliplr(bias_shifted_CI_hi)], ...
         [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    plot(thresholds, bias_shifted_mean, 'k-', 'LineWidth', 1.5);
    ylim([-0.15 0.35]); xlim([0 1]); yline(0, '--r');
    xlabel('V1 bias threshold'); ylabel('HPC bias diff (T1 - T2)');
    title(sprintf('SO phase bin %s', SO_thresholds{npower}));
    set(gca,"TickDir","out",'box','off','Color','none','FontSize',12)

    % B: Proportion vs Bias Diff
    nexttile((npower-1)*3 + 2);
    hold on;
    valid = isfinite(bias_mean) & isfinite(prop_mean);
    fill([bias_CI_lo(valid), fliplr(bias_CI_hi(valid))], ...
         [prop_mean(valid), fliplr(prop_mean(valid))], ...
         colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    plot(bias_mean(valid), prop_mean(valid), '-', 'Color', colour_lines(npower,:), 'LineWidth', 2);
    fill([bias_shifted_CI_lo(valid), fliplr(bias_shifted_CI_hi(valid))], ...
         [prop_shifted_mean(valid), fliplr(prop_shifted_mean(valid))], ...
         [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    plot(bias_shifted_mean(valid), prop_shifted_mean(valid), 'k-', 'LineWidth', 1.5);
    xlim([-0.1 0.35]); xline(0, '--r');
    xlabel('HPC bias diff (T1 - T2)'); ylabel('Proportion of events detected');
    title('Event Proportion vs. Bias Difference');
    set(gca,"TickDir","out",'box','off','Color','none','FontSize',12)

    % C: Proportion CI vs Bias Diff
    nexttile((npower-1)*3 + 3);
    hold on;
    fill([bias_mean(valid), fliplr(bias_mean(valid))], ...
         [prop_CI_lo(valid), fliplr(prop_CI_hi(valid))], ...
         colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    plot(bias_mean(valid), prop_mean(valid), '-', 'Color', colour_lines(npower,:), 'LineWidth', 2);
    fill([bias_shifted_mean(valid), fliplr(bias_shifted_mean(valid))], ...
         [prop_shifted_CI_lo(valid), fliplr(prop_shifted_CI_hi(valid))], ...
         [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    plot(bias_shifted_mean(valid), prop_shifted_mean(valid), 'k-', 'LineWidth', 1.5);
    xlim([-0.1 0.35]); xline(0, '--r');
    xlabel('HPC bias diff (T1 - T2)'); ylabel('Proportion of events detected');
    title('Proportion vs. HPC Bias Difference');
    set(gca,"TickDir","out",'box','off','Color','none','FontSize',12)
    toc
end

clear Fill

fig = figure;
fig.Position = [640 100 2*1100/3 650/2];
% fig.Name = 'KDE bias difference in HPC low vs high SO power (PRE ripple)';
fig.Name = 'KDE bias difference in HPC SO trough synchrony (PRE ripple)';

% Select bins (1 = low, 4 = high)
nexttile;
for npower = [1 2]
    bias_mean = SO_phase_KDE_bias_difference(npower).bias_diff_mean;
    bias_CI_lo = SO_phase_KDE_bias_difference(npower).bias_diff_CI(1,:);
    bias_CI_hi = SO_phase_KDE_bias_difference(npower).bias_diff_CI(2,:);
    thresholds = SO_phase_KDE_bias_difference(npower).thresholds;

    hold on;
    x2 = [thresholds, fliplr(thresholds)];
    y2 = [bias_CI_lo, fliplr(bias_CI_hi)];
    Fill(npower) = fill(x2, y2, colour_lines(npower,:), ...
        'EdgeColor', 'none', 'FaceAlpha', 0.3);
    plot(thresholds, bias_mean, 'Color', colour_lines(npower,:), 'LineWidth', 2);

    xlabel('HPC Bias threshold');
    ylabel('V1 bias diff (T1 - T2)');
    set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12);
    ylim([-0.1 0.35]);
    % ylim([0 1])
end

bias_mean = SO_phase_KDE_bias_difference(2).bias_diff_shifted_mean;
bias_CI_lo = SO_phase_KDE_bias_difference(2).bias_diff_shifted_CI(1,:);
bias_CI_hi = SO_phase_KDE_bias_difference(2).bias_diff_shifted_CI(2,:)

hold on;
x2 = [thresholds, fliplr(thresholds)];
y2 = [bias_CI_lo, fliplr(bias_CI_hi)];
Fill(end +1) = fill(x2, y2, 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
plot(thresholds, bias_mean, 'Color','k', 'LineWidth', 2);

yline(0, '--r');
legend(Fill([1 2 3]), {'SO trough not sync', 'SO trough sync','Shuffled'}, 'box', 'off');

nexttile;
for npower = [1 2]
    bias_mean = SO_phase_KDE_bias_difference(npower).bias_diff_mean;
    bias_CI_lo = SO_phase_KDE_bias_difference(npower).bias_diff_CI(1,:);
    bias_CI_hi = SO_phase_KDE_bias_difference(npower).bias_diff_CI(2,:);
    prop_mean = SO_phase_KDE_bias_difference(npower).prop_mean;

    hold on;
    y2 = [prop_mean, fliplr(prop_mean)];
    x2 = [bias_CI_lo, fliplr(bias_CI_hi)];
    Fill(npower) = fill(x2, y2, colour_lines(npower,:), ...
        'EdgeColor', 'none', 'FaceAlpha', 0.3);
    plot(bias_mean, prop_mean, 'Color', colour_lines(npower,:), 'LineWidth', 2);

    xlabel('V1 bias diff (T1 - T2)');
    ylabel('Proportion of events detected');
    set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12);
    xlim([-0.1 0.35]);
    ylim([0 1])
end


bias_mean = SO_phase_KDE_bias_difference(2).bias_diff_shifted_mean;
bias_CI_lo = SO_phase_KDE_bias_difference(2).bias_diff_shifted_CI(1,:);
bias_CI_hi = SO_phase_KDE_bias_difference(2).bias_diff_shifted_CI(2,:)
prop_mean = SO_phase_KDE_bias_difference(2).prop_shifted_mean;

hold on;
y2 = [prop_mean, fliplr(prop_mean)];
x2 = [bias_CI_lo, fliplr(bias_CI_hi)];
Fill(end + 1) = fill(x2, y2, 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
plot(bias_mean, prop_mean, 'Color','k', 'LineWidth', 2);


xline(0, '--r');
legend(Fill([1 2 3]), {'SO trough not sync', 'SO trough sync','Shuffled'}, 'box', 'off');

% Save results
save('SO_trough_KDE_bias_difference_based_on_PRE_V1_bias.mat', 'SO_phase_KDE_bias_difference');




%%%%% AUC mean + CI bar plot
% Plot layout
fig = figure;
fig.Position = [640 100 281 325]
fig.Name = 'KDE bias V1 AUC SO trough synchrony (PRE)';
data = SO_phase_KDE_bias_difference;
n_bins = length(data);
bar_width = 0.3;      % Width of the bars
group_offset = 0.15;    % Distance from the center integer (half the gap between bars)
hold on;
clear BAR
for i = 1:4
    if i > length(data)
        x_shuf = i;
        y_shuf = data(1).AUC_mean_shuffled;
        bar(x_shuf, y_shuf, bar_width, ...
            'FaceColor', 'k', ...
            'FaceAlpha', 0.15, ...
            'EdgeColor', 'none');
    else


        % --- 1. Plot Shuffled Data (Left Bar) ---
        x_shuf = i - group_offset;
        y_shuf = data(i).AUC_mean_shuffled;

        % Calculate Error Deltas (Errorbar requires length relative to mean, not absolute values)
        % CI is [lower, upper]
        neg_err_shuf = y_shuf - data(i).AUC_CI_shuffled(1);
        pos_err_shuf = data(i).AUC_CI_shuffled(2) - y_shuf;

        % Plot Bar
        if i == 1
            BAR(1) = bar(x_shuf, y_shuf, bar_width, ...
                'FaceColor', 'k', ...
                'FaceAlpha', 0.15, ...
                'EdgeColor', 'none');
        else
            bar(x_shuf, y_shuf, bar_width, ...
                'FaceColor', 'k', ...
                'FaceAlpha', 0.15, ...
                'EdgeColor', 'none');
        end

        % Plot Error Bar
        E = errorbar(x_shuf, y_shuf, neg_err_shuf, pos_err_shuf, ...
            'Color', 'k', 'LineWidth', 1.5, 'CapSize', 8, 'LineStyle', 'none');


        % --- 2. Plot Real Data (Right Bar) ---
        x_real = i + group_offset;
        y_real = data(i).AUC_mean;

        neg_err_real = y_real - data(i).AUC_CI(1);
        pos_err_real = data(i).AUC_CI(2) - y_real;

        % Get specific color for this power bin
        % Assuming colour_lines is size [n_bins x 3]
        this_color = colour_lines(i, :);

        % Plot Bar
        BAR(i+1) = bar(x_real, y_real, bar_width, ...
            'FaceColor', this_color, ...
            'FaceAlpha', 0.3, ...
            'EdgeColor', 'none');

        errorbar(x_real, y_real, neg_err_real, pos_err_real, ...
            'Color', 'k', 'LineWidth', 1.5, 'CapSize', 8, 'LineStyle', 'none');
    end
end

hold off;
% Set X-ticks to be centered on the groups
set(gca, 'XTick', 1:nBins);
xlim([0.5, 4 + 0.5]);

% Labels
ylabel('V1 bias AUC');
xlabel('Power Bins');
legend([BAR(1:end)],{'Shuffled','SO trough unilateral','SO trough bilateral'},'box','off')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE bias difference based on V1 bias'),[])



%% log odds AUC with SO trough + spindle trough or spindle peak (V1→HPC)
% Spindle power binning across both probes
% all_spindle_power = mean(ripple_info.spindle_amplitude, 1);  % avg of probe 1 and 2
% SO_thresholds = prctile(ripple_info.SO_phase, 0:99.9/4:99.9);
SO_thresholds = {'SO trough spindle peak','SO trough spindle trough'};
nBins = 2;

bins_to_use = bin_centers>0 & bin_centers<0.1;
bins_to_select = bin_centers>0 & bin_centers<0.2;
bins_to_use_shifted = bin_centers > -1 & bin_centers < -0.9;
nBoot = 1000;

SO_phase_KDE_bias_difference = struct;

colour_lines = [ ...
    241, 182, 218;   % original end (lightest)
    231, 41, 138    % original start (darkest)
] / 256;



fig = figure;
fig.Position = [640 100 1100 650];
fig.Name = 'KDE bias difference in HPC with SO + spindle trough synchrony';
tiledlayout(nBins, 3, 'TileSpacing', 'compact');

for npower = 1:nBins
    tic
    % Select events in current spindle bin
    % power_index = ripple_info.spindle_amplitude(:,1) > spindle_thresholds(npower,1) & ...
    %               ripple_info.spindle_amplitude(:,1) <= spindle_thresholds(npower+1,1) |...
    %               ripple_info.spindle_amplitude(:,2) > spindle_thresholds(npower,2) & ...
    %               ripple_info.spindle_amplitude(:,2) <= spindle_thresholds(npower+1,2);

    event_index =1:length(z_bias);
    % event_index = find(power_index);

    mean_bias = mean(z_bias_V1(bins_to_select, event_index), 'omitnan');
    % mean_bias_shifted = mean(z_bias(bins_to_use_shifted, event_index), 'omitnan');
    mean_bias_V1 = mean(z_bias(bins_to_use, event_index), 'omitnan');
    selected_events = length(mean_bias);

    thresholds = prctile(abs(mean_bias), 0:10:100);
    thresholds = thresholds(1:end-1);
    nThresh = length(thresholds);

    bias_diff_boot = NaN(nBoot, nThresh);
    prop_events_boot = NaN(nBoot, nThresh);
    bias_diff_shifted_boot = NaN(nBoot, nThresh);
    prop_events_shifted_boot = NaN(nBoot, nThresh);


    % event_phase = ripple_info.SO_phase';


    parfor iBoot = 1:nBoot
        s = RandStream('philox4x32_10', 'Seed', iBoot);
        idx = randi(s, selected_events, selected_events, 1);

        true_idx = find(event_index);

        bb = mean_bias(idx);
        bb_shift = mean_bias;
        boot_V1 = mean_bias_V1(idx);

        diff_tmp = NaN(1, nThresh);
        prop_tmp = NaN(1, nThresh);
        diff_tmp_shifted = NaN(1, nThresh);
        prop_tmp_shifted = NaN(1, nThresh);
        event_phase = ripple_info.SO_phase(true_idx(idx),:)';
        event_phase_shifted = ripple_info.SO_phase(true_idx,:)';

        spindle_phase = ripple_info.spindle_phase(true_idx(idx),:)';
        spindle_phase_shifted = ripple_info.spindle_phase(true_idx,:)';

        for i = 1:nThresh
            th = thresholds(i);

            t1 = bb >= th;
            t2 = bb <= -th;

            if npower == 1 % if trough not sync


                t1 = t1 & (spindle_phase(2,:) > -pi/2 & spindle_phase(2,:) < pi/2) & (event_phase(2,:) >= -pi & event_phase(2,:) < -pi/2 | event_phase(2,:) > pi/2 & event_phase(2,:) <= pi);
                t2 = t2 & (spindle_phase(1,:) > -pi/2 & spindle_phase(1,:) < pi/2) & (event_phase(1,:) >= -pi & event_phase(1,:) < -pi/2 | event_phase(1,:) > pi/2 & event_phase(1,:) <= pi);

                t1_V1 = boot_V1(t1);
                t2_V1 = boot_V1(t2);
                if any(t1) && any(t2)
                    diff_tmp(i) = mean(t1_V1) - mean(t2_V1);
                end

                total_events = mean([sum((spindle_phase(2,:) > -pi/2 & spindle_phase(2,:) < pi/2) & (event_phase(2,:) >= -pi & event_phase(2,:) < -pi/2 | event_phase(2,:) > pi/2 & event_phase(2,:) <= pi)) ...
                    sum((spindle_phase(1,:) > -pi/2 & spindle_phase(1,:) < pi/2) & (event_phase(1,:) >= -pi & event_phase(1,:) < -pi/2 | event_phase(1,:) > pi/2 & event_phase(1,:) <= pi))]);

                prop_tmp(i) = (sum(t1) + sum(t2)) / total_events;



                t1s = bb_shift >= th;
                t2s = bb_shift <= -th;
                t1s = t1s & (spindle_phase_shifted(2,:) > -pi/2 & spindle_phase_shifted(2,:) < pi/2) & (event_phase_shifted(2,:) >= -pi & event_phase_shifted(2,:) < -pi/2 | event_phase_shifted(2,:) > pi/2 & event_phase_shifted(2,:) <= pi);
                t2s = t2s & (spindle_phase_shifted(1,:) > -pi/2 & spindle_phase_shifted(1,:) < pi/2) & (event_phase_shifted(1,:) >= -pi & event_phase_shifted(1,:) < -pi/2 | event_phase_shifted(1,:) > pi/2 & event_phase_shifted(1,:) <= pi);

                t1_V1 = boot_V1(t1s);
                t2_V1 = boot_V1(t2s);

                if any(t1s) && any(t2s)
                    diff_tmp_shifted(i) = mean(t1_V1) - mean(t2_V1);
                end
                prop_tmp_shifted(i) = (sum(t1s) + sum(t2s)) / total_events;

            elseif npower == 2 % if trough sync

                t1 = t1 & (spindle_phase(2,:) >= -pi & spindle_phase(2,:) <= -pi/2 | spindle_phase(2,:) >= pi/2 & spindle_phase(2,:) <= pi) & (event_phase(2,:) >= -pi & event_phase(2,:) <= -pi/2 | event_phase(2,:) >= pi/2 & event_phase(2,:) <= pi);
                t2 = t2 & (spindle_phase(1,:) >= -pi & spindle_phase(1,:) <= -pi/2 | spindle_phase(1,:) >= pi/2 & spindle_phase(1,:) <= pi) & (event_phase(1,:) >= -pi & event_phase(1,:) <= -pi/2 | event_phase(1,:) >= pi/2 & event_phase(1,:) <= pi);

                t1_V1 = boot_V1(t1);
                t2_V1 = boot_V1(t2);
                if any(t1) && any(t2)
                    diff_tmp(i) = mean(t1_V1) - mean(t2_V1);
                end

                total_events = mean([sum((spindle_phase(2,:) >= -pi & spindle_phase(2,:) <= -pi/2 | spindle_phase(2,:) >= pi/2 & spindle_phase(2,:) <= pi) & (event_phase(2,:) >= -pi & event_phase(2,:) <= -pi/2 | event_phase(2,:) >= pi/2 & event_phase(2,:) <= pi)) ...
                    sum((spindle_phase(1,:) >= -pi & spindle_phase(1,:) <= -pi/2 | spindle_phase(1,:) >= pi/2 & spindle_phase(1,:) <= pi) & (event_phase(1,:) >= -pi & event_phase(1,:) <= -pi/2 | event_phase(1,:) >= pi/2 & event_phase(1,:) <= pi))]);

                prop_tmp(i) = (sum(t1) + sum(t2)) / total_events;

                t1s = bb_shift >= th;
                t2s = bb_shift <= -th;
                t1s = t1s & (spindle_phase_shifted(2,:) >= -pi & spindle_phase_shifted(2,:) <= -pi/2 | spindle_phase_shifted(2,:) >= pi/2 & spindle_phase_shifted(2,:) <= pi) & (event_phase_shifted(2,:) >= -pi & event_phase_shifted(2,:) <= -pi/2 | event_phase_shifted(2,:) >= pi/2 & event_phase_shifted(2,:) <= pi);
                t2s = t2s & (spindle_phase_shifted(1,:) >= -pi & spindle_phase_shifted(1,:) <= -pi/2 | spindle_phase_shifted(1,:) >= pi/2 & spindle_phase_shifted(1,:) <= pi) & (event_phase_shifted(1,:) >= -pi & event_phase_shifted(1,:) <= -pi/2 | event_phase_shifted(1,:) >= pi/2 & event_phase_shifted(1,:) <= pi);


                t1_V1 = boot_V1(t1s);
                t2_V1 = boot_V1(t2s);

                if any(t1s) && any(t2s)
                    diff_tmp_shifted(i) = mean(t1_V1) - mean(t2_V1);
                end
                prop_tmp_shifted(i) = (sum(t1s) + sum(t2s)) / total_events;

            end



        end

        bias_diff_boot(iBoot, :) = diff_tmp;
        prop_events_boot(iBoot, :) = prop_tmp;
        bias_diff_shifted_boot(iBoot, :) = diff_tmp_shifted;
        prop_events_shifted_boot(iBoot, :) = prop_tmp_shifted;
    end

    % Stats
    bias_mean = mean(bias_diff_boot, 1, 'omitnan');
    bias_CI_lo = prctile(bias_diff_boot, 2.5, 1);
    bias_CI_hi = prctile(bias_diff_boot, 97.5, 1);
    prop_mean = mean(prop_events_boot, 1, 'omitnan');
    prop_CI_lo = prctile(prop_events_boot, 2.5, 1);
    prop_CI_hi = prctile(prop_events_boot, 97.5, 1);
    bias_shifted_mean = mean(bias_diff_shifted_boot, 1, 'omitnan');
    bias_shifted_CI_lo = prctile(bias_diff_shifted_boot, 2.5, 1);
    bias_shifted_CI_hi = prctile(bias_diff_shifted_boot, 97.5, 1);
    prop_shifted_mean = mean(prop_events_shifted_boot, 1, 'omitnan');
    prop_shifted_CI_lo = prctile(prop_events_shifted_boot, 2.5, 1);
    prop_shifted_CI_hi = prctile(prop_events_shifted_boot, 97.5, 1);

    % Store
    SO_phase_KDE_bias_difference(npower).phase_range{npower} = ...
        SO_thresholds{npower};
    SO_phase_KDE_bias_difference(npower).bias_diff_mean = bias_mean;
    SO_phase_KDE_bias_difference(npower).bias_diff_CI = [bias_CI_lo; bias_CI_hi];
    SO_phase_KDE_bias_difference(npower).prop_mean = prop_mean;
    SO_phase_KDE_bias_difference(npower).prop_CI = [prop_CI_lo; prop_CI_hi];
    SO_phase_KDE_bias_difference(npower).thresholds = thresholds;
    SO_phase_KDE_bias_difference(npower).bias_diff_shifted_mean = bias_shifted_mean;
    SO_phase_KDE_bias_difference(npower).bias_diff_shifted_CI = [bias_shifted_CI_lo; bias_shifted_CI_hi];
    SO_phase_KDE_bias_difference(npower).prop_shifted_mean = prop_shifted_mean;
    SO_phase_KDE_bias_difference(npower).prop_shifted_CI = [prop_shifted_CI_lo; prop_shifted_CI_hi];
    SO_phase_KDE_bias_difference(npower).total_events = total_events;

    % store AUC
    auc_boot = (trapz(thresholds, bias_diff_boot') / (max(thresholds)-min(thresholds)))';
    auc_shift_boot = (trapz(thresholds, bias_diff_shifted_boot') / (max(thresholds)-min(thresholds)))';

    SO_phase_KDE_bias_difference(npower).AUC_mean = mean(auc_boot, 'omitnan');
    SO_phase_KDE_bias_difference(npower).AUC_CI = prctile(auc_boot, [2.5 97.5]);
    SO_phase_KDE_bias_difference(npower).AUC_mean_shuffled = mean(auc_shift_boot, 'omitnan');
    SO_phase_KDE_bias_difference(npower).AUC_CI_shuffled = prctile(auc_shift_boot, [2.5 97.5]);

    % ----------- PLOTS (A, B, C) ----------------

    % A: Bias vs. Threshold
    nexttile((npower-1)*3 + 1);
    hold on;
    fill([thresholds, fliplr(thresholds)], [bias_CI_lo, fliplr(bias_CI_hi)], ...
        colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    plot(thresholds, bias_mean, 'Color', colour_lines(npower,:), 'LineWidth', 2);
    fill([thresholds, fliplr(thresholds)], ...
         [bias_shifted_CI_lo, fliplr(bias_shifted_CI_hi)], ...
         [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    plot(thresholds, bias_shifted_mean, 'k-', 'LineWidth', 1.5);
    ylim([-0.15 0.35]); xlim([0 1]); yline(0, '--r');
    xlabel('V1 bias threshold'); ylabel('HPC bias diff (T1 - T2)');
    title(sprintf('SO phase bin %s', SO_thresholds{npower}));
    set(gca,"TickDir","out",'box','off','Color','none','FontSize',12)

    % B: Proportion vs Bias Diff
    nexttile((npower-1)*3 + 2);
    hold on;
    valid = isfinite(bias_mean) & isfinite(prop_mean);
    fill([bias_CI_lo(valid), fliplr(bias_CI_hi(valid))], ...
         [prop_mean(valid), fliplr(prop_mean(valid))], ...
         colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    plot(bias_mean(valid), prop_mean(valid), '-', 'Color', colour_lines(npower,:), 'LineWidth', 2);
    fill([bias_shifted_CI_lo(valid), fliplr(bias_shifted_CI_hi(valid))], ...
         [prop_shifted_mean(valid), fliplr(prop_shifted_mean(valid))], ...
         [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    plot(bias_shifted_mean(valid), prop_shifted_mean(valid), 'k-', 'LineWidth', 1.5);
    xlim([-0.1 0.35]); xline(0, '--r');
    xlabel('HPC bias diff (T1 - T2)'); ylabel('Proportion of events detected');
    title('Event Proportion vs. Bias Difference');
    set(gca,"TickDir","out",'box','off','Color','none','FontSize',12)

    % C: Proportion CI vs Bias Diff
    nexttile((npower-1)*3 + 3);
    hold on;
    fill([bias_mean(valid), fliplr(bias_mean(valid))], ...
         [prop_CI_lo(valid), fliplr(prop_CI_hi(valid))], ...
         colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    plot(bias_mean(valid), prop_mean(valid), '-', 'Color', colour_lines(npower,:), 'LineWidth', 2);
    fill([bias_shifted_mean(valid), fliplr(bias_shifted_mean(valid))], ...
         [prop_shifted_CI_lo(valid), fliplr(prop_shifted_CI_hi(valid))], ...
         [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    plot(bias_shifted_mean(valid), prop_shifted_mean(valid), 'k-', 'LineWidth', 1.5);
    xlim([-0.1 0.35]); xline(0, '--r');
    xlabel('HPC bias diff (T1 - T2)'); ylabel('Proportion of events detected');
    title('Proportion vs. HPC Bias Difference');
    set(gca,"TickDir","out",'box','off','Color','none','FontSize',12)
    toc
end

clear Fill

fig = figure;
fig.Position = [640 100 2*1100/3 650/2];
% fig.Name = 'KDE bias difference in HPC low vs high SO power (PRE ripple)';
fig.Name = 'KDE bias difference in HPC SO + spindle trough synchrony';

% Select bins (1 = low, 4 = high)
nexttile;
for npower = [1 2]
    bias_mean = SO_phase_KDE_bias_difference(npower).bias_diff_mean;
    bias_CI_lo = SO_phase_KDE_bias_difference(npower).bias_diff_CI(1,:);
    bias_CI_hi = SO_phase_KDE_bias_difference(npower).bias_diff_CI(2,:);
    thresholds = SO_phase_KDE_bias_difference(npower).thresholds;

    hold on;
    x2 = [thresholds, fliplr(thresholds)];
    y2 = [bias_CI_lo, fliplr(bias_CI_hi)];
    Fill(npower) = fill(x2, y2, colour_lines(npower,:), ...
        'EdgeColor', 'none', 'FaceAlpha', 0.3);
    plot(thresholds, bias_mean, 'Color', colour_lines(npower,:), 'LineWidth', 2);

    xlabel('HPC Bias threshold');
    ylabel('V1 bias diff (T1 - T2)');
    set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12);
    ylim([-0.1 0.35]);
    % ylim([0 1])
end

bias_mean = SO_phase_KDE_bias_difference(2).bias_diff_shifted_mean;
bias_CI_lo = SO_phase_KDE_bias_difference(2).bias_diff_shifted_CI(1,:);
bias_CI_hi = SO_phase_KDE_bias_difference(2).bias_diff_shifted_CI(2,:)

hold on;
x2 = [thresholds, fliplr(thresholds)];
y2 = [bias_CI_lo, fliplr(bias_CI_hi)];
Fill(end +1) = fill(x2, y2, 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
plot(thresholds, bias_mean, 'Color','k', 'LineWidth', 2);

yline(0, '--r');
legend(Fill([1 2 3]), {'SO trough not sync', 'SO trough sync','Shuffled'}, 'box', 'off');

nexttile;
for npower = [1 2]
    bias_mean = SO_phase_KDE_bias_difference(npower).bias_diff_mean;
    bias_CI_lo = SO_phase_KDE_bias_difference(npower).bias_diff_CI(1,:);
    bias_CI_hi = SO_phase_KDE_bias_difference(npower).bias_diff_CI(2,:);
    prop_mean = SO_phase_KDE_bias_difference(npower).prop_mean;

    hold on;
    y2 = [prop_mean, fliplr(prop_mean)];
    x2 = [bias_CI_lo, fliplr(bias_CI_hi)];
    Fill(npower) = fill(x2, y2, colour_lines(npower,:), ...
        'EdgeColor', 'none', 'FaceAlpha', 0.3);
    plot(bias_mean, prop_mean, 'Color', colour_lines(npower,:), 'LineWidth', 2);

    xlabel('V1 bias diff (T1 - T2)');
    ylabel('Proportion of events detected');
    set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12);
    xlim([-0.1 0.35]);
    ylim([0 1])
end


bias_mean = SO_phase_KDE_bias_difference(2).bias_diff_shifted_mean;
bias_CI_lo = SO_phase_KDE_bias_difference(2).bias_diff_shifted_CI(1,:);
bias_CI_hi = SO_phase_KDE_bias_difference(2).bias_diff_shifted_CI(2,:)
prop_mean = SO_phase_KDE_bias_difference(2).prop_shifted_mean;

hold on;
y2 = [prop_mean, fliplr(prop_mean)];
x2 = [bias_CI_lo, fliplr(bias_CI_hi)];
Fill(end + 1) = fill(x2, y2, 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
plot(bias_mean, prop_mean, 'Color','k', 'LineWidth', 2);


xline(0, '--r');
legend(Fill([1 2 3]), {'SO trough not sync', 'SO trough sync','Shuffled'}, 'box', 'off');

% Save results
save(fullfile(analysis_folder,'V1-HPC sleep reactivation','SO_spindle_trough_KDE_bias_difference_based_on_V1_bias.mat'), 'SO_phase_KDE_bias_difference');
load(fullfile(analysis_folder,'V1-HPC sleep reactivation','SO_spindle_trough_KDE_bias_difference_based_on_V1_bias.mat'), 'SO_phase_KDE_bias_difference');




%%%%% AUC mean + CI bar plot
% Plot layout
fig = figure;
fig.Position = [640 100 281 325]
fig.Name = 'KDE bias V1 AUC SO + spindle trough synchrony';
data = SO_phase_KDE_bias_difference;
n_bins = length(data);
bar_width = 0.3;      % Width of the bars
group_offset = 0.15;    % Distance from the center integer (half the gap between bars)
hold on;
clear BAR
for i = 1:4
    if i > length(data)
        x_shuf = i;
        y_shuf = data(1).AUC_mean_shuffled;
        bar(x_shuf, y_shuf, bar_width, ...
            'FaceColor', 'k', ...
            'FaceAlpha', 0.15, ...
            'EdgeColor', 'none');
    else


        % --- 1. Plot Shuffled Data (Left Bar) ---
        x_shuf = i - group_offset;
        y_shuf = data(i).AUC_mean_shuffled;

        % Calculate Error Deltas (Errorbar requires length relative to mean, not absolute values)
        % CI is [lower, upper]
        neg_err_shuf = y_shuf - data(i).AUC_CI_shuffled(1);
        pos_err_shuf = data(i).AUC_CI_shuffled(2) - y_shuf;

        % Plot Bar
        if i == 1
            BAR(1) = bar(x_shuf, y_shuf, bar_width, ...
                'FaceColor', 'k', ...
                'FaceAlpha', 0.15, ...
                'EdgeColor', 'none');
        else
            bar(x_shuf, y_shuf, bar_width, ...
                'FaceColor', 'k', ...
                'FaceAlpha', 0.15, ...
                'EdgeColor', 'none');
        end

        % Plot Error Bar
        E = errorbar(x_shuf, y_shuf, neg_err_shuf, pos_err_shuf, ...
            'Color', 'k', 'LineWidth', 1.5, 'CapSize', 8, 'LineStyle', 'none');


        % --- 2. Plot Real Data (Right Bar) ---
        x_real = i + group_offset;
        y_real = data(i).AUC_mean;

        neg_err_real = y_real - data(i).AUC_CI(1);
        pos_err_real = data(i).AUC_CI(2) - y_real;

        % Get specific color for this power bin
        % Assuming colour_lines is size [n_bins x 3]
        this_color = colour_lines(i, :);

        % Plot Bar
        BAR(i+1) = bar(x_real, y_real, bar_width, ...
            'FaceColor', this_color, ...
            'FaceAlpha', 0.3, ...
            'EdgeColor', 'none');

        errorbar(x_real, y_real, neg_err_real, pos_err_real, ...
            'Color', 'k', 'LineWidth', 1.5, 'CapSize', 8, 'LineStyle', 'none');
    end
end

hold off;
% Set X-ticks to be centered on the groups
set(gca, 'XTick', 1:nBins);
xlim([0.5, 4 + 0.5]);

% Labels
ylabel('V1 bias AUC');
xlabel('Power Bins');
legend([BAR(1:end)],{'Shuffled','SO trough spindle peak','SO trough spindle trough'},'box','off')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE bias difference based on V1 bias'),[])



%%%%%% log odds AUC with SO trough + spindle trough or spindle peak (V1→HPC) (PRE)
% Spindle power binning across both probes
% all_spindle_power = mean(ripple_info.spindle_amplitude, 1);  % avg of probe 1 and 2
% SO_thresholds = prctile(ripple_info.SO_phase, 0:99.9/4:99.9);
SO_thresholds = {'SO trough spindle peak','SO trough spindle trough'};
nBins = 2;

bins_to_use = bin_centers>0 & bin_centers<0.1;
bins_to_select = bin_centers>-0.2 & bin_centers<0;
bins_to_use_shifted = bin_centers > -1 & bin_centers < -0.9;
nBoot = 1000;

SO_phase_KDE_bias_difference = struct;

colour_lines = [ ...
    241, 182, 218;   % original end (lightest)
    231, 41, 138    % original start (darkest)
] / 256;



fig = figure;
fig.Position = [640 100 1100 650];
fig.Name = 'KDE bias difference in HPC with SO + spindle trough synchrony (PRE ripple)';
tiledlayout(nBins, 3, 'TileSpacing', 'compact');

for npower = 1:nBins
    tic
    % Select events in current spindle bin
    % power_index = ripple_info.spindle_amplitude(:,1) > spindle_thresholds(npower,1) & ...
    %               ripple_info.spindle_amplitude(:,1) <= spindle_thresholds(npower+1,1) |...
    %               ripple_info.spindle_amplitude(:,2) > spindle_thresholds(npower,2) & ...
    %               ripple_info.spindle_amplitude(:,2) <= spindle_thresholds(npower+1,2);

    event_index =1:length(z_bias);
    % event_index = find(power_index);

    mean_bias = mean(z_bias_V1(bins_to_select, event_index), 'omitnan');
    % mean_bias_shifted = mean(z_bias(bins_to_use_shifted, event_index), 'omitnan');
    mean_bias_V1 = mean(z_bias(bins_to_use, event_index), 'omitnan');
    selected_events = length(mean_bias);

    thresholds = prctile(abs(mean_bias), 0:10:100);
    thresholds = thresholds(1:end-1);
    nThresh = length(thresholds);

    bias_diff_boot = NaN(nBoot, nThresh);
    prop_events_boot = NaN(nBoot, nThresh);
    bias_diff_shifted_boot = NaN(nBoot, nThresh);
    prop_events_shifted_boot = NaN(nBoot, nThresh);


    % event_phase = ripple_info.SO_phase';


    parfor iBoot = 1:nBoot
        s = RandStream('philox4x32_10', 'Seed', iBoot);
        idx = randi(s, selected_events, selected_events, 1);

        true_idx = find(event_index);

        bb = mean_bias(idx);
        bb_shift = mean_bias;
        boot_V1 = mean_bias_V1(idx);

        diff_tmp = NaN(1, nThresh);
        prop_tmp = NaN(1, nThresh);
        diff_tmp_shifted = NaN(1, nThresh);
        prop_tmp_shifted = NaN(1, nThresh);
        event_phase = ripple_info.SO_phase(true_idx(idx),:)';
        event_phase_shifted = ripple_info.SO_phase(true_idx,:)';

        spindle_phase = ripple_info.spindle_phase(true_idx(idx),:)';
        spindle_phase_shifted = ripple_info.spindle_phase(true_idx,:)';

        for i = 1:nThresh
            th = thresholds(i);

            t1 = bb >= th;
            t2 = bb <= -th;

            if npower == 1 % if trough not sync


                t1 = t1 & (spindle_phase(2,:) > -pi/2 & spindle_phase(2,:) < pi/2) & (event_phase(2,:) >= -pi & event_phase(2,:) < -pi/2 | event_phase(2,:) > pi/2 & event_phase(2,:) <= pi);
                t2 = t2 & (spindle_phase(1,:) > -pi/2 & spindle_phase(1,:) < pi/2) & (event_phase(1,:) >= -pi & event_phase(1,:) < -pi/2 | event_phase(1,:) > pi/2 & event_phase(1,:) <= pi);

                t1_V1 = boot_V1(t1);
                t2_V1 = boot_V1(t2);
                if any(t1) && any(t2)
                    diff_tmp(i) = mean(t1_V1) - mean(t2_V1);
                end

                total_events = mean([sum((event_phase(2,:) >= -pi/2 & event_phase(2,:) <= pi/2) & (event_phase(1,:) >= -pi & event_phase(1,:) <= -pi/2 | event_phase(1,:) >= pi/2 & event_phase(1,:) <= pi)) ...
                    sum((event_phase(1,:) >= -pi/2 & event_phase(1,:) <= pi/2) & (event_phase(2,:) >= -pi & event_phase(2,:) <= -pi/2 | event_phase(2,:) >= pi/2 & event_phase(2,:) <= pi))]);

                prop_tmp(i) = (sum(t1) + sum(t2)) / total_events;



                t1s = bb_shift >= th;
                t2s = bb_shift <= -th;
                t1s = t1s & (spindle_phase_shifted(2,:) > -pi/2 & spindle_phase_shifted(2,:) < pi/2) & (event_phase_shifted(2,:) >= -pi & event_phase_shifted(2,:) < -pi/2 | event_phase_shifted(2,:) > pi/2 & event_phase_shifted(2,:) <= pi);
                t2s = t2s & (spindle_phase_shifted(1,:) > -pi/2 & spindle_phase_shifted(1,:) < pi/2) & (event_phase_shifted(1,:) >= -pi & event_phase_shifted(1,:) < -pi/2 | event_phase_shifted(1,:) > pi/2 & event_phase_shifted(1,:) <= pi);

                t1_V1 = boot_V1(t1s);
                t2_V1 = boot_V1(t2s);

                if any(t1s) && any(t2s)
                    diff_tmp_shifted(i) = mean(t1_V1) - mean(t2_V1);
                end
                prop_tmp_shifted(i) = (sum(t1s) + sum(t2s)) / total_events;

            elseif npower == 2 % if trough sync

                t1 = t1 & (spindle_phase(2,:) >= -pi & spindle_phase(2,:) <= -pi/2 | spindle_phase(2,:) >= pi/2 & spindle_phase(2,:) <= pi) & (event_phase(2,:) >= -pi & event_phase(2,:) <= -pi/2 | event_phase(2,:) >= pi/2 & event_phase(2,:) <= pi);
                t2 = t2 & (spindle_phase(1,:) >= -pi & spindle_phase(1,:) <= -pi/2 | spindle_phase(1,:) >= pi/2 & spindle_phase(1,:) <= pi) & (event_phase(1,:) >= -pi & event_phase(1,:) <= -pi/2 | event_phase(1,:) >= pi/2 & event_phase(1,:) <= pi);

                t1_V1 = boot_V1(t1);
                t2_V1 = boot_V1(t2);
                if any(t1) && any(t2)
                    diff_tmp(i) = mean(t1_V1) - mean(t2_V1);
                end

                total_events = sum((event_phase(2,:) >= -pi & event_phase(2,:) <= -pi/2 | event_phase(2,:) >= pi/2 & event_phase(2,:) <= pi) & (event_phase(1,:) >= -pi & event_phase(1,:) <= -pi/2 | event_phase(1,:) >= pi/2 & event_phase(1,:) <= pi));

                prop_tmp(i) = (sum(t1) + sum(t2)) / total_events;

                t1s = bb_shift >= th;
                t2s = bb_shift <= -th;
                t1s = t1s & (spindle_phase_shifted(2,:) >= -pi & spindle_phase_shifted(2,:) <= -pi/2 | spindle_phase_shifted(2,:) >= pi/2 & spindle_phase_shifted(2,:) <= pi) & (event_phase_shifted(2,:) >= -pi & event_phase_shifted(2,:) <= -pi/2 | event_phase_shifted(2,:) >= pi/2 & event_phase_shifted(2,:) <= pi);
                t2s = t2s & (spindle_phase_shifted(1,:) >= -pi & spindle_phase_shifted(1,:) <= -pi/2 | spindle_phase_shifted(1,:) >= pi/2 & spindle_phase_shifted(1,:) <= pi) & (event_phase_shifted(1,:) >= -pi & event_phase_shifted(1,:) <= -pi/2 | event_phase_shifted(1,:) >= pi/2 & event_phase_shifted(1,:) <= pi);


                t1_V1 = boot_V1(t1s);
                t2_V1 = boot_V1(t2s);

                if any(t1s) && any(t2s)
                    diff_tmp_shifted(i) = mean(t1_V1) - mean(t2_V1);
                end
                prop_tmp_shifted(i) = (sum(t1s) + sum(t2s)) / total_events;

            end



        end

        bias_diff_boot(iBoot, :) = diff_tmp;
        prop_events_boot(iBoot, :) = prop_tmp;
        bias_diff_shifted_boot(iBoot, :) = diff_tmp_shifted;
        prop_events_shifted_boot(iBoot, :) = prop_tmp_shifted;
    end

    % Stats
    bias_mean = mean(bias_diff_boot, 1, 'omitnan');
    bias_CI_lo = prctile(bias_diff_boot, 2.5, 1);
    bias_CI_hi = prctile(bias_diff_boot, 97.5, 1);
    prop_mean = mean(prop_events_boot, 1, 'omitnan');
    prop_CI_lo = prctile(prop_events_boot, 2.5, 1);
    prop_CI_hi = prctile(prop_events_boot, 97.5, 1);
    bias_shifted_mean = mean(bias_diff_shifted_boot, 1, 'omitnan');
    bias_shifted_CI_lo = prctile(bias_diff_shifted_boot, 2.5, 1);
    bias_shifted_CI_hi = prctile(bias_diff_shifted_boot, 97.5, 1);
    prop_shifted_mean = mean(prop_events_shifted_boot, 1, 'omitnan');
    prop_shifted_CI_lo = prctile(prop_events_shifted_boot, 2.5, 1);
    prop_shifted_CI_hi = prctile(prop_events_shifted_boot, 97.5, 1);

    % Store
    SO_phase_KDE_bias_difference(npower).phase_range{npower} = ...
        SO_thresholds{npower};
    SO_phase_KDE_bias_difference(npower).bias_diff_mean = bias_mean;
    SO_phase_KDE_bias_difference(npower).bias_diff_CI = [bias_CI_lo; bias_CI_hi];
    SO_phase_KDE_bias_difference(npower).prop_mean = prop_mean;
    SO_phase_KDE_bias_difference(npower).prop_CI = [prop_CI_lo; prop_CI_hi];
    SO_phase_KDE_bias_difference(npower).thresholds = thresholds;
    SO_phase_KDE_bias_difference(npower).bias_diff_shifted_mean = bias_shifted_mean;
    SO_phase_KDE_bias_difference(npower).bias_diff_shifted_CI = [bias_shifted_CI_lo; bias_shifted_CI_hi];
    SO_phase_KDE_bias_difference(npower).prop_shifted_mean = prop_shifted_mean;
    SO_phase_KDE_bias_difference(npower).prop_shifted_CI = [prop_shifted_CI_lo; prop_shifted_CI_hi];
    SO_phase_KDE_bias_difference(npower).total_events = total_events;

    % store AUC
    auc_boot = (trapz(thresholds, bias_diff_boot') / (max(thresholds)-min(thresholds)))';
    auc_shift_boot = (trapz(thresholds, bias_diff_shifted_boot') / (max(thresholds)-min(thresholds)))';

    SO_phase_KDE_bias_difference(npower).AUC_mean = mean(auc_boot, 'omitnan');
    SO_phase_KDE_bias_difference(npower).AUC_CI = prctile(auc_boot, [2.5 97.5]);
    SO_phase_KDE_bias_difference(npower).AUC_mean_shuffled = mean(auc_shift_boot, 'omitnan');
    SO_phase_KDE_bias_difference(npower).AUC_CI_shuffled = prctile(auc_shift_boot, [2.5 97.5]);

    % ----------- PLOTS (A, B, C) ----------------

    % A: Bias vs. Threshold
    nexttile((npower-1)*3 + 1);
    hold on;
    fill([thresholds, fliplr(thresholds)], [bias_CI_lo, fliplr(bias_CI_hi)], ...
        colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    plot(thresholds, bias_mean, 'Color', colour_lines(npower,:), 'LineWidth', 2);
    fill([thresholds, fliplr(thresholds)], ...
         [bias_shifted_CI_lo, fliplr(bias_shifted_CI_hi)], ...
         [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    plot(thresholds, bias_shifted_mean, 'k-', 'LineWidth', 1.5);
    ylim([-0.15 0.35]); xlim([0 1]); yline(0, '--r');
    xlabel('V1 bias threshold'); ylabel('HPC bias diff (T1 - T2)');
    title(sprintf('SO phase bin %s', SO_thresholds{npower}));
    set(gca,"TickDir","out",'box','off','Color','none','FontSize',12)

    % B: Proportion vs Bias Diff
    nexttile((npower-1)*3 + 2);
    hold on;
    valid = isfinite(bias_mean) & isfinite(prop_mean);
    fill([bias_CI_lo(valid), fliplr(bias_CI_hi(valid))], ...
         [prop_mean(valid), fliplr(prop_mean(valid))], ...
         colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    plot(bias_mean(valid), prop_mean(valid), '-', 'Color', colour_lines(npower,:), 'LineWidth', 2);
    fill([bias_shifted_CI_lo(valid), fliplr(bias_shifted_CI_hi(valid))], ...
         [prop_shifted_mean(valid), fliplr(prop_shifted_mean(valid))], ...
         [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    plot(bias_shifted_mean(valid), prop_shifted_mean(valid), 'k-', 'LineWidth', 1.5);
    xlim([-0.1 0.35]); xline(0, '--r');
    xlabel('HPC bias diff (T1 - T2)'); ylabel('Proportion of events detected');
    title('Event Proportion vs. Bias Difference');
    set(gca,"TickDir","out",'box','off','Color','none','FontSize',12)

    % C: Proportion CI vs Bias Diff
    nexttile((npower-1)*3 + 3);
    hold on;
    fill([bias_mean(valid), fliplr(bias_mean(valid))], ...
         [prop_CI_lo(valid), fliplr(prop_CI_hi(valid))], ...
         colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    plot(bias_mean(valid), prop_mean(valid), '-', 'Color', colour_lines(npower,:), 'LineWidth', 2);
    fill([bias_shifted_mean(valid), fliplr(bias_shifted_mean(valid))], ...
         [prop_shifted_CI_lo(valid), fliplr(prop_shifted_CI_hi(valid))], ...
         [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    plot(bias_shifted_mean(valid), prop_shifted_mean(valid), 'k-', 'LineWidth', 1.5);
    xlim([-0.1 0.35]); xline(0, '--r');
    xlabel('HPC bias diff (T1 - T2)'); ylabel('Proportion of events detected');
    title('Proportion vs. HPC Bias Difference');
    set(gca,"TickDir","out",'box','off','Color','none','FontSize',12)
    toc
end

clear Fill

fig = figure;
fig.Position = [640 100 2*1100/3 650/2];
% fig.Name = 'KDE bias difference in HPC low vs high SO power (PRE ripple)';
fig.Name = 'KDE bias difference in HPC SO + spindle trough synchrony (PRE ripple)';

% Select bins (1 = low, 4 = high)
nexttile;
for npower = [1 2]
    bias_mean = SO_phase_KDE_bias_difference(npower).bias_diff_mean;
    bias_CI_lo = SO_phase_KDE_bias_difference(npower).bias_diff_CI(1,:);
    bias_CI_hi = SO_phase_KDE_bias_difference(npower).bias_diff_CI(2,:);
    thresholds = SO_phase_KDE_bias_difference(npower).thresholds;

    hold on;
    x2 = [thresholds, fliplr(thresholds)];
    y2 = [bias_CI_lo, fliplr(bias_CI_hi)];
    Fill(npower) = fill(x2, y2, colour_lines(npower,:), ...
        'EdgeColor', 'none', 'FaceAlpha', 0.3);
    plot(thresholds, bias_mean, 'Color', colour_lines(npower,:), 'LineWidth', 2);

    xlabel('HPC Bias threshold');
    ylabel('V1 bias diff (T1 - T2)');
    set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12);
    ylim([-0.1 0.35]);
    % ylim([0 1])
end

bias_mean = SO_phase_KDE_bias_difference(2).bias_diff_shifted_mean;
bias_CI_lo = SO_phase_KDE_bias_difference(2).bias_diff_shifted_CI(1,:);
bias_CI_hi = SO_phase_KDE_bias_difference(2).bias_diff_shifted_CI(2,:)

hold on;
x2 = [thresholds, fliplr(thresholds)];
y2 = [bias_CI_lo, fliplr(bias_CI_hi)];
Fill(end +1) = fill(x2, y2, 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
plot(thresholds, bias_mean, 'Color','k', 'LineWidth', 2);

yline(0, '--r');
legend(Fill([1 2 3]), {'SO trough not sync', 'SO trough sync','Shuffled'}, 'box', 'off');

nexttile;
for npower = [1 2]
    bias_mean = SO_phase_KDE_bias_difference(npower).bias_diff_mean;
    bias_CI_lo = SO_phase_KDE_bias_difference(npower).bias_diff_CI(1,:);
    bias_CI_hi = SO_phase_KDE_bias_difference(npower).bias_diff_CI(2,:);
    prop_mean = SO_phase_KDE_bias_difference(npower).prop_mean;

    hold on;
    y2 = [prop_mean, fliplr(prop_mean)];
    x2 = [bias_CI_lo, fliplr(bias_CI_hi)];
    Fill(npower) = fill(x2, y2, colour_lines(npower,:), ...
        'EdgeColor', 'none', 'FaceAlpha', 0.3);
    plot(bias_mean, prop_mean, 'Color', colour_lines(npower,:), 'LineWidth', 2);

    xlabel('V1 bias diff (T1 - T2)');
    ylabel('Proportion of events detected');
    set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12);
    xlim([-0.1 0.35]);
    ylim([0 1])
end


bias_mean = SO_phase_KDE_bias_difference(2).bias_diff_shifted_mean;
bias_CI_lo = SO_phase_KDE_bias_difference(2).bias_diff_shifted_CI(1,:);
bias_CI_hi = SO_phase_KDE_bias_difference(2).bias_diff_shifted_CI(2,:)
prop_mean = SO_phase_KDE_bias_difference(2).prop_shifted_mean;

hold on;
y2 = [prop_mean, fliplr(prop_mean)];
x2 = [bias_CI_lo, fliplr(bias_CI_hi)];
Fill(end + 1) = fill(x2, y2, 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
plot(bias_mean, prop_mean, 'Color','k', 'LineWidth', 2);


xline(0, '--r');
legend(Fill([1 2 3]), {'SO trough not sync', 'SO trough sync','Shuffled'}, 'box', 'off');

% Save results
save(fullfile(analysis_folder,'V1-HPC sleep reactivation','SO_spindle_trough_KDE_bias_difference_based_on_PRE_V1_bias.mat'), 'SO_phase_KDE_bias_difference');
load(fullfile(analysis_folder,'V1-HPC sleep reactivation','SO_spindle_trough_KDE_bias_difference_based_on_PRE_V1_bias.mat'), 'SO_phase_KDE_bias_difference');



%%%%% AUC mean + CI bar plot
% Plot layout
fig = figure;
fig.Position = [640 100 281 325]
fig.Name = 'KDE bias V1 AUC SO + spindle trough synchrony (PRE)';
data = SO_phase_KDE_bias_difference;
n_bins = length(data);
bar_width = 0.3;      % Width of the bars
group_offset = 0.15;    % Distance from the center integer (half the gap between bars)
hold on;
clear BAR
for i = 1:4
    if i > length(data)
        x_shuf = i;
        y_shuf = data(1).AUC_mean_shuffled;
        bar(x_shuf, y_shuf, bar_width, ...
            'FaceColor', 'k', ...
            'FaceAlpha', 0.15, ...
            'EdgeColor', 'none');
    else


        % --- 1. Plot Shuffled Data (Left Bar) ---
        x_shuf = i - group_offset;
        y_shuf = data(i).AUC_mean_shuffled;

        % Calculate Error Deltas (Errorbar requires length relative to mean, not absolute values)
        % CI is [lower, upper]
        neg_err_shuf = y_shuf - data(i).AUC_CI_shuffled(1);
        pos_err_shuf = data(i).AUC_CI_shuffled(2) - y_shuf;

        % Plot Bar
        if i == 1
            BAR(1) = bar(x_shuf, y_shuf, bar_width, ...
                'FaceColor', 'k', ...
                'FaceAlpha', 0.15, ...
                'EdgeColor', 'none');
        else
            bar(x_shuf, y_shuf, bar_width, ...
                'FaceColor', 'k', ...
                'FaceAlpha', 0.15, ...
                'EdgeColor', 'none');
        end

        % Plot Error Bar
        E = errorbar(x_shuf, y_shuf, neg_err_shuf, pos_err_shuf, ...
            'Color', 'k', 'LineWidth', 1.5, 'CapSize', 8, 'LineStyle', 'none');


        % --- 2. Plot Real Data (Right Bar) ---
        x_real = i + group_offset;
        y_real = data(i).AUC_mean;

        neg_err_real = y_real - data(i).AUC_CI(1);
        pos_err_real = data(i).AUC_CI(2) - y_real;

        % Get specific color for this power bin
        % Assuming colour_lines is size [n_bins x 3]
        this_color = colour_lines(i, :);

        % Plot Bar
        BAR(i+1) = bar(x_real, y_real, bar_width, ...
            'FaceColor', this_color, ...
            'FaceAlpha', 0.3, ...
            'EdgeColor', 'none');

        errorbar(x_real, y_real, neg_err_real, pos_err_real, ...
            'Color', 'k', 'LineWidth', 1.5, 'CapSize', 8, 'LineStyle', 'none');
    end
end

hold off;
% Set X-ticks to be centered on the groups
set(gca, 'XTick', 1:nBins);
xlim([0.5, 4 + 0.5]);

% Labels
ylabel('V1 bias AUC');
xlabel('Power Bins');
legend([BAR(1:end)],{'Shuffled','SO trough spindle peak','SO trough spindle trough'},'box','off')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE bias difference based on V1 bias'),[])






%% log odds AUC with SO trough + different spindle power (V1→HPC)
% Spindle power binning across both probes
% all_spindle_power = mean(ripple_info.spindle_amplitude, 1);  % avg of probe 1 and 2
% SO_thresholds = prctile(ripple_info.SO_phase, 0:99.9/4:99.9);

SO_thresholds = {'SO trough 0-25 spindle','SO trough 25-50 spindle','SO trough 50-75 spindle','SO trough 75-100 spindle'};
spindle_thresholds = prctile(ripple_info.spindle_amplitude, 0:99.9/4:99.9);

% SO_thresholds = {'SO trough low spindle','SO trough high spindle'};
% spindle_thresholds = prctile(ripple_info.spindle_amplitude, 0:99.9/2:99.9);

% spindle_thresholds1 = prctile(ripple_info.spindle_amplitude(ripple_info.spindle_amplitude(:,1)>0,1), 0:99.9/2:99.9);
% spindle_thresholds2 = prctile(ripple_info.spindle_amplitude(ripple_info.spindle_amplitude(:,2)>0,2),  0:99.9/2:99.9);
% spindle_thresholds = [spindle_thresholds1;spindle_thresholds2]';
nBins =size(spindle_thresholds,1)-1;

% nBins = 2;

bins_to_use = bin_centers>0 & bin_centers<0.1;
bins_to_select = bin_centers>0 & bin_centers<0.15;
bins_to_use_shifted = bin_centers > -1 & bin_centers < -0.9;
nBoot = 1000;

SO_phase_KDE_bias_difference = struct;
% 
% colour_lines = [ ...
%     241, 182, 218;   % original end (lightest)
%     231, 41, 138    % original start (darkest)
% ] / 256;

colour_lines = [ ...
    241, 182, 218;
    226, 132, 187;
    212,  78, 156;
    231,  41, 138] / 256;

fig = figure;
fig.Position = [640 100 1100 2*650];
% fig.Position = [640 100 1100 650];
fig.Name = 'KDE bias difference in HPC with SO trough with different spindle power';
tiledlayout(nBins, 3, 'TileSpacing', 'compact');

for npower = 1:nBins
    tic
    % Select events in current spindle bin
    % power_index = ripple_info.spindle_amplitude(:,1) > spindle_thresholds(npower,1) & ...
    %               ripple_info.spindle_amplitude(:,1) <= spindle_thresholds(npower+1,1) |...
    %               ripple_info.spindle_amplitude(:,2) > spindle_thresholds(npower,2) & ...
    %               ripple_info.spindle_amplitude(:,2) <= spindle_thresholds(npower+1,2);

    event_index =1:length(z_bias);
    % event_index = find(power_index);

    mean_bias = mean(z_bias_V1(bins_to_select, event_index), 'omitnan');
    % mean_bias_shifted = mean(z_bias(bins_to_use_shifted, event_index), 'omitnan');
    mean_bias_V1 = mean(z_bias(bins_to_use, event_index), 'omitnan');
    selected_events = length(mean_bias);

    thresholds = prctile(abs(mean_bias), 0:10:100);
    thresholds = thresholds(1:end-1);
    nThresh = length(thresholds);

    bias_diff_boot = NaN(nBoot, nThresh);
    prop_events_boot = NaN(nBoot, nThresh);
    bias_diff_shifted_boot = NaN(nBoot, nThresh);
    prop_events_shifted_boot = NaN(nBoot, nThresh);


    % event_phase = ripple_info.SO_phase';


    parfor iBoot = 1:nBoot
        s = RandStream('philox4x32_10', 'Seed', iBoot);
        idx = randi(s, selected_events, selected_events, 1);

        true_idx = find(event_index);

        bb = mean_bias(idx);
        bb_shift = mean_bias;
        boot_V1 = mean_bias_V1(idx);

        diff_tmp = NaN(1, nThresh);
        prop_tmp = NaN(1, nThresh);
        diff_tmp_shifted = NaN(1, nThresh);
        prop_tmp_shifted = NaN(1, nThresh);
        event_phase = ripple_info.SO_phase(true_idx(idx),:)';
        event_phase_shifted = ripple_info.SO_phase(true_idx,:)';

        for i = 1:nThresh
            th = thresholds(i);

            t1 = bb >= th;
            t2 = bb <= -th;

            % SO trough with different spindle powers
            t1 = t1 & (ripple_info.spindle_amplitude(true_idx(idx),2) > spindle_thresholds(npower,2) & ...
                ripple_info.spindle_amplitude(true_idx(idx),2) <= spindle_thresholds(npower+1,2))' & ...
                (event_phase(2,:) >= -pi & event_phase(2,:) < -pi/2 | event_phase(2,:) > pi/2 & event_phase(2,:) <= pi);

            t2 = t2 & (ripple_info.spindle_amplitude(true_idx(idx),1) > spindle_thresholds(npower,1) & ...
                ripple_info.spindle_amplitude(true_idx(idx),1) <= spindle_thresholds(npower+1,1))' & ...
                (event_phase(1,:) >= -pi & event_phase(1,:) < -pi/2 | event_phase(1,:) > pi/2 & event_phase(1,:) <= pi);

            t1_V1 = boot_V1(t1);
            t2_V1 = boot_V1(t2);
            if any(t1) && any(t2)
                diff_tmp(i) = mean(t1_V1) - mean(t2_V1);
            end

            total_events = round(mean([sum((ripple_info.spindle_amplitude(:,1) > spindle_thresholds(npower,1) & ...
                ripple_info.spindle_amplitude(:,1) <= spindle_thresholds(npower+1,1))' & ...
                (event_phase(1,:) >= -pi & event_phase(1,:) < -pi/2 | event_phase(1,:) > pi/2 & event_phase(1,:) <= pi))
                sum((ripple_info.spindle_amplitude(:,2) > spindle_thresholds(npower,2) & ...
                ripple_info.spindle_amplitude(:,2) <= spindle_thresholds(npower+1,2))' & ...
                (event_phase(2,:) >= -pi & event_phase(2,:) < -pi/2 | event_phase(2,:) > pi/2 & event_phase(2,:) <= pi))]));

            prop_tmp(i) = (sum(t1) + sum(t2)) / total_events;

            t1s = bb_shift >= th;
            t2s = bb_shift <= -th;

            t1s = t1s & (ripple_info.spindle_amplitude(true_idx,2) > spindle_thresholds(npower,2) & ...
                ripple_info.spindle_amplitude(true_idx,2) <= spindle_thresholds(npower+1,2))' & ...
                (event_phase_shifted(2,:) >= -pi & event_phase_shifted(2,:) < -pi/2 | event_phase_shifted(2,:) > pi/2 & event_phase_shifted(2,:) <= pi);

            t2s = t2s & (ripple_info.spindle_amplitude(true_idx,1) > spindle_thresholds(npower,1) & ...
                ripple_info.spindle_amplitude(true_idx,1) <= spindle_thresholds(npower+1,1))' & ...
                (event_phase_shifted(1,:) >= -pi & event_phase_shifted(1,:) < -pi/2 | event_phase_shifted(1,:) > pi/2 & event_phase_shifted(1,:) <= pi);
            t1_V1 = boot_V1(t1s);
            t2_V1 = boot_V1(t2s);

            if any(t1s) && any(t2s)
                diff_tmp_shifted(i) = mean(t1_V1) - mean(t2_V1);
            end
            prop_tmp_shifted(i) = (sum(t1s) + sum(t2s)) / total_events;

        end

        bias_diff_boot(iBoot, :) = diff_tmp;
        prop_events_boot(iBoot, :) = prop_tmp;
        bias_diff_shifted_boot(iBoot, :) = diff_tmp_shifted;
        prop_events_shifted_boot(iBoot, :) = prop_tmp_shifted;
    end

    % Stats
    bias_mean = mean(bias_diff_boot, 1, 'omitnan');
    bias_CI_lo = prctile(bias_diff_boot, 2.5, 1);
    bias_CI_hi = prctile(bias_diff_boot, 97.5, 1);
    prop_mean = mean(prop_events_boot, 1, 'omitnan');
    prop_CI_lo = prctile(prop_events_boot, 2.5, 1);
    prop_CI_hi = prctile(prop_events_boot, 97.5, 1);
    bias_shifted_mean = mean(bias_diff_shifted_boot, 1, 'omitnan');
    bias_shifted_CI_lo = prctile(bias_diff_shifted_boot, 2.5, 1);
    bias_shifted_CI_hi = prctile(bias_diff_shifted_boot, 97.5, 1);
    prop_shifted_mean = mean(prop_events_shifted_boot, 1, 'omitnan');
    prop_shifted_CI_lo = prctile(prop_events_shifted_boot, 2.5, 1);
    prop_shifted_CI_hi = prctile(prop_events_shifted_boot, 97.5, 1);

    % Store
    SO_phase_KDE_bias_difference(npower).phase_range{npower} = ...
        SO_thresholds{npower};
    SO_phase_KDE_bias_difference(npower).bias_diff_mean = bias_mean;
    SO_phase_KDE_bias_difference(npower).bias_diff_CI = [bias_CI_lo; bias_CI_hi];
    SO_phase_KDE_bias_difference(npower).prop_mean = prop_mean;
    SO_phase_KDE_bias_difference(npower).prop_CI = [prop_CI_lo; prop_CI_hi];
    SO_phase_KDE_bias_difference(npower).thresholds = thresholds;
    SO_phase_KDE_bias_difference(npower).bias_diff_shifted_mean = bias_shifted_mean;
    SO_phase_KDE_bias_difference(npower).bias_diff_shifted_CI = [bias_shifted_CI_lo; bias_shifted_CI_hi];
    SO_phase_KDE_bias_difference(npower).prop_shifted_mean = prop_shifted_mean;
    SO_phase_KDE_bias_difference(npower).prop_shifted_CI = [prop_shifted_CI_lo; prop_shifted_CI_hi];
    SO_phase_KDE_bias_difference(npower).total_events = total_events;

    % store AUC
    auc_boot = (trapz(thresholds, bias_diff_boot') / (max(thresholds)-min(thresholds)))';
    auc_shift_boot = (trapz(thresholds, bias_diff_shifted_boot') / (max(thresholds)-min(thresholds)))';

    SO_phase_KDE_bias_difference(npower).AUC_mean = mean(auc_boot, 'omitnan');
    SO_phase_KDE_bias_difference(npower).AUC_CI = prctile(auc_boot, [2.5 97.5]);
    SO_phase_KDE_bias_difference(npower).AUC_mean_shuffled = mean(auc_shift_boot, 'omitnan');
    SO_phase_KDE_bias_difference(npower).AUC_CI_shuffled = prctile(auc_shift_boot, [2.5 97.5]);

    % ----------- PLOTS (A, B, C) ----------------

    % A: Bias vs. Threshold
    nexttile((npower-1)*3 + 1);
    hold on;
    fill([thresholds, fliplr(thresholds)], [bias_CI_lo, fliplr(bias_CI_hi)], ...
        colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    plot(thresholds, bias_mean, 'Color', colour_lines(npower,:), 'LineWidth', 2);
    fill([thresholds, fliplr(thresholds)], ...
         [bias_shifted_CI_lo, fliplr(bias_shifted_CI_hi)], ...
         [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    plot(thresholds, bias_shifted_mean, 'k-', 'LineWidth', 1.5);
    ylim([-0.15 0.35]); xlim([0 1]); yline(0, '--r');
    xlabel('V1 bias threshold'); ylabel('HPC bias diff (T1 - T2)');
    title(sprintf('SO phase bin %s', SO_thresholds{npower}));
    set(gca,"TickDir","out",'box','off','Color','none','FontSize',12)

    % B: Proportion vs Bias Diff
    nexttile((npower-1)*3 + 2);
    hold on;
    valid = isfinite(bias_mean) & isfinite(prop_mean);
    fill([bias_CI_lo(valid), fliplr(bias_CI_hi(valid))], ...
         [prop_mean(valid), fliplr(prop_mean(valid))], ...
         colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    plot(bias_mean(valid), prop_mean(valid), '-', 'Color', colour_lines(npower,:), 'LineWidth', 2);
    fill([bias_shifted_CI_lo(valid), fliplr(bias_shifted_CI_hi(valid))], ...
         [prop_shifted_mean(valid), fliplr(prop_shifted_mean(valid))], ...
         [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    plot(bias_shifted_mean(valid), prop_shifted_mean(valid), 'k-', 'LineWidth', 1.5);
    xlim([-0.1 0.35]); xline(0, '--r');
    xlabel('HPC bias diff (T1 - T2)'); ylabel('Proportion of events detected');
    title('Event Proportion vs. Bias Difference');
    set(gca,"TickDir","out",'box','off','Color','none','FontSize',12)

    % C: Proportion CI vs Bias Diff
    nexttile((npower-1)*3 + 3);
    hold on;
    fill([bias_mean(valid), fliplr(bias_mean(valid))], ...
         [prop_CI_lo(valid), fliplr(prop_CI_hi(valid))], ...
         colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    plot(bias_mean(valid), prop_mean(valid), '-', 'Color', colour_lines(npower,:), 'LineWidth', 2);
    fill([bias_shifted_mean(valid), fliplr(bias_shifted_mean(valid))], ...
         [prop_shifted_CI_lo(valid), fliplr(prop_shifted_CI_hi(valid))], ...
         [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    plot(bias_shifted_mean(valid), prop_shifted_mean(valid), 'k-', 'LineWidth', 1.5);
    xlim([-0.1 0.35]); xline(0, '--r');
    xlabel('HPC bias diff (T1 - T2)'); ylabel('Proportion of events detected');
    title('Proportion vs. HPC Bias Difference');
    set(gca,"TickDir","out",'box','off','Color','none','FontSize',12)
    toc
end

clear Fill

fig = figure;
% fig.Position = [640 100 2*1100/3 650/2];
fig.Position = [640 100 2*1100/3 650/2];
% fig.Name = 'KDE bias difference in HPC low vs high SO power (PRE ripple)';
fig.Name = 'KDE bias difference in HPC SO trough with different spindle power';

% Select bins (1 = low, 4 = high)
nexttile;
for npower = [1 4]
    bias_mean = SO_phase_KDE_bias_difference(npower).bias_diff_mean;
    bias_CI_lo = SO_phase_KDE_bias_difference(npower).bias_diff_CI(1,:);
    bias_CI_hi = SO_phase_KDE_bias_difference(npower).bias_diff_CI(2,:);
    thresholds = SO_phase_KDE_bias_difference(npower).thresholds;

    hold on;
    x2 = [thresholds, fliplr(thresholds)];
    y2 = [bias_CI_lo, fliplr(bias_CI_hi)];
    Fill(npower) = fill(x2, y2, colour_lines(npower,:), ...
        'EdgeColor', 'none', 'FaceAlpha', 0.3);
    plot(thresholds, bias_mean, 'Color', colour_lines(npower,:), 'LineWidth', 2);

    xlabel('HPC Bias threshold');
    ylabel('V1 bias diff (T1 - T2)');
    set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12);
    ylim([-0.1 0.35]);
    % ylim([0 1])
end

bias_mean = SO_phase_KDE_bias_difference(npower).bias_diff_shifted_mean;
bias_CI_lo = SO_phase_KDE_bias_difference(npower).bias_diff_shifted_CI(1,:);
bias_CI_hi = SO_phase_KDE_bias_difference(npower).bias_diff_shifted_CI(2,:)

hold on;
x2 = [thresholds, fliplr(thresholds)];
y2 = [bias_CI_lo, fliplr(bias_CI_hi)];
Fill(end +1) = fill(x2, y2, 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
plot(thresholds, bias_mean, 'Color','k', 'LineWidth', 2);

yline(0, '--r');
legend(Fill([1 4 5]), {'SO trough low spindle', 'SO trough high spindle','Shuffled'}, 'box', 'off');
% legend(Fill([1 2 3]), {'SO trough low spindle', 'SO trough high spindle','Shuffled'}, 'box', 'off');

nexttile;
for npower = [1 4]
    bias_mean = SO_phase_KDE_bias_difference(npower).bias_diff_mean;
    bias_CI_lo = SO_phase_KDE_bias_difference(npower).bias_diff_CI(1,:);
    bias_CI_hi = SO_phase_KDE_bias_difference(npower).bias_diff_CI(2,:);
    prop_mean = SO_phase_KDE_bias_difference(npower).prop_mean;

    hold on;
    y2 = [prop_mean, fliplr(prop_mean)];
    x2 = [bias_CI_lo, fliplr(bias_CI_hi)];
    Fill(npower) = fill(x2, y2, colour_lines(npower,:), ...
        'EdgeColor', 'none', 'FaceAlpha', 0.3);
    plot(bias_mean, prop_mean, 'Color', colour_lines(npower,:), 'LineWidth', 2);

    xlabel('V1 bias diff (T1 - T2)');
    ylabel('Proportion of events detected');
    set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12);
    xlim([-0.1 0.35]);
    ylim([0 1])
end


bias_mean = SO_phase_KDE_bias_difference(2).bias_diff_shifted_mean;
bias_CI_lo = SO_phase_KDE_bias_difference(2).bias_diff_shifted_CI(1,:);
bias_CI_hi = SO_phase_KDE_bias_difference(2).bias_diff_shifted_CI(2,:)
prop_mean = SO_phase_KDE_bias_difference(2).prop_shifted_mean;

hold on;
y2 = [prop_mean, fliplr(prop_mean)];
x2 = [bias_CI_lo, fliplr(bias_CI_hi)];
Fill(end + 1) = fill(x2, y2, 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
plot(bias_mean, prop_mean, 'Color','k', 'LineWidth', 2);

xline(0, '--r');
% legend(Fill([1 4 5]), {'SO trough low spindle', 'SO trough high spindle','Shuffled'}, 'box', 'off');


% Save results
save(fullfile(analysis_folder,'V1-HPC sleep reactivation','SO_trough_spindle_power_KDE_bias_difference_based_on_V1_bias.mat'), 'SO_phase_KDE_bias_difference');
load(fullfile(analysis_folder,'V1-HPC sleep reactivation','SO_trough_spindle_power_KDE_bias_difference_based_on_V1_bias.mat'), 'SO_phase_KDE_bias_difference');




%%%%% AUC mean + CI bar plot
% Plot layout
fig = figure;
fig.Position = [640 100 281 325]
fig.Name = 'KDE bias V1 AUC SO trough with different spindle power';
data = SO_phase_KDE_bias_difference;
n_bins = length(data);
bar_width = 0.3;      % Width of the bars
group_offset = 0.15;    % Distance from the center integer (half the gap between bars)
hold on;
clear BAR
for i = 1:4
    if i > length(data)
        x_shuf = i;
        y_shuf = data(1).AUC_mean_shuffled;
        bar(x_shuf, y_shuf, bar_width, ...
            'FaceColor', 'k', ...
            'FaceAlpha', 0.15, ...
            'EdgeColor', 'none');
    else


        % --- 1. Plot Shuffled Data (Left Bar) ---
        x_shuf = i - group_offset;
        y_shuf = data(i).AUC_mean_shuffled;

        % Calculate Error Deltas (Errorbar requires length relative to mean, not absolute values)
        % CI is [lower, upper]
        neg_err_shuf = y_shuf - data(i).AUC_CI_shuffled(1);
        pos_err_shuf = data(i).AUC_CI_shuffled(2) - y_shuf;

        % Plot Bar
        if i == 1
            BAR(1) = bar(x_shuf, y_shuf, bar_width, ...
                'FaceColor', 'k', ...
                'FaceAlpha', 0.15, ...
                'EdgeColor', 'none');
        else
            bar(x_shuf, y_shuf, bar_width, ...
                'FaceColor', 'k', ...
                'FaceAlpha', 0.15, ...
                'EdgeColor', 'none');
        end

        % Plot Error Bar
        E = errorbar(x_shuf, y_shuf, neg_err_shuf, pos_err_shuf, ...
            'Color', 'k', 'LineWidth', 1.5, 'CapSize', 8, 'LineStyle', 'none');


        % --- 2. Plot Real Data (Right Bar) ---
        x_real = i + group_offset;
        y_real = data(i).AUC_mean;

        neg_err_real = y_real - data(i).AUC_CI(1);
        pos_err_real = data(i).AUC_CI(2) - y_real;

        % Get specific color for this power bin
        % Assuming colour_lines is size [n_bins x 3]
        this_color = colour_lines(i, :);

        % Plot Bar
        BAR(i+1) = bar(x_real, y_real, bar_width, ...
            'FaceColor', this_color, ...
            'FaceAlpha', 0.3, ...
            'EdgeColor', 'none');

        errorbar(x_real, y_real, neg_err_real, pos_err_real, ...
            'Color', 'k', 'LineWidth', 1.5, 'CapSize', 8, 'LineStyle', 'none');
    end
end

hold off;
% Set X-ticks to be centered on the groups
set(gca, 'XTick', 1:nBins);
xlim([0.5, 4 + 0.5]);

% Labels
ylabel('V1 bias AUC');
xlabel('Power Bins');
legend([BAR(1:end)],{'Shuffled','SO trough low spindle','SO trough high spindle'},'box','off')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE bias difference based on V1 bias'),[])


%%


%% log odds AUC with SO trough + different spindle power (V1→HPC)
% Spindle power binning across both probes
% all_spindle_power = mean(ripple_info.spindle_amplitude, 1);  % avg of probe 1 and 2
% SO_thresholds = prctile(ripple_info.SO_phase, 0:99.9/4:99.9);

SO_thresholds = {'SO trough 0-25 spindle','SO trough 25-50 spindle','SO trough 50-75 spindle','SO trough 75-100 spindle'};
spindle_thresholds = prctile(ripple_info.spindle_amplitude, 0:99.9/4:99.9);

% SO_thresholds = {'SO trough low spindle','SO trough high spindle'};
% spindle_thresholds = prctile(ripple_info.spindle_amplitude, 0:99.9/2:99.9);

% spindle_thresholds1 = prctile(ripple_info.spindle_amplitude(ripple_info.spindle_amplitude(:,1)>0,1), 0:99.9/2:99.9);
% spindle_thresholds2 = prctile(ripple_info.spindle_amplitude(ripple_info.spindle_amplitude(:,2)>0,2),  0:99.9/2:99.9);
% spindle_thresholds = [spindle_thresholds1;spindle_thresholds2]';
nBins =size(spindle_thresholds,1)-1;

% nBins = 2;

bins_to_use = bin_centers>0 & bin_centers<0.1;
bins_to_select = bin_centers>-0.15 & bin_centers<0;
bins_to_use_shifted = bin_centers > -1 & bin_centers < -0.9;
nBoot = 1000;

SO_phase_KDE_bias_difference = struct;
% 
% colour_lines = [ ...
%     241, 182, 218;   % original end (lightest)
%     231, 41, 138    % original start (darkest)
% ] / 256;

colour_lines = [ ...
    241, 182, 218;
    226, 132, 187;
    212,  78, 156;
    231,  41, 138] / 256;

fig = figure;
fig.Position = [640 100 1100 2*650];
% fig.Position = [640 100 1100 650];
fig.Name = 'KDE bias difference in HPC with SO trough with different spindle power (PRE)';
tiledlayout(nBins, 3, 'TileSpacing', 'compact');

for npower = 1:nBins
    tic
    % Select events in current spindle bin
    % power_index = ripple_info.spindle_amplitude(:,1) > spindle_thresholds(npower,1) & ...
    %               ripple_info.spindle_amplitude(:,1) <= spindle_thresholds(npower+1,1) |...
    %               ripple_info.spindle_amplitude(:,2) > spindle_thresholds(npower,2) & ...
    %               ripple_info.spindle_amplitude(:,2) <= spindle_thresholds(npower+1,2);

    event_index =1:length(z_bias);
    % event_index = find(power_index);

    mean_bias = mean(z_bias_V1(bins_to_select, event_index), 'omitnan');
    % mean_bias_shifted = mean(z_bias(bins_to_use_shifted, event_index), 'omitnan');
    mean_bias_V1 = mean(z_bias(bins_to_use, event_index), 'omitnan');
    selected_events = length(mean_bias);

    thresholds = prctile(abs(mean_bias), 0:10:100);
    thresholds = thresholds(1:end-1);
    nThresh = length(thresholds);

    bias_diff_boot = NaN(nBoot, nThresh);
    prop_events_boot = NaN(nBoot, nThresh);
    bias_diff_shifted_boot = NaN(nBoot, nThresh);
    prop_events_shifted_boot = NaN(nBoot, nThresh);


    % event_phase = ripple_info.SO_phase';


    parfor iBoot = 1:nBoot
        s = RandStream('philox4x32_10', 'Seed', iBoot);
        idx = randi(s, selected_events, selected_events, 1);

        true_idx = find(event_index);

        bb = mean_bias(idx);
        bb_shift = mean_bias;
        boot_V1 = mean_bias_V1(idx);

        diff_tmp = NaN(1, nThresh);
        prop_tmp = NaN(1, nThresh);
        diff_tmp_shifted = NaN(1, nThresh);
        prop_tmp_shifted = NaN(1, nThresh);
        event_phase = ripple_info.SO_phase(true_idx(idx),:)';
        event_phase_shifted = ripple_info.SO_phase(true_idx,:)';

        for i = 1:nThresh
            th = thresholds(i);

            t1 = bb >= th;
            t2 = bb <= -th;
% boot_HPC = [];
            % SO trough with different spindle powers
            t1 = t1 & (ripple_info.spindle_amplitude(true_idx(idx),2) > spindle_thresholds(npower,2) & ...
                ripple_info.spindle_amplitude(true_idx(idx),2) <= spindle_thresholds(npower+1,2))' & ...
                (event_phase(2,:) >= -pi & event_phase(2,:) < -pi/2 | event_phase(2,:) > pi/2 & event_phase(2,:) <= pi);

            t2 = t2 & (ripple_info.spindle_amplitude(true_idx(idx),1) > spindle_thresholds(npower,1) & ...
                ripple_info.spindle_amplitude(true_idx(idx),1) <= spindle_thresholds(npower+1,1))' & ...
                (event_phase(1,:) >= -pi & event_phase(1,:) < -pi/2 | event_phase(1,:) > pi/2 & event_phase(1,:) <= pi);

            t1_V1 = boot_V1(t1);
            t2_V1 = boot_V1(t2);
            if any(t1) && any(t2)
                diff_tmp(i) = mean(t1_V1) - mean(t2_V1);
            end

            total_events = round(mean([sum((ripple_info.spindle_amplitude(:,1) > spindle_thresholds(npower,1) & ...
                ripple_info.spindle_amplitude(:,1) <= spindle_thresholds(npower+1,1))' & ...
                (event_phase(1,:) >= -pi & event_phase(1,:) < -pi/2 | event_phase(1,:) > pi/2 & event_phase(1,:) <= pi))
                sum((ripple_info.spindle_amplitude(:,2) > spindle_thresholds(npower,2) & ...
                ripple_info.spindle_amplitude(:,2) <= spindle_thresholds(npower+1,2))' & ...
                (event_phase(2,:) >= -pi & event_phase(2,:) < -pi/2 | event_phase(2,:) > pi/2 & event_phase(2,:) <= pi))]));

            prop_tmp(i) = (sum(t1) + sum(t2)) / total_events;

            t1s = bb_shift >= th;
            t2s = bb_shift <= -th;

            t1s = t1s & (ripple_info.spindle_amplitude(true_idx,2) > spindle_thresholds(npower,2) & ...
                ripple_info.spindle_amplitude(true_idx,2) <= spindle_thresholds(npower+1,2))' & ...
                (event_phase_shifted(2,:) >= -pi & event_phase_shifted(2,:) < -pi/2 | event_phase_shifted(2,:) > pi/2 & event_phase_shifted(2,:) <= pi);

            t2s = t2s & (ripple_info.spindle_amplitude(true_idx,1) > spindle_thresholds(npower,1) & ...
                ripple_info.spindle_amplitude(true_idx,1) <= spindle_thresholds(npower+1,1))' & ...
                (event_phase_shifted(1,:) >= -pi & event_phase_shifted(1,:) < -pi/2 | event_phase_shifted(1,:) > pi/2 & event_phase_shifted(1,:) <= pi);
            t1_V1 = boot_V1(t1s);
            t2_V1 = boot_V1(t2s);

            if any(t1s) && any(t2s)
                diff_tmp_shifted(i) = mean(t1_V1) - mean(t2_V1);
            end
            prop_tmp_shifted(i) = (sum(t1s) + sum(t2s)) / total_events;

        end

        bias_diff_boot(iBoot, :) = diff_tmp;
        prop_events_boot(iBoot, :) = prop_tmp;
        bias_diff_shifted_boot(iBoot, :) = diff_tmp_shifted;
        prop_events_shifted_boot(iBoot, :) = prop_tmp_shifted;
    end

    % Stats
    bias_mean = mean(bias_diff_boot, 1, 'omitnan');
    bias_CI_lo = prctile(bias_diff_boot, 2.5, 1);
    bias_CI_hi = prctile(bias_diff_boot, 97.5, 1);
    prop_mean = mean(prop_events_boot, 1, 'omitnan');
    prop_CI_lo = prctile(prop_events_boot, 2.5, 1);
    prop_CI_hi = prctile(prop_events_boot, 97.5, 1);
    bias_shifted_mean = mean(bias_diff_shifted_boot, 1, 'omitnan');
    bias_shifted_CI_lo = prctile(bias_diff_shifted_boot, 2.5, 1);
    bias_shifted_CI_hi = prctile(bias_diff_shifted_boot, 97.5, 1);
    prop_shifted_mean = mean(prop_events_shifted_boot, 1, 'omitnan');
    prop_shifted_CI_lo = prctile(prop_events_shifted_boot, 2.5, 1);
    prop_shifted_CI_hi = prctile(prop_events_shifted_boot, 97.5, 1);

    % Store
    SO_phase_KDE_bias_difference(npower).phase_range{npower} = ...
        SO_thresholds{npower};
    SO_phase_KDE_bias_difference(npower).bias_diff_mean = bias_mean;
    SO_phase_KDE_bias_difference(npower).bias_diff_CI = [bias_CI_lo; bias_CI_hi];
    SO_phase_KDE_bias_difference(npower).prop_mean = prop_mean;
    SO_phase_KDE_bias_difference(npower).prop_CI = [prop_CI_lo; prop_CI_hi];
    SO_phase_KDE_bias_difference(npower).thresholds = thresholds;
    SO_phase_KDE_bias_difference(npower).bias_diff_shifted_mean = bias_shifted_mean;
    SO_phase_KDE_bias_difference(npower).bias_diff_shifted_CI = [bias_shifted_CI_lo; bias_shifted_CI_hi];
    SO_phase_KDE_bias_difference(npower).prop_shifted_mean = prop_shifted_mean;
    SO_phase_KDE_bias_difference(npower).prop_shifted_CI = [prop_shifted_CI_lo; prop_shifted_CI_hi];
    SO_phase_KDE_bias_difference(npower).total_events = total_events;

    % store AUC
    auc_boot = (trapz(thresholds, bias_diff_boot') / (max(thresholds)-min(thresholds)))';
    auc_shift_boot = (trapz(thresholds, bias_diff_shifted_boot') / (max(thresholds)-min(thresholds)))';

    SO_phase_KDE_bias_difference(npower).AUC_mean = mean(auc_boot, 'omitnan');
    SO_phase_KDE_bias_difference(npower).AUC_CI = prctile(auc_boot, [2.5 97.5]);
    SO_phase_KDE_bias_difference(npower).AUC_mean_shuffled = mean(auc_shift_boot, 'omitnan');
    SO_phase_KDE_bias_difference(npower).AUC_CI_shuffled = prctile(auc_shift_boot, [2.5 97.5]);

    % ----------- PLOTS (A, B, C) ----------------

    % A: Bias vs. Threshold
    nexttile((npower-1)*3 + 1);
    hold on;
    fill([thresholds, fliplr(thresholds)], [bias_CI_lo, fliplr(bias_CI_hi)], ...
        colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    plot(thresholds, bias_mean, 'Color', colour_lines(npower,:), 'LineWidth', 2);
    fill([thresholds, fliplr(thresholds)], ...
         [bias_shifted_CI_lo, fliplr(bias_shifted_CI_hi)], ...
         [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    plot(thresholds, bias_shifted_mean, 'k-', 'LineWidth', 1.5);
    ylim([-0.15 0.35]); xlim([0 1]); yline(0, '--r');
    xlabel('V1 bias threshold'); ylabel('HPC bias diff (T1 - T2)');
    title(sprintf('SO phase bin %s', SO_thresholds{npower}));
    set(gca,"TickDir","out",'box','off','Color','none','FontSize',12)

    % B: Proportion vs Bias Diff
    nexttile((npower-1)*3 + 2);
    hold on;
    valid = isfinite(bias_mean) & isfinite(prop_mean);
    fill([bias_CI_lo(valid), fliplr(bias_CI_hi(valid))], ...
         [prop_mean(valid), fliplr(prop_mean(valid))], ...
         colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    plot(bias_mean(valid), prop_mean(valid), '-', 'Color', colour_lines(npower,:), 'LineWidth', 2);
    fill([bias_shifted_CI_lo(valid), fliplr(bias_shifted_CI_hi(valid))], ...
         [prop_shifted_mean(valid), fliplr(prop_shifted_mean(valid))], ...
         [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    plot(bias_shifted_mean(valid), prop_shifted_mean(valid), 'k-', 'LineWidth', 1.5);
    xlim([-0.1 0.35]); xline(0, '--r');
    xlabel('HPC bias diff (T1 - T2)'); ylabel('Proportion of events detected');
    title('Event Proportion vs. Bias Difference');
    set(gca,"TickDir","out",'box','off','Color','none','FontSize',12)

    % C: Proportion CI vs Bias Diff
    nexttile((npower-1)*3 + 3);
    hold on;
    fill([bias_mean(valid), fliplr(bias_mean(valid))], ...
         [prop_CI_lo(valid), fliplr(prop_CI_hi(valid))], ...
         colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    plot(bias_mean(valid), prop_mean(valid), '-', 'Color', colour_lines(npower,:), 'LineWidth', 2);
    fill([bias_shifted_mean(valid), fliplr(bias_shifted_mean(valid))], ...
         [prop_shifted_CI_lo(valid), fliplr(prop_shifted_CI_hi(valid))], ...
         [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    plot(bias_shifted_mean(valid), prop_shifted_mean(valid), 'k-', 'LineWidth', 1.5);
    xlim([-0.1 0.35]); xline(0, '--r');
    xlabel('HPC bias diff (T1 - T2)'); ylabel('Proportion of events detected');
    title('Proportion vs. HPC Bias Difference');
    set(gca,"TickDir","out",'box','off','Color','none','FontSize',12)
    toc
end

clear Fill

fig = figure;
fig.Position = [640 100 2*1100/3 650/2];
% fig.Name = 'KDE bias difference in HPC low vs high SO power (PRE ripple)';
fig.Name = 'KDE bias difference in HPC SO trough with different spindle power (PRE ripple)';

% Select bins (1 = low, 4 = high)
nexttile;
for npower = [1 4]
    bias_mean = SO_phase_KDE_bias_difference(npower).bias_diff_mean;
    bias_CI_lo = SO_phase_KDE_bias_difference(npower).bias_diff_CI(1,:);
    bias_CI_hi = SO_phase_KDE_bias_difference(npower).bias_diff_CI(2,:);
    thresholds = SO_phase_KDE_bias_difference(npower).thresholds;

    hold on;
    x2 = [thresholds, fliplr(thresholds)];
    y2 = [bias_CI_lo, fliplr(bias_CI_hi)];
    Fill(npower) = fill(x2, y2, colour_lines(npower,:), ...
        'EdgeColor', 'none', 'FaceAlpha', 0.3);
    plot(thresholds, bias_mean, 'Color', colour_lines(npower,:), 'LineWidth', 2);

    xlabel('HPC Bias threshold');
    ylabel('V1 bias diff (T1 - T2)');
    set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12);
    ylim([-0.1 0.35]);
    % ylim([0 1])
end

bias_mean = SO_phase_KDE_bias_difference(npower).bias_diff_shifted_mean;
bias_CI_lo = SO_phase_KDE_bias_difference(npower).bias_diff_shifted_CI(1,:);
bias_CI_hi = SO_phase_KDE_bias_difference(npower).bias_diff_shifted_CI(2,:)

hold on;
x2 = [thresholds, fliplr(thresholds)];
y2 = [bias_CI_lo, fliplr(bias_CI_hi)];
Fill(end +1) = fill(x2, y2, 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
plot(thresholds, bias_mean, 'Color','k', 'LineWidth', 2);

yline(0, '--r');
legend(Fill([1 4 5]), {'SO trough low spindle', 'SO trough high spindle','Shuffled'}, 'box', 'off');
% legend(Fill([1 2 3]), {'SO trough low spindle', 'SO trough high spindle','Shuffled'}, 'box', 'off');

nexttile;
for npower = [1 4]
    bias_mean = SO_phase_KDE_bias_difference(npower).bias_diff_mean;
    bias_CI_lo = SO_phase_KDE_bias_difference(npower).bias_diff_CI(1,:);
    bias_CI_hi = SO_phase_KDE_bias_difference(npower).bias_diff_CI(2,:);
    prop_mean = SO_phase_KDE_bias_difference(npower).prop_mean;

    hold on;
    y2 = [prop_mean, fliplr(prop_mean)];
    x2 = [bias_CI_lo, fliplr(bias_CI_hi)];
    Fill(npower) = fill(x2, y2, colour_lines(npower,:), ...
        'EdgeColor', 'none', 'FaceAlpha', 0.3);
    plot(bias_mean, prop_mean, 'Color', colour_lines(npower,:), 'LineWidth', 2);

    xlabel('V1 bias diff (T1 - T2)');
    ylabel('Proportion of events detected');
    set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12);
    xlim([-0.1 0.35]);
    ylim([0 1])
end


bias_mean = SO_phase_KDE_bias_difference(2).bias_diff_shifted_mean;
bias_CI_lo = SO_phase_KDE_bias_difference(2).bias_diff_shifted_CI(1,:);
bias_CI_hi = SO_phase_KDE_bias_difference(2).bias_diff_shifted_CI(2,:)
prop_mean = SO_phase_KDE_bias_difference(2).prop_shifted_mean;

hold on;
y2 = [prop_mean, fliplr(prop_mean)];
x2 = [bias_CI_lo, fliplr(bias_CI_hi)];
Fill(end + 1) = fill(x2, y2, 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
plot(bias_mean, prop_mean, 'Color','k', 'LineWidth', 2);

xline(0, '--r');
legend(Fill([1 4 5]), {'SO trough low spindle', 'SO trough high spindle','Shuffled'}, 'box', 'off');
% legend(Fill([1 2 3]), {'SO trough low spindle', 'SO trough high spindle','Shuffled'}, 'box', 'off');


% Save results
save(fullfile(analysis_folder,'V1-HPC sleep reactivation','SO_trough_spindle_power_KDE_bias_difference_based_on_PRE_V1_bias.mat'), 'SO_phase_KDE_bias_difference');
load(fullfile(analysis_folder,'V1-HPC sleep reactivation','SO_trough_spindle_power_KDE_bias_difference_based_on_PRE_V1_bias.mat'), 'SO_phase_KDE_bias_difference');




%%%%% AUC mean + CI bar plot
% Plot layout
fig = figure;
fig.Position = [640 100 281 325]
fig.Name = 'KDE bias V1 AUC SO trough with different spindle power (PRE)';
data = SO_phase_KDE_bias_difference;
n_bins = length(data);
bar_width = 0.3;      % Width of the bars
group_offset = 0.15;    % Distance from the center integer (half the gap between bars)
hold on;
clear BAR
for i = 1:4
    if i > length(data)
        x_shuf = i;
        y_shuf = data(1).AUC_mean_shuffled;
        bar(x_shuf, y_shuf, bar_width, ...
            'FaceColor', 'k', ...
            'FaceAlpha', 0.15, ...
            'EdgeColor', 'none');
    else


        % --- 1. Plot Shuffled Data (Left Bar) ---
        x_shuf = i - group_offset;
        y_shuf = data(i).AUC_mean_shuffled;

        % Calculate Error Deltas (Errorbar requires length relative to mean, not absolute values)
        % CI is [lower, upper]
        neg_err_shuf = y_shuf - data(i).AUC_CI_shuffled(1);
        pos_err_shuf = data(i).AUC_CI_shuffled(2) - y_shuf;

        % Plot Bar
        if i == 1
            BAR(1) = bar(x_shuf, y_shuf, bar_width, ...
                'FaceColor', 'k', ...
                'FaceAlpha', 0.15, ...
                'EdgeColor', 'none');
        else
            bar(x_shuf, y_shuf, bar_width, ...
                'FaceColor', 'k', ...
                'FaceAlpha', 0.15, ...
                'EdgeColor', 'none');
        end

        % Plot Error Bar
        E = errorbar(x_shuf, y_shuf, neg_err_shuf, pos_err_shuf, ...
            'Color', 'k', 'LineWidth', 1.5, 'CapSize', 8, 'LineStyle', 'none');


        % --- 2. Plot Real Data (Right Bar) ---
        x_real = i + group_offset;
        y_real = data(i).AUC_mean;

        neg_err_real = y_real - data(i).AUC_CI(1);
        pos_err_real = data(i).AUC_CI(2) - y_real;

        % Get specific color for this power bin
        % Assuming colour_lines is size [n_bins x 3]
        this_color = colour_lines(i, :);

        % Plot Bar
        BAR(i+1) = bar(x_real, y_real, bar_width, ...
            'FaceColor', this_color, ...
            'FaceAlpha', 0.3, ...
            'EdgeColor', 'none');

        errorbar(x_real, y_real, neg_err_real, pos_err_real, ...
            'Color', 'k', 'LineWidth', 1.5, 'CapSize', 8, 'LineStyle', 'none');
    end
end

hold off;
% Set X-ticks to be centered on the groups
set(gca, 'XTick', 1:nBins);
xlim([0.5, 4 + 0.5]);

% Labels
ylabel('V1 bias AUC');
xlabel('Power Bins');
legend([BAR(1:end)],{'Shuffled','SO trough low spindle','SO trough high spindle'},'box','off')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation','KDE bias difference based on V1 bias'),[])



%%






%%%%%%%%%%% UP and DOWN
nBins = 2;

bins_to_use = bin_centers>0 & bin_centers<0.1;
bins_to_use_shifted = bin_centers>-1 & bin_centers<-0.9;
nBoot = 1000;

% Storage structure
early_late_UP_KDE_bias_difference = struct;
index = {(ripple_info.early_UP_index == 1),(ripple_info.late_UP_index == 1)};
title_names = {'Early UP','Late UP'};

sampled_events = sum(ripple_info.UP_END_index == 1);



colour_lines = [ ...
    241, 182, 218;   % original end (lightest)
    231, 41, 138    % original start (darkest)
] / 256;


% Plot layout
fig = figure;
fig.Position = [640 100 1100 650]
fig.Name = 'KDE bias difference in V1 with ripples during early and late UP ';
tiledlayout(nBins, 3, 'TileSpacing', 'compact');


for npower = 1:nBins
    event_index = (singlet_index + index{npower}) > 1;
    mean_bias = mean(z_bias(bins_to_use, event_index), 'omitnan');
    mean_bias_shifted = mean(z_bias(bins_to_use_shifted, event_index), 'omitnan');
    mean_bias_V1 = mean(z_bias_V1(bins_to_use, event_index), 'omitnan');

    total_events = length(mean_bias);

    % Thresholds for bias
    thresholds = prctile(abs(mean_bias),0:10:100);
    thresholds = thresholds(1:end-1);
    nThresh = length(thresholds);

    % Bootstrap storage
    bias_diff_boot = NaN(nBoot, nThresh);
    prop_events_boot = NaN(nBoot, nThresh);
    % bias_diff_shifted_boot = NaN(nBoot, nThresh);
    % prop_events_shifted_boot = NaN(nBoot, nThresh);

    parfor iBoot = 1:nBoot
        s = RandStream('philox4x32_10', 'Seed', iBoot);
        idx = randi(s, total_events, sampled_events, 1);

        
        % boot_bias_shifted = mean_bias_shifted(idx);
        boot_bias = mean_bias(idx);
        boot_V1 = mean_bias_V1(idx);

        diff_tmp = NaN(1, nThresh);
        prop_tmp = NaN(1, nThresh);

        for i = 1:nThresh
            th = thresholds(i);
            t1 = boot_bias >= th;
            t2 = boot_bias <= -th;

            t1_V1 = boot_V1(t1);
            t2_V1 = boot_V1(t2);

            if ~isempty(t1_V1) && ~isempty(t2_V1)
                diff_tmp(i) = mean(t1_V1) - mean(t2_V1);
            end

            prop_tmp(i) = (sum(t1) + sum(t2)) / sampled_events;
        end

        bias_diff_boot(iBoot, :) = diff_tmp;
        prop_events_boot(iBoot, :) = prop_tmp;
        % 
        % diff_tmp_shifted = NaN(1, nThresh);
        % prop_tmp_shifted = NaN(1, nThresh);
        % 
        % for i = 1:nThresh
        %     th = thresholds(i);
        %     t1 = boot_bias_shifted >= th;
        %     t2 = boot_bias_shifted <= -th;
        % 
        %     t1_V1 = boot_V1(t1);
        %     t2_V1 = boot_V1(t2);
        % 
        %     if ~isempty(t1_V1) && ~isempty(t2_V1)
        %         diff_tmp_shifted(i) = mean(t1_V1) - mean(t2_V1);
        %     end
        % 
        %     prop_tmp_shifted(i) = (sum(t1) + sum(t2)) / total_events;
        % end
        % 
        % bias_diff_shifted_boot(iBoot, :) = diff_tmp_shifted;
        % prop_events_shifted_boot(iBoot, :) = prop_tmp_shifted;
    end

    % Compute stats
    bias_mean = mean(bias_diff_boot, 1, 'omitnan');
    bias_CI_lo = prctile(bias_diff_boot, 2.5, 1);
    bias_CI_hi = prctile(bias_diff_boot, 97.5, 1);

    prop_mean = mean(prop_events_boot, 1, 'omitnan');
    prop_CI_lo = prctile(prop_events_boot, 2.5, 1);
    prop_CI_hi = prctile(prop_events_boot, 97.5, 1);

    % Store results
    early_late_UP_KDE_bias_difference(npower).power_range = [power_thresholds(npower), power_thresholds(npower+1)];
    early_late_UP_KDE_bias_difference(npower).bias_diff_mean = bias_mean;
    early_late_UP_KDE_bias_difference(npower).bias_diff_CI = [bias_CI_lo; bias_CI_hi];
    early_late_UP_KDE_bias_difference(npower).prop_mean = prop_mean;
    early_late_UP_KDE_bias_difference(npower).prop_CI = [prop_CI_lo; prop_CI_hi];
    early_late_UP_KDE_bias_difference(npower).thresholds = thresholds;


    % % Compute stats for shuffled (shifted) bias
    % bias_shifted_mean = mean(bias_diff_shifted_boot, 1, 'omitnan');
    % bias_shifted_CI_lo = prctile(bias_diff_shifted_boot, 2.5, 1);
    % bias_shifted_CI_hi = prctile(bias_diff_shifted_boot, 97.5, 1);
    % 
    % prop_shifted_mean = mean(prop_events_shifted_boot, 1, 'omitnan');
    % prop_shifted_CI_lo = prctile(prop_events_shifted_boot, 2.5, 1);
    % prop_shifted_CI_hi = prctile(prop_events_shifted_boot, 97.5, 1);
    % 
    % % Store shifted (shuffled) results
    % ripple_power_bias_difference(npower).bias_diff_shifted_mean = bias_shifted_mean;
    % ripple_power_bias_difference(npower).bias_diff_shifted_CI = [bias_shifted_CI_lo; bias_shifted_CI_hi];
    % ripple_power_bias_difference(npower).prop_shifted_mean = prop_shifted_mean;
    % ripple_power_bias_difference(npower).prop_shifted_CI = [prop_shifted_CI_lo; prop_shifted_CI_hi];

    % ---- Plot A: Bias difference vs. threshold ----
    nexttile((npower-1)*3 + 1);
    hold on;

    % Real
    x2 = [thresholds, fliplr(thresholds)];
    y2 = [bias_CI_lo, fliplr(bias_CI_hi)];
    fill(x2, y2, colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    plot(thresholds, bias_mean, 'Color', colour_lines(npower,:), 'LineWidth', 2);
    % 
    % % Time-shifted
    % y_shift_lo = bias_shifted_CI_lo;
    % y_shift_hi = bias_shifted_CI_hi;
    % y2s = [y_shift_lo, fliplr(y_shift_hi)];
    % fill(x2, y2s, [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    % plot(thresholds, bias_shifted_mean, 'k-', 'LineWidth', 1.5);

    ylim([-0.15 0.25])
    xlim([0 1.4])
    yline(0,'--r')
    xlabel('HPC bias threshold');
    ylabel('V1 bias diff (T1 - T2)');
    title(sprintf('%s', title_names{npower}));
%     grid on;
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    % ---- Plot B: Proportion vs. Bias Difference (X = bias, Y = proportion) ----
    nexttile((npower-1)*3 + 2);
    hold on;
    valid_idx = isfinite(bias_mean) & isfinite(prop_mean);
    x_vals = bias_mean(valid_idx);
    y_vals = prop_mean(valid_idx);
    x_lo = bias_CI_lo(valid_idx);
    x_hi = bias_CI_hi(valid_idx);
    x_shade = [x_lo, fliplr(x_hi)];
    y_shade = [y_vals, fliplr(y_vals)];

    % Real
    fill(x_shade, y_shade, colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    plot(x_vals, y_vals, '-', 'Color', colour_lines(npower,:), 'LineWidth', 2);

    % % Time-shifted
    % x_vals_shift = bias_shifted_mean(valid_idx);
    % y_vals_shift = prop_shifted_mean(valid_idx);
    % x_lo_s = bias_shifted_CI_lo(valid_idx);
    % x_hi_s = bias_shifted_CI_hi(valid_idx);
    % x_shade_s = [x_lo_s, fliplr(x_hi_s)];
    % y_shade_s = [y_vals_shift, fliplr(y_vals_shift)];

    % fill(x_shade_s, y_shade_s, [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    % plot(x_vals_shift, y_vals_shift, 'k-', 'LineWidth', 1.5);

    xlim([-0.1 0.25])
    xline(0,'--r')
    xlabel('V1 bias diff (T1 - T2)');
    ylabel('Proportion of events detected');
    title('Event Proportion vs. Bias Difference');
%     grid on;
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

    % ---- Plot C: V1 Bias Difference vs. Proportion, shaded CI on Y ----
    nexttile((npower-1)*3 + 3);
    hold on;

    % Real
    valid_idx = isfinite(bias_mean) & isfinite(prop_mean);
    x_vals = bias_mean(valid_idx);
    y_vals = prop_mean(valid_idx);
    y_lo = prop_CI_lo(valid_idx);
    y_hi = prop_CI_hi(valid_idx);
    x_shade = [x_vals, fliplr(x_vals)];
    y_shade = [y_lo, fliplr(y_hi)];

    fill(x_shade, y_shade, colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    plot(x_vals, y_vals, '-', 'Color', colour_lines(npower,:), 'LineWidth', 2);

    % % Time-shifted
    % x_vals_shift = bias_shifted_mean(valid_idx);
    % y_vals_shift = prop_shifted_mean(valid_idx);
    % y_lo_s = prop_shifted_CI_lo(valid_idx);
    % y_hi_s = prop_shifted_CI_hi(valid_idx);
    % x_shade_s = [x_vals_shift, fliplr(x_vals_shift)];
    % y_shade_s = [y_lo_s, fliplr(y_hi_s)];
    % 
    % fill(x_shade_s, y_shade_s, [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    % plot(x_vals_shift, y_vals_shift, 'k-', 'LineWidth', 1.5);

    xlim([-0.1 0.25])
    xline(0,'--r')
    xlabel('V1 bias diff (T1 - T2)');
    ylabel('Proportion of events detected');
    title('Proportion vs. V1 Bias Difference');
%     grid on;
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

end

% Plot layout
fig = figure;
fig.Position = [640 100 2*1100/3 650/2]
fig.Name = 'KDE bias difference in V1 with ripples during early vs late UP ';
% tiledlayout(nBins, 3, 'TileSpacing', 'compact');

colour_lines = [ ...
    241, 182, 218;   % original end (lightest)
    231, 41, 138    % original start (darkest)
] / 256;


nexttile
for npower = [1 2]
    bias_mean = early_late_UP_KDE_bias_difference(npower).bias_diff_mean;
    bias_CI_lo = early_late_UP_KDE_bias_difference(npower).bias_diff_CI(1,:);
    bias_CI_hi = early_late_UP_KDE_bias_difference(npower).bias_diff_CI(2,:)

    hold on;
    x2 = [thresholds, fliplr(thresholds)];
    y2 = [bias_CI_lo, fliplr(bias_CI_hi)];
    Fill(npower) = fill(x2, y2, colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    plot(thresholds, bias_mean, 'Color',colour_lines(npower,:), 'LineWidth', 2);

    xlabel('HPC Bias threshold');
    ylabel('V1 bias diff (T1 - T2)');
    %     title(sprintf('Power bin %d: %.2f–%.2f', npower, power_thresholds(npower), power_thresholds(npower+1)));
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
    ylim([-0.1 0.25])
    %     grid on;
end
yline(0,'--r')
legend(Fill(1:2) ,{'early UP','late UP'},'box','off')

nexttile
for npower = [1 2]
    bias_mean = early_late_UP_KDE_bias_difference(npower).bias_diff_mean;
    bias_CI_lo = early_late_UP_KDE_bias_difference(npower).bias_diff_CI(1,:);
    bias_CI_hi = early_late_UP_KDE_bias_difference(npower).bias_diff_CI(2,:)
    prop_mean = early_late_UP_KDE_bias_difference(npower).prop_mean;

    hold on;
    y2 = [prop_mean, fliplr(prop_mean)];
    x2 = [bias_CI_lo, fliplr(bias_CI_hi)];
    Fill(npower) = fill(x2, y2, colour_lines(npower,:), 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    plot(bias_mean, prop_mean, 'Color',colour_lines(npower,:), 'LineWidth', 2);

    xlabel('HPC Bias threshold');
    ylabel('Proportion of events detected');
    %     title(sprintf('Power bin %d: %.2f–%.2f', npower, power_thresholds(npower), power_thresholds(npower+1)));
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
        xlim([-0.1 0.25])
    %     grid on;
end
xline(0,'--r')
legend(Fill(1:2) ,{'early UP','late UP'},'box','off')

save('early_late_UP_KDE_bias_difference.mat', 'early_late_UP_KDE_bias_difference');
save_all_figures(fullfile(analysis_folder,'V1-HPC sleep reactivation'),[])



