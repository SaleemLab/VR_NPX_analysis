%% Loading extracted SO ripples spindles info for plotting and analysis 
addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\OneDrive\Documents\GitHub\VR_NPX_analysis'))


if exist('D:\corticohippocampal_replay')>0
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

load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_whole.mat'));

load(fullfile(analysis_folder,'KDE_reactivation_DOWN_all_POST.mat'),'KDE_reactivation_DOWN_all')
load(fullfile(analysis_folder,'KDE_reactivation_UP_all_POST.mat'),'KDE_reactivation_UP_all')
load(fullfile(analysis_folder,'KDE_reactivation_all_POST.mat'),'KDE_reactivation_all')
load(fullfile(analysis_folder,'KDE_reactivation_V1_DOWN_all_POST.mat'),'KDE_reactivation_V1_DOWN_all')
load(fullfile(analysis_folder,'KDE_reactivation_V1_UP_all_POST.mat'),'KDE_reactivation_V1_UP_all')
load(fullfile(analysis_folder,'KDE_reactivation_V1_all_POST.mat'),'KDE_reactivation_V1_all')




%%%%% Grabbing bias at DOWN-UP and UP-DOWN transition

all_sessions = max(slow_waves_all(1).DOWN_session_count);
sessions_to_process = 1:all_sessions;

% Initialize output variables
% Initialize output
UP_bias_V1_all = cell(1, 2);
DOWN_bias_V1_all = cell(1, 2);

% Parameters
psth_window = [-1 1];
psth_step = 0.01;  % 10 ms bins
nbins = round(diff(psth_window) / psth_step);  % 200 bins
half_bins = nbins / 2;

for nprobe = 1:2
    up_bias_all = [];
    down_bias_all = [];

    for nsession = 1:length(KDE_reactivation_V1_UP_all(nprobe).bias)
        % === Index setup ===
        up_all = find(slow_waves_all(nprobe).UP_session_count == nsession);
        UP_event_index = intersect(up_all, probability(nprobe).UP_all_index);

        temp_index = find(slow_waves_all(nprobe).DOWN_session_count == nsession);
        [~, ia, ib] = intersect(slow_waves_all(nprobe).DOWN_ints(temp_index, 2), ...
                                slow_waves_all(nprobe).UP_ints(UP_event_index, 1));

        up_ints = slow_waves_all(nprobe).UP_ints(up_all, :);
        previous_DOWN_event_index = ia;
        UP_event_index = find(ismember(up_all, UP_event_index));

        down_all = find(slow_waves_all(nprobe).DOWN_session_count == nsession);
        DOWN_event_index = intersect(down_all, probability(nprobe).DOWN_all_index);
        DOWN_event_index = find(ismember(down_all, DOWN_event_index));
        down_ints = slow_waves_all(nprobe).DOWN_ints(down_all, :);

        % Bias data
        bias_up = KDE_reactivation_V1_UP_all(nprobe).bias{nsession};
        bias_down = KDE_reactivation_V1_DOWN_all(nprobe).bias{nsession};
        time_up = KDE_reactivation_V1_UP_all(nprobe).event_bins{nsession}(:,1)';
        time_down = KDE_reactivation_V1_DOWN_all(nprobe).event_bins{nsession}(:,1)';
        event_id_up = KDE_reactivation_V1_UP_all(nprobe).event_id{nsession};
        event_id_down = KDE_reactivation_V1_DOWN_all(nprobe).event_id{nsession};

        %% === Handle UP-centered events ===
        for i = 1:length(UP_event_index)
            up_bias = nan(1, nbins);

            %% [-1,0]s from previous DOWN (aligned to END)
            if i <= length(previous_DOWN_event_index)
                eid = previous_DOWN_event_index(i);
                event_end = down_ints(eid, 2);
                rel_time = time_down - event_end;

                idx = find(event_id_down == eid & rel_time >= -1 & rel_time < 0);
                rel_times = rel_time(idx);
                bin_idx = round((rel_times + 1) / psth_step) + 1;

                valid = bin_idx > 0 & bin_idx <= half_bins;
                up_bias(bin_idx(valid)) = bias_down(idx(valid));
            end

            %% [0,1]s from current UP (aligned to START)
            eid = UP_event_index(i);
            event_start = up_ints(eid, 1);
            rel_time = time_up - event_start;

            idx = find(event_id_up == eid & rel_time >= 0 & rel_time < 1);
            rel_times = rel_time(idx);
            bin_idx = round(rel_times / psth_step) + half_bins;

            valid = bin_idx > half_bins & bin_idx <= nbins;
            up_bias(bin_idx(valid)) = bias_up(idx(valid));

            up_bias_all = [up_bias_all; up_bias];
        end

        %% === Handle DOWN-centered events ===
        for i = 1:length(DOWN_event_index)
            down_bias = nan(1, nbins);

            %% [-1,0]s from previous UP (aligned to END)
            if i <= length(UP_event_index)
                eid = UP_event_index(i);
                event_end = up_ints(eid, 2);
                rel_time = time_up - event_end;

                idx = find(event_id_up == eid & rel_time >= -1 & rel_time < 0);
                rel_times = rel_time(idx);
                bin_idx = round((rel_times + 1) / psth_step) + 1;

                valid = bin_idx > 0 & bin_idx <= half_bins;
                down_bias(bin_idx(valid)) = bias_up(idx(valid));
            end

            %% [0,1]s from current DOWN (aligned to START)
            eid = DOWN_event_index(i);
            event_start = down_ints(eid, 1);
            rel_time = time_down - event_start;

            idx = find(event_id_down == eid & rel_time >= 0 & rel_time < 1);
            rel_times = rel_time(idx);
            bin_idx = round(rel_times / psth_step) + half_bins;

            valid = bin_idx > half_bins & bin_idx <= nbins;
            down_bias(bin_idx(valid)) = bias_down(idx(valid));

            down_bias_all = [down_bias_all; down_bias];
        end
    end

    UP_bias_V1_all{nprobe} = up_bias_all;
    DOWN_bias_V1_all{nprobe} = down_bias_all;
end

KDE_reactivation_PSTH.V1_UP = [UP_bias_V1_all{1} UP_bias_V1_all{2}];
KDE_reactivation_PSTH.V1_DOWN = [DOWN_bias_V1_all{1} DOWN_bias_V1_all{2}];

% load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_markov_normalised.mat'));
% load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_markov.mat'));