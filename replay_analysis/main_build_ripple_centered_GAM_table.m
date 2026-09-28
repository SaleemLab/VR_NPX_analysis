% MAIN_BUILD_RIPPLE_CENTERED_GAM_TABLE
% Ripple-centered counterpart to MAIN_BUILD_UP_DOWN_INFO_GAM_TABLE: instead of
% one row per UP state (summarizing only the first/last ripple), this produces
% one row per ripple event, tagged with which UP state it belongs to, its
% serial position within that UP (1st, 2nd, 3rd, ... and from the end), and
% its own MUA/reactivation log-odds bias, alongside the UP-level outcome
% features (early/late/nextUP log-odds, next/previous DOWN duration, SO power,
% bilateral lag) broadcast from its parent UP state.
%
% This is the data needed for within-UP, all-ripples analyses (e.g. does
% content match or ripple order predict which ripple ends the UP / predicts
% next UP content, "evidence accumulation" across successive ripples in a UP).
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

%% 1. Load data
load(fullfile(analysis_folder, 'slow_waves_all_POST.mat'));
load(fullfile(analysis_folder, 'ripples_all_POST.mat'));
load(fullfile(analysis_folder, 'V1-HPC sleep interaction', 'SO_ripples_probability_whole_combined.mat'));
probability_psth_whole = probability;
clear probability
load(fullfile(analysis_folder, 'V1-HPC sleep reactivation', 'KDE_reactivation_ripples_PSTH.mat'));
load(fullfile(analysis_folder, 'V1-HPC sleep reactivation', 'KDE_reactivation_PSTH.mat'));

sessions_to_process = 1:max(slow_waves_all(1).UP_session_count);

%% 2. SO peak magnitude at UP->DOWN and DOWN->UP transitions (next/previous DOWN "SO power")
SWpeakmag_UD = [];
SWpeakmag_DU = [];

for nprobe = 1:2
    for nsession = 1:max(slow_waves_all(1).UP_session_count)
        C = intersect(find(slow_waves_all(nprobe).DOWN_session_count == nsession), probability_psth_whole(nprobe).DOWN_all_index);
        SWpeakmag_UD = [SWpeakmag_UD; slow_waves_all(nprobe).SWpeakmag(C)];

        [~, ia] = intersect(slow_waves_all(nprobe).DOWN_ints(slow_waves_all(nprobe).DOWN_session_count == nsession, 2), ...
            slow_waves_all(nprobe).UP_ints(intersect(find(slow_waves_all(nprobe).UP_session_count == nsession), probability_psth_whole(nprobe).UP_all_index), 1));
        temp = find(slow_waves_all(nprobe).DOWN_session_count == nsession);
        SWpeakmag_DU = [SWpeakmag_DU; slow_waves_all(nprobe).SWpeakmag(temp(ia))];
    end
end

%% 3. Per-probe absolute-time UP/DOWN/ripple intervals, spike times, and prev/next DOWN linkage
UP_ints = []; DOWN_ints = []; ripple_ints = []; ripple_peaktimes = [];
V1_MUA_spiketimes = []; HC_MUA_spiketimes = [];
prev_down_idx = []; next_down_idx = [];

for nprobe = 1:2
    V1_MUA_spiketimes{nprobe} = [];
    HC_MUA_spiketimes{nprobe} = [];

    UP_ints{nprobe}   = slow_waves_all(nprobe).UP_ints;
    DOWN_ints{nprobe} = slow_waves_all(nprobe).DOWN_ints;
    ripple_peaktimes{nprobe} = ripples_all(nprobe).peaktimes(ripples_all(nprobe).SWS_index == 1);
    ripple_ints{nprobe}      = [ripples_all(nprobe).onset(ripples_all(nprobe).SWS_index == 1), ...
                                 ripples_all(nprobe).offset(ripples_all(nprobe).SWS_index == 1)];

    for nsession = 1:length(sessions_to_process)
        sess_val = sessions_to_process(nsession);

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

    [is_prev, prev_down_idx{nprobe}] = ismembertol(UP_ints{nprobe}(:,1), DOWN_ints{nprobe}(:,2), 1e-10);
    prev_down_idx{nprobe}(~is_prev) = nan;

    [is_next, next_down_idx{nprobe}] = ismembertol(UP_ints{nprobe}(:,2), DOWN_ints{nprobe}(:,1), 1e-10);
    next_down_idx{nprobe}(~is_next) = nan;
end

%% 4. Assemble merged_event_info across probes (selected UP/DOWN events only)
merged_event_info.UP_ints   = [UP_ints{1}(probability_psth_whole(1).UP_all_index, :); UP_ints{2}(probability_psth_whole(2).UP_all_index, :)];
merged_event_info.DOWN_ints = [DOWN_ints{1}(probability_psth_whole(1).DOWN_all_index, :); DOWN_ints{2}(probability_psth_whole(2).DOWN_all_index, :)];

merged_event_info.UP_hemisphere_id   = [ones(length(probability_psth_whole(1).UP_all_index), 1); 2*ones(length(probability_psth_whole(2).UP_all_index), 1)];
merged_event_info.DOWN_hemisphere_id = [ones(length(probability_psth_whole(1).DOWN_all_index), 1); 2*ones(length(probability_psth_whole(2).DOWN_all_index), 1)];

merged_event_info.ripples_peaktimes     = [ripple_peaktimes{1}; ripple_peaktimes{2}];
merged_event_info.ripples_ints          = [ripple_ints{1}; ripple_ints{2}];
merged_event_info.ripples_hemisphere_id = [ones(length(ripple_peaktimes{1}), 1); ones(length(ripple_peaktimes{2}), 1) * 2];

%% 5. Ipsi/contra UP and DOWN overlap detection and lags (for UD_lags/DU_lags)
event_types = {'UP', 'DOWN'};
for n = 1:2
    L_idx = find(merged_event_info.(sprintf('%s_hemisphere_id', event_types{n})) == 1);
    R_idx = find(merged_event_info.(sprintf('%s_hemisphere_id', event_types{n})) == 2);
    L_ints = merged_event_info.(sprintf('%s_ints', event_types{n}))(L_idx, :);
    R_ints = merged_event_info.(sprintf('%s_ints', event_types{n}))(R_idx, :);

    windows_threshold = [0.05, 0.1, 0.2];

    for ngroup = 1:4
        L_overlap_idx{ngroup} = []; R_overlap_idx{ngroup} = [];
        L_lags{ngroup} = []; R_lags{ngroup} = [];
        all_overlap_idx{ngroup} = []; all_lags{ngroup} = [];

        for iL = 1:size(L_ints, 1)
            L_start = L_ints(iL, 1);
            if ngroup == 4
                L_end = L_ints(iL, 2);
                is_overlap = (R_ints(:,1) <= L_end) & (R_ints(:,2) >= L_start);
            else
                L_end = L_ints(iL, 1) + windows_threshold(ngroup);
                is_overlap = (R_ints(:,1) <= L_end) & (R_ints(:,1) + windows_threshold(ngroup) >= L_start);
            end

            if any(is_overlap)
                overlapping_R = find(is_overlap);
                for j = 1:length(overlapping_R)
                    iR = overlapping_R(j);
                    L_overlap_idx{ngroup}(end+1) = L_idx(iL);
                    R_overlap_idx{ngroup}(end+1) = R_idx(iR);
                    L_lags{ngroup}(end+1) = L_ints(iL,1) - R_ints(iR,1);
                    R_lags{ngroup}(end+1) = R_ints(iR,1) - L_ints(iL,1);
                end
            end
        end

        lags = [L_lags{ngroup} R_lags{ngroup}];
        all_overlap_idx{ngroup} = [L_overlap_idx{ngroup} R_overlap_idx{ngroup}];
        all_lags{ngroup} = lags;
    end

    merged_event_info.(sprintf('%s_overlap_idx_all', event_types{n})) = all_overlap_idx;
    merged_event_info.(sprintf('%s_lags_all', event_types{n})) = all_lags;
end
clear L_overlap_idx R_overlap_idx L_lags R_lags all_overlap_idx all_lags

%% 6. Merge bilateral (L/R) duplicate ripple detections into single events
[event_ids_first, event_ids_second] = merge_bilateral_ripple_events(merged_event_info.ripples_hemisphere_id, ...
    merged_event_info.ripples_peaktimes, 0.05);

merged_event_info.ripples_hemisphere_id = merged_event_info.ripples_hemisphere_id(event_ids_first);
merged_event_info.ripples_peaktimes     = merged_event_info.ripples_peaktimes(event_ids_first, :);
merged_event_info.ripples_ints          = merged_event_info.ripples_ints(event_ids_first, :);

ripplePower = [ripples_all(1).peak_zscore(ripples_all(1).SWS_index); ripples_all(2).peak_zscore(ripples_all(2).SWS_index)];
merged_event_info.ripples_power = mean([ripplePower(event_ids_first), ripplePower(event_ids_second)], 2);

%% 7. Session and subject metadata
UP_session_count = [slow_waves_all(1).UP_session_count(probability_psth_whole(1).UP_all_index); ...
                     slow_waves_all(2).UP_session_count(probability_psth_whole(2).UP_all_index)];
subject_id = str2double(cellstr(slow_waves_all(1).subject(UP_session_count, end-1:end)));
[~, ~, subject_id] = unique(subject_id);

merged_event_info.session_id = UP_session_count;
merged_event_info.subject_id = subject_id;

%% 8. Per-ripple HC and V1 reactivation log-odds bias (post-ripple and pre-ripple windows)
timebin = 0.01;
time_windows = [-1 1];
bin_edges = time_windows(1):timebin:time_windows(2);
bin_centers = bin_edges(1:end-1) + timebin/2;

z_bias    = KDE_reactivation_ripples_PSTH.HPC_z_logodds_ripples';
z_bias_V1 = KDE_reactivation_ripples_PSTH.V1_z_logodds_ripples';

z_bias1 = z_bias(isfinite(z_bias));
z_bias(z_bias >= inf)  = prctile(z_bias1, 99.5);
z_bias(z_bias <= -inf) = prctile(z_bias1, 0.5);

z_bias1 = z_bias(isfinite(z_bias_V1));
z_bias_V1(z_bias_V1 >= inf)  = prctile(z_bias1, 99.5);
z_bias_V1(z_bias_V1 <= -inf) = prctile(z_bias1, 0.5);

z_bias    = z_bias    + KDE_reactivation_ripples_PSTH.nan_mask';
z_bias_V1 = z_bias_V1 + KDE_reactivation_ripples_PSTH.nan_mask';

z_bias    = z_bias(:, event_ids_first);
z_bias_V1 = z_bias_V1(:, event_ids_first);

ripple_HC_logodds     = mean(z_bias(bin_centers>0 & bin_centers<0.1, :), 'omitnan')';
ripple_HC_logodds_PRE = mean(z_bias(bin_centers>-0.1 & bin_centers<0, :), 'omitnan')';
ripple_V1_logodds     = mean(z_bias_V1(bin_centers>0 & bin_centers<0.1, :), 'omitnan')';
ripple_V1_logodds_PRE = mean(z_bias_V1(bin_centers>-0.1 & bin_centers<0, :), 'omitnan')';

%% 9. UP/DOWN-epoch reactivation log-odds PSTHs (per probe, inf-clipped)
UP_HPC_log_odds = cell(1,2); UP_V1_log_odds = cell(1,2);
DOWN_HPC_log_odds = cell(1,2); DOWN_V1_log_odds = cell(1,2);

for nprobe = 1:2
    UP_V1_log_odds{nprobe} = KDE_reactivation_PSTH(nprobe).V1_UP_log_odds;
    temp1 = UP_V1_log_odds{nprobe}(isfinite(UP_V1_log_odds{nprobe}));
    UP_V1_log_odds{nprobe}(UP_V1_log_odds{nprobe} >= inf)  = prctile(temp1, 99.5);
    UP_V1_log_odds{nprobe}(UP_V1_log_odds{nprobe} <= -inf) = prctile(temp1, 0.5);

    UP_HPC_log_odds{nprobe} = KDE_reactivation_PSTH(nprobe).HPC_UP_log_odds;
    temp1 = UP_HPC_log_odds{nprobe}(isfinite(UP_HPC_log_odds{nprobe}));
    UP_HPC_log_odds{nprobe}(UP_HPC_log_odds{nprobe} >= inf)  = prctile(temp1, 99.5);
    UP_HPC_log_odds{nprobe}(UP_HPC_log_odds{nprobe} <= -inf) = prctile(temp1, 0.5);

    DOWN_V1_log_odds{nprobe} = KDE_reactivation_PSTH(nprobe).V1_DOWN_log_odds;
    temp1 = DOWN_V1_log_odds{nprobe}(isfinite(DOWN_V1_log_odds{nprobe}));
    DOWN_V1_log_odds{nprobe}(DOWN_V1_log_odds{nprobe} >= inf)  = prctile(temp1, 99.5);
    DOWN_V1_log_odds{nprobe}(DOWN_V1_log_odds{nprobe} <= -inf) = prctile(temp1, 0.5);

    DOWN_HPC_log_odds{nprobe} = KDE_reactivation_PSTH(nprobe).HPC_DOWN_log_odds;
    temp1 = DOWN_HPC_log_odds{nprobe}(isfinite(DOWN_HPC_log_odds{nprobe}));
    DOWN_HPC_log_odds{nprobe}(DOWN_HPC_log_odds{nprobe} >= inf)  = prctile(temp1, 99.5);
    DOWN_HPC_log_odds{nprobe}(DOWN_HPC_log_odds{nprobe} <= -inf) = prctile(temp1, 0.5);
end

%% 10. Previous/next DOWN duration (per UP event)
nUP = size(merged_event_info.UP_ints, 1);
previous_DOWN_duration = nan(nUP, 1);
next_DOWN_duration     = nan(nUP, 1);

row = 0;
for nprobe = 1:2
    UP_idx = probability_psth_whole(nprobe).UP_all_index;
    n = length(UP_idx);
    r = row+1:row+n;

    pdi = prev_down_idx{nprobe}(UP_idx);
    ndi = next_down_idx{nprobe}(UP_idx);
    valid_p = ~isnan(pdi);
    valid_n = ~isnan(ndi);

    prev_dur = nan(n,1); next_dur = nan(n,1);
    prev_dur(valid_p) = DOWN_ints{nprobe}(pdi(valid_p),2) - DOWN_ints{nprobe}(pdi(valid_p),1);
    next_dur(valid_n) = DOWN_ints{nprobe}(ndi(valid_n),2) - DOWN_ints{nprobe}(ndi(valid_n),1);

    previous_DOWN_duration(r) = prev_dur;
    next_DOWN_duration(r) = next_dur;
    row = row + n;
end

%% 11. Early/late UP and nextUP reactivation log-odds (100ms windows, per UP event)
time_window = 0.1;
tvec = KDE_reactivation_PSTH(1).tvec;
tidx_late  = tvec > -time_window & tvec < 0;
tidx_early = tvec > 0 & tvec < time_window;

late_UP_log_odds = nan(nUP,1); late_UP_V1_log_odds = nan(nUP,1);
early_UP_log_odds = nan(nUP,1); early_UP_V1_log_odds = nan(nUP,1);
next_early_UP_log_odds = nan(nUP,1); next_early_UP_V1_log_odds = nan(nUP,1);

row = 0;
for nprobe = 1:2
    UP_idx = probability_psth_whole(nprobe).UP_all_index;
    n = length(UP_idx);
    r = row+1:row+n;

    ndi = next_down_idx{nprobe}(UP_idx);
    valid_n = ~isnan(ndi);
    late_UP_log_odds(r(valid_n))     = mean(DOWN_HPC_log_odds{nprobe}(ndi(valid_n), tidx_late), 2, 'omitnan');
    late_UP_V1_log_odds(r(valid_n))  = mean(DOWN_V1_log_odds{nprobe}(ndi(valid_n), tidx_late), 2, 'omitnan');

    early_UP_log_odds(r)    = mean(UP_HPC_log_odds{nprobe}(UP_idx, tidx_early), 2, 'omitnan');
    early_UP_V1_log_odds(r) = mean(UP_V1_log_odds{nprobe}(UP_idx, tidx_early), 2, 'omitnan');

    next_idx = UP_idx + 1;
    valid_next = next_idx <= size(UP_HPC_log_odds{nprobe}, 1);
    safe_next_idx = next_idx; safe_next_idx(~valid_next) = 1;
    tmp1 = mean(UP_HPC_log_odds{nprobe}(safe_next_idx, tidx_early), 2, 'omitnan'); tmp1(~valid_next) = nan;
    tmp2 = mean(UP_V1_log_odds{nprobe}(safe_next_idx, tidx_early), 2, 'omitnan');  tmp2(~valid_next) = nan;
    next_early_UP_log_odds(r)    = tmp1;
    next_early_UP_V1_log_odds(r) = tmp2;

    row = row + n;
end

%% 12. Ipsi/contra UP-DOWN lags (relies on UP_all_index/DOWN_all_index sharing per-cycle pairing/order)
UD_lags = nan(nUP,1);
UD_lags(merged_event_info.DOWN_overlap_idx_all{end}) = abs(merged_event_info.DOWN_lags_all{end});

DU_lags = nan(nUP,1);
DU_lags(merged_event_info.UP_overlap_idx_all{end}) = abs(merged_event_info.UP_lags_all{end});

%% 13. Build the ripple-centered table
fprintf('Building ripple-centered GAM table...\n');
ripple_tbl = build_ripple_centered_GAM_table(merged_event_info, V1_MUA_spiketimes, HC_MUA_spiketimes, ...
    ripple_HC_logodds, ripple_HC_logodds_PRE, ripple_V1_logodds, ripple_V1_logodds_PRE, ...
    early_UP_log_odds, early_UP_V1_log_odds, ...
    late_UP_log_odds, late_UP_V1_log_odds, ...
    next_early_UP_log_odds, next_early_UP_V1_log_odds, ...
    next_DOWN_duration, previous_DOWN_duration, ...
    SWpeakmag_UD, SWpeakmag_DU, ...
    UD_lags, DU_lags);

my_folder = 'C:/Users/masah/Documents/GitHub/VR_NPX_analysis/UP_DOWN_ripple_GAM';
if ~exist(my_folder, 'dir')
    mkdir(my_folder);
end
writetable(ripple_tbl, fullfile(my_folder, 'ripple_centered_GAM.csv'));

fprintf('Ripple-centered GAM table built and saved: %d ripple rows across %d UP events.\n', ...
    height(ripple_tbl), nUP);
