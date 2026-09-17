% MAIN_BUILD_UP_DOWN_INFO_GAM_TABLE
% Builds the UP_DOWN_info struct and the per-UP-event GAM table
% (UP_DOWN_info_GAM.csv) used for R survival/GAM analysis.
%
% MUA-derived features (last ripple / past ripples / non-ripple MUA,
% first vs second half of UP MUA, etc.) are produced entirely by
% BUILD_LAST_RIPPLE_TO_DOWN_SUMMARY_TABLE (from raw spike times), so that
% logic is not duplicated here. This script adds the complementary,
% non-MUA information: ripple identity/timing, reactivation log-odds bias,
% periripple spindle amplitude (via BUILD_RIPPLE_LOG_ODDS_SPINDLE_TABLE),
% slow-oscillation peak/delta features, ipsi/contra UP-DOWN lags, and
% previous/next DOWN duration.
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
load(fullfile(analysis_folder, 'V1-HPC sleep interaction', 'SO_spindles_probability_whole.mat'));
spindle_probability_psth_whole = probability;
load(fullfile(analysis_folder, 'V1-HPC sleep reactivation', 'KDE_reactivation_ripples_PSTH.mat'));
load(fullfile(analysis_folder, 'V1-HPC sleep reactivation', 'KDE_reactivation_PSTH.mat'));


clear probability
analysis_folder = 'P:\corticohippocampal_replay';
load(fullfile(analysis_folder, 'periripple_LFP_info_V1.mat'));
sessions_to_process = 1:max(slow_waves_all(1).UP_session_count);

%% 2. Reference shank per session/probe (for cortical SO peak z-score lookup)
cortex_ref_shank = [];
for nsession = 1:max(ripples_all(1).session_count)
    for probe_no = 1:2
        cortex_ref_shank(nsession, probe_no) = find(slow_waves_all(probe_no).shank_id{nsession} == ...
            slow_waves_all(probe_no).shank{nsession}(slow_waves_all(probe_no).channel{nsession} == slow_waves_all(probe_no).best_channel(nsession)) ...
            & slow_waves_all(probe_no).probe_hemisphere{nsession} == probe_no);
    end
end

%% 3. Ipsi/contra SO delta z-score and peak magnitude at UP-DOWN and DOWN-UP transitions
ipsi_Delta_peaks_zscore_UD = [];
ipsi_Delta_peaks_zscore_DU = [];
SWpeakmag_UD = [];
SWpeakmag_DU = [];
contra_Delta_peaks_zscore_UD = [];
contra_Delta_peaks_zscore_DU = [];

for nprobe = 1:2
    for nsession = 1:max(slow_waves_all(1).UP_session_count)
        [C, ia, ib] = intersect(find(slow_waves_all(nprobe).DOWN_session_count == nsession), probability_psth_whole(nprobe).DOWN_all_index);

        ipsi_Delta_peaks_zscore_UD = [ipsi_Delta_peaks_zscore_UD; ...
            slow_waves_all(nprobe).DOWN_peaks_zscore{nsession}(cortex_ref_shank(nsession, nprobe), ia)'];
        SWpeakmag_UD = [SWpeakmag_UD; slow_waves_all(nprobe).SWpeakmag(C)];

        [C, ia, ib] = intersect(slow_waves_all(nprobe).DOWN_ints(slow_waves_all(nprobe).DOWN_session_count == nsession, 2), ...
            slow_waves_all(nprobe).UP_ints(intersect(find(slow_waves_all(nprobe).UP_session_count == nsession), probability_psth_whole(nprobe).UP_all_index), 1));

        ipsi_Delta_peaks_zscore_DU = [ipsi_Delta_peaks_zscore_DU; slow_waves_all(nprobe).DOWN_peaks_zscore{nsession}(cortex_ref_shank(nsession, nprobe), ia)'];
        temp = find(slow_waves_all(nprobe).DOWN_session_count == nsession);
        SWpeakmag_DU = [SWpeakmag_DU; slow_waves_all(nprobe).SWpeakmag(temp(ia))];

        [C, ia, ib] = intersect(find(slow_waves_all(nprobe).DOWN_session_count == nsession), probability_psth_whole(nprobe).DOWN_all_index);
        contra_Delta_peaks_zscore_UD = [contra_Delta_peaks_zscore_UD; ...
            slow_waves_all(nprobe).DOWN_peaks_zscore{nsession}(cortex_ref_shank(nsession, abs(nprobe - 3)), ia)'];

        [C, ia, ib] = intersect(slow_waves_all(nprobe).DOWN_ints(slow_waves_all(nprobe).DOWN_session_count == nsession, 2), ...
            slow_waves_all(nprobe).UP_ints(intersect(find(slow_waves_all(nprobe).UP_session_count == nsession), probability_psth_whole(nprobe).UP_all_index), 1));
        contra_Delta_peaks_zscore_DU = [contra_Delta_peaks_zscore_DU; slow_waves_all(nprobe).DOWN_peaks_zscore{nsession}(cortex_ref_shank(nsession, abs(nprobe - 3)), ia)'];
    end
end

%% 4. Per-probe absolute-time UP/DOWN/ripple intervals, spike times, and prev/next DOWN linkage
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

    nUPprobe = length(UP_ints{nprobe});
    [is_prev, prev_down_idx{nprobe}] = ismembertol(UP_ints{nprobe}(:,1), DOWN_ints{nprobe}(:,2), 1e-10);
    prev_down_idx{nprobe}(~is_prev) = nan;

    [is_next, next_down_idx{nprobe}] = ismembertol(UP_ints{nprobe}(:,2), DOWN_ints{nprobe}(:,1), 1e-10);
    next_down_idx{nprobe}(~is_next) = nan;
end

%% 5. Assemble merged_event_info across probes (selected UP/DOWN events only)
merged_event_info.UP_ints   = [UP_ints{1}(probability_psth_whole(1).UP_all_index, :); UP_ints{2}(probability_psth_whole(2).UP_all_index, :)];
merged_event_info.DOWN_ints = [DOWN_ints{1}(probability_psth_whole(1).DOWN_all_index, :); DOWN_ints{2}(probability_psth_whole(2).DOWN_all_index, :)];

merged_event_info.UP_hemisphere_id   = [ones(length(probability_psth_whole(1).UP_all_index), 1); 2*ones(length(probability_psth_whole(2).UP_all_index), 1)];
merged_event_info.DOWN_hemisphere_id = [ones(length(probability_psth_whole(1).DOWN_all_index), 1); 2*ones(length(probability_psth_whole(2).DOWN_all_index), 1)];

merged_event_info.ripples_peaktimes     = [ripple_peaktimes{1}; ripple_peaktimes{2}];
merged_event_info.ripples_ints          = [ripple_ints{1}; ripple_ints{2}];
merged_event_info.ripples_hemisphere_id = [ones(length(ripple_peaktimes{1}), 1); 2*ones(length(ripple_peaktimes{2}), 1)];

% Onset of the next UP state on the same probe/hemisphere (NaN at session end), in merged row order
UP_next_onset = [];
for nprobe = 1:2
    UP_idx = probability_psth_whole(nprobe).UP_all_index;
    next_onset_probe = nan(length(UP_idx), 1);
    next_onset_probe(1:end-1) = UP_ints{nprobe}(UP_idx(2:end), 1);
    UP_next_onset = [UP_next_onset; next_onset_probe];
end

%% 6. Ipsi/contra UP and DOWN overlap detection and lags (for UD_lags/DU_lags)
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
        all_overlap_idx{ngroup} = []; overlap_idx{ngroup} = []; non_overlap_idx{ngroup} = [];
        all_lags{ngroup} = [];

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
        [overlap_idx{ngroup}, ~] = unique(all_overlap_idx{ngroup});
        all_lags{ngroup} = lags;

        nonoverlap_L_idx = setdiff(L_idx, unique(L_overlap_idx{ngroup}));
        nonoverlap_R_idx = setdiff(R_idx, unique(R_overlap_idx{ngroup}));
        non_overlap_idx{ngroup} = unique([nonoverlap_L_idx; nonoverlap_R_idx]);
    end

    merged_event_info.(sprintf('%s_overlap_idx_all', event_types{n})) = all_overlap_idx;
    merged_event_info.(sprintf('%s_non_overlap_idx', event_types{n})) = non_overlap_idx;
    merged_event_info.(sprintf('%s_overlap_idx', event_types{n})) = overlap_idx;
    merged_event_info.(sprintf('%s_lags_all', event_types{n})) = all_lags;
end
clear L_overlap_idx R_overlap_idx L_lags R_lags all_overlap_idx overlap_idx non_overlap_idx all_lags

%% 7. Merge bilateral (L/R) duplicate ripple detections into single events
[event_ids_first, event_ids_second] = merge_bilateral_ripple_events(merged_event_info.ripples_hemisphere_id, ...
    merged_event_info.ripples_peaktimes, 0.05);

merged_event_info.ripples_hemisphere_id = merged_event_info.ripples_hemisphere_id(event_ids_first);
merged_event_info.ripples_peaktimes     = merged_event_info.ripples_peaktimes(event_ids_first, :);
merged_event_info.ripples_ints          = merged_event_info.ripples_ints(event_ids_first, :);

ripplePower = [ripples_all(1).peak_zscore(ripples_all(1).SWS_index); ripples_all(2).peak_zscore(ripples_all(2).SWS_index)];
merged_event_info.ripples_power = mean([ripplePower(event_ids_first), ripplePower(event_ids_second)], 2);

%% 8. Session and subject metadata
UP_session_count = [slow_waves_all(1).UP_session_count(probability_psth_whole(1).UP_all_index); ...
                     slow_waves_all(2).UP_session_count(probability_psth_whole(2).UP_all_index)];
subject_id = str2double(cellstr(slow_waves_all(1).subject(UP_session_count, end-1:end)));
[~, ~, subject_id] = unique(subject_id);

merged_event_info.session_id = UP_session_count;
merged_event_info.subject_id = subject_id;

%% 9. Periripple spindle amplitude (deduplicated to merged ripple order)
spindle_amplitude1 = [periripple_LFP_info_V1(1).spindle_amplitude{1}(:, ripples_all(1).SWS_index==1), periripple_LFP_info_V1(2).spindle_amplitude{1}(:, ripples_all(2).SWS_index==1)];
spindle_amplitude2 = [periripple_LFP_info_V1(1).spindle_amplitude{2}(:, ripples_all(1).SWS_index==1), periripple_LFP_info_V1(2).spindle_amplitude{2}(:, ripples_all(2).SWS_index==1)];
spindle_amplitude = nan([size(spindle_amplitude1), 2]);
spindle_amplitude(:,:,1) = spindle_amplitude1;
spindle_amplitude(:,:,2) = spindle_amplitude2;
spindle_amplitude_temporal = spindle_amplitude(:, event_ids_first, :);
spindle_tvec = periripple_LFP_info_V1(1).tvec;
clear spindle_amplitude spindle_amplitude1 spindle_amplitude2

%% 10. Per-ripple reactivation log-odds bias (post-ripple and pre-ripple window means)
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

mean_z_bias         = mean(z_bias(bin_centers>0 & bin_centers<0.1, :), 'omitnan');
mean_z_bias_PRE     = mean(z_bias(bin_centers>-0.1 & bin_centers<0, :), 'omitnan');
mean_z_bias_V1      = mean(z_bias_V1(bin_centers>0 & bin_centers<0.1, :), 'omitnan');
mean_z_bias_V1_PRE  = mean(z_bias_V1(bin_centers>-0.1 & bin_centers<0, :), 'omitnan');

%% 11. UP/DOWN-epoch reactivation log-odds PSTHs (per probe, inf-clipped)
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

%% 12. MUA-derived per-UP features (last ripple / past ripples / non-ripple / half-UP splits)
fprintf('Building MUA-derived last-ripple-to-DOWN features...\n');
% T = build_last_ripple_to_DOWN_summary_table(merged_event_info, V1_MUA_spiketimes, HC_MUA_spiketimes, ...
%     'time_reference', 'offset', 'non_ripple_scope', 'entire_UP');
% T = build_last_ripple_to_DOWN_summary_table(merged_event_info, ...
%     V1_MUA_spiketimes, HC_MUA_spiketimes, ...
%     'time_reference', 'peak', ...
%     'non_ripple_scope', 'prior_to_last_ripple');

T = build_last_ripple_to_DOWN_summary_table(merged_event_info, ...
    V1_MUA_spiketimes, HC_MUA_spiketimes, ...
    'time_reference', 'offset', ...
    'non_ripple_scope', 'prior_to_last_ripple');

%% 13. Non-MUA per-UP features: ripple identity/timing, log-odds bias, spindle amplitude
ripple_tbl = build_ripple_log_odds_table(merged_event_info, ...
    mean_z_bias, mean_z_bias_PRE, mean_z_bias_V1, mean_z_bias_V1_PRE, ...
    spindle_amplitude_temporal, spindle_tvec, UP_next_onset,'time_reference', 'peak');

ripple_tbl = build_ripple_log_odds_table(merged_event_info, ...
    mean_z_bias, mean_z_bias_PRE, mean_z_bias_V1, mean_z_bias_V1_PRE, ...
    spindle_amplitude_temporal, spindle_tvec, UP_next_onset,'time_reference', 'boundary');

% T.UP_to_last_ripple = T.up_duration-T.last_ripple_duration-T.last_ripple_to_UP_term;

% scatter(T.UP_to_last_ripple,ripple_tbl.time_to_last_ripples)

%% 14. Previous/next DOWN duration (per UP event)
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

%% 15. Early/late UP and DOWN reactivation log-odds (100ms windows, per UP event)
time_window = 0.1;
tvec = KDE_reactivation_PSTH(1).tvec;
tidx_late  = tvec > -time_window & tvec < 0;
tidx_early = tvec > 0 & tvec < time_window;

late_UP_log_odds = nan(nUP,1); late_UP_V1_log_odds = nan(nUP,1);
early_UP_log_odds = nan(nUP,1); early_UP_V1_log_odds = nan(nUP,1);
DOWN_log_odds = nan(nUP,1); DOWN_V1_log_odds_vec = nan(nUP,1);
previous_late_UP_log_odds = nan(nUP,1); previous_late_UP_V1_log_odds = nan(nUP,1);
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
    DOWN_log_odds(r(valid_n))        = mean(DOWN_HPC_log_odds{nprobe}(ndi(valid_n), tidx_early), 2, 'omitnan');
    DOWN_V1_log_odds_vec(r(valid_n)) = mean(DOWN_V1_log_odds{nprobe}(ndi(valid_n), tidx_early), 2, 'omitnan');

    early_UP_log_odds(r)    = mean(UP_HPC_log_odds{nprobe}(UP_idx, tidx_early), 2, 'omitnan');
    early_UP_V1_log_odds(r) = mean(UP_V1_log_odds{nprobe}(UP_idx, tidx_early), 2, 'omitnan');

    pdi = prev_down_idx{nprobe}(UP_idx);
    valid_p = ~isnan(pdi);
    previous_late_UP_log_odds(r(valid_p))    = mean(DOWN_HPC_log_odds{nprobe}(pdi(valid_p), tidx_late), 2, 'omitnan');
    previous_late_UP_V1_log_odds(r(valid_p)) = mean(DOWN_V1_log_odds{nprobe}(pdi(valid_p), tidx_late), 2, 'omitnan');

    next_idx = UP_idx + 1;
    valid_next = next_idx <= size(UP_HPC_log_odds{nprobe}, 1);
    safe_next_idx = next_idx; safe_next_idx(~valid_next) = 1;
    tmp1 = mean(UP_HPC_log_odds{nprobe}(safe_next_idx, tidx_early), 2, 'omitnan'); tmp1(~valid_next) = nan;
    tmp2 = mean(UP_V1_log_odds{nprobe}(safe_next_idx, tidx_early), 2, 'omitnan');  tmp2(~valid_next) = nan;
    next_early_UP_log_odds(r)    = tmp1;
    next_early_UP_V1_log_odds(r) = tmp2;

    row = row + n;
end

%% 16. Ipsi/contra UP-DOWN lags (relies on UP_all_index/DOWN_all_index sharing per-cycle pairing/order)
UD_lags = nan(nUP,1); UD_lags_signed = nan(nUP,1);
UD_lags(merged_event_info.DOWN_overlap_idx_all{end})        = abs(merged_event_info.DOWN_lags_all{end});
UD_lags_signed(merged_event_info.DOWN_overlap_idx_all{end})  = merged_event_info.DOWN_lags_all{end};

DU_lags = nan(nUP,1); DU_lags_signed = nan(nUP,1);
DU_lags(merged_event_info.UP_overlap_idx_all{end})        = abs(merged_event_info.UP_lags_all{end});
DU_lags_signed(merged_event_info.UP_overlap_idx_all{end}) = merged_event_info.UP_lags_all{end};

V1_direction = [-1*ones(length(probability_psth_whole(1).UP_all_index),1); ones(length(probability_psth_whole(2).UP_all_index),1)];

%% 17. Assemble UP_DOWN_info
UP_DOWN_info = table2struct(T, 'ToScalar', true);

ripple_fields = ripple_tbl.Properties.VariableNames;
for i = 1:numel(ripple_fields)
    UP_DOWN_info.(ripple_fields{i}) = ripple_tbl.(ripple_fields{i});
end

UP_DOWN_info.previous_DOWN_duration = previous_DOWN_duration;
UP_DOWN_info.next_DOWN_duration     = next_DOWN_duration;

UP_DOWN_info.ipsi_Delta_peaks_zscore_UD   = ipsi_Delta_peaks_zscore_UD;
UP_DOWN_info.ipsi_Delta_peaks_zscore_DU   = ipsi_Delta_peaks_zscore_DU;
UP_DOWN_info.contra_Delta_peaks_zscore_UD = contra_Delta_peaks_zscore_UD;
UP_DOWN_info.contra_Delta_peaks_zscore_DU = contra_Delta_peaks_zscore_DU;
UP_DOWN_info.SWpeakmag_UD = SWpeakmag_UD;
UP_DOWN_info.SWpeakmag_DU = SWpeakmag_DU;

UP_DOWN_info.late_UP_log_odds = late_UP_log_odds;
UP_DOWN_info.late_UP_V1_log_odds = late_UP_V1_log_odds;
UP_DOWN_info.early_UP_log_odds = early_UP_log_odds;
UP_DOWN_info.early_UP_V1_log_odds = early_UP_V1_log_odds;
UP_DOWN_info.DOWN_log_odds = DOWN_log_odds;
UP_DOWN_info.DOWN_V1_log_odds = DOWN_V1_log_odds_vec;
UP_DOWN_info.previous_late_UP_log_odds = previous_late_UP_log_odds;
UP_DOWN_info.previous_late_UP_V1_log_odds = previous_late_UP_V1_log_odds;
UP_DOWN_info.next_early_UP_log_odds = next_early_UP_log_odds;
UP_DOWN_info.next_early_UP_V1_log_odds = next_early_UP_V1_log_odds;

UP_DOWN_info.UD_lags = UD_lags;
UP_DOWN_info.DU_lags = DU_lags;
UP_DOWN_info.UD_lags_signed = UD_lags_signed;
UP_DOWN_info.DU_lags_signed = DU_lags_signed;

UP_DOWN_info.V1_direction = V1_direction;

UP_DOWN_info.ipsi_spindles_UP    = [nansum(spindle_probability_psth_whole(1).L_spindles_UP(:,50:60)')' > 0; nansum(spindle_probability_psth_whole(2).R_spindles_UP(:,50:60)')' > 0];
UP_DOWN_info.contra_spindles_UP  = [nansum(spindle_probability_psth_whole(1).R_spindles_UP(:,50:60)')' > 0; nansum(spindle_probability_psth_whole(2).L_spindles_UP(:,50:60)')' > 0];
UP_DOWN_info.ipsi_spindles_DOWN  = [nansum(spindle_probability_psth_whole(1).L_spindles_DOWN(:,50:60)')' > 0; nansum(spindle_probability_psth_whole(2).R_spindles_DOWN(:,50:60)')' > 0];
UP_DOWN_info.contra_spindles_DOWN = [nansum(spindle_probability_psth_whole(1).R_spindles_DOWN(:,50:60)')' > 0; nansum(spindle_probability_psth_whole(2).L_spindles_DOWN(:,50:60)')' > 0];
UP_DOWN_info1 = UP_DOWN_info;

% load(fullfile(analysis_folder, 'V1-HPC sleep reactivation', 'UP_DOWN_info_100ms.mat'), 'UP_DOWN_info');

save(fullfile(analysis_folder, 'V1-HPC sleep reactivation', 'UP_DOWN_info_100ms_offset.mat'), 'UP_DOWN_info');

% save(fullfile(analysis_folder, 'V1-HPC sleep reactivation', 'UP_DOWN_info_100ms_peak.mat'), 'UP_DOWN_info');


% 
histogram(UP_DOWN_info1.time_from_last_ripples,-0.03:0.005:2,'Normalization','probability');hold on;histogram(UP_DOWN_info.time_from_last_ripples,-0.03:0.005:2,'Normalization','probability')

scatter(UP_DOWN_info1.last_ripples_log_odds,UP_DOWN_info.last_ripples_log_odds_UP)

scatter(UP_DOWN_info1.time_from_last_ripples,UP_DOWN_info.time_from_last_ripples_UP)

%% 18. Build GAM table for R
load(fullfile(analysis_folder, 'V1-HPC sleep reactivation', 'UP_DOWN_info_100ms.mat'), 'UP_DOWN_info');
load(fullfile(analysis_folder, 'V1-HPC sleep reactivation', 'UP_DOWN_info_100ms_offset.mat'), 'UP_DOWN_info');
% load(fullfile(analysis_folder, 'V1-HPC sleep reactivation', 'UP_DOWN_info_100ms_peak.mat'), 'UP_DOWN_info');

% load(fullfile(analysis_folder, 'V1-HPC sleep reactivation', 'UP_DOWN_info_100ms_peak2.mat'), 'UP_DOWN_info');

UP_DOWN_info.time_from_last_ripple_to_next_UP = UP_DOWN_info.next_DOWN_duration + UP_DOWN_info.time_from_last_ripples;

lastRippleNormalisedUP  = UP_DOWN_info.time_to_last_ripples ./ UP_DOWN_info.up_duration;
firstRippleNormalisedUP = UP_DOWN_info.time_to_first_ripples ./ UP_DOWN_info.up_duration;
lastRippleNormalisedUP(lastRippleNormalisedUP > 1) = 1;
firstRippleNormalisedUP(firstRippleNormalisedUP > 1) = 1;
lastRippleNormalisedUP(lastRippleNormalisedUP < 0) = 0;
firstRippleNormalisedUP(firstRippleNormalisedUP < 0) = 0;

tbl = table(...
    UP_DOWN_info.time_to_first_ripples(:), ...
    UP_DOWN_info.time_from_first_ripples(:), ...
    UP_DOWN_info.time_to_last_ripples(:), ...
    UP_DOWN_info.time_from_last_ripples(:), ...
    UP_DOWN_info.time_from_last_ripple_to_next_UP(:), ...
    UP_DOWN_info.last_ripple_duration,...
    UP_DOWN_info.V1_direction(:), ...
    UP_DOWN_info.first_ripples_power(:), ...
    UP_DOWN_info.last_ripples_power(:), ...
    UP_DOWN_info.next_early_UP_V1_log_odds(:), ...
    UP_DOWN_info.early_UP_V1_log_odds(:), ...
    UP_DOWN_info.late_UP_V1_log_odds(:), ...
    UP_DOWN_info.previous_late_UP_V1_log_odds(:), ...
    UP_DOWN_info.first_ripples_V1_log_odds(:), ...
    UP_DOWN_info.first_ripples_PRE_V1_log_odds(:), ...
    UP_DOWN_info.last_ripples_V1_log_odds(:), ...
    UP_DOWN_info.last_ripples_PRE_V1_log_odds(:), ...
    UP_DOWN_info.next_early_UP_log_odds(:), ...
    UP_DOWN_info.early_UP_log_odds(:), ...
    UP_DOWN_info.late_UP_log_odds(:), ...
    UP_DOWN_info.previous_late_UP_log_odds(:), ...
    UP_DOWN_info.first_ripples_log_odds(:), ...
    UP_DOWN_info.first_ripples_PRE_log_odds(:), ...
    UP_DOWN_info.last_ripples_log_odds(:), ...
    UP_DOWN_info.last_ripples_PRE_log_odds(:), ...
    UP_DOWN_info.up_duration(:), ...
    UP_DOWN_info.next_DOWN_duration(:), ...
    UP_DOWN_info.previous_DOWN_duration(:), ...
    UP_DOWN_info.SWpeakmag_UD(:), ...
    UP_DOWN_info.SWpeakmag_DU(:), ...
    UD_lags, DU_lags, UD_lags_signed, DU_lags_signed, ...
    lastRippleNormalisedUP(:), ...
    firstRippleNormalisedUP(:), ...
    UP_DOWN_info.ripple_count(:), ...
    merged_event_info.subject_id(:), ...
    merged_event_info.session_id(:), ...
    'VariableNames', { ...
    'TimeToFirstRipple', 'TimefromFirstRipple', 'TimeToLastRipple', 'TimefromLastRipple', 'TimetoNextUP', 'LastRippleDuration'...
    'V1Track', 'firstRipplePower', 'lastRipplePower', ...
    'nextUPV1', 'earlyUPV1', 'lateUPV1', 'previousUPV1', ...
    'firstRippleV1', 'firstRippleV1PRE', 'lastRippleV1', 'lastRippleV1PRE', ...
    'nextUPHPC', 'earlyUPHPC', 'lateUPHPC', 'previousUPHPC', ...
    'firstRippleHPC', 'firstRippleHPCPRE', 'lastRippleHPC', 'lastRippleHPCPRE', ...
    'UPDuration', 'nextDOWNDuration', 'previousDOWNDuration', 'nextDOWNSOPower', 'previousDOWNSOPower', ...
    'nextDOWNlag', 'UPlag', 'nextDOWNlagSigned', 'UPlagSigned', ...
    'lastRippleNormalisedUP', 'firstRippleNormalisedUP', ...
    'RippleCounts', 'AnimalID', 'SessionID'});

tbl.secondUPcumIpsiV1 = T.second_half_V1_MUA_sum(:);
tbl.firstUPcumIpsiV1  = T.first_half_V1_MUA_sum(:);
tbl.UPcumDiffIpsiV1   = T.second_half_V1_MUA_sum(:) - T.first_half_V1_MUA_sum(:);

tbl.secondUPcumHPC = T.second_half_HPC_MUA_sum(:);
tbl.firstUPcumHPC  = T.first_half_HPC_MUA_sum(:);
tbl.UPcumDiffHPC   = T.second_half_HPC_MUA_sum(:) - T.first_half_HPC_MUA_sum(:);
tbl.normalisedRippleUPcumHPC = (T.first_half_ripple_HPC_MUA_sum(:) + T.second_half_ripple_HPC_MUA_sum(:)) ./ ...
    (T.first_half_ripple_duration(:) + T.second_half_ripple_duration(:));

tbl.lastRippleMUArate=UP_DOWN_info.last_ripple_HPC_MUA_mean;
tbl.firstRippleMUArate=UP_DOWN_info.first_ripple_HPC_MUA_mean;
tbl.lastRippleMUArateV1=UP_DOWN_info.last_ripple_V1_MUA_mean;
tbl.firstRippleMUArateV1=UP_DOWN_info.first_ripple_V1_MUA_mean;

% writetable(tbl, 'UP_DOWN_info_GAM_peak2.csv');

% writetable(tbl, 'UP_DOWN_info_GAM_peak.csv');
writetable(tbl, 'UP_DOWN_info_GAM_offset.csv');


tbl = readtable('UP_DOWN_info_GAM_offset.csv');

tbl = readtable('UP_DOWN_info_GAM_peak.csv');


tbl1 = readtable('UP_DOWN_info_GAM.csv');

histogram(tbl.TimefromLastRipple);hold on;histogram(tbl1.TimefromLastRipple)
histogram(tbl.lastRipplePower);hold on;histogram(tbl1.lastRipplePower)


x = tbl1.TimefromLastRipple(tbl1.TimefromLastRipple>=0);

histogram(log10(x))

;hold on;histogram(tbl1.TimefromLastRipple)

scatter(tbl.TimefromLastRipple,tbl1.TimefromLastRipple)


fprintf('UP_DOWN_info and GAM table built for %d UP events.\n', height(tbl));

%%


x = tbl.lastRippleMUArate;
y = tbl.lastRipplePower;

plot_2D_scatter(x, y, ...
    'title_name','Last ripple MUA rate vs ripple power',...
    'xlabel_text', 'Last Ripple HPC normalised MUA rate', ...
    'ylabel_text', 'Last Ripple power (z)', ...
    'animal_id', tbl.AnimalID, ...       % Vector or cell array of animal IDs
    'session_id', tbl.SessionID, ...     % Vector or cell array of session IDs
    'use_mixed_effects', true, ...     % Enable LME model for stats
    'title_text', 'Density Scatter', ...
    'show_stats', true, ...
    'add_fit_line', true, ...
    'nbins',200,...
    'smooth_color', true,...
    'smooth_sigma', 3,...
    'xscale', 'log10',...
    'yscale', 'log10')

xlim([log10(0.2) log10(1)])
ylim([log10(5) log10(30)])
save_all_figures(output_dir, [],'SVG_option',1);

% scatter(tbl.lastRippleMUArate,tbl.lastRipplePower)