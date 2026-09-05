function T = build_last_ripple_to_DOWN_summary_table(merged_event_info, V1_MUA_spiketimes, HC_MUA_spiketimes, varargin)
% BUILD_LAST_RIPPLE_TO_DOWN_SUMMARY_TABLE Builds a per-UP event summary table
% for analyzing survival time from the last ripple to UP state termination.
%
% Inputs:
%   merged_event_info - Struct containing:
%       .UP_ints              - Nx2 matrix [onset, offset] of UP states (abs time)
%       .UP_hemisphere_id     - Nx1 vector (1=Left V1, 2=Right V1)
%       .ripples_ints         - Mx2 matrix [onset, offset] of ripples (abs time)
%       .ripples_peaktimes    - Mx1 vector of ripple peak times (abs time)
%       .ripples_power        - Mx1 vector of LFP ripple z-score peak powers
%       .session_id           - Nx1 vector of session IDs
%       .subject_id           - (Optional) Nx1 vector of subject IDs
%   V1_MUA_spiketimes - Cell array {1} (Left V1) and {2} (Right V1) of spike times
%   HC_MUA_spiketimes - Cell array {1} (Left HPC) and {2} (Right HPC) of spike times
%
% Options:
%   'time_reference'   - 'offset' (default) or 'peak' for time from last ripple to UP term
%   'non_ripple_scope' - 'entire_UP' (default) or 'prior_to_last_ripple'
%   'time_bin'         - 0.01 (10ms bin size, default)
%   'min_interval_dur' - 1e-5 minimum duration threshold

p = inputParser;
addParameter(p, 'time_reference', 'offset', @(x) ismember(x, {'offset', 'peak'}));
addParameter(p, 'non_ripple_scope', 'entire_UP', @(x) ismember(x, {'entire_UP', 'prior_to_last_ripple'}));
addParameter(p, 'time_bin', 0.01, @isnumeric);
addParameter(p, 'min_interval_dur', 1e-5, @isnumeric);
addParameter(p, 'tau_recency', 0.5, @isnumeric);
parse(p, varargin{:});

time_ref    = p.Results.time_reference;
nr_scope    = p.Results.non_ripple_scope;
time_bin    = p.Results.time_bin;
min_dur     = p.Results.min_interval_dur;
tau_recency = p.Results.tau_recency;

nUP = size(merged_event_info.UP_ints, 1);
ripIntsAbs  = merged_event_info.ripples_ints;
ripPeakAbs  = merged_event_info.ripples_peaktimes;
ripPowerAbs = merged_event_info.ripples_power;

% Session and subject metadata assignment
if isfield(merged_event_info, 'session_id') && ~isempty(merged_event_info.session_id)
    session_id = merged_event_info.session_id;
else
    session_id = floor(merged_event_info.UP_ints(:,1) / 1000000);
end

if isfield(merged_event_info, 'subject_id') && ~isempty(merged_event_info.subject_id)
    subject_id = merged_event_info.subject_id;
else
    subject_id = session_id;
end

%% 1. Pre-compute 10ms binned session-normalized MUA directly from raw spike times
unique_sessions = unique(session_id);
max_sess_id = max(unique_sessions);

norm_v1_L = cell(max_sess_id, 1);
norm_v1_R = cell(max_sess_id, 1);
norm_hc_L = cell(max_sess_id, 1);
norm_hc_R = cell(max_sess_id, 1);

fprintf('Pre-computing session-normalized 10ms MUA directly from spike times...\n');
for sIdx = 1:length(unique_sessions)
    s = unique_sessions(sIdx);
    s_offset = s * 1000000;
    s_next   = (s + 1) * 1000000;
    
    % Filter spike times for this session
    spk_v1_L = V1_MUA_spiketimes{1}(V1_MUA_spiketimes{1} >= s_offset & V1_MUA_spiketimes{1} < s_next) - s_offset;
    spk_v1_R = V1_MUA_spiketimes{2}(V1_MUA_spiketimes{2} >= s_offset & V1_MUA_spiketimes{2} < s_next) - s_offset;
    spk_hc_L = HC_MUA_spiketimes{1}(HC_MUA_spiketimes{1} >= s_offset & HC_MUA_spiketimes{1} < s_next) - s_offset;
    spk_hc_R = HC_MUA_spiketimes{2}(HC_MUA_spiketimes{2} >= s_offset & HC_MUA_spiketimes{2} < s_next) - s_offset;
    
    max_up_t = max(merged_event_info.UP_ints(session_id == s, 2) - s_offset);
    max_spk_t = max([0; spk_v1_L; spk_v1_R; spk_hc_L; spk_hc_R]);
    max_t = max(max_up_t, max_spk_t) + time_bin;
    
    t_edges = 0:time_bin:max_t;
    
    cnt_v1_L = histcounts(spk_v1_L, t_edges);
    cnt_v1_R = histcounts(spk_v1_R, t_edges);
    cnt_hc_L = histcounts(spk_hc_L, t_edges);
    cnt_hc_R = histcounts(spk_hc_R, t_edges);
    
    min_v1_L = min(cnt_v1_L); p99_v1_L = prctile(cnt_v1_L, 99); if p99_v1_L <= min_v1_L, p99_v1_L = min_v1_L + 1; end
    min_v1_R = min(cnt_v1_R); p99_v1_R = prctile(cnt_v1_R, 99); if p99_v1_R <= min_v1_R, p99_v1_R = min_v1_R + 1; end
    min_hc_L = min(cnt_hc_L); p99_hc_L = prctile(cnt_hc_L, 99); if p99_hc_L <= min_hc_L, p99_hc_L = min_hc_L + 1; end
    min_hc_R = min(cnt_hc_R); p99_hc_R = prctile(cnt_hc_R, 99); if p99_hc_R <= min_hc_R, p99_hc_R = min_hc_R + 1; end
    
    norm_v1_L{s} = min(1, max(0, (cnt_v1_L - min_v1_L) / (p99_v1_L - min_v1_L)));
    norm_v1_R{s} = min(1, max(0, (cnt_v1_R - min_v1_R) / (p99_v1_R - min_v1_R)));
    norm_hc_L{s} = min(1, max(0, (cnt_hc_L - min_hc_L) / (p99_hc_L - min_hc_L)));
    norm_hc_R{s} = min(1, max(0, (cnt_hc_R - min_hc_R) / (p99_hc_R - min_hc_R)));
end

%% 2. Pre-allocate table columns
upID_col                   = (1:nUP)';
session_id_col             = session_id;
subject_id_col             = subject_id;
hemisphere_id_col          = merged_event_info.UP_hemisphere_id;
up_duration_col            = zeros(nUP, 1);
ripple_count_col           = zeros(nUP, 1);
last_ripple_to_UP_term_col = nan(nUP, 1);

% Factor 1: Last Ripple Features (Actual unclipped ripple bounds)
last_ripple_duration_col     = nan(nUP, 1);
last_ripple_power_col        = nan(nUP, 1);
last_ripple_interval_col     = nan(nUP, 1);
last_ripple_HPC_MUA_sum_col  = nan(nUP, 1);
last_ripple_HPC_MUA_mean_col = nan(nUP, 1);
last_ripple_V1_MUA_sum_col   = nan(nUP, 1);
last_ripple_V1_MUA_mean_col  = nan(nUP, 1);

% Last Ripple Peak MUA & Half-Duration MUA metrics
last_ripple_HPC_MUA_peak_col                 = nan(nUP, 1);
last_ripple_V1_MUA_peak_col                  = nan(nUP, 1);
last_ripple_first_half_HPC_MUA_sum_col       = nan(nUP, 1);
last_ripple_first_half_HPC_MUA_mean_col      = nan(nUP, 1);
last_ripple_second_half_HPC_MUA_sum_col      = nan(nUP, 1);
last_ripple_second_half_HPC_MUA_mean_col     = nan(nUP, 1);
last_ripple_first_half_V1_MUA_sum_col        = nan(nUP, 1);
last_ripple_first_half_V1_MUA_mean_col       = nan(nUP, 1);
last_ripple_second_half_V1_MUA_sum_col       = nan(nUP, 1);
last_ripple_second_half_V1_MUA_mean_col      = nan(nUP, 1);

% Factor 2: Past Ripples History (Overall)
past_ripples_count_col        = nan(nUP, 1);
past_ripples_duration_col     = nan(nUP, 1);
past_ripples_HPC_MUA_sum_col  = nan(nUP, 1);
past_ripples_HPC_MUA_mean_col = nan(nUP, 1);
past_ripples_V1_MUA_sum_col   = nan(nUP, 1);
past_ripples_V1_MUA_mean_col  = nan(nUP, 1);

% Past Ripples in 1st Half & 2nd Half (EXCLUDING last ripple)
first_half_past_ripples_HPC_MUA_sum_col   = zeros(nUP, 1);
first_half_past_ripples_HPC_MUA_mean_col  = zeros(nUP, 1);
second_half_past_ripples_HPC_MUA_sum_col  = zeros(nUP, 1);
second_half_past_ripples_HPC_MUA_mean_col = zeros(nUP, 1);

first_half_past_ripples_V1_MUA_sum_col    = zeros(nUP, 1);
first_half_past_ripples_V1_MUA_mean_col   = zeros(nUP, 1);
second_half_past_ripples_V1_MUA_sum_col    = zeros(nUP, 1);
second_half_past_ripples_V1_MUA_mean_col   = zeros(nUP, 1);

first_half_past_ripples_duration_col     = zeros(nUP, 1);
second_half_past_ripples_duration_col    = zeros(nUP, 1);

% Recency-Weighted Past Ripples MUA (tau = 0.1s)
recency_weighted_past_HPC_MUA_sum_col  = zeros(nUP, 1);
recency_weighted_past_HPC_MUA_mean_col = zeros(nUP, 1);
recency_weighted_past_V1_MUA_sum_col   = zeros(nUP, 1);
recency_weighted_past_V1_MUA_mean_col  = zeros(nUP, 1);

% Binned Past Ripples MUA by distance to last ripple peak (0-100ms, 100-200ms, 200+ms)
past_ripples_0_100ms_HPC_MUA_sum_col    = zeros(nUP, 1);
past_ripples_0_100ms_HPC_MUA_mean_col   = zeros(nUP, 1);
past_ripples_0_100ms_V1_MUA_sum_col     = zeros(nUP, 1);
past_ripples_0_100ms_V1_MUA_mean_col    = zeros(nUP, 1);

past_ripples_100_200ms_HPC_MUA_sum_col  = zeros(nUP, 1);
past_ripples_100_200ms_HPC_MUA_mean_col = zeros(nUP, 1);
past_ripples_100_200ms_V1_MUA_sum_col   = zeros(nUP, 1);
past_ripples_100_200ms_V1_MUA_mean_col  = zeros(nUP, 1);

past_ripples_200plus_ms_HPC_MUA_sum_col  = zeros(nUP, 1);
past_ripples_200plus_ms_HPC_MUA_mean_col = zeros(nUP, 1);
past_ripples_200plus_ms_V1_MUA_sum_col   = zeros(nUP, 1);
past_ripples_200plus_ms_V1_MUA_mean_col  = zeros(nUP, 1);

% Binned Non-Ripple MUA prior to last ripple peak by distance (0-100ms, 100-200ms, 200+ms)
non_ripple_0_100ms_HPC_MUA_sum_col    = zeros(nUP, 1);
non_ripple_0_100ms_HPC_MUA_mean_col   = zeros(nUP, 1);
non_ripple_0_100ms_V1_MUA_sum_col     = zeros(nUP, 1);
non_ripple_0_100ms_V1_MUA_mean_col    = zeros(nUP, 1);

non_ripple_100_200ms_HPC_MUA_sum_col  = zeros(nUP, 1);
non_ripple_100_200ms_HPC_MUA_mean_col = zeros(nUP, 1);
non_ripple_100_200ms_V1_MUA_sum_col   = zeros(nUP, 1);
non_ripple_100_200ms_V1_MUA_mean_col  = zeros(nUP, 1);

non_ripple_200plus_ms_HPC_MUA_sum_col  = zeros(nUP, 1);
non_ripple_200plus_ms_HPC_MUA_mean_col = zeros(nUP, 1);
non_ripple_200plus_ms_V1_MUA_sum_col   = zeros(nUP, 1);
non_ripple_200plus_ms_V1_MUA_mean_col  = zeros(nUP, 1);

% Factor 3: Non-Ripple MUA
non_ripple_duration_col     = zeros(nUP, 1);
non_ripple_HPC_MUA_sum_col  = zeros(nUP, 1);
non_ripple_HPC_MUA_mean_col = zeros(nUP, 1);
non_ripple_V1_MUA_sum_col   = zeros(nUP, 1);
non_ripple_V1_MUA_mean_col  = zeros(nUP, 1);

% Overall 1st half & 2nd half UP state MUA metrics
first_half_HPC_MUA_sum_col   = zeros(nUP, 1);
first_half_HPC_MUA_mean_col  = zeros(nUP, 1);
second_half_HPC_MUA_sum_col  = zeros(nUP, 1);
second_half_HPC_MUA_mean_col = zeros(nUP, 1);

first_half_V1_MUA_sum_col    = zeros(nUP, 1);
first_half_V1_MUA_mean_col   = zeros(nUP, 1);
second_half_V1_MUA_sum_col   = zeros(nUP, 1);
second_half_V1_MUA_mean_col  = zeros(nUP, 1);

% 1st half & 2nd half ALL RIPPLE MUA metrics (including last ripple)
first_half_ripple_HPC_MUA_sum_col   = zeros(nUP, 1);
first_half_ripple_HPC_MUA_mean_col  = zeros(nUP, 1);
second_half_ripple_HPC_MUA_sum_col  = zeros(nUP, 1);
second_half_ripple_HPC_MUA_mean_col = zeros(nUP, 1);

first_half_ripple_V1_MUA_sum_col    = zeros(nUP, 1);
first_half_ripple_V1_MUA_mean_col   = zeros(nUP, 1);
second_half_ripple_V1_MUA_sum_col    = zeros(nUP, 1);
second_half_ripple_V1_MUA_mean_col   = zeros(nUP, 1);

% 1st half & 2nd half NON-RIPPLE MUA metrics
first_half_non_ripple_HPC_MUA_sum_col   = zeros(nUP, 1);
first_half_non_ripple_HPC_MUA_mean_col  = zeros(nUP, 1);
second_half_non_ripple_HPC_MUA_sum_col  = zeros(nUP, 1);
second_half_non_ripple_HPC_MUA_mean_col = zeros(nUP, 1);

first_half_non_ripple_V1_MUA_sum_col    = zeros(nUP, 1);
first_half_non_ripple_V1_MUA_mean_col   = zeros(nUP, 1);
second_half_non_ripple_V1_MUA_sum_col    = zeros(nUP, 1);
second_half_non_ripple_V1_MUA_mean_col   = zeros(nUP, 1);

% Cumulative Ripple & Non-Ripple Duration in 1st & 2nd half of UP
first_half_ripple_duration_col     = zeros(nUP, 1);
second_half_ripple_duration_col    = zeros(nUP, 1);
first_half_non_ripple_duration_col  = zeros(nUP, 1);
second_half_non_ripple_duration_col = zeros(nUP, 1);

%% 3. Iterate over each UP event
fprintf('Extracting per-UP summary features for last ripple to DOWN survival...\n');
for iUP = 1:nUP
    upOnset  = merged_event_info.UP_ints(iUP, 1);
    upOffset = merged_event_info.UP_ints(iUP, 2);
    upDur    = upOffset - upOnset;
    hemiID   = merged_event_info.UP_hemisphere_id(iUP);
    sessID   = session_id(iUP);
    
    up_duration_col(iUP) = upDur;
    
    s_offset   = sessID * 1000000;
    rel_onset  = upOnset - s_offset;
    rel_offset = upOffset - s_offset;
    
    bin_start = max(1, floor(rel_onset / time_bin) + 1);
    bin_end   = min(length(norm_v1_L{sessID}), ceil(rel_offset / time_bin));
    if bin_start > bin_end, bin_end = bin_start; end
    
    % Mean HPC MUA across probes
    hpcMUA = mean([norm_hc_L{sessID}(bin_start:bin_end); norm_hc_R{sessID}(bin_start:bin_end)], 1, 'omitnan');
    
    % Ipsilateral V1 MUA based on UP event hemisphere_id (1=Left V1, 2=Right V1)
    if hemiID == 1
        v1MUA = norm_v1_L{sessID}(bin_start:bin_end);
    else
        v1MUA = norm_v1_R{sessID}(bin_start:bin_end);
    end
    
    hpcMUA = hpcMUA(:)';
    v1MUA  = v1MUA(:)';
    nBins  = length(hpcMUA);
    binCenters = (0:nBins-1) * time_bin + time_bin/2;
    
    % First Half and Second Half UP state MUA metrics (Overall)
    half_idx = floor(nBins / 2);
    if half_idx < 1, half_idx = 1; end
    
    first_half_HPC_MUA_sum_col(iUP)  = sum(hpcMUA(1:half_idx));
    first_half_HPC_MUA_mean_col(iUP) = mean(hpcMUA(1:half_idx));
    first_half_V1_MUA_sum_col(iUP)   = sum(v1MUA(1:half_idx));
    first_half_V1_MUA_mean_col(iUP)  = mean(v1MUA(1:half_idx));
    
    if nBins > 1
        second_half_HPC_MUA_sum_col(iUP)  = sum(hpcMUA(half_idx+1:end));
        second_half_HPC_MUA_mean_col(iUP) = mean(hpcMUA(half_idx+1:end));
        second_half_V1_MUA_sum_col(iUP)   = sum(v1MUA(half_idx+1:end));
        second_half_V1_MUA_mean_col(iUP)  = mean(v1MUA(half_idx+1:end));
    else
        second_half_HPC_MUA_sum_col(iUP)  = first_half_HPC_MUA_sum_col(iUP);
        second_half_HPC_MUA_mean_col(iUP) = first_half_HPC_MUA_mean_col(iUP);
        second_half_V1_MUA_sum_col(iUP)   = first_half_V1_MUA_sum_col(iUP);
        second_half_V1_MUA_mean_col(iUP)  = first_half_V1_MUA_mean_col(iUP);
    end
    
    % Find ripples overlapping this UP event
    overlapIdx = find(ripIntsAbs(:,2) > upOnset & ripIntsAbs(:,1) < upOffset);
    
    ripAbsOn  = [];
    ripAbsOff = [];
    ripAbsPk  = [];
    ripPow    = [];
    
    if ~isempty(overlapIdx)
        for r = 1:length(overlapIdx)
            idx = overlapIdx(r);
            rOn_abs  = ripIntsAbs(idx, 1);
            rOff_abs = ripIntsAbs(idx, 2);
            rPk_abs  = ripPeakAbs(idx);
            
            if rOff_abs - rOn_abs >= min_dur
                ripAbsOn  = [ripAbsOn; rOn_abs];
                ripAbsOff = [ripAbsOff; rOff_abs];
                ripAbsPk  = [ripAbsPk; rPk_abs];
                ripPow    = [ripPow; ripPowerAbs(idx)];
            end
        end
    end
    
    % Sort ripples by absolute onset time
    if ~isempty(ripAbsOn)
        [ripAbsOn, sIdx] = sort(ripAbsOn);
        ripAbsOff = ripAbsOff(sIdx);
        ripAbsPk  = ripAbsPk(sIdx);
        ripPow    = ripPow(sIdx);
    end
    
    k = length(ripAbsOn);
    ripple_count_col(iUP) = k;
    
    % Calculate cumulative ripple durations in 1st half and 2nd half of UP
    half_t = upDur / 2;
    dur_1st_rip = 0;
    dur_2nd_rip = 0;
    for rIdx = 1:k
        on_rel  = ripAbsOn(rIdx) - upOnset;
        off_rel = ripAbsOff(rIdx) - upOnset;
        
        ov_1st = max(0, min(half_t, off_rel) - max(0, on_rel));
        if ov_1st > 0
            dur_1st_rip = dur_1st_rip + ov_1st;
        end
        
        ov_2nd = max(0, min(upDur, off_rel) - max(half_t, max(0, on_rel)));
        if ov_2nd > 0
            dur_2nd_rip = dur_2nd_rip + ov_2nd;
        end
    end
    
    first_half_ripple_duration_col(iUP)     = dur_1st_rip;
    second_half_ripple_duration_col(iUP)    = dur_2nd_rip;
    first_half_non_ripple_duration_col(iUP)  = max(0, half_t - dur_1st_rip);
    second_half_non_ripple_duration_col(iUP) = max(0, (upDur - half_t) - dur_2nd_rip);
    
    % Create bin mask for ripples inside UP state
    isRipBin = false(1, nBins);
    for rIdx = 1:k
        rOn_rel  = ripAbsOn(rIdx) - upOnset;
        rOff_rel = ripAbsOff(rIdx) - upOnset;
        isRipBin = isRipBin | (binCenters >= rOn_rel & binCenters < rOff_rel);
    end
    
    % Split masks into 1st half and 2nd half ripple & non-ripple bins
    idx_mask = 1:nBins;
    mask_1st_rip    = (idx_mask <= half_idx) & isRipBin;
    mask_1st_nonrip = (idx_mask <= half_idx) & ~isRipBin;
    
    if nBins > 1
        mask_2nd_rip    = (idx_mask > half_idx) & isRipBin;
        mask_2nd_nonrip = (idx_mask > half_idx) & ~isRipBin;
    else
        mask_2nd_rip    = false(1, nBins);
        mask_2nd_nonrip = false(1, nBins);
    end
    
    % 1st Half Ripple MUA (All ripples)
    if any(mask_1st_rip)
        first_half_ripple_HPC_MUA_sum_col(iUP)  = sum(hpcMUA(mask_1st_rip));
        first_half_ripple_HPC_MUA_mean_col(iUP) = mean(hpcMUA(mask_1st_rip));
        first_half_ripple_V1_MUA_sum_col(iUP)   = sum(v1MUA(mask_1st_rip));
        first_half_ripple_V1_MUA_mean_col(iUP)  = mean(v1MUA(mask_1st_rip));
    else
        first_half_ripple_HPC_MUA_sum_col(iUP)  = 0;
        first_half_ripple_HPC_MUA_mean_col(iUP) = 0;
        first_half_ripple_V1_MUA_sum_col(iUP)   = 0;
        first_half_ripple_V1_MUA_mean_col(iUP)  = 0;
    end
    
    % 2nd Half Ripple MUA (All ripples)
    if any(mask_2nd_rip)
        second_half_ripple_HPC_MUA_sum_col(iUP)  = sum(hpcMUA(mask_2nd_rip));
        second_half_ripple_HPC_MUA_mean_col(iUP) = mean(hpcMUA(mask_2nd_rip));
        second_half_ripple_V1_MUA_sum_col(iUP)   = sum(v1MUA(mask_2nd_rip));
        second_half_ripple_V1_MUA_mean_col(iUP)  = mean(v1MUA(mask_2nd_rip));
    else
        second_half_ripple_HPC_MUA_sum_col(iUP)  = 0;
        second_half_ripple_HPC_MUA_mean_col(iUP) = 0;
        second_half_ripple_V1_MUA_sum_col(iUP)   = 0;
        second_half_ripple_V1_MUA_mean_col(iUP)  = 0;
    end
    
    % 1st Half Non-Ripple MUA
    if any(mask_1st_nonrip)
        first_half_non_ripple_HPC_MUA_sum_col(iUP)  = sum(hpcMUA(mask_1st_nonrip));
        first_half_non_ripple_HPC_MUA_mean_col(iUP) = mean(hpcMUA(mask_1st_nonrip));
        first_half_non_ripple_V1_MUA_sum_col(iUP)   = sum(v1MUA(mask_1st_nonrip));
        first_half_non_ripple_V1_MUA_mean_col(iUP)  = mean(v1MUA(mask_1st_nonrip));
    else
        first_half_non_ripple_HPC_MUA_sum_col(iUP)  = 0;
        first_half_non_ripple_HPC_MUA_mean_col(iUP) = 0;
        first_half_non_ripple_V1_MUA_sum_col(iUP)   = 0;
        first_half_non_ripple_V1_MUA_mean_col(iUP)  = 0;
    end
    
    % 2nd Half Non-Ripple MUA
    if any(mask_2nd_nonrip)
        second_half_non_ripple_HPC_MUA_sum_col(iUP)  = sum(hpcMUA(mask_2nd_nonrip));
        second_half_non_ripple_HPC_MUA_mean_col(iUP) = mean(hpcMUA(mask_2nd_nonrip));
        second_half_non_ripple_V1_MUA_sum_col(iUP)   = sum(v1MUA(mask_2nd_nonrip));
        second_half_non_ripple_V1_MUA_mean_col(iUP)  = mean(v1MUA(mask_2nd_nonrip));
    else
        second_half_non_ripple_HPC_MUA_sum_col(iUP)  = 0;
        second_half_non_ripple_HPC_MUA_mean_col(iUP) = 0;
        second_half_non_ripple_V1_MUA_sum_col(iUP)   = 0;
        second_half_non_ripple_V1_MUA_mean_col(iUP)  = 0;
    end
    
    if k == 0
        % No ripples in this UP event
        non_ripple_duration_col(iUP)     = upDur;
        non_ripple_HPC_MUA_sum_col(iUP)  = sum(hpcMUA);
        non_ripple_HPC_MUA_mean_col(iUP) = mean(hpcMUA);
        non_ripple_V1_MUA_sum_col(iUP)   = sum(v1MUA);
        non_ripple_V1_MUA_mean_col(iUP)  = mean(v1MUA);
    else
        % At least 1 ripple present
        last_abs_on  = ripAbsOn(end);
        last_abs_off = ripAbsOff(end);
        last_abs_pk  = ripAbsPk(end);
        
        % Use ACTUAL full unclipped ripple duration
        last_ripple_duration_col(iUP) = last_abs_off - last_abs_on;
        last_ripple_power_col(iUP)    = ripPow(end);
        if k > 1
            last_ripple_interval_col(iUP) = ripAbsPk(end) - ripAbsPk(end-1);
        end
        
        if strcmp(time_ref, 'peak')
            ref_t = last_abs_pk - upOnset;
        else
            ref_t = last_abs_off - upOnset;
        end
        
        last_ripple_to_UP_term_col(iUP) = upDur - ref_t;
        
        % Extract MUA across actual full ripple bounds (relative to session offset)
        rel_rip_on  = last_abs_on - s_offset;
        rel_rip_off = last_abs_off - s_offset;
        
        bStart_rip = max(1, floor(rel_rip_on / time_bin) + 1);
        bEnd_rip   = min(length(norm_v1_L{sessID}), ceil(rel_rip_off / time_bin));
        if bStart_rip > bEnd_rip, bEnd_rip = bStart_rip; end
        
        lastHpcMUA = mean([norm_hc_L{sessID}(bStart_rip:bEnd_rip); norm_hc_R{sessID}(bStart_rip:bEnd_rip)], 1, 'omitnan');
        if hemiID == 1
            lastV1MUA = norm_v1_L{sessID}(bStart_rip:bEnd_rip);
        else
            lastV1MUA = norm_v1_R{sessID}(bStart_rip:bEnd_rip);
        end
        
        last_ripple_HPC_MUA_sum_col(iUP)  = sum(lastHpcMUA);
        last_ripple_HPC_MUA_mean_col(iUP) = mean(lastHpcMUA);
        last_ripple_V1_MUA_sum_col(iUP)   = sum(lastV1MUA);
        last_ripple_V1_MUA_mean_col(iUP)  = mean(lastV1MUA);
        
        % Peak MUA during last ripple (across actual ripple bounds)
        last_ripple_HPC_MUA_peak_col(iUP) = max(lastHpcMUA);
        last_ripple_V1_MUA_peak_col(iUP)  = max(lastV1MUA);
        
        % 1st Half vs 2nd Half of Last Ripple Duration MUA
        nRipBins = length(lastHpcMUA);
        rHalf    = floor(nRipBins / 2);
        if rHalf < 1, rHalf = 1; end
        
        last_ripple_first_half_HPC_MUA_sum_col(iUP)  = sum(lastHpcMUA(1:rHalf));
        last_ripple_first_half_HPC_MUA_mean_col(iUP) = mean(lastHpcMUA(1:rHalf));
        last_ripple_first_half_V1_MUA_sum_col(iUP)   = sum(lastV1MUA(1:rHalf));
        last_ripple_first_half_V1_MUA_mean_col(iUP)  = mean(lastV1MUA(1:rHalf));
        
        if nRipBins > 1
            last_ripple_second_half_HPC_MUA_sum_col(iUP)  = sum(lastHpcMUA(rHalf+1:end));
            last_ripple_second_half_HPC_MUA_mean_col(iUP) = mean(lastHpcMUA(rHalf+1:end));
            last_ripple_second_half_V1_MUA_sum_col(iUP)   = sum(lastV1MUA(rHalf+1:end));
            last_ripple_second_half_V1_MUA_mean_col(iUP)  = mean(lastV1MUA(rHalf+1:end));
        else
            last_ripple_second_half_HPC_MUA_sum_col(iUP)  = last_ripple_first_half_HPC_MUA_sum_col(iUP);
            last_ripple_second_half_HPC_MUA_mean_col(iUP) = last_ripple_first_half_HPC_MUA_mean_col(iUP);
            last_ripple_second_half_V1_MUA_sum_col(iUP)   = last_ripple_first_half_V1_MUA_sum_col(iUP);
            last_ripple_second_half_V1_MUA_mean_col(iUP)  = last_ripple_first_half_V1_MUA_mean_col(iUP);
        end
        
        % Factor 2: Past Ripples History (EXCLUDING last ripple)
        past_ripples_count_col(iUP) = k - 1;
        if k == 1
            past_ripples_duration_col(iUP)     = 0;
            past_ripples_HPC_MUA_sum_col(iUP)  = 0;
            past_ripples_HPC_MUA_mean_col(iUP) = 0;
            past_ripples_V1_MUA_sum_col(iUP)   = 0;
            past_ripples_V1_MUA_mean_col(iUP)  = 0;
        else
            pastAbsOn  = ripAbsOn(1:end-1);
            pastAbsOff = ripAbsOff(1:end-1);
            past_ripples_duration_col(iUP) = sum(pastAbsOff - pastAbsOn);
            
            pastBinMask = false(1, nBins);
            for pIdx = 1:length(pastAbsOn)
                pOn_rel  = pastAbsOn(pIdx) - upOnset;
                pOff_rel = pastAbsOff(pIdx) - upOnset;
                pastBinMask = pastBinMask | (binCenters >= pOn_rel & binCenters < pOff_rel);
            end
            
            if any(pastBinMask)
                past_ripples_HPC_MUA_sum_col(iUP)  = sum(hpcMUA(pastBinMask));
                past_ripples_HPC_MUA_mean_col(iUP) = mean(hpcMUA(pastBinMask));
                past_ripples_V1_MUA_sum_col(iUP)   = sum(v1MUA(pastBinMask));
                past_ripples_V1_MUA_mean_col(iUP)  = mean(v1MUA(pastBinMask));
            else
                past_ripples_HPC_MUA_sum_col(iUP)  = 0;
                past_ripples_HPC_MUA_mean_col(iUP) = 0;
                past_ripples_V1_MUA_sum_col(iUP)   = 0;
                past_ripples_V1_MUA_mean_col(iUP)  = 0;
            end
            
            % Past ripples MUA in 1st half & 2nd half of UP (EXCLUDING last ripple)
            dur_1st_past = 0;
            dur_2nd_past = 0;
            for pIdx = 1:length(pastAbsOn)
                on_rel  = pastAbsOn(pIdx) - upOnset;
                off_rel = pastAbsOff(pIdx) - upOnset;
                
                ov_1st = max(0, min(half_t, off_rel) - max(0, on_rel));
                if ov_1st > 0, dur_1st_past = dur_1st_past + ov_1st; end
                
                ov_2nd = max(0, min(upDur, off_rel) - max(half_t, max(0, on_rel)));
                if ov_2nd > 0, dur_2nd_past = dur_2nd_past + ov_2nd; end
            end
            first_half_past_ripples_duration_col(iUP)  = dur_1st_past;
            second_half_past_ripples_duration_col(iUP) = dur_2nd_past;
            
            mask_1st_past = (idx_mask <= half_idx) & pastBinMask;
            if nBins > 1
                mask_2nd_past = (idx_mask > half_idx) & pastBinMask;
            else
                mask_2nd_past = false(1, nBins);
            end
            
            if any(mask_1st_past)
                first_half_past_ripples_HPC_MUA_sum_col(iUP)  = sum(hpcMUA(mask_1st_past));
                first_half_past_ripples_HPC_MUA_mean_col(iUP) = mean(hpcMUA(mask_1st_past));
                first_half_past_ripples_V1_MUA_sum_col(iUP)   = sum(v1MUA(mask_1st_past));
                first_half_past_ripples_V1_MUA_mean_col(iUP)  = mean(v1MUA(mask_1st_past));
            end
            
            if any(mask_2nd_past)
                second_half_past_ripples_HPC_MUA_sum_col(iUP)  = sum(hpcMUA(mask_2nd_past));
                second_half_past_ripples_HPC_MUA_mean_col(iUP) = mean(hpcMUA(mask_2nd_past));
                second_half_past_ripples_V1_MUA_sum_col(iUP)   = sum(v1MUA(mask_2nd_past));
                second_half_past_ripples_V1_MUA_mean_col(iUP)  = mean(v1MUA(mask_2nd_past));
            end
            
            % Recency-weighted past ripples MUA
            pastAbsPk = ripAbsPk(1:end-1);
            n_past    = length(pastAbsOn);
            rec_sum_hpc  = 0; rec_mean_hpc = 0;
            rec_sum_v1   = 0; rec_mean_v1  = 0;
            for pIdx = 1:n_past
                dist_to_last = ripAbsPk(end) - pastAbsPk(pIdx);
                w = exp(-dist_to_last / tau_recency);
                
                rel_p_on  = pastAbsOn(pIdx) - s_offset;
                rel_p_off = pastAbsOff(pIdx) - s_offset;
                bStart_p = max(1, floor(rel_p_on / time_bin) + 1);
                bEnd_p   = min(length(norm_v1_L{sessID}), ceil(rel_p_off / time_bin));
                if bStart_p > bEnd_p, bEnd_p = bStart_p; end
                
                pHpcMUA = mean([norm_hc_L{sessID}(bStart_p:bEnd_p); norm_hc_R{sessID}(bStart_p:bEnd_p)], 1, 'omitnan');
                if hemiID == 1
                    pV1MUA = norm_v1_L{sessID}(bStart_p:bEnd_p);
                else
                    pV1MUA = norm_v1_R{sessID}(bStart_p:bEnd_p);
                end
                
                rec_sum_hpc  = rec_sum_hpc  + sum(pHpcMUA)  * w;
                rec_sum_v1   = rec_sum_v1   + sum(pV1MUA)   * w;
                rec_mean_hpc = rec_mean_hpc + mean(pHpcMUA) * w;
                rec_mean_v1  = rec_mean_v1  + mean(pV1MUA)  * w;
            end
            recency_weighted_past_HPC_MUA_sum_col(iUP)  = rec_sum_hpc;
            recency_weighted_past_V1_MUA_sum_col(iUP)   = rec_sum_v1;
            recency_weighted_past_HPC_MUA_mean_col(iUP) = rec_mean_hpc / n_past;
            recency_weighted_past_V1_MUA_mean_col(iUP)  = rec_mean_v1 / n_past;
            
            % Past ripples MUA in distance windows relative to last ripple peak (0-100ms, 100-200ms, 200+ms)
            past_mask_0_100   = false(1, nBins);
            past_mask_100_200 = false(1, nBins);
            past_mask_200plus = false(1, nBins);
            
            for pIdx = 1:n_past
                dist_p = ripAbsPk(end) - pastAbsPk(pIdx);
                pOn_rel  = pastAbsOn(pIdx) - upOnset;
                pOff_rel = pastAbsOff(pIdx) - upOnset;
                
                if dist_p >= 0 && dist_p < 0.100
                    past_mask_0_100 = past_mask_0_100 | (binCenters >= pOn_rel & binCenters < pOff_rel);
                elseif dist_p >= 0.100 && dist_p < 0.200
                    past_mask_100_200 = past_mask_100_200 | (binCenters >= pOn_rel & binCenters < pOff_rel);
                elseif dist_p >= 0.200
                    past_mask_200plus = past_mask_200plus | (binCenters >= pOn_rel & binCenters < pOff_rel);
                end
            end
            
            if any(past_mask_0_100)
                past_ripples_0_100ms_HPC_MUA_sum_col(iUP)  = sum(hpcMUA(past_mask_0_100));
                past_ripples_0_100ms_HPC_MUA_mean_col(iUP) = mean(hpcMUA(past_mask_0_100));
                past_ripples_0_100ms_V1_MUA_sum_col(iUP)   = sum(v1MUA(past_mask_0_100));
                past_ripples_0_100ms_V1_MUA_mean_col(iUP)  = mean(v1MUA(past_mask_0_100));
            end
            
            if any(past_mask_100_200)
                past_ripples_100_200ms_HPC_MUA_sum_col(iUP)  = sum(hpcMUA(past_mask_100_200));
                past_ripples_100_200ms_HPC_MUA_mean_col(iUP) = mean(hpcMUA(past_mask_100_200));
                past_ripples_100_200ms_V1_MUA_sum_col(iUP)   = sum(v1MUA(past_mask_100_200));
                past_ripples_100_200ms_V1_MUA_mean_col(iUP)  = mean(v1MUA(past_mask_100_200));
            end
            
            if any(past_mask_200plus)
                past_ripples_200plus_ms_HPC_MUA_sum_col(iUP)  = sum(hpcMUA(past_mask_200plus));
                past_ripples_200plus_ms_HPC_MUA_mean_col(iUP) = mean(hpcMUA(past_mask_200plus));
                past_ripples_200plus_ms_V1_MUA_sum_col(iUP)   = sum(v1MUA(past_mask_200plus));
                past_ripples_200plus_ms_V1_MUA_mean_col(iUP)  = mean(v1MUA(past_mask_200plus));
            end
        end
        
        % Factor 3: Non-Ripple Activity
        if strcmp(nr_scope, 'prior_to_last_ripple')
            eval_window_max = last_abs_on - upOnset;
            allRipOn_rel  = ripAbsOn(1:end-1) - upOnset;
            allRipOff_rel = ripAbsOff(1:end-1) - upOnset;
        else % 'entire_UP'
            eval_window_max = upDur;
            allRipOn_rel  = ripAbsOn - upOnset;
            allRipOff_rel = ripAbsOff - upOnset;
        end
        
        nonRipBinMask = (binCenters <= eval_window_max);
        for rIdx = 1:length(allRipOn_rel)
            nonRipBinMask = nonRipBinMask & ~(binCenters >= allRipOn_rel(rIdx) & binCenters < allRipOff_rel(rIdx));
        end
        
        totalRipDurInWindow = 0;
        for rIdx = 1:length(allRipOn_rel)
            totalRipDurInWindow = totalRipDurInWindow + max(0, min(eval_window_max, allRipOff_rel(rIdx)) - max(0, allRipOn_rel(rIdx)));
        end
        non_ripple_duration_col(iUP) = max(0, eval_window_max - totalRipDurInWindow);
        
        if any(nonRipBinMask)
            non_ripple_HPC_MUA_sum_col(iUP)  = sum(hpcMUA(nonRipBinMask));
            non_ripple_HPC_MUA_mean_col(iUP) = mean(hpcMUA(nonRipBinMask));
            non_ripple_V1_MUA_sum_col(iUP)   = sum(v1MUA(nonRipBinMask));
            non_ripple_V1_MUA_mean_col(iUP)  = mean(v1MUA(nonRipBinMask));
        else
            non_ripple_HPC_MUA_sum_col(iUP)  = 0;
            non_ripple_HPC_MUA_mean_col(iUP) = 0;
            non_ripple_V1_MUA_sum_col(iUP)   = 0;
            non_ripple_V1_MUA_mean_col(iUP)  = 0;
        end
        
        % Non-ripple MUA prior to last ripple peak in distance windows (0-100ms, 100-200ms, 200+ms)
        last_pk_rel = last_abs_pk - upOnset;
        nonRipPriorMask = (binCenters < last_pk_rel) & ~isRipBin;
        dist_bin = last_pk_rel - binCenters;
        
        nonrip_mask_0_100   = nonRipPriorMask & (dist_bin >= 0 & dist_bin < 0.100);
        nonrip_mask_100_200 = nonRipPriorMask & (dist_bin >= 0.100 & dist_bin < 0.200);
        nonrip_mask_200plus = nonRipPriorMask & (dist_bin >= 0.200);
        
        if any(nonrip_mask_0_100)
            non_ripple_0_100ms_HPC_MUA_sum_col(iUP)  = sum(hpcMUA(nonrip_mask_0_100));
            non_ripple_0_100ms_HPC_MUA_mean_col(iUP) = mean(hpcMUA(nonrip_mask_0_100));
            non_ripple_0_100ms_V1_MUA_sum_col(iUP)   = sum(v1MUA(nonrip_mask_0_100));
            non_ripple_0_100ms_V1_MUA_mean_col(iUP)  = mean(v1MUA(nonrip_mask_0_100));
        end
        
        if any(nonrip_mask_100_200)
            non_ripple_100_200ms_HPC_MUA_sum_col(iUP)  = sum(hpcMUA(nonrip_mask_100_200));
            non_ripple_100_200ms_HPC_MUA_mean_col(iUP) = mean(hpcMUA(nonrip_mask_100_200));
            non_ripple_100_200ms_V1_MUA_sum_col(iUP)   = sum(v1MUA(nonrip_mask_100_200));
            non_ripple_100_200ms_V1_MUA_mean_col(iUP)  = mean(v1MUA(nonrip_mask_100_200));
        end
        
        if any(nonrip_mask_200plus)
            non_ripple_200plus_ms_HPC_MUA_sum_col(iUP)  = sum(hpcMUA(nonrip_mask_200plus));
            non_ripple_200plus_ms_HPC_MUA_mean_col(iUP) = mean(hpcMUA(nonrip_mask_200plus));
            non_ripple_200plus_ms_V1_MUA_sum_col(iUP)   = sum(v1MUA(nonrip_mask_200plus));
            non_ripple_200plus_ms_V1_MUA_mean_col(iUP)  = mean(v1MUA(nonrip_mask_200plus));
        end
    end
end

%% 4. Assemble output table
T = table(...
    upID_col, ...
    session_id_col, ...
    subject_id_col, ...
    hemisphere_id_col, ...
    up_duration_col, ...
    ripple_count_col, ...
    last_ripple_to_UP_term_col, ...
    last_ripple_duration_col, ...
    last_ripple_power_col, ...
    last_ripple_interval_col, ...
    last_ripple_HPC_MUA_sum_col, ...
    last_ripple_HPC_MUA_mean_col, ...
    last_ripple_V1_MUA_sum_col, ...
    last_ripple_V1_MUA_mean_col, ...
    last_ripple_HPC_MUA_peak_col, ...
    last_ripple_V1_MUA_peak_col, ...
    last_ripple_first_half_HPC_MUA_sum_col, ...
    last_ripple_first_half_HPC_MUA_mean_col, ...
    last_ripple_second_half_HPC_MUA_sum_col, ...
    last_ripple_second_half_HPC_MUA_mean_col, ...
    last_ripple_first_half_V1_MUA_sum_col, ...
    last_ripple_first_half_V1_MUA_mean_col, ...
    last_ripple_second_half_V1_MUA_sum_col, ...
    last_ripple_second_half_V1_MUA_mean_col, ...
    past_ripples_count_col, ...
    past_ripples_duration_col, ...
    past_ripples_HPC_MUA_sum_col, ...
    past_ripples_HPC_MUA_mean_col, ...
    past_ripples_V1_MUA_sum_col, ...
    past_ripples_V1_MUA_mean_col, ...
    first_half_past_ripples_HPC_MUA_sum_col, ...
    first_half_past_ripples_HPC_MUA_mean_col, ...
    second_half_past_ripples_HPC_MUA_sum_col, ...
    second_half_past_ripples_HPC_MUA_mean_col, ...
    first_half_past_ripples_V1_MUA_sum_col, ...
    first_half_past_ripples_V1_MUA_mean_col, ...
    second_half_past_ripples_V1_MUA_sum_col, ...
    second_half_past_ripples_V1_MUA_mean_col, ...
    first_half_past_ripples_duration_col, ...
    second_half_past_ripples_duration_col, ...
    recency_weighted_past_HPC_MUA_sum_col, ...
    recency_weighted_past_HPC_MUA_mean_col, ...
    recency_weighted_past_V1_MUA_sum_col, ...
    recency_weighted_past_V1_MUA_mean_col, ...
    past_ripples_0_100ms_HPC_MUA_sum_col, ...
    past_ripples_0_100ms_HPC_MUA_mean_col, ...
    past_ripples_0_100ms_V1_MUA_sum_col, ...
    past_ripples_0_100ms_V1_MUA_mean_col, ...
    past_ripples_100_200ms_HPC_MUA_sum_col, ...
    past_ripples_100_200ms_HPC_MUA_mean_col, ...
    past_ripples_100_200ms_V1_MUA_sum_col, ...
    past_ripples_100_200ms_V1_MUA_mean_col, ...
    past_ripples_200plus_ms_HPC_MUA_sum_col, ...
    past_ripples_200plus_ms_HPC_MUA_mean_col, ...
    past_ripples_200plus_ms_V1_MUA_sum_col, ...
    past_ripples_200plus_ms_V1_MUA_mean_col, ...
    non_ripple_0_100ms_HPC_MUA_sum_col, ...
    non_ripple_0_100ms_HPC_MUA_mean_col, ...
    non_ripple_0_100ms_V1_MUA_sum_col, ...
    non_ripple_0_100ms_V1_MUA_mean_col, ...
    non_ripple_100_200ms_HPC_MUA_sum_col, ...
    non_ripple_100_200ms_HPC_MUA_mean_col, ...
    non_ripple_100_200ms_V1_MUA_sum_col, ...
    non_ripple_100_200ms_V1_MUA_mean_col, ...
    non_ripple_200plus_ms_HPC_MUA_sum_col, ...
    non_ripple_200plus_ms_HPC_MUA_mean_col, ...
    non_ripple_200plus_ms_V1_MUA_sum_col, ...
    non_ripple_200plus_ms_V1_MUA_mean_col, ...
    non_ripple_duration_col, ...
    non_ripple_HPC_MUA_sum_col, ...
    non_ripple_HPC_MUA_mean_col, ...
    non_ripple_V1_MUA_sum_col, ...
    non_ripple_V1_MUA_mean_col, ...
    first_half_HPC_MUA_sum_col, ...
    first_half_HPC_MUA_mean_col, ...
    second_half_HPC_MUA_sum_col, ...
    second_half_HPC_MUA_mean_col, ...
    first_half_V1_MUA_sum_col, ...
    first_half_V1_MUA_mean_col, ...
    second_half_V1_MUA_sum_col, ...
    second_half_V1_MUA_mean_col, ...
    first_half_ripple_HPC_MUA_sum_col, ...
    first_half_ripple_HPC_MUA_mean_col, ...
    second_half_ripple_HPC_MUA_sum_col, ...
    second_half_ripple_HPC_MUA_mean_col, ...
    first_half_ripple_V1_MUA_sum_col, ...
    first_half_ripple_V1_MUA_mean_col, ...
    second_half_ripple_V1_MUA_sum_col, ...
    second_half_ripple_V1_MUA_mean_col, ...
    first_half_non_ripple_HPC_MUA_sum_col, ...
    first_half_non_ripple_HPC_MUA_mean_col, ...
    second_half_non_ripple_HPC_MUA_sum_col, ...
    second_half_non_ripple_HPC_MUA_mean_col, ...
    first_half_non_ripple_V1_MUA_sum_col, ...
    first_half_non_ripple_V1_MUA_mean_col, ...
    second_half_non_ripple_V1_MUA_sum_col, ...
    second_half_non_ripple_V1_MUA_mean_col, ...
    first_half_ripple_duration_col, ...
    second_half_ripple_duration_col, ...
    first_half_non_ripple_duration_col, ...
    second_half_non_ripple_duration_col, ...
    'VariableNames', { ...
    'upID', 'session_id', 'subject_id', 'hemisphere_id', ...
    'up_duration', 'ripple_count', 'last_ripple_to_UP_term', ...
    'last_ripple_duration', 'last_ripple_power', 'last_ripple_interval', ...
    'last_ripple_HPC_MUA_sum', 'last_ripple_HPC_MUA_mean', ...
    'last_ripple_V1_MUA_sum', 'last_ripple_V1_MUA_mean', ...
    'last_ripple_HPC_MUA_peak', 'last_ripple_V1_MUA_peak', ...
    'last_ripple_first_half_HPC_MUA_sum', 'last_ripple_first_half_HPC_MUA_mean', ...
    'last_ripple_second_half_HPC_MUA_sum', 'last_ripple_second_half_HPC_MUA_mean', ...
    'last_ripple_first_half_V1_MUA_sum', 'last_ripple_first_half_V1_MUA_mean', ...
    'last_ripple_second_half_V1_MUA_sum', 'last_ripple_second_half_V1_MUA_mean', ...
    'past_ripples_count', 'past_ripples_duration', ...
    'past_ripples_HPC_MUA_sum', 'past_ripples_HPC_MUA_mean', ...
    'past_ripples_V1_MUA_sum', 'past_ripples_V1_MUA_mean', ...
    'first_half_past_ripples_HPC_MUA_sum', 'first_half_past_ripples_HPC_MUA_mean', ...
    'second_half_past_ripples_HPC_MUA_sum', 'second_half_past_ripples_HPC_MUA_mean', ...
    'first_half_past_ripples_V1_MUA_sum', 'first_half_past_ripples_V1_MUA_mean', ...
    'second_half_past_ripples_V1_MUA_sum', 'second_half_past_ripples_V1_MUA_mean', ...
    'first_half_past_ripples_duration', 'second_half_past_ripples_duration', ...
    'recency_weighted_past_HPC_MUA_sum', 'recency_weighted_past_HPC_MUA_mean', ...
    'recency_weighted_past_V1_MUA_sum', 'recency_weighted_past_V1_MUA_mean', ...
    'past_ripples_0_100ms_HPC_MUA_sum', 'past_ripples_0_100ms_HPC_MUA_mean', ...
    'past_ripples_0_100ms_V1_MUA_sum', 'past_ripples_0_100ms_V1_MUA_mean', ...
    'past_ripples_100_200ms_HPC_MUA_sum', 'past_ripples_100_200ms_HPC_MUA_mean', ...
    'past_ripples_100_200ms_V1_MUA_sum', 'past_ripples_100_200ms_V1_MUA_mean', ...
    'past_ripples_200plus_ms_HPC_MUA_sum', 'past_ripples_200plus_ms_HPC_MUA_mean', ...
    'past_ripples_200plus_ms_V1_MUA_sum', 'past_ripples_200plus_ms_V1_MUA_mean', ...
    'non_ripple_0_100ms_HPC_MUA_sum', 'non_ripple_0_100ms_HPC_MUA_mean', ...
    'non_ripple_0_100ms_V1_MUA_sum', 'non_ripple_0_100ms_V1_MUA_mean', ...
    'non_ripple_100_200ms_HPC_MUA_sum', 'non_ripple_100_200ms_HPC_MUA_mean', ...
    'non_ripple_100_200ms_V1_MUA_sum', 'non_ripple_100_200ms_V1_MUA_mean', ...
    'non_ripple_200plus_ms_HPC_MUA_sum', 'non_ripple_200plus_ms_HPC_MUA_mean', ...
    'non_ripple_200plus_ms_V1_MUA_sum', 'non_ripple_200plus_ms_V1_MUA_mean', ...
    'non_ripple_duration', 'non_ripple_HPC_MUA_sum', ...
    'non_ripple_HPC_MUA_mean', 'non_ripple_V1_MUA_sum', ...
    'non_ripple_V1_MUA_mean', ...
    'first_half_HPC_MUA_sum', 'first_half_HPC_MUA_mean', ...
    'second_half_HPC_MUA_sum', 'second_half_HPC_MUA_mean', ...
    'first_half_V1_MUA_sum', 'first_half_V1_MUA_mean', ...
    'second_half_V1_MUA_sum', 'second_half_V1_MUA_mean', ...
    'first_half_ripple_HPC_MUA_sum', 'first_half_ripple_HPC_MUA_mean', ...
    'second_half_ripple_HPC_MUA_sum', 'second_half_ripple_HPC_MUA_mean', ...
    'first_half_ripple_V1_MUA_sum', 'first_half_ripple_V1_MUA_mean', ...
    'second_half_ripple_V1_MUA_sum', 'second_half_ripple_V1_MUA_mean', ...
    'first_half_non_ripple_HPC_MUA_sum', 'first_half_non_ripple_HPC_MUA_mean', ...
    'second_half_non_ripple_HPC_MUA_sum', 'second_half_non_ripple_HPC_MUA_mean', ...
    'first_half_non_ripple_V1_MUA_sum', 'first_half_non_ripple_V1_MUA_mean', ...
    'second_half_non_ripple_V1_MUA_sum', 'second_half_non_ripple_V1_MUA_mean', ...
    'first_half_ripple_duration', 'second_half_ripple_duration', ...
    'first_half_non_ripple_duration', 'second_half_non_ripple_duration'} ...
);

end
