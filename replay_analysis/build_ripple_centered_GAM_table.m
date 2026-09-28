function ripple_tbl = build_ripple_centered_GAM_table(merged_event_info, V1_MUA_spiketimes, HC_MUA_spiketimes, ...
    ripple_HC_logodds, ripple_HC_logodds_PRE, ripple_V1_logodds, ripple_V1_logodds_PRE, ...
    early_UP_HPC_log_odds, early_UP_V1_log_odds, ...
    late_UP_HPC_log_odds, late_UP_V1_log_odds, ...
    next_UP_HPC_log_odds, next_UP_V1_log_odds, ...
    next_DOWN_duration, previous_DOWN_duration, ...
    next_DOWN_SO_power, previous_DOWN_SO_power, ...
    next_DOWN_lag, UP_lag, varargin)
% BUILD_RIPPLE_CENTERED_GAM_TABLE  One row per ripple event (not per UP state),
% for GAMMs that use ALL ripples within a UP state (evidence-accumulation /
% within-UP ripple-order analyses), rather than only the first/last ripple as
% in the UP-level UP_DOWN_info_GAM table.
%
% Each ripple row carries: which UP it belongs to, its serial position within
% that UP (1st, 2nd, 3rd, ... and from the end), timing relative to UP onset/
% offset/previous-ripple/next-ripple, its own HPC/V1 MUA and reactivation
% log-odds bias (post- and pre-ripple windows) and derived V1-HPC coherence,
% plus the UP-level outcome features (early/late/nextUP log-odds, next/
% previous DOWN duration, SO power and bilateral lag) broadcast from its
% parent UP so within-UP ripple order/content can be related to those
% UP-level outcomes.
%
% Inputs:
%   merged_event_info - struct with UP_ints, UP_hemisphere_id, ripples_ints,
%                        ripples_peaktimes, ripples_power, session_id, subject_id
%                        (deduplicated/merged across hemispheres, as built in
%                        main_build_UP_DOWN_info_GAM_table.m)
%   V1_MUA_spiketimes, HC_MUA_spiketimes - Cell {1}=Left, {2}=Right raw spike times
%   ripple_HC_logodds, ripple_HC_logodds_PRE - Mx1, HPC reactivation log-odds bias
%                                              per ripple (post-/pre-ripple windows)
%   ripple_V1_logodds, ripple_V1_logodds_PRE - Mx1, V1 reactivation log-odds bias
%                                              per ripple (post-/pre-ripple windows)
%   early_UP_HPC_log_odds, early_UP_V1_log_odds   - Nx1 per-UP early-UP (0-100ms
%                                                    after UP onset) log-odds
%   late_UP_HPC_log_odds, late_UP_V1_log_odds     - Nx1 per-UP late-UP (100ms
%                                                    before UP termination) log-odds
%   next_UP_HPC_log_odds, next_UP_V1_log_odds     - Nx1 per-UP early-window log-odds
%                                                    of the FOLLOWING UP state
%   next_DOWN_duration, previous_DOWN_duration    - Nx1 per-UP next/previous DOWN duration
%   next_DOWN_SO_power, previous_DOWN_SO_power    - Nx1 per-UP SO peak magnitude at the
%                                                    UP->DOWN / DOWN->UP transition
%   next_DOWN_lag, UP_lag                         - Nx1 per-UP bilateral (ipsi/contra)
%                                                    lag at the UP->DOWN / DOWN->UP transition
%
% Options:
%   'time_bin'         - 0.01 (10ms bin size, default)
%   'min_interval_dur' - 1e-5 minimum ripple duration threshold
%
% Output:
%   ripple_tbl - MATLAB table, one row per ripple event that occurs inside a UP state.

p = inputParser;
addParameter(p, 'time_bin', 0.01, @isnumeric);
addParameter(p, 'min_interval_dur', 1e-5, @isnumeric);
parse(p, varargin{:});

time_bin = p.Results.time_bin;
min_dur  = p.Results.min_interval_dur;

nUP = size(merged_event_info.UP_ints, 1);
ripIntsAbs  = merged_event_info.ripples_ints;
ripPeakAbs  = merged_event_info.ripples_peaktimes;
ripPowerAbs = merged_event_info.ripples_power;

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
% (identical convention to build_last_ripple_to_DOWN_summary_table.m /
% build_UP_counting_process_intervals.m, so ripple MUA rates here are on the
% same scale as the UP-level and interval-level tables)
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

%% 2. Iterate over each UP event and emit one row per ripple inside it
fprintf('Extracting per-ripple features for the ripple-centered GAM table...\n');

upID_list               = cell(nUP, 1);
sessionID_list          = cell(nUP, 1);
animalID_list           = cell(nUP, 1);
hemisphereID_list       = cell(nUP, 1);
UPDuration_list         = cell(nUP, 1);

rippleIndexInUP_list     = cell(nUP, 1);
rippleIndexFromEnd_list  = cell(nUP, 1);
numRipplesInUP_list      = cell(nUP, 1);
isFirstRipple_list       = cell(nUP, 1);
isLastRipple_list        = cell(nUP, 1);

rippleDuration_list         = cell(nUP, 1);
ripplePower_list            = cell(nUP, 1);
timeFromUPonset_list        = cell(nUP, 1);
timeToUPend_list             = cell(nUP, 1);
rippleNormalisedUP_list      = cell(nUP, 1);
timeFromPreviousRipple_list  = cell(nUP, 1);
timeToNextRipple_list        = cell(nUP, 1);

rippleHC_logodds_list      = cell(nUP, 1);
rippleV1_logodds_list      = cell(nUP, 1);
rippleHC_logoddsPRE_list   = cell(nUP, 1);
rippleV1_logoddsPRE_list   = cell(nUP, 1);
rippleCoherence_list       = cell(nUP, 1);
rippleCoherencePRE_list    = cell(nUP, 1);

rippleHPC_MUA_sum_list   = cell(nUP, 1);
rippleHPC_MUA_mean_list  = cell(nUP, 1);
rippleV1_MUA_sum_list    = cell(nUP, 1);
rippleV1_MUA_mean_list   = cell(nUP, 1);

earlyUPHPC_list  = cell(nUP, 1);
earlyUPV1_list   = cell(nUP, 1);
lateUPHPC_list   = cell(nUP, 1);
lateUPV1_list    = cell(nUP, 1);
nextUPHPC_list   = cell(nUP, 1);
nextUPV1_list    = cell(nUP, 1);

nextDOWNDuration_list     = cell(nUP, 1);
previousDOWNDuration_list = cell(nUP, 1);
nextDOWNSOPower_list      = cell(nUP, 1);
previousDOWNSOPower_list  = cell(nUP, 1);
nextDOWNlag_list          = cell(nUP, 1);
UPlag_list                = cell(nUP, 1);

for iUP = 1:nUP
    upOnset  = merged_event_info.UP_ints(iUP, 1);
    upOffset = merged_event_info.UP_ints(iUP, 2);
    upDur    = upOffset - upOnset;
    hemiID   = merged_event_info.UP_hemisphere_id(iUP);
    sessID   = session_id(iUP);
    subjID   = subject_id(iUP);

    overlapIdx = find(ripPeakAbs(:,1) > upOnset & ripPeakAbs(:,1) < upOffset);
    if isempty(overlapIdx)
        continue
    end

    [~, sIdx] = sort(ripIntsAbs(overlapIdx, 1));
    ripples_index = overlapIdx(sIdx);

    % Drop ripples with (near-)zero duration, same threshold as other builders
    keep = (ripIntsAbs(ripples_index,2) - ripIntsAbs(ripples_index,1)) > min_dur;
    ripples_index = ripples_index(keep);
    if isempty(ripples_index)
        continue
    end

    k = length(ripples_index);

    s_offset = sessID * 1000000;

    % Own HPC/V1 MUA for each ripple in this UP (relative to session offset)
    hpcMUAmean = nan(k,1); hpcMUAsum = nan(k,1);
    v1MUAmean  = nan(k,1); v1MUAsum  = nan(k,1);
    for r = 1:k
        rOn_abs  = ripIntsAbs(ripples_index(r), 1);
        rOff_abs = ripIntsAbs(ripples_index(r), 2);

        rel_on  = rOn_abs  - s_offset;
        rel_off = rOff_abs - s_offset;

        bStart = max(1, floor(rel_on / time_bin) + 1);
        bEnd   = min(length(norm_v1_L{sessID}), ceil(rel_off / time_bin));
        if bStart > bEnd, bEnd = bStart; end

        rHpcMUA = mean([norm_hc_L{sessID}(bStart:bEnd); norm_hc_R{sessID}(bStart:bEnd)], 1, 'omitnan');
        if hemiID == 1
            rV1MUA = norm_v1_L{sessID}(bStart:bEnd);
        else
            rV1MUA = norm_v1_R{sessID}(bStart:bEnd);
        end

        hpcMUAsum(r)  = sum(rHpcMUA);
        hpcMUAmean(r) = mean(rHpcMUA);
        v1MUAsum(r)   = sum(rV1MUA);
        v1MUAmean(r)  = mean(rV1MUA);
    end

    % Timing features relative to UP onset/offset and neighboring ripples
    ripOn_rel  = ripIntsAbs(ripples_index,1) - upOnset;
    ripOff_rel = ripIntsAbs(ripples_index,2) - upOnset;
    ripDur     = ripIntsAbs(ripples_index,2) - ripIntsAbs(ripples_index,1);
    ripPk_abs  = ripPeakAbs(ripples_index);

    timeFromUPonset = ripOn_rel;
    timeToUPend     = upDur - ripOff_rel;

    rippleNormalisedUP = ripOn_rel ./ upDur;
    rippleNormalisedUP(rippleNormalisedUP < 0) = 0;
    rippleNormalisedUP(rippleNormalisedUP > 1) = 1;

    timeFromPreviousRipple = nan(k,1);
    timeToNextRipple       = nan(k,1);
    if k > 1
        timeFromPreviousRipple(2:end) = diff(ripPk_abs);
        timeToNextRipple(1:end-1)     = diff(ripPk_abs);
    end

    rippleIndexInUP    = (1:k)';
    rippleIndexFromEnd = (k:-1:1)';

    % Per-ripple reactivation log-odds bias and derived V1-HPC coherence
    hc     = ripple_HC_logodds(ripples_index);
    v1     = ripple_V1_logodds(ripples_index);
    hcPRE  = ripple_HC_logodds_PRE(ripples_index);
    v1PRE  = ripple_V1_logodds_PRE(ripples_index);

    coherence    = sign(v1 .* hc) .* sqrt(abs(v1 .* hc));
    coherencePRE = sign(v1PRE .* hcPRE) .* sqrt(abs(v1PRE .* hcPRE));

    % Store this UP event's ripple rows
    upID_list{iUP}         = repmat(iUP, k, 1);
    sessionID_list{iUP}    = repmat(sessID, k, 1);
    animalID_list{iUP}     = repmat(subjID, k, 1);
    hemisphereID_list{iUP} = repmat(hemiID, k, 1);
    UPDuration_list{iUP}   = repmat(upDur, k, 1);

    rippleIndexInUP_list{iUP}    = rippleIndexInUP;
    rippleIndexFromEnd_list{iUP} = rippleIndexFromEnd;
    numRipplesInUP_list{iUP}     = repmat(k, k, 1);
    isFirstRipple_list{iUP}      = double(rippleIndexInUP == 1);
    isLastRipple_list{iUP}       = double(rippleIndexFromEnd == 1);

    rippleDuration_list{iUP}        = ripDur;
    ripplePower_list{iUP}           = ripPowerAbs(ripples_index);
    timeFromUPonset_list{iUP}       = timeFromUPonset;
    timeToUPend_list{iUP}           = timeToUPend;
    rippleNormalisedUP_list{iUP}    = rippleNormalisedUP;
    timeFromPreviousRipple_list{iUP} = timeFromPreviousRipple;
    timeToNextRipple_list{iUP}       = timeToNextRipple;

    rippleHC_logodds_list{iUP}    = hc;
    rippleV1_logodds_list{iUP}    = v1;
    rippleHC_logoddsPRE_list{iUP} = hcPRE;
    rippleV1_logoddsPRE_list{iUP} = v1PRE;
    rippleCoherence_list{iUP}     = coherence;
    rippleCoherencePRE_list{iUP}  = coherencePRE;

    rippleHPC_MUA_sum_list{iUP}  = hpcMUAsum;
    rippleHPC_MUA_mean_list{iUP} = hpcMUAmean;
    rippleV1_MUA_sum_list{iUP}   = v1MUAsum;
    rippleV1_MUA_mean_list{iUP}  = v1MUAmean;

    earlyUPHPC_list{iUP} = repmat(early_UP_HPC_log_odds(iUP), k, 1);
    earlyUPV1_list{iUP}  = repmat(early_UP_V1_log_odds(iUP),  k, 1);
    lateUPHPC_list{iUP}  = repmat(late_UP_HPC_log_odds(iUP),  k, 1);
    lateUPV1_list{iUP}   = repmat(late_UP_V1_log_odds(iUP),   k, 1);
    nextUPHPC_list{iUP}  = repmat(next_UP_HPC_log_odds(iUP),  k, 1);
    nextUPV1_list{iUP}   = repmat(next_UP_V1_log_odds(iUP),   k, 1);

    nextDOWNDuration_list{iUP}     = repmat(next_DOWN_duration(iUP),     k, 1);
    previousDOWNDuration_list{iUP} = repmat(previous_DOWN_duration(iUP), k, 1);
    nextDOWNSOPower_list{iUP}      = repmat(next_DOWN_SO_power(iUP),     k, 1);
    previousDOWNSOPower_list{iUP}  = repmat(previous_DOWN_SO_power(iUP), k, 1);
    nextDOWNlag_list{iUP}          = repmat(next_DOWN_lag(iUP), k, 1);
    UPlag_list{iUP}                = repmat(UP_lag(iUP),        k, 1);
end

ripple_tbl = table(...
    vertcat(upID_list{:}), ...
    vertcat(sessionID_list{:}), ...
    vertcat(animalID_list{:}), ...
    vertcat(hemisphereID_list{:}), ...
    vertcat(UPDuration_list{:}), ...
    vertcat(rippleIndexInUP_list{:}), ...
    vertcat(rippleIndexFromEnd_list{:}), ...
    vertcat(numRipplesInUP_list{:}), ...
    vertcat(isFirstRipple_list{:}), ...
    vertcat(isLastRipple_list{:}), ...
    vertcat(rippleDuration_list{:}), ...
    vertcat(ripplePower_list{:}), ...
    vertcat(timeFromUPonset_list{:}), ...
    vertcat(timeToUPend_list{:}), ...
    vertcat(rippleNormalisedUP_list{:}), ...
    vertcat(timeFromPreviousRipple_list{:}), ...
    vertcat(timeToNextRipple_list{:}), ...
    vertcat(rippleHC_logodds_list{:}), ...
    vertcat(rippleV1_logodds_list{:}), ...
    vertcat(rippleHC_logoddsPRE_list{:}), ...
    vertcat(rippleV1_logoddsPRE_list{:}), ...
    vertcat(rippleCoherence_list{:}), ...
    vertcat(rippleCoherencePRE_list{:}), ...
    vertcat(rippleHPC_MUA_sum_list{:}), ...
    vertcat(rippleHPC_MUA_mean_list{:}), ...
    vertcat(rippleV1_MUA_sum_list{:}), ...
    vertcat(rippleV1_MUA_mean_list{:}), ...
    vertcat(earlyUPHPC_list{:}), ...
    vertcat(earlyUPV1_list{:}), ...
    vertcat(lateUPHPC_list{:}), ...
    vertcat(lateUPV1_list{:}), ...
    vertcat(nextUPHPC_list{:}), ...
    vertcat(nextUPV1_list{:}), ...
    vertcat(nextDOWNDuration_list{:}), ...
    vertcat(previousDOWNDuration_list{:}), ...
    vertcat(nextDOWNSOPower_list{:}), ...
    vertcat(previousDOWNSOPower_list{:}), ...
    vertcat(nextDOWNlag_list{:}), ...
    vertcat(UPlag_list{:}), ...
    'VariableNames', { ...
    'upID', 'SessionID', 'AnimalID', 'hemisphere_id', 'UPDuration', ...
    'rippleIndexInUP', 'rippleIndexFromEnd', 'numRipplesInUP', 'isFirstRipple', 'isLastRipple', ...
    'rippleDuration', 'ripplePower', ...
    'timeFromUPonset', 'timeToUPend', 'rippleNormalisedUP', ...
    'timeFromPreviousRipple', 'timeToNextRipple', ...
    'rippleHC_logodds', 'rippleV1_logodds', 'rippleHC_logoddsPRE', 'rippleV1_logoddsPRE', ...
    'rippleCoherence', 'rippleCoherencePRE', ...
    'rippleHPC_MUA_sum', 'rippleHPC_MUA_mean', 'rippleV1_MUA_sum', 'rippleV1_MUA_mean', ...
    'earlyUPHPC', 'earlyUPV1', 'lateUPHPC', 'lateUPV1', 'nextUPHPC', 'nextUPV1', ...
    'nextDOWNDuration', 'previousDOWNDuration', 'nextDOWNSOPower', 'previousDOWNSOPower', ...
    'nextDOWNlag', 'UPlag'} ...
    );

end
