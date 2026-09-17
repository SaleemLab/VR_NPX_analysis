function ripple_info_table = build_ripple_log_odds_spindle_table(merged_event_info, ...
    mean_z_bias, mean_z_bias_PRE, mean_z_bias_V1, mean_z_bias_V1_PRE, ...
    spindle_amplitude_temporal, spindle_tvec, UP_next_onset, varargin)
% BUILD_RIPPLE_LOG_ODDS_SPINDLE_TABLE Per-UP-event ripple identity, timing,
% reactivation log-odds bias, and periripple spindle amplitude features.
%
% This complements BUILD_LAST_RIPPLE_TO_DOWN_SUMMARY_TABLE (which covers
% MUA-derived features) with the non-MUA ripple features needed for the
% UP_DOWN_info / GAM table: which ripple was first/last, its timing
% relative to the UP state, reactivation log-odds bias around it, and
% spindle amplitude around the last ripple.
%
% Inputs:
%   merged_event_info - struct with UP_ints, ripples_ints, ripples_peaktimes,
%                        ripples_power (deduplicated/merged across hemispheres)
%   mean_z_bias, mean_z_bias_PRE       - 1xM HPC reactivation log-odds per ripple (post/pre window)
%   mean_z_bias_V1, mean_z_bias_V1_PRE - 1xM V1 reactivation log-odds per ripple (post/pre window)
%   spindle_amplitude_temporal - TxMxR periripple spindle amplitude (time x ripple x region)
%   spindle_tvec               - Tx1 time vector for spindle_amplitude_temporal (relative to ripple onset)
%   UP_next_onset              - Nx1 absolute onset time of the next UP state on the same
%                                 probe/hemisphere (NaN if this is the last UP in its session)
%
% Options:
%   'min_interval_dur'        - 1e-5 minimum ripple duration to count (default)
%   'spindle_pre_ripple_time' - -0.1 s, time point sampled for pre-ripple spindle power (default)

p = inputParser;
addParameter(p, 'min_interval_dur', 0.03, @isnumeric);
addParameter(p, 'spindle_pre_ripple_time', -0.1, @isnumeric);
addParameter(p, 'time_reference', 'boundary', @(x) ismember(x, {'boundary', 'peak'}));

parse(p, varargin{:});
time_ref    = p.Results.time_reference;
min_dur = p.Results.min_interval_dur;

nUP         = size(merged_event_info.UP_ints, 1);
ripIntsAbs  = merged_event_info.ripples_ints;
ripPeakAbs  = merged_event_info.ripples_peaktimes;
ripPowerAbs = merged_event_info.ripples_power;
nRegions    = size(spindle_amplitude_temporal, 3);

[~, spindle_pre_bin] = min(abs(spindle_tvec - p.Results.spindle_pre_ripple_time));

first_ripples_index     = nan(nUP, 1);
last_ripples_index      = nan(nUP, 1);
time_to_first_ripples   = nan(nUP, 1);
time_from_first_ripples = nan(nUP, 1);
time_to_last_ripples    = nan(nUP, 1);
time_from_last_ripples  = nan(nUP, 1);
last_ripples_duration =nan(nUP, 1);
first_ripples_duration =nan(nUP, 1);
first_ripples_power     = nan(nUP, 1);
last_ripples_power      = nan(nUP, 1);

first_ripples_log_odds        = nan(nUP, 1);
last_ripples_log_odds         = nan(nUP, 1);
first_ripples_PRE_log_odds    = nan(nUP, 1);
last_ripples_PRE_log_odds     = nan(nUP, 1);
first_ripples_V1_log_odds     = nan(nUP, 1);
last_ripples_V1_log_odds      = nan(nUP, 1);
first_ripples_PRE_V1_log_odds = nan(nUP, 1);
last_ripples_PRE_V1_log_odds  = nan(nUP, 1);

second_last_ripples_log_odds        = nan(nUP, 1);
second_last_ripples_PRE_V1_log_odds = nan(nUP, 1);
second_last_ripples_V1_log_odds     = nan(nUP, 1);
previous_ripples_log_odds           = nan(nUP, 1);
previous_ripples_PRE_V1_log_odds    = nan(nUP, 1);
previous_ripples_V1_log_odds        = nan(nUP, 1);

last_ripple_next_spindle_diff  = nan(nUP, nRegions);
last_ripple_next_spindle_power = nan(nUP, nRegions);
first_ripple_spindle_power     = nan(nUP, nRegions);
last_ripple_spindle_power      = nan(nUP, nRegions);

fprintf('Extracting per-UP ripple timing, log-odds and spindle features...\n');
for iUP = 1:nUP
    upOnset  = merged_event_info.UP_ints(iUP, 1);
    upOffset = merged_event_info.UP_ints(iUP, 2);

    % overlapIdx = find(ripIntsAbs(:,1) > upOnset & ripIntsAbs(:,1) < upOffset);
    overlapIdx = find(ripPeakAbs(:,1) > upOnset & ripPeakAbs(:,1) < upOffset);

    if isempty(overlapIdx)
        continue
    end

    [~, sIdx] = sort(ripIntsAbs(overlapIdx, 1));
    ripples_index = overlapIdx(sIdx);

    first_ripples_index(iUP) = ripples_index(1);
    last_ripples_index(iUP)  = ripples_index(end);

    if strcmp(time_ref, 'peak')
        time_from_last_ripples(iUP)  = upOffset - ripPeakAbs(ripples_index(end));
        time_from_first_ripples(iUP) = upOffset - ripPeakAbs(ripples_index(1));
        time_to_last_ripples(iUP)    = ripPeakAbs(ripples_index(end)) - upOnset;
        time_to_first_ripples(iUP)   = ripPeakAbs(ripples_index(1)) - upOnset;
    else
        time_from_last_ripples(iUP)  = upOffset - ripIntsAbs(ripples_index(end),2);
        time_from_first_ripples(iUP) = upOffset - ripIntsAbs(ripples_index(1),2);
        time_to_last_ripples(iUP)    = ripIntsAbs(ripples_index(end),1) - upOnset;
        time_to_first_ripples(iUP)   = ripIntsAbs(ripples_index(1),1) - upOnset;
    end



    first_ripples_power(iUP) = ripPowerAbs(ripples_index(1));
    last_ripples_power(iUP)  = ripPowerAbs(ripples_index(end));

    first_ripples_log_odds(iUP)        = mean_z_bias(ripples_index(1));
    last_ripples_log_odds(iUP)         = mean_z_bias(ripples_index(end));
    first_ripples_PRE_log_odds(iUP)    = mean_z_bias_PRE(ripples_index(1));
    last_ripples_PRE_log_odds(iUP)     = mean_z_bias_PRE(ripples_index(end));
    first_ripples_V1_log_odds(iUP)     = mean_z_bias_V1(ripples_index(1));
    last_ripples_V1_log_odds(iUP)      = mean_z_bias_V1(ripples_index(end));
    first_ripples_PRE_V1_log_odds(iUP) = mean_z_bias_V1_PRE(ripples_index(1));
    last_ripples_PRE_V1_log_odds(iUP)  = mean_z_bias_V1_PRE(ripples_index(end));

    if length(ripples_index) > 1
        second_last_ripples_log_odds(iUP)        = mean_z_bias(ripples_index(end-1));
        second_last_ripples_PRE_V1_log_odds(iUP) = mean_z_bias_V1_PRE(ripples_index(end-1));
        second_last_ripples_V1_log_odds(iUP)     = mean_z_bias_V1(ripples_index(end-1));
    end

    % Ripple immediately before this UP's first ripple, in global ripple order
    if ripples_index(1) > 1
        previous_ripples_log_odds(iUP)        = mean_z_bias(ripples_index(1) - 1);
        previous_ripples_PRE_V1_log_odds(iUP) = mean_z_bias_V1_PRE(ripples_index(1) - 1);
        previous_ripples_V1_log_odds(iUP)     = mean_z_bias_V1(ripples_index(1) - 1);
    end

    % Spindle amplitude around the last ripple of this UP
    last_on_abs = ripIntsAbs(ripples_index(end), 1);
    [~, LFP_bin_UP_end] = min(abs(spindle_tvec - (upOffset - last_on_abs)));

    if ~isnan(UP_next_onset(iUP))
        [~, LFP_bin_next_UP] = min(abs(spindle_tvec - (UP_next_onset(iUP) - last_on_abs)));
        last_ripple_next_spindle_power(iUP,:) = squeeze(spindle_amplitude_temporal(LFP_bin_next_UP, ripples_index(end), :))';
        last_ripple_next_spindle_diff(iUP,:)  = squeeze(spindle_amplitude_temporal(LFP_bin_next_UP, ripples_index(end), :) - ...
            spindle_amplitude_temporal(LFP_bin_UP_end, ripples_index(end), :))';
    end

    first_ripple_spindle_power(iUP,:) = squeeze(spindle_amplitude_temporal(spindle_pre_bin, ripples_index(1), :))';
    last_ripple_spindle_power(iUP,:)  = squeeze(spindle_amplitude_temporal(spindle_pre_bin, ripples_index(end), :))';
end


% time_from_last_ripples(time_from_last_ripples<0)=0;
% time_from_first_ripples(time_from_first_ripples<0)=0;
% time_to_last_ripples(time_to_last_ripples<0)=0;
% time_to_first_ripples(time_to_first_ripples<0)=0;

ripple_info_table = table(...
    first_ripples_index, last_ripples_index, ...
    time_to_first_ripples, time_from_first_ripples, ...
    time_to_last_ripples, time_from_last_ripples, ...
    first_ripples_power, last_ripples_power, ...
    first_ripples_log_odds, last_ripples_log_odds, ...
    first_ripples_PRE_log_odds, last_ripples_PRE_log_odds, ...
    first_ripples_V1_log_odds, last_ripples_V1_log_odds, ...
    first_ripples_PRE_V1_log_odds, last_ripples_PRE_V1_log_odds, ...
    second_last_ripples_log_odds, second_last_ripples_PRE_V1_log_odds, second_last_ripples_V1_log_odds, ...
    previous_ripples_log_odds, previous_ripples_PRE_V1_log_odds, previous_ripples_V1_log_odds, ...
    'VariableNames', { ...
    'first_ripples_index', 'last_ripples_index', ...
    'time_to_first_ripples', 'time_from_first_ripples', ...
    'time_to_last_ripples', 'time_from_last_ripples', ...
    'first_ripples_power', 'last_ripples_power', ...
    'first_ripples_log_odds', 'last_ripples_log_odds', ...
    'first_ripples_PRE_log_odds', 'last_ripples_PRE_log_odds', ...
    'first_ripples_V1_log_odds', 'last_ripples_V1_log_odds', ...
    'first_ripples_PRE_V1_log_odds', 'last_ripples_PRE_V1_log_odds', ...
    'second_last_ripples_log_odds', 'second_last_ripples_PRE_V1_log_odds', 'second_last_ripples_V1_log_odds', ...
    'previous_ripples_log_odds', 'previous_ripples_PRE_V1_log_odds', 'previous_ripples_V1_log_odds'} ...
    );

ripple_info_table.last_ripple_next_spindle_diff  = last_ripple_next_spindle_diff;
ripple_info_table.last_ripple_next_spindle_power = last_ripple_next_spindle_power;
ripple_info_table.first_ripple_spindle_power     = first_ripple_spindle_power;
ripple_info_table.last_ripple_spindle_power      = last_ripple_spindle_power;

end
