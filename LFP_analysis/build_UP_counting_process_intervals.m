function T = build_UP_counting_process_intervals(merged_event_info, V1_MUA_spiketimes, HC_MUA_spiketimes, varargin)
% BUILD_UP_COUNTING_PROCESS_INTERVALS Slices UP state events into counting process
% intervals [start, stop) bounded by hippocampal ripple onsets and offsets.
%
% Computes 10ms binned, session-level min to 99th percentile normalized MUA
% directly from raw spike times (V1_MUA_spiketimes and HC_MUA_spiketimes).
%
% Inputs:
%   merged_event_info - Struct containing:
%       .UP_ints              - Nx2 matrix [onset, offset] of UP states (abs time)
%       .UP_hemisphere_id     - Nx1 vector (1=Left, 2=Right)
%       .ripples_ints         - Mx2 matrix [onset, offset] of ripples (abs time)
%       .ripples_power        - Mx1 vector of LFP ripple z-score peak powers
%       .session_id           - Nx1 vector of session IDs
%       .subject_id           - (Optional) Nx1 vector of subject IDs
%   V1_MUA_spiketimes - Cell array {1} (Left V1) and {2} (Right V1) of spike times
%   HC_MUA_spiketimes - Cell array {1} (Left HPC) and {2} (Right HPC) of spike times
%
% Output:
%   T                 - MATLAB table containing interval-level counting process data.

p = inputParser;
addParameter(p, 'time_bin', 0.01, @isnumeric);         % 10 ms time bin
addParameter(p, 'min_interval_dur', 1e-5, @isnumeric); % minimum duration threshold
parse(p, varargin{:});

time_bin = p.Results.time_bin;
min_dur  = p.Results.min_interval_dur;

nUP = size(merged_event_info.UP_ints, 1);
ripIntsAbs  = merged_event_info.ripples_ints;
ripPowerAbs = merged_event_info.ripples_power;
ripPeakAbs  = merged_event_info.ripples_peaktimes;    
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
    
    % Find max timestamp for binning grid
    max_up_t = max(merged_event_info.UP_ints(session_id == s, 2) - s_offset);
    max_spk_t = max([0; spk_v1_L; spk_v1_R; spk_hc_L; spk_hc_R]);
    max_t = max(max_up_t, max_spk_t) + time_bin;
    
    t_edges = 0:time_bin:max_t;
    
    % Bin raw spike counts (10ms bins)
    cnt_v1_L = histcounts(spk_v1_L, t_edges);
    cnt_v1_R = histcounts(spk_v1_R, t_edges);
    cnt_hc_L = histcounts(spk_hc_L, t_edges);
    cnt_hc_R = histcounts(spk_hc_R, t_edges);
    
    % Min-to-99th percentile normalization per session and channel
    min_v1_L = min(cnt_v1_L); p99_v1_L = prctile(cnt_v1_L, 99); if p99_v1_L <= min_v1_L, p99_v1_L = min_v1_L + 1; end
    min_v1_R = min(cnt_v1_R); p99_v1_R = prctile(cnt_v1_R, 99); if p99_v1_R <= min_v1_R, p99_v1_R = min_v1_R + 1; end
    min_hc_L = min(cnt_hc_L); p99_hc_L = prctile(cnt_hc_L, 99); if p99_hc_L <= min_hc_L, p99_hc_L = min_hc_L + 1; end
    min_hc_R = min(cnt_hc_R); p99_hc_R = prctile(cnt_hc_R, 99); if p99_hc_R <= min_hc_R, p99_hc_R = min_hc_R + 1; end
    
    norm_v1_L{s} = min(1, max(0, (cnt_v1_L - min_v1_L) / (p99_v1_L - min_v1_L)));
    norm_v1_R{s} = min(1, max(0, (cnt_v1_R - min_v1_R) / (p99_v1_R - min_v1_R)));
    norm_hc_L{s} = min(1, max(0, (cnt_hc_L - min_hc_L) / (p99_hc_L - min_hc_L)));
    norm_hc_R{s} = min(1, max(0, (cnt_hc_R - min_hc_R) / (p99_hc_R - min_hc_R)));
end

%% 2. Pre-allocate storage for table rows
upID_list          = cell(nUP, 1);
session_id_list    = cell(nUP, 1);
subject_id_list    = cell(nUP, 1);
hemisphere_id_list = cell(nUP, 1);
start_list         = cell(nUP, 1);
stop_list          = cell(nUP, 1);
duration_list      = cell(nUP, 1);
event_list         = cell(nUP, 1);
censoring_list     = cell(nUP, 1);
inRipple_list      = cell(nUP, 1);
ripplePower_list   = cell(nUP, 1);

rippleHPC_MUA_sum_list   = cell(nUP, 1);
rippleHPC_MUA_mean_list  = cell(nUP, 1);
rippleV1_MUA_sum_list    = cell(nUP, 1);
rippleV1_MUA_mean_list   = cell(nUP, 1);

nonRippleHPC_MUA_sum_list  = cell(nUP, 1);
nonRippleHPC_MUA_mean_list = cell(nUP, 1);
nonRippleV1_MUA_sum_list   = cell(nUP, 1);
nonRippleV1_MUA_mean_list  = cell(nUP, 1);

intervalHPC_MUA_sum_list  = cell(nUP, 1);
intervalHPC_MUA_mean_list = cell(nUP, 1);
intervalV1_MUA_sum_list   = cell(nUP, 1);
intervalV1_MUA_mean_list  = cell(nUP, 1);

cumRippleHPC_MUA_list    = cell(nUP, 1);
cumNonRippleHPC_MUA_list = cell(nUP, 1);
cumTotalHPC_MUA_list     = cell(nUP, 1);

cumRippleV1_MUA_list     = cell(nUP, 1);
cumNonRippleV1_MUA_list  = cell(nUP, 1);
cumTotalV1_MUA_list      = cell(nUP, 1);

cumRippleCount_list      = cell(nUP, 1);

cumRippleHPC_MUA_incl_list    = cell(nUP, 1);
cumNonRippleHPC_MUA_incl_list = cell(nUP, 1);
cumTotalHPC_MUA_incl_list     = cell(nUP, 1);
cumRippleV1_MUA_incl_list     = cell(nUP, 1);
cumNonRippleV1_MUA_incl_list  = cell(nUP, 1);
cumTotalV1_MUA_incl_list      = cell(nUP, 1);
cumRippleCount_incl_list      = cell(nUP, 1);

hasRipple_list      = cell(nUP, 1);
numRipplesInUP_list = cell(nUP, 1);

lastRipplePower_list        = cell(nUP, 1);
lastRippleHPC_MUA_sum_list  = cell(nUP, 1);
lastRippleHPC_MUA_mean_list = cell(nUP, 1);
lastRippleV1_MUA_sum_list   = cell(nUP, 1);
lastRippleV1_MUA_mean_list  = cell(nUP, 1);
timeSinceLastRipple_list    = cell(nUP, 1);

%% 3. Slice UP State Events into Time-Varying Intervals
for iUP = 1:nUP
    upOnset  = merged_event_info.UP_ints(iUP, 1);
    upOffset = merged_event_info.UP_ints(iUP, 2);
    upDur    = upOffset - upOnset;
    hemiID   = merged_event_info.UP_hemisphere_id(iUP);
    
    sessID   = session_id(iUP);
    subjID   = subject_id(iUP);
    
    s_offset  = sessID * 1000000;
    rel_onset = upOnset - s_offset;
    rel_offset = upOffset - s_offset;
    
    % Extract 10ms normalized bin indices for this UP event
    bin_start = max(1, floor(rel_onset / time_bin) + 1);
    bin_end   = min(length(norm_v1_L{sessID}), ceil(rel_offset / time_bin));
    
    if bin_start > bin_end
        bin_end = bin_start;
    end
    
    % HPC MUA is averaged across Left and Right HPC
    hpcMUA = mean([norm_hc_L{sessID}(bin_start:bin_end); norm_hc_R{sessID}(bin_start:bin_end)], 1, 'omitnan');
    
    % V1 MUA is ipsilateral to cortex hemisphere of origin
    if hemiID == 1
        v1MUA = norm_v1_L{sessID}(bin_start:bin_end);
    else
        v1MUA = norm_v1_R{sessID}(bin_start:bin_end);
    end
    
    hpcMUA = hpcMUA(:)';
    v1MUA  = v1MUA(:)';
    
    nBins = length(hpcMUA);
    binCenters = (0:nBins-1) * time_bin + time_bin/2;
    
    % Identify ripples overlapping this UP state
    % overlapIdx = find(ripIntsAbs(:,2) > upOnset & ripIntsAbs(:,1) < upOffset);
    overlapIdx = find(ripPeakAbs(:,1) > upOnset & ripPeakAbs(:,1) < upOffset);

    ripRel = [];
    ripPow = [];
    if ~isempty(overlapIdx)
        for r = 1:length(overlapIdx)
            idx = overlapIdx(r);
            rOn  = max(0, ripIntsAbs(idx, 1) - upOnset);
            rOff = min(upDur, ripIntsAbs(idx, 2) - upOnset);
            if rOff - rOn > min_dur
                ripRel = [ripRel; rOn, rOff];
                ripPow = [ripPow; ripPowerAbs(idx)];
            end
        end
    end

    nRipples = size(ripRel, 1);
    hasRip   = double(nRipples > 0);
    
    % Generate cutpoints at 0, UP duration, and ripple onsets/offsets
    if ~isempty(ripRel)
        cuts = unique([0; ripRel(:); upDur]);
    else
        cuts = [0; upDur];
    end
    cuts = sort(cuts);
    
    % Filter cutpoints that are too close
    validCuts = cuts(1);
    for k = 2:length(cuts)
        if (cuts(k) - validCuts(end)) >= min_dur
            validCuts = [validCuts; cuts(k)];
        end
    end
    if validCuts(end) < upDur
        validCuts(end) = upDur;
    end
    cuts = validCuts;
    
    nInt = length(cuts) - 1;
    
    % Pre-allocate interval vectors for this UP event
    i_start    = zeros(nInt, 1);
    i_stop     = zeros(nInt, 1);
    i_dur      = zeros(nInt, 1);
    i_event    = zeros(nInt, 1);
    i_cens     = zeros(nInt, 1);
    i_inRip    = zeros(nInt, 1);
    i_ripPow   = zeros(nInt, 1);
    
    i_ripHpcSum   = zeros(nInt, 1);
    i_ripHpcMean  = zeros(nInt, 1);
    i_ripV1Sum    = zeros(nInt, 1);
    i_ripV1Mean   = zeros(nInt, 1);
    
    i_nonRipHpcSum  = zeros(nInt, 1);
    i_nonRipHpcMean = zeros(nInt, 1);
    i_nonRipV1Sum   = zeros(nInt, 1);
    i_nonRipV1Mean  = zeros(nInt, 1);
    
    i_intHpcSum   = zeros(nInt, 1);
    i_intHpcMean  = zeros(nInt, 1);
    i_intV1Sum    = zeros(nInt, 1);
    i_intV1Mean   = zeros(nInt, 1);
    
    i_cumRipHPC    = zeros(nInt, 1);
    i_cumNonRipHPC = zeros(nInt, 1);
    i_cumTotHPC    = zeros(nInt, 1);
    
    i_cumRipV1     = zeros(nInt, 1);
    i_cumNonRipV1  = zeros(nInt, 1);
    i_cumTotV1     = zeros(nInt, 1);
    
    i_cumRipCount  = zeros(nInt, 1);

    i_cumRipHPC_incl    = zeros(nInt, 1);
    i_cumNonRipHPC_incl = zeros(nInt, 1);
    i_cumTotHPC_incl    = zeros(nInt, 1);
    i_cumRipV1_incl     = zeros(nInt, 1);
    i_cumNonRipV1_incl  = zeros(nInt, 1);
    i_cumTotV1_incl     = zeros(nInt, 1);
    i_cumRipCount_incl  = zeros(nInt, 1);

    i_lastRipPow          = zeros(nInt, 1);
    i_lastRipHpcSum       = zeros(nInt, 1);
    i_lastRipHpcMean      = zeros(nInt, 1);
    i_lastRipV1Sum        = zeros(nInt, 1);
    i_lastRipV1Mean       = zeros(nInt, 1);
    i_timeSinceLastRipple = zeros(nInt, 1);

    runningRipHPC    = 0;
    runningNonRipHPC = 0;
    runningRipV1     = 0;
    runningNonRipV1  = 0;
    runningRipCount  = 0;

    last_rip_pow      = 0;
    last_rip_hpc_sum  = 0;
    last_rip_hpc_mean = 0;
    last_rip_v1_sum   = 0;
    last_rip_v1_mean  = 0;
    last_rip_off      = NaN;
    
    for j = 1:nInt
        t0 = cuts(j);
        t1 = cuts(j+1);
        dt = t1 - t0;
        
        i_start(j) = t0;
        i_stop(j)  = t1;
        i_dur(j)   = dt;
        
        % Check if [t0, t1) is a ripple interval
        isRip = false;
        powVal = 0;
        if ~isempty(ripRel)
            for r = 1:size(ripRel, 1)
                ovStart = max(t0, ripRel(r,1));
                ovStop  = min(t1, ripRel(r,2));
                if (ovStop - ovStart) >= 0.5 * dt || (ovStop - ovStart) >= 0.01
                    isRip = true;
                    powVal = ripPow(r);
                    break;
                end
            end
        end
        
        % Assign lagged cumulative metrics (values prior to t0)
        i_cumRipHPC(j)    = runningRipHPC;
        i_cumNonRipHPC(j) = runningNonRipHPC;
        i_cumTotHPC(j)    = runningRipHPC + runningNonRipHPC;
        
        i_cumRipV1(j)     = runningRipV1;
        i_cumNonRipV1(j)  = runningNonRipV1;
        i_cumTotV1(j)     = runningRipV1 + runningNonRipV1;
        
        i_cumRipCount(j)  = runningRipCount;

        % Identify bins falling inside [t0, t1)
        binMask = (binCenters >= t0 & binCenters < t1);
        if ~any(binMask)
            [~, nearIdx] = min(abs(binCenters - (t0 + t1)/2));
            binMask = false(1, nBins);
            binMask(nearIdx) = true;
        end
        
        hpcSub = hpcMUA(binMask);
        v1Sub  = v1MUA(binMask);
        
        hpcSumVal  = sum(hpcSub);
        hpcMeanVal = mean(hpcSub);
        v1SumVal   = sum(v1Sub);
        v1MeanVal  = mean(v1Sub);
        
        % General interval MUA
        i_intHpcSum(j)  = hpcSumVal;
        i_intHpcMean(j) = hpcMeanVal;
        i_intV1Sum(j)   = v1SumVal;
        i_intV1Mean(j)  = v1MeanVal;
        
        if isRip
            i_inRip(j)     = 1;
            i_ripPow(j)    = powVal;
            
            % Ripple-specific MUA (0 for non-ripple)
            i_ripHpcSum(j)  = hpcSumVal;
            i_ripHpcMean(j) = hpcMeanVal;
            i_ripV1Sum(j)   = v1SumVal;
            i_ripV1Mean(j)  = v1MeanVal;
            
            i_nonRipHpcSum(j)  = 0;
            i_nonRipHpcMean(j) = 0;
            i_nonRipV1Sum(j)   = 0;
            i_nonRipV1Mean(j)  = 0;
            
            runningRipHPC   = runningRipHPC + hpcSumVal;
            runningRipV1    = runningRipV1  + v1SumVal;
            runningRipCount = runningRipCount + 1;

            % Update carried-forward ripple features
            last_rip_pow      = powVal;
            last_rip_hpc_sum  = hpcSumVal;
            last_rip_hpc_mean = hpcMeanVal;
            last_rip_v1_sum   = v1SumVal;
            last_rip_v1_mean  = v1MeanVal;
            
            % Update last ripple offset
            if ~isempty(ripRel)
                for r = 1:size(ripRel, 1)
                    ovStart = max(t0, ripRel(r,1));
                    ovStop  = min(t1, ripRel(r,2));
                    if (ovStop - ovStart) >= 0.5 * dt || (ovStop - ovStart) >= 0.01
                        last_rip_off = ripRel(r, 2);
                        break;
                    end
                end
            end

            i_timeSinceLastRipple(j) = 0;
        else
            i_inRip(j)     = 0;
            i_ripPow(j)    = 0;
            
            i_ripHpcSum(j)  = 0;
            i_ripHpcMean(j) = 0;
            i_ripV1Sum(j)   = 0;
            i_ripV1Mean(j)  = 0;
            
            % Non-ripple specific MUA (0 for ripple)
            i_nonRipHpcSum(j)  = hpcSumVal;
            i_nonRipHpcMean(j) = hpcMeanVal;
            i_nonRipV1Sum(j)   = v1SumVal;
            i_nonRipV1Mean(j)  = v1MeanVal;
            
            runningNonRipHPC = runningNonRipHPC + hpcSumVal;
            runningNonRipV1  = runningNonRipV1  + v1SumVal;

            if isnan(last_rip_off)
                i_timeSinceLastRipple(j) = 0;
            else
                i_timeSinceLastRipple(j) = max(0, t1 - last_rip_off);
            end
        end
        
        % Assign inclusive cumulative metrics (values up to t1, including current interval)
        i_cumRipHPC_incl(j)    = runningRipHPC;
        i_cumNonRipHPC_incl(j) = runningNonRipHPC;
        i_cumTotHPC_incl(j)    = runningRipHPC + runningNonRipHPC;
        
        i_cumRipV1_incl(j)     = runningRipV1;
        i_cumNonRipV1_incl(j)  = runningNonRipV1;
        i_cumTotV1_incl(j)     = runningRipV1 + runningNonRipV1;
        
        i_cumRipCount_incl(j)  = runningRipCount;

        % Assign carried-forward ripple features
        i_lastRipPow(j)     = last_rip_pow;
        i_lastRipHpcSum(j)  = last_rip_hpc_sum;
        i_lastRipHpcMean(j) = last_rip_hpc_mean;
        i_lastRipV1Sum(j)   = last_rip_v1_sum;
        i_lastRipV1Mean(j)  = last_rip_v1_mean;

        % Event status: last interval ends in DOWN transition (event=1, censoring=0)
        if j == nInt
            i_event(j) = 1;
            i_cens(j)  = 0;
        else
            i_event(j) = 0;
            i_cens(j)  = 1;
        end
    end
    
    % Store interval arrays for this UP event
    upID_list{iUP}          = repmat(iUP, nInt, 1);
    session_id_list{iUP}    = repmat(sessID, nInt, 1);
    subject_id_list{iUP}    = repmat(subjID, nInt, 1);
    hemisphere_id_list{iUP} = repmat(hemiID, nInt, 1);
    start_list{iUP}         = i_start;
    stop_list{iUP}          = i_stop;
    duration_list{iUP}      = i_dur;
    event_list{iUP}         = i_event;
    censoring_list{iUP}     = i_cens;
    inRipple_list{iUP}      = i_inRip;
    ripplePower_list{iUP}   = i_ripPow;
    
    rippleHPC_MUA_sum_list{iUP}   = i_ripHpcSum;
    rippleHPC_MUA_mean_list{iUP}  = i_ripHpcMean;
    rippleV1_MUA_sum_list{iUP}    = i_ripV1Sum;
    rippleV1_MUA_mean_list{iUP}   = i_ripV1Mean;
    
    nonRippleHPC_MUA_sum_list{iUP}  = i_nonRipHpcSum;
    nonRippleHPC_MUA_mean_list{iUP} = i_nonRipHpcMean;
    nonRippleV1_MUA_sum_list{iUP}   = i_nonRipV1Sum;
    nonRippleV1_MUA_mean_list{iUP}  = i_nonRipV1Mean;
    
    intervalHPC_MUA_sum_list{iUP}  = i_intHpcSum;
    intervalHPC_MUA_mean_list{iUP} = i_intHpcMean;
    intervalV1_MUA_sum_list{iUP}   = i_intV1Sum;
    intervalV1_MUA_mean_list{iUP}  = i_intV1Mean;
    
    cumRippleHPC_MUA_list{iUP}    = i_cumRipHPC;
    cumNonRippleHPC_MUA_list{iUP} = i_cumNonRipHPC;
    cumTotalHPC_MUA_list{iUP}     = i_cumTotHPC;
    
    cumRippleV1_MUA_list{iUP}     = i_cumRipV1;
    cumNonRippleV1_MUA_list{iUP}  = i_cumNonRipV1;
    cumTotalV1_MUA_list{iUP}      = i_cumTotV1;
    
    cumRippleCount_list{iUP}      = i_cumRipCount;

    cumRippleHPC_MUA_incl_list{iUP}    = i_cumRipHPC_incl;
    cumNonRippleHPC_MUA_incl_list{iUP} = i_cumNonRipHPC_incl;
    cumTotalHPC_MUA_incl_list{iUP}     = i_cumTotHPC_incl;
    
    cumRippleV1_MUA_incl_list{iUP}     = i_cumRipV1_incl;
    cumNonRippleV1_MUA_incl_list{iUP}  = i_cumNonRipV1_incl;
    cumTotalV1_MUA_incl_list{iUP}      = i_cumTotV1_incl;
    
    cumRippleCount_incl_list{iUP}      = i_cumRipCount_incl;

    hasRipple_list{iUP}      = repmat(hasRip, nInt, 1);
    numRipplesInUP_list{iUP} = repmat(nRipples, nInt, 1);

    lastRipplePower_list{iUP}        = i_lastRipPow;
    lastRippleHPC_MUA_sum_list{iUP}  = i_lastRipHpcSum;
    lastRippleHPC_MUA_mean_list{iUP} = i_lastRipHpcMean;
    lastRippleV1_MUA_sum_list{iUP}   = i_lastRipV1Sum;
    lastRippleV1_MUA_mean_list{iUP}  = i_lastRipV1Mean;
    timeSinceLastRipple_list{iUP}    = i_timeSinceLastRipple;
end

% Concatenate all intervals into a single table
T = table(...
    vertcat(upID_list{:}), ...
    vertcat(session_id_list{:}), ...
    vertcat(subject_id_list{:}), ...
    vertcat(hemisphere_id_list{:}), ...
    vertcat(start_list{:}), ...
    vertcat(stop_list{:}), ...
    vertcat(duration_list{:}), ...
    vertcat(event_list{:}), ...
    vertcat(censoring_list{:}), ...
    vertcat(inRipple_list{:}), ...
    vertcat(ripplePower_list{:}), ...
    vertcat(rippleHPC_MUA_sum_list{:}), ...
    vertcat(rippleHPC_MUA_mean_list{:}), ...
    vertcat(rippleV1_MUA_sum_list{:}), ...
    vertcat(rippleV1_MUA_mean_list{:}), ...
    vertcat(nonRippleHPC_MUA_sum_list{:}), ...
    vertcat(nonRippleHPC_MUA_mean_list{:}), ...
    vertcat(nonRippleV1_MUA_sum_list{:}), ...
    vertcat(nonRippleV1_MUA_mean_list{:}), ...
    vertcat(intervalHPC_MUA_sum_list{:}), ...
    vertcat(intervalHPC_MUA_mean_list{:}), ...
    vertcat(intervalV1_MUA_sum_list{:}), ...
    vertcat(intervalV1_MUA_mean_list{:}), ...
    vertcat(cumRippleHPC_MUA_list{:}), ...
    vertcat(cumNonRippleHPC_MUA_list{:}), ...
    vertcat(cumTotalHPC_MUA_list{:}), ...
    vertcat(cumRippleV1_MUA_list{:}), ...
    vertcat(cumNonRippleV1_MUA_list{:}), ...
    vertcat(cumTotalV1_MUA_list{:}), ...
    vertcat(cumRippleCount_list{:}), ...
    vertcat(cumRippleHPC_MUA_incl_list{:}), ...
    vertcat(cumNonRippleHPC_MUA_incl_list{:}), ...
    vertcat(cumTotalHPC_MUA_incl_list{:}), ...
    vertcat(cumRippleV1_MUA_incl_list{:}), ...
    vertcat(cumNonRippleV1_MUA_incl_list{:}), ...
    vertcat(cumTotalV1_MUA_incl_list{:}), ...
    vertcat(cumRippleCount_incl_list{:}), ...
    vertcat(hasRipple_list{:}), ...
    vertcat(numRipplesInUP_list{:}), ...
    vertcat(lastRipplePower_list{:}), ...
    vertcat(lastRippleHPC_MUA_sum_list{:}), ...
    vertcat(lastRippleHPC_MUA_mean_list{:}), ...
    vertcat(lastRippleV1_MUA_sum_list{:}), ...
    vertcat(lastRippleV1_MUA_mean_list{:}), ...
    vertcat(timeSinceLastRipple_list{:}), ...
    'VariableNames', { ...
    'upID', 'session_id', 'subject_id', 'hemisphere_id', ...
    'start', 'stop', 'duration', 'event', 'censoring', ...
    'inRipple', 'ripplePower', ...
    'rippleHPC_MUA_sum', 'rippleHPC_MUA_mean', ...
    'rippleV1_MUA_sum', 'rippleV1_MUA_mean', ...
    'nonRippleHPC_MUA_sum', 'nonRippleHPC_MUA_mean', ...
    'nonRippleV1_MUA_sum', 'nonRippleV1_MUA_mean', ...
    'intervalHPC_MUA_sum', 'intervalHPC_MUA_mean', ...
    'intervalV1_MUA_sum', 'intervalV1_MUA_mean', ...
    'cumRippleHPC_MUA', 'cumNonRippleHPC_MUA', 'cumTotalHPC_MUA', ...
    'cumRippleV1_MUA', 'cumNonRippleV1_MUA', 'cumTotalV1_MUA', ...
    'cumRippleCount', ...
    'cumRippleHPC_MUA_incl', 'cumNonRippleHPC_MUA_incl', 'cumTotalHPC_MUA_incl', ...
    'cumRippleV1_MUA_incl', 'cumNonRippleV1_MUA_incl', 'cumTotalV1_MUA_incl', ...
    'cumRippleCount_incl', ...
    'hasRipple', 'numRipplesInUP', ...
    'lastRipplePower', 'lastRippleHPC_MUA_sum', 'lastRippleHPC_MUA_mean', ...
    'lastRippleV1_MUA_sum', 'lastRippleV1_MUA_mean', 'timeSinceLastRipple'} ...
);

end
