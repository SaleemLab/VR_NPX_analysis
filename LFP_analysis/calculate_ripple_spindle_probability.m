function probability = calculate_ripple_spindle_probability(slow_waves_all, ripples_all, spindles_all, sessions_to_process, varargin)
% Calculates spindle probability relative to L/R ripple onset or peaktimes

p = inputParser;
addParameter(p,'option','absolute',@ischar);
addParameter(p,'shuffle_option','time_circular_shift',@ischar);
addParameter(p,'time_option','peaktimes',@ischar);
addParameter(p,'time_windows',[-0.5 0.5],@isnumeric);
addParameter(p,'time_bin',0.02,@isnumeric);
addParameter(p,'num_bins',20,@isnumeric);
addParameter(p,'duration_threshold',2,@isnumeric);
parse(p,varargin{:})

option = p.Results.option;
time_option = p.Results.time_option;
time_windows = p.Results.time_windows;
time_bin = p.Results.time_bin;
num_bins = p.Results.num_bins;
duration_threshold = p.Results.duration_threshold;
shuffle_option = p.Results.shuffle_option;

timebin_edge = time_windows(1):time_bin:time_windows(end);
bins_centre = timebin_edge(1)+time_bin/2:time_bin:timebin_edge(end)-time_bin/2;

for nprobe = 1:length(spindles_all)
    for hemi = 1:2 % hemi 1 = L ripple, 2 = R ripple
        ripple_index_all = [];
        binnedArray = [];

        for nsession = 1:length(sessions_to_process)
            session_id = sessions_to_process(nsession);

            ripple_index = find(ripples_all(hemi).session_count == session_id & ripples_all(hemi).SWS_index == 1);
            if isempty(ripple_index), continue; end

            ripple_onset = ripples_all(hemi).onset(ripple_index);

            if contains(shuffle_option,'baseline')
                ripple_onset = ripple_onset - 3;
            end

            % Get spindles for this probe
            spindles_index = find(spindles_all(nprobe).session_count == session_id & spindles_all(nprobe).SWS_index == 1);
            if isempty(spindles_index), continue; end

            if contains(time_option,'onset')
                spindle_times = spindles_all(nprobe).onset(spindles_index);
            else
                spindle_times = [spindles_all(nprobe).onset(spindles_index) spindles_all(nprobe).offset(spindles_index)];
            end

            [prob, temp, event_index] = calculate_event_probability(spindle_times, ripple_onset, timebin_edge, 0);

            if ~contains(shuffle_option,'baseline')
                ripple_times = [ripples_all(hemi).onset(ripple_index) ripples_all(hemi).offset(ripple_index)];
                timebin_edges_all = ripple_onset + bins_centre;

                for i = 1:size(ripple_times,1)
                    if i > 1
                        prev_offset = ripple_times(i-1,2);
                        temp(i, timebin_edges_all(i,:) <= prev_offset) = NaN;
                    end
                    if i < size(ripple_times,1)
                        next_onset = ripple_times(i+1,1);
                        temp(i, timebin_edges_all(i,:) >= next_onset) = NaN;
                    end
                end
            end

            binnedArray = [binnedArray; temp];
            ripple_index_all = [ripple_index_all; event_index];
        end

        tag = ["L","R"];
        ripple_label = tag{hemi};
        
        probability(nprobe).([ripple_label '_spindle']) = binnedArray;
        probability(nprobe).([ripple_label '_spindle_index']) = ripple_index_all;
        probability(nprobe).([ripple_label '_ripple_no']) = sum(ripples_all(hemi).SWS_index == 1);

        % Bootstrap
        tempBoot = [];
        for iBoot = 1:1000
            s = RandStream('mrg32k3a','Seed',iBoot);
            event_id = datasample(s,1:size(binnedArray,1),size(binnedArray,1));
            tempBoot(iBoot,:) = sum(binnedArray(event_id,:), 'omitnan') ./ sum(~isnan(binnedArray(event_id,:)));
        end
        probability(nprobe).([ripple_label '_spindle_bootstrap']) = tempBoot;

        % Shuffle
        if contains(shuffle_option,'time_circular_shift')
            tempShuf = [];
            for iBoot = 1:1000
                temp1 = zeros(size(binnedArray));
                parfor event_id = 1:size(binnedArray,1)
                    s = RandStream('mrg32k3a','Seed',100000*iBoot+event_id);
                    shift = datasample(s,1:size(binnedArray,2),1);
                    bins = circshift(1:size(binnedArray,2), shift);
                    temp1(event_id,:) = binnedArray(event_id,bins);
                end
                tempShuf(iBoot,:) = sum(temp1,1,'omitnan') ./ sum(~isnan(temp1));
            end
            probability(nprobe).([ripple_label '_spindle_shuffled']) = tempShuf;
        end
    end
end
end
