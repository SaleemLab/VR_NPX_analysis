function probability = calculate_spindle_ripple_probability(slow_waves_all, ripples_all, spindles_all, sessions_to_process, varargin)
% Calculates ripple probability relative to L/R spindle onset or peaktimes

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

for hemi = 1:2 % hemi 1 = L spindle, 2 = R spindle
    L_binnedArray = [];
    R_binnedArray = [];
    L_ripple_index_all = [];
    R_ripple_index_all = [];

    for nsession = 1:length(sessions_to_process)
        session_id = sessions_to_process(nsession);

        spindle_index = find(spindles_all(hemi).session_count == session_id & spindles_all(hemi).SWS_index == 1);
        if isempty(spindle_index), continue; end

        spindle_times = spindles_all(hemi).onset(spindle_index);

        if contains(shuffle_option,'baseline')
            spindle_times = spindle_times - 3;
        end

        for ripple_hemi = 1:2 % ripple source: 1 = L ripple, 2 = R ripple
            ripples_index = find(ripples_all(ripple_hemi).session_count == session_id & ripples_all(ripple_hemi).SWS_index == 1);
            if isempty(ripples_index), continue; end

            if contains(time_option,'onset')
                ripple_times = ripples_all(ripple_hemi).peaktimes(ripples_index);
            elseif contains(time_option,'peaktimes')
                ripple_times = ripples_all(ripple_hemi).onset(ripples_index);
            else
                ripple_times = [ripples_all(ripple_hemi).onset(ripples_index) ripples_all(ripple_hemi).offset(ripples_index)];
            end

            [prob, temp, event_index] = calculate_event_probability(ripple_times, spindle_times, timebin_edge, 0);

            if ~contains(shuffle_option,'baseline')
                spindle_win = [spindles_all(hemi).onset(spindle_index) spindles_all(hemi).offset(spindle_index)];
                timebin_edges_all = spindle_times + bins_centre;

                for i = 1:size(spindle_win,1)
                    if i > 1
                        prev_offset = spindle_win(i-1,2);
                        temp(i, timebin_edges_all(i,:) <= prev_offset) = NaN;
                    end
                    if i < size(spindle_win,1)
                        next_onset = spindle_win(i+1,1);
                        temp(i, timebin_edges_all(i,:) >= next_onset) = NaN;
                    end
                end
            end

            if ripple_hemi == 1
                L_binnedArray = [L_binnedArray; temp];
                L_ripple_index_all = [L_ripple_index_all; event_index];
            else
                R_binnedArray = [R_binnedArray; temp];
                R_ripple_index_all = [R_ripple_index_all; event_index];
            end
        end
    end

    tag = ["L","R"];
    spindle_label = tag{hemi};

    probability(hemi).L_ripple = L_binnedArray;
    probability(hemi).R_ripple = R_binnedArray;
    probability(hemi).L_ripple_index = L_ripple_index_all;
    probability(hemi).R_ripple_index = R_ripple_index_all;
    probability(hemi).spindle_no = sum(spindles_all(hemi).SWS_index == 1);

    % Bootstrap
    for s = {'L', 'R'}
        side = s{1};
        tempArray = probability(hemi).([side '_ripple']);
        tempBoot = [];
        for iBoot = 1:1000
            s = RandStream('mrg32k3a','Seed',iBoot);
            event_id = datasample(s,1:size(tempArray,1),size(tempArray,1));
            tempBoot(iBoot,:) = sum(tempArray(event_id,:), 'omitnan') ./ sum(~isnan(tempArray(event_id,:)));
        end
        probability(hemi).([side '_ripple_bootstrap']) = tempBoot;
    end

    % Shuffle
    if contains(shuffle_option,'time_circular_shift')
        for s = {'L', 'R'}
            side = s{1};
            tempArray = probability(hemi).([side '_ripple']);
            tempShuf = [];
            for iBoot = 1:1000
                temp1 = zeros(size(tempArray));
                parfor event_id = 1:size(tempArray,1)
                    s = RandStream('mrg32k3a','Seed',100000*iBoot+event_id);
                    shift = datasample(s,1:size(tempArray,2),1);
                    bins = circshift(1:size(tempArray,2), shift);
                    temp1(event_id,:) = tempArray(event_id,bins);
                end
                tempShuf(iBoot,:) = sum(temp1,1,'omitnan') ./ sum(~isnan(temp1));
            end
            probability(hemi).([side '_ripple_shuffled']) = tempShuf;
        end
    end
end
end
