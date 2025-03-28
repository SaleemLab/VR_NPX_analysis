function probability = calculate_UP_DOWN_spindle_probability(slow_waves_all,spindles_all,sessions_to_process,varargin)
% ripples_all = [];
p = inputParser;
addParameter(p,'option','absolute',@ischar);
addParameter(p,'shuffle_option','time_circular_shift',@ischar);

addParameter(p,'time_option','peaktimes',@ischar);
addParameter(p,'time_wondows',[-0.5 0.5],@isnumeric);
addParameter(p,'time_bin',0.02,@isnumeric);
addParameter(p,'num_bins',20,@isnumeric);
addParameter(p,'duration_threshold',2,@isnumeric);

parse(p,varargin{:})
option = p.Results.option;
time_option = p.Results.time_option;
time_wondows = p.Results.time_wondows;
time_bin = p.Results.time_bin;
num_bins = p.Results.num_bins;
duration_threshold = p.Results.duration_threshold;
shuffle_option = p.Results.shuffle_option;

for nprobe = 1:length(slow_waves_all)
    %%%%%%%%%%%%%%% L spindles
    UP_index_all = [];
    DOWN_index_all = [];
    spindle_index_UP_all = [];
    spindle_index_DOWN_all = [];

    binnedArrayUP = [];
    binnedArrayDOWN = [];
    binnedArrayUPShuffled = [];
    binnedArrayDOWNShuffled = [];

    for nsession = 1:length(sessions_to_process)

        % Find UP followed by a DOWN
        %         UP_DOWN_index = find(slow_waves_all(nprobe).UP_session_count == nsession); % Find UP -> DOWN
        UP_index = find(slow_waves_all(nprobe).UP_session_count == sessions_to_process(nsession)); % Find UP this session
        %         UP_index = UP_index(slow_waves_all(nprobe).UP_DOWN_index(UP_DOWN_index,1)); % Find UP followed by a DOWN
        UP_duration = slow_waves_all(nprobe).UP_ints(UP_index,2)-slow_waves_all(nprobe).UP_ints(UP_index,1);
        UP_index = UP_index(UP_duration<=duration_threshold);
        UP_ints = slow_waves_all(nprobe).UP_ints(UP_index,:);

        % Find DOWN index
        DOWN_index = find(slow_waves_all(nprobe).DOWN_session_count == sessions_to_process(nsession)); % Find DOWN this session
        DOWN_ints = slow_waves_all(nprobe).DOWN_ints(DOWN_index,:);
        [C,ia,ib] = intersect(UP_ints(:,2),DOWN_ints(:,1));
        % [C,ia,ib] = intersect(UP_ints(:,2)+0.01,DOWN_ints(:,1));
        UP_index = UP_index(ia);
        DOWN_index = DOWN_index(ib);
        DOWN_index_all = [DOWN_index_all; DOWN_index];

        UP_ints = slow_waves_all(nprobe).UP_ints(UP_index,:);
        UP_index_all = [UP_index_all; UP_index];

        if contains(shuffle_option,'baseline')
            s = RandStream('mrg32k3a','Seed',1); % Set random seed for resampling
            time_jitter = 2 + (2.5 - 2) * rand(s,1, length(UP_index));
            time_jitter = [time_jitter' time_jitter'];
            UP_ints = slow_waves_all(nprobe).UP_ints(UP_index,:)-time_jitter;

            s = RandStream('mrg32k3a','Seed',2); % Set random seed for resampling
            time_jitter = 2 + (2.5 - 2) * rand(s,1, length(DOWN_index));
            time_jitter = [time_jitter' time_jitter'];
            DOWN_ints = slow_waves_all(nprobe).DOWN_ints(DOWN_index,:)-time_jitter;
        else
            UP_ints = slow_waves_all(nprobe).UP_ints(UP_index,:);
            DOWN_ints = slow_waves_all(nprobe).DOWN_ints(DOWN_index,:);
        end

        spindles_index = find(spindles_all(1).session_count == sessions_to_process(nsession)& spindles_all(1).SWS_index == 1);
        spindle_peaktimes = spindles_all(1).SWS_peaktimes;

        if contains(time_option,'peaktimes')
            spindle_times= spindle_peaktimes;
        else
            spindle_times = [spindles_all(1).onset(spindles_index) spindles_all(1).offset(spindles_index)];
        end

        if contains(option,'normalised')
            % spindle probability during normalised UP duration
            [probability(nprobe).L_spindles_DOWN_session(nsession,:),event_index,normalized_duration,temp] = calculate_relative_event_probability(DOWN_ints,spindle_times,num_bins,0);
        else
            [probability(nprobe).L_spindles_DOWN_session(nsession,:),temp,event_index] = calculate_event_probability(spindle_times,DOWN_ints(:,1),time_wondows(1):time_bin:time_wondows(end),0);
        end

        binnedArrayDOWN=[binnedArrayDOWN; temp];
        spindle_index_DOWN_all = [spindle_index_DOWN_all;event_index];

        if contains(option,'normalised')
            [probability(nprobe).L_spindles_UP_session(nsession,:),event_index,normalized_duration,temp] = calculate_relative_event_probability(UP_ints,spindle_times,num_bins,0);
        else
            [probability(nprobe).L_spindles_UP_session(nsession,:),temp,event_index] = calculate_event_probability(spindle_times,UP_ints(:,1),time_wondows(1):time_bin:time_wondows(end),0);
        end

        binnedArrayUP=[binnedArrayUP; temp];
        spindle_index_UP_all = [spindle_index_UP_all;event_index];

    end
    probability(nprobe).UP_all_index = UP_index_all;
    probability(nprobe).DOWN_all_index = DOWN_index_all;
    probability(nprobe).UP_duration = slow_waves_all(nprobe).UP_ints(UP_index_all,2)-slow_waves_all(nprobe).UP_ints(UP_index_all,1);
    probability(nprobe).DOWN_duration = slow_waves_all(nprobe).DOWN_ints(DOWN_index_all,2)-slow_waves_all(nprobe).DOWN_ints(DOWN_index_all,1);

    probability(nprobe).L_spindles_DOWN = binnedArrayDOWN;
    probability(nprobe).L_spindles_UP = binnedArrayUP;
    probability(nprobe).L_spindles_DOWN_index = spindle_index_DOWN_all;
    probability(nprobe).L_spindles_UP_index = spindle_index_UP_all;

    all_spindle_no = sum(spindles_all(1).SWS_index == 1);
    probability(nprobe).L_spindle_no = all_spindle_no;
    all_UP_no = length(UP_index_all);
    all_DOWN_no = length(DOWN_index_all);

    % bootstrap distribution
    tempUP = [];
    tempDOWN = [];
    for iBoot = 1:1000
        s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
        event_id = datasample(s,1:size(binnedArrayUP,1),size(binnedArrayUP,1));
        tempUP(iBoot,:) = sum(binnedArrayUP(event_id,:))/all_UP_no;

        s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
        event_id = datasample(s,1:size(binnedArrayDOWN,1),size(binnedArrayDOWN,1));
        tempDOWN(iBoot,:) = sum(binnedArrayDOWN(event_id,:))/all_DOWN_no;
    end

    probability(nprobe).L_spindles_DOWN_bootstrap = tempDOWN;
    probability(nprobe).L_spindles_UP_bootstrap = tempUP;

    if contains(shuffle_option,'time_circular_shift')
        % timebin circularly shifted
        tempUP = [];
        tempDOWN = [];
        tic
        for iBoot = 1:1000
            temp1 = [];
            temp2 = [];
            parfor event_id = 1:size(binnedArrayUP,1)
                s = RandStream('mrg32k3a','Seed',100000*iBoot+event_id); % Set random seed for resampling
                %             s = RandStream('mrg32k3a','Seed',i+10000*shuffle_options); % Set random seed for resampling
                bins_to_shift = datasample(s,1:1:size(binnedArrayUP,2),1);
                bins= circshift(1:1:size(binnedArrayUP,2),bins_to_shift);
                temp1(event_id,:) = binnedArrayUP(event_id,bins);
                temp2(event_id,:) = binnedArrayDOWN(event_id,bins);
            end

            tempUP(iBoot,:) = sum(temp1,1)/all_UP_no;
            tempDOWN(iBoot,:) = sum(temp2,1)/all_DOWN_no;
        end
        toc

        probability(nprobe).L_spindles_DOWN_shuffled = tempDOWN;
        probability(nprobe).L_spindles_UP_shuffled = tempUP;

    end


    %%%%%%%%%%%%%%% R spindles
    UP_index_all = [];
    DOWN_index_all = [];
    spindle_index_UP_all = [];
    spindle_index_DOWN_all = [];

    binnedArrayUP = [];
    binnedArrayDOWN = [];
    binnedArrayUPShuffled = [];
    binnedArrayDOWNShuffled = [];

    for nsession = 1:length(sessions_to_process)
        % Find UP followed by a DOWN
        UP_index = find(slow_waves_all(nprobe).UP_session_count == sessions_to_process(nsession)); % Find UP this session
        %         UP_index = UP_index(slow_waves_all(nprobe).UP_DOWN_index(UP_DOWN_index,1)); % Find UP followed by a DOWN
        UP_duration = slow_waves_all(nprobe).UP_ints(UP_index,2)-slow_waves_all(nprobe).UP_ints(UP_index,1);
        UP_index = UP_index(UP_duration<=duration_threshold);
        UP_ints = slow_waves_all(nprobe).UP_ints(UP_index,:);

        % Find DOWN index
        DOWN_index = find(slow_waves_all(nprobe).DOWN_session_count == sessions_to_process(nsession)); % Find DOWN this session
        DOWN_ints = slow_waves_all(nprobe).DOWN_ints(DOWN_index,:);
        [C,ia,ib] = intersect(UP_ints(:,2),DOWN_ints(:,1));
        % [C,ia,ib] = intersect(UP_ints(:,2)+0.01,DOWN_ints(:,1));
        UP_index = UP_index(ia);
        DOWN_index = DOWN_index(ib);
        DOWN_index_all = [DOWN_index_all; DOWN_index];

        UP_ints = slow_waves_all(nprobe).UP_ints(UP_index,:);
        UP_index_all = [UP_index_all; UP_index];

        if contains(shuffle_option,'baseline')
            s = RandStream('mrg32k3a','Seed',1); % Set random seed for resampling
            time_jitter = 2 + (2.5 - 2) * rand(s,1, length(UP_index));
            time_jitter = [time_jitter' time_jitter'];
            UP_ints = slow_waves_all(nprobe).UP_ints(UP_index,:)-time_jitter;

            s = RandStream('mrg32k3a','Seed',2); % Set random seed for resampling
            time_jitter = 2 + (2.5 - 2) * rand(s,1, length(DOWN_index));
            time_jitter = [time_jitter' time_jitter'];
            DOWN_ints = slow_waves_all(nprobe).DOWN_ints(DOWN_index,:)-time_jitter;
        else
            UP_ints = slow_waves_all(nprobe).UP_ints(UP_index,:);
            DOWN_ints = slow_waves_all(nprobe).DOWN_ints(DOWN_index,:);
        end

        spindles_index = find(spindles_all(2).session_count == sessions_to_process(nsession)& spindles_all(2).SWS_index == 1);
        spindle_peaktimes = spindles_all(2).SWS_peaktimes;

        if contains(time_option,'peaktimes')
            spindle_times= spindle_peaktimes;
        else
            spindle_times = [spindles_all(2).onset(spindles_index) spindles_all(2).offset(spindles_index)];
        end

        if contains(option,'normalised')
            % Spindles probability during normalised UP duration
            [probability(nprobe).R_spindles_DOWN_session(nsession,:),event_index,normalized_duration,temp] = calculate_relative_event_probability(DOWN_ints,spindle_times,num_bins,0);
        else
            [probability(nprobe).R_spindles_DOWN_session(nsession,:),temp,event_index] = calculate_event_probability(spindle_times,DOWN_ints(:,1),time_wondows(1):time_bin:time_wondows(end),0);
        end

        binnedArrayDOWN=[binnedArrayDOWN; temp];
        spindle_index_DOWN_all = [spindle_index_DOWN_all;event_index];

        if contains(option,'normalised')
            [probability(nprobe).R_spindles_UP_session(nsession,:),event_index,normalized_duration,temp] = calculate_relative_event_probability(UP_ints,spindle_times,num_bins,0);
        else
            [probability(nprobe).R_spindles_UP_session(nsession,:),temp,event_index] = calculate_event_probability(spindle_times,UP_ints(:,1),time_wondows(1):time_bin:time_wondows(end),0);
        end

        binnedArrayUP=[binnedArrayUP; temp];
        spindle_index_UP_all = [spindle_index_UP_all;event_index];

    end

    probability(nprobe).R_spindles_DOWN = binnedArrayDOWN;
    probability(nprobe).R_spindles_UP = binnedArrayUP;
    probability(nprobe).R_spindles_DOWN_index = spindle_index_DOWN_all;
    probability(nprobe).R_spindles_UP_index = spindle_index_UP_all;

    all_spindle_no = sum(spindles_all(2).SWS_index == 1);
    probability(nprobe).R_spindle_no = all_spindle_no;
    all_UP_no = length(UP_index_all);
    all_DOWN_no = length(DOWN_index_all);

    % bootstrap distribution
    tempUP = [];
    tempDOWN = [];
    for iBoot = 1:1000
        s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
        event_id = datasample(s,1:size(binnedArrayUP,1),size(binnedArrayUP,1));
        tempUP(iBoot,:) = sum(binnedArrayUP(event_id,:))/all_UP_no;

        s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
        event_id = datasample(s,1:size(binnedArrayDOWN,1),size(binnedArrayDOWN,1));
        tempDOWN(iBoot,:) = sum(binnedArrayDOWN(event_id,:))/all_DOWN_no;
    end

    probability(nprobe).R_spindles_DOWN_bootstrap = tempDOWN;
    probability(nprobe).R_spindles_UP_bootstrap = tempUP;

    if contains(shuffle_option,'time_circular_shift')
        % timebin circularly shifted
        tempUP = [];
        tempDOWN = [];
        tic
        for iBoot = 1:1000
            temp1 = [];
            temp2 = [];
            parfor event_id = 1:size(binnedArrayUP,1)
                s = RandStream('mrg32k3a','Seed',100000*iBoot+event_id); % Set random seed for resampling
                %             s = RandStream('mrg32k3a','Seed',i+10000*shuffle_options); % Set random seed for resampling
                bins_to_shift = datasample(s,1:1:size(binnedArrayUP,2),1);
                bins= circshift(1:1:size(binnedArrayUP,2),bins_to_shift);
                temp1(event_id,:) = binnedArrayUP(event_id,bins);
                temp2(event_id,:) = binnedArrayDOWN(event_id,bins);
            end

            tempUP(iBoot,:) = sum(temp1,1)/all_UP_no;
            tempDOWN(iBoot,:) = sum(temp2,1)/all_DOWN_no;
        end
        toc
    end
    probability(nprobe).R_spindles_DOWN_shuffled = tempDOWN;
    probability(nprobe).R_spindles_UP_shuffled = tempUP;
end