function probability = calculate_UP_DOWN_ripple_probability(slow_waves_all,ripples_all,sessions_to_process,varargin)

p = inputParser;
addParameter(p,'option','absolute',@ischar);

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


for nprobe = 1:length(slow_waves_all)
    %%%%%%%%%%%%%%% L ripples
    UP_index_all = [];
    DOWN_index_all = [];
    ripple_index_UP_all = [];
    ripple_index_DOWN_all = [];

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

        ripples_index = find(ripples_all(1).session_count == sessions_to_process(nsession)& ripples_all(1).SWS_index == 1);
        ripple_peaktimes = min(ripples_all(1).SWR_peaktimes{sessions_to_process(nsession)}(ripples_all(1).probe_hemisphere{sessions_to_process(nsession)} == 2,ripples_all(1).SWS_index(ripples_all(1).session_count == sessions_to_process(nsession))==1))';

        if contains(time_option,'peaktimes')
            ripple_times= ripple_peaktimes;
        else
            ripple_times = [ripples_all(1).onset(ripples_index) ripples_all(1).offset(ripples_index)];
        end

        if contains(option,'normalised')
            % Ripple probability during normalised UP duration
            [probability(nprobe).L_ripples_DOWN_session(nsession,:),event_index,normalized_duration,temp] = calculate_relative_event_probability(slow_waves_all(nprobe).DOWN_ints(DOWN_index,:),ripple_times,num_bins,0);
        else
            [probability(nprobe).L_ripples_DOWN_session(nsession,:),temp,event_index] = calculate_event_probability(ripple_times,slow_waves_all(nprobe).DOWN_ints(DOWN_index,1),time_wondows(1):time_bin:time_wondows(end),0);
        end

        binnedArrayDOWN=[binnedArrayDOWN; temp];
        ripple_index_DOWN_all = [ripple_index_DOWN_all;event_index];

        if contains(option,'normalised')
            [probability(nprobe).L_ripples_UP_session(nsession,:),event_index,normalized_duration,temp] = calculate_relative_event_probability(UP_ints,ripple_times,num_bins,0);
        else
            [probability(nprobe).L_ripples_UP_session(nsession,:),temp,event_index] = calculate_event_probability(ripple_times,slow_waves_all(nprobe).UP_ints(UP_index,1),time_wondows(1):time_bin:time_wondows(end),0);
        end

        binnedArrayUP=[binnedArrayUP; temp];
        ripple_index_UP_all = [ripple_index_UP_all;event_index];
        
    end
    probability(nprobe).UP_all_index = UP_index_all;
    probability(nprobe).DOWN_all_index = DOWN_index_all;
    probability(nprobe).UP_duration = slow_waves_all(nprobe).UP_ints(UP_index_all,2)-slow_waves_all(nprobe).UP_ints(UP_index_all,1);
    probability(nprobe).DOWN_duration = slow_waves_all(nprobe).DOWN_ints(DOWN_index_all,2)-slow_waves_all(nprobe).DOWN_ints(DOWN_index_all,1);

    probability(nprobe).L_ripples_DOWN = binnedArrayDOWN;
    probability(nprobe).L_ripples_UP = binnedArrayUP;
    probability(nprobe).L_ripples_DOWN_index = ripple_index_DOWN_all;
    probability(nprobe).L_ripples_UP_index = ripple_index_UP_all;
    
    all_ripple_no = sum(ripples_all(1).SWS_index == 1);
    probability(nprobe).L_ripple_no = all_ripple_no;
    
    % bootstrap distribution
    tempUP = [];
    tempDOWN = [];
    for iBoot = 1:1000
        s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
        event_id = datasample(s,1:size(binnedArrayUP,1),size(binnedArrayUP,1));
        tempUP(iBoot,:) = sum(binnedArrayUP(event_id,:))/all_ripple_no;

        s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
        event_id = datasample(s,1:size(binnedArrayDOWN,1),size(binnedArrayDOWN,1));
        tempDOWN(iBoot,:) = sum(binnedArrayDOWN(event_id,:))/all_ripple_no;
    end

    probability(nprobe).L_ripples_DOWN_bootstrap = tempDOWN;
    probability(nprobe).L_ripples_UP_bootstrap = tempUP;

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

        tempUP(iBoot,:) = sum(temp1,1)/all_ripple_no;
        tempDOWN(iBoot,:) = sum(temp2,1)/all_ripple_no;
    end
    toc
    
    probability(nprobe).L_ripples_DOWN_shuffled = tempDOWN;
    probability(nprobe).L_ripples_UP_shuffled = tempUP;




    %%%%%%%%%%%%%%% R ripples
    UP_index_all = [];
    DOWN_index_all = [];

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

        
        ripples_index = find(ripples_all(2).session_count == sessions_to_process(nsession)& ripples_all(2).SWS_index == 1);
        ripple_peaktimes = min(ripples_all(2).SWR_peaktimes{sessions_to_process(nsession)}(ripples_all(2).probe_hemisphere{sessions_to_process(nsession)} == 2,ripples_all(2).SWS_index(ripples_all(2).session_count == sessions_to_process(nsession))==1))';
        
        if contains(time_option,'peaktimes')
            ripple_times= ripple_peaktimes;
        else
            ripple_times = [ripples_all(2).onset(ripples_index) ripples_all(2).offset(ripples_index)];
        end

        if contains(option,'normalised')
            % Ripple probability during normalised UP duration
            [probability(nprobe).R_ripples_DOWN_session(nsession,:),event_index,normalized_duration,temp] = calculate_relative_event_probability(slow_waves_all(nprobe).DOWN_ints(DOWN_index,:),ripple_times,num_bins,0);
        else
            [probability(nprobe).R_ripples_DOWN_session(nsession,:),temp,event_index] = calculate_event_probability(ripple_times,slow_waves_all(nprobe).DOWN_ints(DOWN_index,1),time_wondows(1):time_bin:time_wondows(end),0);
        end

        binnedArrayDOWN=[binnedArrayDOWN; temp];
        ripple_index_DOWN_all = [ripple_index_DOWN_all;event_index];

        if contains(option,'normalised')
            [probability(nprobe).R_ripples_UP_session(nsession,:),event_index,normalized_duration,temp] = calculate_relative_event_probability(UP_ints,ripple_times,num_bins,0);
        else
            [probability(nprobe).R_ripples_UP_session(nsession,:),temp,event_index] = calculate_event_probability(ripple_times,slow_waves_all(nprobe).UP_ints(UP_index,1),time_wondows(1):time_bin:time_wondows(end),0);
        end

        binnedArrayUP=[binnedArrayUP; temp];
        ripple_index_UP_all = [ripple_index_UP_all;event_index];
        
    end

    probability(nprobe).R_ripples_DOWN = binnedArrayDOWN;
    probability(nprobe).R_ripples_UP = binnedArrayUP;
    probability(nprobe).R_ripples_DOWN_index = ripple_index_DOWN_all;
    probability(nprobe).R_ripples_UP_index = ripple_index_UP_all;

    all_ripple_no = sum(ripples_all(2).SWS_index == 1);
    probability(nprobe).R_ripple_no = all_ripple_no;
    
    % bootstrap distribution
    tempUP = [];
    tempDOWN = [];
    for iBoot = 1:1000
        s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
        event_id = datasample(s,1:size(binnedArrayUP,1),size(binnedArrayUP,1));
        tempUP(iBoot,:) = sum(binnedArrayUP(event_id,:))/all_ripple_no;

        s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
        event_id = datasample(s,1:size(binnedArrayDOWN,1),size(binnedArrayDOWN,1));
        tempDOWN(iBoot,:) = sum(binnedArrayDOWN(event_id,:))/all_ripple_no;
    end

    probability(nprobe).R_ripples_DOWN_bootstrap = tempDOWN;
    probability(nprobe).R_ripples_UP_bootstrap = tempUP;

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

        tempUP(iBoot,:) = sum(temp1,1)/all_ripple_no;
        tempDOWN(iBoot,:) = sum(temp2,1)/all_ripple_no;
    end
    toc
    
    probability(nprobe).R_ripples_DOWN_shuffled = tempDOWN;
    probability(nprobe).R_ripples_UP_shuffled = tempUP;
end