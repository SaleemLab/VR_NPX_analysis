function extract_UP_DOWN_ripples_info(slow_waves_all,ripples_all,behavioural_state_merged_all,sessions_to_process,varargin)

p = inputParser;
addParameter(p,'option','UD',@ischar);

% addParameter(p,'time_option','peaktimes',@ischar);
addParameter(p,'time_wondows',[-0.5 0.5],@isnumeric);
addParameter(p,'time_bin',0.02,@isnumeric);
addParameter(p,'num_bins',20,@isnumeric);
addParameter(p,'duration_threshold',2,@isnumeric);

parse(p,varargin{:})
option = p.Results.option;
% time_option = p.Results.time_option;
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

        if contains(options,'UD')
            [C,ia,ib] = intersect(UP_ints(:,2),DOWN_ints(:,1));
        elseif contains(options,'DU')
            [C,ia,ib] = intersect(UP_ints(:,2),DOWN_ints(:,1));
        end
        
        % [C,ia,ib] = intersect(UP_ints(:,2)+0.01,DOWN_ints(:,1));
        UP_index = UP_index(ia);
        DOWN_index = DOWN_index(ib);
        DOWN_index_all = [DOWN_index_all; DOWN_index];

        UP_ints = slow_waves_all(nprobe).UP_ints(UP_index,:);
        UP_index_all = [UP_index_all; UP_index];

        ripples_index = find(ripples_all(1).session_count == sessions_to_process(nsession)& ripples_all(1).SWS_index == 1);
        ripple_peaktimes = min(ripples_all(1).SWR_peaktimes{sessions_to_process(nsession)}(ripples_all(1).probe_hemisphere{sessions_to_process(nsession)} == 2,ripples_all(1).SWS_index(ripples_all(1).session_count == sessions_to_process(nsession))==1))';
        ripple_times = [ripples_all(1).onset(ripples_index) ripples_all(1).offset(ripples_index)];

        calculate_relative_event_probability(slow_waves_all(nprobe).DOWN_ints(DOWN_index,:),ripple_times,num_bins,0);



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

end