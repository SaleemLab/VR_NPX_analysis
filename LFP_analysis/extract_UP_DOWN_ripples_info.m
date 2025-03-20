function extract_UP_DOWN_ripples_info(slow_waves_all,ripples_all,behavioural_state_merged_all,sessions_to_process,varargin)

p = inputParser;
addParameter(p,'option','UD',@ischar);

% addParameter(p,'time_option','peaktimes',@ischar);
addParameter(p,'time_wondows',[-1 1],@isnumeric);
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

    event_info(nprobe).UP_index = [];
    event_info(nprobe).DOWN_index = [];
    event_info(nprobe).UP_duration = [];
    event_info(nprobe).previous_DOWN_duration = [];
    event_info(nprobe).next_DOWN_duration = [];

    event_info(nprobe).DU_slope_HPC = [];
    event_info(nprobe).DU_slope_V1 = [];

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

        if contains(option,'UD')
            [C,ia,ib] = intersect(UP_ints(:,2),DOWN_ints(:,1));
        elseif contains(option,'DU')
            [C,ia,ib] = intersect(UP_ints(:,1),DOWN_ints(:,2));
        end
        
        % [C,ia,ib] = intersect(UP_ints(:,2)+0.01,DOWN_ints(:,1));
        UP_index = UP_index(ia);
        DOWN_index = DOWN_index(ib);
        DOWN_index_all = [DOWN_index_all; DOWN_index];

        UP_ints = slow_waves_all(nprobe).UP_ints(UP_index,:);
        UP_index_all = [UP_index_all; UP_index];

        ripples_index1 = find(ripples_all(1).session_count == sessions_to_process(nsession)& ripples_all(1).SWS_index == 1);
        ripple_peaktimes1 = min(ripples_all(1).SWR_peaktimes{sessions_to_process(nsession)}(ripples_all(1).probe_hemisphere{sessions_to_process(nsession)} == 2,ripples_all(1).SWS_index(ripples_all(1).session_count == sessions_to_process(nsession))==1))';
        ripple_times1 = [ripples_all(1).onset(ripples_index1) ripples_all(1).offset(ripples_index1)];

        ripples_index2 = find(ripples_all(2).session_count == sessions_to_process(nsession)& ripples_all(2).SWS_index == 1);
        ripple_peaktimes2 = min(ripples_all(2).SWR_peaktimes{sessions_to_process(nsession)}(ripples_all(2).probe_hemisphere{sessions_to_process(nsession)} == 2,ripples_all(2).SWS_index(ripples_all(2).session_count == sessions_to_process(nsession))==1))';
        ripple_times2 = [ripples_all(2).onset(ripples_index2) ripples_all(2).offset(ripples_index2)];


        % log event_info
        event_info(nprobe).UP_index = [event_info(nprobe).UP_index UP_index]; 
        event_info(nprobe).DOWN_index = [event_info(nprobe).DOWN_index DOWN_index]; 
        event_info(nprobe).UP_duration = [event_info(nprobe).UP_duration UP_ints(:,2)-UP_ints(:,1)]; 
        event_info(nprobe).previous_DOWN_duration = [event_info(nprobe).previous_DOWN_duration slow_waves_all(nprobe).DOWN_ints(DOWN_index-1,2)-slow_waves_all(nprobe).DOWN_ints(DOWN_index-1,1)]; 
        event_info(nprobe).next_DOWN_duration = [event_info(nprobe).next_DOWN_duration slow_waves_all(nprobe).DOWN_ints(DOWN_index,2)-slow_waves_all(nprobe).DOWN_ints(DOWN_index,1)];


        %%% Get spike counts this session
        NREMInts = behavioural_state_merged_all.SWS{nsession};

        % Define variables
        tvec = [time_bin/2 slow_waves_all(1).V1_MUA_spiketimes{nsession}(end)];
        tvec_interp1 = tvec(1):time_bin:tvec(end);
        tvec_edges = [tvec_interp1(1)-1/(1/mean(diff(tvec_interp1))*2) tvec_interp1+1/(1/mean(diff(tvec_interp1))*2)];
        w = gausswin(0.03*1/mean(diff(tvec_interp1))); % Smoothed with σ = 30 ms
        w = w / sum(w);

        
        DU_slope_HPC = [];
        DU_slope_V1 = [];

        [status,interval,index] = InIntervals(tvec_interp1,NREMInts);

        for mprobe = 1:length(slow_waves_all)

            binnedArrayUPHPC{mprobe} = [];
            binnedArrayDOWNHPC{mprobe} = [];
            binnedArrayRipplesHPC{mprobe} = [];

            binnedArrayUPV1{mprobe} = [];
            binnedArrayDOWNV1{mprobe} = [];
            binnedArrayRipplesV1{mprobe} = [];


            % Overall sleep spikes for zscoring
            spike_times_sleep = slow_waves_all(mprobe).V1_MUA_spiketimes{nsession};
            % spike_speed =  interp1(tvec,Behaviour.mobility,spike_times,'nearest');
            % spike_times_sleep = spike_times(spike_speed < 1);

            % imagesc(filtfilt(w,1,temp));
            V1_spike_counts{mprobe} = filtfilt(w,1,histcounts(spike_times_sleep,tvec_edges));

            spike_times_sleep = slow_waves_all(mprobe).HPC_MUA_spiketimes{nsession};
            % spike_speed =  interp1(tvec,Behaviour.mobility,spike_times,'nearest');
            % spike_times_sleep = spike_times(spike_speed < 1);
            HPC_spike_counts{mprobe} = filtfilt(w,1,histcounts(spike_times_sleep,tvec_edges));

            %%%% Grab perievent HPC MUA activity
            [psth, bins, ~, ~, ~, temp] = psthAndBA(slow_waves_all(mprobe).HPC_MUA_spiketimes{nsession}, slow_waves_all(nprobe).UP_ints(UP_index,1), time_wondows, time_bin);
            temp = (temp-mean(V1_spike_counts{mprobe}(status)))./std(V1_spike_counts{mprobe}(status));% zscore relative to spike count during sleep
            binnedArrayUPHPC{mprobe} = filtfilt(w,1,temp);

            % [psth, bins, ~, ~, ~, temp] = psthAndBA(slow_waves_all(mprobe).V1_MUA_spiketimes{nsession}, ripple_times(nprobe).DOWN_ints(DOWN_index,1), time_wondows, time_bin);
            % temp = (temp-mean(V1_spike_counts{mprobe}(status)))./std(V1_spike_counts{mprobe}(status));% zscore relative to spike count during sleep
            % binnedArrayDOWNHPC{mprobe} = filtfilt(w,1,temp);            

            %%%% Grab perievent V1 MUA activity
            [psth, bins, ~, ~, ~, temp] = psthAndBA(slow_waves_all(mprobe).V1_MUA_spiketimes{nsession}, slow_waves_all(nprobe).UP_ints(UP_index,1), time_wondows, time_bin);
            temp = (temp-mean(V1_spike_counts{mprobe}(status)))./std(V1_spike_counts{mprobe}(status));% zscore relative to spike count during sleep
            binnedArrayUPV1{mprobe} = filtfilt(w,1,temp);

            DU_slope_HPC(mprobe,:) = mean(diff(binnedArrayUPHPC{mprobe}(:,bins >= -0.02 & bins <= 0.04)'));
            DU_slope_V1(mprobe,:) = mean(diff(binnedArrayUPV1{mprobe}(:,bins >= -0.02 & bins <= 0.04)'));
        end

        % log D-U transition MUA spike slope
        event_info(nprobe).DU_slope_HPC = [event_info(nprobe).DU_slope_HPC DU_slope_HPC];
        event_info(nprobe).DU_slope_V1 = [event_info(nprobe).DU_slope_V1 DU_slope_V1];
        
        % UP times
        [index,event_index,UP_time_index] = InIntervals(tvec_interp1, slow_waves_all(nprobe).UP_ints(UP_index,:));



        [~,event_index,normalized_duration,temp] = ...
            calculate_relative_event_probability(slow_waves_all(nprobe).DOWN_ints(DOWN_index,:),ripple_times,num_bins,0);

        for nevent = 1:length(UP_index)
            UP_index(nevent)

            slow_waves_all(nprobe).UP_ints(UP_index(nevent),:)
        end

        [~,event_index,normalized_duration,temp] = ...
            calculate_relative_event_probability(slow_waves_all(nprobe).DOWN_ints(DOWN_index,:),ripple_times,num_bins,0);

        [~,event_index,normalized_duration,temp] = ...
            calculate_relative_event_probability(slow_waves_all(nprobe).UP_ints(UP_index,:),ripple_times,num_bins,0);

        [psth, bins, ~, ~, ~, temp] = psthAndBA(slow_waves_all(mprobe).V1_MUA_spiketimes{nsession}, slow_waves_all(nprobe).UP_ints(UP_index,1), time_wondows, timebin_size);
        temp = (temp-mean(V1_spike_counts{mprobe}(status)))./std(V1_spike_counts{mprobe}(status));% zscore relative to spike count during sleep
        temp = filtfilt(w,1,temp);
        binnedArrayUPV1{mprobe} = [binnedArrayUPV1{mprobe}; filtfilt(w,1,temp)];


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