function PSTH = calculate_UP_DOWN_relative_PSTH(slow_waves_all,ripples_all,behavioural_state_merged_all,sessions_to_process,varargin)

p = inputParser;
addParameter(p,'option','MUA',@ischar);

addParameter(p,'time_option','absolute',@ischar);
addParameter(p,'time_wondows',[-1 1],@isnumeric);
addParameter(p,'time_bin',0.01,@isnumeric);
addParameter(p,'num_bins',20,@isnumeric);
addParameter(p,'duration_threshold',2,@isnumeric);

parse(p,varargin{:})
option = p.Results.option;
time_option = p.Results.time_option;
time_wondows = p.Results.time_wondows;
time_bin = p.Results.time_bin;
num_bins = p.Results.num_bins;
duration_threshold = p.Results.duration_threshold;
shuffle_options = 0;

for nprobe = 1:length(slow_waves_all)


    for mprobe = 1:length(slow_waves_all)

        binnedArrayUPHPC{mprobe} = [];
        binnedArrayDOWNHPC{mprobe} = [];
        binnedArrayRipplesHPC{mprobe} = [];

        binnedArrayUPV1{mprobe} = [];
        binnedArrayDOWNV1{mprobe} = [];
        binnedArrayRipplesV1{mprobe} = [];
    end

    %%%%%%%%%%%%%%% L ripples
    UP_index_all = [];
    DOWN_index_all = [];
    ripples_index_all = [];
    ripple_index_UP_all = [];
    ripple_index_DOWN_all = [];

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

        ripples_index = find(ripples_all(nprobe).session_count == sessions_to_process(nsession)& ripples_all(nprobe).SWS_index == 1);
        ripple_peaktimes = min(ripples_all(nprobe).SWR_peaktimes{sessions_to_process(nsession)}(ripples_all(nprobe).probe_hemisphere{sessions_to_process(nsession)} == 2,ripples_all(nprobe).SWS_index(ripples_all(nprobe).session_count == sessions_to_process(nsession))==1))';
        % if contains(time_option,'peaktimes')
        ripple_times= ripple_peaktimes;

        ripples_index_all = [ripples_index_all; ripples_index];
        % else
        %     ripple_times = [ripples_all(1).onset(ripples_index) ripples_all(1).offset(ripples_index)];
        % end

        %%% Get spike counts this session
        NREMInts = behavioural_state_merged_all.SWS{nsession};

        % Define variables
        timebin_size = 0.01;
        tvec = [timebin_size/2 slow_waves_all(1).V1_MUA_spiketimes{nsession}(end)];
        tvec_interp1 = tvec(1):timebin_size:tvec(end);
        tvec_edges = [tvec_interp1(1)-1/(1/mean(diff(tvec_interp1))*2) tvec_interp1+1/(1/mean(diff(tvec_interp1))*2)];
        w = gausswin(0.03*1/mean(diff(tvec_interp1))); % Smoothed with σ = 30 ms
        w = w / sum(w);

        [status,interval,index] = InIntervals(tvec_interp1,NREMInts);

        % Mainly used for zscore
        for mprobe = 1:length(slow_waves_all)
            spike_times_sleep = slow_waves_all(mprobe).V1_MUA_spiketimes{nsession};
            % spike_speed =  interp1(tvec,Behaviour.mobility,spike_times,'nearest');
            % spike_times_sleep = spike_times(spike_speed < 1);

            % imagesc(filtfilt(w,1,temp));
            V1_spike_counts{mprobe} = filtfilt(w,1,histcounts(spike_times_sleep,tvec_edges));

            spike_times_sleep = slow_waves_all(mprobe).HPC_MUA_spiketimes{nsession};
            % spike_speed =  interp1(tvec,Behaviour.mobility,spike_times,'nearest');
            % spike_times_sleep = spike_times(spike_speed < 1);
            HPC_spike_counts{mprobe} = filtfilt(w,1,histcounts(spike_times_sleep,tvec_edges));


            if contains(option,'MUA')
                [probabilities,event_index,normalized_duration,binnedArray] = calculate_relative_event_probability(slow_waves_all(nprobe).UP_ints(UP_index,:), slow_waves_all(mprobe).V1_MUA_spiketimes{nsession},num_bins,shuffle_options);
                temp = reshape(binnedArray,1,[]);
                temp = (binnedArray-mean(temp))./std(temp);% zscore relative to spike count during sleep
                binnedArrayUPV1{mprobe} = [binnedArrayUPV1{mprobe}; temp];

                [probabilities,event_index,normalized_duration,binnedArray] = calculate_relative_event_probability(slow_waves_all(nprobe).DOWN_ints(DOWN_index,:), slow_waves_all(mprobe).V1_MUA_spiketimes{nsession},num_bins,shuffle_options);
                temp = reshape(binnedArray,1,[]);
                temp = (binnedArray-mean(temp))./std(temp);% zscore relative to spike count during sleep
                binnedArrayDOWNV1{mprobe} = [binnedArrayDOWNV1{mprobe}; temp];

                [probabilities,event_index,normalized_duration,binnedArray] = calculate_relative_event_probability(slow_waves_all(nprobe).UP_ints(UP_index,:), slow_waves_all(mprobe).HPC_MUA_spiketimes{nsession},num_bins,shuffle_options);
                temp = reshape(binnedArray,1,[]);
                temp = (binnedArray-mean(temp))./std(temp);% zscore relative to spike count during sleep
                binnedArrayUPHPC{mprobe} = [binnedArrayUPHPC{mprobe}; temp];

                [probabilities,event_index,normalized_duration,binnedArray] = calculate_relative_event_probability(slow_waves_all(nprobe).DOWN_ints(DOWN_index,:), slow_waves_all(mprobe).HPC_MUA_spiketimes{nsession},num_bins,shuffle_options);
                temp = reshape(binnedArray,1,[]);
                temp = (binnedArray-mean(temp))./std(temp);% zscore relative to spike count during sleep
                binnedArrayDOWNHPC{mprobe} = [binnedArrayDOWNHPC{mprobe}; temp];

            elseif contains(option,'SUA')

            end


        end
    end

    PSTH(nprobe).timebins = bins;

    PSTH(nprobe).V1_UP = binnedArrayUPV1{1};
    PSTH(nprobe).V1_DOWN = binnedArrayDOWNV1{1};

    PSTH(nprobe).R_V1_UP = binnedArrayUPV1{2};
    PSTH(nprobe).R_V1_DOWN = binnedArrayDOWNV1{2};

    PSTH(nprobe).L_HPC_UP = binnedArrayUPHPC{1};
    PSTH(nprobe).L_HPC_DOWN = binnedArrayDOWNHPC{1};

    PSTH(nprobe).R_HPC_UP = binnedArrayUPHPC{2};
    PSTH(nprobe).R_HPC_DOWN = binnedArrayDOWNHPC{2};

end