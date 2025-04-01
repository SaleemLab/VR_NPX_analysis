function PSTH = calculate_UP_DOWN_relative_PSTH(slow_waves_all,ripples_all,behavioural_state_merged_all,sessions_to_process,varargin)

p = inputParser;
addParameter(p,'option','MUA',@ischar);

addParameter(p,'time_option','absolute',@ischar);
addParameter(p,'time_windows',[-1 1],@isnumeric);
addParameter(p,'time_bin',0.01,@isnumeric);
addParameter(p,'num_bins',20,@isnumeric);
addParameter(p,'duration_threshold',2,@isnumeric);
addParameter(p,'shuffle_option','no',@ischar);

parse(p,varargin{:})
option = p.Results.option;
time_option = p.Results.time_option;
time_windows = p.Results.time_windows;
time_bin = p.Results.time_bin;
num_bins = p.Results.num_bins;
duration_threshold = p.Results.duration_threshold;
shuffle_option= p.Results.shuffle_option;

% timebin_edge = time_windows(1):time_bin:time_windows(end);
% bins_centre = timebin_edge(1)+time_bin/2:time_bin:timebin_edge(end)-time_bin/2;


for nprobe = 1:length(slow_waves_all)


    for mprobe = 1:length(slow_waves_all)

        binnedArrayUPHPC{mprobe} = [];
        binnedArrayDOWNHPC{mprobe} = [];
        binnedArrayRipplesHPC{mprobe} = [];
        binnedArraySpindlesHPC{mprobe} = [];

        binnedArrayUPV1{mprobe} = [];
        binnedArrayDOWNV1{mprobe} = [];
        binnedArrayRipplesV1{mprobe} = [];
        binnedArraySpindlesV1{mprobe} = [];
    end

    %%%%%%%%%%%%%%% L ripples
    UP_index_all = [];
    DOWN_index_all = [];
    ripples_index_all = [];
    spindles_index_all = [];

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
        ripple_peaktimes = min(ripples_all(nprobe).SWR_peaktimes{sessions_to_process(nsession)}(ripples_all(nprobe).probe_hemisphere{sessions_to_process(nsession)} == nprobe,ripples_all(nprobe).SWS_index(ripples_all(nprobe).session_count == sessions_to_process(nsession))==1))';
        % if contains(time_option,'peaktimes')
        ripples_index_all = [ripples_index_all; ripples_index];

        % spindles_index = find(spindles_all(1).session_count == sessions_to_process(nsession)& spindles_all(1).SWS_index == 1);
        % spindle_onset = spindles_all(1).onset(spindles_index);
        % 
        % spindles_index_all = [spindles_index_all; ripples_index];

        if contains(shuffle_option,'baseline')
            % s = RandStream('mrg32k3a','Seed',1); % Set random seed for resampling
            % time_jitter = 2 + (2.5 - 2) * rand(s,1, length(UP_index));
            % time_jitter = [time_jitter' time_jitter'];
            time_jitter = [3*ones(1,length(UP_index))' 3*ones(1,length(UP_index))'];
            UP_ints = slow_waves_all(nprobe).UP_ints(UP_index,:)-time_jitter;

            % s = RandStream('mrg32k3a','Seed',2); % Set random seed for resampling
            % time_jitter = 2 + (2.5 - 2) * rand(s,1, length(DOWN_index));
            % time_jitter = [time_jitter' time_jitter'];
            time_jitter = [3*ones(1,length(UP_index))' 3*ones(1,length(UP_index))'];
            DOWN_ints = slow_waves_all(nprobe).DOWN_ints(DOWN_index,:)-time_jitter;

            time_jitter = [3*ones(1,length(ripples_index))'];
            ripple_peaktimes = ripple_peaktimes-time_jitter;

            % time_jitter = [3*ones(1,length(spindles_index))'];
            % spindle_onset = spindle_onset-time_jitter;
        else
            UP_ints = slow_waves_all(nprobe).UP_ints(UP_index,:);
            DOWN_ints = slow_waves_all(nprobe).DOWN_ints(DOWN_index,:);
            ripple_peaktimes= ripple_peaktimes;
            % spindle_onset = spindle_onset;
        end
        
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

            tic
            if contains(option,'MUA')
                [~,~,~,binnedArray] = calculate_relative_event_probability(slow_waves_all(nprobe).UP_ints(UP_index,:), slow_waves_all(mprobe).V1_MUA_spiketimes{nsession},num_bins,[]);
                temp = reshape(binnedArray,1,[]);
                temp = (binnedArray-mean(temp))./std(temp);% zscore relative to spike count during sleep
                binnedArrayUPV1{mprobe} = [binnedArrayUPV1{mprobe}; temp];

                [~,~,~,binnedArray] = calculate_relative_event_probability(slow_waves_all(nprobe).DOWN_ints(DOWN_index,:), slow_waves_all(mprobe).V1_MUA_spiketimes{nsession},num_bins,[]);
                temp = reshape(binnedArray,1,[]);
                temp = (binnedArray-mean(temp))./std(temp);% zscore relative to spike count during sleep
                binnedArrayDOWNV1{mprobe} = [binnedArrayDOWNV1{mprobe}; temp];

                [~,~,~,binnedArray] = calculate_relative_event_probability(slow_waves_all(nprobe).UP_ints(UP_index,:), slow_waves_all(mprobe).HPC_MUA_spiketimes{nsession},num_bins,[]);
                temp = reshape(binnedArray,1,[]);
                temp = (binnedArray-mean(temp))./std(temp);% zscore relative to spike count during sleep
                binnedArrayUPHPC{mprobe} = [binnedArrayUPHPC{mprobe}; temp];

                [~,~,~,binnedArray] = calculate_relative_event_probability(slow_waves_all(nprobe).DOWN_ints(DOWN_index,:), slow_waves_all(mprobe).HPC_MUA_spiketimes{nsession},num_bins,[]);
                temp = reshape(binnedArray,1,[]);
                temp = (binnedArray-mean(temp))./std(temp);% zscore relative to spike count during sleep
                binnedArrayDOWNHPC{mprobe} = [binnedArrayDOWNHPC{mprobe}; temp];

            elseif contains(option,'SUA')

            end
            toc

        end
    end

    % PSTH(nprobe).timebins = bins;

    PSTH(nprobe).V1_UP = binnedArrayUPV1{1};
    PSTH(nprobe).V1_DOWN = binnedArrayDOWNV1{1};

    PSTH(nprobe).R_V1_UP = binnedArrayUPV1{2};
    PSTH(nprobe).R_V1_DOWN = binnedArrayDOWNV1{2};

    PSTH(nprobe).L_HPC_UP = binnedArrayUPHPC{1};
    PSTH(nprobe).L_HPC_DOWN = binnedArrayDOWNHPC{1};

    PSTH(nprobe).R_HPC_UP = binnedArrayUPHPC{2};
    PSTH(nprobe).R_HPC_DOWN = binnedArrayDOWNHPC{2};

end