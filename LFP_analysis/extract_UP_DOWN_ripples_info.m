function event_info = extract_UP_DOWN_ripples_info(slow_waves_all,ripples_all,behavioural_state_merged_all,sessions_to_process,varargin)

p = inputParser;
addParameter(p,'option','UD',@ischar);

% addParameter(p,'time_option','peaktimes',@ischar);
addParameter(p,'time_wondows',[-1 1],@isnumeric);
addParameter(p,'time_bin',0.01,@isnumeric);
addParameter(p,'num_bins',20,@isnumeric);
addParameter(p,'duration_threshold',2,@isnumeric);
addParameter(p,'nprobe',1,@isnumeric);

parse(p,varargin{:})
option = p.Results.option;
% time_option = p.Results.time_option;
time_wondows = p.Results.time_wondows;
time_bin = p.Results.time_bin;
num_bins = p.Results.num_bins;
duration_threshold = p.Results.duration_threshold;
nprobe =  p.Results.nprobe;

% for nprobe = 1:length(slow_waves_all)

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

    % event_info(nprobe).last_L_ripple_index = [];
    % event_info(nprobe).early_L_ripple_index = [];
    % event_info(nprobe).first_L_ripple_index = [];
    % event_info(nprobe).late_L_ripple_index = [];
    event_info(nprobe).L_ripple_normalised_UP_duration = [];
    event_info(nprobe).L_ripple_normalised_DOWN_duration = [];

    % event_info(nprobe).last_R_ripple_index = [];
    % event_info(nprobe).early_R_ripple_index = [];
    % event_info(nprobe).first_R_ripple_index = [];
    % event_info(nprobe).late_R_ripple_index = [];
    event_info(nprobe).R_ripple_normalised_UP_duration = [];
    event_info(nprobe).R_ripple_normalised_DOWN_duration = [];

    event_info(nprobe).L_ripple_zscore_UP=[];
    event_info(nprobe).L_ripple_zscore_DOWN=[];
    event_info(nprobe).R_ripple_zscore_UP=[];
    event_info(nprobe).R_ripple_zscore_DOWN=[];

    event_info(nprobe).L_ripple_cumulative_duration_UP=[];
    event_info(nprobe).L_ripple_cumulative_duration_DOWN=[];
    event_info(nprobe).R_ripple_cumulative_duration_UP=[];
    event_info(nprobe).R_ripple_cumulative_duration_DOWN=[];

    event_info(nprobe).L_ripple_V1_MUA_peak_UP = [];
    event_info(nprobe).L_ripple_HPC_MUA_peak_UP = [];

    event_info(nprobe).R_ripple_V1_MUA_peak_UP = [];
    event_info(nprobe).R_ripple_HPC_MUA_peak_UP = [];


    event_info(nprobe).L_V1_cumulative_activity_UP= [];
    event_info(nprobe).R_V1_cumulative_activity_UP= [];
    event_info(nprobe).L_HPC_cumulative_activity_UP= [];
    event_info(nprobe).R_HPC_cumulative_activity_UP= [];

    event_info(nprobe).L_ripple_V1_MUA_peak_DOWN = [];
    event_info(nprobe).R_ripple_V1_MUA_peak_DOWN = [];
    event_info(nprobe).L_ripple_HPC_MUA_peak_DOWN = [];
    event_info(nprobe).R_ripple_HPC_MUA_peak_DOWN =[];

    event_info(nprobe).L_V1_cumulative_activity_DOWN = [];
    event_info(nprobe).R_V1_cumulative_activity_DOWN = [];
    event_info(nprobe).L_HPC_cumulative_activity_DOWN = [];
    event_info(nprobe).R_HPC_cumulative_activity_DOWN = [];

    event_info(nprobe).L_ripple_V1_MUA_cumulative_UP = [];
    event_info(nprobe).R_ripple_V1_MUA_cumulative_UP = [];
    event_info(nprobe).L_ripple_HPC_MUA_cumulative_UP = [];
    event_info(nprobe).R_ripple_HPC_MUA_cumulative_UP = [];

    event_info(nprobe).L_ripple_V1_MUA_cumulative_DOWN = [];
    event_info(nprobe).R_ripple_V1_MUA_cumulative_DOWN = [];
    event_info(nprobe).L_ripple_HPC_MUA_cumulative_DOWN = [];
    event_info(nprobe).R_ripple_HPC_MUA_cumulative_DOWN = [];


    for nsession = 1:length(sessions_to_process)
        tic
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

        % if contains(option,'UD')
        [C,ia,ib] = intersect(UP_ints(:,2),DOWN_ints(:,1));
        % elseif contains(option,'DU')
        [C,ia,ib] = intersect(UP_ints(:,1),DOWN_ints(:,2));
        % end

        % [C,ia,ib] = intersect(UP_ints(:,2)+0.01,DOWN_ints(:,1));
        UP_index = UP_index(ia);
        DOWN_index = DOWN_index(ib);
        % DOWN_index_all = [DOWN_index_all; DOWN_index];

        UP_ints = slow_waves_all(nprobe).UP_ints(UP_index,:);
        DOWN_ints = slow_waves_all(nprobe).DOWN_ints(DOWN_index,:);
        % UP_index_all = [UP_index_all; UP_index];

        ripple_index1 = find(ripples_all(1).session_count == sessions_to_process(nsession)& ripples_all(1).SWS_index == 1);
        ripple_peaktimes1 = min(ripples_all(1).SWR_peaktimes{sessions_to_process(nsession)}(ripples_all(1).probe_hemisphere{sessions_to_process(nsession)} == 1,ripples_all(1).SWS_index(ripples_all(1).session_count == sessions_to_process(nsession))==1))';
        ripple_times1 = [ripples_all(1).onset(ripple_index1) ripples_all(1).offset(ripple_index1)];

        ripple_index2 = find(ripples_all(2).session_count == sessions_to_process(nsession)& ripples_all(2).SWS_index == 1);
        ripple_peaktimes2 = min(ripples_all(2).SWR_peaktimes{sessions_to_process(nsession)}(ripples_all(2).probe_hemisphere{sessions_to_process(nsession)} == 2,ripples_all(2).SWS_index(ripples_all(2).session_count == sessions_to_process(nsession))==1))';
        ripple_times2 = [ripples_all(2).onset(ripple_index2) ripples_all(2).offset(ripple_index2)];


        % log event_info
        event_info(nprobe).UP_index = [event_info(nprobe).UP_index; UP_index]; % total UP index
        event_info(nprobe).DOWN_index = [event_info(nprobe).DOWN_index; DOWN_index]; % total DOWN index
        event_info(nprobe).UP_duration = [event_info(nprobe).UP_duration; UP_ints(:,2)-UP_ints(:,1)];
        event_info(nprobe).previous_DOWN_duration = [event_info(nprobe).previous_DOWN_duration; slow_waves_all(nprobe).DOWN_ints(DOWN_index-1,2)-slow_waves_all(nprobe).DOWN_ints(DOWN_index-1,1)];
        event_info(nprobe).next_DOWN_duration = [event_info(nprobe).next_DOWN_duration; slow_waves_all(nprobe).DOWN_ints(DOWN_index,2)-slow_waves_all(nprobe).DOWN_ints(DOWN_index,1)];


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

        [NREM_status,interval,~] = InIntervals(tvec_interp1,NREMInts);

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
            V1_spike_counts{mprobe}(V1_spike_counts{mprobe}>prctile(V1_spike_counts{mprobe}(NREM_status),99.5)) = prctile(V1_spike_counts{mprobe}(NREM_status),99.5); % Cap at 99.5% to minimise the influene of outlier

            spike_times_sleep = slow_waves_all(mprobe).HPC_MUA_spiketimes{nsession};
            % spike_speed =  interp1(tvec,Behaviour.mobility,spike_times,'nearest');
            % spike_times_sleep = spike_times(spike_speed < 1);
            HPC_spike_counts{mprobe} = filtfilt(w,1,histcounts(spike_times_sleep,tvec_edges));
            HPC_spike_counts{mprobe}(HPC_spike_counts{mprobe}>prctile(HPC_spike_counts{mprobe}(NREM_status),99.5)) = prctile(HPC_spike_counts{mprobe}(NREM_status),99.5); % Cap at 99.5% to minimise the influene of outlier

            %%%% Grab perievent HPC MUA activity
            [psth, bins, ~, ~, ~, temp] = psthAndBA(slow_waves_all(mprobe).HPC_MUA_spiketimes{nsession}, slow_waves_all(nprobe).UP_ints(UP_index,1), time_wondows, time_bin);
            temp = (temp-mean(V1_spike_counts{mprobe}(NREM_status)))./std(V1_spike_counts{mprobe}(NREM_status));% zscore relative to spike count during sleep
            binnedArrayUPHPC{mprobe} = filtfilt(w,1,temp);


            %%%% Grab perievent V1 MUA activity
            [psth, bins, ~, ~, ~, temp] = psthAndBA(slow_waves_all(mprobe).V1_MUA_spiketimes{nsession}, slow_waves_all(nprobe).UP_ints(UP_index,1), time_wondows, time_bin);
            temp = (temp-mean(V1_spike_counts{mprobe}(NREM_status)))./std(V1_spike_counts{mprobe}(NREM_status));% zscore relative to spike count during sleep
            binnedArrayUPV1{mprobe} = filtfilt(w,1,temp);

            DU_slope_HPC(mprobe,:) = mean(binnedArrayUPV1{mprobe}(:,bins >0 & bins<=0.02),2)- mean(binnedArrayUPV1{mprobe}(:,bins >=-0.02 & bins<0),2);
            DU_slope_V1(mprobe,:) = mean(binnedArrayUPV1{mprobe}(:,bins >0 & bins<=0.02),2)- mean(binnedArrayUPV1{mprobe}(:,bins >=-0.02 & bins<0),2);
        end



        % grab D-U transition MUA spike slope
        event_info(nprobe).DU_slope_HPC = [event_info(nprobe).DU_slope_HPC DU_slope_HPC];
        event_info(nprobe).DU_slope_V1 = [event_info(nprobe).DU_slope_V1 DU_slope_V1];

        % UP spike count index
        [UP_status,UP_event_index,UP_time_index] = InIntervals(tvec_interp1, slow_waves_all(nprobe).UP_ints(UP_index,:));
        %
        [DOWN_status,DOWN_event_index,DOWN_time_index] = InIntervals(tvec_interp1, slow_waves_all(nprobe).DOWN_ints(DOWN_index,:));

        %%%%% Ripple distribution during UP
        % (event_index -> 1 is ripple index, 2 is UP index, 3 is normalized duration)
        [~,event_index1,~,temp] = ...
            calculate_relative_event_probability(slow_waves_all(nprobe).UP_ints(UP_index,:),ripple_peaktimes1,num_bins,0);

        [~,event_index2,~,temp] = ...
            calculate_relative_event_probability(slow_waves_all(nprobe).UP_ints(UP_index,:),ripple_peaktimes2,num_bins,0);

        if ~isempty(event_index1)
            event_info(nprobe).L_ripple_normalised_UP_duration = [event_info(nprobe).L_ripple_normalised_UP_duration;   [ripple_index1(event_index1(:,1)) UP_index(event_index1(:,2)) event_index1(:,3)]];
        else
            event_index1 = [0 0 0];
        end

        if ~isempty(event_index2)
            event_info(nprobe).R_ripple_normalised_UP_duration = [event_info(nprobe).R_ripple_normalised_UP_duration;   [ripple_index2(event_index2(:,1)) UP_index(event_index2(:,2)) event_index2(:,3)]];
        else
            event_index2 = [0 0 0];
        end
        for nevent = 1:length(UP_index)
            % This event ripple
            this_event_index = find(event_index1(:,2) == nevent);

            if ~isempty(this_event_index)
                event_info(nprobe).L_ripple_zscore_UP = [event_info(nprobe).L_ripple_zscore_UP; ripples_all(1).peak_zscore(ripple_index1(event_index1(this_event_index,1)))];

                event_info(nprobe).L_ripple_cumulative_duration_UP = [event_info(nprobe).L_ripple_cumulative_duration_UP...
                    sum(ripples_all(1).offset(ripple_index1(event_index1(this_event_index,1))) - ripples_all(1).onset(ripple_index1(event_index1(this_event_index,1))))];


                for i = 1:length(this_event_index)

                    tidx = tvec_interp1 >= ripples_all(1).onset(ripple_index1(event_index1(this_event_index(i),1))) & tvec_interp1 <= ripples_all(1).offset(ripple_index1(event_index1(this_event_index(i),1)));

                    event_info(nprobe).L_ripple_V1_MUA_peak_UP = [event_info(nprobe).L_ripple_V1_MUA_peak_UP; ...
                        max(V1_spike_counts{1}(tidx)-mean(V1_spike_counts{1}(NREM_status)))./std(V1_spike_counts{1}(NREM_status))...
                        max(V1_spike_counts{2}(tidx)-mean(V1_spike_counts{2}(NREM_status)))./std(V1_spike_counts{2}(NREM_status))];

                    event_info(nprobe).L_ripple_HPC_MUA_peak_UP = [event_info(nprobe).L_ripple_HPC_MUA_peak_UP; ...
                        max(HPC_spike_counts{1}(tidx)-mean(HPC_spike_counts{1}(NREM_status)))./std(HPC_spike_counts{1}(NREM_status))...
                        max(HPC_spike_counts{2}(tidx)-mean(HPC_spike_counts{2}(NREM_status)))./std(HPC_spike_counts{2}(NREM_status))];

                    event_info(nprobe).L_ripple_V1_MUA_cumulative_UP = [event_info(nprobe).L_ripple_V1_MUA_cumulative_UP; ...
                        sum((V1_spike_counts{1}(tidx)-min(V1_spike_counts{1}(NREM_status)))./(max(V1_spike_counts{1}(NREM_status))-min(V1_spike_counts{1}(NREM_status))))...
                        sum((V1_spike_counts{2}(tidx)-min(V1_spike_counts{2}(NREM_status)))./(max(V1_spike_counts{2}(NREM_status))-min(V1_spike_counts{2}(NREM_status))))];

                    event_info(nprobe).L_ripple_HPC_MUA_cumulative_UP = [event_info(nprobe).L_ripple_HPC_MUA_cumulative_UP; ...
                        sum((HPC_spike_counts{1}(tidx)-min(HPC_spike_counts{1}(NREM_status)))./(max(HPC_spike_counts{1}(NREM_status))-min(HPC_spike_counts{1}(NREM_status))))...
                        sum((HPC_spike_counts{2}(tidx)-min(HPC_spike_counts{2}(NREM_status)))./(max(HPC_spike_counts{2}(NREM_status))-min(HPC_spike_counts{2}(NREM_status))))];


                end

            end

            this_event_index = find(event_index2(:,2) == nevent);

            if ~isempty(this_event_index)
                event_info(nprobe).R_ripple_zscore_UP = [event_info(nprobe).R_ripple_zscore_UP; ripples_all(2).peak_zscore(ripple_index2(event_index2(this_event_index,1)))];

                event_info(nprobe).R_ripple_cumulative_duration_UP = [event_info(nprobe).R_ripple_cumulative_duration_UP...
                    sum(ripples_all(2).offset(ripple_index2(event_index2(this_event_index,1))) - ripples_all(2).onset(ripple_index2(event_index2(this_event_index,1))))];

                for i = 1:length(this_event_index)
                    tidx = tvec_interp1 >= ripples_all(2).onset(ripple_index2(event_index2(this_event_index(i),1))) & tvec_interp1 <= ripples_all(2).offset(ripple_index2(event_index2(this_event_index(i),1)));

                    event_info(nprobe).R_ripple_V1_MUA_peak_UP = [event_info(nprobe).R_ripple_V1_MUA_peak_UP; ...
                        max(V1_spike_counts{1}(tidx)-mean(V1_spike_counts{1}(NREM_status)))./std(V1_spike_counts{1}(NREM_status))...
                        max(V1_spike_counts{2}(tidx)-mean(V1_spike_counts{2}(NREM_status)))./std(V1_spike_counts{2}(NREM_status))];

                    event_info(nprobe).R_ripple_HPC_MUA_peak_UP = [event_info(nprobe).R_ripple_HPC_MUA_peak_UP; ...
                        max(HPC_spike_counts{1}(tidx)-mean(HPC_spike_counts{1}(NREM_status)))./std(HPC_spike_counts{1}(NREM_status))...
                        max(HPC_spike_counts{2}(tidx)-mean(HPC_spike_counts{2}(NREM_status)))./std(HPC_spike_counts{2}(NREM_status))];

                    event_info(nprobe).R_ripple_V1_MUA_cumulative_UP = [event_info(nprobe).R_ripple_V1_MUA_cumulative_UP; ...
                        sum((V1_spike_counts{1}(tidx)-min(V1_spike_counts{1}(NREM_status)))./(max(V1_spike_counts{1}(NREM_status))-min(V1_spike_counts{1}(NREM_status))))...
                        sum((V1_spike_counts{2}(tidx)-min(V1_spike_counts{2}(NREM_status)))./(max(V1_spike_counts{2}(NREM_status))-min(V1_spike_counts{2}(NREM_status))))];

                    event_info(nprobe).R_ripple_HPC_MUA_cumulative_UP = [event_info(nprobe).R_ripple_HPC_MUA_cumulative_UP; ...
                        sum((HPC_spike_counts{1}(tidx)-min(HPC_spike_counts{1}(NREM_status)))./(max(HPC_spike_counts{1}(NREM_status))-min(HPC_spike_counts{1}(NREM_status))))...
                        sum((HPC_spike_counts{2}(tidx)-min(HPC_spike_counts{2}(NREM_status)))./(max(HPC_spike_counts{2}(NREM_status))-min(HPC_spike_counts{2}(NREM_status))))];
                end
            end

            % min-max normalised cumulative spiking

            event_info(nprobe).L_V1_cumulative_activity_UP= [event_info(nprobe).L_V1_cumulative_activity_UP sum((V1_spike_counts{1}(UP_event_index == nevent)-min(V1_spike_counts{1}(NREM_status)))./(max(V1_spike_counts{1}(NREM_status))-min(V1_spike_counts{1}(NREM_status))))];

            event_info(nprobe).R_V1_cumulative_activity_UP= [event_info(nprobe).R_V1_cumulative_activity_UP sum(V1_spike_counts{2}(UP_event_index == nevent)-min(V1_spike_counts{2}(NREM_status)))./(max(V1_spike_counts{1}(NREM_status))-min(V1_spike_counts{2}(NREM_status)))];

            event_info(nprobe).L_HPC_cumulative_activity_UP= [event_info(nprobe).L_HPC_cumulative_activity_UP sum((HPC_spike_counts{1}(UP_event_index == nevent)-min(HPC_spike_counts{1}(NREM_status))))./(max(V1_spike_counts{1}(NREM_status))-min(HPC_spike_counts{1}(NREM_status)))];

            event_info(nprobe).R_HPC_cumulative_activity_UP= [event_info(nprobe).R_HPC_cumulative_activity_UP sum(HPC_spike_counts{2}(UP_event_index == nevent)-min(HPC_spike_counts{2}(NREM_status)))./(max(V1_spike_counts{1}(NREM_status))-min(HPC_spike_counts{2}(NREM_status)))];

        end

        % event_info(nprobe).L_ripple_zscore_UP=[];
        % event_info(nprobe).L_ripple_zscore_DOWN=[];
        % event_info(nprobe).R_ripple_zscore_UP=[];
        % event_info(nprobe).R_ripple_zscore_DOWN=[];
        %
        % event_info(nprobe).L_ripple_cumulative_duration_UP=[];
        % event_info(nprobe).L_ripple_cumulative_duration_DOWN=[];
        % event_info(nprobe).R_ripple_cumulative_duration_UP=[];
        % event_info(nprobe).R_ripple_cumulative_duration_DOWN=[];


        %%%%% Ripple distribution during DOWN
        [DOWN_status,DOWN_event_index,DOWN_time_index] = InIntervals(tvec_interp1, slow_waves_all(nprobe).DOWN_ints(DOWN_index,:));


        % Ripple distribution during DOWN
        % (event_index -> 1 is ripple index, 2 is UP index, 3 is normalized duration)
        [~,event_index1,~,temp] = ...
            calculate_relative_event_probability(slow_waves_all(nprobe).DOWN_ints(DOWN_index,:),ripple_peaktimes1,num_bins,0);

        [~,event_index2,~,temp] = ...
            calculate_relative_event_probability(slow_waves_all(nprobe).DOWN_ints(DOWN_index,:),ripple_peaktimes2,num_bins,0);


        if ~isempty(event_index1)
            event_info(nprobe).L_ripple_normalised_DOWN_duration = [event_info(nprobe).L_ripple_normalised_DOWN_duration;   [ripple_index1(event_index1(:,1)) DOWN_index(event_index1(:,2)) event_index1(:,3)]];
        else
            event_index1 = [0 0 0];
        end

        if ~isempty(event_index2)
            event_info(nprobe).R_ripple_normalised_DOWN_duration = [event_info(nprobe).R_ripple_normalised_DOWN_duration;   [ripple_index2(event_index2(:,1)) DOWN_index(event_index2(:,2)) event_index2(:,3)]];
        else
            event_index2 = [0 0 0];
        end

        for nevent = 1:length(DOWN_index)
            % This event ripple
            this_event_index = find(event_index1(:,2) == nevent);

            if ~isempty(this_event_index)
                event_info(nprobe).L_ripple_zscore_DOWN = [event_info(nprobe).L_ripple_zscore_DOWN; ripples_all(1).peak_zscore(ripple_index1(event_index1(this_event_index,1)))];

                event_info(nprobe).L_ripple_cumulative_duration_DOWN = [event_info(nprobe).L_ripple_cumulative_duration_DOWN...
                    sum(ripples_all(1).offset(ripple_index1(event_index1(this_event_index,1))) - ripples_all(1).onset(ripple_index1(event_index1(this_event_index,1))))];


                for i = 1:length(this_event_index)

                    tidx = tvec_interp1 >= ripples_all(1).onset(ripple_index1(event_index1(this_event_index(i),1))) & tvec_interp1 <= ripples_all(1).offset(ripple_index1(event_index1(this_event_index(i),1)));

                    event_info(nprobe).L_ripple_V1_MUA_peak_DOWN = [event_info(nprobe).L_ripple_V1_MUA_peak_DOWN; ...
                        max(V1_spike_counts{1}(tidx)-mean(V1_spike_counts{1}(NREM_status)))./std(V1_spike_counts{1}(NREM_status))...
                        max(V1_spike_counts{2}(tidx)-mean(V1_spike_counts{2}(NREM_status)))./std(V1_spike_counts{2}(NREM_status))];

                    event_info(nprobe).L_ripple_HPC_MUA_peak_DOWN = [event_info(nprobe).L_ripple_HPC_MUA_peak_DOWN; ...
                        max(HPC_spike_counts{1}(tidx)-mean(HPC_spike_counts{1}(NREM_status)))./std(HPC_spike_counts{1}(NREM_status))...
                        max(HPC_spike_counts{2}(tidx)-mean(HPC_spike_counts{2}(NREM_status)))./std(HPC_spike_counts{2}(NREM_status))];

                    event_info(nprobe).L_ripple_V1_MUA_cumulative_DOWN = [event_info(nprobe).L_ripple_V1_MUA_cumulative_DOWN; ...
                        sum((V1_spike_counts{1}(tidx)-min(V1_spike_counts{1}(NREM_status)))./(max(V1_spike_counts{1}(NREM_status))-min(V1_spike_counts{1}(NREM_status))))...
                        sum((V1_spike_counts{2}(tidx)-min(V1_spike_counts{2}(NREM_status)))./(max(V1_spike_counts{2}(NREM_status))-min(V1_spike_counts{2}(NREM_status))))];

                    event_info(nprobe).L_ripple_HPC_MUA_cumulative_DOWN = [event_info(nprobe).L_ripple_HPC_MUA_cumulative_DOWN; ...
                        sum((HPC_spike_counts{1}(tidx)-min(HPC_spike_counts{1}(NREM_status)))./(max(HPC_spike_counts{1}(NREM_status))-min(HPC_spike_counts{1}(NREM_status))))...
                        sum((HPC_spike_counts{2}(tidx)-min(HPC_spike_counts{2}(NREM_status)))./(max(HPC_spike_counts{2}(NREM_status))-min(HPC_spike_counts{2}(NREM_status))))];

                end

            end

            this_event_index = find(event_index2(:,2) == nevent);

            if ~isempty(this_event_index)
                event_info(nprobe).R_ripple_zscore_DOWN = [event_info(nprobe).R_ripple_zscore_DOWN; ripples_all(2).peak_zscore(ripple_index2(event_index2(this_event_index,1)))];

                event_info(nprobe).R_ripple_cumulative_duration_DOWN = [event_info(nprobe).R_ripple_cumulative_duration_DOWN...
                    sum(ripples_all(2).offset(ripple_index2(event_index2(this_event_index,1))) - ripples_all(2).onset(ripple_index2(event_index2(this_event_index,1))))];

                for i = 1:length(this_event_index)
                    tidx = tvec_interp1 >= ripples_all(2).onset(ripple_index2(event_index2(this_event_index(i),1))) & tvec_interp1 <= ripples_all(2).offset(ripple_index2(event_index2(this_event_index(i),1)));

                    event_info(nprobe).R_ripple_V1_MUA_peak_DOWN = [event_info(nprobe).R_ripple_V1_MUA_peak_DOWN; ...
                        max(V1_spike_counts{1}(tidx)-mean(V1_spike_counts{1}(NREM_status)))./std(V1_spike_counts{1}(NREM_status))...
                        max(V1_spike_counts{2}(tidx)-mean(V1_spike_counts{2}(NREM_status)))./std(V1_spike_counts{2}(NREM_status))];

                    event_info(nprobe).R_ripple_HPC_MUA_peak_DOWN = [event_info(nprobe).R_ripple_HPC_MUA_peak_DOWN; ...
                        max(HPC_spike_counts{1}(tidx)-mean(HPC_spike_counts{1}(NREM_status)))./std(HPC_spike_counts{1}(NREM_status))...
                        max(HPC_spike_counts{2}(tidx)-mean(HPC_spike_counts{2}(NREM_status)))./std(HPC_spike_counts{2}(NREM_status))];

                    event_info(nprobe).R_ripple_V1_MUA_cumulative_DOWN = [event_info(nprobe).R_ripple_V1_MUA_cumulative_DOWN; ...
                        sum((V1_spike_counts{1}(tidx)-min(V1_spike_counts{1}(NREM_status)))./(max(V1_spike_counts{1}(NREM_status))-min(V1_spike_counts{1}(NREM_status))))...
                        sum((V1_spike_counts{2}(tidx)-min(V1_spike_counts{2}(NREM_status)))./(max(V1_spike_counts{2}(NREM_status))-min(V1_spike_counts{2}(NREM_status))))];

                    event_info(nprobe).R_ripple_HPC_MUA_cumulative_DOWN = [event_info(nprobe).R_ripple_HPC_MUA_cumulative_DOWN; ...
                        sum((HPC_spike_counts{1}(tidx)-min(HPC_spike_counts{1}(NREM_status)))./(max(HPC_spike_counts{1}(NREM_status))-min(HPC_spike_counts{1}(NREM_status))))...
                        sum((HPC_spike_counts{2}(tidx)-min(HPC_spike_counts{2}(NREM_status)))./(max(HPC_spike_counts{2}(NREM_status))-min(HPC_spike_counts{2}(NREM_status))))];
                end
            end

            % min-max normalised cumulative spiking

             event_info(nprobe).L_V1_cumulative_activity_DOWN= [event_info(nprobe).L_V1_cumulative_activity_DOWN sum((V1_spike_counts{1}(DOWN_event_index == nevent)-min(V1_spike_counts{1}(NREM_status)))./(max(V1_spike_counts{1}(NREM_status))-min(V1_spike_counts{1}(NREM_status))))];

            event_info(nprobe).R_V1_cumulative_activity_DOWN= [event_info(nprobe).R_V1_cumulative_activity_DOWN sum(V1_spike_counts{2}(DOWN_event_index == nevent)-min(V1_spike_counts{2}(NREM_status)))./(max(V1_spike_counts{1}(NREM_status))-min(V1_spike_counts{2}(NREM_status)))];

            event_info(nprobe).L_HPC_cumulative_activity_DOWN= [event_info(nprobe).L_HPC_cumulative_activity_DOWN sum((HPC_spike_counts{1}(DOWN_event_index == nevent)-min(HPC_spike_counts{1}(NREM_status))))./(max(V1_spike_counts{1}(NREM_status))-min(HPC_spike_counts{1}(NREM_status)))];

            event_info(nprobe).R_HPC_cumulative_activity_DOWN= [event_info(nprobe).R_HPC_cumulative_activity_DOWN sum(HPC_spike_counts{2}(DOWN_event_index == nevent)-min(HPC_spike_counts{2}(NREM_status)))./(max(V1_spike_counts{1}(NREM_status))-min(HPC_spike_counts{2}(NREM_status)))];

        end
        toc
    end

% end

