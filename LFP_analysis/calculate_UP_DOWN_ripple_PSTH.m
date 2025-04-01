function PSTH = calculate_UP_DOWN_ripple_PSTH(slow_waves_all,ripples_all,behavioural_state_merged_all,sessions_to_process,varargin)

p = inputParser;
addParameter(p,'option','MUA',@ischar);

addParameter(p,'time_option','absolute',@ischar);
addParameter(p,'time_windows',[-1 1],@isnumeric);
addParameter(p,'time_bin',0.01,@isnumeric);
addParameter(p,'num_bins',20,@isnumeric);
addParameter(p,'duration_threshold',2,@isnumeric);
addParameter(p,'shuffle_option','time_circular_shift',@ischar);
% addParameter(p,'shuffle_option','time_circular_shift',@ischar);

parse(p,varargin{:})
option = p.Results.option;
time_option = p.Results.time_option;
time_windows = p.Results.time_windows;
time_bin = p.Results.time_bin;
num_bins = p.Results.num_bins;
duration_threshold = p.Results.duration_threshold;
shuffle_option= p.Results.shuffle_option;

timebin_edge = time_windows(1):time_bin:time_windows(end);
bins_centre = timebin_edge(1)+time_bin/2:time_bin:timebin_edge(end)-time_bin/2;


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
        ripple_peaktimes = min(ripples_all(nprobe).SWR_peaktimes{sessions_to_process(nsession)}(ripples_all(nprobe).probe_hemisphere{sessions_to_process(nsession)} == nprobe,ripples_all(nprobe).SWS_index(ripples_all(nprobe).session_count == sessions_to_process(nsession))==1))';
        % if contains(time_option,'peaktimes')
        ripple_times= ripple_peaktimes;

        ripples_index_all = [ripples_index_all; ripples_index];

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

        else
            UP_ints = slow_waves_all(nprobe).UP_ints(UP_index,:);
            DOWN_ints = slow_waves_all(nprobe).DOWN_ints(DOWN_index,:);
            ripple_peaktimes= ripple_peaktimes;
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


            if contains(option,'MUA')

                % V1 MUA UP
                [psth, bins, ~, ~, ~, temp] = psthAndBA(slow_waves_all(mprobe).V1_MUA_spiketimes{nsession}, UP_ints(:,1), time_windows, timebin_size);

                if ~contains(shuffle_option,'baseline')
                    timebin_edges_all = UP_ints(:,1) + bins_centre;  % Absolute times of peri-event window

                    for i = 1:size(UP_ints,1)

                        timebin_edges_all(i,:);
                        % Previous UP (skip if this is the first UP)
                        if i > 1
                            prev_offset = UP_ints(i-1,2);
                            % Find peri-time indices within the previous UP state
                            mask_prev =  timebin_edges_all(i,:) <= prev_offset;
                            temp(i, mask_prev) = NaN;
                        end

                        % Next UP (skip if this is the last UP)
                        if i < size(UP_ints,1)
                            next_onset = UP_ints(i+1,1);
                            % Find peri-time indices within the next UP state
                            mask_next = timebin_edges_all(i,:) >= next_onset;
                            temp(i, mask_next) = NaN;
                        end
                    end
                end

                temp = (temp-mean(V1_spike_counts{mprobe}(status)))./std(V1_spike_counts{mprobe}(status));% zscore relative to spike count during sleep
                binnedArrayUPV1{mprobe} = [binnedArrayUPV1{mprobe}; filtfilt(w,1,temp)];



               % V1 MUA DOWN
                [psth, bins, ~, ~, ~, temp] = psthAndBA(slow_waves_all(mprobe).V1_MUA_spiketimes{nsession}, DOWN_ints(:,1), time_windows, timebin_size);

                if ~contains(shuffle_option,'baseline')
                    timebin_edges_all = DOWN_ints(:,1) + bins_centre;  % Absolute times of peri-event window
                    for i = 1:size(DOWN_ints,1)

                        timebin_edges_all(i,:);
                        % Previous DOWN (skip if this is the first UP)
                        if i > 1
                            prev_offset = DOWN_ints(i-1,2);
                            % Find peri-time indices within the previous UP state
                            mask_prev =  timebin_edges_all(i,:) <= prev_offset;
                            temp(i, mask_prev) = NaN;
                        end

                        % Next DOWN (skip if this is the last UP)
                        if i < size(DOWN_ints,1)
                            next_onset = DOWN_ints(i+1,1);
                            % Find peri-time indices within the next UP state
                            mask_next = timebin_edges_all(i,:) >= next_onset;
                            temp(i, mask_next) = NaN;
                        end
                    end
                end

                temp = (temp-mean(V1_spike_counts{mprobe}(status)))./std(V1_spike_counts{mprobe}(status));% zscore relative to spike count during sleep
                binnedArrayDOWNV1{mprobe} = [binnedArrayDOWNV1{mprobe}; filtfilt(w,1,temp)];


               % V1 MUA ripples
                [psth, bins, ~, ~, ~, temp] = psthAndBA(slow_waves_all(mprobe).V1_MUA_spiketimes{nsession}, ripple_peaktimes, time_windows, timebin_size);

                if ~contains(shuffle_option,'baseline')
                    ripple_times = [ripples_all(nprobe).onset(ripples_index) ripples_all(nprobe).offset(ripples_index)];
                    timebin_edges_all = ripple_peaktimes + bins_centre;  % Absolute times of peri-event window
                    for i = 1:size(ripple_peaktimes,1)
                        % Previous DOWN (skip if this is the first UP)
                        if i > 1
                            prev_offset = ripple_times(i-1,2);
                            % Find peri-time indices within the previous UP state
                            mask_prev =  timebin_edges_all(i,:) <= prev_offset;
                            temp(i, mask_prev) = NaN;
                        end

                        % Next DOWN (skip if this is the last UP)
                        if i < size(ripple_times,1)
                            next_onset = ripple_times(i+1,1);
                            % Find peri-time indices within the next UP state
                            mask_next = timebin_edges_all(i,:) >= next_onset;
                            temp(i, mask_next) = NaN;
                        end
                    end
                end

                temp = (temp-mean(V1_spike_counts{mprobe}(status)))./std(V1_spike_counts{mprobe}(status));% zscore relative to spike count during sleep
                binnedArrayRipplesV1{mprobe} = [binnedArrayRipplesV1{mprobe}; filtfilt(w,1,temp)];


                %%%%% HPC MUA
                % HPC MUA UP
                [psth, bins, ~, ~, ~, temp] = psthAndBA(slow_waves_all(mprobe).HPC_MUA_spiketimes{nsession}, UP_ints(:,1), time_windows, timebin_size);

                if ~contains(shuffle_option,'baseline')
                    timebin_edges_all = UP_ints(:,1) + bins_centre;  % Absolute times of peri-event window

                    for i = 1:size(UP_ints,1)

                        timebin_edges_all(i,:);
                        % Previous UP (skip if this is the first UP)
                        if i > 1
                            prev_offset = UP_ints(i-1,2);
                            % Find peri-time indices within the previous UP state
                            mask_prev =  timebin_edges_all(i,:) <= prev_offset;
                            temp(i, mask_prev) = NaN;
                        end

                        % Next UP (skip if this is the last UP)
                        if i < size(UP_ints,1)
                            next_onset = UP_ints(i+1,1);
                            % Find peri-time indices within the next UP state
                            mask_next = timebin_edges_all(i,:) >= next_onset;
                            temp(i, mask_next) = NaN;
                        end
                    end
                end

                temp = (temp-mean(HPC_spike_counts{mprobe}(status)))./std(HPC_spike_counts{mprobe}(status));% zscore relative to spike count during sleep
                binnedArrayUPHPC{mprobe} = [binnedArrayUPHPC{mprobe}; filtfilt(w,1,temp)];

                


                % HPC MUA DOWN
                [psth, bins, ~, ~, ~, temp] = psthAndBA(slow_waves_all(mprobe).HPC_MUA_spiketimes{nsession}, DOWN_ints(:,1), time_windows, timebin_size);

                if ~contains(shuffle_option,'baseline')
                    timebin_edges_all = DOWN_ints(:,1) + bins_centre;  % Absolute times of peri-event window
                    for i = 1:size(DOWN_ints,1)

                        timebin_edges_all(i,:);
                        % Previous DOWN (skip if this is the first UP)
                        if i > 1
                            prev_offset = DOWN_ints(i-1,2);
                            % Find peri-time indices within the previous UP state
                            mask_prev =  timebin_edges_all(i,:) <= prev_offset;
                            temp(i, mask_prev) = NaN;
                        end

                        % Next DOWN (skip if this is the last UP)
                        if i < size(DOWN_ints,1)
                            next_onset = DOWN_ints(i+1,1);
                            % Find peri-time indices within the next UP state
                            mask_next = timebin_edges_all(i,:) >= next_onset;
                            temp(i, mask_next) = NaN;
                        end
                    end
                end

                temp = (temp-mean(HPC_spike_counts{mprobe}(status)))./std(HPC_spike_counts{mprobe}(status));% zscore relative to spike count during sleep
                binnedArrayDOWNHPC{mprobe} = [binnedArrayDOWNHPC{mprobe}; filtfilt(w,1,temp)];



                % HPC MUA ripple
                [psth, bins, ~, ~, ~, temp] = psthAndBA(slow_waves_all(mprobe).HPC_MUA_spiketimes{nsession}, ripple_peaktimes, time_windows, timebin_size);
                temp = (temp-mean(HPC_spike_counts{mprobe}(status)))./std(HPC_spike_counts{mprobe}(status));% zscore relative to spike count during sleep
                binnedArrayRipplesHPC{mprobe} = [binnedArrayRipplesHPC{mprobe}; filtfilt(w,1,temp)];

                if ~contains(shuffle_option,'baseline')
                    ripple_times = [ripples_all(nprobe).onset(ripples_index) ripples_all(nprobe).offset(ripples_index)];
                    timebin_edges_all = ripple_peaktimes + bins_centre;  % Absolute times of peri-event window
                    for i = 1:size(ripple_peaktimes,1)
                        % Previous DOWN (skip if this is the first UP)
                        if i > 1
                            prev_offset = ripple_times(i-1,2);
                            % Find peri-time indices within the previous UP state
                            mask_prev =  timebin_edges_all(i,:) <= prev_offset;
                            temp(i, mask_prev) = NaN;
                        end

                        % Next DOWN (skip if this is the last UP)
                        if i < size(ripple_times,1)
                            next_onset = ripple_times(i+1,1);
                            % Find peri-time indices within the next UP state
                            mask_next = timebin_edges_all(i,:) >= next_onset;
                            temp(i, mask_next) = NaN;
                        end
                    end
                end
            elseif contains(option,'SUA')

            end


        end
    end

    PSTH(nprobe).timebins = bins;

    PSTH(nprobe).L_V1_UP = binnedArrayUPV1{1};
    PSTH(nprobe).L_V1_DOWN = binnedArrayDOWNV1{1};
    PSTH(nprobe).L_V1_ripples = binnedArrayRipplesV1{1};

    PSTH(nprobe).R_V1_UP = binnedArrayUPV1{2};
    PSTH(nprobe).R_V1_DOWN = binnedArrayDOWNV1{2};
    PSTH(nprobe).R_V1_ripples = binnedArrayRipplesV1{2};

    PSTH(nprobe).L_HPC_UP = binnedArrayUPHPC{1};
    PSTH(nprobe).L_HPC_DOWN = binnedArrayDOWNHPC{1};
    PSTH(nprobe).L_HPC_ripples = binnedArrayRipplesHPC{1};

    PSTH(nprobe).R_HPC_UP = binnedArrayUPHPC{2};
    PSTH(nprobe).R_HPC_DOWN = binnedArrayDOWNHPC{2};
    PSTH(nprobe).R_HPC_ripples = binnedArrayRipplesHPC{2};

    % V1 L
    % bootstrap distribution
    tempUP = [];
    tempDOWN = [];
    tempRipples = [];
    for iBoot = 1:1000
        s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
        event_id = datasample(s,1:size(binnedArrayUPV1{1},1),size(binnedArrayUPV1{1},1));
        tempUP(iBoot,:) = mean(binnedArrayUPV1{1}(event_id,:),'omitnan');

        s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
        event_id = datasample(s,1:size(binnedArrayDOWNV1{1},1),size(binnedArrayDOWNV1{1},1));
        tempDOWN(iBoot,:) = mean(binnedArrayDOWNV1{1}(event_id,:),'omitnan');

        s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
        event_id = datasample(s,1:size(binnedArrayRipplesV1{1},1),size(binnedArrayRipplesV1{1},1));
        tempRipples(iBoot,:) = mean(binnedArrayRipplesV1{1}(event_id,:),'omitnan');
    end

    PSTH(nprobe).L_V1_UP_bootstrap = tempUP;
    PSTH(nprobe).L_V1_DOWN_bootstrap = tempDOWN;
    PSTH(nprobe).L_V1_ripples_bootstrap = tempRipples;

    if contains(shuffle_option,'time_circular_shift')
        % timebin circularly shifted
        tempUP = [];
        tempDOWN = [];
        tempRipples = [];
        disp('V1 L spike shift shuffle')
        tic
        for iBoot = 1:1000
            temp1 = [];
            temp2 = [];
            parfor event_id = 1:size(binnedArrayUPV1{1},1)
                s = RandStream('mrg32k3a','Seed',100000*iBoot+event_id); % Set random seed for resampling
                %             s = RandStream('mrg32k3a','Seed',i+10000*shuffle_options); % Set random seed for resampling
                bins_to_shift = datasample(s,1:1:size(binnedArrayUPV1{1},2),1);
                bins= circshift(1:1:size(binnedArrayUPV1{1},2),bins_to_shift);
                temp1(event_id,:) = binnedArrayUPV1{1}(event_id,bins);
                temp2(event_id,:) = binnedArrayDOWNV1{1}(event_id,bins);
            end
            tempUP(iBoot,:) = mean(temp1,1,'omitnan');
            tempDOWN(iBoot,:) = mean(temp2,1,'omitnan');

            temp1 = [];
            parfor event_id = 1:size(binnedArrayRipplesV1{1},1)
                s = RandStream('mrg32k3a','Seed',100000*iBoot+event_id); % Set random seed for resampling
                %             s = RandStream('mrg32k3a','Seed',i+10000*shuffle_options); % Set random seed for resampling
                bins_to_shift = datasample(s,1:1:size(binnedArrayRipplesV1{1},2),1);
                bins= circshift(1:1:size(binnedArrayRipplesV1{1},2),bins_to_shift);
                temp1(event_id,:) = binnedArrayRipplesV1{1}(event_id,bins);
            end
            tempRipples(iBoot,:) = mean(temp1,1,'omitnan');

        end
        toc

        PSTH(nprobe).L_V1_UP_shuffled = tempUP;
        PSTH(nprobe).L_V1_DOWN_shuffled = tempDOWN;
        PSTH(nprobe).L_V1_ripples_shuffled = tempRipples;
    end


    % V1 R
    % bootstrap distribution
    tempUP = [];
    tempDOWN = [];
    tempRipples = [];
    for iBoot = 1:1000
        s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
        event_id = datasample(s,1:size(binnedArrayUPV1{2},1),size(binnedArrayUPV1{2},1));
        tempUP(iBoot,:) = mean(binnedArrayUPV1{2}(event_id,:),'omitnan');

        s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
        event_id = datasample(s,1:size(binnedArrayDOWNV1{2},1),size(binnedArrayDOWNV1{2},1));
        tempDOWN(iBoot,:) = mean(binnedArrayDOWNV1{2}(event_id,:),'omitnan');

        s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
        event_id = datasample(s,1:size(binnedArrayRipplesV1{2},1),size(binnedArrayRipplesV1{2},1));
        tempRipples(iBoot,:) = mean(binnedArrayRipplesV1{2}(event_id,:),'omitnan');
    end

    PSTH(nprobe).R_V1_UP_bootstrap = tempUP;
    PSTH(nprobe).R_V1_DOWN_bootstrap = tempDOWN;
    PSTH(nprobe).R_V1_ripples_bootstrap = tempRipples;

    if contains(shuffle_option,'time_circular_shift')
        % timebin circularly shifted
        tempUP = [];
        tempDOWN = [];
        tempRipples = [];
        disp('V1 R spike shift shuffle')
        tic
        for iBoot = 1:1000
            temp1 = [];
            temp2 = [];
            parfor event_id = 1:size(binnedArrayUPV1{2},1)
                s = RandStream('mrg32k3a','Seed',100000*iBoot+event_id); % Set random seed for resampling
                %             s = RandStream('mrg32k3a','Seed',i+10000*shuffle_options); % Set random seed for resampling
                bins_to_shift = datasample(s,1:1:size(binnedArrayUPV1{2},2),1);
                bins= circshift(1:1:size(binnedArrayUPV1{2},2),bins_to_shift);
                temp1(event_id,:) = binnedArrayUPV1{2}(event_id,bins);
                temp2(event_id,:) = binnedArrayDOWNV1{2}(event_id,bins);
            end
            tempUP(iBoot,:) = mean(temp1,1,'omitnan');
            tempDOWN(iBoot,:) = mean(temp2,1,'omitnan');

            temp1 = [];
            parfor event_id = 1:size(binnedArrayRipplesV1{2},1)
                s = RandStream('mrg32k3a','Seed',100000*iBoot+event_id); % Set random seed for resampling
                %             s = RandStream('mrg32k3a','Seed',i+10000*shuffle_options); % Set random seed for resampling
                bins_to_shift = datasample(s,1:1:size(binnedArrayRipplesV1{2},2),1);
                bins= circshift(1:1:size(binnedArrayRipplesV1{2},2),bins_to_shift);
                temp1(event_id,:) = binnedArrayRipplesV1{2}(event_id,bins);
            end
            tempRipples(iBoot,:) = mean(temp1,1,'omitnan');

        end
        toc

        PSTH(nprobe).R_V1_UP_shuffled = tempUP;
        PSTH(nprobe).R_V1_DOWN_shuffled = tempDOWN;
        PSTH(nprobe).R_V1_ripples_shuffled = tempRipples;
    end



    % HPC L
    % bootstrap distribution
    tempUP = [];
    tempDOWN = [];
    tempRipples = [];
    for iBoot = 1:1000
        s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
        event_id = datasample(s,1:size(binnedArrayUPHPC{1},1),size(binnedArrayUPHPC{1},1));
        tempUP(iBoot,:) = mean(binnedArrayUPHPC{1}(event_id,:),'omitnan');

        s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
        event_id = datasample(s,1:size(binnedArrayDOWNHPC{1},1),size(binnedArrayDOWNHPC{1},1));
        tempDOWN(iBoot,:) = mean(binnedArrayDOWNHPC{1}(event_id,:),'omitnan');

        s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
        event_id = datasample(s,1:size(binnedArrayRipplesHPC{1},1),size(binnedArrayRipplesHPC{1},1));
        tempRipples(iBoot,:) = mean(binnedArrayRipplesHPC{1}(event_id,:),'omitnan');
    end

    PSTH(nprobe).L_HPC_UP_bootstrap = tempUP;
    PSTH(nprobe).L_HPC_DOWN_bootstrap = tempDOWN;
    PSTH(nprobe).L_HPC_ripples_bootstrap = tempRipples;

    % timebin circularly shifted
    if contains(shuffle_option,'time_circular_shift')
        tempUP = [];
        tempDOWN = [];
        tempRipples = [];
        disp('V1 L spike shift shuffle')
        tic
        for iBoot = 1:1000
            temp1 = [];
            temp2 = [];
            parfor event_id = 1:size(binnedArrayUPHPC{1},1)
                s = RandStream('mrg32k3a','Seed',100000*iBoot+event_id); % Set random seed for resampling
                %             s = RandStream('mrg32k3a','Seed',i+10000*shuffle_options); % Set random seed for resampling
                bins_to_shift = datasample(s,1:1:size(binnedArrayUPHPC{1},2),1);
                bins= circshift(1:1:size(binnedArrayUPHPC{1},2),bins_to_shift);
                temp1(event_id,:) = binnedArrayUPHPC{1}(event_id,bins);
                temp2(event_id,:) = binnedArrayDOWNHPC{1}(event_id,bins);
            end
            tempUP(iBoot,:) = mean(temp1,1,'omitnan');
            tempDOWN(iBoot,:) = mean(temp2,1,'omitnan');

            temp1 = [];
            parfor event_id = 1:size(binnedArrayRipplesHPC{1},1)
                s = RandStream('mrg32k3a','Seed',100000*iBoot+event_id); % Set random seed for resampling
                %             s = RandStream('mrg32k3a','Seed',i+10000*shuffle_options); % Set random seed for resampling
                bins_to_shift = datasample(s,1:1:size(binnedArrayRipplesHPC{1},2),1);
                bins= circshift(1:1:size(binnedArrayRipplesHPC{1},2),bins_to_shift);
                temp1(event_id,:) = binnedArrayRipplesHPC{1}(event_id,bins);
            end
            tempRipples(iBoot,:) = mean(temp1,1,'omitnan');

        end
        toc

        PSTH(nprobe).L_HPC_UP_shuffled = tempUP;
        PSTH(nprobe).L_HPC_DOWN_shuffled = tempDOWN;
        PSTH(nprobe).L_HPC_ripples_shuffled = tempRipples;
    end


    % HPC R
    % bootstrap distribution
    tempUP = [];
    tempDOWN = [];
    tempRipples = [];
    for iBoot = 1:1000
        s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
        event_id = datasample(s,1:size(binnedArrayUPHPC{2},1),size(binnedArrayUPHPC{2},1));
        tempUP(iBoot,:) = mean(binnedArrayUPHPC{2}(event_id,:),'omitnan');

        s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
        event_id = datasample(s,1:size(binnedArrayDOWNHPC{2},1),size(binnedArrayDOWNHPC{2},1));
        tempDOWN(iBoot,:) = mean(binnedArrayDOWNHPC{2}(event_id,:),'omitnan');

        s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
        event_id = datasample(s,1:size(binnedArrayRipplesHPC{2},1),size(binnedArrayRipplesHPC{2},1));
        tempRipples(iBoot,:) = mean(binnedArrayRipplesHPC{2}(event_id,:),'omitnan');
    end

    PSTH(nprobe).R_HPC_UP_bootstrap = tempUP;
    PSTH(nprobe).R_HPC_DOWN_bootstrap = tempDOWN;
    PSTH(nprobe).R_HPC_ripples_bootstrap = tempRipples;

    if contains(shuffle_option,'time_circular_shift')
        % timebin circularly shifted
        tempUP = [];
        tempDOWN = [];
        tempRipples = [];
        disp('HPC R spike shift shuffle')
        tic
        for iBoot = 1:1000
            temp1 = [];
            temp2 = [];
            parfor event_id = 1:size(binnedArrayUPHPC{2},1)
                s = RandStream('mrg32k3a','Seed',100000*iBoot+event_id); % Set random seed for resampling
                %             s = RandStream('mrg32k3a','Seed',i+10000*shuffle_options); % Set random seed for resampling
                bins_to_shift = datasample(s,1:1:size(binnedArrayUPHPC{2},2),1);
                bins= circshift(1:1:size(binnedArrayUPHPC{2},2),bins_to_shift);
                temp1(event_id,:) = binnedArrayUPHPC{2}(event_id,bins);
                temp2(event_id,:) = binnedArrayDOWNHPC{2}(event_id,bins);
            end
            tempUP(iBoot,:) = mean(temp1,1,'omitnan');
            tempDOWN(iBoot,:) = mean(temp2,1,'omitnan');

            temp1 = [];
            parfor event_id = 1:size(binnedArrayRipplesHPC{2},1)
                s = RandStream('mrg32k3a','Seed',100000*iBoot+event_id); % Set random seed for resampling
                %             s = RandStream('mrg32k3a','Seed',i+10000*shuffle_options); % Set random seed for resampling
                bins_to_shift = datasample(s,1:1:size(binnedArrayRipplesHPC{2},2),1);
                bins= circshift(1:1:size(binnedArrayRipplesHPC{2},2),bins_to_shift);
                temp1(event_id,:) = binnedArrayRipplesHPC{2}(event_id,bins);
            end
            tempRipples(iBoot,:) = mean(temp1,1,'omitnan');

        end
        toc

        PSTH(nprobe).R_HPC_UP_shuffled = tempUP;
        PSTH(nprobe).R_HPC_DOWN_shuffled = tempDOWN;
        PSTH(nprobe).R_HPC_ripples_shuffled = tempRipples;
    end
end