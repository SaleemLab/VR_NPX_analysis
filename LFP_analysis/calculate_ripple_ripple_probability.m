function probability = calculate_ripple_ripple_probability(slow_waves_all,sessions_to_process,varargin)
% ripples_all = [];
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
ripples_all = [];
for nprobe = 1:length(slow_waves_all)
    UP_index_all = [];
    DOWN_index_all = [];
    ripple_index_UP_all = [];
    ripple_index_DOWN_all = [];
    ripple_index_all=[];

    binnedArray = [];
    binnedArrayShuffled = [];

    for nsession = 1:length(sessions_to_process)

       ripples_index = find(ripples_all(1).session_count == sessions_to_process(nsession)& ripples_all(1).SWS_index == 1);
        ripple_peaktimes = min(ripples_all(1).SWR_peaktimes{sessions_to_process(nsession)}(ripples_all(1).probe_hemisphere{sessions_to_process(nsession)} == 1,ripples_all(1).SWS_index(ripples_all(1).session_count == sessions_to_process(nsession))==1))';

        ripple_peaktimes1 = min(ripples_all(nprobe).SWR_peaktimes{sessions_to_process(nsession)}(ripples_all(nprobe).probe_hemisphere{sessions_to_process(nsession)} == nprobe,ripples_all(nprobe).SWS_index(ripples_all(nprobe).session_count == sessions_to_process(nsession))==1))';
        % ripple_peaktimes1 = min(ripples_all(nprobe).SWR_peaktimes{sessions_to_process(nsession)}(ripples_all(nprobe).probe_hemisphere{sessions_to_process(nsession)} == 2,ripples_all(nprobe).SWS_index(ripples_all(nprobe).session_count == sessions_to_process(nsession))==1))';
        % if contains(time_option,'peaktimes')
        % else
            % ripple_times = [ripples_all(1).onset(ripples_index) ripples_all(1).offset(ripples_index)];
        % end
        if contains(time_option,'peaktimes')
            ripple_times1= ripple_peaktimes1;
        else
            ripple_times1 = [ripples_all(nprobe).onset(ripples_index) ripples_all(nprobe).offset(ripples_index)];
        end

        if contains(shuffle_option,'baseline')
            % s = RandStream('mrg32k3a','Seed',1); % Set random seed for resampling
            % time_jitter = 2 + (2.5 - 2) * rand(s,1, length(UP_index));
            % time_jitter = [time_jitter' time_jitter'];
            time_jitter = [3*ones(1,length(ripples_index))' 3*ones(1,length(ripples_index))'];

            % if contains(time_option,'peaktimes')
            time_jitter = 3*ones(1,length(ripples_index))';
            % else
            %     time_jitter = [3*ones(1,length(ripples_index))' 3*ones(1,length(ripples_index))'];
            % end

            ripple_peaktimes = ripple_peaktimes-time_jitter;
        end


        % Probability of L ripple
        [probability(nprobe).L_ripples_session(nsession,:),temp,event_index] = calculate_event_probability(ripple_times1,ripple_peaktimes,time_windows(1):time_bin:time_windows(end),0);
        % end


        if ~contains(shuffle_option,'baseline')
            ripple_times = [ripples_all(1).onset(ripples_index) ripples_all(1).offset(ripples_index)];
            timebin_edges_all = ripple_times + bins_centre;  % Absolute times of peri-event window

            for i = 1:size(ripple_times,1)

                timebin_edges_all(i,:);
                % Previous DOWN (skip if this is the first DOWN)
                if i > 1
                    prev_offset = ripple_times(i-1,2);
                    % Find peri-time indices within the previous DOWN state
                    mask_prev =  timebin_edges_all(i,:) <= prev_offset;
                    temp(i, mask_prev) = NaN;
                end

                % Next DOWN (skip if this is the last DOWN)
                if i < size(DOWN_ints,1)
                    next_onset = ripple_times(i+1,1);
                    % Find peri-time indices within the next DOWN state
                    mask_next = timebin_edges_all(i,:) >= next_onset;
                    temp(i, mask_next) = NaN;
                end
            end
        end

        binnedArray=[binnedArray; temp];
        ripple_index_all = [ripple_index_all;event_index];
    end

    probability(nprobe).L_ripples = binnedArray;
    probability(nprobe).L_ripples_index = ripple_index_all;% probability of ripples relative to Left ripple onset
    
    all_ripple_no = sum(ripples_all(1).SWS_index == 1);
    probability(nprobe).L_ripple_no = all_ripple_no;


    % bootstrap distribution
    tempUP = [];
    for iBoot = 1:1000
        s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
        event_id = datasample(s,1:size(binnedArray,1),size(binnedArray,1));
        tempUP(iBoot,:) =  sum(binnedArray(event_id,:),'omitnan')./sum(~isnan(binnedArray(event_id,:)));
    end

    probability(nprobe).L_ripples_bootstrap = tempUP;

    if contains(shuffle_option,'time_circular_shift')

        % timebin circularly shifted
        tempUP = [];
        tic
        for iBoot = 1:1000
            temp1 = [];
            temp2 = [];
            parfor event_id = 1:size(binnedArray,1)
                s = RandStream('mrg32k3a','Seed',100000*iBoot+event_id); % Set random seed for resampling
                %             s = RandStream('mrg32k3a','Seed',i+10000*shuffle_options); % Set random seed for resampling
                bins_to_shift = datasample(s,1:1:size(binnedArray,2),1);
                bins= circshift(1:1:size(binnedArray,2),bins_to_shift);
                temp1(event_id,:) = binnedArray(event_id,bins);
            end

            tempUP(iBoot,:) = sum(temp1,1,'omitnan')./sum(~isnan(temp1));

        end
        toc

        probability(nprobe).L_ripples_shuffled = tempUP;
    end



    
    %%%%%%% R Ripple 
    ripple_index_all=[];
    binnedArray = [];

    for nsession = 1:length(sessions_to_process)

        ripples_index = find(ripples_all(2).session_count == sessions_to_process(nsession)& ripples_all(2).SWS_index == 1);
        ripple_peaktimes = min(ripples_all(2).SWR_peaktimes{sessions_to_process(nsession)}(ripples_all(2).probe_hemisphere{sessions_to_process(nsession)} == 2,ripples_all(2).SWS_index(ripples_all(2).session_count == sessions_to_process(nsession))==1))';

        ripple_peaktimes1 = min(ripples_all(nprobe).SWR_peaktimes{sessions_to_process(nsession)}(ripples_all(nprobe).probe_hemisphere{sessions_to_process(nsession)} == nprobe,ripples_all(nprobe).SWS_index(ripples_all(nprobe).session_count == sessions_to_process(nsession))==1))';
        % ripple_peaktimes1 = min(ripples_all(nprobe).SWR_peaktimes{sessions_to_process(nsession)}(ripples_all(nprobe).probe_hemisphere{sessions_to_process(nsession)} == 2,ripples_all(nprobe).SWS_index(ripples_all(nprobe).session_count == sessions_to_process(nsession))==1))';
        % if contains(time_option,'peaktimes')
        % else
        % ripple_times = [ripples_all(1).onset(ripples_index) ripples_all(1).offset(ripples_index)];
        % end
        if contains(time_option,'peaktimes')
            ripple_times1= ripple_peaktimes1;
        else
            ripple_times1 = [ripples_all(nprobe).onset(ripples_index) ripples_all(nprobe).offset(ripples_index)];
        end

        if contains(shuffle_option,'baseline')
            % s = RandStream('mrg32k3a','Seed',1); % Set random seed for resampling
            % time_jitter = 2 + (2.5 - 2) * rand(s,1, length(UP_index));
            % time_jitter = [time_jitter' time_jitter'];
            time_jitter = [3*ones(1,length(ripples_index))' 3*ones(1,length(ripples_index))'];

            % if contains(time_option,'peaktimes')
            time_jitter = 3*ones(1,length(ripples_index))';
            % else
            %     time_jitter = [3*ones(1,length(ripples_index))' 3*ones(1,length(ripples_index))'];
            % end

            ripple_peaktimes = ripple_peaktimes-time_jitter;
        end


        % Probability of R ripple
        [probability(nprobe).R_ripples_session(nsession,:),temp,event_index] = calculate_event_probability(ripple_times1,ripple_peaktimes,time_windows(1):time_bin:time_windows(end),0);
        % end


        if ~contains(shuffle_option,'baseline')
            ripple_times = [ripples_all(2).onset(ripples_index) ripples_all(2).offset(ripples_index)];
            timebin_edges_all = ripple_times + bins_centre;  % Absolute times of peri-event window

            for i = 1:size(ripple_times,1)

                timebin_edges_all(i,:);
                % Previous DOWN (skip if this is the first DOWN)
                if i > 1
                    prev_offset = ripple_times(i-1,2);
                    % Find peri-time indices within the previous DOWN state
                    mask_prev =  timebin_edges_all(i,:) <= prev_offset;
                    temp(i, mask_prev) = NaN;
                end

                % Next DOWN (skip if this is the last DOWN)
                if i < size(DOWN_ints,1)
                    next_onset = ripple_times(i+1,1);
                    % Find peri-time indices within the next DOWN state
                    mask_next = timebin_edges_all(i,:) >= next_onset;
                    temp(i, mask_next) = NaN;
                end
            end
        end

        binnedArray=[binnedArray; temp];
        ripple_index_all = [ripple_index_all;event_index];
    end

    probability(nprobe).R_ripples = binnedArray;
    probability(nprobe).R_ripples_index = ripple_index_all;% probability of ripples relative to Left ripple onset

    all_ripple_no = sum(ripples_all(2).SWS_index == 1);
    probability(nprobe).R_ripple_no = all_ripple_no;


    % bootstrap distribution
    tempUP = [];
    for iBoot = 1:1000
        s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
        event_id = datasample(s,1:size(binnedArray,1),size(binnedArray,1));
        tempUP(iBoot,:) =  sum(binnedArray(event_id,:),'omitnan')./sum(~isnan(binnedArray(event_id,:)));
    end

    probability(nprobe).R_ripples_bootstrap = tempUP;

    if contains(shuffle_option,'time_circular_shift')

        % timebin circularly shifted
        tempUP = [];
        tic
        for iBoot = 1:1000
            temp1 = [];
            temp2 = [];
            parfor event_id = 1:size(binnedArray,1)
                s = RandStream('mrg32k3a','Seed',100000*iBoot+event_id); % Set random seed for resampling
                %             s = RandStream('mrg32k3a','Seed',i+10000*shuffle_options); % Set random seed for resampling
                bins_to_shift = datasample(s,1:1:size(binnedArray,2),1);
                bins= circshift(1:1:size(binnedArray,2),bins_to_shift);
                temp1(event_id,:) = binnedArray(event_id,bins);
            end

            tempUP(iBoot,:) = sum(temp1,1,'omitnan')./sum(~isnan(temp1));

        end
        toc

        probability(nprobe).R_ripples_shuffled = tempUP;
    end


end