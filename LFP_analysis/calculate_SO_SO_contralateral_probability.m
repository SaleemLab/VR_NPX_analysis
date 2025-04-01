function probability = calculate_SO_SO_contralateral_probability(slow_waves_all,sessions_to_process,varargin)
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

for nprobe = 1:length(slow_waves_all)
    if nprobe == 1
        mprobe = 2;
    else
        mprobe = 1;
    end
    %%%%%%%%%%%%%%% L ripples
    UP_index_all = [];
    UP_index_all1 = [];

    DOWN_index_all = [];
    DOWN_index_all1 = [];


    binnedArrayUU = [];
    binnedArrayDD = [];
    binnedArrayUD = [];
    binnedArrayDU = [];

    for nsession = 1:length(sessions_to_process)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ipsilateral UP and DOWN
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


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% contralateral UP and DOWN
        % Find UP followed by a DOWN 
        %         UP_DOWN_index = find(slow_waves_all(nprobe).UP_session_count == nsession); % Find UP -> DOWN
        UP_index1 = find(slow_waves_all(mprobe).UP_session_count == sessions_to_process(nsession)); % Find UP this session
        %         UP_index = UP_index(slow_waves_all(nprobe).UP_DOWN_index(UP_DOWN_index,1)); % Find UP followed by a DOWN
        UP_duration1 = slow_waves_all(mprobe).UP_ints(UP_index1,2)-slow_waves_all(mprobe).UP_ints(UP_index1,1);
        UP_index1 = UP_index1(UP_duration1<=duration_threshold);
        UP_ints1 = slow_waves_all(mprobe).UP_ints(UP_index1,:);


        % Find DOWN index
        DOWN_index1 = find(slow_waves_all(mprobe).DOWN_session_count == sessions_to_process(nsession)); % Find DOWN this session
        DOWN_ints1 = slow_waves_all(mprobe).DOWN_ints(DOWN_index1,:);
        [C,ia,ib] = intersect(UP_ints1(:,2),DOWN_ints1(:,1));
        % [C,ia,ib] = intersect(UP_ints(:,2)+0.01,DOWN_ints(:,1));
        UP_index1 = UP_index1(ia);
        DOWN_index1 = DOWN_index1(ib);
        
        DOWN_ints1 = slow_waves_all(mprobe).DOWN_ints(DOWN_index1,:);
        DOWN_index_all1 = [DOWN_index_all1; DOWN_index1];

        UP_ints1 = slow_waves_all(mprobe).UP_ints(UP_index1,:);
        UP_index_all1 = [UP_index_all1; UP_index1];


        if contains(shuffle_option,'baseline')
            % s = RandStream('mrg32k3a','Seed',1); % Set random seed for resampling
            % time_jitter = 2 + (2.5 - 2) * rand(s,1, length(UP_index));
            % time_jitter = [time_jitter' time_jitter'];
            time_jitter = [3*ones(1,length(UP_index))' 3*ones(1,length(UP_index))'];
            UP_ints = slow_waves_all(nprobe).UP_ints(UP_index,:)-time_jitter;

            % s = RandStream('mrg32k3a','Seed',2); % Set random seed for resampling
            % time_jitter = 2 + (2.5 - 2) * rand(s,1, length(DOWN_index));
            % time_jitter = [time_jitter' time_jitter'];
            time_jitter = [3*ones(1,length(DOWN_index))' 3*ones(1,length(DOWN_index))'];
            DOWN_ints = slow_waves_all(nprobe).DOWN_ints(DOWN_index,:)-time_jitter;
        else
            UP_ints = slow_waves_all(nprobe).UP_ints(UP_index,:);
            DOWN_ints = slow_waves_all(nprobe).DOWN_ints(DOWN_index,:);
        end

        if contains(time_option,'onset')
            time_index= 1;
        else
            time_index= 1:2;
        end

         % Probability of UP During U-D
        if contains(option,'normalised')
            [~,~,~,temp] = calculate_relative_event_probability(DOWN_ints,UP_ints1(:,time_index),num_bins,0);
        else
            [~,temp,~] = calculate_event_probability(UP_ints1(:,time_index),DOWN_ints(:,1),time_windows(1):time_bin:time_windows(end),0);

            timebin_edges_all = DOWN_ints(:,1) + bins_centre;  % Absolute times of peri-event window

            if ~contains(shuffle_option,'baseline')
                for i = 1:size(DOWN_ints,1)

                    timebin_edges_all(i,:);
                    % Previous DOWN (skip if this is the first DOWN)
                    if i > 1
                        prev_offset = DOWN_ints(i-1,2);
                        % Find peri-time indices within the previous DOWN state
                        mask_prev =  timebin_edges_all(i,:) <= prev_offset;
                        temp(i, mask_prev) = NaN;
                    end

                    % Next DOWN (skip if this is the last DOWN)
                    if i < size(DOWN_ints,1)
                        next_onset = DOWN_ints(i+1,1);
                        % Find peri-time indices within the next DOWN state
                        mask_next = timebin_edges_all(i,:) >= next_onset;
                        temp(i, mask_next) = NaN;
                    end
                end
            end
        end

        binnedArrayUD=[binnedArrayUD; temp];

        % Probability of DOWN During D-U
        if contains(option,'normalised')
            [~,~,~,temp] = calculate_relative_event_probability(UP_ints,DOWN_ints1(:,time_index),num_bins,0);
        else
            [~,temp,~] = calculate_event_probability(DOWN_ints1(:,time_index),UP_ints(:,1),time_windows(1):time_bin:time_windows(end),0);

            timebin_edges_all = UP_ints(:,1) + bins_centre;  % Absolute times of peri-event window

            if ~contains(shuffle_option,'baseline')
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
        end

        binnedArrayDU=[binnedArrayDU; temp];

        % Probability of UP During D-U
        if contains(option,'normalised')
            [~,~,~,temp] = calculate_relative_event_probability(UP_ints,UP_ints1(:,time_index),num_bins,0);
        else
            [~,temp,~] = calculate_event_probability(UP_ints1(:,time_index),UP_ints(:,1),time_windows(1):time_bin:time_windows(end),0);

            timebin_edges_all = UP_ints(:,1) + bins_centre;  % Absolute times of peri-event window

            if ~contains(shuffle_option,'baseline')
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
        end

        binnedArrayUU=[binnedArrayUU; temp];
        
        % Probability of DOWN During U-D
        if contains(option,'normalised')
            [~,~,~,temp] = calculate_relative_event_probability(DOWN_ints,DOWN_ints1(:,time_index),num_bins,0);
        else
            [~,temp,~] = calculate_event_probability(DOWN_ints1(:,time_index),DOWN_ints(:,1),time_windows(1):time_bin:time_windows(end),0);

            timebin_edges_all = DOWN_ints(:,1) + bins_centre;  % Absolute times of peri-event window

            if ~contains(shuffle_option,'baseline')
                for i = 1:size(DOWN_ints,1)

                    timebin_edges_all(i,:);
                    % Previous UP (skip if this is the first UP)
                    if i > 1
                        prev_offset = DOWN_ints(i-1,2);
                        % Find peri-time indices within the previous UP state
                        mask_prev =  timebin_edges_all(i,:) <= prev_offset;
                        temp(i, mask_prev) = NaN;
                    end

                    % Next UP (skip if this is the last UP)
                    if i < size(DOWN_ints,1)
                        next_onset = DOWN_ints(i+1,1);
                        % Find peri-time indices within the next UP state
                        mask_next = timebin_edges_all(i,:) >= next_onset;
                        temp(i, mask_next) = NaN;
                    end
                end
            end
        end

        binnedArrayDD=[binnedArrayDD; temp];
    end
    probability(nprobe).UP_all_index = UP_index_all;
    probability(nprobe).DOWN_all_index = DOWN_index_all;
    probability(nprobe).UP_DOWN = binnedArrayUD; % UP during D-U
    probability(nprobe).DOWN_UP = binnedArrayDU; % DOWN during U-D
    probability(nprobe).UP_UP = binnedArrayUU; % UP during D-U
    probability(nprobe).DOWN_DOWN = binnedArrayDD; % DOWN during U-D

    all_spindle_no = sum(spindles_all(1).SWS_index == 1);
    probability(nprobe).L_spindle_no = all_spindle_no;
    all_UP_no = length(UP_index_all);
    all_DOWN_no = length(DOWN_index_all);

    % bootstrap distribution
    tempUD = [];
    tempDU = [];
    tempUU = [];
    tempDD = [];
    for iBoot = 1:1000
        s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
        event_id = datasample(s,1:size(probability(nprobe).UP_DOWN,1),size(probability(nprobe).UP_DOWN,1));
        tempUD(iBoot,:) = sum(probability(nprobe).UP_DOWN(event_id,:),'omitnan')./sum(~isnan(probability(nprobe).UP_DOWN(event_id,:)));

        s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
        event_id = datasample(s,1:size(probability(nprobe).DOWN_UP,1),size(probability(nprobe).DOWN_UP,1));
        tempDU(iBoot,:) = sum(probability(nprobe).DOWN_UP(event_id,:),'omitnan')./sum(~isnan(probability(nprobe).DOWN_UP(event_id,:)));

        s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
        event_id = datasample(s,1:size(probability(nprobe).UP_UP,1),size(probability(nprobe).UP_UP,1));
        tempUU(iBoot,:) = sum(probability(nprobe).UP_UP(event_id,:),'omitnan')./sum(~isnan(probability(nprobe).UP_UP(event_id,:)));

        s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
        event_id = datasample(s,1:size(probability(nprobe).DOWN_DOWN,1),size(probability(nprobe).DOWN_DOWN,1));
        tempDD(iBoot,:) = sum(probability(nprobe).DOWN_DOWN(event_id,:),'omitnan')./sum(~isnan(probability(nprobe).DOWN_DOWN(event_id,:)));
    end

    probability(nprobe).UP_DOWN_bootstrap = tempUD; % UP during D-U
    probability(nprobe).DOWN_UP_bootstrap = tempDU; % DOWN during U-D
    probability(nprobe).UP_UP_bootstrap = tempUU; % UP during D-U
    probability(nprobe).DOWN_DOWN_bootstrap = tempDD; % UP during D-U

    if contains(shuffle_option,'time_circular_shift')
        % timebin circularly shifted
        tempUD = [];
        tempDU = [];
        tempUU = [];
        tempDD = [];
        tic
        for iBoot = 1:1000
            temp1 = [];
            temp2 = [];
            temp3 = [];
            temp4 = [];
            parfor event_id = 1:size(probability(nprobe).UP_DOWN,1)
                s = RandStream('mrg32k3a','Seed',100000*iBoot+event_id); % Set random seed for resampling
                %             s = RandStream('mrg32k3a','Seed',i+10000*shuffle_options); % Set random seed for resampling
                bins_to_shift = datasample(s,1:1:size(probability(nprobe).UP_DOWN,2),1);
                bins= circshift(1:1:size(probability(nprobe).UP_DOWN,2),bins_to_shift);

                temp1(event_id,:) = probability(nprobe).UP_DOWN(event_id,bins);
                temp2(event_id,:) = probability(nprobe).DOWN_UP(event_id,bins);
                temp3(event_id,:) = probability(nprobe).UP_UP(event_id,bins);
                temp4(event_id,:) = probability(nprobe).DOWN_DOWN(event_id,bins);
            end

            tempUD(iBoot,:) = sum(temp1,1,'omitnan')./sum(~isnan(temp1));
            tempDU(iBoot,:) = sum(temp2,1,'omitnan')./sum(~isnan(temp2));
            tempUU(iBoot,:) = sum(temp3,1,'omitnan')./sum(~isnan(temp3));
            tempDD(iBoot,:) = sum(temp4,1,'omitnan')./sum(~isnan(temp4));
        end
        toc
        probability(nprobe).UP_DOWN_shuffled = tempUD; % UP during D-U
        probability(nprobe).DOWN_UP_shuffled = tempDU; % DOWN during U-D
        probability(nprobe).UP_UP_shuffled = tempUU; % UP during D-U
        probability(nprobe).DOWN_DOWN_shuffled = tempDD; % UP during D-U
    end
end