%% read training data
function behaviour = read_training_data(path,date)
% date = ['20',date(5:6),'-',date(3:4),'-',date(1:2)];
reward_path = dir(fullfile(path,['*Reward',date,'*.csv']));
lick_performance_path = dir(fullfile(path,['*Lick_Performance',date,'*.csv']));
wheel_log_path = dir(fullfile(path,['*WheelLog',date,'*.csv']));
trial_path = dir(fullfile(path,['*Trial_info',date,'*.csv']));

if ~contains(reward_path(1).name,'Manual')
    reward = readtable([path,reward_path(1).name]);
else
    reward = readtable([path,reward_path(2).name]);
end
lick_performance = readtable([path,lick_performance_path(1).name]);
if ~isempty(wheel_log_path)
    wheel_log = readtable([path,wheel_log_path(1).name]);
end

if ~isempty(trial_path)
    trial_info = readtable([path,trial_path(1).name]);
end

%% wheel log
if ~isempty(wheel_log_path)
    %% wheel frame
    behaviour.wheel_frame = table2array(wheel_log(:,"Var1"));

    %% wheel Time
    behaviour.wheel_time = table2array(wheel_log(:,"Var2"));

    %% wheel computer time
    behaviour.wheel_computer_time = table2cell(wheel_log(:,"Var3"));

    %% wheel raw input
    behaviour.wheel_raw_input = table2cell(wheel_log(:,"Var4"));
    behaviour.wheel_raw_input = cellfun(@(x) x(8:end), behaviour.wheel_raw_input, 'UniformOutput',false);
    behaviour.wheel_raw_input = cellfun(@str2double,behaviour.wheel_raw_input);

    %% wheel left lick
    behaviour.wheel_lick_L = table2cell(wheel_log(:,"Var5"));
    behaviour.wheel_lick_L = cellfun(@(x) x(7:end), behaviour.wheel_lick_L, 'UniformOutput',false);
    behaviour.wheel_lick_L = cellfun(@str2double,behaviour.wheel_lick_L);

    %% wheel right lick
    behaviour.wheel_lick_R = table2cell(wheel_log(:,"Var6"));
    behaviour.wheel_lick_R = cellfun(@(x) x(7:end), behaviour.wheel_lick_R, 'UniformOutput',false);
    behaviour.wheel_lick_R = cellfun(@str2double,behaviour.wheel_lick_R);

    %% wheel quad state
    behaviour.wheel_quad_state = table2cell(wheel_log(:,"Var7"));
    behaviour.wheel_quad_state = cellfun(@(x) x(11:end), behaviour.wheel_quad_state, 'UniformOutput',false);
    behaviour.wheel_quad_state = cellfun(@str2double,behaviour.wheel_quad_state);
    %% wheel camera frame count
    behaviour.wheel_eye_frame_count = table2cell(wheel_log(:,"Var8"));
    behaviour.wheel_eye_frame_count = cellfun(@(x) x(15:end), behaviour.wheel_eye_frame_count, 'UniformOutput',false);
    behaviour.wheel_eye_frame_count = cellfun(@str2double,behaviour.wheel_eye_frame_count);

    behaviour.wheel_reference_timestamp = table2cell(wheel_log(:,"Var9"));
    behaviour.wheel_reference_timestamp = cellfun(@(x) x(6:end), behaviour.wheel_reference_timestamp, 'UniformOutput',false);
    behaviour.wheel_reference_timestamp = cellfun(@str2double,behaviour.wheel_reference_timestamp);
   
    %% wheel timestamp of serial string read (time of day, total millisecond)
    behaviour.wheel_arduino_read_time = table2cell(wheel_log(:,"Var9"));
    behaviour.wheel_arduino_read_time = cellfun(@(x) x(6:end), behaviour.wheel_arduino_read_time, 'UniformOutput',false);
    behaviour.wheel_arduino_read_time = cellfun(@str2double,behaviour.wheel_arduino_read_time);

    %% wheel position
    behaviour.wheel_position = table2cell(wheel_log(:,"Var10"));
    behaviour.wheel_position = cellfun(@(x) x(10:end), behaviour.wheel_position, 'UniformOutput',false);
    behaviour.wheel_position = cellfun(@(x) x(1:end-1), behaviour.wheel_position, 'UniformOutput',false);
    behaviour.wheel_position = cellfun(@str2double,behaviour.wheel_position);
end

%% reward type, 1 = passive, 10 = active, non if fail the trial (not licking for active, aborted both the same)

behaviour.reward_type = table2cell(reward(:,"Var1"));
behaviour.reward_type = cellfun(@(x) x(2:end), behaviour.reward_type, 'UniformOutput',false);
behaviour.reward_type = cellfun(@str2double,behaviour.reward_type);

% If first datapoint is due to initialisation 
if behaviour.reward_type(1) == 0
    behaviour.reward_type(1) = [];
    reward(1,:) = [];
end
    
%%  Reward Delivery Time Based on Eye Frame Count, or Time

if lick_performance_path.datenum > datenum('14-June-2023 13:00:00')
    reward_delivery_time = table2array(reward(:,"Var2")); % May want to simplify this in the future
    
    for n = 1:length(reward_delivery_time)
        start_index = find(reward_delivery_time(n) == behaviour.wheel_reference_timestamp);
        
        if isempty(start_index) % in rare cases where a frame is dropped, use next frame to get the lick time
            for count = 1:10
                start_index = find(reward_delivery_time(n)+count == behaviour.wheel_reference_timestamp);
                if ~isempty(start_index)
                    break
                end
            end
        end
        
        behaviour.reward_delivery_time(n,1) = behaviour.wheel_time(start_index);
    end
else
    behaviour.reward_delivery_eye_frame_count = table2array(reward(:,"Var2"));
    
    for n = 1:length(behaviour.reward_delivery_eye_frame_count)
        start_index = find(behaviour.reward_delivery_eye_frame_count(n) == behaviour.wheel_eye_frame_count);
        
        if isempty(start_index) % in rare cases where a frame is dropped, use next frame to get the lick time
            for count = 1:10
                start_index = find(reward_delivery_time(n)+count == behaviour.wheel_eye_frame_count);
                if ~isempty(start_index)
                    break
                end
            end
        end
        behaviour.reward_delivery_time(n,1) = behaviour.wheel_time(start_index(1));
    end
end

%% Track 1 count
behaviour.track1_count = table2array(reward(:,"Var3"));

%% Track 2 count
behaviour.track2_count = table2array(reward(:,"Var4"));

%% Track 1 performance
behaviour.track1_performance = table2array(reward(:,"Var5"));

%% Track 2 performance
behaviour.track2_performance = table2array(reward(:,"Var6"));

%% reward delivery position
behaviour.reward_position = table2cell(reward(:,"Var7"));
behaviour.reward_position = cellfun(@(x) x(1:end-1), behaviour.reward_position, 'UniformOutput',false);
behaviour.reward_position = cellfun(@str2double,behaviour.reward_position);

%% LickState
behaviour.lick_state = table2cell(lick_performance(:,"Var1"));
behaviour.lick_state = cellfun(@(x) x(2:end), behaviour.lick_state, 'UniformOutput',false);
behaviour.lick_state = cellfun(@str2double,behaviour.lick_state);
if behaviour.lick_state(1) == 0
    behaviour.lick_state(1) = [];
    lick_performance(1,:) = [];
end

%% Onset Eye Frame Count of Lick

if lick_performance_path.datenum > datenum('14-June-2023 13:00:00')
    behaviour.lick_eye_time =  table2array(lick_performance(:,"Var2"));
    for n = 1:length(behaviour.lick_eye_time)
        start_index = find(behaviour.lick_eye_time(n) == behaviour.wheel_arduino_read_time);

        if isempty(start_index) % in rare cases where a frame is dropped, use next frame to get the lick time
            for count = 1:10
                start_index = find(behaviour.lick_eye_time(n)+count == behaviour.wheel_arduino_read_time);
                if ~isempty(start_index)
                    break
                end
            end
        end
        behaviour.lick_time(n,1) = behaviour.wheel_time(start_index(1));
    end
else
    behaviour.lick_eye_frame_count =  table2array(lick_performance(:,"Var2"));
    
    for n = 1:length(behaviour.lick_eye_frame_count)
        start_index = find(behaviour.lick_eye_frame_count(n) == behaviour.wheel_eye_frame_count);

        if isempty(start_index) % in rare cases where a frame is dropped, use next frame to get the lick time
            for count = 1:10
                start_index = find(behaviour.lick_eye_frame_count(n)+count == behaviour.wheel_eye_frame_count);
                if ~isempty(start_index)
                    break
                end
            end
        end
        behaviour.lick_time(n,1) = behaviour.wheel_time(start_index(1));
    end
end

%% Track ID respective to lick
behaviour.lick_track_ID = table2array(lick_performance(:,"Var3"));

%% lick track 1 lap count
behaviour.lick_track1_count = table2array(lick_performance(:,"Var4"));

%% lick track 2 lap count
behaviour.lick_track2_count = table2array(lick_performance(:,"Var5"));

%% reward delivery position
behaviour.lick_position = table2cell(lick_performance(:,"Var6"));
behaviour.lick_position = cellfun(@(x) x(1:end-1), behaviour.lick_position, 'UniformOutput',false);
behaviour.lick_position = cellfun(@str2double,behaviour.lick_position);


if  ~isempty(trial_path)
% Trial info is added 
    behaviour.track_ID_all = table2cell(trial_info(:,"Var1"));
    behaviour.track_ID_all = cellfun(@(x) x(2:end), behaviour.track_ID_all, 'UniformOutput',false);
    behaviour.track_ID_all = cellfun(@str2double,behaviour.track_ID_all);

    behaviour.lap_ID_all = zeros(length(behaviour.track_ID_all),1);
    behaviour.lap_ID_all(find(behaviour.track_ID_all == 1)) = 1:sum(behaviour.track_ID_all == 1);
    behaviour.lap_ID_all(find(behaviour.track_ID_all == 2)) = 1:sum(behaviour.track_ID_all == 2);

    if trial_path.datenum < datenum('10-June-2023 05:00:00')
        %% Start Time eyeframe count
        behaviour.start_time_eye_frame_count = table2array(trial_info(:,"Var2"));
        behaviour.start_time_eye_frame_count(behaviour.track_ID_all==0) = [];

        for n = 1:length(behaviour.start_time_eye_frame_count)
            start_index = find(behaviour.start_time_eye_frame_count(n) == behaviour.wheel_eye_frame_count);

            if isempty(start_index) % in rare cases where a frame is dropped, use next frame to get the lick time
                for count = 1:10
                    start_index = find(behaviour.start_time_eye_frame_count(n)+count == behaviour.wheel_eye_frame_count);
                    if ~isempty(start_index)
                        break
                    end
                end
            end
            behaviour.start_time_all(n,1) = behaviour.wheel_time(start_index(1));
        end

    elseif trial_path.datenum > datenum('10-June-2023 05:00:00')
        %% Start Time based on wheel_arduino_read_time
        behaviour.start_time = table2array(trial_info(:,"Var2"));
        behaviour.start_time(behaviour.track_ID_all==0) = [];

        for n = 1:length(behaviour.start_time)
            start_index = find(behaviour.start_time(n) == behaviour.wheel_arduino_read_time);

            if isempty(start_index) % in rare cases where a frame is dropped, use next frame to get the lick time
                for count = 1:10
                    start_index = find(behaviour.start_time(n)+count == behaviour.wheel_arduino_read_time);
                    if ~isempty(start_index)
                        break
                    end
                end
            end
            behaviour.start_time_all(n,1) = behaviour.wheel_time(start_index(1));
        end

    else
        behaviour.start_time_all = table2array(trial_info(:,"Var2"));
    end


    behaviour.lap_ID_all(behaviour.track_ID_all==0) = [];
    behaviour.track_ID_all(behaviour.track_ID_all==0) = [];
    
else
    behaviour.start_time_all = [];
end

end