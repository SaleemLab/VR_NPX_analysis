% function [stimulusData,eyeData,wheelData,photodiodeData,stimTimes] = importAndAlignBonsaiLogs(EPHYS_DATAPATH,TRIALDATA_DATAPATH, PERIPHERALS_DATAPATH,EYEDATA_DATAPATH)
% 
% Program to import bonsai logs and eye data and align to ePhys data
% 
% inputs : EPHYS_DATAPATH, TRIALDATA_DATAPATH, PERIPHERALS_DATAPATH, EYEDATA_DATAPATH
% outputs: stimulusData,eyeData,wheelData,photodiodeData (photdiode data derived from the
%               peripherals dataset)
%
% History:  SGS 24/2/2022 ported from import_process_BonsaiData
%           SGS 28/2/2022 added stimulus data parsing
%           SGS 04/4/2022 shifted imec extract to here and allowed function to either call precomputed
%                           asyncpulse times or load anew from LF (long) or AP (very long) imec data
%           SGS 14/4/2022 moved to using toolbox from Jennifer Collonel
%                       (SpikeGLX_Datafile_Tools) to read the async pulse from the ephys file
%           FR 09/06/2022 changed how PD was being read due to stimuli
%           adjustment issues
% Dependencies:
%   import_bonsai_peripherals
%   import_bonsai_eyedata
%   import_BonVisionParams
%   import_processWheelEphys
%   alignBonsaiAndEphysData
%
%   Uses the following libraries, which are included in the folder
%   https://billkarsh.github.io/SpikeGLX/ - SpikeGLX_Datafile_Tools
%   
function [behaviour,position] = import_and_align_Masa_VR_Bonsai(StimulusName,options)

% Program to import bonsai logs and eye data and align to ePhys data and
% align the data based on the delay between quad and photodioide signal

% Find Bonsai files
bonsai_files_names = options(1).bonsai_files_names;
% bonsai_files_names = {bonsai_files.name};

photodiode_path = bonsai_files_names(contains(bonsai_files_names,'PDAsync'));
peripheral_path = bonsai_files_names(contains(bonsai_files_names,'WheelLog'));
% eyedata_path = bonsai_files_names(contains(bonsai_files_names,'EyeLog'));
trial_path = bonsai_files_names(contains(bonsai_files_names,'Trial_info'));
reward_path = bonsai_files_names(contains(bonsai_files_names,'Reward'));
lick_performance_path = bonsai_files_names(contains(bonsai_files_names,'Lick_Performance'));
BONSAI_DATAPATH = options(1).BONSAI_DATAPATH;


DIR = dir(fullfile(options(1).EPHYS_DATAPATH,'*syncpulseTimes.mat'))
if ~isempty(DIR) % everything sync to probe 1
    load(fullfile(options(1).EPHYS_DATAPATH,DIR.name))
    syncTimes_ephys = syncTimes_ephys.on; % use upswings currently
else
    [AP_FILE,LF_FILE] = findImecBinFile(options(1).EPHYS_DATAPATH);
    %     FILE_TO_USE = AP_FILE;
    binpath = fullfile(options(1).EPHYS_DATAPATH,AP_FILE);
    %     binpath = fullfile(options(1).EPHYS_DATAPATH,LF_FILE);
    syncTimes_ephys = SGLXextractSyncPulseFromBinFile(binpath);
    parseDate = date;
    [~,fname] = fileparts(options(1).EPHYS_DATAPATH);
    save(fullfile(options(1).EPHYS_DATAPATH,[fname,'_syncpulseTimes.mat']),'syncTimes_ephys','parseDate');
    syncTimes_ephys = syncTimes_ephys.on; % use upswings currently
end


% Import wheel data, eye data, position data and photodiodide data 
% and synchronise them to async pulse 
peripherals = import_bonsai_peripherals(fullfile([BONSAI_DATAPATH,'\',char(peripheral_path)]));
[peripherals] = alignBonsaiToEphysSyncTimes(peripherals,syncTimes_ephys);

% eyeData = import_bonsai_eyedata(fullfile([BONSAI_DATAPATH,'\',char(EYEDATA_DATAPATH)]));
% [eyeData] = alignBonsaiToEphysSyncTimes(eyeData,syncTimes_ephys);
eyeData = [];

if exist('photodiode_path','var') && ~isempty(photodiode_path)
    [photodiode, photodiode_sync] = import_bonsai_photodiode(fullfile([BONSAI_DATAPATH,'\',char(photodiode_path)]));
    photodiode.Sync = photodiode.S;
    photodiode.Time = photodiode.T;
    photodiode.Photodiode = photodiode.P;
    photodiode = rmfield(photodiode,{'T','S','P'});
    [photodiode] = alignBonsaiToEphysSyncTimes(photodiode,syncTimes_ephys);

else
    photodiode = []; photodiode_sync = [];
end

reward = readtable(fullfile([BONSAI_DATAPATH,'\',char(reward_path)]));
lick_performance = readtable(fullfile([BONSAI_DATAPATH,'\',char(lick_performance_path)]));


if ~isempty(trial_path)
    trial_info = readtable(fullfile([BONSAI_DATAPATH,'\',char(trial_path)]));
end


%% Correct sglxtime based on photodiode and quadstate

photodiode.Photodiode_smoothed = smoothdata(photodiode.Photodiode,'movmedian',10);
pd_ON_OFF = photodiode.Photodiode_smoothed >= 80; % find ON and OFF
pd_ON = find(diff(pd_ON_OFF) == 1)+1; % Find ON
pd_OFF = find(diff(pd_ON_OFF) == -1)+1; % Find OFF

if length(pd_ON) == length(pd_OFF)
    block_length = abs(pd_ON-pd_OFF);
elseif length(pd_ON)>length(pd_OFF)
    for i = 1:length(pd_OFF)
        block_length(i) = abs(pd_ON(i) - pd_OFF(i));
    end
else
    for i = 1:length(pd_ON)
        block_length(i) = abs(pd_ON(i) - pd_OFF(i));
    end
end
blocks_ind = find(block_length>80);% sample rate of photodiode is 3000 per second, 50 samples per frame (60 frame per second), 
% photodiode.binary_pd = zeros(length(photodiode.Photodiode),1);
% for nblock = 1:length(blocks_ind)
%     photodiode.binary_pd(pd_ON(blocks_ind(nblock)):pd_OFF(blocks_ind(nblock))) = 1;
% end

idx_trial_start = find(diff(peripherals.QuadState)==1) + 1;
idx_trial_end = find(diff(peripherals.QuadState)==-1) + 1;
frame_trial_start = peripherals.FrameNumber(idx_trial_start);
frame_trial_end = peripherals.FrameNumber(idx_trial_end);
temp_tbl.start_time = peripherals.sglxTime(idx_trial_start);
temp_tbl.end_time = peripherals.sglxTime(idx_trial_end);
nTrials = length(temp_tbl.start_time);

% peripherals.sglxTime(idx_trial_start);
% peripherals.sglxTime(idx_trial_end);
% photodiode.sglxTime(pd_ON(blocks_ind));
% photodiode.sglxTime(pd_OFF(blocks_ind));

% if difference between photodiode and quadstate is less than 10%, it is probably more or less acruate
if abs(length(blocks_ind)-length(idx_trial_start))/length(idx_trial_start) <    0.1

    disp('Aligning bonsai timestamp based on the delay between photodiode and quad state (may take 10 mins)')
    pdstart = [];
    tic
    % photodiode_index = 1:length(photodiode.FrameNumber);
    % blocks_index = 1:length(blocks_ind);
    H = waitbar(0,'Finding photodiode ON edges');

    for itrial =1:length(idx_trial_start)
        waitbar(itrial/length(idx_trial_start),H)
        % find frame where Bonsai asked to change quad
        frameidx_start = find(photodiode.FrameNumber <= frame_trial_start(itrial),1,'last');

        % find next index where the photodiode detected a quad change
        temp_idx = find(pd_ON(blocks_ind) > frameidx_start,1,'first');
        idx_start = pd_ON(blocks_ind(temp_idx))-1;
        %
        %     % find frame where Bonsai asked to change quad
        %     frameidx_start = photodiode_index(photodiode.FrameNumber <= frame_trial_start(itrial));
        %     frameidx_start = frameidx_start(end);
        %
        %     % find next index where the photodiode detected a quad change
        %     temp_idx = blocks_index(pd_ON(blocks_ind) > frameidx_start);
        %
        %     if ~isempty(temp_idx)
        %         temp_idx = temp_idx(1);
        %     end
        %     idx_start = pd_ON(blocks_ind(temp_idx))-1;

        % convert to spike GLX time
        if ~isempty(idx_start)
            pdstart(itrial) = photodiode.sglxTime(idx_start);
        else
            pdstart(itrial) = nan;
        end
    end
    toc

    pdend = []
    tic
    H = waitbar(0,'Finding photodiode OFF edges');
    for itrial =1:length(idx_trial_end)
        waitbar(itrial/length(idx_trial_end),H)

        % find frame where Bonsai asked to change quad
        frameidx_end = find(photodiode.FrameNumber <= frame_trial_end(itrial),1,'last');

        % find next index where the photodiode detected a quad change
        temp_idx = find(pd_OFF(blocks_ind) > frameidx_end,1,'first');
        idx_end = pd_OFF(blocks_ind(temp_idx))-1;

        %     % find frame where Bonsai asked to change quad
        %     frameidx_end = photodiode_index(photodiode.FrameNumber <= frame_trial_end(itrial));
        %     frameidx_end = frameidx_end(end);
        %
        %     % find next index where the photodiode detected a quad change
        %     temp_idx = blocks_index(pd_ON(blocks_ind) > frameidx_end);
        %
        %     if ~isempty(temp_idx)
        %         temp_idx = temp_idx(1);
        %     end
        %
        %     idx_end = pd_ON(blocks_ind(temp_idx))-1;

        % convert to spike GLX time
        if ~isempty(idx_end)
            pdend(itrial) = photodiode.sglxTime(idx_end);
        else
            pdend(itrial) = nan;
        end
    end
    toc

    photodiodeData.stim_on.sglxTime = pdstart';
    photodiodeData.stim_off.sglxTime = pdend';
else
    options(1).photodiode_failure = 1;
end


if contains(StimulusName,'Masa2tracks')
    StimulusType = 'Masa2tracks';
else
     StimulusType = 'Track';
end

if isfield(options,'photodiode_failure') 
     peripherals.corrected_sglxTime = peripherals.sglxTime;
     disp('photodiode signal failed. Bonsai data not corrected')
else
    switch StimulusType
        case 'replay_Masa2tracks'
            [peripherals] = alignBonsaiToPhotodiode(peripherals,sort([photodiodeData.stim_on.sglxTime; photodiodeData.stim_off.sglxTime]),'replay');
        case 'Masa2tracks' % For now, just subtract from the average delay (the variance is quite small (e.g. 0.5s to 0.6s))
            [peripherals] = alignBonsaiToPhotodiode(peripherals,sort([photodiodeData.stim_on.sglxTime; photodiodeData.stim_off.sglxTime]),[]);
        case 'Track'
            [peripherals] = alignBonsaiToPhotodiode(peripherals,sort([photodiodeData.stim_on.sglxTime; photodiodeData.stim_off.sglxTime]),[]);
    end
end

% figure;
% plot(photodiode.sglxTime(1:20000),photodiode.Photodiode_smoothed(1:20000))
% hold on
% plot(peripherals.corrected_sglxTime(1:300),peripherals.QuadState(1:300)*100)
% plot(photodiode.sglxTime(1:20000),photodiode.binary_pd(1:20000)*120)


% raw_LFP(nchannel,:) = interp1(peripherals.sglxTime,raw_LFP(nchannel,:),probe_1_tvec(1:2:end),'linear'); %probe 1 LFP time to match probe 2 LFP time


%% wheel frame
behaviour.wheel_frame = peripherals.FrameNumber;

%% wheel Time
behaviour.wheel_time = peripherals.FrameTime;

%% wheel computer time
behaviour.wheel_computer_time = peripherals.FrameComputerDateVec;

%% wheel spikeGLX time (Use this for ephys)
behaviour.sglxTime = peripherals.corrected_sglxTime;
behaviour.original_sglxTime = peripherals.sglxTime;

%% wheel raw input
behaviour.wheel_raw_input = peripherals.Wheel;

%% wheel left lick
behaviour.wheel_lick_L = peripherals.LickL;

%% wheel right lick
behaviour.wheel_lick_R = peripherals.LickR;

%% wheel quad state
behaviour.wheel_quad_state = peripherals.QuadState;

%% wheel camera frame count
behaviour.wheel_eye_frame_count = peripherals.EyeFrameCount;
behaviour.wheel_reference_timestamp = peripherals.Time;

%% wheel timestamp of serial string read (time of day, total millisecond)
behaviour.wheel_arduino_read_time = peripherals.Time;

%% wheel position and speed
behaviour.wheel_position = peripherals.Position;

% %% wheel Photodiode
% behaviour.Photodiode = peripherals.Photodiode;

% Calculate wheel speed at each time point
tick_to_cm_conversion = 3.1415*20/1024; % radius 20 cm and one circle is 1024 ticks;
speed = [0; diff(peripherals.Wheel*tick_to_cm_conversion)];
speed(speed<-100) = 0;
speed(speed>100) = 0;
behaviour.speed = speed/mean(diff(behaviour.sglxTime));


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

DIR = dir(fullfile([BONSAI_DATAPATH,'\',char(reward_path)]))

if DIR.datenum > datenum('14-June-2023 13:00:00')
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
        
        behaviour.reward_delivery_time(n,1) = behaviour.sglxTime(start_index);
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
        behaviour.reward_delivery_time(n,1) = behaviour.sglxTime(start_index(1));
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

if DIR.datenum > datenum('14-June-2023 13:00:00')
    behaviour.lick_eye_time =  table2array(lick_performance(:,"Var2"));

    for n = 1:length(behaviour.lick_eye_time)
        start_index = find(behaviour.lick_eye_time(n) == behaviour.wheel_reference_timestamp);

        if isempty(start_index) % in rare cases where a frame is dropped, use next frame to get the lick time
            for count = 1:10
                start_index = find(behaviour.lick_eye_time(n)+count == behaviour.wheel_reference_timestamp);
                if ~isempty(start_index)
                    break
                end
            end
        end

        behaviour.lick_time(n,1) = behaviour.sglxTime(start_index);
    end

else
    behaviour.lick_eye_frame_count =  table2array(lick_performance(:,"Var2"));
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


if  ~isempty(trial_path) & contains(StimulusName,'RUN')

%     if table2array(trial_info(1,"Var5")) >100
%         trial_info(1,:) =  [];
%     end

    behaviour.track_ID_all = table2cell(trial_info(:,"Var1"));
    behaviour.track_ID_all = cellfun(@(x) x(2:end), behaviour.track_ID_all, 'UniformOutput',false);
    behaviour.track_ID_all = cellfun(@str2double,behaviour.track_ID_all);

    behaviour.lap_ID_all = zeros(length(behaviour.track_ID_all),1);
    behaviour.lap_ID_all(find(behaviour.track_ID_all == 1)) = 1:sum(behaviour.track_ID_all == 1);
    behaviour.lap_ID_all(find(behaviour.track_ID_all == 2)) = 1:sum(behaviour.track_ID_all == 2);
    if size(trial_info,2) > 5
        behaviour.trial_type = table2cell(trial_info(:,"Var6"));
        behaviour.trial_type = cellfun(@(x) x(1:1), behaviour.trial_type, 'UniformOutput',false);
        behaviour.trial_type = cellfun(@str2double,behaviour.trial_type); % 1 is active only and 2 is hybrid
    end
    if DIR.datenum > datenum('14-June-2023 13:00:00')
        %% Start Time eyeframe count
        behaviour.start_time_all = table2array(trial_info(:,"Var2"));
        for n = 1:length(behaviour.start_time_all)
            start_index = find(behaviour.start_time_all(n) == behaviour.wheel_reference_timestamp);

            if isempty(start_index) % in rare cases where a frame is dropped, use next frame to get the lick time
                for count = 1:10
                    start_index = find(behaviour.start_time_all(n)+count == behaviour.wheel_reference_timestamp);
                    if ~isempty(start_index)
                        break
                    end
                end
            end
            behaviour.start_time_all(n,1) = behaviour.sglxTime(start_index(1));
        end

    else
        behaviour.start_time_all = table2array(trial_info(:,"Var2"));
    end
else
    behaviour.start_time_all = [];
end

if ~isfield(options,'photodiode_failure')
    behaviour.pd_on.sglxTime = pdstart';
    behaviour.pd_off.sglxTime = pdend';
end

%% Convert behaviour data into position data for extracting laps

data.linear(1).linear = nan(1,length(behaviour.wheel_position));
data.linear(1).linear(find(behaviour.wheel_position <= -990)) = abs(behaviour.wheel_position(find(behaviour.wheel_position <= -990))+1140);
data.linear(1).linear(data.linear(1).linear>140) = 140;


behaviour.track1_position = data.linear(1).linear;
data.x = nan(1,length(behaviour.wheel_position));
data.x(~isnan(data.linear(1).linear)) = data.linear(1).linear(~isnan(data.linear(1).linear));

if sum(find(behaviour.wheel_position >= -141 & behaviour.wheel_position < 0)) > 0 % If there is track 2
    data.linear(2).linear = nan(1,length(behaviour.wheel_position));
    data.linear(2).linear(find(behaviour.wheel_position >= -141 & behaviour.wheel_position < 10)) ...
        = abs(behaviour.wheel_position(find(behaviour.wheel_position >= -141 & behaviour.wheel_position < 10))+140);
    data.linear(2).linear(data.linear(2).linear>140) = 140;
    behaviour.track2_position = data.linear(2).linear;
    data.x(~isnan(data.linear(2).linear)) = data.linear(2).linear(~isnan(data.linear(2).linear));
end

data.t = behaviour.sglxTime';
data.v = behaviour.speed';
% data.v(2:end) = diff(data.x);
data.v_cm = data.v;
position = data;

end