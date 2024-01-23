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
function   [behaviour,task_info,peripherals] = import_and_align_Masa_VR_Bonsai(StimulusName,options)

% Program to import bonsai logs and eye data and align to ePhys data and
% align the data based on the delay between quad and photodioide signal
% Everything is resampled to 60Hz

% Find Bonsai files
bonsai_files_names = options(1).bonsai_files_names;
% bonsai_files_names = {bonsai_files.name};

photodiode_path = bonsai_files_names(contains(bonsai_files_names,'PDAsync'));
peripheral_path = bonsai_files_names(contains(bonsai_files_names,'WheelLog'));
trial_path = bonsai_files_names(contains(bonsai_files_names,'Trial_info'));
reward_path = bonsai_files_names(contains(bonsai_files_names,'Reward'));
lick_performance_path = bonsai_files_names(contains(bonsai_files_names,'Lick_Performance'));

% DLC Eyedata pupil size
DLC_EYEDATA_DATAPATH = char(fullfile(options.BONSAI_DATAPATH,options.bonsai_files_names(contains(options.bonsai_files_names,'DLC'))));% new eye log
FACEDATA_DATAPATH = char(fullfile(options.BONSAI_DATAPATH,options.bonsai_files_names(contains(options.bonsai_files_names,'proc'))));% face data;

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

if ~isempty(DLC_EYEDATA_DATAPATH)
    eye_data_raw = readmatrix(DLC_EYEDATA_DATAPATH);
    [pupil_ellipse,new_tracking_points_coordinates] = eye_data_conversion(eye_data_raw);
end

if ~isempty(FACEDATA_DATAPATH)
    facedata = load(FACEDATA_DATAPATH); % load face energy data
end

%%%%%%
%%%%%%
% For now, many reduntant variables are saved. Will remove soon
% Only resampled at 60Hz would be saved
behaviour = [];

%% wheel frame
% behaviour.wheel_frame_count = resample(peripherals.FrameNumber,peripherals.corrected_sglxTime,60)';
behaviour.wheel_frame_count = interp1(peripherals.corrected_sglxTime,peripherals.FrameNumber,peripherals.corrected_sglxTime(1):1/60:peripherals.corrected_sglxTime(end),'previous');
%% wheel Time
behaviour.wheel_time = interp1(peripherals.corrected_sglxTime,peripherals.FrameTime,peripherals.corrected_sglxTime(1):1/60:peripherals.corrected_sglxTime(end),'linear');

%% wheel computer time
% This is the total milisecond of the day from Bonsai. This info is also
% saved in other meta file such as trial info which can be used to find when lick happens and when each lap starts  
behaviour.computer_timevec = interp1(peripherals.corrected_sglxTime,peripherals.Time,peripherals.corrected_sglxTime(1):1/60:peripherals.corrected_sglxTime(end),'linear');

%% wheel spikeGLX time (Use this for ephys)
behaviour.sglxTime = interp1(peripherals.corrected_sglxTime,peripherals.corrected_sglxTime,peripherals.corrected_sglxTime(1):1/60:peripherals.corrected_sglxTime(end),'linear');

behaviour.sglxTime_uncorrected = interp1(peripherals.corrected_sglxTime,peripherals.sglxTime,peripherals.corrected_sglxTime(1):1/60:peripherals.corrected_sglxTime(end),'linear');
behaviour.tvec = behaviour.sglxTime;% standard timevec used for all behaviour related analysis

%% wheel raw input
behaviour.wheel_raw_input = interp1(peripherals.corrected_sglxTime,peripherals.Wheel,peripherals.corrected_sglxTime(1):1/60:peripherals.corrected_sglxTime(end),'linear');

% Not needed
% %% wheel left lick
% behaviour.wheel_lick_L = peripherals.LickL;
% 
% %% wheel right lick
% behaviour.wheel_lick_R = peripherals.LickR;

% %% wheel quad state
% behaviour.wheel_quad_state = peripherals.QuadState;

%% eye data
behaviour.camera_frame_count = interp1(peripherals.sglxTime,peripherals.EyeFrameCount,peripherals.sglxTime(1):1/60:peripherals.sglxTime(end),'previous');
behaviour.eye_coordinates  = new_tracking_points_coordinates(behaviour.camera_frame_count+1,:);% bonsai eye camera frame and DLC frame is 0 based hence adding 1
behaviour.pupil_size =  pupil_ellipse(behaviour.camera_frame_count+1,6)'; % pupil ellipse area 
behaviour.pupil_movement_angle =  pupil_ellipse(behaviour.camera_frame_count+1,8)';% pupil movement direction (in degree) from centre of the previous fitted ellipse to the centre of current fitted ellipse
behaviour.pupil_movement_distance =  pupil_ellipse(behaviour.camera_frame_count+1,7)'; % pupil movement distance from centre of the previous fitted ellipse to the centre of current fitted ellipse 
behaviour.pupil_ellipse =  pupil_ellipse(behaviour.camera_frame_count+1,1:5)'; % good for plotting pupil ellipse



%% facemap - face motion energy and SVD
behaviour.face_motion_enegy = facedata.motion_1(behaviour.camera_frame_count+1); % total motion energt
behaviour.face_motion_SVD = facedata.motSVD_1(behaviour.camera_frame_count+1,1:100); % 1st 100 SVD face energy variable
behaviour.face_motion_mask = facedata.motMask_1(:,1:100); % Mask (spatial component) 1st 100 

%% wheel timestamp of serial string read (time of day, total millisecond)
% This is subject to delete as well, because in the future we will have one
% timestamp (probbaly arduino logged timestamp) as the reference behaviour timestamp
% (which is then aligned to spikeglx time and then corrected based on photodiode)
% behaviour.wheel_arduino_read_time = resample(peripherals.Time,peripherals.corrected_sglxTime,60);

%% 'Real' wheel position and real speed
behaviour.wheel_position = interp1(peripherals.corrected_sglxTime,peripherals.Position,peripherals.corrected_sglxTime(1):1/60:peripherals.corrected_sglxTime(end),'linear');

% Real speed based on wheel raw input
% For actual speed related analysis, This should be used.

tick_to_cm_conversion = 3.1415*20/1024; % radius 20 cm and one circle is 1024 ticks;
speed = [0 diff(behaviour.wheel_raw_input*tick_to_cm_conversion)];
speed(speed<-100) = 0;% big change in speed often due to teleportation or wheel tick resetting
speed(speed>100) = 0;
behaviour.speed = speed./([0 diff(behaviour.sglxTime)]);% 

%% Virtual position and virtual speed
% Convert raw wheel wheel position into virtual position (0 - 140cm for both tracks)
% Position for grey screen period outside of running blocks will be nan
% Position at the end of each track (also looks grey screen) is still
% logged

% Add track 1 position
behaviour.position = nan(1,length(behaviour.wheel_position));
behaviour.track_ID = zeros(1,length(behaviour.wheel_position));
behaviour.position(find(behaviour.wheel_position <= -990)) = abs(behaviour.wheel_position(find(behaviour.wheel_position <= -990))+1140);
behaviour.track_ID((find(behaviour.wheel_position <= -990))) = 1;

% Add track 2 position
behaviour.position(find(behaviour.wheel_position >= -141 & behaviour.wheel_position < 10)) ...
    = abs(behaviour.wheel_position(find(behaviour.wheel_position >= -141 & behaviour.wheel_position < 10))+140);
behaviour.track_ID(find(behaviour.wheel_position >= -141 & behaviour.wheel_position < 10)) = 2;

behaviour.position(behaviour.position>140) = 140;

% %% wheel Photodiode
% behaviour.Photodiode = peripherals.Photodiode;

%% task_info
% Unlike behavour, task_info will save all the meta info about the task
% such as lap start time, lick time, reward delivery time in discrete time
% rather than continuous timeseries.

%% Reward type
% 1 = T1 passive, 2 = T2 passive, 10 = T1 active and 20 = T2 active 
task_info = [];
task_info.reward_type = table2cell(reward(:,"Var1"));
task_info.reward_type = cellfun(@(x) x(2:end), task_info.reward_type, 'UniformOutput',false);
task_info.reward_type = cellfun(@str2double,task_info.reward_type);

% If first datapoint is due to initialisation 
if task_info.reward_type(1) == 0
    task_info.reward_type(1) = [];
    reward(1,:) = [];
end

%%  Reward Delivery Time Based on Eye Frame Count, or Time

DIR = dir(fullfile([BONSAI_DATAPATH,'\',char(reward_path)]));
task_info.reward_delivery_time = [];
task_info.reward_delivery_time_original = [];

if DIR.datenum > datenum('14-June-2023 13:00:00') & DIR.datenum < datenum('01-Jan-2024 13:00:00')
    reward_delivery_time = table2array(reward(:,"Var2")); % May want to simplify this in the future
    
    % Initially search for reward delivery time using original spikeglx timestamp (not 60Hz resampled one)
    for n = 1:length(reward_delivery_time)
        start_index = find(reward_delivery_time(n) == peripherals.Time);
        
%         start_index = interp1(behaviour.computer_timevec,behaviour.computer_timevec,reward_delivery_time(n),'nearest') == behaviour.computer_timevec;

        if isempty(start_index) % in rare cases where a frame is dropped, use next frame to get the lick time
            for count = 1:10
                start_index = find(reward_delivery_time(n)+count == peripherals.Time);
                if ~isempty(start_index)
                    break
                end
            end
        end
        % Original reward time
        task_info.reward_delivery_time_original(n,1) = peripherals.corrected_sglxTime(start_index);
        % Reward time basewd on resampled 60Hz timestamp
        task_info.reward_delivery_time(n,1) = interp1(behaviour.sglxTime,behaviour.sglxTime,task_info.reward_delivery_time_original(n,1),'nearest');
    end
elseif DIR.datenum < datenum('14-June-2023 13:00:00')
    % Before 14/06/2023 During training, camera/eye frame count was used to find
    % when reward delivered.
    reward_delivery_camera_frame_count = table2array(reward(:,"Var2"));
    
    for n = 1:length(reward_delivery_camera_frame_count)
        start_index = find(reward_delivery_camera_frame_count(n) == peripherals.EyeFrameCount);
        
        if isempty(start_index) % in rare cases where a frame is dropped, use next frame to get the lick time
            for count = 1:10
                start_index = find(reward_delivery_camera_frame_count(n)+count == peripherals.EyeFrameCount);
                if ~isempty(start_index)
                    break
                end
            end
        end

        % Original reward time
        task_info.reward_delivery_time_original(n,1) = peripherals.corrected_sglxTime(start_index);
        % Reward time basewd on resampled 60Hz timestamp
        task_info.reward_delivery_time(n,1) = interp1(behaviour.sglxTime,behaviour.sglxTime,task_info.reward_delivery_time_original(n,1),'nearest');
    end
end

%% Track 1 count
% For each lap, what is the current cumulative track 1 count
% For example, for 11th lap, it may be 3rd lap of track 2 but with track 1 count
% of 8
task_info.reward_track1_count = table2array(reward(:,"Var3"));

%% Track 2 count
% For each lap, what is the current cumulative track 2 count
task_info.reward_track2_count = table2array(reward(:,"Var4"));


%% Track 1 performance
% Probably not needed
task_info.track1_performance = table2array(reward(:,"Var5"));

%% Track 2 performance
task_info.track2_performance = table2array(reward(:,"Var6"));

%% reward delivery position
task_info.reward_position = table2cell(reward(:,"Var7"));
task_info.reward_position = cellfun(@(x) x(1:end-1), task_info.reward_position, 'UniformOutput',false);
task_info.reward_position = cellfun(@str2double,task_info.reward_position);

%% Lick State
% 1 is left and 2 is right
task_info.lick_state = table2cell(lick_performance(:,"Var1"));
task_info.lick_state = cellfun(@(x) x(2:end), task_info.lick_state, 'UniformOutput',false);
task_info.lick_state = cellfun(@str2double,task_info.lick_state);

if isempty(task_info.lick_state) == 1
    task_info.lick_state = [];
elseif task_info.lick_state(1) == 0
    task_info.lick_state(1) = [];
    lick_performance(1,:) = [];
end

%% Lick time
task_info.lick_time_original = [];
task_info.lick_time = [];

if DIR.datenum > datenum('14-June-2023 13:00:00') & DIR.datenum < datenum('01-Jan-2024 13:00:00')
    lick_time =  table2array(lick_performance(:,"Var2"));

    for n = 1:length(lick_time)
        start_index = find(lick_time(n) == peripherals.Time);

        if isempty(start_index) % in rare cases where a frame is dropped, use next frame to get the lick time
            for count = 1:10
                start_index = find(lick_time(n)+count == peripherals.Time);
                if ~isempty(start_index)
                    break
                end
            end
        end
        
        task_info.lick_time_original(n,1) = peripherals.corrected_sglxTime(start_index);
        % Reward time basewd on resampled 60Hz timestamp
        task_info.lick_time(n,1) = interp1(behaviour.sglxTime,behaviour.sglxTime,task_info.lick_time_original(n,1),'nearest');
    end

elseif DIR.datenum < datenum('14-June-2023 13:00:00')

    % Before 14/06/2023 During training, camera/eye frame count was used to find
    % when licked.
    lick_camera_frame_count = table2array(reward(:,"Var2"));

    for n = 1:length(lick_camera_frame_count)
        start_index = find(lick_camera_frame_count(n) == peripherals.EyeFrameCount);

        if isempty(start_index) % in rare cases where a frame is dropped, use next frame to get the lick time
            for count = 1:10
                start_index = find(lick_camera_frame_count(n)+count == peripherals.EyeFrameCount);
                if ~isempty(start_index)
                    break
                end
            end
        end
        task_info.lick_time_original(n,1) = peripherals.corrected_sglxTime(start_index);
        % Reward time basewd on resampled 60Hz timestamp
        task_info.lick_time(n,1) = interp1(behaviour.sglxTime,behaviour.sglxTime,task_info.reward_delivery_time_original(n,1),'nearest');
    end
end


%% Track ID respective to lick
task_info.lick_track_ID = table2array(lick_performance(:,"Var3"));

%% lick track 1 lap count
% redundant info
task_info.lick_track1_count = table2array(lick_performance(:,"Var4"));

%% lick track 2 lap count
% redundant info
task_info.lick_track2_count = table2array(lick_performance(:,"Var5"));

%% reward delivery position
task_info.lick_position = table2cell(lick_performance(:,"Var6"));
task_info.lick_position = cellfun(@(x) x(1:end-1), task_info.lick_position, 'UniformOutput',false);
task_info.lick_position = cellfun(@str2double,task_info.lick_position);

%% Lick count time series

% behaviour.lick_count = []; 


% Based on task info lick time (only after task started)
% t_edges = behaviour.sglxTime(1)-(behaviour.sglxTime(2) - behaviour.sglxTime(1))/2:1/60:behaviour.sglxTime(end)+(behaviour.sglxTime(end) - behaviour.sglxTime(end-1))/2;
% for lick_id = 1:max(task_info.lick_state) % left = 1 and right = 2
%      [N,edges,bin] = histcounts(task_info.lick_time(find(task_info.lick_state==lick_id)),t_edges);
%      behaviour.lick_count(lick_id,:) = N;
% end

% Based on Peripheral-logged lick
behaviour.lick_count(1,:) = interp1(peripherals.sglxTime,[0; diff(peripherals.LickL)],peripherals.sglxTime(1):1/60:peripherals.sglxTime(end),'previous');
behaviour.lick_count(2,:) = interp1(peripherals.sglxTime,[0; diff(peripherals.LickR)],peripherals.sglxTime(1):1/60:peripherals.sglxTime(end),'previous');

%% Lap info 

if  ~isempty(trial_path) & contains(StimulusName,'RUN')

%     if table2array(trial_info(1,"Var5")) >100
%         trial_info(1,:) =  [];
%     end
    
    % Track ID of each lap 
    task_info.track_ID_all = table2cell(trial_info(:,"Var1"));
    task_info.track_ID_all = cellfun(@(x) x(2:end), task_info.track_ID_all, 'UniformOutput',false);
    task_info.track_ID_all = cellfun(@str2double,task_info.track_ID_all);
    
    % Block ID of each lap
    block_transition = [1; diff(task_info.track_ID_all)];
    block_transition(block_transition~=0) = 1;
    task_info.block_ID_all = cumsum(block_transition);

    % Lap ID (Track-specific lap ID) 
    % E.g. 9th lap may be 1st lap on track 2.
    task_info.lap_ID_all = zeros(length(task_info.track_ID_all),1);
    task_info.lap_ID_all(find(task_info.track_ID_all == 1)) = 1:sum(task_info.track_ID_all == 1);
    task_info.lap_ID_all(find(task_info.track_ID_all == 2)) = 1:sum(task_info.track_ID_all == 2);

    % Trial type for each lap
    if size(trial_info,2) > 5
        task_info.trial_type = table2cell(trial_info(:,"Var6"));
        task_info.trial_type = cellfun(@(x) x(1:1), task_info.trial_type, 'UniformOutput',false);
        task_info.trial_type = cellfun(@str2double,task_info.trial_type); % 1 is active only and 2 is hybrid
    end

    if DIR.datenum > datenum('14-June-2023 13:00:00') & DIR.datenum < datenum('01-Jan-2024 13:00:00')
        %% Start Time eyeframe count
        start_time_all = table2array(trial_info(:,"Var2"));
        for n = 1:length(start_time_all)
            start_index = find(start_time_all(n) == peripherals.Time);

            if isempty(start_index) % in rare cases where a frame is dropped, use next frame to get the lick time
                for count = 1:10
                    start_index = find(start_time_all(n)+count == peripherals.Time);
                    if ~isempty(start_index)
                        break
                    end
                end
            end

            task_info.start_time_all_original(n,1) = peripherals.corrected_sglxTime(start_index);
            % Reward time basewd on resampled 60Hz timestamp
            task_info.start_time_all(n,1) = interp1(behaviour.sglxTime,behaviour.sglxTime,task_info.start_time_all_original(n,1),'nearest');
        end

    else
        task_info.start_time_all = table2array(trial_info(:,"Var2"));
    end

else
    task_info.start_time_all = [];
end


%% Extract laps info
task_info.complete_laps_id = [];
task_info.aborted_laps_id = [];

if ~isempty(task_info.start_time_all)

    start_indices= [];

    x = behaviour.position;  % position data
    t= behaviour.sglxTime; % time

    for nlap = 1:length(task_info.track_ID_all)
        start_indices(nlap) = find(t == task_info.start_time_all(nlap));
    end

    for nlap = 1:length(start_indices)

        if nlap < length(start_indices)
            current_lap_x = x(start_indices(nlap):start_indices(nlap+1));
            current_lap_t = t(start_indices(nlap):start_indices(nlap+1));
        else
            current_lap_x = x(start_indices(nlap):end);
            current_lap_t = t(start_indices(nlap):end);
        end


        if length(current_lap_x) > 1 && length(current_lap_x(~isnan(current_lap_x))) >1% Only if getting more than 1 datapoint (maybe noise)
            on_track_x = current_lap_x(~isnan(current_lap_x));
            on_track_t = current_lap_t(~isnan(current_lap_x));

            if sum(on_track_x==0)>0
                start_position = find(on_track_x == 0);
                on_track_x = on_track_x(start_position(1):end);
                on_track_t = on_track_t(start_position(1):end);
            end

            if on_track_x(end) ~= on_track_x(end-1)
                on_track_x(end) = on_track_x(end-1);
                on_track_t(end) = on_track_t(end-1);
            end

            [last_position last_position_index] = max(on_track_x);
%             task_info.end_time_all(nlap) = on_track_t(last_position_index);

            if last_position >= 139 % sometimes last lap ends before 140cm
                task_info.end_time_all(nlap) = on_track_t(last_position_index); % End time in terms of reaching end of track
                task_info.complete_laps_id = [task_info.complete_laps_id nlap];% Lap id here is all laps (not track-specific id)

            else
                % If never reached the end, then end time is the end of the
                % entire lap (jumped to next lap after reaching time threshold)
                task_info.end_time_all(nlap) = on_track_t(end);
                task_info.aborted_laps_id = [task_info.aborted_laps_id nlap]; % Lap id here is all laps (not track-specific id)
            end
            %                 lap_count = lap_count + 1;
        end
    end
end


if ~isfield(options,'photodiode_failure')
    task_info.pd_on.sglxTime = pdstart';
    task_info.pd_off.sglxTime = pdend';
end
% extract_laps_masa(1,behaviour,position)

end