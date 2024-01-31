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
function   [behaviour,task_info,peripherals] = import_and_align_visual_stimuli_Bonsai(StimulusName,options)


EPHYS_DATAPATH = options.EPHYS_DATAPATH;

if isfield(options,'bonsai_files_names')
    PERIPHERALS_DATAPATH = char(fullfile(options.BONSAI_DATAPATH,options.bonsai_files_names(contains(options.bonsai_files_names,'WheelLog'))));
    TRIALDATA_DATAPATH = char(fullfile(options.BONSAI_DATAPATH,options.bonsai_files_names(contains(options.bonsai_files_names,'Trial'))));
    PHOTODIODE_DATAPATH = char(fullfile(options.BONSAI_DATAPATH,options.bonsai_files_names(contains(options.bonsai_files_names,'PDA'))));
    EYEDATA_DATAPATH = char(fullfile(options.BONSAI_DATAPATH,options.bonsai_files_names(contains(options.bonsai_files_names,'EyeLog'))));% Old eye log
    DLC_EYEDATA_DATAPATH = char(fullfile(options.BONSAI_DATAPATH,options.bonsai_files_names(contains(options.bonsai_files_names,'DLC'))));% new eye log
    FACEDATA_DATAPATH = char(fullfile(options.BONSAI_DATAPATH,options.bonsai_files_names(contains(options.bonsai_files_names,'proc'))));% face data;
else
    PERIPHERALS_DATAPATH = options.PERIPHERALS_DATAPATH;
    EYEDATA_DATAPATH = options.EYEDATA_DATAPATH;
    TRIALDATA_DATAPATH = options.TRIALDATA_DATAPATH;
end

if contains(options.bonsai_files_names{1,1},'SparseNoise_fullscreen')
    TRIALDATA_DATAPATH = char(fullfile(options.BONSAI_DATAPATH,options.bonsai_files_names(contains(options.bonsai_files_names,'.bin'))));
end


if contains(options.StimulusName,'Checkerboard')
    options.paradigm = 'photodiode_timestamp';
else
    options.paradigm = 'masa';
end

%%%%%
% Set defaults
if ~isfield(options,'verboseFlag') || isempty(options.verboseFlag)
    verboseFlag = 1;
else
    verboseFlag = options.verboseFlag;
end

if ~isfield(options,'otherMetrics') || isempty(options.otherMetrics)
    otherMetrics = {'wheelspeed','pupilarea'};
else
    otherMetrics = options.otherMetrics;
end

%%%%%
% Step 1: Extract relevant data from wheel, pupil and stimulus files and
% sync different timebases
if verboseFlag
    fprintf('\n\tImporting stimulus, peripherals, eye and photodiode data from \n\t\t%s and \n\t\t%s and \n\t\t%s and \n\t\t%s...',TRIALDATA_DATAPATH,PERIPHERALS_DATAPATH,EYEDATA_DATAPATH,PHOTODIODE_DATAPATH)
end
tstart = tic;

[stimData,eyeData,peripherals,photodiodeData,stimTimes] = import_and_align_bonsai_logs(EPHYS_DATAPATH,TRIALDATA_DATAPATH,PERIPHERALS_DATAPATH,EYEDATA_DATAPATH,PHOTODIODE_DATAPATH,options);
% [stimData,eyeData,peripherals,photodiodeData,stimTimes] = import_and_align_bonsai_logs(EPHYS_DATAPATH,TRIALDATA_DATAPATH,PERIPHERALS_DATAPATH,EYEDATA_DATAPATH,[],options);


if ~isempty(DLC_EYEDATA_DATAPATH)
    eye_data_raw = readmatrix(DLC_EYEDATA_DATAPATH);
    [pupil_ellipse,new_tracking_points_coordinates] = eye_data_conversion(eye_data_raw);
end

if ~isempty(FACEDATA_DATAPATH)
    facedata = load(FACEDATA_DATAPATH); % load face energy data
    if isfield(facedata,'motion_1')
        facedata = [];
    end
end

basicDataLoadTime = toc(tstart);

if verboseFlag
    fprintf('\n\t\tTook %3.1fs to load Bonsai wheel and eye data',basicDataLoadTime)
end

%%%%%%
%%%%%%
% For now, many reduntant variables are saved. Will remove soon
% Only resampled at 60Hz would be saved
behaviour = [];

%% wheel frame
% behaviour.wheel_frame_count = resample(peripherals.FrameNumber,peripherals.sglxTime,60)';
behaviour.wheel_frame_count = interp1(peripherals.sglxTime,peripherals.FrameNumber,peripherals.sglxTime(1):1/60:peripherals.sglxTime(end),'previous');
%% wheel Time
behaviour.wheel_time = interp1(peripherals.sglxTime,peripherals.FrameTime,peripherals.sglxTime(1):1/60:peripherals.sglxTime(end),'linear');

%% wheel computer time
% This is the total milisecond of the day from Bonsai. This info is also
% saved in other meta file such as trial info which can be used to find when lick happens and when each lap starts  
behaviour.computer_timevec = interp1(peripherals.sglxTime,peripherals.Time,peripherals.sglxTime(1):1/60:peripherals.sglxTime(end),'linear');

%% wheel spikeGLX time (Use this for ephys)
% behaviour.sglxTime = resample(peripherals.sglxTime,peripherals.sglxTime,60)';
behaviour.sglxTime = interp1(peripherals.sglxTime,peripherals.sglxTime,peripherals.sglxTime(1):1/60:peripherals.sglxTime(end),'linear');

% behaviour.sglxTime_uncorrected = resample(peripherals.sglxTime,peripherals.corrected_sglxTime,60)';
behaviour.tvec = behaviour.sglxTime;% standard timevec used for all behaviour related analysis

%% eye camera frame count
% behaviour.camera_frame_count = resample(peripherals.EyeFrameCount,peripherals.sglxTime,60)';
behaviour.camera_frame_count = interp1(peripherals.sglxTime,peripherals.EyeFrameCount,peripherals.sglxTime(1):1/60:peripherals.sglxTime(end),'previous');
% behaviour.wheel_reference_timestamp = peripherals.Time;% probably not needed....

%% eye data
if ~isempty(DLC_EYEDATA_DATAPATH)
    behaviour.camera_frame_count = interp1(peripherals.sglxTime,peripherals.EyeFrameCount,peripherals.sglxTime(1):1/60:peripherals.sglxTime(end),'previous');
    behaviour.eye_coordinates  = new_tracking_points_coordinates(behaviour.camera_frame_count+1,:);% bonsai eye camera frame and DLC frame is 0 based hence adding 1
    behaviour.pupil_size =  pupil_ellipse(behaviour.camera_frame_count+1,6)';
    behaviour.pupil_movement_angle =  pupil_ellipse(behaviour.camera_frame_count+1,8)';
    behaviour.pupil_movement_distance =  pupil_ellipse(behaviour.camera_frame_count+1,7)';
end

%% facemap - face motion energy and SVD
if ~isempty(FACEDATA_DATAPATH)
    behaviour.face_motion_enegy = facedata.motion_1(behaviour.camera_frame_count+1); % total motion energt
    behaviour.face_motion_SVD = facedata.motSVD_1(behaviour.camera_frame_count+1,1:100); % 1st 100 SVD face energy variable (temporal component)
    behaviour.face_motion_mask = facedata.motMask_1(:,1:100); % Mask (spatial component) 1st 100
end
%% wheel raw input
behaviour.wheel_raw_input = interp1(peripherals.sglxTime,peripherals.Wheel,peripherals.sglxTime(1):1/60:peripherals.sglxTime(end),'linear');

% Not needed
% %% wheel left lick
% behaviour.wheel_lick_L = peripherals.LickL;
% 
% %% wheel right lick
% behaviour.wheel_lick_R = peripherals.LickR;

% %% wheel quad state
% behaviour.wheel_quad_state = peripherals.QuadState;

%% Peripheral-logged lick count
behaviour.lick_count = [];
behaviour.lick_count(1,:) = interp1(peripherals.sglxTime,[0; diff(peripherals.LickL)],peripherals.sglxTime(1):1/60:peripherals.sglxTime(end),'previous');
behaviour.lick_count(2,:) = interp1(peripherals.sglxTime,[0; diff(peripherals.LickR)],peripherals.sglxTime(1):1/60:peripherals.sglxTime(end),'previous');

%% wheel timestamp of serial string read (time of day, total millisecond)
% This is subject to delete as well, because in the future we will have one
% timestamp (probbaly arduino logged timestamp) as the reference behaviour timestamp
% (which is then aligned to spikeglx time and then corrected based on photodiode)
% behaviour.wheel_arduino_read_time = resample(peripherals.Time,peripherals.corrected_sglxTime,60);

%% Wheel speed
% Real speed based on wheel raw input
% For actual speed related analysis, This should be used.

tick_to_cm_conversion = 3.1415*20/1024; % radius 20 cm and one circle is 1024 ticks;
speed = [0 diff(behaviour.wheel_raw_input*tick_to_cm_conversion)];
speed(speed<-100) = 0;% big change in speed often due to teleportation or wheel tick resetting
speed(speed>100) = 0;
behaviour.speed = speed./([0 diff(behaviour.sglxTime)]);% 

%% task_info
% Unlike behavour, task_info will save all the meta info about the task
% such as lap start time, lick time, reward delivery time in discrete time
% rather than continuous timeseries.

task_info = [];
% Task info about the onset of stimuli and/or stimuli info
task_info.stim_onset = stimTimes;% Stimulus onset time

if ~isempty(stimData.Properties.VariableNames)
    if contains(stimData.Properties.VariableNames,'stim_matrix') % For sparse noise
        task_info.stim_matrix = stimData.stim_matrix; %
    end
end

switch(StimulusName)
    case {'StaticGratings_short'}
        task_info.stim_orientation = readmatrix('X:\ibn-vision\CODE\DEV\BONSAI\Diao\dome_dual_DT\Grating_trials_short.CSV');
    case {'StaticGratings_long'}
        task_info.stim_orientation = readmatrix('X:\ibn-vision\CODE\DEV\BONSAI\Diao\dome_dual_DT\Grating_trials.CSV');
    case {'SparseNoise_fullscreen'}
        if contains(stimData.Properties.VariableNames,'stim_matrix')  % For sparse noise
            task_info.stim_matrix = stimData.stim_matrix; %
        end
end

if  ~isempty(photodiodeData)
    task_info.pd_on.sglxTime = photodiodeData.stim_on.sglxTime';
    task_info.pd_off.sglxTime = photodiodeData.stim_off.sglxTime';
    if isfield(photodiodeData,'photodiode_failure')
        behaviour.photodiode_failure = photodiodeData.photodiode_failure;
    end
else % if photodiode empty
    behaviour.photodiode_failure = 1;
end

end