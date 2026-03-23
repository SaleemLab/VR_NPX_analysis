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
function   [behaviour,task_info,peripherals] = import_and_align_visual_stimuli_Bonsai_Ellie(StimulusName,options)


EPHYS_DATAPATH = options.EPHYS_DATAPATH;

if isfield(options,'bonsai_files_names')
    PERIPHERALS_DATAPATH = char(fullfile(options.BONSAI_DATAPATH,options.bonsai_files_names(contains(options.bonsai_files_names,'Video'))));
    TRIALDATA_DATAPATH = char(fullfile(options.BONSAI_DATAPATH,options.bonsai_files_names(contains(options.bonsai_files_names,'StimulusParams'))));
    PHOTODIODE_DATAPATH = char(fullfile(options.BONSAI_DATAPATH,options.bonsai_files_names(contains(options.bonsai_files_names,'Photodiode'))));
    EYEDATA_DATAPATH = char(fullfile(options.BONSAI_DATAPATH,options.bonsai_files_names(contains(options.bonsai_files_names,'Video'))));% Old eye log
    DLC_EYEDATA_DATAPATH = char(fullfile(options.BONSAI_DATAPATH,options.bonsai_files_names(contains(options.bonsai_files_names,'DLC'))));% new eye log
    FACEDATA_DATAPATH = char(fullfile(options.BONSAI_DATAPATH,options.bonsai_files_names(contains(options.bonsai_files_names,'proc'))));% face data;
else
    PERIPHERALS_DATAPATH = options.PERIPHERALS_DATAPATH;
    EYEDATA_DATAPATH = options.EYEDATA_DATAPATH;
    TRIALDATA_DATAPATH = options.TRIALDATA_DATAPATH;
end

if contains(options.bonsai_files_names{1,1},'SparseNoise')
    TRIALDATA_DATAPATH = char(fullfile(options.BONSAI_DATAPATH,options.bonsai_files_names(contains(options.bonsai_files_names,'.bin'))));
end


if contains(options.StimulusName,'Checkerboard')|contains(options.StimulusName,'FullScreenFlash')
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
    %otherMetrics = {'wheelspeed','pupilarea'};
    otherMetrics = {'pupilarea'};
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

[stimData,eyeData,peripherals,photodiodeData,stimTimes] = import_and_align_bonsai_logs_Ellie(EPHYS_DATAPATH,TRIALDATA_DATAPATH,PERIPHERALS_DATAPATH,EYEDATA_DATAPATH,PHOTODIODE_DATAPATH,options);
% [stimData,eyeData,peripherals,photodiodeData,stimTimes] = import_and_align_bonsai_logs(EPHYS_DATAPATH,TRIALDATA_DATAPATH,PERIPHERALS_DATAPATH,EYEDATA_DATAPATH,[],options);


if ~isempty(DLC_EYEDATA_DATAPATH)
    eye_data_raw = readmatrix(DLC_EYEDATA_DATAPATH);
    [pupil_ellipse,new_tracking_points_coordinates] = eye_data_conversion(eye_data_raw);
end

if ~isempty(FACEDATA_DATAPATH)
    facedata = load(FACEDATA_DATAPATH); % load face energy data
    if ~isfield(facedata,'motion_1')
        facedata = [];
        FACEDATA_DATAPATH = [];
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


if length(peripherals.sglxTime) > length(unique(peripherals.sglxTime))
    [C,ia,ic] = unique(peripherals.sglxTime);
    unique_id = find(diff(ic)==0);
    peripherals.sglxTime(unique_id+1)=nan;% make repeated values -> nan
end

all_fields = fieldnames(peripherals);

temp = peripherals.sglxTime;
for i = 1:length(all_fields)
    if length(peripherals.(all_fields{i})) == length(temp)
        peripherals.(all_fields{i})(isnan(temp))=[];
    end
end

if length(eyeData.sglxTime) > length(unique(eyeData.sglxTime))
    [C,ia,ic] = unique(eyeData.sglxTime);
    unique_id = find(diff(ic)==0);
    eyeData.sglxTime(unique_id+1)=nan;% make repeated values -> nan
end

all_fields = fieldnames(eyeData);
temp = eyeData.sglxTime;

for i = 1:length(all_fields)
    if length(eyeData.(all_fields{i})) == length(temp)
        eyeData.(all_fields{i})(isnan(temp))=[];
    end
end




%% wheel frame
% behaviour.wheel_frame_count = resample(peripherals.FrameNumber,peripherals.sglxTime,60)';
%behaviour.wheel_frame_count = interp1(peripherals.sglxTime,peripherals.RenderFrameCount,peripherals.sglxTime(1):1/60:peripherals.sglxTime(end),'previous');

%% wheel spikeGLX time (Use this for ephys)
% behaviour.sglxTime = resample(peripherals.sglxTime,peripherals.sglxTime,60)';
behaviour.sglxTime = interp1(peripherals.sglxTime,peripherals.sglxTime,peripherals.sglxTime(1):1/60:peripherals.sglxTime(end),'linear');

% behaviour.sglxTime_uncorrected = resample(peripherals.sglxTime,peripherals.corrected_sglxTime,60)';
behaviour.tvec = behaviour.sglxTime;% standard timevec used for all behaviour related analysis

%% eye camera frame count
% behaviour.camera_frame_count = resample(peripherals.EyeFrameCount,peripherals.sglxTime,60)';
behaviour.camera_frame_count = interp1(eyeData.sglxTime,eyeData.VideoFrameCount,peripherals.sglxTime(1):1/60:peripherals.sglxTime(end),'previous');
% behaviour.wheel_reference_timestamp = peripherals.Time;% probably not needed....
behaviour.arduino_time_sixtyHz = interp1(eyeData.sglxTime,eyeData.ArduinoTime,peripherals.sglxTime(1):1/60:peripherals.sglxTime(end),'linear');
behaviour.sync_pulse_sixtyHz = interp1(eyeData.sglxTime,eyeData.LEDAsync,peripherals.sglxTime(1):1/60:peripherals.sglxTime(end),'linear');
behaviour.pupil_X_sixtyHz = interp1(eyeData.sglxTime,eyeData.PupilCentroid_X,peripherals.sglxTime(1):1/60:peripherals.sglxTime(end),'linear');
behaviour.pupil_Y_sixtyHz = interp1(eyeData.sglxTime,eyeData.PupilCentroid_Y,peripherals.sglxTime(1):1/60:peripherals.sglxTime(end),'linear');
behaviour.pupil_area_sixtyHz = interp1(eyeData.sglxTime,eyeData.PupilArea,peripherals.sglxTime(1):1/60:peripherals.sglxTime(end),'linear');
behaviour.pupil_majlen_sixtyHz = interp1(eyeData.sglxTime,eyeData.PupilMajLength,peripherals.sglxTime(1):1/60:peripherals.sglxTime(end),'linear');
behaviour.pupil_minlen_sixtyHz = interp1(eyeData.sglxTime,eyeData.PupilMinLength,peripherals.sglxTime(1):1/60:peripherals.sglxTime(end),'linear');


%% eye data
if ~isempty(DLC_EYEDATA_DATAPATH)
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
%behaviour.wheel_raw_input = interp1(peripherals.sglxTime,peripherals.Wheel,peripherals.sglxTime(1):1/60:peripherals.sglxTime(end),'linear');

% Not needed
% %% wheel left lick
% behaviour.wheel_lick_L = peripherals.LickL;
% 
% %% wheel right lick
% behaviour.wheel_lick_R = peripherals.LickR;

% %% wheel quad state
% behaviour.wheel_quad_state = peripherals.QuadState;

%% Wheel speed
% Real speed based on wheel raw input
% For actual speed related analysis, This should be used.

%tick_to_cm_conversion = 3.1415*20/1024; % radius 20 cm and one circle is 1024 ticks;
%speed = [0 diff(behaviour.wheel_raw_input*tick_to_cm_conversion)];
%speed(speed<-100) = 0;% big change in speed often due to teleportation or wheel tick resetting
%speed(speed>100) = 0;
%behaviour.speed = speed./([0 diff(behaviour.sglxTime)]);% 

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
    case {'TRAIN', 'TRAIN_1', 'TRAIN_2', 'OP_Tuning', 'Direction_Tuning', 'DCBA', 'OMIT', 'E_CD', 'D_CD', 'D___', 'A___',...
            'A_1000ms', 'A_1000ms_1', 'A_1000ms_2', 'A_1000ms_25pc', 'A_1000ms_50pc', 'A_1000ms_75pc', 'E_1000ms', 'TRAIN250',...
            'A_50ms', 'A_500ms', 'A_200ms', 'A_300ms', 'OMIT50grey', 'ADCD', 'lowcontB', 'lowcontDsubbingB', 'lowcontTRAIN',...
            'F_150ms', 'F_150ms_1', 'F_150ms_2', 'F_1000ms'}
        task_info.stim_delay = stimData.Delay;
        task_info.stim_contrast = stimData.Contrast;
        task_info.stim_orientation = stimData.Orientation;

    case {'GAVNIK_ABCD', 'GAVNIK_ABCD_1', 'GAVNIK_ABCD_2', 'GAVNIK DCBA', 'GAVNIK_A_CD', 'GAVNIK_E_CD', 'GAVNIK_D_CD',...
            'GAVNIK200_ABCD', 'GAVNIK200_A_CD', 'GAVNIK200_E_CD', 'GAVNIK250_ABCD', 'GAVNIK250_ABCD_1', 'GAVNIK250_ABCD_2',...
            'GAVNIK250_A_CD', 'GAVNIK250_E_CD', 'GAVNIK250 DCBA',...
            'GAVNIK_D___', 'GAVNIK_A___'}
        task_info.stim_orientation = stimData.Orientation;
        
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
    if isfield(photodiodeData,'warning')
        task_info.warning = photodiodeData.warning;
    end
else % if photodiode empty
    behaviour.photodiode_failure = 1;
end

end