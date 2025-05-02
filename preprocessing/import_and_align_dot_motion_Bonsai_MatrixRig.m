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
function   [behaviour,task_info,peripherals] = import_and_align_dot_motion_Bonsai_MatrixRig(options)


EPHYS_DATAPATH = options.EPHYS_DATAPATH;

if isfield(options,'bonsai_files_names')
    PERIPHERALS_DATAPATH = char(fullfile(options.BONSAI_DATAPATH,options.bonsai_files_names(contains(options.bonsai_files_names,'Wheel'))));
    TRIALDATA_DATAPATH = char(fullfile(options.BONSAI_DATAPATH,options.bonsai_files_names(contains(options.bonsai_files_names,'Trial'))));
    PHOTODIODE_DATAPATH = char(fullfile(options.BONSAI_DATAPATH,options.bonsai_files_names(contains(options.bonsai_files_names,'Photodiode'))));
    EYEDATA_DATAPATH = char(fullfile(options.BONSAI_DATAPATH,options.bonsai_files_names(contains(options.bonsai_files_names,'Video'))));% Old eye log
    DLC_EYEDATA_DATAPATH = char(fullfile(options.BONSAI_DATAPATH,options.bonsai_files_names(contains(options.bonsai_files_names,'DLC'))));% new eye log
    FACEDATA_DATAPATH = char(fullfile(options.BONSAI_DATAPATH,options.bonsai_files_names(contains(options.bonsai_files_names,'proc'))));% face data;
else
    PERIPHERALS_DATAPATH = options.PERIPHERALS_DATAPATH;
    EYEDATA_DATAPATH = options.EYEDATA_DATAPATH;
    TRIALDATA_DATAPATH = options.TRIALDATA_DATAPATH;
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

%[stimData,eyeData,peripherals,photodiodeData,stimTimes] = import_and_align_bonsai_logs_MatrixRig(EPHYS_DATAPATH,TRIALDATA_DATAPATH,PERIPHERALS_DATAPATH,EYEDATA_DATAPATH,PHOTODIODE_DATAPATH,options);
% [stimData,eyeData,peripherals,photodiodeData,stimTimes] = import_and_align_bonsai_logs(EPHYS_DATAPATH,TRIALDATA_DATAPATH,PERIPHERALS_DATAPATH,EYEDATA_DATAPATH,[],options);
% Step 1: Get stimulus times and data parameters
stimData = [];



% Step 2: import behavioural / eye etc

% Wheel data
peripherals= [];
thisTable = readtable(PERIPHERALS_DATAPATH,'Delimiter',{',','=','{','}','(',')'});
thisTable = rmmissing(thisTable,1); % Remove nan due to delimiter

for n = thisTable.Properties.VariableNames
    peripherals.(thisTable.Properties.VariableNames{n}) = thisTable.(thisTable.Properties.VariableNames{n});
end


% Video data
eyeData= [];
thisTable = readtable(EYEDATA_DATAPATH,'Delimiter',{',','=','{','}','(',')'});
thisTable = rmmissing(thisTable,1); % Remove nan due to delimiter

for n = thisTable.Properties.VariableNames
    eyeData.(thisTable.Properties.VariableNames{n}) = thisTable.(thisTable.Properties.VariableNames{n});
end

% Photodiode
if exist('PHOTODIODE_DATAPATH','var') && ~isempty(PHOTODIODE_DATAPATH)
    thisTable = readtable(PHOTODIODE_DATAPATH,'Delimiter',{',','=','{','}','(',')'});
    thisTable = rmmissing(thisTable,1); % Remove nan due to delimiter

    for n = thisTable.Properties.VariableNames
        photodiode.(thisTable.Properties.VariableNames{n}) = thisTable.(thisTable.Properties.VariableNames{n});
    end

else
    photodiode = []; photodiode_sync = [];
end

% Step 3: obtain imec sync pulse from precomputed or imec file
Nidq = [];
td = dir(fullfile(EPHYS_DATAPATH,'..','*NidqTimes*'));
if ~isempty(td)
    load(fullfile(EPHYS_DATAPATH,'..',td(1).name),'Nidq');
    syncTimes_ephys = Nidq.bonsai_sync_on;% use upswings currently
    syncPulse_ephys = Nidq.bonsai_sync;
    syncPulse_ephysTimes = Nidq.sglxTime;
else
    error('nidq and ephys sync pulse extraction and alignment not done!')
    return
end


% Step 4: align nidq and Bonsai logged data

%%%%%%
% For photodiode and convert to spike GLX timeThe
if sum(photodiode.AsyncPulse>0) <= 0 % no sync pulse detected so align photodiode to nidq in sglx time and then align peripehrals to photodiode sglxtime
    'warning!!! sync pulse missing'



        baseline_corrected_PD_nidq = Nidq.photodiode-movmean(Nidq.photodiode,120);
        initial_centroids = [ prctile(baseline_corrected_PD_nidq,10);prctile(baseline_corrected_PD_nidq,90)];
        PD_nidq_kmean = kmeans(baseline_corrected_PD_nidq',2,'Start', initial_centroids)-1;
        PD_nidq_onset = double([0;diff(PD_nidq_kmean)] > 0);

        baseline_corrected_PD_bonsai = photodiode.Photodiode-movmean(photodiode.Photodiode,120);
        initial_centroids = [ prctile(baseline_corrected_PD_bonsai,10);prctile(baseline_corrected_PD_bonsai,90)];
        PD_bonsai_kmean = kmeans(baseline_corrected_PD_bonsai,2,'Start', initial_centroids)-1;
        PD_bonsai_onset = double([0;diff(PD_bonsai_kmean)] > 0);
        
        photodiode.sglxTime = align_A_time_to_B_time(PD_bonsai_onset,photodiode.BonsaiTime,PD_nidq_onset,Nidq.sglxTime);


else
    if exist('photodiode','var') && ~isempty(photodiode)
        photodiode.Sync =  photodiode.AsyncPulse;
        photodiode.Time = photodiode.ArduinoTime;
    %     photodiode.Time = photodiode.BonsaiTime;
    
        [photodiode] = alignBonsaiToEphysSyncTimes(photodiode,syncTimes_ephys);
    %     [photodiode] = alignBonsaiToEphysSyncPulse(photodiode,syncPulse_ephys,syncPulse_ephysTimes);
        
    else
        photodiode = [];
    end
end


%%%%%%
% For wheel data and eyedaya and convert to spike GLX time-base
[C,ia,ic]=unique(photodiode.ArduinoTime);
peripherals.sglxTime = interp1(photodiode.ArduinoTime(ia),photodiode.sglxTime(ia),peripherals.ArduinoTime,'previous');
eyeData.sglxTime = interp1(photodiode.ArduinoTime(ia),photodiode.sglxTime(ia),eyeData.ArduinoTime,'previous');
wheelData = peripherals;


% extract stimuli time based on photodiode timestamp

photodiode.Photodiode_smoothed = smoothdata(photodiode.Photodiode,'movmedian',5);

pd_ON_OFF = photodiode.Photodiode_smoothed >= mean(photodiode.Photodiode); % find ON and OFF
pd_ON = find(diff(pd_ON_OFF) == 1); % Find ON
pd_OFF = find(diff(pd_ON_OFF) == -1); % Find OFF

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

blocks_ind = find(block_length>10); % sample rate of photodiode is 1000 per second so a block size should be quite large

   % plot(photodiode.sglxTime,photodiode.Photodiode_smoothed); hold on; scatter(photodiode.sglxTime(pd_OFF(blocks_ind)),50,'b');
   % scatter(photodiode.sglxTime(pd_ON(blocks_ind)),50,'r');

   % plot(Nidq.sglxTime,Nidq.photodiode/100);
photodiodeData = [];
photodiodeData.stim_on.sglxTime = photodiode.sglxTime(pd_ON(blocks_ind)');
photodiodeData.stim_off.sglxTime = photodiode.sglxTime(pd_OFF(blocks_ind)');

if photodiodeData.stim_off.sglxTime(end) > Nidq.sglxTime(end)
    % if last visual stimuli is later than last timestamp of the recording
    % then it is probably because of wrong alignment (happens normally only during short sessions such as checkerboard)
    % in this rare case use Nidq photodiode directly

    Nidq.photodiode_smoothed = smoothdata(Nidq.photodiode,'movmedian',50);

    pd_ON_OFF = Nidq.photodiode_smoothed >= mean(Nidq.photodiode); % find ON and OFF
    pd_ON = find(diff(pd_ON_OFF) == 1); % Find ON
    pd_OFF = find(diff(pd_ON_OFF) == -1); % Find OFF

    if length(pd_ON) == length(pd_OFF)
        block_length = abs(pd_ON-pd_OFF);
    elseif length(pd_ON)>length(pd_OFF)
        block_length = [];
        for i = 1:length(pd_OFF)
            block_length(i) = abs(pd_ON(i) - pd_OFF(i));
        end
    else
        for i = 1:length(pd_ON)
            block_length(i) = abs(pd_ON(i) - pd_OFF(i));
        end
    end

    blocks_ind = find(block_length>10); % sample rate of photodiode is 1000 per second so a block size should be quite large

    % plot(photodiode.sglxTime,photodiode.Photodiode_smoothed); hold on; scatter(photodiode.sglxTime(pd_OFF(blocks_ind)),50,'b');
    % scatter(photodiode.sglxTime(pd_ON(blocks_ind)),50,'r');

    % plot(Nidq.sglxTime,Nidq.photodiode/100);
    photodiodeData = [];
    photodiodeData.stim_on.sglxTime = Nidq.sglxTime(pd_ON(blocks_ind))';
    photodiodeData.stim_off.sglxTime = Nidq.sglxTime(pd_OFF(blocks_ind))';
    photodiodeData.warning = 'Nidq photodiode';
    disp('Nidq photodiode directly due to bad alignment between Bonsai and Ephys')
end


% Step 5: process wheel data (skipping this just save all peripheral data)
% wheelData.pos = peripherals.Wheel;
% wheelData.Time = peripherals.sglxTime;
% wheelData = import_processWheelEphys(wheelData,'gaussian',12);
% wheelData.sglxTime = wheelData.Time;

% Step 6: extract stimulus times

     stimTimes = photodiodeData.stim_on.sglxTime;

fprintf('\n'); warning('Havent yet enabled checking of number of stimulus expected vs number of photodiode changes'); fprintf('\n');

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
behaviour.wheel_frame_count = interp1(peripherals.sglxTime,peripherals.RenderFrameCount,peripherals.sglxTime(1):1/60:peripherals.sglxTime(end),'previous');

%% wheel spikeGLX time (Use this for ephys)
% behaviour.sglxTime = resample(peripherals.sglxTime,peripherals.sglxTime,60)';
behaviour.sglxTime = interp1(peripherals.sglxTime,peripherals.sglxTime,peripherals.sglxTime(1):1/60:peripherals.sglxTime(end),'linear');

% behaviour.sglxTime_uncorrected = resample(peripherals.sglxTime,peripherals.corrected_sglxTime,60)';
behaviour.tvec = behaviour.sglxTime;% standard timevec used for all behaviour related analysis

%% eye camera frame count
% behaviour.camera_frame_count = resample(peripherals.EyeFrameCount,peripherals.sglxTime,60)';
behaviour.camera_frame_count = interp1(eyeData.sglxTime,eyeData.VideoFrameCount,peripherals.sglxTime(1):1/60:peripherals.sglxTime(end),'previous');
% behaviour.wheel_reference_timestamp = peripherals.Time;% probably not needed....

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
behaviour.wheel_raw_input = interp1(peripherals.sglxTime,peripherals.Wheel,peripherals.sglxTime(1):1/60:peripherals.sglxTime(end),'linear');

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
trial_param = readtable(TRIALDATA_DATAPATH);
task_info.dot_velocity = trial_param.VelX1(1:length(stimTimes));
end