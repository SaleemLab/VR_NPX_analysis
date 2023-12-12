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
function [MousePos,peripherals,photodiodeData,eyeData] = import_and_align_Masa_VR_Bonsai_old(StimulusName,options)

% Program to import bonsai logs and eye data and align to ePhys data and
% align the data based on the delay between quad and photodioide signal

% Find Bonsai files
bonsai_files = dir(fullfile(options.BONSAI_DATAPATH,[StimulusName, '*']));
bonsai_files_names = {bonsai_files.name};

PHOTODIODE_DATAPATH = bonsai_files_names(contains(bonsai_files_names,'PDAsync'));
PERIPHERALS_DATAPATH = bonsai_files_names(contains(bonsai_files_names,'WheelLog'));
EYEDATA_DATAPATH = bonsai_files_names(contains(bonsai_files_names,'EyeLog'));
STIMULI_DATAPATH = bonsai_files_names(contains(bonsai_files_names,'MousePos'));
BONSAI_DATAPATH = options.BONSAI_DATAPATH;

DIR = dir(fullfile(options.EPHYS_DATAPATH,'*syncpulseTimes.mat'))

if ~isempty(DIR)
    load(fullfile(options.EPHYS_DATAPATH,DIR.name))
    syncTimes_ephys = syncTimes_ephys.on; % use upswings currently
else
    [AP_FILE,LF_FILE] = findImecBinFile(options.EPHYS_DATAPATH);
%     FILE_TO_USE = AP_FILE;
    binpath = fullfile(options.EPHYS_DATAPATH,AP_FILE);
    %     binpath = fullfile(options.EPHYS_DATAPATH,LF_FILE);
    syncTimes_ephys = SGLXextractSyncPulseFromBinFile(binpath);
    parseDate = date;
    [~,fname] = fileparts(options.EPHYS_DATAPATH);
    save(fullfile(options.EPHYS_DATAPATH,[fname,'_syncpulseTimes.mat']),'syncTimes_ephys','parseDate');
    syncTimes_ephys = syncTimes_ephys.on; % use upswings currently
end

% Import wheel data, eye data, position data and photodiodide data 
% and synchronise them to async pulse 
peripherals = import_bonsai_peripherals(fullfile([BONSAI_DATAPATH,'\',char(PERIPHERALS_DATAPATH)]));
[peripherals] = alignBonsaiToEphysSyncTimes(peripherals,syncTimes_ephys);

% eyeData = import_bonsai_eyedata(fullfile([BONSAI_DATAPATH,'\',char(EYEDATA_DATAPATH)]));
% [eyeData] = alignBonsaiToEphysSyncTimes(eyeData,syncTimes_ephys);
eyeData = [];

if exist('PHOTODIODE_DATAPATH','var') && ~isempty(PHOTODIODE_DATAPATH)
    [photodiode, photodiode_sync] = import_bonsai_photodiode(fullfile([BONSAI_DATAPATH,'\',char(PHOTODIODE_DATAPATH)]));
    [photodiode] = alignBonsaiToEphysSyncTimes(photodiode,syncTimes_ephys);
else
    photodiode = []; photodiode_sync = [];
end


% Load Position data (may put these into Peripheral in the future)
temp = readtable(fullfile([BONSAI_DATAPATH,'\',char(STIMULI_DATAPATH)]));
tic
b=regexp(table2array(temp(:,4)),'\d+(\.)?(\d+)?','match');
MousePos.pos = str2double([b{:}]);
MousePos.quad=table2array(temp(:,6));
MousePos.Sync=table2array(temp(:,5));
b=regexp(table2array(temp(:,7)),'\d+(\.)?(\d+)?','match');
MousePos.Time = str2double([b{:}])';
MousePos.pos = MousePos.pos(1:length(MousePos.Time))';
toc
[MousePos] = alignBonsaiToEphysSyncTimes(MousePos,syncTimes_ephys);

% Calculate wheel speed at each time point
tick_to_cm_conversion = 3.1415*20/1024; % radius 20 cm and one circle is 1024 ticks;
speed = [0; diff(peripherals.Wheel*tick_to_cm_conversion)];
speed(speed<-100) = 0;
speed(speed>100) = 0;
peripherals.speed = speed;

% 
idx_trial_start = find(diff(peripherals.QuadState)==1) + 1;
idx_trial_end = find(diff(peripherals.QuadState)==-1) + 1;
frame_trial_start = peripherals.FrameNumber(idx_trial_start);
frame_trial_end = peripherals.FrameNumber(idx_trial_end);
temp_tbl.start_time = peripherals.sglxTime(idx_trial_start);
temp_tbl.end_time = peripherals.sglxTime(idx_trial_end);
nTrials = length(temp_tbl.start_time);
% get trial start/end times from photodiode signal
pd_thresh_up = 400;
pd_thresh_down = 100;
for itrial =1:length(temp_tbl.start_time);
    % find frame where Bonsai asked to change quad
    frameidx_start = find(photodiode.FrameNumber <= frame_trial_start(itrial),1,'last');
    % find next index where the photodiode detected a quad change
    temp_idx = find(photodiode.Photodiode(frameidx_start:end)>pd_thresh_up,1,'first');
    idx_start = frameidx_start+temp_idx-1;
    % convert to spike GLX time
    pdstart(itrial) = photodiode.sglxTime(idx_start);
end

for itrial =1:length(temp_tbl.end_time)
    % to get photodiode down we need to do some median filtering to
    % remove some awkward down swings when the quad is white
    frameidx_end = find(photodiode.FrameNumber >= frame_trial_end(itrial),1,'first');
    temp_idx = find(smoothdata(photodiode.Photodiode(frameidx_end:end),'movmedian', 5)<pd_thresh_down,1,'first');
    idx_end = frameidx_end+temp_idx-1;
    pdend(itrial) = photodiode.sglxTime(idx_end);
end

photodiodeData.stim_on.sglxTime = pdstart';
photodiodeData.stim_off.sglxTime = pdend';

switch StimulusName
    case 'replay_Masa2tracks'
        [MousePos] = alignBonsaiToPhotodiode(MousePos,sort([photodiodeData.stim_on.sglxTime; photodiodeData.stim_off.sglxTime]),'replay');
    case 'Masa2tracks'
        [MousePos] = alignBonsaiToPhotodiode(MousePos,sort([photodiodeData.stim_on.sglxTime; photodiodeData.stim_off.sglxTime]),[]);
end

