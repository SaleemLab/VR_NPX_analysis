% function [stimulusData,eyeData,wheelData,photodiodeData,stimTimes] = import_and_align_bonsai_logs(EPHYS_DATAPATH,TRIALDATA_DATAPATH, PERIPHERALS_DATAPATH,EYEDATA_DATAPATH)
%
% Program to import bonsai logs and eye data and align to ePhys data
% Modified from [stimulusData,eyeData,wheelData,photodiodeData,stimTimes] = importAndAlignBonsaiLogs(EPHYS_DATAPATH,TRIALDATA_DATAPATH, PERIPHERALS_DATAPATH,EYEDATA_DATAPATH)
% Now called import_and_align_bonsai_logs (just for differentiation)
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
%           MT 19/01/2024 wheelData is now just peripherals
% Dependencies:
%   import_bonsai_peripherals
%   import_bonsai_eyedata
%   import_BonVisionParams
%   import_processWheelEphys (currently disabled)
%   alignBonsaiAndEphysData
%
%   Uses the following libraries, which are included in the folder
%   https://billkarsh.github.io/SpikeGLX/ - SpikeGLX_Datafile_Tools
%
function [stimData,eyeData,wheelData,photodiodeData,stimTimes] = import_and_align_bonsai_logs_MatrixRig(EPHYS_DATAPATH,TRIALDATA_DATAPATH,PERIPHERALS_DATAPATH,EYEDATA_DATAPATH,PHOTODIODE_DATAPATH,options)

% Step 1: Get stimulus times and data parameters
if ~strcmp(options.paradigm,'photodiode_timestamp')
    [stimData,StimulusName] = import_BonVisionParams(TRIALDATA_DATAPATH);
else % If based on photodiode timestamp in case trial info not saved properly
    % Parse for one of the following stimulus names (you can add to the list here, and to the switch case below)
    specialFileTypes = {'StaticGratings_short','StaticGratings_long','StaticGratings','DriftingGratings','SparseNoise','SparseNoise_fullscreen',...
        'StaticTF','Grey','GratingsPhase','OriAdapt','PosAdapt','DriftingTF','Checkerboard','FullScreenFlash'};
    % See if there is a match within the file name that indicates that this is
    % from the special list
    StimulusTypeMatcher = zeros(1,length(specialFileTypes));
    for tt = 1:length(specialFileTypes)
        StimulusTypeMatcher(tt) = contains(TRIALDATA_DATAPATH,specialFileTypes{tt});
    end
    if length(find(StimulusTypeMatcher)) == 1
        StimulusName = specialFileTypes{find(StimulusTypeMatcher)};
    elseif contains(filepath,'SparseNoise_fullscreen') == 1 % If equal SparseNoise_fullscreen
        StimulusName = 'SparseNoise_fullscreen';
    end

    stimData = [];
    stimData.Properties.VariableNames = [];
end

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
else
    error('nidq and ephys sync pulse extraction and alignment not done!')
    return
end


% Step 4: align nidq and Bonsai logged data

%%%%%%
% For photodiode and convert to spike GLX time-base
if exist('photodiode','var') && ~isempty(photodiode)
    photodiode.Sync =  photodiode.AsyncPulse;
    photodiode.Time = photodiode.ArduinoTime;

    [photodiode] = alignBonsaiToEphysSyncTimes(photodiode,syncTimes_ephys);
else
    photodiode = [];
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

%    plot(photodiode.Photodiode_smoothed); hold on; scatter(pd_OFF(blocks_ind),50);
photodiodeData = [];
photodiodeData.stim_on.sglxTime = photodiode.sglxTime(pd_ON(blocks_ind)');
photodiodeData.stim_off.sglxTime = photodiode.sglxTime(pd_OFF(blocks_ind)');

% Step 5: process wheel data (skipping this just save all peripheral data)
% wheelData.pos = peripherals.Wheel;
% wheelData.Time = peripherals.sglxTime;
% wheelData = import_processWheelEphys(wheelData,'gaussian',12);
% wheelData.sglxTime = wheelData.Time;

% Step 6: extract stimulus times
switch(lower(StimulusName))
    case {'sparsenoise','sparsenoise_fullscreen','checkerboard'}
        stimTimes = sort(cat(1,photodiodeData.stim_on.sglxTime,photodiodeData.stim_off.sglxTime),'ascend');
        %         fprintf('\n\tUsing both photodiode upward and downward transitions as timestamps for stimulus timing')
        if contains(StimulusName,'SparseNoise')
            
            average_stim_duration = (stimTimes(end)-stimTimes(1))/size(stimData,1);
            stim_id = (find(diff(stimTimes)>average_stim_duration*2));

            if length(stim_id) == 1 % if happended once we can remove all the stimuli stime between these two points (happens rarely)
                stimData(stim_id:stim_id+size(stimData,1)-length(stimTimes)-1,:)=[];
            elseif length(stim_id) <4 % if happended few times we can remove all the stimuli stime between these two points (happens rarely)
                stimData(stim_id(1):end,:)=[];
                stimTimes(stim_id:end)=[];
            end

        end

    case {'staticgratings','staticgratings_short','staticgratings_long'}
        stimTimes = sort(cat(1,photodiodeData.stim_on.sglxTime),'ascend');
    case 'grey'
        stimTimes = photodiodeData.stim_on.sglxTime;
    case {'oriadapt','posadapt'}
        %         PD_timestamps = zeros(1,length(photodiodeData.stim_on.bonsaiTime)+length(photodiodeData.stim_off.bonsaiTime));
        %         PD_timestamps(1:2:length(photodiodeData.stim_off.bonsaiTime)*2) = photodiodeData.stim_off.bonsaiTime;
        %         PD_timestamps(2:2:length(photodiodeData.stim_on.bonsaiTime)*2) = photodiodeData.stim_on.bonsaiTime;
        %
        %         OE_ts = zeros(1,length(photodiodeData.stim_on.sglxTime)+length(photodiodeData.stim_off.sglxTime));
        %         OE_ts(1:2:length(photodiodeData.stim_off.sglxTime)*2) = photodiodeData.stim_off.sglxTime;
        %         OE_ts(2:2:length(photodiodeData.stim_on.sglxTime)*2) = photodiodeData.stim_on.sglxTime;

        PD_timestamps = photodiodeData.stim_on.bonsaiTime;

        OE_ts = photodiodeData.stim_on.sglxTime;

        % Do the correction for every photodiode reversal
        stimTimes = stimData.Time./1000;

        for i=1:length(PD_timestamps)-1
            bonsai_ts = PD_timestamps(i);
            oe_ts = OE_ts(i);
            offset(i) = bonsai_ts-oe_ts;
            if i==1
                inds = find(stimData.Time./1000<=PD_timestamps(i));
                stimTimes(inds) = stimTimes(inds)-offset(i);
            end
            inds = find(stimData.Time./1000>=PD_timestamps(i) & stimData.Time./1000<PD_timestamps(i+1));
            stimTimes(inds) = stimTimes(inds)-offset(i);
        end
        i = length(PD_timestamps);
        bonsai_ts = PD_timestamps(i);
        oe_ts = OE_ts(i);
        offset(i) = bonsai_ts-oe_ts;
        inds = find(stimData.Time./1000>PD_timestamps(i));
        stimTimes(inds) = stimTimes(inds)-offset(i);
        stimData.stimTimes = stimTimes;


    otherwise
        % Either use the stimulus times from the stimulus log or photodiode
        %
        % For use of stimulus logged times...
        % stimTimes = stimData.sglxTime;
        % fprintf('\n\tUsing stimulus logged times as timestamps for stimulus timing')
        %
        % For use times from the photodiode upswing
        stimTimes = photodiodeData.stim_on.sglxTime;
        fprintf('\n\tUsing photodiode upward transitions as timestamps for stimulus timing')
end
fprintf('\n'); warning('Havent yet enabled checking of number of stimulus expected vs number of photodiode changes'); fprintf('\n');

