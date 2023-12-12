% function [OE,EyeDat] = parseRawToCollated_spikes(thisAnimal,thisSession,thisFileNum,options)
%
% Inputs
%       thisAnimal  [= 'M19062']
%       thisSession [= '20190911']
%       thisFileNum [= 0]
%       options -
%           options.BonsaiPath = fullfile(serverPath,'DATA','SUBJECTS',upper(thisAnimal),'Extras',upper(thisSession));
%           options.conToppath = fullfile(serverPath,'DATA','SUBJECTS',upper(thisAnimal),'Ephys',upper(thisSession));
%               NB if following fields exist they will be used to directly
%                   access the eye data and OE file
%               options.thisOEfileName - full path to OE file
%               options.thisEYEfileName - fullpath to EYE file
%           options.sync_input = 'ch';
%           options.sync_channel = 20;
%           options.wheel_input = 'ch';
%           options.wheel_channels = [22:23];
%           options.includeWheel = 1;
%           options.UseOneSignal = 0; % 0 or 1
%           options.chanNumbersToInclude = 9;
%           options.EyeSR = 31;
%           options.stim_dur = 0.5; % For PD calculations only
%           options.photo_br = 1;
%           options.NewOESampleRate = 1000;
%
% Dependencies
%       load_eyeWheelPdData
%       readCONandPDandWheelNew
%       calculateSpeedFromWheelData
%       EyePreprocess
%       ts_EyePD
%
% History
%   SGS 8th March 2020 Wrote it adapted from SRP_RunCONtoMAT_EPhys
%   SGS 21st April 2020 Wrote it adapted from parseRawToCollated and  SparseNoiseSyncR
%
function [OE,EyeDat,options] = parseRawToCollated(thisAnimal,thisSession,thisFileNum,options)

if ~exist('thisAnimal', 'var')
    thisAnimal = 'M19075';
end
if ~exist('thisSession', 'var')
    thisSession = '20190617';
end
if ~exist('thisFileNum', 'var')
    thisFileNum = 0;
end

% Default options for debugging
if ~exist('options', 'var')
    if ispc
        if exist('X:\ibn-vision','dir')==7
            serverPath = 'X:\ibn-vision\';
        else
            serverPath = 'X:\';
        end
    elseif ismac
        serverPath = '/Users/s.solomon/Filestore/Research2/ibn-vision';
    end
    
    % Temporary (for M19062-20190911)
    options.BonsaiPath = fullfile(serverPath,'DATA','SUBJECTS',upper(thisAnimal),'Extras',upper(thisSession));
    options.conToppath = fullfile(serverPath,'DATA','SUBJECTS',upper(thisAnimal),'Ephys',upper(thisSession));
    options.sync_input = 'ch';
    options.sync_channel = 20;
    options.wheel_input = 'ch';
    options.wheel_channels = 22:23;
    options.includeWheel = 1;
    options.UseOneSignal = 0; % 0 or 1
    options.chanNumbersToInclude = 9;
    options.EyeSR = 31;
    options.stim_dur = 0.5; % For PD calculations only
    options.photo_br = 1;
    options.NewOESampleRate = 1000;
end

%%%%%%%
% Find Ephys folder matching with thisSession
% to get the ephys, photodiode and wheel data fromthe open ephys file
try
    % User may have specified OE file directly as part of options structure
    if isfield(options, 'thisOEfileName') && ~isempty(options.thisOEfileName)
        thisOEfileName = fullfile(options.conToppath,options.thisOEfileName);
    else
        thisOEfile = fullfile(options.conToppath, [thisAnimal,'_*_', num2str(thisFileNum)]);
        folder_dir = dir(thisOEfile);
        folder_name = folder_dir.name;
        thisOEfileName = fullfile(options.conToppath,folder_name);
    end
    [OELFPData,OEttl,tOEwheelData, fileinfo] = readCONandPDandWheelNew(thisOEfileName,options);
    mainSampleRate = fileinfo.data_info.header.sampleRate;
%     [OELFPData,OEttl,tOEwheelData, fileinfo] = readSpikesAndPDandWheel(thisOEfileName,options);
    
    
    OE.mainSampleRate = mainSampleRate;
    OE.OEFileName = thisOEfileName;
    OE.fileinfo = fileinfo;
    if options.includeWheel==1
        %Calculate speed
        OEwheelData = calculateSpeedFromWheelData(tOEwheelData(2,:),tOEwheelData(1,:), options);
    else
        OEwheelData = [];
    end
catch
    if options.logErrorsFlag
        fprintf(options.errorlogfile,'\n%s - %s - %s\n\n',thisAnimal,thisSession,date);
        fprintf(options.errorlogfile,'\n\t error in parseRawToCollated:readCONandPDandWheelNew');
    end
end

%%%%%
% Load the eye, wheel and Pd data from the bonsai file
% Find the relevant Bonsai file Read the file and return a table containing stim and photodiode params
try
    thisEYEfileName = [];
    if isfield(options, 'thisEYEfileName') && ~isempty(options.thisEYEfileName)
        thisEYEfileName = fullfile(options.BonsaiPath,options.thisEYEfileName);
    else
        thisBonsaiFile = fullfile(options.BonsaiPath, [thisAnimal,'_eyeWheelPdData_', num2str(thisFileNum),'_*']);
        folder_dir = dir(thisBonsaiFile);
        if isempty(folder_dir)
            thisBonsaiFile = fullfile(options.BonsaiPath, [thisAnimal,'_eyeWheelPdData_*']);
            folder_dir = dir(thisBonsaiFile);
        end
        if ~isempty(folder_dir)
            folder_name = folder_dir.name;
            thisEYEfileName = fullfile(options.BonsaiPath,folder_name);
        end
    end
    if ~isempty(thisEYEfileName)
        EyeWhPd_dat = load_eyeWheelPdData(thisEYEfileName);
        
        % Do some preprocessing
        [EyeWhPd_dat, EyeTrackerParams] = EyePreprocess(EyeWhPd_dat);
        videoRate = 1/(mean(diff(EyeWhPd_dat.pdMsSinceStartOfDay))/1000);
        
        EyeDat.SampleRate = videoRate;
        EyeDat.EyeFileName = thisEYEfileName;
        EyeDat.EyeTrackerParams = EyeTrackerParams;
        % Get the photodiode change based on the bonsai file
        if options.photo_br == 0
            switch options.paradigm
                case 'srp'
                    options.target_upPhases = length(OEttl.upPhases);
            end
        end
        [EyeupPhases, EyedownPhases] = ts_EyePD(EyeWhPd_dat,options);
        EyeDat.upPhases = EyeupPhases;
        
        % Store some stuff
        if options.photo_br == 0
            switch options.paradigm
                case 'srp'
                    % Otherwise it is very unstable, and we should just use
                    % the upPhases + the stimulus duration
                    EyeDat.downPhases = EyeDat.upPhases+round(EyeDat.SampleRate*options.stim_dur); %
                otherwise
                    EyeDat.downPhases = [];
            end
        elseif options.photo_br == 1
            EyeDat.downPhases = EyedownPhases;
        end
        
        % Store traces of PD around phases
        nSamples = fix(0.1*EyeDat.SampleRate);
        tcol = repmat(-nSamples:nSamples,length(EyeupPhases),1);
        trow = repmat(EyeupPhases(:),1,size(tcol,2));
        indices = trow+tcol;
        indices(indices>length(EyeWhPd_dat.pdOutput)) = length(EyeWhPd_dat.pdOutput);
        EyeDat.upPhases_traces = EyeWhPd_dat.pdOutput(indices);
        
        if ~isempty(EyeDat.downPhases)
            tcol = repmat(-nSamples:nSamples,length(EyeDat.downPhases),1);
            trow = repmat(EyeDat.downPhases(:),1,size(tcol,2));
            indices = trow+tcol;
            indices(indices>length(EyeWhPd_dat.pdOutput)) = length(EyeWhPd_dat.pdOutput);
            EyeDat.downPhases_traces = EyeWhPd_dat.pdOutput(indices);
        else
            EyeDat.downPhases_traces = [];
        end
        
    else
        warning('No bonsai file for eye data (%s): not processing eye data',thisBonsaiFile)
        EyeDat = []; % No save file for this session M20029-20200211
    end
catch
    if options.logErrorsFlag
        fprintf(options.errorlogfile,'\n%s - %s - %s\n\n',thisAnimal,thisSession,date);
        fprintf(options.errorlogfile,'\n\t error in parseRawToCollated:ts_EyePD');
    end
end
% Downsample OE data and wheel data if necessary
if isfield(options,'NewOESampleRate') && ~isempty(options.NewOESampleRate) && options.NewOESampleRate < mainSampleRate
    [OELFPData,OEttl,OEwheelData] = downSampleLFPAndWheel(OELFPData,OEttl,OEwheelData,mainSampleRate,options);
    OE.SampleRate = options.NewOESampleRate;
else
    OE.SampleRate = mainSampleRate;
end
%%%%%%%%%
% Store traces of PD around phases
try
    nSamples = fix(0.1*OE.SampleRate);
    tcol = repmat(-nSamples:nSamples,length(OEttl.upPhases),1);
    trow = repmat(OEttl.upPhases(:),1,size(tcol,2));
    tIndices = trow+tcol;
    tIndices(tIndices>length(OEttl.ttlData)) = length(OEttl.ttlData);
    OEttl.upPhases_traces = OEttl.ttlData(tIndices);
    tcol = repmat(-nSamples:nSamples,length(OEttl.downPhases),1);
    trow = repmat(OEttl.downPhases(:),1,size(tcol,2));
    tIndices = trow+tcol;
    tIndices(tIndices>length(OEttl.ttlData)) = length(OEttl.ttlData);
    OEttl.downPhases_traces = OEttl.ttlData(tIndices);
    OEttl.ttlData = [];
catch
    if options.logErrorsFlag
        fprintf(options.errorlogfile,'\n%s - %s - %s\n\n',thisAnimal,thisSession,date);
        fprintf(options.errorlogfile,'\n\t error in storing PD traces');
    end
end

% % % Downsample eye video data if necessary [SGS: STILL TO DO]
% % if ~isempty(EyeDat) && isfield(options,'NewOESampleRate') && ~isempty(options.NewOESampleRate) && options.NewOESampleRate < videoRate
% %     error('Downsampling of video data not yet done')
% %     OEdecimationRate = videoRate/options.NewOESampleRate;
% %     if rem(OEdecimationRate,1)    % if not an integer decimation rate
% %         error('Target sample rate is %3.1f but needs to be an integer multiple of main videoRate rate of %3.1f Hz',options.NewOESampleRate,videoRate)
% %     end
% % end

%%%%%%%%%
% Collapse up and down phases to get ttl.timestamps
try
    if ~isempty(EyeDat.downPhases)
        Eye_phases = zeros(1,length(EyeDat.upPhases)+length(EyeDat.downPhases));
        Eye_phases(1:2:length(EyeDat.upPhases)*2) = EyeDat.upPhases;
        Eye_phases(2:2:length(EyeDat.downPhases)*2) = EyeDat.downPhases;
    else
        Eye_phases = EyeDat.upPhases;
    end
    EyeDat.ttl_timestamps = Eye_phases/EyeDat.SampleRate;
catch
    if options.logErrorsFlag
        fprintf(options.errorlogfile,'\n%s - %s - %s\n\n',thisAnimal,thisSession,date);
        fprintf(options.errorlogfile,'\n\t error in collapsing Eye_phases');
    end
end
try
    all_phases = zeros(1,length(OEttl.upPhases)+length(OEttl.downPhases));
    all_phases(1:2:length(OEttl.upPhases)*2) = OEttl.upPhases;
    all_phases(2:2:length(OEttl.downPhases)*2) = OEttl.downPhases;
    OEttl.ttl_timestamps = all_phases./OE.SampleRate;
catch
    if options.logErrorsFlag
        fprintf(options.errorlogfile,'\n%s - %s - %s\n\n',thisAnimal,thisSession,date);
        fprintf(options.errorlogfile,'\n\t error in collapsing OE all_phases');
    end
end

%%%%%%%%%
% Align timing of OE data and Bonsai data
try
    if ~isempty(EyeWhPd_dat)
        % Check to see that the same number of PD are in EYE and OE
        % Find scale and shift factors required to align the two recordings
        % Shift
        eye_first_ts = EyeupPhases(1)/EyeDat.SampleRate;
        oe_first_ts = OEttl.upPhases(1)/OE.SampleRate;
        EYE_OFFSET = oe_first_ts-eye_first_ts;
        %         EYE_OFFSET = OEttl.ttl_timestamps(1)-EyeDat.ttl_timestamps(1);
        fprintf('\nEYE recordings are %3.4fs different to OE recordings',EYE_OFFSET)
        % Offset the eye data samples
        EyeDat.Offset_reOE_applied = EYE_OFFSET;
        EyeDat.ttl_timestamps = EyeDat.ttl_timestamps+EYE_OFFSET;
        EyeDat.data = EyeWhPd_dat(:,{'EyeX_deg','EyeY_deg','EyeArea','EyeArea_mm2','valid'});
        EyeDat.data.timestamps = EyeWhPd_dat.eyeMsSinceStartOfDay/1000+EYE_OFFSET;
        
        % Scale
        if length(EyeDat.ttl_timestamps) ~= length(OEttl.ttl_timestamps)
            if length(EyeDat.ttl_timestamps) > length(OEttl.ttl_timestamps)
                wording = 'MORE';
            else
                wording = 'LESS';
            end
            warning('%s PD reversals in EYE (%01d) than OE (%01d)', wording,length(EyeDat.ttl_timestamps),length(OEttl.ttl_timestamps))
        end
        EYE_SCALE = range(EyeDat.ttl_timestamps)/range(OEttl.ttl_timestamps);
        if abs(EYE_SCALE-1) > 0.001
            if EYE_SCALE-1 > 0
                wording = 'LONGER';
            else
                wording = 'SHORTER';
            end
            warning('Timescale in EYE is %s than OE', wording)
        end
        EyeDat.Scale_reOE_unapplied = EYE_SCALE;
        % Could scale but not currently (seems to do more harm than good)
        % eg.    Eye_ttl.ttl_timestamps = (EYE_SCALE*(Eye_ttl.ttl_timestamps-Eye_ttl.ttl_timestamps(1))+Eye_ttl.ttl_timestamps(1))+EYE_OFFSET;
    end
catch
    if options.logErrorsFlag
        fprintf(options.errorlogfile,'\n%s - %s - %s\n\n',thisAnimal,thisSession,date);
        fprintf(options.errorlogfile,'\n\t error in aligning EYE to OE');
    end
end

% Return
OE.OELFPData = OELFPData;
OE.OEttl = OEttl;
if options.includeWheel==1
    OE.OEwheelData = OEwheelData;
end
OE.Animal = thisAnimal;
OE.Session = thisSession;
OE.Filenum = thisFileNum;
if ~exist('EyeDat','var')
    EyeDat = [];
else
    EyeDat.Animal = thisAnimal;
    EyeDat.Session = thisSession;
    EyeDat.Filenum = thisFileNum;
end