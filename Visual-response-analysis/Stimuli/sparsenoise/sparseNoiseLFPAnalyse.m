% Program to take in Bonsai, EPhys and Eye data and organise it appropriately to make RF maps
% Inputs - animal, session, filenum
%            - options, including
% Outputs
function [initMap, lfp_data,wheel_data,eye_data] = sparseNoiseLFPAnalyse(thisAnimal,thisSession,thisFileNum,options)

%%%%%%%%%
% Dummy options for loading data (options from getParametersForSRPexperiments for this M20037/20200221)
if ~exist('options','var')
    options.photo_br = 1; % 0 or 1
    options.sync_input = 'adc';
    options.sync_channel = 1;
    options.wheel_input = 'adc';
    options.wheel_channels = 3:4;
    options.includeWheel = 1;
    options.UseOneSignal = 0;
    options.chanNumbersToInclude = 9;
    options.EyeSR = 31;
    options.baselineTimeForPDInS = 5;
    options.logErrorsFlag = 0;
end
%%%%%%%%%

%%%%%%%%%
% Other defaults
if nargin < 4
    thisAnimal = 'M19165';
    thisSession = '191219';
    thisFileNum = 6;
end

if~isfield(options,'serverPath')
    if ispc
        if exist('X:\ibn-vision','dir')==7
            options.serverPath = 'X:\ibn-vision\';
        else
            options.serverPath = 'X:\';
        end
    elseif ismac
%         options.serverPath = '/Users/s.solomon/Desktop/TEMP/';
        options.serverPath = '/Users/s.solomon/Filestore/Research2/ibn-vision/';
    end
end
if~isfield(options,'plotflag')
    options.plotflag = 1;
end
%%%%%%%%%

%%%%%%%%%
% Some defaults including known properties of the sparse noise stimulus (should be logged somewhere?)
if~isfield(options,'stim_dur')
    options.stim_dur = 0.5;
end
if~isfield(options,'grid_size')
    options.grid_size = [10 10]; % 10x10 x y
end
if~isfield(options,'mapSampleRate')
    options.mapSampleRate = 60; % Hz
end
%%%%%%%%%

%%%%%%%%%
% Set up some paths if not provided
if ~isfield(options,'BonsaiPath')
    options.BonsaiPath = fullfile(options.serverPath,'DATA','SUBJECTS',upper(thisAnimal),'Extras',upper(thisSession));
end
if ~isfield(options,'conToppath')
    options.conToppath = fullfile(options.serverPath,'DATA','SUBJECTS',upper(thisAnimal),'Ephys',upper(thisSession));
end
%%%%%%%%%

%%%%%
% Process Open Ephys file
% For LFP...
options.NewOESampleRate = options.mapSampleRate;
[OE,EyeDat,options] = parseRawToCollated(thisAnimal,thisSession,thisFileNum,options);
ttltimestamps = union(OE.OEttl.upPhases, OE.OEttl.downPhases);

%%%%%
% Process bonvision file containing stimulus information
% Read in the bonsai stimulus log
thisBonsaiFile = fullfile(options.BonsaiPath, [ thisAnimal,'_SparseNoise_Log_', '*']);
folder_dir = dir(thisBonsaiFile);
if isempty(folder_dir) % Not always with animal prefix
    thisBonsaiFile = fullfile(options.BonsaiPath, [ 'SparseNoise_Log_', '*']);
    folder_dir = dir(thisBonsaiFile);
end
if isempty(folder_dir) % Not always the same file root as above
    thisBonsaiFile = fullfile(options.BonsaiPath, [ thisAnimal,'_quaddata', '*']); % May need to be able to include other options here eg including filenumber
    folder_dir = dir(thisBonsaiFile);
end
folder_name = folder_dir.name;
thisBonsaiFileName = fullfile(options.BonsaiPath,folder_name);
fileID=fopen(thisBonsaiFileName);
thisBinFile=fread(fileID);
fclose(fileID);

% Translate stimulus into -1:1 scale
stim_matrix = zeros(1,length(thisBinFile));
stim_matrix(thisBinFile==0)=-1;
stim_matrix(thisBinFile==255)=1;
stim_matrix(thisBinFile==128)=0;
% Make a NxM grid from the stimulus log
stim_matrix = reshape(stim_matrix, [options.grid_size(1), options.grid_size(2), length(thisBinFile)/options.grid_size(1)/options.grid_size(2)]);
stim_matrix = stim_matrix(:,:,1:end-1); % The last upswing should be ignored

% Compare 
if length(ttltimestamps)<size(stim_matrix,3)
    warning('Number of ttls (%d) LESS than number of stimulus instances (%d).',length(ttltimestamps),size(stim_matrix,3));
    stim_matrix = stim_matrix(:,:,1:length(ttltimestamps));
elseif length(ttltimestamps)>size(stim_matrix,3)
    warning('Number of ttls (%d) MORE than number of stimulus instances (%d).',length(ttltimestamps),size(stim_matrix,3));
    ttltimestamps = ttltimestamps(1:size(stim_matrix,3));
end

% Get the lfp and wheel signal for each stimulus presentation
% - samples to keep after ttl pulse (equals to stimulus duration)
samples_to_keep =  options.stim_dur*options.mapSampleRate;
lfp_data = nan(length(ttltimestamps),samples_to_keep);
wheel_data = lfp_data;
for i=1:length(ttltimestamps)
    lfp_data(i,:) = OE.OELFPData(1,ttltimestamps(i):ttltimestamps(i)+samples_to_keep-1);
    wheel_data(i,:) = OE.OEwheelData.speedSmoothed(1,ttltimestamps(i):ttltimestamps(i)+samples_to_keep-1);
end
wheel_data = nansum(wheel_data,2);

% Get EYE Area data (average over trial)
if ~isempty(EyeDat) && isfield(EyeDat,'upPhases')
    eyettltimestamps = union(EyeDat.upPhases, EyeDat.downPhases);
    if length(ttltimestamps)~=length(eyettltimestamps)
        warning('Different  number of LFP ttls and EYE ttls - ignoring EYE')
        eye_data = NaN(length(ttltimestamps),1);
    else
        eyesamples_to_keep = fix(options.stim_dur*EyeDat.SampleRate);
        eye_data = nan(length(eyettltimestamps),eyesamples_to_keep);
        eye_valid = eye_data;
        for i=1:length(eyettltimestamps)
            eye_data(i,:) = EyeDat.data.EyeArea_mm2(eyettltimestamps(i):eyettltimestamps(i)+eyesamples_to_keep-1)';
            eye_valid(i,:) = EyeDat.data.valid(eyettltimestamps(i):eyettltimestamps(i)+eyesamples_to_keep-1)';
        end
        eye_valid(eye_valid == 0) = NaN;
        eye_data = eye_data.*eye_valid;
        eye_data  = nansum(eye_data,2);
    end
else
    eye_data = [];
end
%%%%%
% Call mapping progamme
options.figName = sprintf('%s :: %s (file %01d)',    thisAnimal, thisSession,  thisFileNum);
initMap = sparseNoiseAnalysis(stim_matrix,lfp_data,wheel_data,eye_data,options);
