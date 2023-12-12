% Stefano LFP sparse noise maps
% SGS 24th August 2020
% Following specifics adjusted for SparseNoise in M19165 Session 191219

thisAnimal = 'M19165';
thisSession = '191219';
thisFileNum = 6;

% Some defaults
options.photo_br = 0; % 0 (pulse width modulated) or 1 (constant)
options.sync_input = 'adc';
options.sync_channel = 4; % photodiode channel
options.OE_pD_thresholds = [0.08 0.2]; % lower and upper thresholds - lower currently not working
options.wheel_input = 'adc';
options.wheel_channels = 5:6; 
options.includeWheel = 1;
options.UseOneSignal = 0;
options.chanNumbersToInclude = 2;
options.EyeSR = 31;
options.logErrorsFlag = 0;

options.plotflag = 1;

options.stim_dur = 0.1;
options.grid_size = [8 8]; % 10x10 x y
options.mapSampleRate = 120; % Hz

options.mapMethod = 'mean';

%%%%%%%%%
% Set up some paths
if ispc
    if exist('X:\ibn-vision','dir')==7
        options.serverPath = 'X:\ibn-vision\';
    else
        options.serverPath = 'X:\';
    end
elseif ismac
    options.serverPath = '/Users/s.solomon/Filestore/Research2/ibn-vision/';
end
options.BonsaiPath = fullfile(options.serverPath,'DATA','SUBJECTS',upper(thisAnimal),'VRBehaviour',upper(thisSession),'SparseNoise');
options.conToppath = fullfile(options.serverPath,'DATA','SUBJECTS',upper(thisAnimal),'ePhys',upper(thisSession));


[initMap, lfp_data,wheel_data,eye_data] = sparseNoiseLFPAnalyse(thisAnimal,thisSession,thisFileNum,options);