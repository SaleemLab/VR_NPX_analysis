% Programme to estiamte receptive fields from sparse noise in NP1 / NP2
% data
% Dependencies: NPXAnalysis2022;
%               Visual-response-analysis/Stimuli/sparsenoise

% First load the data
if ismac
    ROOTPATH = '/Users/s.solomon/Filestore/Research2/ibn-vision';
else
    ROOTPATH = 'X:\ibn-vision';
end

% Specify file
SUBJECT = 'M22008';%'M21200';
SESSION = '20220408';%'20211206';
StimulusName = 'SparseNoise';
gs = [0];
nChannelsToBin = 24;    % Bin firing rates across channels for RF maps
channelRange = [1 384]; % Look at these channels

% Some defaults
options.importMode = 'KS'; % LF or MUA or KS
options.BinWidth = 1/60; % resolution (in s) of output resps (e.g. 1/60)
options.AnalysisTimeWindow = [0 0.5];% two-element vector specifying time window around stim-on (e.g. [-0.25 1.25])
options.ks_unitType = 'good'; % 'mua', 'good' or ''

% Get correct paths
DATAPATH = fullfile(ROOTPATH,'DATA','SUBJECTS');
EPHYS_DATAPATH = fullfile(DATAPATH,SUBJECT,'ephys',SESSION);
BONSAI_DATAPATH = fullfile(DATAPATH,SUBJECT,'bonsai',SESSION);
% If Kilosort 
switch(options.importMode)
    case 'KS'
        options.KS_DATAPATH = fullfile(EPHYS_DATAPATH,'kilosort');
end

% Step 2: Find csv files associated with desired stimulus
[BehaviourDataFiles,EyeDataFiles,TrialParamsFiles] = getBonsaiFileNames(StimulusName,BONSAI_DATAPATH);

% Set bonsai data paths
PERIPHERALS_DATAPATH = fullfile(BONSAI_DATAPATH, BehaviourDataFiles{1});
EYEDATA_DATAPATH = fullfile(BONSAI_DATAPATH, EyeDataFiles{1});
TRIALDATA_DATAPATH = fullfile(BONSAI_DATAPATH,TrialParamsFiles{1});
options.PERIPHERALS_DATAPATH = PERIPHERALS_DATAPATH; % full path to BehaviourData bonsai logfile
options.EYEDATA_DATAPATH = EYEDATA_DATAPATH; % full path to EyeTrackingData bonsai logfile
options.TRIALDATA_DATAPATH = TRIALDATA_DATAPATH; % full path to TrialData bonsai logfile

% Set ephys data path
options.gFileNum = gs(1);
folderName = findGFolder(EPHYS_DATAPATH,options.gFileNum);
T_EPHYS_DATAPATH = fullfile(EPHYS_DATAPATH,folderName,[folderName,'_imec0']);
options.EPHYS_DATAPATH = T_EPHYS_DATAPATH; % full path to ephys recording
options.KS_CATGT_FNAME = 'CatGT_M22008_20220408.log';

% Extract data
[resps,otherData,stimData,~,~,~,timeVector,options] = extractAndCollateNPData(options);

%%
% Calculate sparsenoise RFs
stim_matrix = cat(3,stimData.stim_matrix{:}); % N rows x M cols x nFrames
sn_options.grid_size = [size(stim_matrix,1) size(stim_matrix,2)];
sn_options.mapSampleRate = 60; % Hz
sn_options.mapsToShow = {'linear','black','white','contrast'};
sn_options.mapMethod = 'fitlm'; % fitlm mean
sn_options.framesToShow = [1 4 6 8 10 12];  % at 60 Hz would be 8 ms, 40, 72 etc
sn_options.plotflag = 1;

% Extract wheel and eye data
wheel_data = squeeze(nanmean(otherData(1,:,:),2))';
eye_data = squeeze(nanmean(otherData(2,:,:),2))';
% Choose a channel...or two
channelBins = channelRange(1):nChannelsToBin:channelRange(2);
for theseChannels = 1:length(channelBins)-1
    % Get average across channels (zscore?)
    switch(options.importMode)
        case 'KS'
            tp = options.peakChannel >= channelBins(theseChannels) & options.peakChannel < channelBins(theseChannels+1);
        otherwise
            tp = channelBins(theseChannels):channelBins(theseChannels+1);            
    end
    lfp_data = mean(resps(tp,:,:),1); % returns 1 x N time bins x nFrames
    lfp_data = squeeze(lfp_data)';   % sparsenoise needs frames x N time bins

    sn_options.figName = sprintf('%s :: %s (file %01d) Channels %01d:%01d',    SUBJECT, SESSION,  gs(1), channelBins(theseChannels),channelBins(theseChannels+1));
    initMap = sparseNoiseAnalysis(stim_matrix,lfp_data,wheel_data,eye_data,sn_options); % check to make sure this is the correct one - should be in Visua;-response-analysis/Stimuli/sparsenoise
    drawnow;
end
