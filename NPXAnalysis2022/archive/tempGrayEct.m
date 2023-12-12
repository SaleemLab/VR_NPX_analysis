% Note that:
%       within each SUBJECT there are 'bonsai' and 'ephys' folders
%       within each of those there are date folders (for recording sessions)
%       in ephys/Date/ folder there will be N folders, for each recording,
%               indicated by SUBJECT_DATE_gN
%           in which there is a folder '...gN_imec0'
%               in which is .bin & .meta files for both 'ap' (30kHz?) and 'lf' (2500 Hz) 
%               e.g. M22003_20220221_g0_t0.imec0.lf.bin and M22003_20220221_g0_t0.imec0.ap.bin
%           There will, if sorted, also be a 'kilosort' folder, which will
%           have the results of sorting on the concatenated data files
%               At the minimum this will require the following file types:
%
%
%       in bonsai/Date/ folder there will be 3 .csv files/ recording
%                 e.g
%                 'DriftingGratings__EyeTrackingData_Session_2022-02-21T15_13_35'
%                     There are 9 columns, which are ....
%                           Centroid.X	Centroid.Y	Orientation
%                           MajorAxisLength	MinorAxisLength	Area 
%                           Time pd/async(?) FrameNumber
%                 'DriftingGratings_BehaviourData2022-02-21T15_13_35' 
%                       There are 5 columns 
%                           Wheel	Sync	EyeFrameCount	Photodiode	Time
%                 'DriftingGratings_TrialParams2022-02-21T15_13_52'
%                       There are 2 columns (for a drifting grating stimulus anyway)
%                           DriftingGratingOrientation, Time    
%       
if ismac
    ROOTPATH = '/Users/s.solomon/Filestore/Research2/ibn-vision';
else
    ROOTPATH = 'X:\ibn-vision';
end

%%%%%
% For Grey
%
% SUBJECT = 'M22002';
% SESSION = '20220216';
% gs = [2];
%
% SUBJECT = 'M22003';
% SESSION = '20220221';
% gs = [5];   % Note that for some reason this looks like it is g=6 in the session, but best correlation is with g = 5

% SUBJECT = 'M22003';
% SESSION = '20220222';
% gs = [5];       % Note that appears as g = 4 but is actually g = 5
%
% SUBJECT = 'M22003';
% SESSION = '20220223';
% gs = [5];
% 
% SUBJECT = 'M22004';
% SESSION = '20220323';
% gs = [3];
% 
SUBJECT = 'M22004';
SESSION = '20220324';
gs = [5];
%
% SUBJECT = 'M22004';
% SESSION = '20220325';
% gs = [5];
%%%%%%

% Example:
StimulusName = 'Grey';
NP2_h0_STRING = ''; % '' if NP1, 'h0_' if NP2
options.importMode = 'KS'; % 'LF' or 'MUA' or 'KS'
options.BinWidth = 1; % resolution (in s) of output resps (e.g. 1/60)
options.AnalysisTimeWindow = [];%[-0.25 1.5];% two-element vector specifying time window around stim-on (e.g. [-0.25 1.25])
options.ks_unitType = ''; % 'mua', 'good' or ''

% Step 1: get correct paths
DATAPATH = fullfile(ROOTPATH,'DATA','SUBJECTS');
EPHYS_DATAPATH = fullfile(DATAPATH,SUBJECT,'ephys',SESSION);
BONSAI_DATAPATH = fullfile(DATAPATH,SUBJECT,'bonsai',SESSION);
% If Kilosort 
switch(options.importMode)
    case 'KS'
        options.KS_DATAPATH = fullfile(EPHYS_DATAPATH,'kilosort3');
end

% Step 2: Find csv files associated with desired stimulus
[BehaviourDataFiles,EyeDataFiles,TrialParamsFiles] = getBonsaiFileNames(StimulusName,BONSAI_DATAPATH);

% Step 3: We want to match these files with the relevant ePhys file
% Not done yet

% Step 4: Call parsing routine
for thisFile = 1%:length(BehaviourDataFiles)
    % Set bonsai data paths
    PERIPHERALS_DATAPATH = fullfile(BONSAI_DATAPATH, BehaviourDataFiles{thisFile});
    EYEDATA_DATAPATH = fullfile(BONSAI_DATAPATH, EyeDataFiles{thisFile});
    TRIALDATA_DATAPATH = fullfile(BONSAI_DATAPATH,TrialParamsFiles{thisFile});
    options.PERIPHERALS_DATAPATH = PERIPHERALS_DATAPATH; % full path to BehaviourData bonsai logfile
    options.EYEDATA_DATAPATH = EYEDATA_DATAPATH; % full path to EyeTrackingData bonsai logfile
    options.TRIALDATA_DATAPATH = TRIALDATA_DATAPATH; % full path to TrialData bonsai logfile
    
    % Set ephys data path (temporary, we want to automate this sometime)
    options.gFileNum = gs(thisFile);
    T_EPHYS_DATAPATH = fullfile(EPHYS_DATAPATH,sprintf('%s_%s_%sg%01d',SUBJECT,SESSION,NP2_h0_STRING,gs(thisFile)),sprintf('%s_%s_%sg%01d_imec0',SUBJECT,SESSION,NP2_h0_STRING,gs(thisFile)));
    T_EPHYS_DATAPATH = fullfile(T_EPHYS_DATAPATH,sprintf('%s_%s_%sg%01d_t0.imec0',SUBJECT,SESSION,NP2_h0_STRING,gs(thisFile)));
    options.EPHYS_DATAPATH = T_EPHYS_DATAPATH; % full path to ephys recording
    
    [tresps,totherData,tstimData,~,~,~,ttimeVector,options] = extractAndCollateNPData(options);
    if thisFile == 1
        resps = tresps;
        otherData = totherData;
    else
        resps = cat(3,resps,tresps);
        otherData = cat(3,otherData,totherData);
    end
end


% Temporary example of how to calculate correlations during grey
% First plot the data as a function of time
% normresps = tresps./repmat(nanmean(tresps,2),1,size(tresps,2));
% Can contrain to epochs of stationary 
tp = find(abs(totherData(1,:))<0.001);
% tp = 1:size(totherData,2);
ttresps= tresps(:,tp);
% Zscore data (need to transpose because z-score works across columns)
ttresps = zscore(ttresps')';
% Combine units in same channel
newttimeVector = ttimeVector(tp);
ttotherData = totherData(:,tp);
newttresps = zeros(384,size(ttresps,2));
for thisChannel = 1:384
    tp2 = find(options.peakChannel >= thisChannel-1 & options.peakChannel <= thisChannel+1);
    newttresps(thisChannel,:) = sum(ttresps(tp2,:),1);
end

figure;
subplot(2,1,1);
h1 = imagesc(newttimeVector,1:size(newttresps,1),newttresps,[0 3]);
set(gca,'YDir','normal')
subplot(4,1,3);
h2 = plotyy(newttimeVector,ttotherData(1,:),newttimeVector,ttotherData(2,:)); hold on
%plot(ttimeVector(tp),0,'ro')
linkaxes([h1 h2(1) h2(2)],'x')
set(gca,'XLim',[min(ttimeVector) max(ttimeVector)])

% Now calculate correlation matrix!
% - as a function of time
covresps = cov(newttresps);
figure; 
imagesc(covresps,[-1 1])
set(gca,'YDir','normal')
axis square

% - as a function of unit
covresps = cov(newttresps');
figure; 
imagesc(covresps,[-3 3])
set(gca,'YDir','normal')
axis square
