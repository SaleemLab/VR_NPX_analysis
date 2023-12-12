% Programe to illustrate how to find data files and load in data in
% template format

%%%%%%%%%
% If all goes well, you should only need to vary the next few lines
% Example:
SUBJECT = 'M22030';%'M21200';
SESSION = '20220617';%'20211206';
StimulusName = 'PosAdapt';
options.paradigm = 'ad'; % That's required in the importAndAlignBonsaiLogs 
options.importMode = 'KS'; % LF or MUA or KS
options.BinWidth = 1/60; % resolution (in s) of output resps (e.g. 1/60)
options.AnalysisTimeWindow = [0 0.1];% two-element vector specifying time window around stim-on (e.g. [-0.25 1.25])
options.ks_unitType = 'good'; % 'mua', 'good' or '' if using KS
%%%%%%%%%

%%%%%%%%%
% Step 1: Set up paths
if ismac
    ROOTPATH = '/Users/s.solomon/Filestore/Research2/ibn-vision';
else
    ROOTPATH = 'X:\ibn-vision';
end
DATAPATH = fullfile(ROOTPATH,'DATA','SUBJECTS');
EPHYS_DATAPATH = fullfile(DATAPATH,SUBJECT,'ephys',SESSION);
BONSAI_DATAPATH = fullfile(DATAPATH,SUBJECT,'bonsai',SESSION);

%%%%%%%%%
% Step 2: Find csv files associated with desired stimulus
fileTable = getBoundGAndBonsaiFileNames(SUBJECT,SESSION,DATAPATH);
theseFiles = find(contains(fileTable.StimulusName,StimulusName));
BehaviourDataFiles = fileTable.WheelLog(theseFiles);
EyeDataFiles = fileTable.EyeLog(theseFiles);
TrialParamsFiles = fileTable.Trial(theseFiles);
PDFiles = fileTable.PD(theseFiles);
gNums = fileTable.gNumber(theseFiles);
if ~isempty(PDFiles)
    options.PD_FLAG = 1;    % Have saved the Photodiode output   
else
    options.PD_FLAG = 0;
end
% If Kilosort 
switch(options.importMode)
    case 'KS'
        options.KS_DATAPATH = fullfile(EPHYS_DATAPATH,'kilosort');
        options.KS_CATGT_FNAME = ['CatGT_', SUBJECT, '_', SESSION, '.log'];
end
%%%%%%%%%
% Step 3: Load data
for thisFile = 1:length(BehaviourDataFiles)
    % Set bonsai data paths
    options.PERIPHERALS_DATAPATH = fullfile(BONSAI_DATAPATH, BehaviourDataFiles{thisFile});
    options.EYEDATA_DATAPATH = fullfile(BONSAI_DATAPATH, EyeDataFiles{thisFile});
    options.TRIALDATA_DATAPATH = fullfile(BONSAI_DATAPATH,TrialParamsFiles{thisFile});
    if options.PD_FLAG
        options.PHOTODIODE_DATAPATH = fullfile(BONSAI_DATAPATH,PDFiles{thisFile});
    end

    % Set ephys data path (temporary, we want to automate this sometime)
    folderName = findGFolder(EPHYS_DATAPATH,gNums(thisFile));
    T_EPHYS_DATAPATH = fullfile(EPHYS_DATAPATH,folderName,[folderName,'_imec0']);
    options.EPHYS_DATAPATH = T_EPHYS_DATAPATH; % full path to ephys recording
    options.gFileNum = gNums(thisFile);
    [tresps,totherData,tstimData,~,~,~,ttimeVector] = extractAndCollateNPData(options);
    if thisFile == 1
        resps = tresps;
        otherData = totherData;
    else
        resps = cat(3,resps,tresps);
        otherData = cat(3,otherData,totherData);
    end
end
