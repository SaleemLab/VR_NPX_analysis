%%%%
% Define some paths
if ispc
    if logical(exist('X:\ibn-vision','file'))
        options.ibnvisionPath = 'X:\ibn-vision';
        options.figurePosVec = [685,42,681,642];
%         options.figurePosVec = [1921,71,1680,933];
    else
        options.ibnvisionPath = 'X:\';
        options.figurePosVec  = [962,42,958,954];
    end
elseif ismac
%     options.ibnvisionPath = '/Users/s.solomon/Filestore/Research2/ibn-vision';
    options.ibnvisionPath ='/Users/s.solomon/Desktop/TEMP';
end
options.dataPath = fullfile(options.ibnvisionPath,'DATA','SUBJECTS');

% Define the bleach correction method
options.bleachCorrMethod = 'entire';
options.plotBleachCorrFlag = 0;
% Define some other options
options.plotWheelFlag = 0;
options.plotPdDistFlag = 0;
options.plotDfFlag = 0;
options.trialType = 'sparse_noise';
options.pdMethod = 'GrayHasNoWhite';
% options.fig1Handle = fig1Handle;
options.figs2plot = [1:6];
options.plotAllAnimalsFlag = 0;
options.plotAllDatesFlag = 1;
options.correctAcqDelayFlag = 1;
options.saveFiguresFlag = 0;
options.SkipPSTHsFlag = 0;
options.analyseEWPflag= 0;
options.fitGaussesFlag = 1;
options.analyseRunningFlag = 0;

% ANIMALS_OI = {'M19019','M19053','M19054','M19055','M19126','M19127','M19128','M19129','M19130','M19131'};
ANIMALS_OI = {'M19127'};


% read the data table
METADATAFILE = fullfile(options.ibnvisionPath,'DATA','DATABASES',...
    'DB-Miniscope','TW','TWLookUpTables_Miniscope_Vis_Stim.xlsx');
META_TABLE = readtable(METADATAFILE,'ReadRowNames', 0,'ReadVariableNames', 1,'Sheet','Everything');

% Get sparse noise data
SN_SESSIONS_INDEX = contains(META_TABLE.SPARSE_NOISE,'Y');

ANIMAL_OI_INDEX = zeros(size(META_TABLE,1));
DATE_OI_INDEX = ANIMAL_OI_INDEX;

if options.plotAllAnimalsFlag
    ANIMAL_OI_INDEX = ones(size(META_TABLE,1));
else 
    for a = 1:length(ANIMALS_OI)
        ANIMAL_OI_INDEX = ANIMAL_OI_INDEX + contains(META_TABLE.ANIMAL,ANIMALS_OI{a});
    end
end

if options.plotAllDatesFlag
    DATE_OI_INDEX = ANIMAL_OI_INDEX;
else
    DATE_OI_INDEX = META_TABLE.DATE == 191130;
end

SN_SESSIONS = find(all([SN_SESSIONS_INDEX,ANIMAL_OI_INDEX,DATE_OI_INDEX],2));

close all
for currSessCounter = 2:length(SN_SESSIONS)
    tic
    thisData = META_TABLE(SN_SESSIONS(currSessCounter),:);
    doiFiles = fullfile(options.dataPath,thisData.ANIMAL{1},'MiniscopeData',num2str(thisData.DATE(1)));
    options.rigID = thisData.RIG{1};
    options.currAnimal = str2num(thisData.ANIMAL{1}(2:end));
    options.currDate = thisData.DATE(1);
    
    if options.saveFiguresFlag
        folderName = fullfile(options.ibnvisionPath,'DATA\PROJECTS\SC_Miniscope_Imaging\Sparse_Noise_Experiments',...
            num2str(options.currAnimal), num2str(options.currDate));
        if exist(folderName)
            continue
        end
    end
    
    % First get data
    outdat = readInHeadRestrainedDataNew(doiFiles,1,options);
    % save([doiFile,'\SLAprocessedData'],'outSLA')
    %     drawnow
    if isempty(outdat)
        continue
    end
    close(figure(1))   
    
    if options.saveFiguresFlag
        saveAllFigs(folderName)
    end
    
    currSessCounter/length(SN_SESSIONS)
    toc
    close all
    animalsList(currSessCounter) = str2num(thisData.ANIMAL{1}(2:end));
    datesList(currSessCounter) = thisData.DATE(1);
    outdatStore(currSessCounter) = outdat;
end
