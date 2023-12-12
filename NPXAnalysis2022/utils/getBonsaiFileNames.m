% Find relevant files in bonsai directory for this stimulus
%
% History
%  24/2/2022 SGS Wrote it - ported from batch_processing.m
%  27/2/2022 SGS updated to make sure multiple files are aligned in list
%               correctly, by date in filename
%  03/4/2022 SGS updated to include PDAnalog files
%  07/4/2022 SGS Updated to accommodate new log style
% 10/01/2023 MT include sparse noise full screen version

function [BehaviourDataFiles,EyeDataFiles,TrialParamsFiles,PDFiles,PDSyncFiles] = getBonsaiFileNames(StimulusName,BONSAI_DATAPATH)
% Find relevant fles in bonsai folder
bonsai_files = dir(fullfile(BONSAI_DATAPATH,['*', StimulusName, '*']));
bonsai_files_names = {bonsai_files.name};
% Switch depending on whether file types are new (>7th April 2022) or old (<7th April 2022)
BehaviourDataFiles1 = bonsai_files_names(contains(bonsai_files_names,'BehaviourData'));
BehaviourDataFiles2 = bonsai_files_names(contains(bonsai_files_names,'WheelLog'));
if ~isempty(BehaviourDataFiles1) & isempty(BehaviourDataFiles2)
    fileStyle = 'old';
elseif isempty(BehaviourDataFiles1) & ~isempty(BehaviourDataFiles2)
    fileStyle = 'new';
else
    error('Not sure which type of files to find')
end
switch(fileStyle)
    case 'old'
        % Identify file types
        BehaviourDataFiles = bonsai_files_names(contains(bonsai_files_names,'BehaviourData'));
        EyeDataFiles = bonsai_files_names(contains(bonsai_files_names,'EyeTrackingData'));
        TrialParamsFiles = bonsai_files_names(logical([contains(bonsai_files_names,'TrialParams') + contains(bonsai_files_names,'SparseNoise_Log')]));
        
        % Check to see if equivalent number of files present
        if length(BehaviourDataFiles) ~= length(EyeDataFiles)  || length(BehaviourDataFiles) ~= length(TrialParamsFiles)  || length(EyeDataFiles) ~= length(TrialParamsFiles)
            error('Found unequal number of behavioural, eye and trial files for target')
        end
        
        % Get timestamps
        BehaviourDataFilesTimes = cellfun(@(x) getTimeStampsFromFileName(x),BehaviourDataFiles);
        EyeDataFilesTimes = cellfun(@(x) getTimeStampsFromFileName(x),EyeDataFiles);
        TrialParamsFilesTimes = cellfun(@(x) getTimeStampsFromFileName(x),TrialParamsFiles);
        
        % Sort in ascending time
        BehaviourDataFiles = sortData(BehaviourDataFiles,BehaviourDataFilesTimes);
        EyeDataFiles = sortData(EyeDataFiles,EyeDataFilesTimes);
        TrialParamsFiles = sortData(TrialParamsFiles,TrialParamsFilesTimes);
        
        if nargout > 3
            PDFiles = bonsai_files_names(contains(bonsai_files_names,'PDAnalog'));
            PDSyncFiles = bonsai_files_names(contains(bonsai_files_names,'SyncPulseAnalog'));
            PDFilesTimes = cellfun(@(x) getTimeStampsFromFileName(x),PDFiles);
            PDSyncFilesTimes = cellfun(@(x) getTimeStampsFromFileName(x),PDSyncFiles);
            PDFiles = sortData(PDFiles,PDFilesTimes);
            PDSyncFiles = sortData(PDSyncFiles,PDSyncFilesTimes);
        end
        
    case 'new'
        % Identify file types
        BehaviourDataFiles = bonsai_files_names(contains(bonsai_files_names,'WheelLog'));
        EyeDataFiles = bonsai_files_names(contains(bonsai_files_names,'EyeLog'));
        TrialParamsFiles = bonsai_files_names(contains(bonsai_files_names,'Trial'));
        TrialParamsFiles2 = bonsai_files_names(contains(bonsai_files_names,'SparseNoise_Log'));
        TrialParamsFiles3 = bonsai_files_names(contains(bonsai_files_names,'SparseNoise_fullscreen_Log'));
        
        if ~isempty(TrialParamsFiles2)
            TrialParamsFiles = TrialParamsFiles2;
        elseif ~isempty(TrialParamsFiles3) %If fullscreen,choose fullscreen sparsenoise
            TrialParamsFiles = TrialParamsFiles3;

        end
        
        % Check to see if equivalent number of files present
        if length(BehaviourDataFiles) ~= length(EyeDataFiles)  || length(BehaviourDataFiles) ~= length(TrialParamsFiles)  || length(EyeDataFiles) ~= length(TrialParamsFiles)
            error('Found unequal number of behavioural, eye and trial files for target')
        end
        
        % Get timestamps
        BehaviourDataFilesTimes = cellfun(@(x) getTimeStampsFromFileName(x),BehaviourDataFiles);
        EyeDataFilesTimes = cellfun(@(x) getTimeStampsFromFileName(x),EyeDataFiles);
        TrialParamsFilesTimes = cellfun(@(x) getTimeStampsFromFileName(x),TrialParamsFiles);
        
        % Sort in ascending time
        BehaviourDataFiles = sortData(BehaviourDataFiles,BehaviourDataFilesTimes);
        EyeDataFiles = sortData(EyeDataFiles,EyeDataFilesTimes);
        TrialParamsFiles = sortData(TrialParamsFiles,TrialParamsFilesTimes);
        
        if nargout > 3
            PDFiles = bonsai_files_names(contains(bonsai_files_names,'PDAsync'));
            PDFilesTimes = cellfun(@(x) getTimeStampsFromFileName(x),PDFiles);
            PDFiles = sortData(PDFiles,PDFilesTimes);
            PDSyncFiles = [];
        end
end


end


function sortedData = sortData(data,sorter)
[~,ind] = sort(sorter,'ascend');
sortedData = data(ind);
end
