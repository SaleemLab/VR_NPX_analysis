% Programme to return bound g-files and bonsai log files for given
% recording

% FR Edited on 16/05/2022 to deal with error on loading data from
% recordings from before we started saving AnalogPD

function fileTable = getBoundGAndBonsaiFileNames(SUBJECT,SESSION,DATAPATH)

% First, check to see if there is a 'processed' folder in the subject root
BONSAI_DATAPATH = fullfile(DATAPATH,SUBJECT,'bonsai',SESSION);
PROCESSED_DATAPATH = fullfile(DATAPATH,SUBJECT,'processed',SESSION);

% Check if 'processed' folder exists
if ~exist(fullfile(DATAPATH,SUBJECT,'processed'),'dir')
    [SUCCESS,MESSAGE,MESSAGEID] = mkdir(fullfile(DATAPATH,SUBJECT),'processed');
end
% Check if SESSION folder exists
if ~exist(PROCESSED_DATAPATH,'dir')
    [SUCCESS,MESSAGE,MESSAGEID] = mkdir(PROCESSED_DATAPATH);
end
% Define summary file name
SummaryFileName = sprintf('%s_%s_SummaryFile.csv',SUBJECT,SESSION);
if exist(fullfile(PROCESSED_DATAPATH,SummaryFileName),'file')
    % If there is log file ... read table in
    fileTable = readtable(fullfile(PROCESSED_DATAPATH,SummaryFileName));
else
    % If there is not a log file ... build it and save
    td1 = dir(fullfile(BONSAI_DATAPATH,'*WheelLog*.csv')); % >7th April 2022
    td2 = dir(fullfile(BONSAI_DATAPATH,'*BehaviourData*.csv')); % For < 7th April 2022
    td = cat(1,td1,td2);
    tx = {}; [tx{1:length(td)}] = deal(td.name);
    allStimulusNames = strtok(tx,'_'); % Find the stimuli in the folder
    allStimulusFilesTimes = cellfun(@(x) getTimeStampsFromFileName(x),tx);
    [~,ind] = sort(allStimulusFilesTimes,'ascend');
    allStimulusNames = allStimulusNames(ind);
    
    % For each stimulus
    % Make N by 2 matrix of fieldname + value type
    variable_names_types = [["StimulusName", "string"]; ...
        ["gNumber", "double"]; ...
        ["Trial", "string"]; ...
        ["WheelLog", "string"]; ...
        ["EyeLog", "string"]; ...
        ["PD", "string"]; ...
        ["Comments", "string"]];
    % Make table using fieldnames & value types from above
    fileTable = table('Size',[0,size(variable_names_types,1)],...
        'VariableNames', variable_names_types(:,1),...
        'VariableTypes', variable_names_types(:,2));
    
    warning off
    thisRow = 0;
    for thisStimulus = 1:length(allStimulusNames)
        % Retrieve this stimulus
        StimulusName = allStimulusNames{thisStimulus};
        % Find relevant log files and construct
        [BehaviourDataFiles,EyeDataFiles,TrialParamsFiles,PDFiles] = getBonsaiFileNames(StimulusName,BONSAI_DATAPATH);
        % Build table
        for thisFile = 1:length(BehaviourDataFiles)
            thisRow = thisRow+1;
            fileTable.StimulusName(thisRow) = StimulusName;
            fileTable.gNumber(thisRow) = NaN;
            fileTable.WheelLog(thisRow) = BehaviourDataFiles{thisFile};
            fileTable.EyeLog(thisRow) = EyeDataFiles{thisFile};
            fileTable.Trial(thisRow) = TrialParamsFiles{thisFile};
            if ~isempty(PDFiles)                
                fileTable.PD(thisRow) = PDFiles{thisFile};
            else
                fileTable.PD(thisRow) = 'Missing';
            end
        end
    end
    warning on
    % Write table to file
    writetable(fileTable,fullfile(PROCESSED_DATAPATH,SummaryFileName),'WriteVariableNames',true,'Delimiter',',')
    
    fprintf('\nPlease open the following file and enter g-numbers to proceed:\n\t%s\n',fullfile(PROCESSED_DATAPATH,SummaryFileName))
end