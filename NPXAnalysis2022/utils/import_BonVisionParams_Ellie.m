% function [datatable,StimulusName] = import_BonVisionParams(filepath,varargin)
%   Function to collate different import styles for different types of
%   BonVision stimulus parameter logging
% 
% Inputs: 
%   filepath (full path to Bonsai log file). If this includes a
%       'StimulusName' from the list below, that is all you need to do
%
%   Variable input arguments in standard Matlab pairwise notation for additional arguments 
%       e.g. ...,'startRow',10,'endRow',20) - this will result in declaring
%       a variable with name startRow and value 10, and a variable with name endRow and value 20
%    
%       One of these arguments can be 'StimulusType' if the csv filename isn't
%       explicit or needs to be overwritten. e.g. include the pair ...,'StimulusType','StaticGrating') 
%
% Currently supports 'StimulusName' of: 
%       [Standard] which includes 'SpeedTextures', 'Asym', 'SpeedDiscrim', 'Textures2D', 'SpeedTuning','Speed2D','FastSpeed'
%       [StaticGratings] for 'StaticGratings' 
%
% Example call:
%   datatable = import_BonVisionParams( ...
%                   '/Users/s.solomon/Filestore/Research2/ibn-vision/DATA/SUBJECTS/M22003/bonsai/20220223/StaticGratings_TrialParams2022-02-23T14_36_17.csv', ...
%                   'StimulusName','StaticGratings');
%
% History: SGS Wrote it, collating the various functions included below (in
%   each case the original function is in utils and is named 'import_Bonsai*' while the realted function
%   below is 'import_BonVision*'
%          FR Added DriftingTF as a case on 10/06/2022
%
function [datatable,StimulusName] = import_BonVisionParams_Ellie(filepath,varargin)

% Check to see if variable exists
if ~exist('filepath','var') || isempty(filepath)
    error('No filepath specified')
end

% Parse for one of the following stimulus names (you can add to the list here, and to the switch case below)
specialFileTypes = {'StaticGratings','DriftingGratings','SparseNoise_fullscreen','StaticTF','Grey','GratingsPhase', ...
        'OriAdapt','PosAdapt','DriftingTF','Checkerboard', ... 
        'SparseNoise','OP_Tuning','Direction_Tuning', ...
        'TRAIN', 'DCBA', 'OMIT', 'E_CD', 'ADCD', 'lowcontB', 'lowcontDsubbingB', 'lowcontTRAIN', ... 
        'GAVNIK_ABCD', 'GAVNIK_A_CD', 'GAVNIK_E_CD', 'GAVNIK DCBA', 'GAVNIK_D___', 'GAVNIK_A___', 'GAVNIK_ABCD_1',...
        'GAVNIK_ABCD_2', 'GAVNIK_D_CD', 'D___', 'A___', 'TRAIN_1', 'TRAIN_2', 'D_CD',...
        'A_50ms', 'A_500ms', 'A_200ms', 'A_300ms', 'OMIT50grey',...
        'A_1000ms', 'A_1000ms_1', 'A_1000ms_2', 'A_1000ms_25pc', 'A_1000ms_50pc', 'A_1000ms_75pc', 'E_1000ms', 'TRAIN250',...
        'GAVNIK200_ABCD', 'GAVNIK200_A_CD', 'GAVNIK200_E_CD', 'GAVNIK250_ABCD', 'GAVNIK250_ABCD_1', 'GAVNIK250_ABCD_2',...
        'GAVNIK250_A_CD', 'GAVNIK250_E_CD', 'GAVNIK250 DCBA', 'F_150ms', 'F_150ms_1', 'F_150ms_2', 'F_1000ms'};

% See if there is a match within the file name that indicates that this is
% from the special list
StimulusTypeMatcher = zeros(1,length(specialFileTypes));
[~, filename, ~] = fileparts(filepath);

% Extract tag between subject ID and _StimulusParams
stimPattern = regexp(filename, '^[^_]+_(.*)_StimulusParams', 'tokens');
if ~isempty(stimPattern)
    stimulusTag = strtrim(stimPattern{1}{1});
else
    stimulusTag = '';
end

for tt = 1:length(specialFileTypes)
    StimulusTypeMatcher(tt) = strcmp(stimulusTag, specialFileTypes{tt});
end

if length(find(StimulusTypeMatcher)) == 1
    StimulusName = specialFileTypes{find(StimulusTypeMatcher)};
elseif contains(filepath,'SparseNoise_fullscreen') == 1 % If equal SparseNoise_fullscreen
    StimulusName = 'SparseNoise_fullscreen';
end

% Process other potential inputs
if exist('varargin','var')
    % If additional variables are passed in, use standard pairwise format
    nVargin = length(varargin)/2;
    for thisVar = 1:nVargin
        % Declare each variable with the parameter value passed through
        eval(sprintf('%s = varargin{%01d};',varargin{(thisVar-1)*2+1},(thisVar-1)*2+2));
    end    
end
if ~exist('StimulusName','var') || isempty(StimulusName)
    StimulusName = 'Standard';
    warning('No StimulusType specified in filepath: using default script for loading Bonsai stimulus parameters')
end

dateNames = regexp(filepath,'(?<year>\d+)-(?<month>\d+)-(?<day>\d+)','match');
if etime(datevec(dateNames{1}),[2022 04 07 00 00 00]) < 0 % Happened before April 7th
    loadflag = 'old';
else
    loadflag = 'new';
end

switch(StimulusName)
    case 'Standard'
        switch(loadflag)
            case 'old'
                datatable = import_BonVisionParamsStandard_trials(filepath);
            case 'new'
                datatable = import_BonVisionParamsStandard_New(filepath);
        end
    case {'StaticGratings','DriftingGratings','GratingsPhase'}
        switch(loadflag)
            case 'old'
                datatable = import_BonVisionParamsStaticGratings_trials(filepath);
            case 'new'
                thisTable = readtable(filepath,'Delimiter',{',','{','=','}','(',')'});
                varCols = [1 2 3 4];    % FrameNumber and FrameTime
                varNames = {'FrameNumber','FrameTime','FrameComputerDateVec','Orientation'};
                datatable = thisTable(:,varCols);
                datatable.Properties.VariableNames = varNames;
                datatable.FrameComputerDateVec = datevec(datatable{:,'FrameComputerDateVec'},'yyyy-mm-ddTHH:MM:SS.FFF');
        end
    case 'SparseNoise' 
        if ~exist('grid_size','var')
            dateNames = regexp(filepath,'(?<year>\d+)-(?<month>\d+)-(?<day>\d+)','match');
            if etime(datevec(dateNames{1}),[2024 04 01 00 00 00]) < 0 % Happened before April 1st
                grid_size = [8,8];
                warning('Using default grid_size of 8x8 for SparseNoise')
            else
                grid_size = [9 16]; % Ellie will use [9 16]
            end
        end
        datatable = import_BonVisionParamsSparseNoise_bin(filepath,grid_size);
        
    case 'SparseNoise_fullscreen'
        grid_size = [8 16];
        datatable = import_BonVisionParamsSparseNoise_bin(filepath,grid_size);

    case {'TRAIN', 'DCBA', 'OMIT', 'E_CD', 'ADCD', 'lowcontB', 'lowcontDsubbingB', 'lowcontTRAIN', ...
            'OP_Tuning', 'Direction_Tuning', 'D___', 'A___', 'TRAIN_1', 'TRAIN_2', 'D_CD', ...
            'A_50ms', 'A_500ms', 'A_200ms', 'A_300ms', 'OMIT50grey'...
            'A_1000ms', 'A_1000ms_1', 'A_1000ms_2', 'A_1000ms_25pc', 'A_1000ms_50pc', 'A_1000ms_75pc', 'E_1000ms', 'TRAIN250',...
            'F_150ms', 'F_150ms_1', 'F_150ms_2', 'F_1000ms'}
        thisTable = readtable(filepath,'Delimiter',{',','{','=','}','(',')'});
        varCols = [4 14 20 22];    % % Extract the relevant columns (4 contains delay, 14 contains Contrast
        % 20 contains the grating Orientation/Direction, 22 contains the ArduinoTime of presentation)
        varNames = {'Delay','Contrast','Orientation','Time'};
        datatable = thisTable(:,varCols);
        datatable.Properties.VariableNames = varNames;
        % Remove the last row (the last 'stimulus') because the PD doesn't go off after this stimulus 
        % EB 25/2/2026 - I have amended the script so that a
        % dummy pd_off is padded in to allow the final pd_on
        %datatable(end, :) = [];
        %datatable.FrameComputerDateVec = datevec(datatable{:,'FrameComputerDateVec'},'yyyy-mm-ddTHH:MM:SS.FFF');
        
    case {'GAVNIK_ABCD', 'GAVNIK DCBA', 'GAVNIK_ABCD_1', 'GAVNIK_ABCD_2',...
            'GAVNIK200_ABCD', 'GAVNIK250_ABCD','GAVNIK250_ABCD_1', 'GAVNIK250_ABCD_2', 'GAVNIK250 DCBA'} % Orientation value is "grey" for interstimulus grey screen.
        thisTable = readtable(filepath,'Delimiter',{',','{','=','}','(',')'});
        varCols = [1 3];    % % Extract the relevant columns (1 contains Orienation, 3 contains the ArduinoTime of presentation)
        varNames = {'Orientation','Time'};
        datatable = thisTable(:,varCols);
        datatable.Properties.VariableNames = varNames;
        %datatable.FrameComputerDateVec = datevec(datatable{:,'FrameComputerDateVec'},'yyyy-mm-ddTHH:MM:SS.FFF');
        
        % Remove rows where Orientation is NaN
        datatable = datatable(~isnan(datatable.Orientation), :);
        
        
    case {'GAVNIK_A_CD', 'GAVNIK_E_CD', 'GAVNIK_D_CD', 'GAVNIK200_A_CD', 'GAVNIK200_E_CD',...
            'GAVNIK250_A_CD', 'GAVNIK250_E_CD'} % Orientation value is "grey" for interstimulus grey screen.
        % Orientation is "omission" for omitted second element
        thisTable = readtable(filepath,'Delimiter',{',','{','=','}','(',')'});
        varCols = [1 3];    % % Extract the relevant columns (1 contains Orienation, 3 contains the ArduinoTime of presentation)
        varNames = {'Orientation','Time'};
        datatable = thisTable(:,varCols);
        datatable.Properties.VariableNames = varNames;
        %datatable.FrameComputerDateVec = datevec(datatable{:,'FrameComputerDateVec'},'yyyy-mm-ddTHH:MM:SS.FFF');
        
        % Deal with NaNs in Orientation column
        stim_orientations = datatable.Orientation;
        original_orientations = stim_orientations;  % Preserve original before modifying

        % Find indices where NaN occurs
        nan_indices = isnan(stim_orientations);
        % Create a logical mask to keep track of which values to remove
        remove_mask = false(size(stim_orientations));

        % Loop through stim_orientations to modify NaNs
        for i = 1:length(stim_orientations)
            if nan_indices(i)
                if i > 3 && isnan(original_orientations(i-3))
                    remove_mask(i) = true;  % Remove NaN if 3 steps earlier was also NaN (i.e. Remove NaN where it was "grey"
                else
                    stim_orientations(i) = 0;  % Otherwise, replace NaN with 0 (i.e. where NaN was "omission" of second element)
                end
            end
        end

        % Apply the cleaned orientations back to the table
        datatable.Orientation = stim_orientations;
        % Remove rows marked for removal
        datatable(remove_mask, :) = [];
    
    case {'GAVNIK_A___', 'GAVNIK_D___'}  % Single-element GAVNIK stimuli (A or D alone) 
        % as there are hardly any numbers in the logging file, Matlab treats all entries as strings.
        % Read in the stimulus CSV file
        thisTable = readtable(filepath,'Delimiter',{',','{','=','}','(',')'});
        varCols = [1 3];    % Columns 1 = Orientation, 3 = ArduinoTime
        varNames = {'Orientation','Time'};
        datatable = thisTable(:,varCols);
        datatable.Properties.VariableNames = varNames;

        % Replace NaN orientations (BonVision logs omissions as NaN)
        stim_orientations = datatable.Orientation;

        % Convert Orientation entries to numeric values
        numeric_orientations = [];
        remove_mask = false(length(stim_orientations),1);   %% Track which rows to remove
        for i = 1:length(stim_orientations)
            s = stim_orientations{i};

            % Try to convert to number
            num = str2double(s);
            
            if strcmpi(s, 'grey')
                remove_mask(i) = true;              %% mark for removal
                continue;                           %% skip adding anything for "grey"
            elseif ~isnan(num)
                % If it's a valid number (e.g., '10', '45', etc.), use it directly
                numeric_orientations(end+1) = num;
            elseif strcmpi(s, 'Omission') % Give NaN values (logged as "omission") in stim_orientations a value of 0
                numeric_orientations(end+1) = 0;  % Use 0 for omission
            else
                warning('Unrecognized stimulus orientation: %s', s);
            end
        end

        % Store cleaned numeric orientations in the table
        datatable(remove_mask, :) = []; %% remove “grey” rows
        stim_orientations = numeric_orientations';
        datatable.Orientation = stim_orientations;

        
    case {'StaticTF','DriftingTF'}
        switch(loadflag)
            case 'old'
                thisTable = readtable(filepath,'Delimiter',{',','{','=','}','(',')'});
                datatable = thisTable(:,{'Var1','Var2','Var4','Var5'});
                datatable.Properties.VariableNames = {'FrameNumber','FrameTime','StimIndex','Time'};
            case 'new'
                thisTable = readtable(filepath,'Delimiter',{',','{','=','}','(',')'});
                varCols = [1 2 3 5 6];    % FrameNumber and FrameTime
                varNames = {'FrameNumber','FrameTime','FrameComputerDateVec','StimIndex','Time'};
                datatable = thisTable(:,varCols);
                datatable.Properties.VariableNames = varNames;
                datatable.FrameComputerDateVec = datevec(datatable{:,'FrameComputerDateVec'},'yyyy-mm-ddTHH:MM:SS.FFF');
        end
    case 'Grey'
        % April 7th 2022
        datatable = table();
        datatable.FrameNumber = NaN;
        datatable.FrameTime = 0;
        datatable.StimIndex = NaN;
        datatable.Time = NaN;
    case {'OriAdapt','PosAdapt'}
        % AP 20/4/2022
        ta = readtable(filepath,'Delimiter',{',','{','=','}','(',')'});
        varCols = [1 2 3];    % FrameNumber and FrameTime
        ta = ta(:,varCols);
        %[a] = cellfun(@(x) textscan(x,'%s','Delimiter',':'),thisTable.Value,'UniformOutput',false);
        tb = readtable(filepath,'Delimiter',{',','{','=','}','(',')',':'});
        ta = horzcat(ta,tb(:,{'Var7','Var8'}));
        ta.Var9 = ta.Var7; ta.Var10 = ta.Var8;
        datatable = horzcat(ta(1:2:end,[1:3,5]),ta(2:2:end,7));
        if strcmp(StimulusName,'OriAdapt')
            varNames = {'FrameNumber','FrameTime','FrameComputerDateVec','Orientation','Contrast'};
        else
            varNames = {'FrameNumber','FrameTime','FrameComputerDateVec','Position','Contrast'};
        end
        datatable.Properties.VariableNames = varNames;
        datatable.FrameComputerDateVec = datevec(datatable{:,'FrameComputerDateVec'},'yyyy-mm-ddTHH:MM:SS.FFF');
        timevec = datatable.FrameComputerDateVec;
        timevecDayStart = timevec;
        timevecDayStart(:,4:end) = 0; % set to beginning of day
        datatable.Time = etime(timevec,timevecDayStart)*1000;
end

end

function datatable = import_BonVisionParamsStaticGratings_trials_New(filepath)
thisTable = readtable(filepath,'Delimiter',{',','{','=','}','(',')'});
varCols = [1 2 3 4];    % FrameNumber and FrameTime
varNames = {'FrameNumber','FrameTime','FrameComputerDateVec','Orientation'};
datatable = thisTable(:,varCols);
datatable.Properties.VariableNames = varNames;
datatable.FrameComputerDateVec = datevec(datatable{:,'FrameComputerDateVec'},'yyyy-mm-ddTHH:MM:SS.FFF');

end

function datatable = import_BonVisionParamsStandard_New(filepath)
% Post April 7th 2022
thisTable = readtable(filepath,'Delimiter',{',','{','=','}','(',')'});
% Find strings in the table and read them in, and assign the
% next variable column as the data
varCols = [1 2 3];    % FrameNumber and FrameTime
varNames = {'FrameNumber','FrameTime','FrameComputerDateVec'};
td = table2cell(thisTable(1,:));
for kk=1:length(td)
    if ischar(td{kk}) && ~isempty(td{kk}) && isempty(str2num(td{kk}(1)) == 2)
        varNames = cat(2, varNames, td{kk});    % Add the variable column name
        varCols = cat(2, varCols, kk+1);        % Add the variable column
    end
end
datatable = thisTable(:,varCols);
datatable.Properties.VariableNames = varNames;
% Convert computer time to
datatable.FrameComputerDateVec = datevec(datatable{:,'FrameComputerDateVec'},'yyyy-mm-ddTHH:MM:SS.FFF');
end

function datatable = import_BonVisionParamsSparseNoise_bin(filename,grid_size)
% import SparseNoise binary file
fileID=fopen(filename);
thisBinFile=fread(fileID);
fclose(fileID);

% Translate stimulus into -1:1 scale
stim_matrix = zeros(1,length(thisBinFile));
stim_matrix(thisBinFile==0)=-1;
stim_matrix(thisBinFile==255)=1;
stim_matrix(thisBinFile==128)=0;

% Make a NxM grid from the stimulus log
stim_matrix = reshape(stim_matrix, [grid_size(1), grid_size(2), length(thisBinFile)/grid_size(1)/grid_size(2)]);
stim_matrix = stim_matrix(:,:,1:end-1); % The last 'stimulus'

% Make it a table
datatable = table;
for thisTrial = 1:size(stim_matrix,3)
    datatable.stim_matrix{thisTrial,1} = squeeze(stim_matrix(:,:,thisTrial));
end
end

function datatable = import_BonVisionParamsStandard_trials(filename)
% import dot corridor trial params file
%filename = 'X:\ibn-vision\DATA\SUBJECTS\M21192\bonsai\20210817\SpeedTextures_TrialParams2021-08-17T14_39_35.csv'
datatable_orig = readtable(filename);
nVars = width(datatable_orig);
nTrials = height(datatable_orig);

datatable = table;

for ivar = 1:nVars-1 % last var is just timestamp
    varStr = ['Var', num2str(ivar)];
    str  = char(datatable_orig.(varStr)(1));
    idx1 = find(isletter(str)); idx1 = idx1(1); % find first letter
    idx2 = strfind(str, '=');
    newVarName = str(idx1:idx2-1);
    
    newVals = nan(nTrials,1);
    for itrial = 1:nTrials
        str = char(datatable_orig.(varStr){itrial});
        idx = strfind(str,'=');
        if isempty(idx)
            newstr = str;
        else
            newstr = str(idx+1:end);
        end
        sscanf(newstr, '%g', 1);
        newVals(itrial) = sscanf(newstr, '%g', 1);
    end
    datatable.(newVarName) = newVals;    
end
end

function datatable = import_BonVisionParamsStaticGratings_trials(filename)
% import static or drifting gratings stimulus data file
% SGS 24/2/2022
datatable_orig = readtable(filename);
datatable = datatable_orig(2:end,:); % First line is the onset of grey screen
datatable = table2cell(datatable);
datatable(:,1) = cellfun(@(x) str2double(x),datatable(:,1),'UniformOutput',false);
datatable = cell2table(datatable);
datatable.Properties.VariableNames{1} = 'Orientation';
datatable.Properties.VariableNames{2} = 'Time';         % Changeed from 'Timestamp' to 'Time' o 10th April 2022

end

% function trialParams_table = import_FastSpeed_TrialParams(filename, startRow, endRow)
% %IMPORTFILE Import numeric data from a text file as a matrix.
% %   FASTSPEEDPARAMS20210817T153401 = IMPORTFILE(FILENAME) Reads data from
% %   text file FILENAME for the default selection.
% %
% %   FASTSPEEDPARAMS20210817T153401 = IMPORTFILE(FILENAME, STARTROW, ENDROW)
% %   Reads data from rows STARTROW through ENDROW of text file FILENAME.
% %
% % Example:
% %   FastSpeedParams20210817T153401 = importfile('FastSpeedParams2021-08-17T15_34_01.csv', 1, 1807);
% %
% %    See also TEXTSCAN.
% 
% % Auto-generated by MATLAB on 2021/08/21 17:10:22
% 
% %% Initialize variables.
% delimiter = ',';
% if nargin<=2
%     startRow = 2;
%     endRow = inf;
% end
% 
% %% Read columns of data as text:
% % For more information, see the TEXTSCAN documentation.
% formatSpec = '%s%s%s%s%s%[^\n\r]';
% 
% %% Open the text file.
% fileID = fopen(filename,'r');
% 
% %% Read columns of data according to the format.
% % This call is based on the structure of the file used to generate this
% % code. If an error occurs for a different file, try regenerating the code
% % from the Import Tool.
% dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
% for block=2:length(startRow)
%     frewind(fileID);
%     dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
%     for col=1:length(dataArray)
%         dataArray{col} = [dataArray{col};dataArrayBlock{col}];
%     end
% end
% 
% %% Close the text file.
% fclose(fileID);
% 
% %% Convert the contents of columns containing numeric text to numbers.
% % Replace non-numeric text with NaN.
% raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
% for col=1:length(dataArray)-1
%     raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
% end
% numericData = NaN(size(dataArray{1},1),size(dataArray,2));
% 
% for col=[1,2,3,4,5]
%     % Converts text in the input cell array to numbers. Replaced non-numeric
%     % text with NaN.
%     rawData = dataArray{col};
%     for row=1:size(rawData, 1)
%         % Create a regular expression to detect and remove non-numeric prefixes and
%         % suffixes.
%         regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
%         try
%             result = regexp(rawData(row), regexstr, 'names');
%             numbers = result.numbers;
%             
%             % Detected commas in non-thousand locations.
%             invalidThousandsSeparator = false;
%             if numbers.contains(',')
%                 thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
%                 if isempty(regexp(numbers, thousandsRegExp, 'once'))
%                     numbers = NaN;
%                     invalidThousandsSeparator = true;
%                 end
%             end
%             % Convert numeric text to numbers.
%             if ~invalidThousandsSeparator
%                 numbers = textscan(char(strrep(numbers, ',', '')), '%f');
%                 numericData(row, col) = numbers{1};
%                 raw{row, col} = numbers{1};
%             end
%         catch
%             raw{row, col} = rawData{row};
%         end
%     end
% end
% 
% 
% %% Replace non-numeric cells with NaN
% R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
% raw(R) = {NaN}; % Replace non-numeric cells
% 
% %% Create output variable
% trialParams_table = table;
% trialParams_table.PhotodiodeSign = cell2mat(raw(:, 1));
% trialParams_table.Vel = cell2mat(raw(:, 2));
% trialParams_table.Dist_lower = cell2mat(raw(:, 3));
% trialParams_table.Dist_upper = cell2mat(raw(:, 4));
% trialParams_table.VarName5 = cell2mat(raw(:, 5));
% end

