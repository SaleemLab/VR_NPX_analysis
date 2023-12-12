% Function to import bnsai logs of eye data
% Revised old version, now using generic bonsai import function
% SGS 03/04/2022

function [EyeData] = import_eye_data(filename, startRow,endRow)
dateNames = regexp(filename,'(?<year>\d+)-(?<month>\d+)-(?<day>\d+)','match');
if etime(datevec(dateNames{1}),[2022 04 07 00 00 00]) < 0 % Happened before April 7th
    loadflag = 'old';
else
    loadflag = 'new';
end
switch(loadflag)
    case 'new'    % After 7th April 2022
        % Read as table
        thisTable = readtable(filename,'Delimiter',{',','{','=','}'});
        if isempty(thisTable)
            EyeData = [];
            disp('No eye data found')
            return
        end
        varCols = [1 2 3];    % FrameNumber and FrameTime
        varNames = {'FrameNumber','FrameTime','FrameComputerDateVec'};
        td = table2cell(thisTable(1,:));
        for kk=1:length(td)
            if ischar(td{kk}) && ~isempty(td{kk}) && isempty(str2num(td{kk}(1)) == 2)
                varNames = cat(2, varNames, td{kk});    % Add the variable column name
                varCols = cat(2, varCols, kk+1);        % Add the variable column
            end
        end
        EyeData = thisTable(:,varCols);
        EyeData.Properties.VariableNames = varNames;
        % Convert computer time to
        EyeData.FrameComputerDateVec = datevec(EyeData{:,'FrameComputerDateVec'},'yyyy-mm-ddTHH:MM:SS.FFF');
        EyeData.CentroidX = EyeData.X;
        EyeData.CentroidY = EyeData.Y;
        EyeData = table2struct(EyeData,'ToScalar',true);
        EyeData = rmfield(EyeData,'X');
        EyeData = rmfield(EyeData,'Y');
        
    case 'old'
        % Before 7th April 2022
        % Initialize variables
        delimiter = ',';
        if nargin < 2
            startRow = 2;
            endRow = inf;
        end
        
        % Import
        data = import_bonsai_file(filename,delimiter,startRow, endRow);
        % %
        % Allocate imported array to column variable names
        EyeData.CentroidX = data.Item1_Centroid_X;
        EyeData.CentroidY = data.Item1_Centroid_Y;
        EyeData.Orientation = data.Item1_Orientation;
        EyeData.MajorAxisLength = data.Item1_MajorAxisLength;
        EyeData.MinorAxisLength = data.Item1_MinorAxisLength;
        EyeData.Area = data.Item1_Area;
        EyeData.Time = data.Item2_Item1;
        EyeData.Sync = data.Item2_Item2;
        EyeData.EyeFrameCount = data.Item3;
end