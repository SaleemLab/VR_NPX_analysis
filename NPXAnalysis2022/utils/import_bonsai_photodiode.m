% Funciton to import bonsai analog photodiode data
% SGS 03/04/2022
function [photodiode,photodiode_sync] = import_bonsai_photodiode(pd_filename,sync_filename)
dateNames = regexp(pd_filename,'(?<year>\d+)-(?<month>\d+)-(?<day>\d+)','match');
if etime(datevec(dateNames{1}),[2022 04 07 00 00 00]) < 0 % Happened before April 7th
    loadflag = 'old';
else
    loadflag = 'new';
end
switch(loadflag)
    case 'new' % From 7th April 2022
        % Read as table
        thisTable = readtable(pd_filename,'Delimiter',{',','(',')','{','=','}'});
%         pd_filename = 'Z:\ibn-vision\DATA\SUBJECTS\M23017\stimuli\20230629\SparseNoise_fullscreen_PDAsync2023-06-29T18_06_18.csv';
        varCols = [1 2 3];    % FrameNumber and FrameTime
        varNames = {'FrameNumber','FrameTime','FrameComputerDateVec'};
        td = table2cell(thisTable(1,:));
        for kk=1:length(td)
            if ischar(td{kk}) && ~isempty(td{kk}) && isempty(str2num(td{kk}(1)) == 2)
                varNames = cat(2, varNames, td{kk});    % Add the variable column name
                varCols = cat(2, varCols, kk+1);        % Add the variable column
            end
        end
        index = find(contains(varNames, 'SyncPulse', 'IgnoreCase', true));
        if ~isempty(index)
            varNames{index} = 'Sync';
        end

        if length(varCols) <= 3 % Hard coded for now, for visual stimuli (checkerboard and sparsenoise and static grating) PDAsync structure missing 'variable title', added manually
            varCols = [varCols 5 6 7]; % where 5 is photodiode signal, 6 is Sync pulse and 7 is timestamp
            varNames = [varNames,'Photodiode','Sync','Time'];
        end
        photodiode = thisTable(:,varCols);
        photodiode.Properties.VariableNames = varNames;
        % Convert computer time to
        photodiode.FrameComputerDateVec = datevec(photodiode{:,'FrameComputerDateVec'},'yyyy-mm-ddTHH:MM:SS.FFF');
        photodiode = table2struct(photodiode,'ToScalar',true);
        photodiode_sync = [];
        
        % datenum(thisTable{:,'ComputerTimestamp'},'yyyy-mm-ddTHH:MM:SS')
        % Following code fragment takes a long tie because it is in a
        % frameloop and there are hundreds of thousands of frames
        % Should be a better way of doing this
        if 0
            % Can go through and find the individual frames, find the real timestamp at the frame transition, then add the
            % difference from that timestamp to the Frametime
            frames = unique(photodiode.FrameNumber);
            photodiode.TimeReFrame = zeros(size(photodiode.FrameTime));
            for thisFrame = 1:length(frames)
                theseFrames = find(photodiode.FrameNumber == frames(thisFrame));
                photodiode.TimeReFrame(theseFrames) = photodiode.FrameTime(theseFrames)+photodiode.Time(theseFrames)-photodiode.Time(theseFrames(1));
            end
        end

%         if exist(photodiode.P) 
%             photodiode.P
%             photodiode.P = 
%             photodiode.Time = photodiode.T;
%         end
        
    case 'old' % Before 7th April 2022
        delimiter = ',';
        if nargin<=2
            startRow = 2;
            endRow = inf;
        end
        
        % Use generic import for photodiode
        photodiode = import_bonsai_file(pd_filename, delimiter, startRow, endRow);
        % Would like to use generic import for sync but currently encoded as
        % Boolean TRUE FALSE
        % Open the photodiode text file.
        fileID = fopen(sync_filename,'r');
        % Format for each line of text:
        formatSpec = '%f%s%[^\n\r]';
        dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
        for block=2:length(startRow)
            frewind(fileID);
            dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
            for col=1:length(dataArray)
                dataArray{col} = [dataArray{col};dataArrayBlock{col}];
            end
        end
        % Close the text file.
        fclose(fileID);
        % Reassign
        photodiode_sync.Time = dataArray{1};
        Sync = nominal(dataArray{2});
        photodiode_sync.Sync = zeros(size(Sync));
        photodiode_sync.Sync(Sync == nominal('True')) = 1;
        
end