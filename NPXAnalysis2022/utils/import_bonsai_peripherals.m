% Modified by M.T on 01/07 to also extract virtual position and asynchrouns
% pulse and photodiode signal. (photodiode signal is sampled at 60Hz in
% wheel log, so usually we should't use it)

function [peripherals] = import_bonsai_peripherals(filename,startRow,endRow)
% Scan date in filepath to determine which import method to use

dateNames = regexp(filename,'(?<year>\d+)-(?<month>\d+)-(?<day>\d+)','match');
if etime(datevec(dateNames{1}),[2022 04 07 00 00 00]) < 0 % Happened before April 7th
    loadflag = 'old';
else
    loadflag = 'new';
end
switch(loadflag)
    case 'new'
        % Read as table - after 7th April 2022
%         thisTable = readtable(filename,'Delimiter',{',','{','=','}'});
        thisTable = readtable(filename,'Delimiter',{',','=','{','}','(',')'});
        thisTable = rmmissing(thisTable,2); % Remove nan due to delimiter


        varCols = [1 2 3];    % FrameNumber and FrameTime
        varNames = {'FrameNumber','FrameTime','FrameComputerDateVec'};
        td = table2cell(thisTable(1,:));
        for kk=1:length(td)
            if ischar(td{kk}) && ~isempty(td{kk}) && isempty(str2num(td{kk}(1)) == 2)
                varNames = cat(2, varNames, td{kk});    % Add the variable column name
                varCols = cat(2, varCols, kk+1);        % Add the variable column
            end
        end
        peripherals = thisTable(:,varCols);
        peripherals.Properties.VariableNames = varNames;
        % Convert computer time to
        peripherals.FrameComputerDateVec = datevec(peripherals{:,'FrameComputerDateVec'},'yyyy-mm-ddTHH:MM:SS.FFF');
        peripherals = table2struct(peripherals,'ToScalar',true);

        % Add asynch pulse and photodiode
        if isfield(peripherals,'SyncPulse') % if SyncPulse (no VR Position)
            peripherals.Sync = peripherals.SyncPulse;
            peripherals = rmfield(peripherals,'SyncPulse');
            peripherals.Photodiode = table2array(thisTable(:,end));% Last column is photodiode
        elseif isfield(peripherals,'Position') % If VR position

            if size(thisTable,2)==20
                peripherals.Sync = table2array(thisTable(:,end-2));
                peripherals.Photodiode = table2array(thisTable(:,end-1));
                peripherals.VideoSync = table2array(thisTable(:,end));
            else
                peripherals.Sync = table2array(thisTable(:,end-1));
                peripherals.Photodiode = table2array(thisTable(:,end));
                peripherals.VideoSync = [];
            end
        end



    case 'old'
        % Before 7th April 2022
        % Initialize variables.
        delimiter = ',';
        if nargin<=2
            startRow = 2;
            endRow = inf;
        end
        
        % Replaced previous code with generic bonsai import:
        peripherals = import_bonsai_file(filename, delimiter, startRow, endRow);
        
        % If a field called Quadstate is not present, assign it the value of Photodiode
        if ~isfield(peripherals,'QuadState')
            if isfield(peripherals,'Photodiode')
                peripherals.QuadState = peripherals.Photodiode;
            else
                peripherals.Photodiode = [];
                peripherals.QuadState = [];
            end
        end
        
        
        %%%%%%%
        % OLD - hard code columns
        % % % switch(length(dataArray))
        % % %     case 6
        % % %         Pre 22/3/2022 when didn't store photdiode
        % % %         Wheel	Sync EyeFrameCount 'Photodiode' Time
        % % %         Wheel = dataArray{:, 1};
        % % %         Sync = dataArray{:, 2};
        % % %         EyeFrameCount = dataArray{:, 3};
        % % %         QuadState = dataArray{:, 4}; % note called photodiode in .csv file but actually quadstate
        % % %         Time = dataArray{:, 5};
        % % %         Photodiode = QuadState;
        % % %     case 7
        % % %         Post 22/3/2022 when stored photodiode as well as sync pulse
        % % %         Wheel	Photodiode	Sync	QuadState	EyeFrameCount	Time
        % % %         Wheel = dataArray{:, 1};
        % % %         Photodiode = dataArray{:, 2};
        % % %         Sync = dataArray{:, 3};
        % % %         QuadState = dataArray{:, 4};
        % % %         EyeFrameCount = dataArray{:, 5};
        % % %         Time = dataArray{:,6};
        % % %     case 8
        % % %         2021 NP2 when stored repcount as well as photodiode as well as sync pulse
        % % %         Wheel	Sync	Lick	EyeFrameCount	'Photodiode'	RepCounter	Time
        % % %         Wheel = dataArray{:, 1};
        % % %         Photodiode = dataArray{:, 2};
        % % %         Sync = dataArray{:, 3};
        % % %         QuadState = dataArray{:, 4};
        % % %         EyeFrameCount = dataArray{:, 5};
        % % %         Time = dataArray{:,6};
        % % %
        % % % end
        % % %
        % % % generate output struct
        % % % peripherals.Wheel = Wheel;
        % % % peripherals.Sync = Sync;
        % % % peripherals.EyeFrameCount = EyeFrameCount;
        % % % peripherals.Photodiode = Photodiode;
        % % % peripherals.QuadState = QuadState;
        % % % peripherals.Time = Time;
end




