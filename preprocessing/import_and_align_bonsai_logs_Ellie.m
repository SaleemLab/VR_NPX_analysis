% function [stimulusData,eyeData,wheelData,photodiodeData,stimTimes] = import_and_align_bonsai_logs(EPHYS_DATAPATH,TRIALDATA_DATAPATH, PERIPHERALS_DATAPATH,EYEDATA_DATAPATH)
%
% Program to import bonsai logs and eye data and align to ePhys data
% Modified from [stimulusData,eyeData,wheelData,photodiodeData,stimTimes] = importAndAlignBonsaiLogs(EPHYS_DATAPATH,TRIALDATA_DATAPATH, PERIPHERALS_DATAPATH,EYEDATA_DATAPATH)
% Now called import_and_align_bonsai_logs (just for differentiation)
% inputs : EPHYS_DATAPATH, TRIALDATA_DATAPATH, PERIPHERALS_DATAPATH, EYEDATA_DATAPATH
% outputs: stimulusData,eyeData,wheelData,photodiodeData (photdiode data derived from the
%               peripherals dataset)
%
% History:  SGS 24/2/2022 ported from import_process_BonsaiData
%           SGS 28/2/2022 added stimulus data parsing
%           SGS 04/4/2022 shifted imec extract to here and allowed function to either call precomputed
%                           asyncpulse times or load anew from LF (long) or AP (very long) imec data
%           SGS 14/4/2022 moved to using toolbox from Jennifer Collonel
%                       (SpikeGLX_Datafile_Tools) to read the async pulse from the ephys file
%           FR 09/06/2022 changed how PD was being read due to stimuli
%           adjustment issues
%           MT 19/01/2024 wheelData is now just peripherals
% Dependencies:
%   import_bonsai_peripherals
%   import_bonsai_eyedata
%   import_BonVisionParams
%   import_processWheelEphys (currently disabled)
%   alignBonsaiAndEphysData
%
%   Uses the following libraries, which are included in the folder
%   https://billkarsh.github.io/SpikeGLX/ - SpikeGLX_Datafile_Tools
%
function [stimData,eyeData,peripherals,photodiodeData,stimTimes] = import_and_align_bonsai_logs_Ellie(EPHYS_DATAPATH,TRIALDATA_DATAPATH,PERIPHERALS_DATAPATH,EYEDATA_DATAPATH,PHOTODIODE_DATAPATH,options)

% Step 1: Get stimulus times and data parameters
if ~strcmp(options.paradigm,'photodiode_timestamp')
    [stimData,StimulusName] = import_BonVisionParams_Ellie(TRIALDATA_DATAPATH);
else % If based on photodiode timestamp in case trial info not saved properly
    % Parse for one of the following stimulus names (you can add to the list here, and to the switch case below)
    specialFileTypes = {'StaticGratings_short','StaticGratings_long','StaticGratings','DriftingGratings','SparseNoise','SparseNoise_fullscreen',...
        'StaticTF','Grey','GratingsPhase','OriAdapt','PosAdapt','DriftingTF','Checkerboard','FullScreenFlash', 'OP_Tuning',...
        'TRAIN', 'DCBA', 'OMIT', 'E_CD', 'ADCD', 'lowcontB', 'lowcontDsubbingB', 'lowcontTRAIN', 'GAVNIK_ABCD', 'GAVNIK_A_CD',...
        'GAVNIK_E_CD', 'GAVNIK DCBA', 'GAVNIK_D___', 'GAVNIK_A___', 'GAVNIK_ABCD_1', 'GAVNIK_ABCD_2', 'GAVNIK_D_CD',...
        'TRAIN_1', 'TRAIN_2', 'D___', 'A___', 'D_CD',...
        'A_50ms', 'A_500ms', 'A_200ms', 'A_300ms', 'OMIT50grey', 'Direction_Tuning'...
        'A_1000ms', 'A_1000ms_1', 'A_1000ms_2', 'A_1000ms_25pc', 'A_1000ms_50pc', 'A_1000ms_75pc', 'E_1000ms', 'TRAIN250',...
        'GAVNIK200_ABCD', 'GAVNIK200_A_CD', 'GAVNIK200_E_CD', 'GAVNIK250_ABCD'};

    % See if there is a match within the file name that indicates that this is
    % from the special list
    StimulusTypeMatcher = zeros(1,length(specialFileTypes));
    [~, filename, ~] = fileparts(TRIALDATA_DATAPATH);
    
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

    stimData = [];
    stimData.Properties.VariableNames = [];
end

% Step 2: import behavioural / eye etc

% peripherals
peripherals= [];
thisTable = readtable(PERIPHERALS_DATAPATH,'Delimiter',{',','=','{','}','(',')'});
thisTable = rmmissing(thisTable,1); % Remove nan due to delimiter

for n = thisTable.Properties.VariableNames
    peripherals.(thisTable.Properties.VariableNames{n}) = thisTable.(thisTable.Properties.VariableNames{n});
end


% Video data
eyeData= [];
thisTable = readtable(EYEDATA_DATAPATH,'Delimiter',{',','=','{','}','(',')'});
thisTable = rmmissing(thisTable,1); % Remove nan due to delimiter

for n = thisTable.Properties.VariableNames
    eyeData.(thisTable.Properties.VariableNames{n}) = thisTable.(thisTable.Properties.VariableNames{n});
end

% Photodiode
if exist('PHOTODIODE_DATAPATH','var') && ~isempty(PHOTODIODE_DATAPATH)
    thisTable = readtable(PHOTODIODE_DATAPATH,'Delimiter',{',','=','{','}','(',')'});
    thisTable = rmmissing(thisTable,1); % Remove nan due to delimiter

    for n = thisTable.Properties.VariableNames
        photodiode.(thisTable.Properties.VariableNames{n}) = thisTable.(thisTable.Properties.VariableNames{n});
    end

else
    photodiode = []; photodiode_sync = [];
end

% Step 3: obtain imec sync pulse from precomputed or imec file
Nidq = [];
td = dir(fullfile(EPHYS_DATAPATH,'..','*NidqTimes*'));
if ~isempty(td)
    load(fullfile(EPHYS_DATAPATH,'..',td(1).name),'Nidq');
    syncTimes_ephys = Nidq.bonsai_sync_on;% use upswings currently
    syncPulse_ephys = Nidq.bonsai_sync;
    syncPulse_ephysTimes = Nidq.sglxTime;
else
    error('nidq and ephys sync pulse extraction and alignment not done!')
    return
end


% Step 4: align nidq and Bonsai logged data

%%%%%%
% For photodiode and convert to spike GLX timeThe
if exist('photodiode','var') && ~isempty(photodiode)
    photodiode.Sync =  photodiode.AsyncPulse;
    photodiode.Time = photodiode.ArduinoTime;
%     photodiode.Time = photodiode.BonsaiTime;

    [photodiode] = alignBonsaiToEphysSyncTimes(photodiode,syncTimes_ephys);
%     [photodiode] = alignBonsaiToEphysSyncPulse(photodiode,syncPulse_ephys,syncPulse_ephysTimes);
    
else
    photodiode = [];
end


%%%%%%
% For wheel data and eyedaya and convert to spike GLX time-base
[C,ia,ic]=unique(photodiode.ArduinoTime);
peripherals.sglxTime = interp1(photodiode.ArduinoTime(ia),photodiode.sglxTime(ia),peripherals.ArduinoTime,'previous');
eyeData.sglxTime = interp1(photodiode.ArduinoTime(ia),photodiode.sglxTime(ia),eyeData.ArduinoTime,'previous');
%wheelData = peripherals;

disp(['StimulusName: ', StimulusName]);
switch(StimulusName)
    case {'GAVNIK_ABCD', 'GAVNIK_A_CD', 'GAVNIK_E_CD', 'GAVNIK DCBA', 'GAVNIK_A___', 'GAVNIK_D___', 'GAVNIK_ABCD_1',...
            'GAVNIK_ABCD_2', 'GAVNIK_D_CD', 'GAVNIK200_ABCD', 'GAVNIK200_A_CD', 'GAVNIK200_E_CD', 'GAVNIK250_ABCD',} % for GAVNIK stimuli, the PD quad begins grey, goes white for first stim A, black for second B,
            % white for third C, black for fourth D and then grey for the greyscreen period between the sequences of four stimuli
        
        % Compute the trailing moving median (of the current value and the four values prior)
        photodiode_median = movmedian(photodiode.Photodiode, [4 0]);

        % Initialize transition markers
        pd_ON_OFF = zeros(size(photodiode.Photodiode));

        % Threshold conditions
        %up_condition_1 = (photodiode_median < 280 & photodiode_median > 220);  
        %up_condition_2 = (photodiode_median < 110);
        %down_condition = (photodiode_median > 380);

        % Cluster the photodiode signal into 3 levels (black, grey, white)
        [idx, centers] = kmeans(double(photodiode.Photodiode(:)), 3);  
        centers = sort(centers);  % ensure order: low < mid < high
        centers = centers(:)';    % force row vector [1×3]

        low_center  = centers(1);
        mid_center  = centers(2);
        high_center = centers(3);

        % Assign photodiode_median to the closest cluster
        [~, cluster_id] = min(abs(photodiode_median(:) - centers), [], 2);
        
        % Define conditions by cluster ID
        up_condition_1 = (cluster_id == 2);  % grey/mid
        up_condition_2 = (cluster_id == 1);  % black/low
        down_condition = (cluster_id == 3);  % white/high

        % Track the last detected transition
        last_transition = -1;  % -1 means no transition detected yet
        state = 0;  % - 0 = off state, 1 = on state

        % Loop through the signal to detect transitions
        for i = 1:length(photodiode.Photodiode)-50 % Ensure index is within bounds
            
            % Define tolerance (you can tweak this)
            tol = 0.3 * (high_center - mid_center);

        % Check upward transitions (only if last transition was NOT up)
            if (last_transition ~= 1) && ((up_condition_1(i) && all(photodiode.Photodiode(i+4:i+5) > (high_center - tol))) || ...% if coming from grey quad and the next value and three subsequent values exceed 280
                                  (up_condition_2(i) && all(photodiode.Photodiode(i+7:i+9) > (high_center - tol)))) % if coming from black quad and the values of positions 4-7 ahead exceed 350 (i.e. when the quad is going to white but not grey, an upswing will be detected).
                pd_ON_OFF(i) = 1;
                state = 1;          % Set state to ON
                last_transition = 1;  % Mark last transition as UP
                disp(['UP Transition at ', num2str(i)]);  % Debugging output
            end

            % Keep the state as 1 until a downswing is detected
            if state == 1
                pd_ON_OFF(i) = 1;  % Keep state as 1 (ON)
            end

            % Check downward transitions (only if last transition was NOT down)
            if (last_transition ~= 0) && (down_condition(i) && all(photodiode.Photodiode(i:i+1) < (high_center - tol)))
                pd_ON_OFF(i) = 0;
                state = 0;          % Set state to OFF
                last_transition = 0;  % Mark last transition as DOWN
                disp(['DOWN Transition at ', num2str(i)]);  % Debugging output
            end
    
        end

        % Extract ON and OFF transition indices
        pd_ON = find(diff(pd_ON_OFF) == 1);  % Find the indices of ON transitions
        pd_OFF = find(diff(pd_ON_OFF) == -1);  % Find the indices of OFF transitions
        %pd_ON_times = Nidq.sglxTime(pd_ON);
        %pd_OFF_times = Nidq.sglxTime(pd_OFF);

    otherwise
        photodiode.Photodiode_smoothed = smoothdata(photodiode.Photodiode,'movmedian',5);

        pd_ON_OFF = photodiode.Photodiode_smoothed >= mean(photodiode.Photodiode); % find ON and OFF
        pd_ON = find(diff(pd_ON_OFF) == 1); % Find ON
        pd_OFF = find(diff(pd_ON_OFF) == -1); % Find OFF

        % Remove any photodiode transitions that occur >5 seconds after the last transition (PD fluctuations can occur when turning off Bonsai)
        % Convert indices to SpikeGLX time
        pd_ON_times = photodiode.sglxTime(pd_ON);
        pd_OFF_times = photodiode.sglxTime(pd_OFF);
        %pd_ON_times = Nidq.sglxTime(pd_ON);
        %pd_OFF_times = Nidq.sglxTime(pd_OFF);

        % Find time differences between consecutive transitions
        pd_ON_diffs = diff(pd_ON_times);
        pd_OFF_diffs = diff(pd_OFF_times);

        %%%%% Keep only transitions where the time difference is <= 5 seconds
        valid_ON_idx = [true; pd_ON_diffs <= 5];  % Keep first transition and filter the rest
        valid_OFF_idx = [true; pd_OFF_diffs <= 5];
        %valid_ON_idx = [true, pd_ON_diffs <= 5];  % Keep first transition and filter the rest
        %valid_OFF_idx = [true, pd_OFF_diffs <= 5];
        % Find the last valid ON transition
        last_valid_ON = find(valid_ON_idx, 1, 'last');

        % Remove any transitions occurring more than 5 seconds after the last valid ON
        valid_ON_idx((last_valid_ON + 1):end) = false;
        valid_OFF_idx((last_valid_ON + 1):end) = false;

        % Apply filtering
        pd_ON = pd_ON(valid_ON_idx);
        pd_OFF = pd_OFF(valid_OFF_idx);

        % Convert back to SpikeGLX time
        pd_ON_times = photodiode.sglxTime(pd_ON);
        pd_OFF_times = photodiode.sglxTime(pd_OFF);
        %pd_ON_times = Nidq.sglxTime(pd_ON);
        %pd_OFF_times = Nidq.sglxTime(pd_OFF);

end


% Now calculate the duration (in number of samples taken during the presentation) of each ON presentation
if length(pd_ON) == length(pd_OFF)
    block_length = abs(pd_ON-pd_OFF);
elseif length(pd_ON)>length(pd_OFF)
    for i = 1:length(pd_OFF)
        block_length(i) = abs(pd_ON(i) - pd_OFF(i));
    end
else
    for i = 1:length(pd_ON)
        block_length(i) = abs(pd_ON(i) - pd_OFF(i));
    end
end

blocks_ind = find(block_length>10); % sample rate of photodiode is 1000 per second so a block size should be quite large

% Remove the last pd_ON transition if there is a mismatch (since pd_ON is one greater than pd_OFF)
if length(pd_ON) > length(pd_OFF)
    pd_ON = pd_ON(1:end-1); % Remove the last ON transition
end

% Ensure blocks_ind is valid within the length of pd_OFF (now pd_ON and pd_OFF are aligned)
blocks_ind = blocks_ind(1:min(end, length(pd_OFF))); % Trim blocks_ind to match pd_OFF length

   % plot(photodiode.sglxTime,photodiode.Photodiode_smoothed); hold on; scatter(photodiode.sglxTime(pd_OFF(blocks_ind)),50,'b');
   % scatter(photodiode.sglxTime(pd_ON(blocks_ind)),50,'r');

   % plot(Nidq.sglxTime,Nidq.photodiode/100);
photodiodeData = [];
photodiodeData.stim_on.sglxTime = photodiode.sglxTime(pd_ON(blocks_ind)');
photodiodeData.stim_off.sglxTime = photodiode.sglxTime(pd_OFF(blocks_ind)');

whos photodiodeData photodiodeData.stim_off
class(photodiodeData.stim_off.sglxTime), size(photodiodeData.stim_off.sglxTime)

whos Nidq
class(Nidq.sglxTime), size(Nidq.sglxTime)


if photodiodeData.stim_off.sglxTime(end) > Nidq.sglxTime(end)
    % if last visual stimuli is later than last timestamp of the recording
    % then it is probably because of wrong alignment (happens normally only during short sessions such as checkerboard)
    % in this rare case use Nidq photodiode directly

    Nidq.photodiode_smoothed = smoothdata(Nidq.photodiode,'movmedian',50);

    pd_ON_OFF = Nidq.photodiode_smoothed >= mean(Nidq.photodiode); % find ON and OFF
    pd_ON = find(diff(pd_ON_OFF) == 1); % Find ON
    pd_OFF = find(diff(pd_ON_OFF) == -1); % Find OFF

    if length(pd_ON) == length(pd_OFF)
        block_length = abs(pd_ON-pd_OFF);
    elseif length(pd_ON)>length(pd_OFF)
        for i = 1:length(pd_OFF)
            block_length(i) = abs(pd_ON(i) - pd_OFF(i));
        end
    else
        for i = 1:length(pd_ON)
            block_length(i) = abs(pd_ON(i) - pd_OFF(i));
        end
    end

    blocks_ind = find(block_length>10); % sample rate of photodiode is 1000 per second so a block size should be quite large

    % plot(photodiode.sglxTime,photodiode.Photodiode_smoothed); hold on; scatter(photodiode.sglxTime(pd_OFF(blocks_ind)),50,'b');
    % scatter(photodiode.sglxTime(pd_ON(blocks_ind)),50,'r');

    % plot(Nidq.sglxTime,Nidq.photodiode/100);
    photodiodeData = [];
    photodiodeData.stim_on.sglxTime = Nidq.sglxTime(pd_ON(blocks_ind))';
    photodiodeData.stim_off.sglxTime = Nidq.sglxTime(pd_OFF(blocks_ind))';
    photodiodeData.warning = 'Nidq photodiode';
    disp('Nidq photodiode directly due to bad alignment between Bonsai and Ephys')
end


% Step 5: process wheel data (skipping this just save all peripheral data)
% wheelData.pos = peripherals.Wheel;
% wheelData.Time = peripherals.sglxTime;
% wheelData = import_processWheelEphys(wheelData,'gaussian',12);
% wheelData.sglxTime = wheelData.Time;

% Step 6: extract stimulus times
disp(['StimulusName: ', StimulusName]);
switch(StimulusName)
    case {'sparsenoise', 'SparseNoise_1', 'SparseNoise_2', 'sparsenoise_fullscreen','checkerboard'}
        stimTimes = sort(cat(1,photodiodeData.stim_on.sglxTime,photodiodeData.stim_off.sglxTime),'ascend');
        %         fprintf('\n\tUsing both photodiode upward and downward transitions as timestamps for stimulus timing')
        if contains(StimulusName,'SparseNoise')
            
            average_stim_duration = (stimTimes(end)-stimTimes(1))/size(stimData,1);
            stim_id = (find(diff(stimTimes)>average_stim_duration*2));

            if isempty(stim_id)
                
            elseif length(stim_id) == 1 % if happended once we can remove all the stimuli stime between these two points (happens rarely)
                stimData(stim_id:stim_id+size(stimData,1)-length(stimTimes)-1,:)=[];
            elseif length(stim_id) <4 % if happended few times we can remove all the stimuli stime between these two points (happens rarely)
                stimData(stim_id(1):end,:)=[];
                stimTimes(stim_id:end)=[];
            end

        end

    case {'TRAIN', 'OP_Tuning', 'Direction_Tuning', 'DCBA', 'OMIT', 'E_CD', 'D___', 'A___', 'TRAIN_1', 'TRAIN_2', 'D_CD',...
            'A_50ms', 'A_500ms', 'A_200ms', 'A_300ms', 'OMIT50grey', 'ADCD', 'lowcontB', 'lowcontDsubbingB', 'lowcontTRAIN', 'staticgratings','staticgratings_short','staticgratings_long',...
            'A_1000ms', 'A_1000ms_1', 'A_1000ms_2', 'A_1000ms_25pc', 'A_1000ms_50pc', 'A_1000ms_75pc', 'E_1000ms', 'TRAIN250'}
        stimTimes = sort(cat(1,photodiodeData.stim_on.sglxTime),'ascend'); %cases where a stimulus is on the screen only when the quad is white
    
    case {'GAVNIK_ABCD', 'GAVNIK_A_CD', 'GAVNIK_E_CD', 'GAVNIK DCBA', 'GAVNIK_D___', 'GAVNIK_A___', 'GAVNIK_ABCD_1', 'GAVNIK_ABCD_2',...
            'GAVNIK_D_CD', 'GAVNIK200_ABCD', 'GAVNIK200_A_CD', 'GAVNIK200_E_CD', 'GAVNIK250_ABCD'}
        stimTimes = sort(cat(1,photodiodeData.stim_on.sglxTime,photodiodeData.stim_off.sglxTime),'ascend'); % Sort the concatenated array in ascending (i.e. chronological) order
    
    case 'grey'
        stimTimes = photodiodeData.stim_on.sglxTime;
    case {'oriadapt','posadapt'}
        %         PD_timestamps = zeros(1,length(photodiodeData.stim_on.bonsaiTime)+length(photodiodeData.stim_off.bonsaiTime));
        %         PD_timestamps(1:2:length(photodiodeData.stim_off.bonsaiTime)*2) = photodiodeData.stim_off.bonsaiTime;
        %         PD_timestamps(2:2:length(photodiodeData.stim_on.bonsaiTime)*2) = photodiodeData.stim_on.bonsaiTime;
        %
        %         OE_ts = zeros(1,length(photodiodeData.stim_on.sglxTime)+length(photodiodeData.stim_off.sglxTime));
        %         OE_ts(1:2:length(photodiodeData.stim_off.sglxTime)*2) = photodiodeData.stim_off.sglxTime;
        %         OE_ts(2:2:length(photodiodeData.stim_on.sglxTime)*2) = photodiodeData.stim_on.sglxTime;

        PD_timestamps = photodiodeData.stim_on.bonsaiTime;

        OE_ts = photodiodeData.stim_on.sglxTime;

        % Do the correction for every photodiode reversal
        stimTimes = stimData.Time./1000;

        for i=1:length(PD_timestamps)-1
            bonsai_ts = PD_timestamps(i);
            oe_ts = OE_ts(i);
            offset(i) = bonsai_ts-oe_ts;
            if i==1
                inds = find(stimData.Time./1000<=PD_timestamps(i));
                stimTimes(inds) = stimTimes(inds)-offset(i);
            end
            inds = find(stimData.Time./1000>=PD_timestamps(i) & stimData.Time./1000<PD_timestamps(i+1));
            stimTimes(inds) = stimTimes(inds)-offset(i);
        end
        i = length(PD_timestamps);
        bonsai_ts = PD_timestamps(i);
        oe_ts = OE_ts(i);
        offset(i) = bonsai_ts-oe_ts;
        inds = find(stimData.Time./1000>PD_timestamps(i));
        stimTimes(inds) = stimTimes(inds)-offset(i);
        stimData.stimTimes = stimTimes;


    otherwise
        % Either use the stimulus times from the stimulus log or photodiode
        %
        % For use of stimulus logged times...
        % stimTimes = stimData.sglxTime;
        % fprintf('\n\tUsing stimulus logged times as timestamps for stimulus timing')
        %
        % For use times from the photodiode upswing
        stimTimes = photodiodeData.stim_on.sglxTime;
        fprintf('\n\tUsing photodiode upward transitions as timestamps for stimulus timing')
end
fprintf('\n'); warning('Havent yet enabled checking of number of stimulus expected vs number of photodiode changes'); fprintf('\n');

