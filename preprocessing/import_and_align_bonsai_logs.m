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
function [stimData,eyeData,wheelData,photodiodeData,stimTimes] = import_and_align_bonsai_logs(EPHYS_DATAPATH,TRIALDATA_DATAPATH,PERIPHERALS_DATAPATH,EYEDATA_DATAPATH,PHOTODIODE_DATAPATH,options)

% Step 1: Get stimulus times and data parameters
if ~strcmp(options.paradigm,'photodiode_timestamp')
    [stimData,StimulusName] = import_BonVisionParams(TRIALDATA_DATAPATH);
else % If based on photodiode timestamp in case trial info not saved properly
    % Parse for one of the following stimulus names (you can add to the list here, and to the switch case below)
    specialFileTypes = {'StaticGratings_short','StaticGratings_long','StaticGratings','DriftingGratings','SparseNoise','SparseNoise_fullscreen','StaticTF','Grey','GratingsPhase','OriAdapt','PosAdapt','DriftingTF','Checkerboard'};
    % See if there is a match within the file name that indicates that this is
    % from the special list
    StimulusTypeMatcher = zeros(1,length(specialFileTypes));
    for tt = 1:length(specialFileTypes)
        StimulusTypeMatcher(tt) = contains(TRIALDATA_DATAPATH,specialFileTypes{tt});
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
peripherals = import_bonsai_peripherals(PERIPHERALS_DATAPATH);
eyeData = import_bonsai_eyedata(EYEDATA_DATAPATH);
if exist('PHOTODIODE_DATAPATH','var') && ~isempty(PHOTODIODE_DATAPATH)
    [photodiode, photodiode_sync] = import_bonsai_photodiode(PHOTODIODE_DATAPATH);
else
    photodiode = []; photodiode_sync = [];
end

% Step 3: obtain imec sync pulse from precomputed or imec file
td = dir(fullfile(EPHYS_DATAPATH,'*syncpulseTimes*'));
if ~isempty(td)
    load(fullfile(EPHYS_DATAPATH,td(1).name),'syncTimes_ephys');
    syncTimes_ephys = syncTimes_ephys.on;% use upswings currently
else
    %%%%% Old - use Neuropixels memmap
    % %     imec = Neuropixel.ImecDataset(EPHYS_DATAPATH);
    % %     if imec.syncInLFFile
    % %         mm = imec.memmapLF_full(); % get LF data
    % %         vec = mm.Data.x(imec.syncChannelIndex,:); % get the sync channel from LFP data
    % %         vec_idx = find(diff(vec)>=0.5); % Find upswings in sync pulse
    % %         syncTimes_ephys.on = (vec_idx+1)./imec.fsLF; % add 1 to vec_idx to compensate for diff and convert from samples to s (fsLF is sample rate)
    % %         vec_idx = find(diff(vec)<=-0.5); % Find downswings in sync pulse
    % %         syncTimes_ephys.off = (vec_idx+1)./imec.fsLF; % add 1 to vec_idx to compensate for diff and convert from samples to s (fsLF is sample rate)
    % %     elseif imec.syncInAPFile
    % %         fprintf('\n'); warning('Extracting sync-pulses from AP file...this takes about 5 minutes'); fprintf('\n');
    % %         syncTimes_ephys = extractSyncPulseFromAPFile(imec);
    % %     else
    % %         error('This programme relies on sync source')
    % %     end
    
    %%%%% New - use SpikeGLX_Datafile_Tools
    % If there is an LF file use that, otherwise AP file
    [ap_file,lf_file] = findImecBinFile(EPHYS_DATAPATH);
    if ~isempty(lf_file)
        binpath = fullfile(EPHYS_DATAPATH,lf_file);
    else
        fprintf('\n'); warning('Extracting sync-pulses from AP file...this takes about 5 minutes'); fprintf('\n');
        binpath = fullfile(EPHYS_DATAPATH,ap_file);
    end
    syncTimes_ephys = SGLXextractSyncPulseFromBinFile(binpath);
    parseDate = date;
    [~,fname] = fileparts(EPHYS_DATAPATH);
    save(fullfile(EPHYS_DATAPATH,[fname,'_syncpulseTimes.mat']),'syncTimes_ephys','parseDate');
    syncTimes_ephys = syncTimes_ephys.on; % use upswings currently
end


% Step 4: align ephys and Bonsai logged data
[stimData,peripherals,eyeData,photodiode] = alignBonsaiAndEphysData(syncTimes_ephys,stimData,peripherals,eyeData,photodiode);
wheelData = peripherals;

if strcmp(options.paradigm,'photodiode_timestamp')
    % extract stimuli time based on photodiode timestamp

    photodiode.Photodiode_smoothed = smoothdata(photodiode.Photodiode,'movmedian',5);

    pd_ON_OFF = photodiode.Photodiode_smoothed >= 80; % find ON and OFF
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

    blocks_ind = find(block_length>200); % sample rate of photodiode is 3000 per second so a block size should be quite large

    %    plot(photodiode.Photodiode_smoothed); hold on; scatter(pd_OFF(blocks_ind),50);
    photodiodeData = [];
    photodiodeData.stim_on.sglxTime = photodiode.sglxTime (pd_ON(blocks_ind)');
    photodiodeData.stim_off.sglxTime = photodiode.sglxTime(pd_OFF(blocks_ind)');
end

% Step 5: get trial start and end times from photodiode QuadState (should
% replace with actual Photodiode at some point)
if isempty(photodiode)
    fprintf('\n'); warning('Stimulus times derived from peripherals.QuadState'); fprintf('\n');
    idx_trial_start = find(diff(peripherals.QuadState)>=0.5) + 1;
    idx_trial_end = find(diff(peripherals.QuadState)<=-0.5) + 1;
    photodiodeData.stim_on.bonsaiTime = peripherals.Time(idx_trial_start)./1000; % convert to seconds
    photodiodeData.stim_off.bonsaiTime = peripherals.Time(idx_trial_end)./1000;
    photodiodeData.stim_on.sglxTime = peripherals.sglxTime(idx_trial_start);
    photodiodeData.stim_off.sglxTime = peripherals.sglxTime(idx_trial_end);
elseif options.paradigm==1
    ttl = convertPDtoTimeStamps_FR(photodiode.Photodiode,400,120,1000,options);
    photodiodeData.stim_on.bonsaiTime = photodiode.Time(ttl.upPhases)./1000;
    photodiodeData.stim_on.sglxTime = photodiode.sglxTime(ttl.upPhases);
elseif strcmp(options.paradigm,'sn')
    ttl = convertPDtoTimeStamps_FR(photodiode.Photodiode,400,120,1000,options);
    photodiodeData.stim_on.bonsaiTime = photodiode.Time(ttl.upPhases)./1000;
    photodiodeData.stim_on.sglxTime = photodiode.sglxTime(ttl.upPhases);
    photodiodeData.stim_off.bonsaiTime = photodiode.Time(ttl.downPhases)./1000;
    photodiodeData.stim_off.sglxTime = photodiode.sglxTime(ttl.downPhases);
elseif strcmp(options.paradigm,'ad')
    %Amalia included this. May not work for other paradigms.
    options.photo_br = 0;
    ttl = convertPDtoTimeStamps(photodiode.Photodiode,300,120,1000,options);
    photodiodeData.stim_on.bonsaiTime = photodiode.Time(ttl.upPhases)./1000;
    photodiodeData.stim_on.sglxTime = photodiode.sglxTime(ttl.upPhases);

elseif strcmp(options.paradigm,'masa')
    idx_trial_start = find(diff(peripherals.QuadState)==1) + 1;
    idx_trial_end = find(diff(peripherals.QuadState)==-1) + 1;
    frame_trial_start = peripherals.FrameNumber(idx_trial_start);
    frame_trial_end = peripherals.FrameNumber(idx_trial_end);
    temp_tbl.start_time = peripherals.sglxTime(idx_trial_start);
    temp_tbl.end_time = peripherals.sglxTime(idx_trial_end);
    nTrials = length(temp_tbl.start_time);
    % get trial start/end times from photodiode signal
    pd_thresh_up = 400;
    pd_thresh_down = 100;
    for itrial =1:nTrials
        % find frame where Bonsai asked to change quad
        frameidx_start = find(photodiode.FrameNumber <= frame_trial_start(itrial),1,'last');
        % find next index where the photodiode detected a quad change
        temp_idx = find(photodiode.Photodiode(frameidx_start:end)>pd_thresh_up,1,'first');
        idx_start = frameidx_start+temp_idx-1;
        % convert to spike GLX time
        pdstart(itrial) = photodiode.sglxTime(idx_start);
        % to get photodiode down we need to do some median filtering to
        % remove some awkward down swings when the quad is white
        frameidx_end = find(photodiode.FrameNumber >= frame_trial_end(itrial),1,'first');
        temp_idx = find(smoothdata(photodiode.Photodiode(frameidx_end:end),'movmedian', 5)<pd_thresh_down,1,'first');
        idx_end = frameidx_end+temp_idx-1;
        pdend(itrial) = photodiode.sglxTime(idx_end);
    end

    photodiodeData.stim_on.sglxTime = pdstart';
    photodiodeData.stim_off.sglxTime = pdend';

elseif strcmp(options.paradigm,'SG')
        %hot fix of the static grating quadstate issue
     
    idx_trial_start = find(diff(peripherals.QuadState)==1) + 1;
    idx_trial_end = find(diff(peripherals.QuadState)==-1) + 1;
    idx_trial_diff = diff(idx_trial_start);
    idx_trial_diff_end = diff(idx_trial_end);
    idx_trial_start_tmp = zeros(size(stimData,1),1);
    idx_trial_end_tmp = zeros(size(stimData,1),1);
    idx_count = 0;
    idx_error = [];
    for IDX = 1:size(idx_trial_start,1)-1
        idx_count = idx_count+1;
        idx_trial_start_tmp(idx_count) = idx_trial_start(IDX);
        idx_trial_end_tmp(idx_count) = idx_trial_end(IDX);
        if idx_trial_diff(IDX) > 150
            idx_count = idx_count+1;
            idx_error = [idx_error;idx_count-1;idx_count];
            idx_trial_start_tmp(idx_count)=idx_trial_start(IDX)+120;
            idx_trial_end_tmp(idx_count)=idx_trial_end(IDX)+120;
        end
        
    end
    idx_trial_start_tmp(idx_count+1) = idx_trial_start(IDX+1);
    idx_trial_end_tmp(idx_count+1) = idx_trial_end(IDX+1);
    idx_trial_start = idx_trial_start_tmp;
    idx_trial_end = idx_trial_end_tmp;
    frame_trial_start = peripherals.FrameNumber(idx_trial_start);
    frame_trial_end = peripherals.FrameNumber(idx_trial_end);
    temp_tbl.start_time = peripherals.sglxTime(idx_trial_start);
    temp_tbl.end_time = peripherals.sglxTime(idx_trial_end);
    nTrials = length(temp_tbl.start_time);
    % get trial start/end times from photodiode signal
    pd_thresh_up = 400;
    pd_thresh_down = 100;
    for itrial =1:nTrials
        % find frame where Bonsai asked to change quad
        frameidx_start = find(photodiode.FrameNumber <= frame_trial_start(itrial),1,'last');
        % find next index where the photodiode detected a quad change
        temp_idx = find(photodiode.Photodiode(frameidx_start:end)>pd_thresh_up,1,'first');
        idx_start = frameidx_start+temp_idx-1;
        % convert to spike GLX time
        pdstart(itrial) = photodiode.sglxTime(idx_start);
        % to get photodiode down we need to do some median filtering to
        % remove some awkward down swings when the quad is white
        frameidx_end = find(photodiode.FrameNumber >= frame_trial_end(itrial),1,'first');
        temp_idx = find(smoothdata(photodiode.Photodiode(frameidx_end:end),'movmedian', 5)<pd_thresh_down,1,'first');
        idx_end = frameidx_end+temp_idx-1;
        pdend(itrial) = photodiode.sglxTime(idx_end);
    end

    photodiodeData.stim_on.sglxTime = pdstart';
    photodiodeData.stim_off.sglxTime = pdend';
    wheelData.staticgrating_idx_error = idx_error;
end

% Step 5: process wheel data (skipping this just save all peripheral data)
% wheelData.pos = peripherals.Wheel;
% wheelData.Time = peripherals.sglxTime;
% wheelData = import_processWheelEphys(wheelData,'gaussian',12);
% wheelData.sglxTime = wheelData.Time;

% Step 6: extract stimulus times
switch(lower(StimulusName))
    case {'sparsenoise','sparsenoise_fullscreen','checkerboard'}
        stimTimes = sort(cat(1,photodiodeData.stim_on.sglxTime,photodiodeData.stim_off.sglxTime),'ascend');
%         fprintf('\n\tUsing both photodiode upward and downward transitions as timestamps for stimulus timing')
    case {'staticgratings','staticgratings_short','staticgratings_long'}
        stimTimes = sort(cat(1,photodiodeData.stim_on.sglxTime),'ascend');
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

