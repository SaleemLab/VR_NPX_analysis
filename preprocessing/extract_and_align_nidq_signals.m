function extract_and_align_nidq_signals(options)
% This code extracts  nidq signals then aligns nidq signals to Ephys signal
DIR = dir(fullfile(options.EPHYS_DATAPATH,'..','*_NidqTimes.mat'))
if isempty(DIR)
    
    td = dir(fullfile(options.EPHYS_DATAPATH,'..','*nidq.bin'));
    nidq_bin = fullfile(td.folder,td.name);
    td = dir(fullfile(options.EPHYS_DATAPATH,'..','*nidq.meta'));
    nidq_meta = fullfile(td.folder,td.name);
    
    binpath=nidq_bin;
    % Extract path and filename
    [fpath,binName,extension] = fileparts(binpath);
    binName = [binName,extension];
    % Parse the corresponding metafile
    meta = ReadMeta(binpath);
    nChan = str2double(meta.nSavedChans);
    nSamp = floor(10.0 * SampRate(meta)); % Read in 10s interval chunks - could improve this
    nFileSamp = str2double(meta.fileSizeBytes) / (2 * nChan);
    
    % For a digital channel: read this digital word dw in the saved file
    % (1-based). For imec data there is never more than one saved digital word.
    dw = 1;
    % Read these lines in dw (0-based).
    % For 3B2 imec data: the sync pulse is stored in line 6.
    % May be 1 or more line indices.
    dLineList = [0,1,6];
    
    
    Nidq.tvec = (1:nFileSamp)./SampRate(meta); % nidq timevec
    Nidq.bonsai_sync = zeros(1,nFileSamp);% Sync pulse to sync Bonsai and nidq/Ephys
    Nidq.nidq_sync = zeros(1,nFileSamp); % Sync pulse to sync nidq and Ephys (last Sync channel in ephys bin)
    Nidq.photodiode = zeros(1,nFileSamp);% photodiode
    Nidq.wheel_forward = zeros(1,nFileSamp);% Wheel ticks forward
    Nidq.wheel_backward = zeros(1,nFileSamp);% Wheel ticks backward
    Nidq.lick_left = zeros(1,nFileSamp);%
    Nidq.lick_right = zeros(1,nFileSamp);%
    Nidq.lick_left = zeros(1,nFileSamp);%
    Nidq.lick_right = zeros(1,nFileSamp);%
    
    nIter = ceil(nFileSamp/nSamp);
    for thisFilePart = 0:nIter-1 %load 10 seconds segment of the data
        thisSampleStart = thisFilePart*nSamp;
        dataArray = ReadBin(thisSampleStart, nSamp, meta, binName, fpath);
        
        Nidq.bonsai_sync(:,thisSampleStart+1:thisSampleStart+size(dataArray,2)) = dataArray(1,:);% Sync pulse to sync Bonsai and nidq/Ephys
        Nidq.nidq_sync(:,thisSampleStart+1:thisSampleStart+size(dataArray,2)) = dataArray(end,:);% % Sync pulse to sync nidq and Ephys (last Sync channel in ephys bin)
        Nidq.photodiode(:,thisSampleStart+1:thisSampleStart+size(dataArray,2)) = dataArray(2,:);% photodiode
        Nidq.wheel_forward(:,thisSampleStart+1:thisSampleStart+size(dataArray,2)) = dataArray(3,:);% Wheel ticks forward
        Nidq.wheel_backward(:,thisSampleStart+1:thisSampleStart+size(dataArray,2)) = dataArray(4,:);% Wheel ticks backward
        Nidq.reward_left(:,thisSampleStart+1:thisSampleStart+size(dataArray,2)) = dataArray(7,:);%
        Nidq.reward_right(:,thisSampleStart+1:thisSampleStart+size(dataArray,2)) = dataArray(8,:);%
        Nidq.lick_left(:,thisSampleStart+1:thisSampleStart+size(dataArray,2)) = dataArray(5,:);%
        Nidq.lick_right(:,thisSampleStart+1:thisSampleStart+size(dataArray,2)) = dataArray(6,:);%
    end
    
    
    %%% Regular sync pulse (ephys and Nidq)
    % Find upswings in sync pulse
    vec_idx = find(diff(Nidq.nidq_sync)>=mean(Nidq.nidq_sync));
    % add 1 to vec_idx to compensate for diff and convert from samples to s
    Nidq.nidq_sync_on = (vec_idx+1)./SampRate(meta);
    % Find downswings in sync pulse
    vec_idx = find(diff(Nidq.nidq_sync)<=-mean(Nidq.nidq_sync));
    % add 1 to vec_idx to compensate for diff and convert from samples to s
    Nidq.nidq_sync_off = (vec_idx+1)./SampRate(meta);
    
    %%% Extract Ephys Sync pulse from last channel (ephys and nidq)
    fname = erase(td.name,'.nidq.meta');
    DIR = dir(fullfile(options.EPHYS_DATAPATH,'..','*syncpulseTimes.mat'))
    %     DIR = [];
    if ~isempty(DIR)
        load(fullfile(options.EPHYS_DATAPATH,'..',DIR.name))
    else
        [ap_file,lf_file] = findImecBinFile(options.EPHYS_DATAPATH);
        if ~isempty(lf_file)
            binpath = fullfile(options.EPHYS_DATAPATH,lf_file);
        else
            fprintf('\n'); warning('Extracting sync-pulses from AP file...this takes about 5 minutes or more'); fprintf('\n');
            binpath = fullfile(options.EPHYS_DATAPATH,ap_file);
        end
        syncTimes_ephys = SGLXextractSyncPulseFromBinFile(binpath);
        parseDate = date;
        save(fullfile(options.EPHYS_DATAPATH,'..',[fname,'_syncpulseTimes.mat']),'syncTimes_ephys','parseDate','-v7.3'); % used to be saved inside the ephys probe folder, now just same place where niqd is saved
    end
    
    
    %  [~,fname] = fileparts(options.EPHYS_DATAPATH);
    % save(fullfile(options.EPHYS_DATAPATH,[fname,'_syncpulseTimes.mat']),'syncTimes_ephys','parseDate');
    [Nidq] = alignNiqdToEphysSyncTimes(Nidq,syncTimes_ephys.on);
    
    
    %%%
    % Store the spikeglx based event times
    %%% Regular sync pulse (ephys and Nidq)
    % Find upswings in sync pulse
    vec_idx = find(diff(Nidq.nidq_sync)>=mean(Nidq.nidq_sync));
    % add 1 to vec_idx to compensate for diff and convert from samples to s
    Nidq.nidq_sync_on = Nidq.sglxTime(vec_idx+1);
    % Find downswings in sync pulse
    vec_idx = find(diff(Nidq.nidq_sync)<=-mean(Nidq.nidq_sync));
    % add 1 to vec_idx to compensate for diff and convert from samples to s
    Nidq.nidq_sync_off = Nidq.sglxTime(vec_idx+1);
    
    %%% Irregular sync pulse (bonsai and Nidq)
    % Find upswings in sync pulse
    vec_idx = find(diff(Nidq.bonsai_sync)>=mean(Nidq.bonsai_sync));
    % add 1 to vec_idx to compensate for diff and convert from samples to s
    Nidq.bonsai_sync_on = Nidq.sglxTime(vec_idx+1);
    % Find downswings in sync pulse
    vec_idx = find(diff(Nidq.bonsai_sync)<=-mean(Nidq.bonsai_sync));
    % add 1 to vec_idx to compensate for diff and convert from samples to s
    Nidq.bonsai_sync_off = Nidq.sglxTime(vec_idx+1);
    
    
    %%% Photodiode (bonsai and Nidq)
    % Find upswings in sync pulse
    vec_idx = find(Nidq.photodiode>=mean(Nidq.photodiode));
    Nidq.photodiode_on = Nidq.sglxTime(vec_idx);
    % Find downswings in sync pulse
    vec_idx = find(Nidq.photodiode<=mean(Nidq.photodiode));
    Nidq.photodiode_off = Nidq.sglxTime(vec_idx);
    parseDate = date;
    % [~,fname] = fileparts(options.EPHYS_DATAPATH);
    
    
    %%% Lick Left (bonsai and Nidq) (for lick it is positive going to 0 due to blockage of IR)
    % Find upswings in sync pulse
    vec_idx = find(diff(Nidq.lick_left)<=-10000);
    % add 1 to vec_idx to compensate for diff and convert from samples to s
    Nidq.lick_left_on = Nidq.sglxTime(vec_idx+1);
    % Find downswings in sync pulse
    vec_idx = find(diff(Nidq.lick_left)>=10000);
    % add 1 to vec_idx to compensate for diff and convert from samples to s
    Nidq.lick_left_off = Nidq.sglxTime(vec_idx+1);
    
    %%% Lick Right (bonsai and Nidq)
    % Find upswings in sync pulse
    vec_idx = find(diff(Nidq.lick_right)<=-10000);
    % add 1 to vec_idx to compensate for diff and convert from samples to s
    Nidq.lick_right_on = Nidq.sglxTime(vec_idx+1);
    % Find downswings in sync pulse
    vec_idx = find(diff(Nidq.lick_right)>=10000);
    % add 1 to vec_idx to compensate for diff and convert from samples to s
    Nidq.lick_right_off = Nidq.sglxTime(vec_idx+1);
    
    
    %%% Reward Left (bonsai and Nidq)
    % Find upswings in sync pulse
    vec_idx = find(diff(Nidq.reward_left)>=10000);
    % add 1 to vec_idx to compensate for diff and convert from samples to s
    Nidq.reward_left_on = Nidq.sglxTime(vec_idx+1);
    % Find downswings in sync pulse
    vec_idx = find(diff(Nidq.reward_left)<=-10000);
    % add 1 to vec_idx to compensate for diff and convert from samples to s
    Nidq.reward_left_off = Nidq.sglxTime(vec_idx+1);
    
    %%% Reward Left (bonsai and Nidq)
    % Find upswings in sync pulse
    vec_idx = find(diff(Nidq.reward_right)>=10000);
    % add 1 to vec_idx to compensate for diff and convert from samples to s
    Nidq.reward_right_on = Nidq.sglxTime(vec_idx+1);
    % Find downswings in sync pulse
    vec_idx = find(diff(Nidq.reward_right)<=-10000);
    % add 1 to vec_idx to compensate for diff and convert from samples to s
    Nidq.reward_right_off = Nidq.sglxTime(vec_idx+1);
    parseDate = date;
    % [~,fname] = fileparts(options.EPHYS_DATAPATH);
    
    save(fullfile(options.EPHYS_DATAPATH,'..',[fname,'_NidqTimes.mat']),'Nidq','parseDate','-v7.3');
end

end


