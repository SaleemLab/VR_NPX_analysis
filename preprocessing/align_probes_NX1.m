% function [stimData, peripherals,eyeData] = alignBonsaiAndEphysData(imec,stimData,peripherals,eyeData)
% Function to align probe 1 and probe 2 
% 1. First TTL ON recorded by probe 1 and probe 2 AP sample band to determine the offset
% 2. Use last and first TTL ON to determine the scaling factor.
% -3. calculate the probe 2 AP sample number in terms of probe 1 AP sample number
% Automatically, LFP band sample can be corrected because it is percisely
% acquired every 12th spike band sample 

% https://open-ephys.github.io/gui-docs/Tutorials/Data-Synchronization.html
% This is the putative, crude alignment method. But since the average
% offset between two probes is between 1 and 2 ms, this method may be
% sufficient. 

function [stimData,peripherals,eyeData,photodiode] = align_probes_NX1(options)

for nprobe = 1:length(options)
    tic
    DIR = dir(fullfile(options(nprobe).EPHYS_DATAPATH,'*syncpulseTimes.mat'))
%     DIR = [];
    [~,fname] = fileparts(options(nprobe).EPHYS_DATAPATH);

    if ~isempty(DIR)
        load(fullfile(options(nprobe).EPHYS_DATAPATH,DIR.name))
        
%         syncTimes_ephys = syncTimes_ephys.on; % use upswings currently
    else
        [AP_FILE,LF_FILE] = findImecBinFile(options(nprobe).EPHYS_DATAPATH);
        %     FILE_TO_USE = AP_FILE;
        binpath = fullfile(options(nprobe).EPHYS_DATAPATH,AP_FILE);
        %     binpath = fullfile(options.EPHYS_DATAPATH,LF_FILE);
        syncTimes_ephys = SGLXextractSyncPulseFromBinFile(binpath);
        parseDate = date;
        
        save(fullfile(options(nprobe).EPHYS_DATAPATH,[fname,'_syncpulseTimes.mat']),'syncTimes_ephys','parseDate');
%         syncTimes_ephys = syncTimes_ephys.on; % use upswings currently
    end

    syncTimes_ephys_multi(nprobe).probe = syncTimes_ephys;
    toc
end
% sync_offset = syncTimes_ephys_multi(1).probe.on - syncTimes_ephys_multi(2).probe.on;
% sync_offset(1);

vec_idx1 = find(diff(syncTimes_ephys_multi(1).probe.Sync)>=0.5);
vec_idx2 = find(diff(syncTimes_ephys_multi(2).probe.Sync)>=0.5);
scaling_factor = (vec_idx1(end)-vec_idx1(1))/(vec_idx2(end)-vec_idx2(1));

% aligned_AP_sample_number = 1:length(syncTimes_ephys_multi(2).probe.Sync);
aligned_AP_sample_number= round(([1:length(syncTimes_ephys_multi(2).probe.Sync)] - vec_idx2(1)) * scaling_factor + vec_idx1(1));

% Probe 2 AP data sample in terms of Probe 1 AP data
save(fullfile(options(2).EPHYS_DATAPATH,[fname,'_aligned_AP_sample_number.mat']),'scaling_factor','aligned_AP_sample_number','parseDate');

end
