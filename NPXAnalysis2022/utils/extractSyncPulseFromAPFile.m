% Program to try to extract the sync pulse efficiently from AP files in NP1 or NP2 
% data sets
%  Inputs - imec file structure, sync pulse channel number (default is that in the AP file but cna be overridden)
%  Output - times of onset and offset of the sync pulse
% History: SGS 5th April 2022 Wrote it
%
function syncTimes_ephys = extractSyncPulseFromAPFile(imec, syncChannelNum)
if ~exist('syncChannelNum','var') || isempty(syncChannelNum)
    syncChannelNum = imec.syncChannelId;
    syncChannelIndex = imec.syncChannelIndex;
else
    syncChannelIndex = find(imec.channelIds == syncChannelNum);
end
% The AP file is stored in the variable imec.pathAP
% We should be able to extract this from the memmap by indexing the
% appropriate value; we only need to extract a subset of values
downsampleRate = fix(imec.fsAP/1000);    % should be 30x - ie transform AP sampling frequency to 1000 Hz
mmap = imec.memmapAP_full(); % get AP data map

% **Note that resultant indexes will be at the downsampled rate!
value = mmap.Data.x(syncChannelIndex, 1:downsampleRate:size(mmap.Data.x,2)); % extract at downsampled rate 
% Find upswings in sync pulse
vec_idx = find(diff(value)>=0.5);
% add 1 to vec_idx to compensate for diff and convert from samples to s 
syncTimes_ephys.on = (vec_idx+1)./(imec.fsAP/downsampleRate); 
% Find downswings in sync pulse
vec_idx = find(diff(value)<=-0.5); 
% add 1 to vec_idx to compensate for diff and convert from samples to s 
syncTimes_ephys.off = (vec_idx+1)./(imec.fsAP/downsampleRate); 

