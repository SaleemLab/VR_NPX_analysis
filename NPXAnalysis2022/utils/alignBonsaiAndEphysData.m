% function [stimData, peripherals,eyeData] = alignBonsaiAndEphysData(imec,stimData,peripherals,eyeData)
%   Function to align the bonsai stimulus, peripherals and pupil recordings with the
%   ephys recordings. 
%
% This function uses the LFP channels for NP1 data. Would like to make it
%   generic and able to use full spectrum data from NP2 files.
%
% Inputs: syncTimes_ephys (from imec file), stimData, peripherals and
%               eyeData structures (and photodiode structures for newer files)
% Outputs: 'sglxTime' fields are added to the stimData, peripherals, and eyeData structures
%               (and photodiode structures) which are otherwise unchanged
%
% History: SGS March 2022: Wrote it
%          SGS 4/4/2022:   Shifted imec call to importAndAlignBonsaiLogs
%          SGS 10/4/2022:  Added photodiode 
%          SGS 12/4/2022:  Restructured so only one transparent method, and
%                          included initial parse to improve xcorr; also
%                          forced interp1 to linear and allowed
%                          extrapolation
function [stimData,peripherals,eyeData,photodiode] = alignBonsaiAndEphysData(syncTimes_ephys,stimData,peripherals,eyeData,photodiode)

%%%%%%
% For wheel data and convert to spike GLX time-base
[peripherals] = alignBonsaiToEphysSyncTimes(peripherals,syncTimes_ephys);

if ~isempty(eyeData)
    [eyeData] = alignBonsaiToEphysSyncTimes(eyeData,syncTimes_ephys);
end

%%%%%%
% For photodiode and convert to spike GLX time-base
if exist('photodiode','var') && ~isempty(photodiode)
    [photodiode] = alignBonsaiToEphysSyncTimes(photodiode,syncTimes_ephys);
else 
    photodiode = [];
end

%%%%%
% For stimData table - if it has a 'Timestamp' field (loosely defined)
ts_column = find(strcmp(stimData.Properties.VariableNames,'Time')); % replaced 'timestamp' with 'time' 10th April 2022
if ~isempty(ts_column)
    td = cat(1,stimData{:,ts_column(1)});
    if ~isempty(photodiode) % If photodiode is present take from that....
        tp = find(~isnan(photodiode.Time) & ~isnan(photodiode.sglxTime));
        t1 = unique(photodiode.Time(tp));
        t2 = unique(photodiode.sglxTime(tp));
    else % ...otherwise take from eye
        tp = find(~isnan(eyeData.Time) & ~isnan(eyeData.sglxTime));
        t1 = eyeData.Time(tp);
        t2 = eyeData.sglxTime(tp);
    end
    stimTimes = interp1(t1(:)', t2(:)',td(:)')';
    stimData.sglxTime = stimTimes;
end
end

function [bonsai_data] = alignBonsaiToEphysSyncTimes(bonsai_data,syncTimes_ephys)

bonsai_idx = find(diff(bonsai_data.Sync)==1);
syncTimes_bonsai = bonsai_data.Time(bonsai_idx+1)./1000; % add 1 to idx to compensate for diff and convert from ms to s
% Especially in early recordings in a session there may be significant
% delay between the start of Bonsai measurements and the start of ephys
% measurements. This can make xcorr less reliable as there are a large
% number of potential correspondences
% We therefore first attempt to align the overall sync pulse traces ...
[r, lags] = xcorr(diff(syncTimes_bonsai), diff(syncTimes_ephys));
[~, joint_idx] = max(r);
best_lag = lags(joint_idx);
% ...then remove the syncpulses from bonsai trace that occur before the
% optimum lag
if best_lag > 0
    syncTimes_bonsai = syncTimes_bonsai(syncTimes_bonsai>syncTimes_bonsai(best_lag));
end
% ...and rerun the xcorr
[r, lags] = xcorr(diff(syncTimes_bonsai), diff(syncTimes_ephys));
[~, joint_idx] = max(r);
best_lag = lags(joint_idx);

% check which async pulse sequence is available first (ie. which was turned
% on first)
if best_lag < 0
    nSyncOffset = -best_lag+1;
    t_npix = syncTimes_ephys(nSyncOffset:end); % sync sglx times
    if length(t_npix) > length(syncTimes_bonsai)
        t_npix = t_npix(1:length(syncTimes_bonsai));
    end
    t_bonsai = syncTimes_bonsai(1:numel(t_npix)); %sync bonsai times
else
    nSyncOffset = best_lag+1;
    t_npix = syncTimes_ephys(1:end);
    if length(t_npix) > length(syncTimes_bonsai)
        t_npix = t_npix(1:length(syncTimes_bonsai));
    end
    t_bonsai = syncTimes_bonsai(nSyncOffset:numel(t_npix)+nSyncOffset-1);
end
bonsai_data.sglxTime = interp1(t_bonsai, t_npix, bonsai_data.Time./1000,'linear','extrap');

end
