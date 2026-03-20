% Originally based on alignBonsaiAndEphysData
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


% Now input rather than syncTimes_ephys (swing time) -> now use syncPulse_ephys (sync pulse trace) and
% syncPulse_ephysTimes (tvec)
%function [stimData,peripherals,eyeData,photodiode] = alignAllBonsaiAndEphysData(syncPulse_ephys,syncPulse_ephysTimes,stimData,peripherals,eyeData,photodiode)

% 
function [stimData,peripherals,eyeData,photodiode] = alignAllBonsaiAndEphysData(syncPulse_ephys,syncPulse_ephysTimes,stimData,peripherals,eyeData,photodiode)

%%%%%%
% For wheel data and convert to spike GLX time-base
% [peripherals] = alignBonsaiToEphysSyncTimes(peripherals,syncTimes_ephys);
[peripherals] = alignBonsaiToEphysSyncPulse(peripherals,syncPulse_ephys,syncPulse_ephysTimes);
if ~isempty(eyeData)
%     [eyeData] = alignBonsaiToEphysSyncTimes(eyeData,syncTimes_ephys);
    [eyeData] = alignBonsaiToEphysSyncPulse(eyeData,syncPulse_ephys,syncPulse_ephysTimes);
end

%%%%%%
% For photodiode and convert to spike GLX time-base
if exist('photodiode','var') && ~isempty(photodiode)
%     [photodiode] = alignBonsaiToEphysSyncTimes(photodiode,syncTimes_ephys);
    [photodiode] = alignBonsaiToEphysSyncPulse(photodiode,syncPulse_ephys,syncPulse_ephysTimes);
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
