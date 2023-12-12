% Program to align two sync pulses from logs (or from log and npix)
%   In this function, the timebase in log1 (log1.Time)
%   is shifted to the timebase in log2 (log2.Time),
%   and the resultant timebase is appended to log1 as log1.realignedTime
%
%   If you would like to align log timebase to npix timebase, you should therefore
%   include the npix as the second input
% SGS 03/04/2022
function [log1,log2] = alignSyncPulses(log1,log2)

% Find upswings in sync pulses in log1 and log2
idx1 = find(diff(log1.Sync)==1);
syncTimes1 = log1.Time(idx1+1); % add 1 to idx to compensate for diff 

idx2 = find(diff(log2.Sync)==1);
syncTimes2 = log2.Time(idx1+1); % add 1 to idx to compensate for diff 

% Find correlation between two 
[r, lags] = xcorr(diff(syncTimes1), diff(syncTimes2));
[~, joint_idx] = max(r);
best_lag = lags(joint_idx);

% check which async pulse sequence is available first (ie. which was turned
% on first)
if best_lag < 0
    nSyncOffset = -best_lag+1;
    t2 = syncTimes2(nSyncOffset:end); % sync sglx times
    if length(t2) > length(syncTimes1)
        t2 = t2(1:length(syncTimes1));
    end
    t1 = syncTimes1(1:numel(t2)); %sync bonsai times
else
    nSyncOffset = best_lag+1;
    t2 = syncTimes2(1:end);
    if length(t2) > length(syncTimes1)
        t2 = t2(1:length(syncTimes1));
    end
    t1 = syncTimes1(nSyncOffset:numel(t2)+nSyncOffset-1);
end
log1.realignedTime = interp1(t1, t2, log1.Time,'linear','extrap');
