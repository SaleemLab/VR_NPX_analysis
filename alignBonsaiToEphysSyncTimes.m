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