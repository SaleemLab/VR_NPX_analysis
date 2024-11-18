function [Nidq_data] = alignNiqdToEphysSyncTimes(Nidq_data,syncTimes_ephys)

% used to align nidq to ephys regular sync pulse
syncTimes_nidq = Nidq_data.nidq_sync_on;

% Especially in early recordings in a session there may be significant
% delay between the start of Bonsai measurements and the start of ephys
% measurements. This can make xcorr less reliable as there are a large
% number of potential correspondences
% We therefore first attempt to align the overall sync pulse traces ...
[r, lags] = xcorr(diff(syncTimes_nidq), diff(syncTimes_ephys));
[~, joint_idx] = max(r);
best_lag = lags(joint_idx);
% ...then remove the syncpulses from bonsai trace that occur before the
% optimum lag
if best_lag > 0
    syncTimes_nidq = syncTimes_nidq(syncTimes_nidq>syncTimes_nidq(best_lag));
end
% ...and rerun the xcorr
[r, lags] = xcorr(diff(syncTimes_nidq), diff(syncTimes_ephys));
[~, joint_idx] = max(r);
best_lag = lags(joint_idx);

% check which async pulse sequence is available first (ie. which was turned
% on first)
if best_lag < 0
    nSyncOffset = -best_lag+1;
    t_npix = syncTimes_ephys(nSyncOffset:end); % sync sglx times
    if length(t_npix) > length(syncTimes_nidq)
        t_npix = t_npix(1:length(syncTimes_nidq));
    end
    t_bonsai = syncTimes_nidq(1:numel(t_npix)); %sync bonsai times
else
    nSyncOffset = best_lag+1;
    t_npix = syncTimes_ephys(1:end);
    if length(t_npix) > length(syncTimes_nidq)
        t_npix = t_npix(1:length(syncTimes_nidq));
    end
    t_bonsai = syncTimes_nidq(nSyncOffset:numel(t_npix)+nSyncOffset-1);
end
Nidq_data.sglxTime = interp1(t_bonsai, t_npix, Nidq_data.tvec,'linear','extrap');

end