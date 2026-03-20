function [bonsai_data] = alignBonsaiToEphysSyncPulse(bonsai_data,syncPulse_ephys,syncPulse_ephysTimes)

[bonsai_tvec,ia,ic] = unique((bonsai_data.Time-bonsai_data.Time(1))./1000);
ephys_SR = 1./mean((diff(syncPulse_ephysTimes)));
% Resample bonsai_data.Sync to match syncTimes_ephys
bonsai_interp = interp1(bonsai_tvec, bonsai_data.Sync(ia), syncPulse_ephysTimes, 'previous','extrap');
bonsai_interp(isnan(bonsai_interp)) = 0;

[r, lags] = xcorr(normalize(bonsai_interp), normalize(syncPulse_ephys));
[~, joint_idx] = max(r);
best_lag = lags(joint_idx);


bonsai_data.sglxTime = (bonsai_data.Time-bonsai_data.Time(1))./1000-best_lag*(1/ephys_SR);
% figure;plot(bonsai_data.sglxTime,bonsai_data.Photodiode);hold on;plot(bonsai_data.Time./1000-(bonsai_data.Time(1)./1000-bonsai_data.sglxTime(1)),bonsai_data.Photodiode)
end