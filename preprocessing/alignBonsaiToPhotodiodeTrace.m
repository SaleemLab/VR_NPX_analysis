function [bonsai_data] = alignBonsaiToPhotodiodeTrace(bonsai_data,photodiode_signal,photodiode_tvec)
[bonsai_tvec,ia,ic] = unique(bonsai_data.sglxTime);
ephys_SR = 1./mean((diff(photodiode_tvec)));
% Resample bonsai_data.Sync to match syncTimes_ephys
bonsai_interp = interp1(bonsai_tvec, bonsai_data.Photodiode(ia), photodiode_tvec, 'previous','extrap');
bonsai_interp(isnan(bonsai_interp)) = 0;

[r, lags] = xcorr(normalize(bonsai_interp), normalize(photodiode_signal));
[~, joint_idx] = max(r);
best_lag = lags(joint_idx);


bonsai_data.corrected_sglxTime = (bonsai_data.sglxTime)-best_lag*(1/ephys_SR);
% figure;plot(bonsai_data.sglxTime,bonsai_data.Photodiode);hold on;plot(bonsai_data.Time./1000-(bonsai_data.Time(1)./1000-bonsai_data.sglxTime(1)),bonsai_data.Photodiode)
end