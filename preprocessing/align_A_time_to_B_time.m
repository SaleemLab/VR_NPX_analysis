function A_time_in_B = align_A_time_to_B_time(A_signal,A_time,B_signal,B_time)
% whole signal cross correlation between A and B but align A to B signals
% in B signal time


[A_tvec,ia,ic] = unique((A_time-A_time(1)));
B_SR = 1./mean((diff(B_time)));
% Resample bonsai_data.Sync to match syncTimes_ephys
A_interp = interp1(A_tvec, A_signal(ia), B_time, 'nearest');
A_interp(isnan(A_interp)) = 0;





[r, lags] = xcorr(normalize(A_interp), normalize(B_signal));
[~, joint_idx] = max(r);
best_lag = lags(joint_idx);


A_time_in_B = (A_time-A_time(1))-best_lag*(1/B_SR);
% figure;plot(bonsai_data.sglxTime,bonsai_data.Photodiode);hold on;plot(bonsai_data.Time./1000-(bonsai_data.Time(1)./1000-bonsai_data.sglxTime(1)),bonsai_data.Photodiode)
end