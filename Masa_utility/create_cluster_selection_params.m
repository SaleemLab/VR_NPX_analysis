function params = create_cluster_selection_params

% default values of the params
params.isi_violations_ratio = @(x) x<=0.1;
params.amplitude_cutoff = @(x) x<=0.1; %0.01 if strict and removes lots of units
params.amplitude_median = @(x) x>50; %IBL 50
params.sliding_rp_violation = @(x) x<=0.1; % chosen as 10% at IBL
params.drift_ptp = @(x) x<5;
params.drift_std = @(x) x<1;
params.drift_mad = @(x) x<1;
params.num_negative_peaks = @(x) x<=1;
params.num_positive_peaks = @(x) x<=2;
params.peak_to_valley = @(x) x<=0.0008 & x>= 0.0002;
params.amplitude_cv_median = @(x) x<=0.7;
params.amplitude_cv_range = @(x) x<=0.7;

end