function params = create_cluster_selection_params(varargin)

% Default values
p = inputParser;
addParameter(p,'sorting_option','spikeinterface',@isstr) % Select channels for analysis (default is all the channles or at least all channels loaded (e.g. only from one cloumn))

% assign parameters (either defaults or given)
parse(p,varargin{:});
sorting_option = p.Results.sorting_option;

if contains(sorting_option,'spikeinterface')
    % default values of the params
%     params.isi_violations_ratio = @(x) x<=0.1;
    params.amplitude_cutoff = @(x) x<=0.1; %0.01 if strict and removes lots of units
    params.amplitude_median = @(x) abs(x)>50; %IBL 50
    params.sliding_rp_violation = @(x) x<=0.1; % chosen as 10% at IBL
    % params.drift_ptp = @(x) x<5;
    % params.drift_std = @(x) x<1;
    % params.drift_mad = @(x) x<1;
    params.num_negative_peaks = @(x) x<=1;
    params.num_positive_peaks = @(x) x<=2;
    params.peak_to_valley = @(x) x<=0.0008 & x>= 0.0002;
    params.amplitude_cv_median = @(x) x<=0.7;
    params.amplitude_cv_range = @(x) x<=0.7;

elseif contains(sorting_option,'old')
    params.amplitude_cutoff = @(x) x<=0.1; %0.01 if strict and removes lots of units
    params.amplitude = @(x) x>50; %IBL 50
    params.isi_viol = @(x) x<=0.1;
    params.duration = @(x) x<=0.8 & x>= 0.2;
elseif contains(sorting_option,'masa') 
    params.amplitude_cutoff = @(x) x<=0.1; %0.01 if strict and removes lots of units
    params.amplitude_median = @(x) abs(x)>20; %IBL 50 but bombcell 20
    params.sliding_rp_violation = @(x) x<=0.1; % chosen as 10% at IBL
    params.num_negative_peaks = @(x) x<=1;
    params.num_positive_peaks = @(x) x<=2;
    params.peak_to_valley = @(x) x<=0.0008 & x>= 0.0001; % Peters et al 2022
    params.amplitude_cv_median = @(x) x<=0.7;
    params.amplitude_cv_range = @(x) x<=0.7;
    params.firing_rate = @(x) x>= 0.01;% Due to sleep/wave difference, we
    %can have silent neuron during sleep or awake
    % params.presence_ratio = @(x) x>= 0.3;% Due to sleep/wave difference,
    %we can have silent neuron during sleep or awake during sleep or awake
    params.sd_ratio = @(x) x<=4;
    params.snr = @(x) x>= 2;
    params.half_width = @(x) x<= 0.0003;
else

end

end