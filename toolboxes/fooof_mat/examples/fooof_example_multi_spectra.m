%% FOOOF Matlab Wrapper Example - Multiple PSDs
%
% This example computes example power spectra models on a group of
% power spectra, and prints out the results.
%

%% Run Example
pyenv("ExecutionMode","OutOfProcess")
% Load data
load('data/ch_dat_one.mat');
load('data/ch_dat_two.mat');

% Combine into a multi-channel data matrix
chs_dat = [ch_dat_one; ch_dat_two]';

% Calculate power spectra with Welch's method
[psds, freqs] = pwelch(chs_dat, 500, [], [], s_rate);

% Transpose, to make inputs row vectors
freqs = freqs';

% FOOOF settings
settings = struct();
f_range = [1, 30];

% Run FOOOF across a group of power spectra
fooof_results = fooof_group(freqs, psds, f_range, settings);

% Check out the FOOOF Results
fooof_results
fooof_results(2).aperiodic_params(2)
%aperiodic_params_: a list of aperiodic parameters, stored as [Offset, (Knee), Exponent]
%peak_params_: all periodic parameters, where each row is a peak, as [CF, PW, BW]
%r_squared_: the r-squared of the model, as compared to the original data
%error_: the error of the model, as compared to the original data