function [lags, corr_mean, corr_CI, shuffle_mean, shuffle_CI] = ...
    xcorr_bias(HPC_bias, V1_bias, max_lag_bins, num_shuffles, num_boot, varargin)
%XCORR_BIAS Fast cross-correlation between two neural bias signals (e.g., HPC and V1)
%
%   [lags, corr_mean, corr_CI, shuffle_mean, shuffle_CI] = ...
%       xcorr_bias(HPC_bias, V1_bias, max_lag_bins, num_shuffles, num_boot, 'use_gpu', true)
%
%   Computes lagged Pearson cross-correlations between z-scored bias time series across ripple events.
%   Fast implementation using dot products and optional GPU acceleration. Includes bootstrap and shuffle CIs.
%
%   Inputs:
%     HPC_bias      - [nEvents x nTimeBins] matrix of z-scored bias from HPC
%     V1_bias       - [nEvents x nTimeBins] matrix of z-scored bias from V1
%     max_lag_bins  - Maximum lag in bins (positive and negative)
%     num_shuffles  - Number of circular shuffle iterations for null distribution
%     num_boot      - Number of bootstrap resampling iterations
%     'use_gpu'     - (Optional) Boolean flag to use GPU if available [default: false]
%
%   Outputs:
%     lags          - Vector of lag values (in bins), e.g., -max_lag:1:max_lag
%     corr_mean     - Mean real cross-correlation across trials (1 x num_lags)
%     corr_CI       - [2 x num_lags] Bootstrap-based 95% confidence interval (rows: [2.5th, 97.5th] percentiles)
%     shuffle_mean  - Mean cross-correlation from shuffled data
%     shuffle_CI    - [2 x num_lags] 95% CI of the shuffle distribution

% Masahiro Takigawa 2024
% ------------------------------
% Parse optional input arguments
% ------------------------------
p = inputParser;
addParameter(p, 'use_gpu', false, @(x) islogical(x) || isnumeric(x));
parse(p, varargin{:});
use_gpu = logical(p.Results.use_gpu);

% ------------------------------
% Setup
% ------------------------------
[N_events, T] = size(HPC_bias);
lags = -max_lag_bins:max_lag_bins;
num_lags = length(lags);

% ------------------------------
% Z-score each trial (row-wise)
% ------------------------------
HPC_z = zscore(HPC_bias, 0, 2);
V1_z = zscore(V1_bias, 0, 2);

if use_gpu
    HPC_z = gpuArray(HPC_z);
    V1_z = gpuArray(V1_z);
end

% ------------------------------
% Real cross-correlation: trial-wise dot product across lags
% ------------------------------
corr_mat = zeros(N_events, num_lags);
for li = 1:num_lags
    lag = lags(li);
    if lag < 0
        h = HPC_z(:, 1:end+lag);
        v = V1_z(:, 1-lag:end);
    elseif lag > 0
        h = HPC_z(:, 1+lag:end);
        v = V1_z(:, 1:end-lag);
    else
        h = HPC_z;
        v = V1_z;
    end
    % Pearson via dot product (z-scored data)
    corr_mat(:, li) = sum(h .* v, 2) ./ (size(h,2) - 1);
end

corr_mean = mean(corr_mat, 1);
if use_gpu
    corr_mean = gather(corr_mean);  % bring back to CPU
end

% ------------------------------
% Bootstrap confidence intervals across events
% ------------------------------
boot_corrs = zeros(num_boot, num_lags);
parfor b = 1:num_boot
    idx = randsample(N_events, N_events, true);  % bootstrap resample
    boot_corrs(b,:) = mean(corr_mat(idx, :), 1);
end
corr_CI = prctile(boot_corrs, [2.5 97.5], 1);

% ------------------------------
% Shuffle control: circularly shift V1 bias within each trial
% ------------------------------
shuffle_corrs = zeros(num_shuffles, num_lags);
parfor s = 1:num_shuffles
    V1_shifted = zeros(size(V1_z));
    for r = 1:N_events
        shift = randi(T);  % circular shift amount
        V1_shifted(r, :) = circshift(V1_z(r,:), [0, shift]);
    end

    temp_corr = zeros(1, num_lags);
    for li = 1:num_lags
        lag = lags(li);
        if lag < 0
            h = HPC_z(:, 1:end+lag);
            v = V1_shifted(:, 1-lag:end);
        elseif lag > 0
            h = HPC_z(:, 1+lag:end);
            v = V1_shifted(:, 1:end-lag);
        else
            h = HPC_z;
            v = V1_shifted;
        end
        temp_corr(li) = mean(sum(h .* v, 2) ./ (size(h,2) - 1));
    end
    shuffle_corrs(s, :) = temp_corr;
end

shuffle_mean = mean(shuffle_corrs, 1);
shuffle_CI = prctile(shuffle_corrs, [2.5 97.5], 1);
end
