function [mean_boot, LCI, UCI,temp] = calculate_bootstrap_CI(data_matrix, index,varargin)
% SETTINGS
p = inputParser;
addParameter(p, 'nSamples',  size(data_matrix, 1));
addParameter(p, 'nBoot', 1000);
addParameter(p, 'nseed', 0);
parse(p, varargin{:});

nSamples = p.Results.nSamples;
nBoot = p.Results.nBoot;
nseed = p.Results.nseed;

temp = zeros(nBoot, size(data_matrix,2));

data_matrix1 = data_matrix(isfinite(data_matrix));
data_matrix(data_matrix>=inf) = prctile(data_matrix1,99.5);
data_matrix(data_matrix<=-inf) = prctile(data_matrix1,0.5);

parfor iBoot = 1:nBoot
    s = RandStream('mrg32k3a','Seed',iBoot+nseed);
    resample_idx = datasample(s, index, nSamples);
    temp(iBoot,:) = mean(data_matrix(resample_idx,:), 'omitnan');
end
LCI = prctile(temp, 2.5);
UCI = prctile(temp, 97.5);
mean_boot = mean(temp,1);
end