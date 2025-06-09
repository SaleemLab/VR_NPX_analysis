function [mean_boot, LCI, UCI] = calculate_bootstrap_CI(data_matrix, index)
    nBoot = 1000;
    nSamples = size(data_matrix, 1);
    temp = zeros(nBoot, size(data_matrix,2));
    parfor iBoot = 1:nBoot
        s = RandStream('mrg32k3a','Seed',iBoot);
        resample_idx = datasample(s, index, nSamples);
        temp(iBoot,:) = mean(data_matrix(resample_idx,:), 'omitnan');
    end
    LCI = prctile(temp, 2.5);
    UCI = prctile(temp, 97.5);
    mean_boot = mean(temp,1);
end