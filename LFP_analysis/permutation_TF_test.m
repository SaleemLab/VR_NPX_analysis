function [tf_mean, pvals] = permutation_TF_test(data, nperm, alpha)
    % data: [freq x time x event]
    data = permute(data, [3 1 2]);  % [event x freq x time]
    [nEvent, nFreq, nTime] = size(data);
    tf_mean = squeeze(mean(data, 1));  % True mean [freq x time]

    null_dist = zeros(nperm, nFreq, nTime, 'single');

    for p = 1:nperm
        s = RandStream('philox4x32_10', 'Seed', p);
        flip_signs = (rand(s,nEvent,1) > 0.5)*2 - 1;  % ±1
        flip_mat = reshape(flip_signs, [nEvent,1,1]);
        perm_data = squeeze(mean(data .* flip_mat, 1,'omitnan'));
        null_dist(p,:,:) = perm_data;
    end

    tf_mean_expanded = reshape(tf_mean, [1, size(tf_mean,1), size(tf_mean,2)]);

    pvals = squeeze(mean(abs(null_dist) > abs(tf_mean_expanded), 1,'omitnan'));
    % sig_mask = double(pvals < alpha);
    % sig_mask(sig_mask==0) = nan;
end