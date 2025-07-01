function [tf_mean, pvals, sig_mask,null_dist] = permutation_TF_test(data, nperm, varargin)
p = inputParser;
addParameter(p, 'shuffle_option',1);
addParameter(p, 'alpha',0.05);
addParameter(p, 'data2',[]);
parse(p, varargin{:});
shuffle_option = p.Results.shuffle_option;
alpha = p.Results.alpha;

data = permute(data, [3 1 2]);  % [event x freq x time]
[nEvent, nFreq, nTime] = size(data);

null_dist = zeros(nperm, nFreq, nTime, 'single');

if shuffle_option == 1
    tf_mean = squeeze(mean(data, 1));  % True mean [freq x time]

    for p = 1:nperm
        s = RandStream('philox4x32_10', 'Seed', p);
        % Sign-flip
        flip_signs = (rand(s,nEvent,1) > 0.5)*2 - 1;
        flip_mat = reshape(flip_signs, [nEvent,1,1]);
        null_dist(p,:,:) = mean(data .* flip_mat, 1, 'omitnan');
    end
elseif shuffle_option == 2 & ~isempty(data2)
    complex_diff = exp(1i * (data - data2));  % element-wise phase diff
    tf_mean = squeeze(abs(mean(complex_diff, 1, 'omitnan')));  % [freq x time]

    % Circular shift on data2 (second region), keep data fixed
    for p = 1:nperm
        shifted_data2 = nan(size(data2), 'like', data2);

        s = RandStream('philox4x32_10', 'Seed', p);  % consistent seed

        for e = 1:nEvent
            shift_amt = randi(s, [1, nTime]);  % circular shift amount
            shifted_data2(e,:,:) = circshift(data2(e,:,:), [0 0 shift_amt]);
        end

        % Compute PLV between original data and shifted data2
        complex_diff = exp(1i * (data - shifted_data2));  % element-wise phase diff
        perm_plv = squeeze(abs(mean(complex_diff, 1, 'omitnan')));  % [freq x time]
        null_dist(p,:,:) = perm_plv;
    end
end


tf_mean_expanded = reshape(tf_mean, [1, size(tf_mean,1), size(tf_mean,2)]);
abs_data = abs(tf_mean_expanded);
abs_null = abs(null_dist);
threshold = prctile(abs_null, 100*(1-alpha), 1);

raw_sig = abs_data > threshold;
real_cluster_stat = sum(raw_sig(:));
null_cluster_stats = squeeze(sum(sum(abs_null > threshold, 2), 3));
cluster_p = mean(null_cluster_stats >= real_cluster_stat);

if cluster_p < alpha
    sig_mask = squeeze(raw_sig);
else
    sig_mask = zeros(nFreq, nTime);
end

pvals = squeeze(mean(abs_null > abs_data, 1, 'omitnan'));
end
