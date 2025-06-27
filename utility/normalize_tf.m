function S_out = normalize_tf(S, baseline_idx, method)
% [freq x time x trials]
    switch method
        case 'dB'
            baseline_mean = mean(S(:, baseline_idx, :), 2, 'omitnan');  % [freq x 1 x trials]
            S_out = 10 * log10(bsxfun(@rdivide, S, baseline_mean));
        case 'z'
            baseline_mean = mean(S(:, baseline_idx, :), 2, 'omitnan');
            baseline_std  = std(S(:, baseline_idx, :), 0, 2, 'omitnan');
            S_out = bsxfun(@rdivide, ...
                    bsxfun(@minus, S, baseline_mean), ...
                    baseline_std);

        case 'z_all'
            % Reshape to [freq x (time * event)]
            S_reshaped = reshape(S, size(S,1), []);

            % Compute mean and std per frequency across all bins
            mu = mean(S_reshaped, 2, 'omitnan');        % [freq x 1]
            sigma = std(S_reshaped, 0, 2, 'omitnan');   % [freq x 1]

            % Z-score normalization per frequency
            S_out = bsxfun(@rdivide, bsxfun(@minus, S, mu), sigma);  % [freq x time x event]


        otherwise
            error('Unknown normalization method: %s', method);
    end
end
