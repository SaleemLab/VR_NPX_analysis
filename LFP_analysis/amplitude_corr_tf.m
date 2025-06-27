function amp_corr = amplitude_corr_tf(amp1, amp2, method)
    % amp1, amp2: [freq x time x trials]
    % method: 'spearman' or 'pearson'
    
    if nargin < 3
        method = 'spearman';  % default
    end

    [nFreq, nTime, nTrials] = size(amp1);
    
    % Reshape to 2D: [freq*time x trials]
    A1 = reshape(amp1, [], nTrials);  % [nFreq*nTime x trials]
    A2 = reshape(amp2, [], nTrials);  % same size

    % Remove rows with too many NaNs
    valid_mask = sum(isnan(A1) | isnan(A2), 2) < (nTrials - 5);  % at least 5 trials

    % Initialize output
    amp_corr = nan(nFreq * nTime, 1);

    % Apply correlation
    if strcmpi(method, 'spearman')
        % Rank transform
        A1 = tiedrank(A1')';  % transpose to [trials x N], rank, transpose back
        A2 = tiedrank(A2')';
    end

    % Normalize (zero-mean, unit std)
    A1z = (A1 - mean(A1, 2, 'omitnan')) ./ std(A1, 0, 2, 'omitnan');
    A2z = (A2 - mean(A2, 2, 'omitnan')) ./ std(A2, 0, 2, 'omitnan');

    % Correlation = dot product of z-scored vectors / (n-1)
    dotprod = sum(A1z .* A2z, 2, 'omitnan');
    amp_corr(valid_mask) = dotprod(valid_mask) ./ (nTrials - 1);

    % Reshape back to [freq x time]
    amp_corr = reshape(amp_corr, nFreq, nTime);
end
