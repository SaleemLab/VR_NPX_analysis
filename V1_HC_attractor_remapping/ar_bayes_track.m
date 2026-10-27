function [LO, post] = ar_bayes_track(counts, maps, tau, method, rate_floor)
% Track log odds with position marginalised (flat prior over the 2 x nPos
% states); same likelihoods as nature_neuroscience_addon_analysis/bayes_decode_two_track.m
%   counts: cells x bins;  maps: cells x states (Hz);  tau: bin duration (s,
%   scalar or 1 x bins; Poisson only);  method: 'multinomial' or 'poisson'
% LO: bins x 1, log P(T1) - log P(T2);  post: states x bins (single)
maps = max(maps, rate_floor);
switch method
    case 'poisson'
        tau = tau(:)';
        if isscalar(tau), tau = tau * ones(1, size(counts, 2)); end
        logL = log(maps)' * counts - sum(maps, 1)' * tau;
    case 'multinomial'
        logL = log(maps ./ sum(maps, 1))' * counts;
end
nPos = size(maps, 2) / 2;
LO = logsumexp(logL(1:nPos, :)) - logsumexp(logL(nPos + 1:end, :));
LO = LO(:);
if nargout > 1
    logL = logL - max(logL, [], 1);
    post = exp(logL);
    post = single(post ./ sum(post, 1));
end
end

function y = logsumexp(x)
m = max(x, [], 1);
y = m + log(sum(exp(x - m), 1));
end
