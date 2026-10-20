function posterior = bayes_decode_two_track(counts, rate_maps, time_bin, method, rate_floor)
% Bayesian position decoding over the joint state space of both tracks.
%   counts    : nCells x nBins spike counts
%   rate_maps : nCells x nStates rate maps in Hz (states = [Track 1 bins, Track 2 bins])
%   time_bin  : bin duration in s (scalar or 1 x nBins)
%   method    : 'poisson'     - P(n|x) = prod_i Poisson(n_i; tau*f_i(x))
%               'multinomial' - P(n|x,N) = Multinomial(n; N, p(x)), p_i(x) = f_i(x)/sum_j f_j(x)
%   rate_floor: minimum rate (Hz) applied to the maps before decoding (default 0.1)
% posterior: nStates x nBins, each column sums to 1 over both tracks (flat prior).
%
% Multinomial is the Poisson likelihood conditioned on the total spike count N:
%   log P_poisson(n|x) = sum_i n_i log p_i(x) + [N log F(x) - tau F(x)] + const,
% with F(x) = sum_i f_i(x). Dropping the bracketed population-rate term makes the
% decoder depend only on which cells fire relative to each other, so a global gain
% change (all rates x c, e.g. sleep vs RUN) leaves the posterior unchanged. Bins with
% no spikes have a flat multinomial posterior.

if nargin < 5 || isempty(rate_floor), rate_floor = 0.1; end
rate_maps = max(rate_maps, rate_floor);
counts = full(double(counts));

switch lower(method)
    case 'poisson'
        tau = time_bin(:)';
        if isscalar(tau), tau = tau * ones(1, size(counts, 2)); end
        logL = log(rate_maps)' * counts - sum(rate_maps, 1)' * tau;   % nStates x nBins
    case 'multinomial'
        p = rate_maps ./ sum(rate_maps, 1);
        logL = log(p)' * counts;
    otherwise
        error('Unknown method %s', method);
end

logL = logL - max(logL, [], 1);
posterior = exp(logL);
posterior = posterior ./ sum(posterior, 1);
end
