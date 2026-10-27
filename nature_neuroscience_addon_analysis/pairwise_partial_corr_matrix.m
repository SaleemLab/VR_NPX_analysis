function [r, c] = pairwise_partial_corr_matrix(X, Y, Z1, Z2)
% Partial correlation between every column of X (events x nV) and every column of Y
% (events x nH), controlling for two covariates that may differ per cell:
%   Z1 (events x nV): covariate of the X cell (e.g. leave-one-out V1 MUA)
%   Z2 (events x nH): covariate of the Y cell (e.g. leave-one-out HC MUA)
% For pair (i, j) this equals corr of the residuals of X(:, i) and Y(:, j) after
% least-squares regression on [1 Z1(:, i) Z2(:, j)] (as partial_corr in
% pairwise_V1_HC_coupling_POST), computed in closed form for all pairs at once.
% Pairs with a singular covariance return NaN.
% Second output c: the partial covariance (covariance of the two residuals, normalised by
% the number of events), in spikes^2 per window.
n = size(X, 1);
X = X - mean(X, 1); Y = Y - mean(Y, 1); Z1 = Z1 - mean(Z1, 1); Z2 = Z2 - mean(Z2, 1);
Sxx = sum(X .^ 2, 1)' / n;          % nV x 1
Syy = sum(Y .^ 2, 1) / n;           % 1 x nH
Sxy = X' * Y / n;                   % nV x nH
p = sum(X .* Z1, 1)' / n;           % cov(x, z1), nV x 1
q = X' * Z2 / n;                    % cov(x, z2), nV x nH
s = Z1' * Y / n;                    % cov(y, z1), nV x nH
t = sum(Y .* Z2, 1) / n;            % cov(y, z2), 1 x nH
a = sum(Z1 .^ 2, 1)' / n;           % var(z1), nV x 1
d = sum(Z2 .^ 2, 1) / n;            % var(z2), 1 x nH
b = Z1' * Z2 / n;                   % cov(z1, z2), nV x nH
det_z = a .* d - b .^ 2;
Sxx_z = Sxx - (p .^ 2 .* d - 2 * p .* q .* b + q .^ 2 .* a) ./ det_z;
Syy_z = Syy - (s .^ 2 .* d - 2 * s .* t .* b + t .^ 2 .* a) ./ det_z;
Sxy_z = Sxy - (p .* s .* d - p .* t .* b - q .* s .* b + q .* t .* a) ./ det_z;
r = Sxy_z ./ sqrt(Sxx_z .* Syy_z);
bad = ~isfinite(r) | det_z <= eps * (a .* d) | Sxx_z <= 0 | Syy_z <= 0;
r(bad) = NaN;
c = Sxy_z;
c(bad) = NaN;
end
