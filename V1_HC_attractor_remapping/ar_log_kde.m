function lp = ar_log_kde(train, test, h)
% Log density of a Gaussian product-kernel KDE (bandwidth h in every
% dimension), the same estimate as mvksdensity(train, test, 'Bandwidth', h) but
% computed in log space, so test points far from the training cloud give finite
% log densities instead of 0 (and log odds of +-Inf or NaN).
n = size(train, 1); d = size(train, 2);
lp = NaN(size(test, 1), 1);
const = -log(n) - d * log(h) - d / 2 * log(2 * pi);
chunk = 2000;
for i = 1:chunk:size(test, 1)
    idx = i:min(i + chunk - 1, size(test, 1));
    A = -pdist2(test(idx, :), train, 'squaredeuclidean') / (2 * h^2);
    m = max(A, [], 2);
    lp(idx) = m + log(sum(exp(A - m), 2)) + const;
end
end
