function [p2, p1] = paired_perm_test(d, tail)
% Exact sign-flip permutation test on mean(d), enumerating all 2^n sign patterns.
% d: paired differences (finite). tail: 'right' (mean(d) > 0) or 'left'.
% Patterns are split into low bits (enumerated as a vector) and high bits (looped)
% so n up to ~25 stays cheap in memory.
a = abs(d(:)); n = numel(a);
T_obs = sum(d);                          % sum is equivalent to mean for fixed n
tol = 1e-12 * max(1, sum(a));
k = min(n, 16);
low = (1 - 2 * (dec2bin(0:2^k - 1, k) == '1')) * a(1:k);
high_bits = n - k;
n_abs = 0; n_ge = 0; n_le = 0;
for h = 0:2^high_bits - 1
    if high_bits > 0
        s = 1 - 2 * (dec2bin(h, high_bits) == '1');
        T = low + s * a(k + 1:end);
    else
        T = low;
    end
    n_abs = n_abs + sum(abs(T) >= abs(T_obs) - tol);
    n_ge = n_ge + sum(T >= T_obs - tol);
    n_le = n_le + sum(T <= T_obs + tol);
end
N = 2^n;
p2 = n_abs / N;
if strcmp(tail, 'right')
    p1 = n_ge / N;
else
    p1 = n_le / N;
end
end
