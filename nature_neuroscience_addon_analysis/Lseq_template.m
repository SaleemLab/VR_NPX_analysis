function T = Lseq_template(M)
% Ordered firing probability matrix for L_seq (Huh et al. 2026, Nat Commun 17:8304)
% from ordered spike-pair counts accumulated over active-phase trials.
%   M : nCells x nCells, M(j,i) = N_{j->i}, spike pairs in which target i fires after
%       reference j (sum of Lseq_pair_counts over trials)
% For each target cell i:
%   <q_ij> = (N_{j->i} - N_{i->j}) / N_ij, N_ij = N_{j->i} + N_{i->j}
%   prior (pooled over reference cells j, weighted by N_ij):
%     Q_i0     = sum_j N_ij <q_ij> / sum_j N_ij
%     sigma_i2 = per-pair-outcome variance of the +1/-1 ordering outcomes around Q_i0
%     sigma_i02= N_ij-weighted variance of <q_ij> around Q_i0 (between-pair spread)
%   posterior mean Qbar_ij = (N_ij/sigma_i02 <q_ij> + Q_i0/sigma_i2) / (N_ij/sigma_i02 + 1/sigma_i2)
%   (variances as in the authors' figshare code, Make_pMat_f1.m)
%   GLM: logit P(C_j->C_i) = c_i Qbar_ij + c_i0, fit as a binomial GLM of N_{j->i}
%        out of N_ij on Qbar_ij (the authors use a probit link)
% Pairs without any spike pair (N_ij = 0) get P = 0.5 (W = 0), as the authors' code.
% Pairwise probabilities are normalised so P(j->i) + P(i->j) = 1 and bounded to
% [0.01 0.99]. With the large pair counts of RUN laps the shrinkage is negligible and
% Qbar ~ <q>.
% Output T: P (P(a,b) = probability that b fires after a), W = log P - log P'
% (antisymmetric log odds used by Lseq_score), N, q, Qbar, c (nCells x 2: c_i0, c_i).

nC = size(M, 1);
M(1:nC + 1:end) = 0;
N = M + M';
q = (M - M') ./ N;                     % q(j,i) = <q_ij>: target i, reference j
q(N == 0) = NaN;
Qbar = NaN(nC);
Praw = NaN(nC);
c = NaN(nC, 2);
for i = 1:nC
    j = setdiff(find(N(:, i) > 0), i);
    if isempty(j), continue; end
    n = N(j, i); qi = q(j, i);
    Q0 = sum(n .* qi) / sum(n);
    s2 = (sum(M(j, i)) * (1 - Q0)^2 + sum(M(i, j)) * (1 + Q0)^2) / sum(n);
    s02 = max(sum(n .* (qi - Q0).^2) / sum(n), 1e-6);
    s2 = max(s2, 1e-6);
    Qbar(j, i) = (n / s02 .* qi + Q0 / s2) ./ (n / s02 + 1 / s2);
    if numel(j) >= 3 && std(Qbar(j, i)) > 0
        b = logistic_fit(Qbar(j, i), M(j, i), n);
        c(i, :) = b';
        Praw(j, i) = 1 ./ (1 + exp(-(b(1) + b(2) * Qbar(j, i))));
    else
        Praw(j, i) = (1 + Qbar(j, i)) / 2;
    end
end
Praw(isnan(Praw)) = 0.5;
P = Praw ./ (Praw + Praw');
P = min(max(P, 0.01), 0.99);
P(1:nC + 1:end) = 0.5;
T.P = P;
T.W = log(P) - log(P');
T.N = N;
T.q = q;
T.Qbar = Qbar;
T.c = c;
end

function b = logistic_fit(x, k, n)
% Binomial GLM (logit link) of k successes out of n on x by Newton-Raphson (IRLS);
% same estimates as glmfit(x, [k n], 'binomial'), much faster for two parameters.
X = [ones(numel(x), 1) x(:)];
k = k(:); n = n(:);
b = zeros(2, 1);
softplus = @(e) max(e, 0) + log1p(exp(-abs(e)));
dev = @(b) -sum(k .* (X * b) - n .* softplus(X * b));
for it = 1:100
    mu = 1 ./ (1 + exp(-X * b));
    w = n .* mu .* (1 - mu);
    step = (X' * (X .* w) + 1e-10 * eye(2)) \ (X' * (k - n .* mu));
    d0 = dev(b);
    while dev(b + step) > d0 + 1e-12 && norm(step) > 1e-12   % step halving
        step = step / 2;
    end
    b = b + step;
    if max(abs(step)) < 1e-9, break; end
end
end
