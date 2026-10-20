function S = Lseq_score(t, cells, Ws, nShuffle, rs)
% Log sequence likelihood ratio L_seq of one candidate event (Huh et al. 2026,
% Nat Commun 17:8304), k = 2 segmentation:
%   log P2*(R|S) = 1/2 [ sum_k log P(C_k -> C_k+1) + sum_{k<l} log P(C_k -> C_l) ]
%   L_seq = log P2*(forward) - log P2*(reverse)
% = 1/2 sum_ab (Ccons_ab + Mall_ab) W_ab, with W = log P - log P' (Lseq_template),
% Ccons = consecutive-spike transition counts and Mall = all ordered spike-pair counts
% of the event (same-cell pairs ignored; they cancel between forward and reverse).
% Pairs of spikes at exactly the same time are excluded from both terms (as the
% authors' code).
% L_seq is divided by the number of contributing spike-pair terms, so it is the mean
% log likelihood ratio per pair (positive = forward, negative = reverse).
%   t, cells : spike times and cell indices of the event's spikes
%   Ws       : cell array of W matrices (templates) scored with the same shuffles
%   nShuffle : cell identity shuffles (permutation of the event's participating cells)
%   rs       : RandStream
% S.L, S.z, S.p_fwd, S.p_rev : 1 x numel(Ws); S.n_spikes, S.n_cells, S.n_terms

[t, order] = sort(t(:));
cells = cells(order);
nW = numel(Ws);
S.n_spikes = numel(cells);
[u, ~, ci] = unique(cells);
m = numel(u);
S.n_cells = m;
S.L = NaN(1, nW); S.z = NaN(1, nW); S.p_fwd = NaN(1, nW); S.p_rev = NaN(1, nW);
S.n_terms = 0;
if m < 2, return; end
[~, ~, g] = unique(t);                          % time groups (ties share a group)
G = full(sparse(g, ci, 1, max(g), m));
before = cumsum(G, 1) - G;                      % spikes of each cell strictly earlier
A = before' * G;
step = diff(t) > 0;                             % consecutive spikes not tied
A = A + accumarray([ci([step; false]) ci([false; step])], 1, [m m]);
A(1:m + 1:end) = 0;
A = A / 2;
S.n_terms = sum(A(:));
if S.n_terms == 0, return; end
perm_idx = zeros(nShuffle, m);
for k = 1:nShuffle
    perm_idx(k, :) = randperm(rs, m);
end
for w = 1:nW
    Wu = Ws{w}(u, u);
    S.L(w) = sum(A .* Wu, 'all') / S.n_terms;
    Ls = zeros(nShuffle, 1);
    for k = 1:nShuffle
        pk = perm_idx(k, :);
        Ls(k) = sum(A .* Wu(pk, pk), 'all');
    end
    Ls = Ls / S.n_terms;
    S.z(w) = (S.L(w) - mean(Ls)) / std(Ls);
    S.p_fwd(w) = (1 + sum(Ls >= S.L(w))) / (1 + nShuffle);
    S.p_rev(w) = (1 + sum(Ls <= S.L(w))) / (1 + nShuffle);
end
end
