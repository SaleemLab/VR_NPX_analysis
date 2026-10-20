function M = Lseq_pair_counts(t, cells, nCells)
% Ordered spike-pair counts of one trial for the sequence-likelihood (L_seq) method
% (Huh et al. 2026, Nat Commun 17:8304).
%   t, cells : spike times and cell indices (1..nCells) of one trial
% M(a,b) = number of spike pairs (any lag within the trial) in which cell b fires
% after cell a, i.e. N_{a->b} with reference a and target b. Same-cell pairs are 0.
% Spikes at exactly the same time are ordered neither way (excluded, as the authors'
% code).
if isempty(t)
    M = zeros(nCells);
    return
end
[~, ~, g] = unique(t(:));                       % time groups (ties share a group)
G = full(sparse(g, cells(:), 1, max(g), nCells));
before = cumsum(G, 1) - G;                      % spikes of each cell strictly earlier
M = before' * G;
M(1:nCells + 1:end) = 0;
end
