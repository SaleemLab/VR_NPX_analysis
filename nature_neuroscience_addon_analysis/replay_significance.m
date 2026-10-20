function [sig, p_comb] = replay_significance(E, types, alpha)
% Replay significance of events per track from one or more shuffle types
% (Takigawa et al. 2024 eLife: with several shuffles, the event must pass every one of them;
% the combined p is the largest of their p-values).
% E:     struct of p-values, E.(type) = nEv x 2 (T1, T2); NaN = not scored / not run
% types: cell array of shuffle types to combine, e.g. {'spike', 'field'}
% alpha: scalar or vector of alpha levels
% sig:    nEv x 2 x numel(alpha) logical, p < alpha on that track for every type in types
%         (NaN p -> not significant). An event significant on both tracks is a multi-track
%         event: it counts as a T1 and a T2 event, and once in the proportion of events.
% p_comb: nEv x 2, largest p over types (NaN if any type is NaN)
p_all = zeros([size(E.(types{1})) numel(types)]);
for k = 1:numel(types)
    p_all(:, :, k) = E.(types{k});
end
p_comb = max(p_all, [], 3, 'includenan');
sig = p_comb < reshape(alpha, 1, 1, []);
end
