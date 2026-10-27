function laps = ar_balanced_train_laps(T, b, h, seed)
% Training laps for the decoder of block b: laps outside blocks b-h..b+h, with
% the majority track randomly subsampled to the number of laps of the minority
% track; empty if a track has no laps
laps = find(~ismember(T.block, b - h:b + h));
l1 = laps(T.track(laps) == 1); l2 = laps(T.track(laps) == 2);
n = min(numel(l1), numel(l2));
if n == 0, laps = []; return; end
rs = RandStream('mrg32k3a', 'Seed', seed);
if numel(l1) > n, l1 = sort(l1(randperm(rs, numel(l1), n))); end
if numel(l2) > n, l2 = sort(l2(randperm(rs, numel(l2), n))); end
laps = sort([l1; l2]);
end
