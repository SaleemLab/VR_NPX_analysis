function G = ar_bin_session(Ds, T, P)
% Spike counts and behaviour of one RUN1 session in 100 ms non-overlapping bins
% (P.bin_width), of two kinds, stacked in one list:
%   kind 0: a grid over the whole session (laps and grey screens), starting at
%           ceil(tvec(1) / 0.05) * 0.05 (the same bins as the PLS training bins
%           of the earlier 50 ms-stepped version)
%   kind 1: bins locked to each valid transition's lap end t_off(k-1), with an
%           edge exactly at the event: P.align_window around it
%   kind 2: the same, locked to track onset t_on(k)
% Event-locked bins never straddle a lap end or an onset.
%   Ds : one element of RUN1_decoding_input.mat;  T : its lap table (lap_tables.mat)
% G fields (one row per bin):
%   kind, trans (lap k of the transition for kinds 1/2, else 0), rel (bin index
%   from the event: bin covers event + [rel, rel + 1] * bin_width), t_start,
%   t_centre, counts (nBins x nCells, single), pos (cm at the centre, NaN outside
%   laps), speed (mean cm/s), lap (k if the bin lies fully inside lap k, else 0),
%   segment (lap k whose grey screen or lap contains the centre:
%   (t_off(k-1), t_off(k)]), t_from_on (centre - t_on(segment)),
%   t_from_off_prev (centre - t_off(segment-1), NaN for segment 1),
%   nonoverlap (true for grid bins: used for training and normalisation)

w = P.bin_width;
t0 = ceil(Ds.tvec(1) / 0.05) * 0.05;
t1 = Ds.tvec(end) - w;
grid = (t0:w:t1)';
trans = find(T.valid_transition);
rel = (round(P.align_window(1) / w):round(P.align_window(2) / w) - 1)';
nr = numel(rel);
off_start = T.t_off(trans - 1)' + rel * w;    % nr x nTrans
on_start = T.t_on(trans)' + rel * w;
G.kind = [zeros(numel(grid), 1); ones(numel(off_start), 1); 2 * ones(numel(on_start), 1)];
G.trans = [zeros(numel(grid), 1); repelem(trans, nr); repelem(trans, nr)];
G.rel = [zeros(numel(grid), 1); repmat(rel, numel(trans), 1); repmat(rel, numel(trans), 1)];
G.t_start = [grid; off_start(:); on_start(:)];
G.t_centre = G.t_start + w / 2;
ok = G.t_start >= Ds.tvec(1) & G.t_start + w <= Ds.tvec(end);
G = structfun(@(x) x(ok), G, 'UniformOutput', false);
nBins = numel(G.t_start);

% spike counts: bins of each kind are non-overlapping, so count per kind with
% the sorted bin edges (start, end) of that kind
nCells = numel(Ds.cell_id);
[~, cell_idx] = ismember(Ds.spike_id, Ds.cell_id);
keep = cell_idx > 0;
st = Ds.spike_times(keep); sc = cell_idx(keep);
counts = zeros(nBins, nCells, 'single');
for kd = 0:2
    idx = find(G.kind == kd);
    [s, order] = sort(G.t_start(idx));
    idx = idx(order);
    edges = reshape([s s + w]', [], 1);
    % contiguous bins share an edge: nudge each end edge down by 1 us so that
    % the edge vector is strictly increasing (spikes exactly on an edge are
    % negligible)
    edges(2:2:end) = edges(2:2:end) - 1e-6;
    b = discretize(st, edges);
    real = mod(b, 2) == 1;                      % odd intervals are bins
    bi = (b(real) + 1) / 2;
    counts(idx, :) = single(accumarray([bi sc(real)], 1, [numel(idx) nCells]));
end
G.counts = counts;

% behaviour: mean speed over the bin, position at the centre
sample_speed = Ds.speed; sample_speed(isnan(sample_speed)) = 0;
cs = [0; cumsum(sample_speed)];
i0 = discretize(G.t_start, [-Inf; Ds.tvec; Inf]);       % first sample index >= start (approx)
i1 = discretize(G.t_start + w, [-Inf; Ds.tvec; Inf]);
n = max(i1 - i0, 1);
G.speed = (cs(min(i1, numel(cs))) - cs(min(i0, numel(cs)))) ./ n;
sample = interp1(Ds.tvec, (1:numel(Ds.tvec))', G.t_centre, 'nearest', 'extrap');
G.pos = Ds.position(sample);

nLaps = height(T);
G.lap = zeros(nBins, 1);
for k = 1:nLaps
    G.lap(G.t_start >= T.t_on(k) & G.t_start + w <= T.t_off(k)) = k;
end
G.pos(G.lap == 0) = NaN;
G.segment = discretize(G.t_centre, [-Inf; T.t_off(1:end-1); Inf], 'IncludedEdge', 'right');
G.t_from_on = G.t_centre - T.t_on(G.segment);
prev_off = [NaN; T.t_off(1:end-1)];
G.t_from_off_prev = G.t_centre - prev_off(G.segment);
G.nonoverlap = G.kind == 0;
end
