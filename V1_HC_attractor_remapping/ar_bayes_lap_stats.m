function [occ, cnt] = ar_bayes_lap_stats(Ds, T, P)
% Occupancy (s) and spike counts per lap and 10 cm bin, running samples only
% (in lap, position defined, speed > P.bayes_speed), as in RUN1_session_behaviour.
% occ: nLaps x nPos;  cnt: nCells x nPos x nLaps (cells of Ds.cell_id)
tvec = Ds.tvec;
nLaps = height(T); nPos = numel(P.pos_edges) - 1; nCells = numel(Ds.cell_id);
lap = zeros(size(tvec));
for k = 1:nLaps
    lap(tvec >= T.t_on(k) & tvec < T.t_off(k)) = k;
end
pos_bin = discretize(Ds.position, P.pos_edges);
running = lap > 0 & ~isnan(pos_bin) & Ds.speed > P.bayes_speed;
dt = [diff(tvec); median(diff(tvec))];
occ = accumarray([lap(running) pos_bin(running)], dt(running), [nLaps nPos]);
[~, cell_idx] = ismember(Ds.spike_id, Ds.cell_id);
keep = cell_idx > 0 & Ds.spike_times >= tvec(1) & Ds.spike_times <= tvec(end);
st = Ds.spike_times(keep); sc = cell_idx(keep);
sample = interp1(tvec, (1:numel(tvec))', st, 'nearest');
sl = lap(sample);
ok = running(sample) & sl > 0;
ok(ok) = st(ok) >= T.t_on(sl(ok)) & st(ok) < T.t_off(sl(ok));
cnt = accumarray([sc(ok) pos_bin(sample(ok)) sl(ok)], 1, [nCells nPos nLaps]);
end
