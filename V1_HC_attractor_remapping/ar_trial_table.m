function T = ar_trial_table(x)
% One row per transition from epoch-table rows x (one region and decoder):
% zLO per epoch (end_prev, grey_early, grey_late, onset, onset_100, onset_200,
% onset_late) and, for grey_late and the onset epochs, log spike count (spk_*),
% speed (spd_*) and licks L - R (lick_*), z-scored within session (spikes,
% speed, lap_in_block). Also the late grey of the previous / next transition in
% the same block (grey_late_prevtr / _nexttr, for temporal controls), the track
% signs, track pair and grouping factors for LMEs.
epochs = ["end_prev", "grey_early", "grey_late", "onset", "onset_100", "onset_200", "onset_late"];
with_cov = ["grey_late", "onset", "onset_100", "onset_200", "onset_late"];
keys = unique(x(:, {'session', 'lap'}), 'rows');
T = keys;
d = x(x.epoch == "grey_late", :);
[~, ia, ib] = intersect(keys, d(:, {'session', 'lap'}));
for v = {'animal', 'track', 'prev_track', 'is_switch', 'lap_in_block', 'block', 'grey_dur'}
    T.(v{1}) = NaN(height(T), 1); T.(v{1})(ia) = d.(v{1})(ib);
end
for e = epochs
    d = x(x.epoch == e, :);
    [~, ia, ib] = intersect(keys, d(:, {'session', 'lap'}));
    T.(e) = NaN(height(T), 1); T.(e)(ia) = d.zLO(ib);
    if ismember(e, with_cov)
        s = "_" + e;
        T.("spk" + s) = NaN(height(T), 1); T.("spk" + s)(ia) = log1p(d.spikes_per_win(ib));
        T.("spd" + s) = NaN(height(T), 1); T.("spd" + s)(ia) = d.speed(ib);
        T.("lick" + s) = NaN(height(T), 1); T.("lick" + s)(ia) = d.licks_L(ib) - d.licks_R(ib);
    end
end
T.is_switch = logical(T.is_switch);
T.sgn_next = 2 * (T.track == 1) - 1;
T.sgn_prev = 2 * (T.prev_track == 1) - 1;
T.grey_late_prevtr = NaN(height(T), 1); T.grey_late_nexttr = NaN(height(T), 1);
for i = 1:height(T)
    j = find(T.session == T.session(i) & T.lap == T.lap(i) - 1 & T.block == T.block(i) & ~T.is_switch);
    if ~isempty(j) && ~T.is_switch(i), T.grey_late_prevtr(i) = T.grey_late(j); end
    j = find(T.session == T.session(i) & T.lap == T.lap(i) + 1 & T.block == T.block(i));
    if ~isempty(j), T.grey_late_nexttr(i) = T.grey_late(j); end
end
zs = ["spk_" + with_cov, "spd_" + with_cov, "lap_in_block"];
[g, ~] = findgroups(T.session);
for v = zs
    mu = splitapply(@(z) mean(z, 'omitnan'), T.(v), g);
    sd = splitapply(@(z) std(z, 'omitnan'), T.(v), g);
    T.(v) = (T.(v) - mu(g)) ./ sd(g);
end
T.gl_x_switch = T.grey_late .* T.is_switch;
T.pair = categorical(T.prev_track * 10 + T.track);
T.animal_c = categorical(T.animal);
T.sess_c = categorical(T.session);
T.blk_c = categorical(T.session * 100 + T.block);
end
