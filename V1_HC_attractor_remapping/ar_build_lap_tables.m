function L = ar_build_lap_tables()
% Per-session lap, block and transition table of RUN1, with licks and rewards
% moved onto the ephys clock.
% Inputs (copies in inputs\): RUN1_decoding_input.mat (lap times and tracks in
% ephys time) and behaviour\per_lap_table.mat (Bonsai lap table, licks).
% Bonsai lap start times (ms) are mapped to ephys lap onsets by a linear fit per
% session (all 22 sessions: same tracks lap by lap, residual <= 0.2 s).
% Output: results\data\lap_tables.mat, struct array L (one per session):
%   subject, session, session_index, animal_index, clock_fit, licks (table: t,
%   side 1 = left / 2 = right, lap), laps (table, one row per lap):
%     lap, track, t_on, t_off           track onset / offset (= grey onset), s
%     block, lap_in_block, block_len, laps_to_end (0 = last lap of block)
%     is_switch      first lap of a block other than the first block
%     prev_track     track of lap k-1 (NaN for lap 1)
%     run_len_prev   laps of prev_track in a row up to lap k-1
%     grey_before    t_on(k) - t_off(k-1), s
%     valid_transition  lap >= 2 and grey_before within P.grey_limits
%     is_probe, is_manual, reward_time (ephys s), speed_lap (mean cm/s)

P = ar_setup();
S = load(fullfile(P.inputs, 'RUN1_decoding_input.mat'), 'D'); D = S.D;
B = load(fullfile(P.inputs, 'behaviour', 'per_lap_table.mat'), 'laps', 'licks', 'runs');

L = struct();
for s = 1:numel(D)
    Ds = D(s);
    rid = B.runs.run_id(B.runs.session_id == s & B.runs.phase == "familiar" & B.runs.run_label == "RUN1");
    assert(numel(rid) == 1, 'Session %d: %d RUN1 runs', s, numel(rid));
    bl = B.laps(B.laps.run_id == rid, :);
    nLaps = size(Ds.lap_times, 1);
    assert(height(bl) >= nLaps && isequal(bl.track(1:nLaps), Ds.lap_track(:)), 'Session %d: lap tracks differ', s);
    bl = bl(1:nLaps, :);
    fit = polyfit(bl.t_start, Ds.lap_times(:, 1), 1);
    resid = Ds.lap_times(:, 1) - polyval(fit, bl.t_start);
    assert(max(abs(resid)) < 0.3, 'Session %d: clock fit residual %.2f s', s, max(abs(resid)));

    T = table();
    T.lap = (1:nLaps)';
    T.track = Ds.lap_track(:);
    T.t_on = Ds.lap_times(:, 1);
    T.t_off = Ds.lap_times(:, 2);
    T.block = cumsum([1; diff(T.track) ~= 0]);
    T.lap_in_block = zeros(nLaps, 1); T.block_len = zeros(nLaps, 1);
    for b = 1:max(T.block)
        in = find(T.block == b);
        T.lap_in_block(in) = 1:numel(in);
        T.block_len(in) = numel(in);
    end
    T.laps_to_end = T.block_len - T.lap_in_block;
    T.is_switch = T.lap_in_block == 1 & T.block > 1;
    T.prev_track = [NaN; T.track(1:end-1)];
    T.run_len_prev = [NaN; T.lap_in_block(1:end-1)];
    T.grey_before = [NaN; T.t_on(2:end) - T.t_off(1:end-1)];
    T.valid_transition = T.lap >= 2 & T.grey_before >= P.grey_limits(1) & T.grey_before <= P.grey_limits(2);
    T.is_probe = bl.is_probe > 0;
    T.is_manual = logical(bl.is_manual);
    T.reward_time = polyval(fit, bl.reward_time);
    in_lap = @(k) Ds.tvec >= T.t_on(k) & Ds.tvec < T.t_off(k);
    T.speed_lap = arrayfun(@(k) mean(Ds.speed(in_lap(k)), 'omitnan'), (1:nLaps)');
    assert(isequal(T.block, bl.block_id) && isequal(T.lap_in_block, bl.lap_in_block), 'Session %d: block ids differ', s);

    lk = B.licks(B.licks.run_id == rid, :);
    K = table();
    K.t = polyval(fit, lk.time);
    K.side = lk.side;
    K.lap = lk.lap;

    L(s).subject = Ds.subject;
    L(s).session = Ds.session;
    L(s).session_index = s;
    L(s).animal_index = find(strcmp(P.subjects, Ds.subject));
    L(s).clock_fit = fit;
    L(s).clock_resid_max = max(abs(resid));
    L(s).laps = T;
    L(s).licks = K;
    fprintf('%2d %s %s: %d laps, %d blocks, %d switches, %d/%d valid transitions (grey %.1f-%.1f s), %d licks\n', ...
        s, Ds.subject, Ds.session, nLaps, max(T.block), sum(T.is_switch), sum(T.valid_transition), nLaps - 1, ...
        min(T.grey_before(T.valid_transition)), max(T.grey_before(T.valid_transition)), height(K));
end

out_dir = fullfile(P.results, 'data');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
save(fullfile(out_dir, 'lap_tables.mat'), 'L');
end
