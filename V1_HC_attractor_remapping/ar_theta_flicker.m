function C = ar_theta_flicker(varargin)
% Theta-paced flicker between track representations around track onset, in
% the spirit of Jezek et al. 2011 (Nature, doi 10.1038/nature10439): decode every
% theta cycle and ask whether, on the first lap of a block, cycles keep
% expressing the old track for a few cycles before settling on the new one.
% Theta cycles: HPC MUA phase (ar_theta_phase; phase 0 = minimal population
% firing, cycles between consecutive phase-0 points), 100-200 ms long, with the
% animal running (mean speed > 'speed' cm/s) through the cycle. Cycles are
% numbered from track onset: 1, 2, ... start at or after onset; -1, -2, ... end
% at or before onset; the cycle that straddles onset is 0.
% Per cycle and region (HPC; V1 without P.V1_exclude_subjects):
%   LO     multinomial Bayes track log odds (position marginalised), block-held-
%          out lobo3 decoder with track-balanced training laps (as ar_decode_session)
%   pv_r   Jezek-style population-vector correlation of the cycle's spike counts
%          with each track's rate map at the current position (in-lap cycles)
%   other  decoded as the track other than the upcoming lap's (LO sign), only
%          for cycles with >= 'min_spikes' spikes (else undetermined)
% Summaries: fraction of 'other' cycles and LO toward the current track by cycle
% number, first lap of a block vs laps 6+; label changes in cycles 1-6; the
% first cycle after which the current track holds; trial GLME of other-track
% cycles 1-5 on the late-grey bias toward the other track (x switch), plus HPC-V1
% late-grey coupling and V1 bias (ar_grey_coupling); per-session example plots.
% Output: results\theta_flicker\ (cycles.mat, stats csv, figures, examples\)
p = inputParser;
addParameter(p, 'speed', 5);
addParameter(p, 'min_spikes', 3);
addParameter(p, 'window', [-1.5 1.5]);
addParameter(p, 'examples', true);
addParameter(p, 'summary_only', false);   % reuse the saved cycles.mat
parse(p, varargin{:});
O = p.Results;
P = ar_setup();
Lt = load(fullfile(P.results, 'data', 'lap_tables.mat'), 'L'); L = Lt.L;
out_dir = fullfile(P.results, 'theta_flicker');
if ~exist(fullfile(out_dir, 'examples'), 'dir'), mkdir(fullfile(out_dir, 'examples')); end
if O.summary_only
    S = load(fullfile(out_dir, 'cycles.mat'), 'C'); C = S.C;
    summarise(C, O, P, out_dir);
    if O.examples
        for s = unique(C.session)'
            plot_examples(C(C.session == s & C.region == "HPC", :), L(s), P, out_dir);
        end
    end
    return
end
S = load(fullfile(P.inputs, 'RUN1_decoding_input.mat'), 'D'); D = S.D;
bad = find(ismember(P.subjects, P.V1_exclude_subjects));
nPos = numel(P.pos_edges) - 1;
rows = {};
for s = 1:numel(D)
    Ds = D(s); T = L(s).laps;
    F = load(fullfile(P.results, 'theta_phase', sprintf('session_%02d.mat', s)), 't', 'phase');
    ph = double(F.phase); tt = F.t;
    wrap = find(ph(1:end - 1) > 2 & ph(2:end) < -2) + 1;
    cs = tt(wrap(1:end - 1)) - 0.005; ce = tt(wrap(2:end)) - 0.005;
    dur = ce - cs;
    keep = dur >= 0.1 & dur <= 0.2;
    cs = cs(keep); ce = ce(keep);
    cspd = cycle_mean(Ds.tvec, Ds.speed, cs, ce);
    cpos = interp1(Ds.tvec, Ds.position, (cs + ce) / 2, 'nearest');
    trans = find(T.valid_transition);
    % cycles in the window around each onset
    tr_of = zeros(numel(cs), 1); idx = NaN(numel(cs), 1);
    for i = 1:numel(trans)
        k = trans(i); t0 = T.t_on(k);
        in = (cs + ce) / 2 >= t0 + O.window(1) & (cs + ce) / 2 <= t0 + O.window(2);
        post = in & cs >= t0; pre = in & ce <= t0; mid = in & cs < t0 & ce > t0;
        ip = find(post); idx(ip) = 1:numel(ip);
        iq = find(pre); idx(iq) = -(numel(iq):-1:1);
        idx(mid) = 0;
        tr_of(in) = k;
    end
    use = tr_of > 0 & cspd > O.speed;
    cs = cs(use); ce = ce(use); idx = idx(use); tr_of = tr_of(use); cspd = cspd(use); cpos = cpos(use);
    in_lap = cs >= T.t_on(tr_of) & ce <= T.t_off(tr_of);
    cpos(~in_lap) = NaN;
    [occ, cnt] = ar_bayes_lap_stats(Ds, T, P);
    nBlocks = max(T.block);
    for reg = {'HPC', 'V1'}
        if strcmp(reg{1}, 'V1') && ismember(L(s).animal_index, bad), continue; end
        cells = find(contains(Ds.region, reg{1}));
        counts = cycle_counts(Ds, cells, cs, ce);          % cycles x cells
        LO = NaN(numel(cs), 1); r_cur = LO; r_oth = LO;
        pb = discretize(cpos, P.pos_edges);
        for b = 1:nBlocks
            w = find(T.block(tr_of) == b);
            if isempty(w), continue; end
            train = ismember((1:height(T))', ar_balanced_train_laps(T, b, 1, 1000 * s + 10 * b + 1));
            maps = ar_rate_maps(occ, cnt(cells, :, :), train, T.track, nPos);
            LO(w) = ar_bayes_track(counts(w, :)', maps, 1, 'multinomial', P.rate_floor);
            for j = w(:)'
                if isnan(pb(j)), continue; end
                cur = T.track(tr_of(j)); oth = 3 - cur;
                x = counts(j, :)';
                if sum(x > 0) < 2, continue; end
                r_cur(j) = corr(x, maps(:, (cur - 1) * nPos + pb(j)));
                r_oth(j) = corr(x, maps(:, (oth - 1) * nPos + pb(j)));
            end
        end
        nspk = sum(counts, 2);
        sgn = 2 * (T.track(tr_of) == 1) - 1;
        X = table();
        X.session = repmat(s, numel(cs), 1); X.animal = repmat(L(s).animal_index, numel(cs), 1);
        X.region = repmat(string(reg{1}), numel(cs), 1);
        X.lap = tr_of; X.cycle = idx; X.t_start = cs; X.dur = ce - cs; X.speed = cspd; X.pos = cpos;
        X.track = T.track(tr_of); X.prev_track = T.prev_track(tr_of); X.lap_in_block = T.lap_in_block(tr_of);
        X.is_switch = T.is_switch(tr_of); X.block = T.block(tr_of);
        X.n_spikes = nspk; X.LO = LO; X.LO_cur = LO .* sgn;
        X.other = double(X.LO_cur < 0); X.other(nspk < O.min_spikes | X.LO_cur == 0) = NaN;
        X.r_cur = r_cur; X.r_oth = r_oth;
        rows{end + 1} = X; %#ok<AGROW>
    end
    fprintf('Session %2d: %d cycles around %d onsets\n', s, numel(cs), numel(trans));
end
C = vertcat(rows{:});
save(fullfile(out_dir, 'cycles.mat'), 'C', 'O', '-v7.3');
summarise(C, O, P, out_dir);
if O.examples
    for s = unique(C.session)'
        plot_examples(C(C.session == s & C.region == "HPC", :), L(s), P, out_dir);
    end
end
end

%% ------------------------------------------------------------------------
function m = cycle_mean(tv, v, cs, ce)
v(isnan(v)) = 0;
c = [0; cumsum(v)];
i0 = discretize(cs, [-Inf; tv; Inf]); i1 = discretize(ce, [-Inf; tv; Inf]);
n = max(i1 - i0, 1);
m = (c(min(i1, numel(c))) - c(min(i0, numel(c)))) ./ n;
end

function counts = cycle_counts(Ds, cells, cs, ce)
% spike counts per (non-overlapping) cycle and cell
[~, ci] = ismember(Ds.spike_id, Ds.cell_id(cells));
ok = ci > 0; st = Ds.spike_times(ok); ci = ci(ok);
[s0, order] = sort(cs); e0 = ce(order);
edges = reshape([s0 e0 - 1e-6]', [], 1);
b = discretize(st, edges);
real = mod(b, 2) == 1;
bi = (b(real) + 1) / 2;
c = accumarray([bi ci(real)], 1, [numel(cs) numel(cells)]);
counts = zeros(numel(cs), numel(cells));
counts(order, :) = c;
end

%% ------------------------------------------------------------------------
function summarise(C, O, P, out_dir)
acol = lines(6);
cyc = -6:10;
groups = {'first lap (switch)', @(x) x.lap_in_block == 1 & x.block > 1; 'lap 2', @(x) x.lap_in_block == 2; ...
    'laps 6+', @(x) x.lap_in_block >= 6};
gcol = [0.85 0.2 0.2; 0.95 0.6 0.1; 0.1 0.25 0.6];
stats = {};
fig = figure('Position', [50 50 1800 760], 'Color', 'w');
tl = tiledlayout(2, 4, 'TileSpacing', 'compact');
for r = 1:2
    reg = ["HPC", "V1"]; reg = reg(r);
    x0 = C(C.region == reg, :);
    us = unique(x0.session);
    % A: fraction of cycles decoded as the other track; B: LO toward current
    for mI = 1:2
        nexttile(tl); hold on;
        for g = 1:size(groups, 1)
            M = NaN(numel(us), numel(cyc));
            for i = 1:numel(us)
                xs = x0(x0.session == us(i), :);
                xs = xs(groups{g, 2}(xs), :);
                for c = 1:numel(cyc)
                    xc = xs(xs.cycle == cyc(c), :);
                    if mI == 1, M(i, c) = mean(xc.other, 'omitnan'); else, M(i, c) = mean(xc.LO_cur, 'omitnan'); end
                end
            end
            m = mean(M, 1, 'omitnan'); se = std(M, 0, 1, 'omitnan') ./ sqrt(sum(isfinite(M), 1));
            errorbar(cyc + (g - 2) * 0.1, m, se, 'o-', 'Color', gcol(g, :), 'MarkerFaceColor', gcol(g, :), 'CapSize', 0, 'LineWidth', 1.2);
        end
        xline(0.5, 'k-');
        if mI == 1, yline(0.5, 'k:'); ylabel('Fraction of cycles decoded as other track'); else, yline(0, 'k:'); ylabel('Track log odds toward current track'); end
        xlabel('Theta cycle from onset (1 = first full cycle)');
        title(sprintf('%s: %s', reg, ifelse(mI == 1, 'other-track cycles', 'per-cycle log odds')), 'FontWeight', 'normal');
        if r == 1 && mI == 1, legend(groups(:, 1), 'Box', 'off', 'Location', 'northeast'); end
        set(gca, 'XTick', cyc(1:2:end));
    end
    % C: switch - laps 6+ in other-track fraction, cycles 1-3 and 4-6, per session
    nexttile(tl); hold on;
    wins = {1:3, 4:6, 7:10, -3:-1}; wl = {'cycles 1-3', '4-6', '7-10', 'pre -3..-1'};
    for w = 1:numel(wins)
        v = NaN(numel(us), 1); an = v;
        for i = 1:numel(us)
            xs = x0(x0.session == us(i) & ismember(x0.cycle, wins{w}), :);
            a = mean(xs.other(groups{1, 2}(xs)), 'omitnan'); b = mean(xs.other(groups{3, 2}(xs)), 'omitnan');
            v(i) = a - b; an(i) = xs.animal(1);
        end
        ok = isfinite(v);
        scatter(w + 0.1 * randn(sum(ok), 1), v(ok), 16, acol(an(ok), :), 'filled', 'MarkerFaceAlpha', 0.6);
        errorbar(w + 0.3, mean(v, 'omitnan'), std(v, 'omitnan') / sqrt(sum(ok)), 'ko', 'MarkerFaceColor', 'k', 'CapSize', 0, 'LineWidth', 1.5);
        [ua, ~, ia] = unique(an(ok)); am = accumarray(ia, v(ok), [], @mean);
        pv = signrank(v(ok));
        text(w, max(v(ok)) + 0.02, sprintf('p %.2g\n%d/%d', pv, sum(am > 0), numel(am)), 'HorizontalAlignment', 'center', 'FontSize', 7);
        stats(end + 1, :) = {reg, "other_switch_minus_late", string(wl{w}), mean(v, 'omitnan'), pv, sum(v(ok) > 0), sum(ok), sum(am > 0), numel(am)}; %#ok<AGROW>
    end
    yline(0, 'k:');
    set(gca, 'XTick', 1:numel(wins), 'XTickLabel', wl); xlim([0.5 numel(wins) + 0.5]);
    ylabel('First lap - laps 6+: other-track fraction');
    title(sprintf('%s: excess old-track cycles on first laps', reg), 'FontWeight', 'normal');
    % D: trial model - does late-grey bias toward the other track predict other-track cycles 1-5?
    nexttile(tl); hold on;
    tm = trial_model(x0, P, reg);
    stats = [stats; tm.rows]; %#ok<AGROW>
    names = tm.names; est = tm.est; ci = tm.ci;
    for k = 1:numel(names)
        plot(ci(k, :), [k k], '-', 'Color', [0.5 0.5 0.5]);
        plot(est(k), k, 'ko', 'MarkerFaceColor', 'k');
    end
    xline(0, 'k:');
    set(gca, 'YTick', 1:numel(names), 'YTickLabel', names, 'YDir', 'reverse', 'TickLabelInterpreter', 'none');
    ylim([0.5 numel(names) + 0.5]);
    xlabel('GLME log odds of an other-track cycle (cycles 1-5)');
    title(sprintf('%s: what predicts old-track cycles after onset', reg), 'FontWeight', 'normal');
end
title(tl, sprintf('Theta-cycle decoding around track onset (multinomial Bayes, held-out blocks; running > %d cm/s; >= %d spikes)', O.speed, O.min_spikes));
ar_save_figure(fig, fullfile(out_dir, 'theta_flicker_summary'));
close(fig);
T = cell2table(stats, 'VariableNames', {'region', 'measure', 'term', 'estimate', 'p', 'n_pos', 'n', 'n_anim_pos', 'n_anim'});
writetable(T, fullfile(out_dir, 'theta_flicker_stats.csv'));
end

function tm = trial_model(x0, P, reg)
% Per transition: number of other-track cycles among determined cycles 1-5;
% predictors: late-grey (last 1 s) bias toward the other track (zLO from the
% epoch table), switch, their interaction, HPC-V1 late-grey coupling and V1
% late-grey bias toward the other track (V1 animals), spikes and speed
S = load(fullfile(P.results, 'data', 'epoch_table.mat'), 'E'); E = S.E;
G = load(fullfile(P.results, 'grey_coupling', 'grey_coupling.mat'), 'R'); cpl = G.R.bayes_multinomial.trial_coupling;
x = x0(ismember(x0.cycle, 1:5), :);
[g, keys] = findgroups(x(:, {'session', 'animal', 'lap', 'track', 'is_switch'}));
keys.n_other = splitapply(@(v) sum(v == 1), x.other, g);
keys.n_det = splitapply(@(v) sum(isfinite(v)), x.other, g);
keys.spk = splitapply(@(v) log1p(mean(v)), x.n_spikes, g);
keys.spd = splitapply(@mean, x.speed, g);
keys = keys(keys.n_det > 0, :);
e = E(E.decoder == "bayes_multinomial" & E.epoch == "grey_late", :);
for rr = ["HPC", "V1"]
    er = e(e.region == rr, :);
    [~, ia, ib] = intersect(keys(:, {'session', 'lap'}), er(:, {'session', 'lap'}));
    v = NaN(height(keys), 1); v(ia) = -er.zLO(ib) .* (2 * (keys.track(ia) == 1) - 1);   % toward the other track
    keys.("grey_other_" + rr) = v;
end
[~, ia, ib] = intersect(keys(:, {'session', 'lap'}), cpl(:, {'session', 'lap'}));
keys.coupling = NaN(height(keys), 1); keys.coupling(ia) = cpl.coupling_late(ib);
keys.is_switch = double(keys.is_switch);
keys.animal_c = categorical(keys.animal); keys.sess_c = categorical(keys.session);
for v = ["spk", "spd"]
    keys.(v) = (keys.(v) - mean(keys.(v), 'omitnan')) ./ std(keys.(v), 'omitnan');
end
own = "grey_other_" + reg;
keys.own = keys.(own);
oreg = setdiff(["HPC", "V1"], reg);
keys.oth = keys.("grey_other_" + oreg);
f1 = 'n_other ~ own * is_switch + spk + spd + (1|animal_c) + (1|sess_c)';
k1 = keys(all(isfinite([keys.own keys.spk keys.spd]), 2), :);
m1 = fitglme(k1, f1, 'Distribution', 'Binomial', 'BinomialSize', k1.n_det);
tm.names = {}; tm.est = []; tm.ci = []; tm.rows = {};
tm = collect(tm, m1, {'own', 'is_switch', 'own:is_switch'}, {'late-grey bias to other track', 'first lap of block', 'bias x first lap'}, reg, "glme_main");
% coupling and the other region's bias (conflict hypothesis): V1 sessions only
k2 = keys(all(isfinite([keys.own keys.oth keys.coupling keys.spk keys.spd]), 2), :);
if height(k2) > 50
    f2 = 'n_other ~ own * coupling + oth * coupling + is_switch + spk + spd + (1|animal_c) + (1|sess_c)';
    m2 = fitglme(k2, f2, 'Distribution', 'Binomial', 'BinomialSize', k2.n_det);
    tm = collect(tm, m2, {'coupling', 'own:coupling', 'coupling:oth', 'oth'}, ...
        {'HPC-V1 late-grey coupling', 'own bias x coupling', sprintf('%s bias x coupling', oreg), ...
        sprintf('%s late-grey bias to other', oreg)}, reg, "glme_conflict");
end
end

function tm = collect(tm, m, terms, labels, reg, model)
c = m.Coefficients; ci = coefCI(m);
for k = 1:numel(terms)
    j = strcmp(c.Name, terms{k});
    if ~any(j) && contains(terms{k}, ':')     % fitglme may order interaction terms either way
        ab = split(terms{k}, ':'); j = strcmp(c.Name, [ab{2} ':' ab{1}]);
    end
    if ~any(j), continue; end
    tm.names{end + 1} = labels{k}; tm.est(end + 1, 1) = c.Estimate(j); tm.ci(end + 1, :) = ci(j, :);
    tm.rows(end + 1, :) = {reg, model, string(terms{k}), c.Estimate(j), c.pValue(j), NaN, m.NumObservations, NaN, NaN};
end
end

function v = ifelse(c, a, b)
if c, v = a; else, v = b; end
end

%% ------------------------------------------------------------------------
function plot_examples(x, Ls, P, out_dir)
% Rows: up to 8 first laps of blocks (left) and 8 late laps (right); bars =
% per-cycle log odds toward the lap's track (green toward it, red toward the
% other track; open = < 3 spikes), cycles -6..10 around onset
T = Ls.laps;
sw = unique(x.lap(x.lap_in_block == 1 & x.block > 1));
la = unique(x.lap(x.lap_in_block >= 6));
rs = RandStream('mrg32k3a', 'Seed', x.session(1));
if numel(la) > 8, la = sort(la(randperm(rs, numel(la), 8))); end
sw = sw(1:min(8, numel(sw)));
fig = figure('Position', [50 50 1500 1100], 'Color', 'w');
tl = tiledlayout(8, 2, 'TileSpacing', 'tight', 'Padding', 'compact');
lim = prctile(abs(x.LO_cur), 98);
for row = 1:8
    for colI = 1:2
        nexttile(tl, (row - 1) * 2 + colI); hold on;
        laps = {sw, la}; laps = laps{colI};
        if row > numel(laps), axis off; continue; end
        k = laps(row);
        xk = x(x.lap == k & x.cycle >= -6 & x.cycle <= 10, :);
        for j = 1:height(xk)
            v = max(min(xk.LO_cur(j), lim), -lim);
            c = [0.15 0.6 0.3]; if v < 0, c = [0.85 0.2 0.2]; end
            if xk.n_spikes(j) < 3, fc = 'none'; else, fc = c; end
            bar(xk.cycle(j), v, 0.8, 'FaceColor', fc, 'EdgeColor', c);
        end
        xline(0.5, 'k-'); yline(0, 'k:'); ylim([-lim lim]); xlim([-6.6 10.6]);
        set(gca, 'XTick', -6:2:10, 'TickDir', 'out', 'Box', 'off', 'FontSize', 7);
        if row < min(8, numel(laps)), set(gca, 'XTickLabel', []); end
        trk = {'T1 (L)', 'T2 (R)'};
        title(sprintf('lap %d, %s, lap %d of block', k, trk{T.track(k)}, T.lap_in_block(k)), 'FontWeight', 'normal', 'FontSize', 7);
    end
end
title(tl, sprintf('Session %d (%s %s), HPC per-theta-cycle log odds toward the lap''s track: left = first laps of blocks, right = laps 6+', ...
    x.session(1), Ls.subject, Ls.session), 'FontSize', 9);
xlabel(tl, 'Theta cycle from onset (vertical line = onset; green = current track, red = other; open = < 3 spikes)', 'FontSize', 8);
ar_save_figure(fig, fullfile(out_dir, 'examples', sprintf('session_%02d', x.session(1))));
close(fig);
end
