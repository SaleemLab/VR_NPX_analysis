function C = ar_single_cell()
% Single-cell counterpart of the population track bias around lap transitions.
% Per cell (spatial cells of RUN1_decoding_input.mat):
%   track preference (in-lap rate maps, running > 5 cm/s, all laps):
%     pref = (r_T1 - r_T2) / (r_T1 + r_T2) over the whole track, and pref_start
%     over 0-20 cm
%   per transition: firing rate in end_prev (-1..0 s from lap end), grey_early
%     (0..1 s after lap end), grey_late (-1..0 s before onset), onset (0..0.2 s),
%     onset_late (0.2..1 s)
%   grey preference: d' of grey rates after a T1 lap vs after a T2 lap (all
%     transitions), raw and after removing a slow trend (centred moving mean over
%     31 transitions, ~3 blocks: keeps block alternation, removes slower drift)
%   onset switch effect: for laps of each track X, rate on the first lap of a
%     block minus the mean over laps 6+, signed toward the old track: per cell and
%     track, paired with the cell's preference for the old (previous) track
%   switch step (drift-robust): rate at transition i = 1 (after the first lap of
%     the new block) minus i = 0 (switch transition), signed toward the new track,
%     vs a baseline of consecutive within-block transitions (signed toward the
%     block track); a context effect makes the step, but not the baseline,
%     correlate with the in-track preference
% Statistics: Spearman across cells (pooled, and per session -> signed rank /
% animals), within region (HPC, V1_L, V1_R; V1 excludes P.V1_exclude_subjects).
% Output: results\single_cell\cells.mat (table C, one row per cell) and figures.
P = ar_setup();
S = load(fullfile(P.inputs, 'RUN1_decoding_input.mat'), 'D'); D = S.D;
Lt = load(fullfile(P.results, 'data', 'lap_tables.mat'), 'L'); L = Lt.L;
out_dir = fullfile(P.results, 'single_cell');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
epochs = {'end_prev', 'grey_early', 'grey_late', 'onset', 'onset_late'};
rows = {};
for s = 1:numel(D)
    Ds = D(s); T = L(s).laps;
    nC = numel(Ds.cell_id);
    [pref, pref_start, rate_mean, peak_pos] = track_preference(Ds, T, P);
    trans = find(T.valid_transition);
    iv = zeros(numel(trans), 2, numel(epochs));
    iv(:, :, 1) = [T.t_off(trans - 1) - 1, T.t_off(trans - 1)];
    iv(:, :, 2) = [T.t_off(trans - 1), T.t_off(trans - 1) + 1];
    iv(:, :, 3) = [T.t_on(trans) - 1, T.t_on(trans)];
    iv(:, :, 4) = [T.t_on(trans), T.t_on(trans) + 0.2];
    iv(:, :, 5) = [T.t_on(trans) + 0.2, T.t_on(trans) + 1];
    R = NaN(nC, numel(trans), numel(epochs));
    [~, ci] = ismember(Ds.spike_id, Ds.cell_id);
    for c = 1:nC
        st = Ds.spike_times(ci == c);
        for e = 1:numel(epochs)
            R(c, :, e) = (count_in(st, iv(:, 2, e)) - count_in(st, iv(:, 1, e))) ./ (iv(:, 2, e) - iv(:, 1, e));
        end
    end
    prev = T.prev_track(trans); trk = T.track(trans); lib = T.lap_in_block(trans);
    % switch-locked pairs (transition i = 0 -> i = 1 after each switch) and
    % within-block consecutive pairs (both in the same block, i >= 1), as indices
    % into trans, with the sign of the (new / block) track (+1 = T1)
    [~, pos_in] = ismember((1:height(T))', trans);
    pairs_sw = []; pairs_base = [];
    for k = trans'
        k1 = k + 1;
        if k1 > height(T) || pos_in(k1) == 0, continue; end
        sgn = 2 * (T.track(k) == 1) - 1;
        if T.is_switch(k) && T.block(k1) == T.block(k)
            pairs_sw(end + 1, :) = [pos_in(k) pos_in(k1) sgn]; %#ok<AGROW>
        elseif ~T.is_switch(k) && T.block(k1) == T.block(k)
            pairs_base(end + 1, :) = [pos_in(k) pos_in(k1) sgn]; %#ok<AGROW>
        end
    end
    for c = 1:nC
        row = {s, L(s).animal_index, string(Ds.subject), Ds.cell_id(c), string(Ds.region(c)), rate_mean(c), pref(c), ...
            pref_start(c), peak_pos(c)};
        for e = 1:numel(epochs)
            x = R(c, :, e)';
            xd = x - movmean(x, 31);
            row = [row, {dprime(x(prev == 1), x(prev == 2)), dprime(xd(prev == 1), xd(prev == 2)), ...
                ranksum_p(x(prev == 1), x(prev == 2))}]; %#ok<AGROW>
            % switch effect per track X, signed toward the old track
            sw = NaN(1, 2);
            for X = 1:2
                a = x(trk == X & lib == 1); b = x(trk == X & lib >= 6);
                if numel(a) >= 2 && numel(b) >= 2
                    sw(X) = (mean(a) - mean(b)) / max(sqrt((var(a) + var(b)) / 2), 0.5);
                end
            end
            % switch step (i = 1 minus i = 0) and within-block consecutive baseline,
            % signed toward the new / block track, in units of the cell's SD
            sd = max(std(x, 'omitnan'), 0.1);
            step = mean(pairs_sw(:, 3) .* (x(pairs_sw(:, 2)) - x(pairs_sw(:, 1))), 'omitnan') / sd;
            base = mean(pairs_base(:, 3) .* (x(pairs_base(:, 2)) - x(pairs_base(:, 1))), 'omitnan') / sd;
            row = [row, {sw(1), sw(2), step, base}]; %#ok<AGROW>
        end
        rows(end + 1, :) = row; %#ok<AGROW>
    end
    fprintf('Session %d: %d cells\n', s, nC);
end
names = {'session', 'animal', 'subject', 'cell_id', 'region', 'rate', 'pref', 'pref_start', 'peak_pos'};
for e = 1:numel(epochs)
    names = [names, strcat(epochs{e}, {'_d_prev', '_d_prev_detrend', '_p_prev', '_sw_T1lap', '_sw_T2lap', ...
        '_swstep', '_swbase'})]; %#ok<AGROW>
end
C = cell2table(rows, 'VariableNames', names);
save(fullfile(out_dir, 'cells.mat'), 'C');
plot_single_cell(C, P, epochs, out_dir);
end

function n = count_in(st, t)
% number of spikes before each time in t (binary search on sorted spike times)
n = discretize(t, [-Inf; unique(st(:)); Inf]) - 1;
end

function d = dprime(a, b)
a = a(isfinite(a)); b = b(isfinite(b));
if numel(a) < 3 || numel(b) < 3, d = NaN; return; end
s = sqrt((var(a) + var(b)) / 2);
if s == 0, d = 0; else, d = (mean(a) - mean(b)) / s; end
end

function p = ranksum_p(a, b)
a = a(isfinite(a)); b = b(isfinite(b));
if numel(a) < 3 || numel(b) < 3, p = NaN; return; end
p = ranksum(a, b);
end

function [pref, pref_start, rate_mean, peak_pos] = track_preference(Ds, T, P)
tvec = Ds.tvec; nLaps = height(T); nPos = numel(P.pos_edges) - 1;
lap = zeros(size(tvec));
for k = 1:nLaps, lap(tvec >= T.t_on(k) & tvec < T.t_off(k)) = k; end
pb = discretize(Ds.position, P.pos_edges);
run = lap > 0 & ~isnan(pb) & Ds.speed > P.bayes_speed;
tr = zeros(size(tvec)); tr(lap > 0) = T.track(lap(lap > 0));
dt = [diff(tvec); median(diff(tvec))];
occ = accumarray([tr(run) pb(run)], dt(run), [2 nPos]);
[~, ci] = ismember(Ds.spike_id, Ds.cell_id);
ok = ci > 0 & Ds.spike_times >= tvec(1) & Ds.spike_times <= tvec(end);
smp = interp1(tvec, (1:numel(tvec))', Ds.spike_times(ok), 'nearest');
ci = ci(ok); keep = run(smp);
cnt = accumarray([ci(keep) tr(smp(keep)) pb(smp(keep))], 1, [numel(Ds.cell_id) 2 nPos]);
r_all = sum(cnt, 3) ./ sum(occ, 2)';
start = 1:2;   % 0-20 cm
r_st = sum(cnt(:, :, start), 3) ./ sum(occ(:, start), 2)';
pref = (r_all(:, 1) - r_all(:, 2)) ./ (r_all(:, 1) + r_all(:, 2));
pref_start = (r_st(:, 1) - r_st(:, 2)) ./ (r_st(:, 1) + r_st(:, 2));
rate_mean = mean(r_all, 2);
% peak position (cm, bin centre) on the cell's preferred track
rate = cnt ./ reshape(occ, [1 2 nPos]);
pt = 1 + (pref < 0);
peak_pos = NaN(numel(pref), 1);
for c = 1:numel(pref)
    [~, k] = max(squeeze(rate(c, pt(c), :)));
    peak_pos(c) = (k - 0.5) * (P.pos_edges(2) - P.pos_edges(1));
end
end

function plot_single_cell(C, P, epochs, out_dir)
% Rows: regions. Columns: grey preference vs track preference (grey_early,
% grey_late; raw and detrended), and onset switch effect vs old-track preference.
regions = {'HPC', 'V1_L', 'V1_R'};
bad = find(ismember(P.subjects, P.V1_exclude_subjects));
acol = lines(6);
fig = figure('Position', [50 50 1700 330 * numel(regions)], 'Color', 'w');
tl = tiledlayout(numel(regions), 5, 'TileSpacing', 'compact');
stats = {};
for r = 1:numel(regions)
    x = C(contains(C.region, regions{r}), :);
    if startsWith(regions{r}, 'V1'), x(ismember(x.animal, bad), :) = []; end
    panels = {'grey_early_d_prev', 'grey_late_d_prev', 'grey_late_d_prev_detrend', 'end_prev_d_prev'};
    for k = 1:numel(panels)
        nexttile(tl); hold on;
        y = x.(panels{k});
        scatter(x.pref, y, 8, acol(x.animal, :), 'filled', 'MarkerFaceAlpha', 0.5);
        [rho, p] = corr(x.pref, y, 'Type', 'Spearman', 'Rows', 'complete');
        sess = session_rho(x, 'pref', panels{k});
        stats(end + 1, :) = {string(regions{r}), string(panels{k}), height(x), rho, p, mean(sess.rho, 'omitnan'), ...
            signrank_p(sess.rho), sum(sess.rho > 0), sum(isfinite(sess.rho)), animals_pos(sess)}; %#ok<AGROW>
        xline(0, 'k:'); yline(0, 'k:');
        xlabel('In-track preference (T1 - T2)'); ylabel('d'' grey rate after T1 vs T2 lap');
        title(sprintf('%s %s: rho %.2f (p %.2g); sessions %.2f, p %.2g', regions{r}, panels{k}, rho, p, ...
            mean(sess.rho, 'omitnan'), signrank_p(sess.rho)), 'Interpreter', 'none', 'FontWeight', 'normal', 'FontSize', 8);
    end
    % onset switch effect vs preference for the old track (both lap types pooled)
    nexttile(tl); hold on;
    y = [x.onset_sw_T1lap; x.onset_sw_T2lap];
    old_pref = [-x.pref; x.pref];   % for a T1 lap the old track is T2
    xa = [x; x];
    scatter(old_pref, y, 8, acol(xa.animal, :), 'filled', 'MarkerFaceAlpha', 0.5);
    [rho, p] = corr(old_pref, y, 'Type', 'Spearman', 'Rows', 'complete');
    xa.old_pref = old_pref; xa.sw = y;
    sess = session_rho(xa, 'old_pref', 'sw');
    stats(end + 1, :) = {string(regions{r}), "onset_switch_vs_old_pref", height(xa), rho, p, mean(sess.rho, 'omitnan'), ...
        signrank_p(sess.rho), sum(sess.rho > 0), sum(isfinite(sess.rho)), animals_pos(sess)}; %#ok<AGROW>
    xline(0, 'k:'); yline(0, 'k:');
    xlabel('Preference for the old (previous block) track'); ylabel('Onset rate: first lap - laps 6+ (d'')');
    title(sprintf('%s onset: rho %.2f (p %.2g); sessions p %.2g', regions{r}, rho, p, signrank_p(sess.rho)), ...
        'FontWeight', 'normal', 'FontSize', 8);
end
title(tl, 'Single cells: grey-screen and onset firing vs track preference (colour = animal)');
ar_save_figure(fig, fullfile(out_dir, 'single_cell_summary'), true);

% Switch step vs within-block baseline: session rho between the step and the
% in-track preference (dots = sessions, colour = animal)
meas = {'end_prev', 'grey_early', 'grey_late', 'onset', 'onset_late'};
fig = figure('Position', [50 50 1500 330 * numel(regions)], 'Color', 'w');
tl = tiledlayout(numel(regions), 1, 'TileSpacing', 'compact');
for r = 1:numel(regions)
    x = C(contains(C.region, regions{r}), :);
    if startsWith(regions{r}, 'V1'), x(ismember(x.animal, bad), :) = []; end
    nexttile(tl); hold on;
    for k = 1:numel(meas)
        for b = 1:2
            v = [meas{k} '_' char(string(ifelse(b == 1, 'swstep', 'swbase')))];
            sess = session_rho(x, 'pref', v);
            [rho, p] = corr(x.pref, x.(v), 'Type', 'Spearman', 'Rows', 'complete');
            stats(end + 1, :) = {string(regions{r}), string(v), height(x), rho, p, mean(sess.rho, 'omitnan'), ...
                signrank_p(sess.rho), sum(sess.rho > 0), sum(isfinite(sess.rho)), animals_pos(sess)}; %#ok<AGROW>
            xx = k + (b - 1.5) * 0.35;
            scatter(xx + 0.05 * randn(height(sess), 1), sess.rho, 15, acol(sess.animal, :), 'filled', 'MarkerFaceAlpha', 0.6);
            c = ifelse(b == 1, [0 0 0], [0.6 0.6 0.6]);
            errorbar(xx, mean(sess.rho, 'omitnan'), std(sess.rho, 'omitnan') / sqrt(sum(isfinite(sess.rho))), 'o', ...
                'Color', c, 'MarkerFaceColor', c, 'CapSize', 0, 'LineWidth', 1.5);
            text(xx, 0.6, sprintf('%.2g', signrank_p(sess.rho)), 'HorizontalAlignment', 'center', 'FontSize', 7);
        end
    end
    yline(0, 'k:');
    set(gca, 'XTick', 1:numel(meas), 'XTickLabel', meas, 'TickLabelInterpreter', 'none', 'TickDir', 'out', 'Box', 'off');
    ylabel('Session rho (step vs preference)');
    title(sprintf('%s: switch step i=1 - i=0 (black) vs within-block baseline (grey); signed-rank p above', regions{r}), ...
        'FontWeight', 'normal', 'Interpreter', 'none');
end
ar_save_figure(fig, fullfile(out_dir, 'single_cell_switch_step'));

% Field-location split of the switch step: cells with the preferred-track peak
% at the track start (< 40 cm), middle (40-100 cm) or end (>= 100 cm)
loc = {[0 40], [40 100], [100 140]}; loc_name = {'start', 'middle', 'end'};
for r = 1:numel(regions)
    x0 = C(contains(C.region, regions{r}), :);
    if startsWith(regions{r}, 'V1'), x0(ismember(x0.animal, bad), :) = []; end
    for g = 1:numel(loc)
        x = x0(x0.peak_pos >= loc{g}(1) & x0.peak_pos < loc{g}(2), :);
        for v = {'grey_early_swstep', 'grey_late_swstep', 'grey_early_swbase', 'grey_late_swbase', 'onset_swstep'}
            sess = session_rho(x, 'pref', v{1});
            [rho, p] = corr(x.pref, x.(v{1}), 'Type', 'Spearman', 'Rows', 'complete');
            stats(end + 1, :) = {string(regions{r}), string(v{1}) + "_" + loc_name{g}, height(x), rho, p, ...
                mean(sess.rho, 'omitnan'), signrank_p(sess.rho), sum(sess.rho > 0), sum(isfinite(sess.rho)), animals_pos(sess)}; %#ok<AGROW>
        end
    end
end
T = cell2table(stats, 'VariableNames', {'region', 'measure', 'n_cells', 'rho_pooled', 'p_pooled', ...
    'rho_session_mean', 'p_session_signrank', 'n_sess_pos', 'n_sess', 'animals_pos'});
writetable(T, fullfile(out_dir, 'single_cell_stats.csv'));
end

function v = ifelse(c, a, b)
if c, v = a; else, v = b; end
end

function sess = session_rho(x, a, b)
us = unique(x.session);
rho = NaN(numel(us), 1); an = NaN(numel(us), 1);
for i = 1:numel(us)
    xi = x(x.session == us(i), :);
    if sum(isfinite(xi.(a)) & isfinite(xi.(b))) >= 8
        rho(i) = corr(xi.(a), xi.(b), 'Type', 'Spearman', 'Rows', 'complete');
    end
    an(i) = xi.animal(1);
end
sess = table(us, an, rho, 'VariableNames', {'session', 'animal', 'rho'});
end

function p = signrank_p(x)
x = x(isfinite(x));
if numel(x) < 2, p = NaN; else, p = signrank(x); end
end

function s = animals_pos(sess)
[ua, ~, ia] = unique(sess.animal);
m = accumarray(ia, sess.rho, [], @(v) mean(v, 'omitnan'));
s = sprintf('%d/%d', sum(m > 0), sum(isfinite(m)));
end
