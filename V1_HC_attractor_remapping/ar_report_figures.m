function ar_report_figures()
% Essential figures for the results report (results\report\):
%   fig1_drift       decoders and drift: lobo3 vs lobo5 accuracy; grey bias of
%                    switch vs non-switch transitions (mirror image); drift-
%                    cancelled grey bias across decoders
%   fig3_grey        HPC grey persistence: switch-locked curves (all-laps
%                    decoder), switch step across decoders, single cells by field
%   fig4_onset       first-lap enhancement: onset bias and firing vs lap in
%                    block, switch - late (spike/speed-controlled), V1 single cells
%   fig5_trial       trial-to-trial prediction: LME coefficients across models
%                    and decoders, per-session partial correlation
% Session dots coloured by animal, diamonds = animal means, black = mean +- SEM.
P = ar_setup();
out_dir = fullfile(P.results, 'report');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
S = load(fullfile(P.results, 'data', 'epoch_table.mat'), 'E'); E = S.E;
Lt = load(fullfile(P.results, 'data', 'lap_tables.mat'), 'L'); L = Lt.L;
bad = find(ismember(P.subjects, P.V1_exclude_subjects));
E(startsWith(E.region, "V1") & ismember(E.animal, bad), :) = [];
E.sgn_prev = 2 * (E.prev_track == 1) - 1;
E.sgn_next = 2 * (E.track == 1) - 1;
acol = lines(6);
set(groot, 'defaultAxesFontSize', 10, 'defaultAxesTickDir', 'out', 'defaultAxesBox', 'off');

fig1(E, P, acol, out_dir);
fig3(E, L, P, acol, out_dir);
fig4(E, P, acol, out_dir);
fig5(P, acol, out_dir);
end

%% ------------------------------------------------------------------------
function fig1(E, P, acol, out_dir)
Q = readtable(fullfile(P.results, 'qc', 'session_qc.csv'), 'TextType', 'string');
fig = figure('Position', [50 50 1500 430], 'Color', 'w');
tl = tiledlayout(1, 3, 'TileSpacing', 'loose');
% A: lobo3 vs lobo5 in-lap accuracy
nexttile(tl); hold on;
mk = {'o', 's'}; regs = ["HPC", "V1"];
for r = 1:2
    q = Q(Q.region == regs(r) & Q.decoder == "bayes_multinomial", :);
    if regs(r) == "V1", q(ismember(q.animal, find(ismember(P.subjects, P.V1_exclude_subjects))), :) = []; end
    for i = 1:height(q)
        h = scatter(q.acc_rest_lobo3(i), q.acc_rest_lobo5(i), 36, acol(q.animal(i), :), mk{r}, 'LineWidth', 1.2);
        if r == 1, h.MarkerFaceColor = acol(q.animal(i), :); end
    end
end
plot([0.3 1], [0.3 1], 'k:'); yline(0.5, ':', 'Color', [0.6 0.6 0.6]);
axis square; xlim([0.3 1]); ylim([0.3 1]);
xlabel('Track accuracy, held out b\pm1 (lobo3)'); ylabel('Track accuracy, held out b\pm2 (lobo5)');
title({'A  Training blocks further away in time', 'decode worse (filled HPC, open V1)'}, 'FontWeight', 'normal');
% B: HPC early grey, switch vs non-switch, per decoder
decs = ["bayes_multinomial", "bayes_multinomial_lobo5", "bayes_multinomial_all"];
dlab = {'lobo3', 'lobo5', 'all laps'};
nexttile(tl); hold on;
for d = 1:3
    x = E(E.decoder == decs(d) & E.region == "HPC" & E.epoch == "grey_early", :);
    [sw, ns, an] = switch_nonswitch(x);
    jitter_dots(d - 0.18, ns, an, acol); jitter_dots(d + 0.18, sw, an, acol);
    errorbar(d - 0.18, mean(ns), sem(ns), 'ko', 'MarkerFaceColor', 'k', 'CapSize', 0, 'LineWidth', 1.5);
    errorbar(d + 0.18, mean(sw), sem(sw), 'o', 'Color', [0.85 0.2 0.2], 'MarkerFaceColor', [0.85 0.2 0.2], 'CapSize', 0, 'LineWidth', 1.5);
end
yline(0, 'k:');
set(gca, 'XTick', 1:3, 'XTickLabel', dlab); xlim([0.4 3.6]);
ylabel('HPC early grey: bias toward previous track (z)');
title({'B  Switch (red) and non-switch (black) transitions', 'swap places between lobo3 and lobo5'}, 'FontWeight', 'normal');
% C: drift-cancelled (switch + non-switch) / 2
nexttile(tl); hold on;
conds = {"HPC", "grey_early"; "HPC", "grey_late"; "V1", "grey_early"; "V1", "grey_late"};
for k = 1:size(conds, 1)
    for d = 1:3
        x = E(E.decoder == decs(d) & E.region == conds{k, 1} & E.epoch == conds{k, 2}, :);
        [sw, ns, an] = switch_nonswitch(x);
        v = (sw + ns) / 2;
        xx = k + (d - 2) * 0.25;
        jitter_dots(xx, v, an, acol, 0.04);
        errorbar(xx, mean(v), sem(v), 'ko', 'MarkerFaceColor', 'k', 'CapSize', 0, 'LineWidth', 1.5);
        if d == 1, text(xx, -0.32, dlab{d}, 'FontSize', 7, 'Rotation', 90, 'HorizontalAlignment', 'left'); end
        if d > 1, text(xx, -0.32, dlab{d}, 'FontSize', 7, 'Rotation', 90, 'HorizontalAlignment', 'left'); end
    end
end
yline(0, 'k:');
set(gca, 'XTick', 1:4, 'XTickLabel', {'HPC early', 'HPC late', 'V1 early', 'V1 late'}); xlim([0.5 4.5]); ylim([-0.35 0.6]);
ylabel('Bias toward previous track, (switch + non-switch)/2 (z)');
title({'C  Averaging the two cancels decoder drift:', 'all decoders agree'}, 'FontWeight', 'normal');
ar_save_figure(fig, fullfile(out_dir, 'fig1_drift'));
close(fig);
end

function [sw, ns, an] = switch_nonswitch(x)
% session means of the centred bias toward the previous track for switch and
% non-switch transitions (centre = mean of previous-T1 and previous-T2 means)
us = unique(x.session);
sw = NaN(numel(us), 1); ns = sw; an = sw;
for i = 1:numel(us)
    xs = x(x.session == us(i), :);
    c = mean([mean(xs.zLO(xs.prev_track == 1), 'omitnan') mean(xs.zLO(xs.prev_track == 2), 'omitnan')]);
    tp = (xs.zLO - c) .* xs.sgn_prev;
    sw(i) = mean(tp(xs.is_switch), 'omitnan'); ns(i) = mean(tp(~xs.is_switch), 'omitnan');
    an(i) = xs.animal(1);
end
end

%% ------------------------------------------------------------------------
function fig3(E, L, P, acol, out_dir)
idx = -5:6;
fig = figure('Position', [50 50 1900 430], 'Color', 'w');
tl = tiledlayout(1, 4, 'TileSpacing', 'loose');
% A: switch-locked curves, all-laps decoder (grey is never in its training data)
nexttile(tl); hold on;
ep = ["grey_early", "grey_late"]; cl = [0 0 0; 0.6 0.6 0.6];
for e = 1:2
    x = E(E.decoder == "bayes_multinomial_all" & E.region == "HPC" & E.epoch == ep(e), :);
    [Y, ~] = switch_curves(x, L, idx);
    m = mean(Y, 1, 'omitnan'); se = std(Y, 0, 1, 'omitnan') ./ sqrt(sum(isfinite(Y), 1));
    errorbar(idx + (e - 1.5) * 0.12, m, se, 'o-', 'Color', cl(e, :), 'MarkerFaceColor', cl(e, :), 'CapSize', 0, 'LineWidth', 1.5);
end
xline(0.5, 'k--'); yline(0, 'k:');
xlabel('Lap transition relative to the block switch'); ylabel('HPC bias toward the new track (z)');
legend({'early grey (0-1 s after lap end)', 'late grey (1 s before onset)'}, 'Box', 'off', 'Location', 'northwest');
title({'A  HPC grey flips after one new lap', 'all-laps decoder; 0 = old -> first new lap'}, 'FontWeight', 'normal');
% B: step 1 - 0 for HPC early grey across decoders, animals
nexttile(tl); hold on;
decs = ["bayes_multinomial", "bayes_multinomial_lobo5", "bayes_multinomial_all", "pls", "pls_lobo5", "pls_all"];
dlab = {'Bayes lobo3', 'Bayes lobo5', 'Bayes all', 'PLS lobo3', 'PLS lobo5', 'PLS all'};
for d = 1:numel(decs)
    for e = 1:2
        x = E(E.decoder == decs(d) & E.region == "HPC" & E.epoch == ep(e), :);
        [Y, an] = switch_curves(x, L, idx);
        step = Y(:, idx == 1) - Y(:, idx == 0);
        xx = d + (e - 1.5) * 0.3;
        [ua, ~, ia] = unique(an);
        am = accumarray(ia, step, [], @(v) mean(v, 'omitnan'));
        scatter(xx + zeros(size(am)), am, 30, acol(ua, :), 'd', 'filled');
        errorbar(xx, mean(step, 'omitnan'), sem(step), 'o', 'Color', cl(e, :), 'MarkerFaceColor', cl(e, :), 'CapSize', 0, 'LineWidth', 1.5);
    end
end
yline(0, 'k:');
set(gca, 'XTick', 1:numel(decs), 'XTickLabel', dlab); xtickangle(30); xlim([0.4 numel(decs) + 0.6]);
ylabel('Step: after first new lap - before switch (z)');
title({'B  Step for every decoder', 'black = early, grey = late; diamonds = animals'}, 'FontWeight', 'normal');
% C: single cells, HPC: session rho of step vs track preference by field location
Cs = load(fullfile(P.results, 'single_cell', 'cells.mat'), 'C'); C = Cs.C;
x0 = C(contains(C.region, 'HPC'), :);
nexttile(tl); hold on;
groups = {[0 140], [0 40], [40 100], [100 140]}; glab = {'all', 'start field', 'middle', 'end field'};
for g = 1:numel(groups)
    x = x0(x0.peak_pos >= groups{g}(1) & x0.peak_pos < groups{g}(2) + (g == 1), :);
    for b = 1:2
        v = 'grey_early_swstep'; c = [0 0 0];
        if b == 2, v = 'grey_early_swbase'; c = [0.6 0.6 0.6]; end
        [rho, an] = session_rho(x, 'pref', v);
        xx = g + (b - 1.5) * 0.3;
        jitter_dots(xx, rho, an, acol, 0.05);
        errorbar(xx, mean(rho, 'omitnan'), sem(rho), 'o', 'Color', c, 'MarkerFaceColor', c, 'CapSize', 0, 'LineWidth', 1.5);
    end
end
yline(0, 'k:');
set(gca, 'XTick', 1:4, 'XTickLabel', glab); xlim([0.5 4.5]);
ylabel('Session rho: firing change vs track preference');
title({'C  Single HPC cells (end-of-track cells)', 'black = switch step, grey = baseline'}, 'FontWeight', 'normal');
% D: switch step by time after lap end (visual latency vs persistence),
% all-laps decoder, lap-end-locked 100 ms bins
S = load(fullfile(P.results, 'data', 'epoch_table.mat'), 'TC'); TC = S.TC;
bad = find(ismember(P.subjects, P.V1_exclude_subjects));
wins = {[0 0.2], [0.2 0.5], [0.5 1]}; wl = {'0-0.2 s', '0.2-0.5 s', '0.5-1 s'};
regs = {'HPC', 'V1'}; rc = [0.15 0.45 0.25; 0.45 0.2 0.6];
nexttile(tl); hold on;
for r = 1:2
    sel = find(strcmp({TC.off.region}, regs{r}) & strcmp({TC.off.decoder}, 'bayes_multinomial_all'));
    for w = 1:numel(wins)
        cols = TC.t_off > wins{w}(1) & TC.t_off < wins{w}(2);
        stp = []; an = [];
        for i = sel
            s = TC.off(i).session; D = TC.trans(TC.trans.session == s, :);
            if strcmp(regs{r}, 'V1') && ismember(D.animal(1), bad), continue; end
            T = L(s).laps;
            v = NaN(height(T), 1); v(D.lap) = mean(double(TC.off(i).z(:, cols)), 2, 'omitnan');
            d = [];
            for k0 = find(T.is_switch)'
                if k0 + 1 > height(T) || T.block(k0 + 1) ~= T.block(k0), continue; end
                d(end + 1) = (2 * (T.track(k0) == 1) - 1) * (v(k0 + 1) - v(k0)); %#ok<AGROW>
            end
            stp(end + 1) = mean(d, 'omitnan'); an(end + 1) = D.animal(1); %#ok<AGROW>
        end
        xx = w + (r - 1.5) * 0.3;
        [ua, ~, ia] = unique(an); am = accumarray(ia(:), stp(:), [], @mean);
        scatter(xx + 0.1 + zeros(size(am)), am, 26, acol(ua, :), 'd', 'filled', 'MarkerEdgeColor', 'k');
        errorbar(xx, mean(stp), sem(stp(:)), 'o', 'Color', rc(r, :), 'MarkerFaceColor', rc(r, :), 'CapSize', 0, 'LineWidth', 1.5);
    end
end
yline(0, 'k:');
set(gca, 'XTick', 1:numel(wins), 'XTickLabel', wl); xlim([0.5 numel(wins) + 0.5]);
xlabel('Time after lap end'); ylabel('Switch step toward new track (z)');
title({'D  V1: first 200 ms only (visual latency);', 'HPC persists ~1 s (green HPC, purple V1)'}, 'FontWeight', 'normal');
ar_save_figure(fig, fullfile(out_dir, 'fig3_grey'));
close(fig);
end

function [Y, an] = switch_curves(x, L, idx)
% session x transition-index curves of zLO toward the new track around switches
us = unique(x.session);
Y = NaN(numel(us), numel(idx)); an = NaN(numel(us), 1);
for si = 1:numel(us)
    s = us(si); T = L(s).laps;
    xs = x(x.session == s, :);
    v = NaN(height(T), 1); v(xs.lap) = xs.zLO;
    k0s = find(T.is_switch);
    M = NaN(numel(k0s), numel(idx));
    for j = 1:numel(k0s)
        k0 = k0s(j); b = T.block(k0); sgn = 2 * (T.track(k0) == 1) - 1;
        for ii = 1:numel(idx)
            k = k0 + idx(ii);
            if k < 2 || k > height(T), continue; end
            if idx(ii) < 0 && ~(T.block(k) == b - 1 && T.block(k - 1) == b - 1), continue; end
            if idx(ii) > 0 && ~(T.block(k) == b && T.block(k - 1) == b), continue; end
            M(j, ii) = v(k) * sgn;
        end
    end
    Y(si, :) = mean(M, 1, 'omitnan'); an(si) = xs.animal(1);
end
end

%% ------------------------------------------------------------------------
function fig4(E, P, acol, out_dir)
fig = figure('Position', [50 50 1700 430], 'Color', 'w');
tl = tiledlayout(1, 4, 'TileSpacing', 'loose');
regs = ["HPC", "V1"]; rc = [0.15 0.45 0.25; 0.45 0.2 0.6];
lib_max = 10;
% A, B: onset bias toward the current track and onset firing vs lap in block
for panel = 1:2
    nexttile(tl); hold on;
    for r = 1:2
        x = E(E.decoder == "bayes_multinomial" & E.region == regs(r) & E.epoch == "onset", :);
        x = x(x.block > 1, :);   % first block has no switch
        us = unique(x.session); M = NaN(numel(us), lib_max);
        for i = 1:numel(us)
            xs = x(x.session == us(i), :);
            lib = min(xs.lap_in_block, lib_max);
            if panel == 1
                y = xs.zLO .* xs.sgn_next;
            else
                y = xs.spikes_per_win / mean(xs.spikes_per_win, 'omitnan') * 100;
            end
            M(i, :) = accumarray(lib, y, [lib_max 1], @(v) mean(v, 'omitnan'), NaN)';
        end
        m = mean(M, 1, 'omitnan'); se = std(M, 0, 1, 'omitnan') ./ sqrt(sum(isfinite(M), 1));
        errorbar(1:lib_max, m, se, 'o-', 'Color', rc(r, :), 'MarkerFaceColor', rc(r, :), 'CapSize', 0, 'LineWidth', 1.5);
    end
    set(gca, 'XTick', 1:lib_max, 'XTickLabel', [string(1:lib_max - 1) "10+"]); xlim([0.5 lib_max + 0.5]);
    xlabel('Lap in block (1 = first lap after a switch)');
    if panel == 1
        ylabel('Onset (0-200 ms) bias toward current track (z)');
        title({'A  Onset track coding is strongest on the first lap', 'and declines over the block'}, 'FontWeight', 'normal');
        legend(regs, 'Box', 'off');
    else
        ylabel('Onset firing (% of session mean)');
        title({'B  Firing at onset is also highest on the first lap'}, 'FontWeight', 'normal');
    end
end
% C: switch - laps 6+ with spikes and speed regressed out
nexttile(tl); hold on;
decs = ["bayes_multinomial", "pls"];
for r = 1:2
    for d = 1:2
        x = E(E.decoder == decs(d) & E.region == regs(r) & E.epoch == "onset", :);
        us = unique(x.session); v = NaN(numel(us), 1); an = v;
        for i = 1:numel(us)
            xs = x(x.session == us(i), :);
            c = mean([mean(xs.zLO(xs.prev_track == 1), 'omitnan') mean(xs.zLO(xs.prev_track == 2), 'omitnan')]);
            y = (xs.zLO - c) .* xs.sgn_next;
            X = [ones(height(xs), 1) log1p(xs.spikes_per_win) xs.speed];
            ok = all(isfinite(X), 2) & isfinite(y);
            b = X(ok, :) \ y(ok);
            y = y - X(:, 2:3) * b(2:3);
            v(i) = mean(y(xs.lap_in_block == 1), 'omitnan') - mean(y(xs.lap_in_block >= 6), 'omitnan');
            an(i) = xs.animal(1);
        end
        xx = r + (d - 1.5) * 0.35;
        jitter_dots(xx, v, an, acol, 0.05);
        errorbar(xx, mean(v, 'omitnan'), sem(v), 'ko', 'MarkerFaceColor', 'k', 'CapSize', 0, 'LineWidth', 1.5);
        text(xx, -0.55, char(strrep(decs(d), 'bayes_', '')), 'FontSize', 7, 'HorizontalAlignment', 'center');
    end
end
yline(0, 'k:');
set(gca, 'XTick', 1:2, 'XTickLabel', regs); xlim([0.4 2.6]); ylim([-0.6 1.2]);
ylabel('First lap - laps 6+, onset (z, spikes & speed removed)');
title({'C  The first-lap gain is not just more spikes', 'or faster running'}, 'FontWeight', 'normal');
% D: single cells: onset change on first lap vs preference for the old track
Cs = load(fullfile(P.results, 'single_cell', 'cells.mat'), 'C'); C = Cs.C;
bad = find(ismember(P.subjects, P.V1_exclude_subjects));
nexttile(tl); hold on;
cr = {'HPC', 'V1_L', 'V1_R'};
for k = 1:3
    x = C(contains(C.region, cr{k}), :);
    if startsWith(cr{k}, 'V1'), x(ismember(x.animal, bad), :) = []; end
    xa = [x; x]; xa.old_pref = [-x.pref; x.pref]; xa.sw = [x.onset_sw_T1lap; x.onset_sw_T2lap];
    [rho, an] = session_rho(xa, 'old_pref', 'sw');
    jitter_dots(k, rho, an, acol, 0.08);
    errorbar(k + 0.25, mean(rho, 'omitnan'), sem(rho), 'ko', 'MarkerFaceColor', 'k', 'CapSize', 0, 'LineWidth', 1.5);
end
yline(0, 'k:');
set(gca, 'XTick', 1:3, 'XTickLabel', {'HPC', 'V1 L', 'V1 R'}, 'TickLabelInterpreter', 'none'); xlim([0.5 3.5]);
ylabel('Session rho: first-lap onset change vs old-track preference');
title({'D  V1 cells tuned to the new track fire more', 'on its first lap (negative = new-track cells up)'}, 'FontWeight', 'normal');
ar_save_figure(fig, fullfile(out_dir, 'fig4_onset'));
close(fig);
end

%% ------------------------------------------------------------------------
function fig5(P, acol, out_dir)
C = readtable(fullfile(P.results, 'trial_prediction', 'trial_prediction_stats.csv'), 'TextType', 'string');
R = load(fullfile(P.results, 'trial_prediction', 'trial_prediction.mat'), 'R'); R = R.R;
fig = figure('Position', [50 50 1600 460], 'Color', 'w');
tl = tiledlayout(1, 3, 'TileSpacing', 'loose');
decs = ["bayes_multinomial", "bayes_multinomial_lobo5", "pls", "pls_lobo5"];
dmk = {'o', 's', '^', 'v'};
terms = {"main_onset", "grey_late", 'late grey -> 0-200 ms (main)'; ...
    "main_onset_100", "grey_late", 'late grey -> 0-100 ms'; ...
    "main_onset_200", "grey_late", 'late grey -> 100-200 ms'; ...
    "main_onset_late", "grey_late", 'late grey -> 0.2-1 s'; ...
    "plus_early_end", "grey_late", 'late grey, + early grey & lap end'; ...
    "plus_early_end", "grey_early", 'early grey (same model)'; ...
    "plus_early_end", "end_prev", 'lap end (same model)'; ...
    "neighbour_control", "grey_late_prev_transition", 'late grey of previous trial'; ...
    "neighbour_control", "grey_late_next_transition", 'late grey of next trial'};
regs = ["HPC", "V1"];
for r = 1:2
    nexttile(tl); hold on;
    for k = 1:size(terms, 1)
        for d = 1:numel(decs)
            row = C(C.decoder == decs(d) & C.region == regs(r) & C.model == terms{k, 1} & C.term == terms{k, 2}, :);
            if isempty(row), continue; end
            y = k + (d - 2.5) * 0.15;
            plot([row.ci_lo row.ci_hi], [y y], '-', 'Color', [0.5 0.5 0.5]);
            plot(row.estimate, y, dmk{d}, 'Color', 'k', 'MarkerFaceColor', [1 1 1] * (d > 2) * 0.7);
        end
    end
    % cross-region term
    other = regs(3 - r);
    k = size(terms, 1) + 1;
    for d = 1:numel(decs)
        row = C(C.decoder == decs(d) & C.region == regs(r) & C.model == "cross_from_" + other & C.term == other + "_grey_late", :);
        y = k + (d - 2.5) * 0.15;
        plot([row.ci_lo row.ci_hi], [y y], '-', 'Color', [0.5 0.5 0.5]);
        plot(row.estimate, y, dmk{d}, 'Color', 'k', 'MarkerFaceColor', [1 1 1] * (d > 2) * 0.7);
    end
    xline(0, 'k:');
    set(gca, 'YTick', 1:k, 'YTickLabel', [terms(:, 3); {sprintf('%s late grey (other region)', other)}], 'YDir', 'reverse');
    ylim([0.5 k + 0.5]);
    xlabel(sprintf('LME coefficient on %s onset bias (z per z)', regs(r)));
    title({sprintf('%s  %s onset bias predicted by...', char('A' + r - 1), regs(r)), ...
        'circle/square = Bayes lobo3/lobo5, triangles = PLS'}, 'FontWeight', 'normal');
end
% C: per-session partial Spearman (late grey vs onset), main decoder
nexttile(tl); hold on;
for r = 1:2
    for d = [1 3]
        ps = R.(decs(d)).(regs(r)).session_partial;
        xx = r + (d - 2) * 0.18;
        jitter_dots(xx, ps.rho, ps.animal, acol, 0.04);
        errorbar(xx, mean(ps.rho, 'omitnan'), sem(ps.rho), 'ko', 'MarkerFaceColor', 'k', 'CapSize', 0, 'LineWidth', 1.5);
        text(xx, -0.25, char(strrep(strrep(decs(d), 'bayes_', ''), '_', ' ')), 'FontSize', 7, 'HorizontalAlignment', 'center');
    end
end
yline(0, 'k:');
set(gca, 'XTick', 1:2, 'XTickLabel', regs); xlim([0.5 2.5]); ylim([-0.3 0.8]);
ylabel('Session partial Spearman, late grey vs onset 0-200 ms');
title({'C  Positive in nearly every session and all animals', '(covariates removed; dots = sessions)'}, 'FontWeight', 'normal');
ar_save_figure(fig, fullfile(out_dir, 'fig5_trial'));
close(fig);
end

%% ------------------------------------------------------------------------
function jitter_dots(x, v, an, acol, w)
if nargin < 5, w = 0.06; end
ok = isfinite(v);
scatter(x + w * (rand(sum(ok), 1) - 0.5) * 2, v(ok), 14, acol(an(ok), :), 'filled', 'MarkerFaceAlpha', 0.6);
[ua, ~, ia] = unique(an(ok));
am = accumarray(ia, v(ok), [], @mean);
scatter(x + 0.12 + zeros(size(am)), am, 30, acol(ua, :), 'd', 'filled', 'MarkerEdgeColor', 'k');
end

function s = sem(v)
v = v(isfinite(v));
s = std(v) / sqrt(numel(v));
end

function [rho, an] = session_rho(x, a, b)
us = unique(x.session);
rho = NaN(numel(us), 1); an = NaN(numel(us), 1);
for i = 1:numel(us)
    xi = x(x.session == us(i), :);
    if sum(isfinite(xi.(a)) & isfinite(xi.(b))) >= 8
        rho(i) = corr(xi.(a), xi.(b), 'Type', 'Spearman', 'Rows', 'complete');
    end
    an(i) = xi.animal(1);
end
end
