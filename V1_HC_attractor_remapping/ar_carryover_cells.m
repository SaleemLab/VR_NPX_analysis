function [Cc, Sx, Ab] = ar_carryover_cells(varargin)
% Which cells carry the trial-to-trial carry-over of track bias from the last
% 1 s of grey (late grey) into the first 200 ms of the track (onset)?
% Per session and region (HPC; V1, V1_L, V1_R without P.V1_exclude_subjects),
% spike counts per transition in onset-locked 100 ms bins (ar_bin_session):
% late grey -1..0 s, onset 0..0.2 s, onset_late 0.2..1 s, and a within-grey
% control (-1.2..-0.2 s -> -0.2..0 s). Counts are residualised per cell by the
% mean of the same block x previous x upcoming track group.
%  1. Linear stand-in for the decoder: LO_lin = sum_i w_i n_i, w_i = log rate
%     ratio T1/T2 (ar_cell_properties). Its pre -> post covariance splits exactly
%     into each cell's own persistence (diagonal: w_i^2 cov(n_i pre, n_i post))
%     and co-fluctuation between cells (off-diagonal; same- vs opposite-
%     preference pairs). Validated against the Bayesian carry-over.
%  2. Per cell: carry-over rho_i (Spearman, residual counts late grey vs onset;
%     also onset_late and the within-grey control) and share of the population
%     carry-over k_i = w_i * sum_j w_j cov(n_i pre, n_j post) / total.
%  3. Ablation with the actual decoder (multinomial Bayes, block-held-out lobo3,
%     track-balanced): carry-over rho without a cell group (start / middle / end
%     field, high-rate > 10 Hz) minus rho without the same number of random
%     cells (nRandom draws).
% Output: results\carryover_cells\ (cells.mat, csv, figures)
p = inputParser;
addParameter(p, 'nRandom', 20);
parse(p, varargin{:});
O = p.Results;
P = ar_setup();
S = load(fullfile(P.inputs, 'RUN1_decoding_input.mat'), 'D'); D = S.D;
Lt = load(fullfile(P.results, 'data', 'lap_tables.mat'), 'L'); L = Lt.L;
out_dir = fullfile(P.results, 'carryover_cells');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
bad = find(ismember(P.subjects, P.V1_exclude_subjects));
nPos = numel(P.pos_edges) - 1;
w = P.bin_width;
ep = struct('late', -10:-1, 'onset', 0:1, 'onset_late', 2:9, 'ctrl_pre', -12:-3, 'ctrl_post', -2:-1);
cell_rows = {}; sess_rows = {}; abl_rows = {};
rng(17);
for s = 1:numel(D)
    Ds = D(s); T = L(s).laps;
    trans = find(T.valid_transition); nTr = numel(trans);
    Cp = ar_cell_properties(Ds, T, P);
    Cp.theta_mvl = theta_mvl(Ds, P, s);
    G = ar_bin_session(Ds, T, P);
    % onset-locked bin index per transition and relative bin
    rels = -15:14;
    idx = NaN(nTr, numel(rels));
    r2 = find(G.kind == 2);
    [~, it] = ismember(G.trans(r2), trans); [~, ir] = ismember(G.rel(r2), rels);
    ok = it > 0 & ir > 0;
    idx(sub2ind(size(idx), it(ok), ir(ok))) = r2(ok);
    grp = findgroups(T.block(trans), T.prev_track(trans), T.track(trans));
    % residual counts per epoch (transitions x cells)
    names = fieldnames(ep);
    R = struct();
    for e = 1:numel(names)
        cols = ismember(rels, ep.(names{e}));
        X = zeros(nTr, numel(Cp.cell_id));
        okt = all(~isnan(idx(:, cols)), 2);
        for i = find(okt)'
            X(i, :) = sum(G.counts(idx(i, cols), :), 1);
        end
        X(~okt, :) = NaN;
        R.(names{e}) = residualise(X, grp);
    end
    [occ, cnt] = ar_bayes_lap_stats(Ds, T, P);
    for reg = {'HPC', 'V1', 'V1_L', 'V1_R'}
        if startsWith(reg{1}, 'V1') && ismember(L(s).animal_index, bad), continue; end
        switch reg{1}
            case 'HPC', cells = contains(Cp.region, 'HPC');
            case 'V1', cells = contains(Cp.region, 'V1');
            otherwise, cells = Cp.region == reg{1};
        end
        cells = find(cells);
        if numel(cells) < P.min_cells, continue; end
        wv = Cp.w(cells);
        % 1. linear stand-in and covariance decomposition
        for pr = {{'late', 'onset', 'across'}, {'late', 'onset_late', 'across_late'}, {'ctrl_pre', 'ctrl_post', 'within_grey'}}
            A = R.(pr{1}{1})(:, cells); B = R.(pr{1}{2})(:, cells);
            okr = all(isfinite(A), 2) & all(isfinite(B), 2);
            A = A(okr, :); B = B(okr, :);
            Cab = (A - mean(A))' * (B - mean(B)) / (size(A, 1) - 1);
            total = wv' * Cab * wv;
            diagc = sum(wv.^2 .* diag(Cab));
            WW = wv * wv' .* Cab; WW(logical(eye(numel(wv)))) = 0;
            same = sum(WW(sign(wv) == sign(wv')), 'all'); opp = sum(WW(sign(wv) ~= sign(wv')), 'all');
            rho_lin = corr(A * wv, B * wv, 'Type', 'Spearman');
            % components in correlation units: they add up to the Pearson r of the
            % linear readouts (bounded, unlike fractions of a small total)
            den = sqrt((wv' * cov(A) * wv) * (wv' * cov(B) * wv));
            sess_rows(end + 1, :) = {s, L(s).animal_index, string(reg{1}), string(pr{1}{3}), numel(cells), sum(okr), ...
                rho_lin, total, diagc / total, (same + opp) / total, same / total, opp / total, ...
                total / den, diagc / den, same / den, opp / den}; %#ok<AGROW>
            if strcmp(pr{1}{3}, 'across')
                k = wv .* (Cab * wv) / total;            % share of the population carry-over
                kd = wv.^2 .* diag(Cab) / total;
                rho_i = NaN(numel(cells), 1); rho_l = rho_i; rho_c = rho_i;
                Al = R.onset_late(okr, cells); Ac = R.ctrl_pre(:, cells); Bc = R.ctrl_post(:, cells);
                for j = 1:numel(cells)
                    if std(A(:, j)) > 0 && std(B(:, j)) > 0, rho_i(j) = corr(A(:, j), B(:, j), 'Type', 'Spearman'); end
                    if std(Al(:, j), 'omitnan') > 0, rho_l(j) = corr(A(:, j), Al(:, j), 'Type', 'Spearman', 'Rows', 'complete'); end
                    oc = isfinite(Ac(:, j)) & isfinite(Bc(:, j));
                    if std(Ac(oc, j)) > 0 && std(Bc(oc, j)) > 0, rho_c(j) = corr(Ac(oc, j), Bc(oc, j), 'Type', 'Spearman'); end
                end
                Cx = Cp(cells, :);
                Cx.session = repmat(s, numel(cells), 1); Cx.animal = repmat(L(s).animal_index, numel(cells), 1);
                Cx.group_region = repmat(string(reg{1}), numel(cells), 1);
                Cx.rho_onset = rho_i; Cx.rho_onset_late = rho_l; Cx.rho_within_grey = rho_c;
                Cx.share = k; Cx.share_diag = kd;
                cell_rows{end + 1} = Cx; %#ok<AGROW>
            end
        end
        % 3. ablation with the decoder
        groups = {'start field', Cp.field_group(cells) == "start"; 'middle field', Cp.field_group(cells) == "middle"; ...
            'end field', Cp.field_group(cells) == "end"; 'high rate', Cp.high_rate(cells)};
        bins_late = idx(:, ismember(rels, ep.late)); bins_on = idx(:, ismember(rels, ep.onset));
        full = decoded_carryover(cells, true(numel(cells), 1), G, T, trans, occ, cnt, s, P, nPos, bins_late, bins_on, grp);
        for g = 1:size(groups, 1)
            drop = groups{g, 2};
            if sum(drop) < 2 || sum(~drop) < P.min_cells, continue; end
            rab = decoded_carryover(cells, ~drop, G, T, trans, occ, cnt, s, P, nPos, bins_late, bins_on, grp);
            rnd = NaN(O.nRandom, 1);
            for q = 1:O.nRandom
                keep = true(numel(cells), 1); keep(randperm(numel(cells), sum(drop))) = false;
                rnd(q) = decoded_carryover(cells, keep, G, T, trans, occ, cnt, s, P, nPos, bins_late, bins_on, grp);
            end
            abl_rows(end + 1, :) = {s, L(s).animal_index, string(reg{1}), string(groups{g, 1}), sum(drop), numel(cells), ...
                full, rab, mean(rnd), rab - mean(rnd)}; %#ok<AGROW>
        end
    end
    fprintf('Session %d done\n', s);
end
Cc = vertcat(cell_rows{:});
Sx = cell2table(sess_rows, 'VariableNames', {'session', 'animal', 'region', 'pair', 'n_cells', 'n_trans', 'rho_lin', ...
    'total_cov', 'frac_diag', 'frac_offdiag', 'frac_same_pref_pairs', 'frac_opposite_pref_pairs', ...
    'r_total', 'r_own', 'r_same_pairs', 'r_opp_pairs'});
Ab = cell2table(abl_rows, 'VariableNames', {'session', 'animal', 'region', 'group', 'n_dropped', 'n_cells', 'rho_full', ...
    'rho_without_group', 'rho_without_random', 'effect'});
save(fullfile(out_dir, 'carryover_cells.mat'), 'Cc', 'Sx', 'Ab', '-v7.3');
writetable(Sx, fullfile(out_dir, 'carryover_sessions.csv'));
writetable(Ab, fullfile(out_dir, 'carryover_ablation.csv'));
summarise(Cc, Sx, Ab, P, out_dir);
end

%% ------------------------------------------------------------------------
function Xr = residualise(X, grp)
Xr = X;
for k = unique(grp)'
    r = grp == k;
    Xr(r, :) = X(r, :) - mean(X(r, :), 1, 'omitnan');
end
end

function mvl = theta_mvl(Ds, P, s)
% theta-phase mean vector length per cell while running > 5 cm/s (ar_theta_phase)
f = fullfile(P.results, 'theta_phase', sprintf('session_%02d.mat', s));
mvl = NaN(numel(Ds.cell_id), 1);
if ~exist(f, 'file'), return; end
F = load(f, 't', 'phase');
spd = interp1(Ds.tvec, Ds.speed, F.t, 'linear', 0);
[~, ci] = ismember(Ds.spike_id, Ds.cell_id);
b = discretize(Ds.spike_times, [F.t - 0.005; F.t(end) + 0.005]);
ok = ci > 0 & ~isnan(b);
ok(ok) = spd(b(ok)) > 5;
ph = double(F.phase(b(ok)));
z = accumarray(ci(ok), exp(1i * ph), [numel(Ds.cell_id) 1]);
n = accumarray(ci(ok), 1, [numel(Ds.cell_id) 1]);
mvl = abs(z) ./ n; mvl(n < 50) = NaN;
end

function rho = decoded_carryover(cells, keep, G, T, trans, occ, cnt, s, P, nPos, bins_late, bins_on, grp)
% carry-over of the multinomial Bayes log odds (late grey mean -> onset mean,
% residualised by block x track pair), decoding with the kept cells only
cl = cells(keep);
nTr = numel(trans);
lo_late = NaN(nTr, 1); lo_on = NaN(nTr, 1);
nl = size(bins_late, 2); no = size(bins_on, 2);
for b = unique(T.block(trans))'
    tr = find(T.block(trans) == b & all(~isnan([bins_late bins_on]), 2));
    if isempty(tr), continue; end
    train = ismember((1:height(T))', ar_balanced_train_laps(T, b, 1, 1000 * s + 10 * b + 1));
    maps = ar_rate_maps(occ, cnt(cl, :, :), train, T.track, nPos);
    bins = [bins_late(tr, :) bins_on(tr, :)]';           % (nl + no) x numel(tr)
    lo = ar_bayes_track(double(G.counts(bins(:), cl))', maps, P.bin_width, 'multinomial', P.rate_floor);
    lo = reshape(lo, nl + no, numel(tr));
    lo_late(tr) = mean(lo(1:nl, :), 1); lo_on(tr) = mean(lo(nl + 1:end, :), 1);
end
x = residualise(lo_late, grp); y = residualise(lo_on, grp);
rho = corr(x, y, 'Type', 'Spearman', 'Rows', 'complete');
end

%% ------------------------------------------------------------------------
function summarise(Cc, Sx, Ab, P, out_dir)
acol = lines(6);
regs = ["HPC", "V1", "V1_L", "V1_R"];
stats = {};
fig = figure('Position', [50 50 1900 800], 'Color', 'w');
tl = tiledlayout(2, 4, 'TileSpacing', 'compact');
% A: decomposition of the linear carry-over r into own-cell persistence and
% co-fluctuation of same- / opposite-preference pairs (correlation units; the
% three add up to r_total), across onset (filled) and within grey (open)
nexttile(tl); hold on;
parts = {'r_total', 'r_own', 'r_same_pairs', 'r_opp_pairs'}; pl = {'total', 'own cell', 'same-pref pairs', 'opp-pref pairs'};
for r = 1:2
    for pairI = 1:2
        pn = ["across", "within_grey"]; pn = pn(pairI);
        x = Sx(Sx.region == regs(r) & Sx.pair == pn, :);
        for k = 1:numel(parts)
            xx = (r - 1) * 5 + k + (pairI - 1.5) * 0.3; v = x.(parts{k});
            if pairI == 1, mf = 'k'; else, mf = 'w'; end
            errorbar(xx, mean(v), std(v) / sqrt(numel(v)), 'ko', 'MarkerFaceColor', mf, 'CapSize', 0, 'LineWidth', 1.2);
            stats(end + 1, :) = {regs(r), "r_decomposition_" + pn, string(pl{k}), mean(v), signrank(v), sum(v > 0), numel(v)}; %#ok<AGROW>
        end
    end
end
yline(0, 'k:');
set(gca, 'XTick', [1:4 6:9], 'XTickLabel', ["HPC " + string(pl), "V1 " + string(pl)]); xtickangle(35);
ylabel('Contribution to carry-over r');
title({'A  Own-cell persistence vs co-fluctuation', 'filled = across onset, open = within grey'}, 'FontWeight', 'normal');
% B: linear stand-in vs Bayesian carry-over
nexttile(tl); hold on;
S = load(fullfile(P.results, 'trial_prediction', 'trial_prediction.mat'), 'R'); Rt = S.R;
for r = 1:2
    ps = Rt.bayes_multinomial.(char(regs(r))).session_partial;
    x = Sx(Sx.region == regs(r) & Sx.pair == "across", :);
    [~, ia, ib] = intersect(x.session, ps.session);
    mk = {'o', 's'};
    scatter(ps.rho(ib), x.rho_lin(ia), 30, acol(x.animal(ia), :), mk{r}, 'filled');
    rr = corr(ps.rho(ib), x.rho_lin(ia), 'Rows', 'complete');
    stats(end + 1, :) = {regs(r), "lin_vs_bayes_session_r", "", rr, NaN, NaN, numel(ia)}; %#ok<AGROW>
end
plot([0 0.8], [0 0.8], 'k:');
xlabel('Bayes carry-over (partial rho, section 5)'); ylabel('Linear stand-in carry-over (rho)');
title({'B  The linear stand-in tracks the decoder', 'circles HPC, squares V1'}, 'FontWeight', 'normal');
% C, D: single-cell carry-over by field group and by rate (HPC, V1)
for r = 1:2
    nexttile(tl); hold on;
    x = Cc(Cc.group_region == regs(r), :);
    gl = {'start', 'middle', 'end', 'low rate', 'high rate'};
    sel = {x.field_group == "start" & ~x.high_rate, x.field_group == "middle" & ~x.high_rate, x.field_group == "end" & ~x.high_rate, ...
        ~x.high_rate, x.high_rate};
    for g = 1:numel(gl)
        for mI = 1:2
            v = ifelse(mI == 1, x.rho_onset, x.rho_within_grey);
            sm = session_means(v(sel{g}), x.session(sel{g}));
            col = ifelse(mI == 1, [0 0 0], [0.6 0.6 0.6]);
            errorbar(g + (mI - 1.5) * 0.3, mean(sm, 'omitnan'), std(sm, 'omitnan') / sqrt(sum(isfinite(sm))), 'o', ...
                'Color', col, 'MarkerFaceColor', col, 'CapSize', 0, 'LineWidth', 1.5);
            if mI == 1
                stats(end + 1, :) = {regs(r), "cell_rho_onset", string(gl{g}), mean(sm, 'omitnan'), signrank(sm), sum(sm > 0), sum(isfinite(sm))}; %#ok<AGROW>
            end
        end
    end
    yline(0, 'k:');
    set(gca, 'XTick', 1:numel(gl), 'XTickLabel', gl); xlim([0.5 numel(gl) + 0.5]);
    ylabel('Cell carry-over rho (session mean)');
    title({sprintf('%s  %s single-cell carry-over', char('C' + r - 1), regs(r)), ...
        'black across onset, grey within grey'}, 'FontWeight', 'normal');
end
% E: share of the population carry-over vs share of cells, by group
nexttile(tl); hold on;
gnames = ["start", "middle", "end", "high rate"];
for r = 1:2
    x = Cc(Cc.group_region == regs(r), :);
    us = unique(x.session);
    for g = 1:numel(gnames)
        sh = NaN(numel(us), 1); nf = sh;
        for i = 1:numel(us)
            xs = x(x.session == us(i), :);
            if g < 4, m = xs.field_group == gnames(g); else, m = xs.high_rate; end
            sh(i) = sum(xs.share(m)); nf(i) = mean(m);
        end
        xx = (r - 1) * 5 + g;
        bar(xx - 0.2, mean(nf), 0.35, 'FaceColor', [0.8 0.8 0.8]);
        bar(xx + 0.2, mean(sh, 'omitnan'), 0.35, 'FaceColor', [0.2 0.45 0.75]);
        errorbar(xx + 0.2, mean(sh, 'omitnan'), std(sh, 'omitnan') / sqrt(numel(us)), 'k', 'CapSize', 0);
        stats(end + 1, :) = {regs(r), "share_minus_cell_fraction", gnames(g), mean(sh - nf, 'omitnan'), signrank(sh - nf), sum(sh - nf > 0), numel(us)}; %#ok<AGROW>
    end
end
yline(0, 'k:');
set(gca, 'XTick', [1:4 6:9], 'XTickLabel', ["HPC " + gnames, "V1 " + gnames]); xtickangle(35);
ylabel('Fraction'); title({'E  Share of carry-over (blue) vs share of cells (grey)', 'linear stand-in'}, 'FontWeight', 'normal');
% F: ablation with the decoder
nexttile(tl); hold on;
gl = ["start field", "middle field", "end field", "high rate"];
for r = 1:2
    x = Ab(Ab.region == regs(r), :);
    for g = 1:numel(gl)
        v = x.effect(x.group == gl(g)); an = x.animal(x.group == gl(g));
        if isempty(v), continue; end
        xx = (r - 1) * 5 + g;
        scatter(xx + 0.1 * randn(numel(v), 1), v, 14, acol(an, :), 'filled', 'MarkerFaceAlpha', 0.6);
        errorbar(xx + 0.3, mean(v), std(v) / sqrt(numel(v)), 'ko', 'MarkerFaceColor', 'k', 'CapSize', 0, 'LineWidth', 1.5);
        [ua, ~, ia] = unique(an); am = accumarray(ia, v, [], @mean);
        stats(end + 1, :) = {regs(r), "ablation_effect", gl(g), mean(v), signrank(v), sum(v > 0), numel(v)}; %#ok<AGROW>
        text(xx, max(v) + 0.03, sprintf('%.2g', signrank(v)), 'HorizontalAlignment', 'center', 'FontSize', 7);
        clear ua am
    end
end
yline(0, 'k:');
set(gca, 'XTick', [1:4 6:9], 'XTickLabel', ["HPC -" + gl, "V1 -" + gl]); xtickangle(35);
ylabel('Carry-over rho: without group - without random cells');
title({'F  Ablation with the Bayes decoder', '< 0: the group carries more than its share'}, 'FontWeight', 'normal');
% G: V1 hemispheres, single-cell carry-over and decomposition
nexttile(tl); hold on;
for r = 3:4
    x = Cc(Cc.group_region == regs(r), :);
    sm = session_means(x.rho_onset, x.session);
    y = Sx(Sx.region == regs(r) & Sx.pair == "across", :);
    errorbar(r - 2 - 0.15, mean(sm, 'omitnan'), std(sm, 'omitnan') / sqrt(sum(isfinite(sm))), 'ko', 'MarkerFaceColor', 'k', 'CapSize', 0, 'LineWidth', 1.5);
    errorbar(r - 2 + 0.15, mean(y.rho_lin), std(y.rho_lin) / sqrt(height(y)), 'o', 'Color', [0.2 0.45 0.75], 'MarkerFaceColor', [0.2 0.45 0.75], 'CapSize', 0, 'LineWidth', 1.5);
    stats(end + 1, :) = {regs(r), "cell_rho_onset_all", "", mean(sm, 'omitnan'), signrank(sm), sum(sm > 0), sum(isfinite(sm))}; %#ok<AGROW>
    stats(end + 1, :) = {regs(r), "population_rho_lin", "", mean(y.rho_lin), signrank(y.rho_lin), sum(y.rho_lin > 0), height(y)}; %#ok<AGROW>
end
yline(0, 'k:'); set(gca, 'XTick', 1:2, 'XTickLabel', {'V1 L', 'V1 R'}); xlim([0.5 2.5]);
ylabel('rho'); title({'G  V1 hemispheres', 'black = single cells, blue = population (linear)'}, 'FontWeight', 'normal');
% H: LME of cell carry-over on properties (HPC)
nexttile(tl); hold on;
x = Cc(Cc.group_region == "HPC" & isfinite(Cc.rho_onset), :);
x.abs_pref = abs(x.pref); x.log_rate = log(x.session_rate);
x.field_group = categorical(x.field_group, ["middle", "start", "end"]);
x.animal_c = categorical(x.animal); x.sess_c = categorical(x.session);
lme = fitlme(x, 'rho_onset ~ field_group + log_rate + abs_pref + theta_mvl + (1|animal_c) + (1|sess_c)');
c = lme.Coefficients; ci = coefCI(lme);
terms = c.Name(2:end);
for k = 1:numel(terms)
    plot(ci(k + 1, :), [k k], '-', 'Color', [0.5 0.5 0.5]);
    plot(c.Estimate(k + 1), k, 'ko', 'MarkerFaceColor', 'k');
    stats(end + 1, :) = {"HPC", "cell_rho_lme", string(terms{k}), c.Estimate(k + 1), c.pValue(k + 1), NaN, lme.NumObservations}; %#ok<AGROW>
end
xline(0, 'k:');
set(gca, 'YTick', 1:numel(terms), 'YTickLabel', strrep(terms, '_', ' '), 'YDir', 'reverse'); ylim([0.5 numel(terms) + 0.5]);
xlabel('Effect on cell carry-over rho'); title({'H  HPC cell carry-over vs properties (LME)', 'field vs middle; log rate; |pref|; theta MVL'}, 'FontWeight', 'normal');
title(tl, 'Single-cell decomposition of the late grey -> onset (0-200 ms) carry-over');
ar_save_figure(fig, fullfile(out_dir, 'carryover_cells'));
close(fig);
T = cell2table(stats, 'VariableNames', {'region', 'measure', 'term', 'estimate', 'p', 'n_pos', 'n'});
writetable(T, fullfile(out_dir, 'carryover_stats.csv'));
end

function sm = session_means(v, s)
us = unique(s); sm = NaN(numel(us), 1);
for i = 1:numel(us), sm(i) = mean(v(s == us(i)), 'omitnan'); end
end

function v = ifelse(c, a, b)
if c, v = a; else, v = b; end
end
