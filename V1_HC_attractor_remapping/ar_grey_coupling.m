function R = ar_grey_coupling(varargin)
% Do HPC and V1 track biases co-fluctuate in the grey screen (early: 0.3-1.3 s
% after lap end, excluding V1's offset response; late: last 1 s before onset),
% and does that coupling shape how well the late-grey bias predicts the onset
% bias? Lags are computed within each window only (see ar_grey_lag for lead /
% lag over the whole grey).
% Per session (V1 animals only), event-locked 100 ms zLO bins:
%   residual = zLO minus the mean time course of transitions with the same
%   previous x upcoming track pair (removes shared stereotyped dynamics such as
%   the decay after lap end); 'ctrl' also regresses out each region's log spike
%   count and speed per bin.
%   between-trial r : trial means (block x pair means removed) across trials
%   within-trial r  : bins after removing each trial's mean, concatenated over
%                     trials; lagged versions (lag > 0: V1 later than HPC); null
%                     by pairing HPC and V1 from different trials (same pair)
%   concatenated r  : residual bins without removing trial means (mixes both)
%   per-trial r     : correlation of the 10 trial-demeaned bins (Fisher z)
% Trial models (ar_trial_lme): onset ~ grey_late * coupling (interaction = does
% the pre-onset bias predict the onset bias better on coupled trials?), for HPC
% and V1 onset (0-200, 0-100, 100-200 ms); factors of coupling (LME).
% Output: results\grey_coupling\ (csv, mat, figure)
p = inputParser;
addParameter(p, 'decoders', {'bayes_multinomial', 'pls'});
addParameter(p, 'nShuffle', 200);
parse(p, varargin{:});
O = p.Results;
P = ar_setup();
S = load(fullfile(P.results, 'data', 'epoch_table.mat'), 'E', 'TC'); E = S.E; TC = S.TC;
out_dir = fullfile(P.results, 'grey_coupling');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
bad = find(ismember(P.subjects, P.V1_exclude_subjects));
lags = -5:5;
% early grey starts 0.3 s after lap end: the first ~200 ms carry V1's latency /
% offset response to the screen turning grey (a visual event, not grey activity)
win = struct('late', struct('field', 'on', 'cols', find(TC.t_on > -1 & TC.t_on < 0)), ...
    'early', struct('field', 'off', 'cols', find(TC.t_off > 0.3 & TC.t_off < 1.3)));
wn = fieldnames(win);
rows = {}; srows = {};
R = struct();
rng(11);
for d = 1:numel(O.decoders)
    dec = O.decoders{d};
    trial_c = table();
    for w = 1:numel(wn)
        Wd = win.(wn{w});
        A = TC.(Wd.field);
        sH = find(strcmp({A.region}, 'HPC') & strcmp({A.decoder}, dec));
        for i = 1:numel(sH)
            s = A(sH(i)).session;
            Ds = TC.trans(TC.trans.session == s, :);
            if ismember(Ds.animal(1), bad), continue; end
            iV = find(strcmp({A.region}, 'V1') & strcmp({A.decoder}, dec) & [A.session] == s);
            if isempty(iV), continue; end
            for ctrl = [false true]
                H = residual(A(sH(i)), Wd.cols, Ds, ctrl);
                V = residual(A(iV), Wd.cols, Ds, ctrl);
                st = coupling_stats(H, V, Ds, lags, O.nShuffle);
                srows(end + 1, :) = {string(dec), string(wn{w}), ctrl, s, Ds.animal(1), st.between, st.within, ...
                    st.concat, st.null_mean, st.null_sd, st.lagged}; %#ok<AGROW>
                if ~ctrl
                    tc = table(repmat(s, height(Ds), 1), Ds.lap, st.per_trial, 'VariableNames', {'session', 'lap', ['coupling_' wn{w}]});
                    if w == 1, trial_c = [trial_c; tc]; %#ok<AGROW>
                    else
                        [~, ia, ib] = intersect(trial_c(:, {'session', 'lap'}), tc(:, {'session', 'lap'}));
                        if ~ismember(['coupling_' wn{w}], trial_c.Properties.VariableNames), trial_c.(['coupling_' wn{w}]) = NaN(height(trial_c), 1); end
                        trial_c.(['coupling_' wn{w}])(ia) = tc.(['coupling_' wn{w}])(ib);
                    end
                end
            end
        end
    end
    R.(dec).trial_coupling = trial_c;
    % trial models
    for reg = {'HPC', 'V1'}
        x = E(E.decoder == dec & E.region == reg{1} & ~ismember(E.animal, bad), :);
        T = ar_trial_table(x);
        [~, ia, ib] = intersect(T(:, {'session', 'lap'}), trial_c(:, {'session', 'lap'}));
        T.coupling_late = NaN(height(T), 1); T.coupling_early = NaN(height(T), 1);
        T.coupling_late(ia) = trial_c.coupling_late(ib); T.coupling_early(ia) = trial_c.coupling_early(ib);
        for resp = {'onset', 'onset_100', 'onset_200'}
            [~, st] = ar_trial_lme(T, 'grey_late', {'coupling_late', 'grey_late:coupling_late'}, resp{1});
            rows(end + 1, :) = coef_row(dec, reg{1}, ['interaction_' resp{1}], 'grey_late:coupling_late', ar_lme_coef(st.lme, 'grey_late:coupling_late')); %#ok<AGROW>
            rows(end + 1, :) = coef_row(dec, reg{1}, ['interaction_' resp{1}], 'grey_late', st); %#ok<AGROW>
        end
        % factors of late-grey coupling (from the HPC table; covariates are
        % per region, so V1 spikes are added from the V1 table)
        if strcmp(reg{1}, 'HPC'), TH = T; else, TV = T; end
    end
    [~, ia, ib] = intersect(TH(:, {'session', 'lap'}), TV(:, {'session', 'lap'}));
    F = TH(ia, :);
    F.spk_grey_late_V1 = TV.spk_grey_late(ib);
    F.abs_bias_H = abs(TH.grey_late(ia)); F.abs_bias_V = abs(TV.grey_late(ib));
    F.agree = sign(TH.grey_late(ia)) == sign(TV.grey_late(ib));
    [g, ~] = findgroups(F.session);
    mu = splitapply(@mean, F.lap, g); sd = splitapply(@std, F.lap, g);
F.lap_z = (F.lap - mu(g)) ./ sd(g);
    vars = {'spd_grey_late', 'spk_grey_late', 'spk_grey_late_V1', 'lap_in_block', 'is_switch', 'abs_bias_H', ...
        'abs_bias_V', 'lick_grey_late', 'grey_dur', 'lap_z'};
    F = F(all(isfinite([F.coupling_late F{:, vars}]), 2), :);
    lme = fitlme(F, ['coupling_late ~ ' strjoin(vars, ' + ') ' + (1|animal_c) + (1|sess_c)']);
    for v = vars
        nm = v{1}; if strcmp(nm, 'is_switch'), nm = 'is_switch_1'; end
        rows(end + 1, :) = coef_row(dec, 'HPC_V1', 'coupling_factors', v{1}, ar_lme_coef(lme, nm)); %#ok<AGROW>
    end
    % does mean-level agreement (both biased to the same track) matter?
    TH2 = TH(ia, :); TH2.grey_late_V1 = TV.grey_late(ib);
    [~, st] = ar_trial_lme(TH2, 'grey_late', {'grey_late_V1', 'grey_late:grey_late_V1'});
    rows(end + 1, :) = coef_row(dec, 'HPC', 'agreement', 'grey_late:grey_late_V1', ar_lme_coef(st.lme, 'grey_late:grey_late_V1')); %#ok<AGROW>
end
C = cell2table(rows, 'VariableNames', {'decoder', 'region', 'model', 'term', 'estimate', 'ci_lo', 'ci_hi', 'p', 'n'});
writetable(C, fullfile(out_dir, 'grey_coupling_models.csv'));
Sx = cell2table(srows, 'VariableNames', {'decoder', 'window', 'ctrl', 'session', 'animal', 'r_between', 'r_within', ...
    'r_concat', 'null_mean', 'null_sd', 'r_lagged'});
save(fullfile(out_dir, 'grey_coupling.mat'), 'R', 'C', 'Sx', 'lags', '-v7.3');
writetable(removevars(Sx, 'r_lagged'), fullfile(out_dir, 'grey_coupling_sessions.csv'));
summarise(Sx, C, lags, P, out_dir);
end

%% ------------------------------------------------------------------------
function Z = residual(A, cols, Ds, ctrl)
% bins of the window minus the mean time course of the same track pair;
% optionally regress out log spike count and speed of each bin (session level)
Z = double(A.z(:, cols));
g = findgroups(Ds.prev_track, Ds.track);
for k = unique(g)'
    Z(g == k, :) = Z(g == k, :) - mean(Z(g == k, :), 1, 'omitnan');
end
if ctrl
    X = [ones(numel(Z), 1) log1p(double(reshape(A.spikes(:, cols), [], 1))) double(reshape(A.speed(:, cols), [], 1))];
    y = Z(:);
    ok = all(isfinite(X), 2) & isfinite(y);
    b = X(ok, :) \ y(ok);
    y(ok) = y(ok) - X(ok, :) * b;
    Z = reshape(y, size(Z));
end
end

function st = coupling_stats(H, V, Ds, lags, nShuffle)
% between-trial, within-trial (lagged), concatenated, per-trial coupling
% between-trial: trial means with block x pair means removed
mh = mean(H, 2, 'omitnan'); mv = mean(V, 2, 'omitnan');
g = findgroups(Ds.block, Ds.prev_track, Ds.track);
gm = splitapply(@(x) mean(x, 'omitnan'), mh, g); mh = mh - gm(g);
gm = splitapply(@(x) mean(x, 'omitnan'), mv, g); mv = mv - gm(g);
st.between = corr(mh, mv, 'Rows', 'complete');
% concatenated (residual bins, trial means kept)
st.concat = corr(H(:), V(:), 'Rows', 'complete');
% within-trial
Hc = H - mean(H, 2, 'omitnan'); Vc = V - mean(V, 2, 'omitnan');
st.within = corr(Hc(:), Vc(:), 'Rows', 'complete');
nb = size(H, 2);
st.lagged = NaN(1, numel(lags));
for k = 1:numel(lags)
    L = lags(k);
    if L >= 0, a = Hc(:, 1:nb - L); b = Vc(:, 1 + L:nb);
    else, a = Hc(:, 1 - L:nb); b = Vc(:, 1:nb + L); end
    st.lagged(k) = corr(a(:), b(:), 'Rows', 'complete');
end
% null: V1 from another trial with the same track pair
gp = findgroups(Ds.prev_track, Ds.track);
nul = NaN(nShuffle, 1);
for sh = 1:nShuffle
    perm = (1:size(V, 1))';
    for k = unique(gp)'
        ix = find(gp == k);
        perm(ix) = ix(randperm(numel(ix)));
    end
    Vs = Vc(perm, :);
    nul(sh) = corr(Hc(:), Vs(:), 'Rows', 'complete');
end
st.null_mean = mean(nul); st.null_sd = std(nul);
% per trial
st.per_trial = NaN(size(H, 1), 1);
for i = 1:size(H, 1)
    ok = isfinite(Hc(i, :)) & isfinite(Vc(i, :));
    if sum(ok) >= 6 && std(Hc(i, ok)) > 0 && std(Vc(i, ok)) > 0
        st.per_trial(i) = atanh(min(max(corr(Hc(i, ok)', Vc(i, ok)'), -0.99), 0.99));
    end
end
end

function row = coef_row(dec, reg, model, term, st)
row = {string(dec), string(reg), string(model), string(term), st.est, st.ci(1), st.ci(2), st.p, st.n};
end

%% ------------------------------------------------------------------------
function summarise(Sx, C, lags, P, out_dir)
acol = lines(6);
fig = figure('Position', [50 50 1700 430], 'Color', 'w');
tl = tiledlayout(1, 4, 'TileSpacing', 'loose');
x = Sx(Sx.decoder == "bayes_multinomial" & ~Sx.ctrl, :);
% A: lagged within-trial correlation
nexttile(tl); hold on;
wc = {'late', [0 0 0]; 'early', [0.6 0.6 0.6]};
for w = 1:2
    xw = x(x.window == wc{w, 1}, :);
    L = xw.r_lagged;
    if iscell(L), L = vertcat(L{:}); end
    m = mean(L, 1, 'omitnan'); se = std(L, 0, 1, 'omitnan') / sqrt(size(L, 1));
    errorbar(lags * 100, m, se, 'o-', 'Color', wc{w, 2}, 'MarkerFaceColor', wc{w, 2}, 'CapSize', 0, 'LineWidth', 1.5);
    yline(mean(xw.null_mean), ':', 'Color', wc{w, 2});
end
xline(0, 'k:');
xlabel('Lag (ms; > 0: V1 after HPC)'); ylabel('Within-trial HPC-V1 correlation');
legend({'late grey', 'null', 'early grey', 'null'}, 'Box', 'off');
title({'A  Bin-by-bin co-fluctuation within trials', 'multinomial Bayes, mean +- SEM over sessions'}, 'FontWeight', 'normal');
% B: between / within / concatenated, per session
nexttile(tl); hold on;
meas = {'r_between', 'r_within', 'r_concat'}; ml = {'between', 'within', 'concat.'};
for w = 1:2
    xw = x(x.window == wc{w, 1}, :);
    for k = 1:3
        xx = (w - 1) * 4 + k;
        v = xw.(meas{k});
        scatter(xx + 0.08 * randn(height(xw), 1), v, 18, acol(xw.animal, :), 'filled', 'MarkerFaceAlpha', 0.7);
        errorbar(xx + 0.3, mean(v, 'omitnan'), std(v, 'omitnan') / sqrt(sum(isfinite(v))), 'ko', 'MarkerFaceColor', 'k', 'CapSize', 0, 'LineWidth', 1.5);
    end
end
yline(0, 'k:');
set(gca, 'XTick', [1:3 5:7], 'XTickLabel', ["late " + string(ml), "early " + string(ml)]); xtickangle(30);
ylabel('HPC-V1 correlation');
title({'B  Between trials, within trials, concatenated', 'dots = sessions (colour = animal)'}, 'FontWeight', 'normal');
% C: interaction coefficients
nexttile(tl); hold on;
mods = {'interaction_onset', 'interaction_onset_100', 'interaction_onset_200'};
mlab = {'0-200 ms', '0-100 ms', '100-200 ms'};
regs = {'HPC', 'V1'}; decs = {'bayes_multinomial', 'pls'}; mk = {'o', '^'};
yl = {};
yy = 0;
for r = 1:2
    for m = 1:3
        yy = yy + 1; yl{end + 1} = sprintf('%s onset %s', regs{r}, mlab{m}); %#ok<AGROW>
        for d = 1:2
            row = C(C.decoder == decs{d} & C.region == regs{r} & C.model == mods{m} & C.term == "grey_late:coupling_late", :);
            if isempty(row), continue; end
            y = yy + (d - 1.5) * 0.25;
            plot([row.ci_lo row.ci_hi], [y y], '-', 'Color', [0.5 0.5 0.5]);
            plot(row.estimate, y, mk{d}, 'Color', 'k', 'MarkerFaceColor', [1 1 1] * (d - 1) * 0.7);
        end
    end
end
xline(0, 'k:');
set(gca, 'YTick', 1:yy, 'YTickLabel', yl, 'YDir', 'reverse'); ylim([0.5 yy + 0.5]);
xlabel('late grey x coupling (LME)');
title({'C  Does coupling boost the pre -> post prediction?', 'circle = Bayes, triangle = PLS'}, 'FontWeight', 'normal');
% D: factors of late-grey coupling
nexttile(tl); hold on;
f = C(C.model == "coupling_factors", :);
terms = unique(f.term, 'stable');
for k = 1:numel(terms)
    for d = 1:2
        row = f(f.term == terms{k} & f.decoder == decs{d}, :);
        y = k + (d - 1.5) * 0.25;
        plot([row.ci_lo row.ci_hi], [y y], '-', 'Color', [0.5 0.5 0.5]);
        plot(row.estimate, y, mk{d}, 'Color', 'k', 'MarkerFaceColor', [1 1 1] * (d - 1) * 0.7);
    end
end
xline(0, 'k:');
set(gca, 'YTick', 1:numel(terms), 'YTickLabel', strrep(terms, '_', ' '), 'YDir', 'reverse'); ylim([0.5 numel(terms) + 0.5]);
xlabel('Effect on per-trial late-grey coupling (Fisher z)');
title({'D  What goes with stronger coupling?', 'LME, animal and session random intercepts'}, 'FontWeight', 'normal');
ar_save_figure(fig, fullfile(out_dir, 'grey_coupling'));
close(fig);
end
