function R = ar_trial_prediction(varargin)
% Does the track bias in the last 1 s of grey (grey_late) predict the bias at
% track onset, trial by trial, beyond other factors? Responses: onset (0-200 ms,
% main), onset_100 (0-100 ms), onset_200 (100-200 ms), onset_late (0.2-1 s).
% Signed log odds toward T1 (zLO); LME (ar_trial_lme) on transitions with
% covariates: track pair, lap in block, spikes and speed (main and signed by
% track) and licks in both epochs; random intercepts animal, session, block.
% Variants: + grey_early and end_prev; + grey_late of the previous / next
% transition in the same block (temporal control); grey_late x is_switch; per
% animal; per-session partial Spearman; cross-region (one region's grey_late
% predicting the other's onset beyond its own).
% Time-resolved: correlation across transitions of residual zLO (block x pair
% means removed) between every pre-onset and post-onset 100 ms bin (event-locked
% bins, edges at onset), same transition and next transition in the block.
% Options: 'decoders', 'measure' ('zLO' default, 'LO' = raw log odds; the
% time-resolved part always uses zLO, the only per-bin values stored)
% Output: results\trial_prediction[_rawLO]\ (csv, mat, figures)
p = inputParser;
addParameter(p, 'decoders', {'bayes_multinomial', 'pls', 'bayes_multinomial_lobo5', 'pls_lobo5'});
addParameter(p, 'measure', 'zLO', @(x) ismember(x, {'zLO', 'LO'}));
parse(p, varargin{:});
O = p.Results;
P = ar_setup();
S = load(fullfile(P.results, 'data', 'epoch_table.mat'), 'E', 'TC'); E = S.E; TC = S.TC;
out_dir = fullfile(P.results, 'trial_prediction');
if strcmp(O.measure, 'LO')
    E.zLO = E.LO;
    out_dir = [out_dir '_rawLO'];
end
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
bad = find(ismember(P.subjects, P.V1_exclude_subjects));
regions = {'HPC', 'V1'};
rows = {};
R = struct();
for d = 1:numel(O.decoders)
    dec = O.decoders{d};
    W = struct();
    for r = 1:numel(regions)
        x = E(E.decoder == dec & E.region == regions{r}, :);
        if startsWith(regions{r}, 'V1'), x(ismember(x.animal, bad), :) = []; end
        W.(regions{r}) = ar_trial_table(x);
    end
    for r = 1:numel(regions)
        T = W.(regions{r});
        reg = regions{r};
        for resp = {'onset', 'onset_100', 'onset_200', 'onset_late'}
            [m, st] = ar_trial_lme(T, 'grey_late', {}, resp{1});
            rows(end + 1, :) = coef_row(dec, reg, ['main_' resp{1}], 'grey_late', st); %#ok<AGROW>
            if strcmp(resp{1}, 'onset'), R.(dec).(reg).main = m; end
        end
        [~, st] = ar_trial_lme(T, 'grey_late', {'grey_early', 'end_prev'});
        rows(end + 1, :) = coef_row(dec, reg, 'plus_early_end', 'grey_late', st); %#ok<AGROW>
        rows(end + 1, :) = coef_row(dec, reg, 'plus_early_end', 'grey_early', ar_lme_coef(st.lme, 'grey_early')); %#ok<AGROW>
        rows(end + 1, :) = coef_row(dec, reg, 'plus_early_end', 'end_prev', ar_lme_coef(st.lme, 'end_prev')); %#ok<AGROW>
        [~, st] = ar_trial_lme(T, 'grey_late', {'grey_late_prevtr', 'grey_late_nexttr'});
        rows(end + 1, :) = coef_row(dec, reg, 'neighbour_control', 'grey_late', st); %#ok<AGROW>
        rows(end + 1, :) = coef_row(dec, reg, 'neighbour_control', 'grey_late_prev_transition', ar_lme_coef(st.lme, 'grey_late_prevtr')); %#ok<AGROW>
        rows(end + 1, :) = coef_row(dec, reg, 'neighbour_control', 'grey_late_next_transition', ar_lme_coef(st.lme, 'grey_late_nexttr')); %#ok<AGROW>
        [~, st] = ar_trial_lme(T, 'grey_late', {'gl_x_switch'});
        rows(end + 1, :) = coef_row(dec, reg, 'switch_interaction', 'grey_late_x_switch', ar_lme_coef(st.lme, 'gl_x_switch')); %#ok<AGROW>
        for a = unique(T.animal)'
            [~, st] = ar_trial_lme(T(T.animal == a, :), 'grey_late', {});
            rows(end + 1, :) = coef_row(dec, reg, sprintf('animal_%s', P.subjects{a}), 'grey_late', st); %#ok<AGROW>
        end
        ps = session_partial(T);
        R.(dec).(reg).session_partial = ps;
        [ua2, ~, ia] = unique(ps.animal);
        am = accumarray(ia, ps.rho, [], @(v) mean(v, 'omitnan'));
        rows(end + 1, :) = {string(dec), string(reg), "session_partial_spearman", "grey_late", mean(ps.rho, 'omitnan'), ...
            NaN, NaN, signrank(ps.rho), sum(ps.rho > 0), sum(isfinite(ps.rho)), sum(am > 0), numel(ua2), ...
            string(strtrim(sprintf('%.3f ', am)))}; %#ok<AGROW>
    end
    % cross-region: one region's late grey predicting the other's onset
    keys = intersect(W.HPC(:, {'session', 'lap'}), W.V1(:, {'session', 'lap'}));
    [~, ih] = ismember(keys, W.HPC(:, {'session', 'lap'}));
    [~, iv] = ismember(keys, W.V1(:, {'session', 'lap'}));
    for target = {'HPC', 'V1'}
        other = setdiff({'HPC', 'V1'}, target);
        if strcmp(target{1}, 'HPC'), T = W.HPC(ih, :); T.other_grey_late = W.V1.grey_late(iv);
        else, T = W.V1(iv, :); T.other_grey_late = W.HPC.grey_late(ih); end
        [~, st] = ar_trial_lme(T, 'grey_late', {'other_grey_late'});
        rows(end + 1, :) = coef_row(dec, target{1}, ['cross_from_' other{1}], [other{1} '_grey_late'], ar_lme_coef(st.lme, 'other_grey_late')); %#ok<AGROW>
        rows(end + 1, :) = coef_row(dec, target{1}, ['cross_from_' other{1}], 'own_grey_late', st); %#ok<AGROW>
    end
    if strcmp(O.measure, 'zLO'), R.(dec).xtime = cross_time(TC, dec, regions, bad); end
end
C = cell2table(rows, 'VariableNames', {'decoder', 'region', 'model', 'term', 'estimate', 'ci_lo', 'ci_hi', 'p', ...
    'n_pos', 'n', 'n_anim_pos', 'n_anim', 'animal_values'});
writetable(C, fullfile(out_dir, 'trial_prediction_stats.csv'));
save(fullfile(out_dir, 'trial_prediction.mat'), 'R', 'C', '-v7.3');
if strcmp(O.measure, 'zLO'), plot_xtime(R, O.decoders, regions, out_dir); end
end

%% ------------------------------------------------------------------------
function row = coef_row(dec, reg, model, term, st)
row = {string(dec), string(reg), string(model), string(term), st.est, st.ci(1), st.ci(2), st.p, ...
    NaN, st.n, NaN, NaN, ""};
end

function ps = session_partial(T)
% Spearman between onset (0-200 ms) and grey_late after regressing the
% covariates (block, pair, spikes / speed / licks terms) out of both, per session
us = unique(T.session);
rho = NaN(numel(us), 1); an = NaN(numel(us), 1);
for i = 1:numel(us)
    t = T(T.session == us(i), :);
    X = [dummyvar(categorical(t.block)) dummyvar(t.pair) t.spk_onset t.spk_onset .* t.sgn_next t.spd_onset ...
        t.spd_onset .* t.sgn_next t.spk_grey_late t.spk_grey_late .* t.sgn_prev t.spd_grey_late ...
        t.spd_grey_late .* t.sgn_prev t.lick_onset t.lick_grey_late];
    ok = all(isfinite([X t.onset t.grey_late]), 2);
    an(i) = t.animal(1);
    if sum(ok) < 20, continue; end
    X = X(ok, :); X = X(:, std(X, 0, 1) > 0); X = [ones(sum(ok), 1) X];
    ry = t.onset(ok) - X * (X \ t.onset(ok));
    rx = t.grey_late(ok) - X * (X \ t.grey_late(ok));
    rho(i) = corr(rx, ry, 'Type', 'Spearman');
end
ps = table(us, an, rho, 'VariableNames', {'session', 'animal', 'rho'});
end

%% ------------------------------------------------------------------------
function X = cross_time(TC, dec, regions, bad)
% Residual zLO (onset-locked bins) with the block x pair mean removed per bin;
% Spearman across transitions between each pre-onset and post-onset bin, same
% transition and next transition in the block; curve = correlation of each
% pre-onset bin with the 0-200 ms onset mean. Mean Fisher z over sessions.
t = TC.t_on;
pre = find(t < 0); post = find(t > 0);
on = t > 0 & t < 0.2;
D = TC.trans;
X = struct('t', t, 'pre', t(pre), 'post', t(post));
for r = 1:numel(regions)
    sel = find(strcmp({TC.on.region}, regions{r}) & strcmp({TC.on.decoder}, dec));
    Msame = {}; Mnext = {}; Cs = []; Cn = []; an = []; ss = [];
    for i = 1:numel(sel)
        s = TC.on(sel(i)).session;
        Ds = D(D.session == s, :);
        if startsWith(regions{r}, 'V1') && ismember(Ds.animal(1), bad), continue; end
        z = double(TC.on(sel(i)).z);
        g = findgroups(Ds.block, Ds.prev_track, Ds.track);
        for k = unique(g)'
            z(g == k, :) = z(g == k, :) - mean(z(g == k, :), 1, 'omitnan');
        end
        y = mean(z(:, on), 2, 'omitnan');
        nxt = NaN(height(Ds), 1);
        for j = 1:height(Ds)
            q = find(Ds.lap == Ds.lap(j) + 1 & Ds.block == Ds.block(j), 1);
            if ~isempty(q), nxt(j) = q; end
        end
        okn = ~isnan(nxt);
        Msame{end + 1} = atanh(corr(z(:, pre), z(:, post), 'Type', 'Spearman', 'Rows', 'pairwise')); %#ok<AGROW>
        Mnext{end + 1} = atanh(corr(z(okn, pre), z(nxt(okn), post), 'Type', 'Spearman', 'Rows', 'pairwise')); %#ok<AGROW>
        Cs(end + 1, :) = atanh(corr(z(:, pre), y, 'Type', 'Spearman', 'Rows', 'pairwise'))'; %#ok<AGROW>
        Cn(end + 1, :) = atanh(corr(z(okn, pre), y(nxt(okn)), 'Type', 'Spearman', 'Rows', 'pairwise'))'; %#ok<AGROW>
        an(end + 1) = Ds.animal(1); ss(end + 1) = s; %#ok<AGROW>
    end
    X.(regions{r}).same = mean(cat(3, Msame{:}), 3, 'omitnan');
    X.(regions{r}).next = mean(cat(3, Mnext{:}), 3, 'omitnan');
    X.(regions{r}).curve_same = Cs; X.(regions{r}).curve_next = Cn;
    X.(regions{r}).animal = an; X.(regions{r}).session = ss;
end
end

function plot_xtime(R, decoders, regions, out_dir)
cmap = [linspace(0, 1, 32)' linspace(0, 1, 32)' ones(32, 1); ones(32, 1) linspace(1, 0, 32)' linspace(1, 0, 32)'];
for d = 1:numel(decoders)
    X = R.(decoders{d}).xtime;
    fig = figure('Position', [50 50 1500 380 * numel(regions)], 'Color', 'w');
    tl = tiledlayout(numel(regions), 3, 'TileSpacing', 'compact');
    for r = 1:numel(regions)
        Y = X.(regions{r});
        for k = 1:2
            nexttile(tl);
            M = Y.same; ttl = 'same transition';
            if k == 2, M = Y.next; ttl = 'pre of transition k, post of k+1 (control)'; end
            imagesc(X.post, X.pre, tanh(M), [-0.3 0.3]); axis xy; colorbar; colormap(gca, cmap);
            xlabel('Post-onset bin centre (s)'); ylabel('Pre-onset bin centre (s)');
            title(sprintf('%s: %s', regions{r}, ttl), 'FontWeight', 'normal');
        end
        nexttile(tl); hold on;
        for k = 1:2
            C = Y.curve_same; c = [0 0 0];
            if k == 2, C = Y.curve_next; c = [0.6 0.6 0.6]; end
            m = mean(C, 1, 'omitnan'); se = std(C, 0, 1, 'omitnan') / sqrt(size(C, 1));
            fill([X.pre fliplr(X.pre)], tanh([m + se fliplr(m - se)]), c, 'EdgeColor', 'none', 'FaceAlpha', 0.2);
            plot(X.pre, tanh(m), 'o-', 'Color', c, 'MarkerFaceColor', c, 'MarkerSize', 4, 'LineWidth', 1.5);
        end
        yline(0, 'k:');
        xlabel('Pre-onset bin centre (s)'); ylabel('Spearman r with onset bias (0-200 ms)');
        title(sprintf('%s: black = same transition, grey = next transition', regions{r}), 'FontWeight', 'normal');
        set(gca, 'TickDir', 'out', 'Box', 'off');
    end
    title(tl, sprintf('%s: residual bias before vs after onset, 100 ms bins (mean over sessions)', decoders{d}), 'Interpreter', 'none');
    ar_save_figure(fig, fullfile(out_dir, ['cross_time_' decoders{d}]));
    close(fig);
end
end
