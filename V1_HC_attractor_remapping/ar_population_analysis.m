function ar_population_analysis(varargin)
% Population track bias around lap transitions (epoch_table.mat):
%  1. time courses of the bias toward the upcoming track, aligned on the previous
%     lap's end and on track onset, by position of the upcoming lap in its block
%     (first lap of block = switch, 2, 3-5, 6+)
%  2. epoch summaries: bias toward the previous (block) track in the grey screen,
%     and onset decoding of switch vs late laps
%  3. grey-screen bias vs number of same-track laps just run (run length)
%  4. trial-by-trial chain end_prev -> grey_early -> grey_late -> onset on
%     residuals (zLO minus its mean for the same previous/upcoming track pair;
%     main version also demeaned within block and with spike count / speed
%     regressed out), and HPC-V1 co-variation of the residuals within each epoch
% Bias measures: zLO centred per session, region, decoder and epoch (or time
% point) by the average of the previous-T1 and previous-T2 means, then signed
% toward the previous or upcoming track. Every statistic is given per session,
% per animal and as an LME with animal and session random intercepts.
% Options: 'decoders' (default all variants of ar_epoch_table; _all decoders are
%          in-sample for lap epochs, so only their grey epochs are meaningful),
%          'regions' ({'HPC','V1'}),
%          'V1_all_animals' (false: V1 excludes P.V1_exclude_subjects),
%          'out_name' (results subfolder, default 'population'),
%          'measure' ('zLO' default; 'LO' = raw log odds), 'timecourses' (true)
% Output: results\<out_name>\<decoder>\ figures (.fig/.png/.pdf) and stats csv

p = inputParser;
addParameter(p, 'decoders', {'bayes_multinomial', 'bayes_poisson', 'pls', 'bayes_multinomial_lobo5', ...
    'bayes_poisson_lobo5', 'pls_lobo5', 'bayes_multinomial_all', 'pls_all'});
addParameter(p, 'regions', {'HPC', 'V1'});
addParameter(p, 'V1_all_animals', false);
addParameter(p, 'out_name', 'population');
addParameter(p, 'measure', 'zLO', @(x) ismember(x, {'zLO', 'LO'}));
addParameter(p, 'timecourses', true);
parse(p, varargin{:});
O = p.Results;
P = ar_setup();
S = load(fullfile(P.results, 'data', 'epoch_table.mat'), 'E', 'TC');
E = S.E; TC = S.TC;
if strcmp(O.measure, 'LO')
    % raw log odds (log P(T1)/P(T2), mean over the epoch's windows) instead of
    % the session z-scored, clipped zLO; time courses exist only as zLO
    E.zLO = E.LO;
end
E = add_measures(E, P, O);

for d = 1:numel(O.decoders)
    dec = O.decoders{d};
    out_dir = fullfile(P.results, O.out_name, dec);
    if ~exist(out_dir, 'dir'), mkdir(out_dir); end
    Ed = E(E.decoder == dec, :);
    if O.timecourses, plot_timecourses(TC, dec, O.regions, P, O, out_dir); end
    rows = {};
    rows = [rows; switch_locked(Ed, O.regions, out_dir, dec, P)]; %#ok<AGROW>
    rows = [rows; epoch_summary(Ed, O.regions, out_dir, dec)]; %#ok<AGROW>
    rows = [rows; run_length_analysis(Ed, O.regions, out_dir, dec)]; %#ok<AGROW>
    for v = {'resid', 'resid_blk', 'resid_ctrl_blk'}
        rows = [rows; chain_analysis(Ed, O.regions, out_dir, dec, v{1})]; %#ok<AGROW>
    end
    T = cell2table(rows, 'VariableNames', {'decoder', 'region', 'analysis', 'measure', 'mean', 'sem', ...
        'p_signrank', 'n_sess_pos', 'n_sess', 'n_anim_pos', 'n_anim', 'lme_est', 'lme_ci_lo', 'lme_ci_hi', ...
        'lme_p', 'n_trans', 'animal_means'});
    writetable(T, fullfile(out_dir, 'population_stats.csv'));
    close all
end
end

%% ------------------------------------------------------------------------
function E = add_measures(E, P, O)
E.sgn_prev = 2 * (E.prev_track == 1) - 1;
E.sgn_next = 2 * (E.track == 1) - 1;
E.group = repmat("lap6plus", height(E), 1);
E.group(E.lap_in_block == 1) = "switch";
E.group(E.lap_in_block == 2) = "lap2";
E.group(E.lap_in_block >= 3 & E.lap_in_block <= 5) = "lap3_5";
if ~O.V1_all_animals
    bad = find(ismember(P.subjects, P.V1_exclude_subjects));
    E(startsWith(E.region, "V1") & ismember(E.animal, bad), :) = [];
end
[g, ~] = findgroups(E.session, E.region, E.decoder, E.epoch);
c = splitapply(@(z, pt) balanced_mean(z, pt), E.zLO, E.prev_track, g);
E.zc = E.zLO - c(g);
E.toward_prev = E.zc .* E.sgn_prev;
E.toward_next = E.zc .* E.sgn_next;
% residual after removing the mean of the same previous x upcoming track pair
[g2, ~] = findgroups(E.session, E.region, E.decoder, E.epoch, E.prev_track, E.track);
m = splitapply(@(z) mean(z, 'omitnan'), E.zLO, g2);
E.resid = E.zLO - m(g2);
% same, also regressing out log spike count per window and speed within each
% group (log odds magnitude grows with spike count: shared gain or running
% speed fluctuations would otherwise correlate residuals across epochs/regions)
E.resid_ctrl = NaN(height(E), 1);
idx = accumarray(g2, (1:height(E))', [], @(x) {x});
for k = 1:numel(idx)
    i = idx{k};
    X = [ones(numel(i), 1) log1p(E.spikes_per_win(i)) E.speed(i)];
    y = E.zLO(i);
    ok = all(isfinite(X), 2) & isfinite(y);
    if sum(ok) < 6, continue; end
    b = X(ok, :) \ y(ok);
    E.resid_ctrl(i(ok)) = y(ok) - X(ok, :) * b;
end
% within-block versions: each block has its own held-out decoder, whose offset on
% out-of-distribution (grey) activity is shared by all transitions of the block;
% demeaning within block (non-switch transitions, >= 3 per block) removes it
[g3, ~] = findgroups(E.session, E.region, E.decoder, E.epoch, E.block, E.is_switch);
for v = ["resid", "resid_ctrl"]
    m = splitapply(@(z) mean(z, 'omitnan'), E.(v), g3);
    n = splitapply(@(z) sum(isfinite(z)), E.(v), g3);
    x = E.(v) - m(g3);
    x(n(g3) < 3 | E.is_switch) = NaN;
    E.(v + "_blk") = x;
end
% lick side bias toward the previous track (licks L - R signed)
E.lick_toward_prev = (E.licks_L - E.licks_R) .* E.sgn_prev;
% bias toward the upcoming track with log spike count and speed regressed out
% across all transitions of a session (per region, decoder, epoch): tests whether
% switch vs late differences survive matching for firing and running
E.toward_next_ctrl = NaN(height(E), 1);
[g4, ~] = findgroups(E.session, E.region, E.decoder, E.epoch);
idx4 = accumarray(g4, (1:height(E))', [], @(x) {x});
for k = 1:numel(idx4)
    i = idx4{k};
    X = [ones(numel(i), 1) log1p(E.spikes_per_win(i)) E.speed(i)];
    y = E.toward_next(i);
    ok = all(isfinite(X), 2) & isfinite(y);
    if sum(ok) < 10, continue; end
    b = X(ok, :) \ y(ok);
    E.toward_next_ctrl(i(ok)) = y(ok) - X(ok, 2:3) * b(2:3);
end
end

function c = balanced_mean(z, prev_track)
c = mean([mean(z(prev_track == 1), 'omitnan') mean(z(prev_track == 2), 'omitnan')], 'omitnan');
end

%% ------------------------------------------------------------------------
function plot_timecourses(TC, dec, regions, P, O, out_dir)
groups = ["switch", "lap2", "lap3_5", "lap6plus"];
gcol = [0.85 0.2 0.2; 0.95 0.6 0.1; 0.4 0.65 0.9; 0.1 0.25 0.6];
D = TC.trans;
grp = repmat("lap6plus", height(D), 1);
grp(D.lap_in_block == 1) = "switch"; grp(D.lap_in_block == 2) = "lap2";
grp(D.lap_in_block >= 3 & D.lap_in_block <= 5) = "lap3_5";
sgn_next = 2 * (D.track == 1) - 1;
aligns = {'off', 'on'};
fig = figure('Position', [50 50 1500 380 * numel(regions)], 'Color', 'w');
tl = tiledlayout(numel(regions), 3, 'TileSpacing', 'compact');
fig2 = figure('Position', [50 50 1500 380 * numel(regions)], 'Color', 'w');
tl2 = tiledlayout(numel(regions), 3, 'TileSpacing', 'compact');
acol = lines(6);
for r = 1:numel(regions)
    for a = 1:2
        A = TC.(aligns{a});
        tt = TC.(['t_' aligns{a}]);
        sel = find(strcmp({A.region}, regions{r}) & strcmp({A.decoder}, dec));
        sess_curves = NaN(numel(sel), numel(tt), numel(groups));
        sess_ids = NaN(numel(sel), 1);
        for i = 1:numel(sel)
            s = A(sel(i)).session;
            if startsWith(regions{r}, 'V1') && ~O.V1_all_animals && ...
                    ismember(D.animal(find(D.session == s, 1)), find(ismember(P.subjects, P.V1_exclude_subjects)))
                continue
            end
            rows = D.session == s;
            z = double(A(sel(i)).z);
            pt = D.prev_track(rows);
            c = mean([mean(z(pt == 1, :), 1, 'omitnan'); mean(z(pt == 2, :), 1, 'omitnan')], 1);
            tn = (z - c) .* sgn_next(rows);
            for g = 1:numel(groups)
                sess_curves(i, :, g) = mean(tn(grp(rows) == groups(g), :), 1, 'omitnan');
            end
            sess_ids(i) = s;
        end
        keep = ~isnan(sess_ids);
        sess_curves = sess_curves(keep, :, :); sess_ids = sess_ids(keep);
        figure(fig); nexttile(tl, (r - 1) * 3 + a); hold on;
        for g = 1:numel(groups)
            m = mean(sess_curves(:, :, g), 1, 'omitnan');
            se = std(sess_curves(:, :, g), 0, 1, 'omitnan') / sqrt(numel(sess_ids));
            fill([tt fliplr(tt)], [m + se fliplr(m - se)], gcol(g, :), 'EdgeColor', 'none', 'FaceAlpha', 0.2, 'HandleVisibility', 'off');
            plot(tt, m, 'o-', 'Color', gcol(g, :), 'MarkerFaceColor', gcol(g, :), 'MarkerSize', 3, 'LineWidth', 1.5);
        end
        yline(0, 'k:'); xline(0, 'k-'); xlim([-1 1]);
        if a == 1, xlabel('Time from previous lap end (s)'); else, xlabel('Time from track onset (s)'); end
        ylabel('Bias toward upcoming track (z)');
        title(sprintf('%s, n = %d sessions (100 ms bins)', regions{r}, numel(sess_ids)), 'FontWeight', 'normal');
        set(gca, 'TickDir', 'out', 'Box', 'off');
        if r == 1 && a == 1, legend(["first lap (switch)", "lap 2", "laps 3-5", "laps 6+"], 'Box', 'off', 'Location', 'northwest'); end
        % switch minus late, per animal
        diffc = sess_curves(:, :, 1) - sess_curves(:, :, 4);
        an = D.animal(arrayfun(@(s) find(D.session == s, 1), sess_ids));
        figure(fig2); nexttile(tl2, (r - 1) * 3 + a); hold on;
        ua = unique(an);
        for k = 1:numel(ua)
            plot(tt, mean(diffc(an == ua(k), :), 1, 'omitnan'), 'Color', acol(ua(k), :), 'LineWidth', 1);
        end
        plot(tt, mean(diffc, 1, 'omitnan'), 'k', 'LineWidth', 2);
        yline(0, 'k:'); xline(0, 'k-');
        if a == 1, xlabel('Time from previous lap end (s)'); else, xlabel('Time from track onset (s)'); end
        ylabel('Switch - laps 6+ (toward upcoming, z)');
        title(sprintf('%s: per animal (colour), mean (black)', regions{r}), 'FontWeight', 'normal');
        set(gca, 'TickDir', 'out', 'Box', 'off');
        if r == 1 && a == 1, legend([string(P.subjects(ua)) "mean"], 'Box', 'off', 'Location', 'southwest'); end
        if a == 2
            % zoom on onset 0-1.5 s for decoding accuracy build-up
            figure(fig); nexttile(tl, (r - 1) * 3 + 3); hold on;
            zoom = tt >= -0.5 & tt <= 1;
            for g = 1:numel(groups)
                m = mean(sess_curves(:, zoom, g), 1, 'omitnan');
                se = std(sess_curves(:, zoom, g), 0, 1, 'omitnan') / sqrt(numel(sess_ids));
                fill([tt(zoom) fliplr(tt(zoom))], [m + se fliplr(m - se)], gcol(g, :), 'EdgeColor', 'none', 'FaceAlpha', 0.2);
                plot(tt(zoom), m, 'Color', gcol(g, :), 'LineWidth', 1.5);
            end
            yline(0, 'k:'); xline(0, 'k-'); xline(0.2, 'k:');
            xlabel('Time from track onset (s)'); ylabel('Bias toward upcoming track (z)');
            title('Onset (zoom)', 'FontWeight', 'normal'); set(gca, 'TickDir', 'out', 'Box', 'off');
            figure(fig2); nexttile(tl2, (r - 1) * 3 + 3); hold on;
            scatter(an + 0.1 * randn(size(an)), mean(diffc(:, tt > 0 & tt < 0.2), 2), 25, acol(an, :), 'filled');
            scatter(an + 0.3, mean(diffc(:, tt > -1 & tt < 0), 2), 25, acol(an, :));
            yline(0, 'k:');
            set(gca, 'XTick', 1:6, 'XTickLabel', P.subjects, 'TickDir', 'out', 'Box', 'off');
            ylabel('Switch - laps 6+ (z)');
            title('Sessions: filled = onset 0-200 ms, open = grey -1-0 s', 'FontWeight', 'normal');
        end
    end
end
title(tl, sprintf('%s: bias toward the upcoming track by lap in block (mean +- SEM over sessions)', dec), 'Interpreter', 'none');
title(tl2, sprintf('%s: first lap of block minus laps 6+, per animal', dec), 'Interpreter', 'none');
ar_save_figure(fig, fullfile(out_dir, 'timecourse_toward_upcoming'));
ar_save_figure(fig2, fullfile(out_dir, 'timecourse_switch_minus_late_per_animal'));
end

%% ------------------------------------------------------------------------
function rows = switch_locked(Ed, regions, out_dir, dec, P)
% Switch-locked (regression-discontinuity) analysis, robust to slow drift.
% For every block switch (first lap k0 of a new block), transition index
% i = k - k0 (transition k = lap k-1 -> lap k): i <= -1 are inside the old block,
% i = 0 is the switch transition (previous lap old track, upcoming lap new track),
% i >= 1 are inside the new block (previous lap already the new track).
% y = zLO signed toward the NEW track (no centring: T1->T2 and T2->T1 switches
% cancel any decoder offset). Session curves = mean over switches.
%   step     = y(i = 1) - y(i = 0): change after one lap of the new track; in the
%              grey screen both epochs are ~one lap apart and decoded by the same
%              held-out decoder (block of the new track) or by the all-laps decoder
%   baseline = mean of within-block consecutive changes y(i+1) - y(i) for
%              i = -4..-2 and 1..4 (signed the same way): what drift alone gives
%   pre-step = y(i = 0) - y(i = -1): decoder changes here (held-out decoders),
%              the previous track does not
Lt = load(fullfile(P.results, 'data', 'lap_tables.mat'), 'L'); L = Lt.L;
idx = -5:6;
epochs = ["end_prev", "grey_early", "grey_late", "onset", "onset_late"];
rows = {};
acol = lines(6);
fig = figure('Position', [50 50 1800 360 * numel(regions)], 'Color', 'w');
tl = tiledlayout(numel(regions), numel(epochs) + 1, 'TileSpacing', 'compact');
for r = 1:numel(regions)
    step_stats = cell(1, numel(epochs));
    for e = 1:numel(epochs)
        x = Ed(Ed.region == regions{r} & Ed.epoch == epochs(e), :);
        us = unique(x.session);
        Y = NaN(numel(us), numel(idx)); an = NaN(numel(us), 1);
        for si = 1:numel(us)
            s = us(si); T = L(s).laps;
            xs = x(x.session == s, :);
            v = NaN(height(T), 1); v(xs.lap) = xs.zLO;
            k0s = find(T.is_switch);
            M = NaN(numel(k0s), numel(idx));
            for j = 1:numel(k0s)
                k0 = k0s(j); b = T.block(k0);
                sgn = 2 * (T.track(k0) == 1) - 1;
                for ii = 1:numel(idx)
                    k = k0 + idx(ii);
                    if k < 2 || k > height(T), continue; end
                    % transition k must lie in the old block (i < 0: laps k-1, k in b-1)
                    % or the new block (i > 0: laps k-1, k in b)
                    if idx(ii) < 0 && ~(T.block(k) == b - 1 && T.block(k - 1) == b - 1), continue; end
                    if idx(ii) > 0 && ~(T.block(k) == b && T.block(k - 1) == b), continue; end
                    M(j, ii) = v(k) * sgn;
                end
            end
            Y(si, :) = mean(M, 1, 'omitnan');
            an(si) = xs.animal(1);
        end
        % controls (stats only): the same step for transitions without licks in
        % the epoch, and for spikes per window, speed and licks toward the new
        % track (unsigned measures: they should not step with the previous track)
        % (toward_next_ctrl: already signed toward the upcoming = new track for i >= 0)
        ctrl = {'nolick', 'spikes_per_win', 'speed', 'lick_toward_new', 'toward_next_ctrl'};
        for c = 1:numel(ctrl)
            Yc = switch_curves(x, L, idx, ctrl{c});
            stc = vec_stats(Yc(:, idx == 1) - Yc(:, idx == 0), us, an);
            rows(end + 1, :) = stat_row(dec, regions{r}, "switch_step_ctrl_" + ctrl{c}, epochs(e), stc); %#ok<AGROW>
        end
        nexttile(tl); hold on;
        ua = unique(an);
        for k = 1:numel(ua)
            plot(idx, mean(Y(an == ua(k), :), 1, 'omitnan'), '-', 'Color', [acol(ua(k), :) 0.6], 'LineWidth', 0.8);
        end
        m = mean(Y, 1, 'omitnan'); se = std(Y, 0, 1, 'omitnan') ./ sqrt(sum(isfinite(Y), 1));
        errorbar(idx, m, se, 'ko-', 'MarkerFaceColor', 'k', 'CapSize', 0, 'LineWidth', 1.5);
        xline(-0.5, ':', 'Color', [0.5 0.5 0.5]); xline(0.5, 'k--'); yline(0, 'k:');
        xlabel('Transition from switch'); ylabel('Bias toward the new track (z)');
        % stats
        d = @(a, b) Y(:, idx == a) - Y(:, idx == b);
        step = d(1, 0);
        pre = d(0, -1);
        base = mean([d(-3, -4) d(-2, -3) d(-1, -2) d(2, 1) d(3, 2) d(4, 3) d(5, 4)], 2, 'omitnan');
        st = vec_stats(step, us, an); rows(end + 1, :) = stat_row(dec, regions{r}, 'switch_step_1_minus_0', epochs(e), st); %#ok<AGROW>
        st2 = vec_stats(step - base, us, an); rows(end + 1, :) = stat_row(dec, regions{r}, 'switch_step_minus_baseline', epochs(e), st2); %#ok<AGROW>
        st3 = vec_stats(pre, us, an); rows(end + 1, :) = stat_row(dec, regions{r}, 'switch_prestep_0_minus_m1', epochs(e), st3); %#ok<AGROW>
        st4 = vec_stats(base, us, an); rows(end + 1, :) = stat_row(dec, regions{r}, 'switch_baseline_step', epochs(e), st4); %#ok<AGROW>
        for ii = 1:numel(idx)
            sti = vec_stats(Y(:, ii), us, an);
            rows(end + 1, :) = stat_row(dec, regions{r}, "switch_curve_i" + idx(ii), epochs(e), sti); %#ok<AGROW>
        end
        title({sprintf('%s %s', regions{r}, epochs(e)), sprintf('step %.2f (p %.2g, %d/%d anim); vs base p %.2g', st.mean, ...
            st.p_signrank, st.n_anim_pos, st.n_anim, st2.p_signrank)}, 'Interpreter', 'none', 'FontWeight', 'normal', 'FontSize', 8);
        set(gca, 'TickDir', 'out', 'Box', 'off', 'XTick', idx);
        step_stats{e} = st2;
    end
    % summary: step minus baseline per epoch, animals
    nexttile(tl); hold on;
    for e = 1:numel(epochs)
        plot_stat(e, step_stats{e}, acol);
    end
    finish_axes(epochs, 'Step (1 - 0) minus baseline (z)', sprintf('%s: switch step beyond drift', regions{r}));
end
title(tl, {sprintf('%s: switch-locked bias toward the new track (black = mean +- SEM over sessions, colour = animals)', dec), ...
    'Transition i: lap k0+i-1 -> lap k0+i (k0 = first lap of new block). Dashed (0.5): previous track changes; grey dotted (-0.5): held-out decoder changes'}, ...
    'Interpreter', 'none');
ar_save_figure(fig, fullfile(out_dir, 'switch_locked'));
end

function Y = switch_curves(x, L, idx, what)
% Session x transition-index curves around switches (as in switch_locked) for
% 'nolick' (zLO toward the new track, transitions with licks in the epoch
% removed), 'spikes_per_win', 'speed' (unsigned) or 'lick_toward_new'
us = unique(x.session);
Y = NaN(numel(us), numel(idx));
for si = 1:numel(us)
    s = us(si); T = L(s).laps;
    xs = x(x.session == s, :);
    v = NaN(height(T), 1); signed = true;
    switch what
        case 'nolick'
            z = xs.zLO; z(xs.licks_L + xs.licks_R > 0) = NaN; v(xs.lap) = z;
        case 'lick_toward_new'
            v(xs.lap) = xs.licks_L - xs.licks_R;    % x sign of the new track below
        otherwise
            v(xs.lap) = xs.(what); signed = false;
    end
    k0s = find(T.is_switch);
    M = NaN(numel(k0s), numel(idx));
    for j = 1:numel(k0s)
        k0 = k0s(j); b = T.block(k0);
        sgn = 1;
        if signed, sgn = 2 * (T.track(k0) == 1) - 1; end
        for ii = 1:numel(idx)
            k = k0 + idx(ii);
            if k < 2 || k > height(T), continue; end
            if idx(ii) < 0 && ~(T.block(k) == b - 1 && T.block(k - 1) == b - 1), continue; end
            if idx(ii) > 0 && ~(T.block(k) == b && T.block(k - 1) == b), continue; end
            M(j, ii) = v(k) * sgn;
        end
    end
    Y(si, :) = mean(M, 1, 'omitnan');
end
end

function st = vec_stats(v, us, an)
% session values v -> session_vector_stats-like struct (LME across animals)
st = session_vector_stats(v, us, an);
end

%% ------------------------------------------------------------------------
function rows = epoch_summary(Ed, regions, out_dir, dec)
% (a) bias toward the previous track per epoch (all transitions; the animal
%     cannot know a switch is coming, so pre-switch transitions are included)
% (b) bias toward the upcoming track and Bayesian p(true track): switch - late
epochs = ["end_prev", "grey_early", "grey_late", "onset", "onset_late", "lap_rest"];
rows = {};
fig = figure('Position', [50 50 1600 370 * numel(regions)], 'Color', 'w');
tl = tiledlayout(numel(regions), 4, 'TileSpacing', 'compact');
acol = lines(6);
for r = 1:numel(regions)
    Er = Ed(Ed.region == regions{r}, :);
    % (a) toward previous track, by epoch
    nexttile(tl); hold on;
    for e = 1:numel(epochs)
        x = Er(Er.epoch == epochs(e), :);
        st = ar_stats(x.toward_prev, x.session, x.animal);
        rows(end + 1, :) = stat_row(dec, regions{r}, 'toward_prev_all', epochs(e), st); %#ok<AGROW>
        plot_stat(e, st, acol);
    end
    finish_axes(epochs, 'Bias toward previous lap''s track (z)', sprintf('%s: all transitions', regions{r}));
    % non-switch only (previous = upcoming), grey epochs
    nexttile(tl); hold on;
    for e = 1:numel(epochs)
        x = Er(Er.epoch == epochs(e) & ~Er.is_switch, :);
        st = ar_stats(x.toward_prev, x.session, x.animal);
        rows(end + 1, :) = stat_row(dec, regions{r}, 'toward_prev_nonswitch', epochs(e), st); %#ok<AGROW>
        x = Er(Er.epoch == epochs(e) & Er.is_switch, :);
        st2 = ar_stats(x.toward_prev, x.session, x.animal);
        rows(end + 1, :) = stat_row(dec, regions{r}, 'toward_prev_switch', epochs(e), st2); %#ok<AGROW>
        plot_stat(e - 0.15, st, acol); plot_stat(e + 0.15, st2, acol, 'r');
    end
    finish_axes(epochs, 'Bias toward previous lap''s track (z)', 'black = no switch, red = switch (prev = old track)');
    % (b) switch - late, toward upcoming
    nexttile(tl); hold on;
    for e = 1:numel(epochs)
        d = paired_session_diff(Er(Er.epoch == epochs(e), :), 'toward_next');
        rows(end + 1, :) = stat_row(dec, regions{r}, 'toward_next_switch_minus_late', epochs(e), d); %#ok<AGROW>
        plot_stat(e, d, acol);
    end
    finish_axes(epochs, 'Switch - laps 6+: toward upcoming (z)', 'per session (dots, colour = animal)');
    % (c) Bayesian p(true track) at onset epochs: switch vs late
    nexttile(tl); hold on;
    ep2 = ["end_prev", "onset", "onset_late", "lap_rest"];
    for e = 1:numel(ep2)
        x = Er(Er.epoch == ep2(e), :);
        if all(isnan(x.p_true_track)), continue; end
        d = paired_session_diff(x, 'p_true_track');
        rows(end + 1, :) = stat_row(dec, regions{r}, 'p_true_switch_minus_late', ep2(e), d); %#ok<AGROW>
        plot_stat(e, d, acol);
    end
    finish_axes(ep2, 'Switch - laps 6+: p(true track)', 'Bayesian posterior on the lap''s own track');
    % confound checks (stats only): no-lick transitions, and speed / spike count
    % differences between switch and late laps
    for e = 1:numel(epochs)
        x = Er(Er.epoch == epochs(e) & Er.licks_L + Er.licks_R == 0, :);
        st = ar_stats(x.toward_prev, x.session, x.animal);
        rows(end + 1, :) = stat_row(dec, regions{r}, 'toward_prev_all_nolick', epochs(e), st); %#ok<AGROW>
        x = Er(Er.epoch == epochs(e), :);
        rows(end + 1, :) = stat_row(dec, regions{r}, 'toward_next_switch_minus_late_ctrl', epochs(e), ...
            paired_session_diff(x, 'toward_next_ctrl')); %#ok<AGROW>
        rows(end + 1, :) = stat_row(dec, regions{r}, 'speed_switch_minus_late', epochs(e), paired_session_diff(x, 'speed')); %#ok<AGROW>
        rows(end + 1, :) = stat_row(dec, regions{r}, 'spikes_switch_minus_late', epochs(e), paired_session_diff(x, 'spikes_per_win')); %#ok<AGROW>
        rows(end + 1, :) = stat_row(dec, regions{r}, 'lick_toward_prev', epochs(e), ar_stats(x.lick_toward_prev, x.session, x.animal)); %#ok<AGROW>
    end
end
title(tl, sprintf('%s: epoch summaries (bars = mean +- SEM over sessions, LME p in title of points)', dec), 'Interpreter', 'none');
ar_save_figure(fig, fullfile(out_dir, 'epoch_summary'));
end

function d = paired_session_diff(x, var)
% switch minus laps 6+ per session (session means), as an ar_stats-like struct
sw = x(x.group == "switch", :); la = x(x.group == "lap6plus", :);
[us, ~] = unique(x.session);
sm = NaN(numel(us), 1); sa = NaN(numel(us), 1);
for i = 1:numel(us)
    sm(i) = mean(sw.(var)(sw.session == us(i)), 'omitnan') - mean(la.(var)(la.session == us(i)), 'omitnan');
    sa(i) = x.animal(find(x.session == us(i), 1));
end
ok = isfinite(sm);
d = struct();
d.sessions = us(ok); d.session_mean = sm(ok); d.session_animal = sa(ok);
[ua, ~, ia] = unique(d.session_animal);
d.animals = ua; d.animal_mean = accumarray(ia, d.session_mean, [], @mean);
d.mean = mean(d.session_mean); d.sem = std(d.session_mean) / sqrt(sum(ok));
d.p_signrank = signrank(d.session_mean);
d.n_sess_pos = sum(d.session_mean > 0); d.n_sess = sum(ok);
d.n_anim_pos = sum(d.animal_mean > 0); d.n_anim = numel(ua);
% transition-level LME: var ~ is_switch (switch vs laps 6+)
y = [sw.(var); la.(var)];
tbl = table(y, [true(height(sw), 1); false(height(la), 1)], categorical([sw.session; la.session]), ...
    categorical([sw.animal; la.animal]), 'VariableNames', {'y', 'is_switch', 'session', 'animal'});
d.n_trans = sum(isfinite(y));
try
    lme = fitlme(tbl, 'y ~ is_switch + (1|animal) + (1|animal:session)');
    k = strcmp(lme.Coefficients.Name, 'is_switch_1');
    ci = coefCI(lme);
    d.lme_est = lme.Coefficients.Estimate(k); d.lme_p = lme.Coefficients.pValue(k); d.lme_ci = ci(k, :);
catch
    d.lme_est = NaN; d.lme_p = NaN; d.lme_ci = [NaN NaN];
end
end

%% ------------------------------------------------------------------------
function rows = run_length_analysis(Ed, regions, out_dir, dec)
% Grey-screen bias toward the previous track vs how many laps of that track the
% animal has just run (run_len_prev); pre-switch vs other transitions at matched
% run length (the animal cannot know a switch is coming: should not differ)
rows = {};
epochs = ["end_prev", "grey_early", "grey_late", "onset"];
rl_bins = {1, 2, 3:4, 5:6, 7:8, 9:13};
rl_lab = {'1', '2', '3-4', '5-6', '7-8', '9+'};
fig = figure('Position', [50 50 1600 370 * numel(regions)], 'Color', 'w');
tl = tiledlayout(numel(regions), numel(epochs), 'TileSpacing', 'compact');
for r = 1:numel(regions)
    for e = 1:numel(epochs)
        x = Ed(Ed.region == regions{r} & Ed.epoch == epochs(e), :);
        nexttile(tl); hold on;
        for pre_sw = [false true]
            xs = x(x.is_switch == pre_sw, :);
            [us] = unique(x.session);
            M = NaN(numel(us), numel(rl_bins));
            for i = 1:numel(us)
                for b = 1:numel(rl_bins)
                    M(i, b) = mean(xs.toward_prev(xs.session == us(i) & ismember(xs.run_len_prev, rl_bins{b})), 'omitnan');
                end
            end
            n = sum(~isnan(M), 1);
            m = mean(M, 1, 'omitnan'); se = std(M, 0, 1, 'omitnan') ./ sqrt(n);
            if pre_sw, c = [0.85 0.2 0.2]; else, c = [0.1 0.1 0.1]; end
            errorbar((1:numel(rl_bins)) + 0.1 * pre_sw, m, se, 'o-', 'Color', c, 'MarkerFaceColor', c, 'CapSize', 0);
        end
        % slope of bias on log run length, non-switch transitions only
        xn = x(~x.is_switch, :);
        cov = table(log(xn.run_len_prev) - mean(log(xn.run_len_prev)), 'VariableNames', {'log_run'});
        st = ar_stats(xn.toward_prev, xn.session, xn.animal, 'covariates', cov);
        slope = NaN; ps = NaN;
        if ~isempty(st.lme)
            k = strcmp(st.lme.Coefficients.Name, 'log_run');
            slope = st.lme.Coefficients.Estimate(k); ps = st.lme.Coefficients.pValue(k);
            ci = coefCI(st.lme);
            stl = st; stl.lme_est = slope; stl.lme_p = ps; stl.lme_ci = ci(k, :);
            stl.mean = NaN; stl.sem = NaN; stl.p_signrank = NaN;
            rows(end + 1, :) = stat_row(dec, regions{r}, 'toward_prev_vs_log_run_length', epochs(e), stl); %#ok<AGROW>
        end
        yline(0, 'k:');
        set(gca, 'XTick', 1:numel(rl_bins), 'XTickLabel', rl_lab, 'TickDir', 'out', 'Box', 'off');
        xlabel('Laps of previous track in a row'); ylabel('Bias toward previous track (z)');
        title(sprintf('%s %s: slope %.3f / log lap, p %.3g', regions{r}, epochs(e), slope, ps), 'FontWeight', 'normal', 'Interpreter', 'none');
    end
end
title(tl, sprintf('%s: bias toward previous track vs run length (black = next lap same track, red = next lap switches)', dec), 'Interpreter', 'none');
ar_save_figure(fig, fullfile(out_dir, 'run_length'));
end

%% ------------------------------------------------------------------------
function rows = chain_analysis(Ed, regions, out_dir, dec, var)
% Trial-by-trial: Spearman correlation across transitions (within session) of
% residual zLO between epochs, and between HPC and V1 in the same epoch.
% Residuals remove the mean of each previous x upcoming track pair, so block
% structure alone cannot create correlations ('resid'); '_ctrl' also regresses
% out log spike count and speed within each pair; '_blk' demeans within block
% (removes the per-block decoder offset; non-switch transitions only).
rows = {};
Ed.resid = Ed.(var);
epochs = ["end_prev", "grey_early", "grey_late", "onset", "onset_late"];
nE = numel(epochs);
fig = figure('Position', [50 50 1400 420], 'Color', 'w');
tl = tiledlayout(1, numel(regions) + 1, 'TileSpacing', 'compact');
acol = lines(6);
for r = 1:numel(regions)
    x = Ed(Ed.region == regions{r}, {'session', 'animal', 'lap', 'epoch', 'resid'});
    x.epoch = categorical(x.epoch);
    W = unstack(x, 'resid', 'epoch');
    us = unique(W.session);
    Rm = NaN(nE, nE, numel(us)); an = NaN(numel(us), 1);
    for i = 1:numel(us)
        w = W(W.session == us(i), :);
        A = w{:, cellstr(epochs)};
        Rm(:, :, i) = corr(A, 'Type', 'Spearman', 'Rows', 'pairwise');
        an(i) = w.animal(1);
    end
    Zm = atanh(min(max(Rm, -0.999), 0.999));
    nexttile(tl);
    imagesc(tanh(mean(Zm, 3, 'omitnan')), [-0.3 0.3]); axis square; colorbar;
    colormap(gca, redblue());
    set(gca, 'XTick', 1:nE, 'XTickLabel', epochs, 'YTick', 1:nE, 'YTickLabel', epochs, 'TickLabelInterpreter', 'none');
    title(sprintf('%s: residual correlation (mean over %d sessions)', regions{r}, numel(us)), 'FontWeight', 'normal');
    for a = 1:nE - 1
        b = a + 1;
        z = squeeze(Zm(a, b, :));
        st = session_vector_stats(z, us, an);
        rows(end + 1, :) = stat_row(dec, regions{r}, var + "_corr_" + epochs(a), epochs(b), st); %#ok<AGROW>
        text(b, a, sprintf('%.2f\np=%.2g', tanh(st.mean), st.p_signrank), 'HorizontalAlignment', 'center', 'FontSize', 7);
    end
    z = squeeze(Zm(1, 4, :)); st = session_vector_stats(z, us, an);
    rows(end + 1, :) = stat_row(dec, regions{r}, var + "_corr_end_prev", "onset", st); %#ok<AGROW>
end
% HPC vs V1 residuals, same epoch
if numel(regions) >= 2
    nexttile(tl); hold on;
    for e = 1:nE
        xh = Ed(Ed.region == regions{1} & Ed.epoch == epochs(e), {'session', 'animal', 'lap', 'resid'});
        xv = Ed(Ed.region == regions{2} & Ed.epoch == epochs(e), {'session', 'lap', 'resid'});
        xh.Properties.VariableNames{'resid'} = 'r1'; xv.Properties.VariableNames{'resid'} = 'r2';
        J = innerjoin(xh, xv, 'Keys', {'session', 'lap'});
        us = unique(J.session); z = NaN(numel(us), 1); an = NaN(numel(us), 1);
        for i = 1:numel(us)
            j = J(J.session == us(i), :);
            z(i) = atanh(corr(j.r1, j.r2, 'Type', 'Spearman', 'Rows', 'complete'));
            an(i) = j.animal(1);
        end
        st = session_vector_stats(z, us, an);
        rows(end + 1, :) = stat_row(dec, 'HPC_V1', var + "_corr_HPC_V1", epochs(e), st); %#ok<AGROW>
        plot_stat(e, st, acol);
    end
    finish_axes(epochs, 'HPC-V1 residual correlation (Fisher z)', 'Trial-by-trial HPC-V1 co-variation');
end
title(tl, sprintf('%s: trial-by-trial chain of track bias residuals (%s)', dec, var), 'Interpreter', 'none');
ar_save_figure(fig, fullfile(out_dir, ['chain_correlations_' var]));
end

function st = session_vector_stats(z, us, an)
ok = isfinite(z);
z = z(ok); us = us(ok); an = an(ok);
st.sessions = us; st.session_mean = z; st.session_animal = an;
[ua, ~, ia] = unique(an);
st.animals = ua; st.animal_mean = accumarray(ia, z, [], @mean);
st.mean = mean(z); st.sem = std(z) / sqrt(numel(z));
st.p_signrank = signrank(z);
st.n_sess_pos = sum(z > 0); st.n_sess = numel(z);
st.n_anim_pos = sum(st.animal_mean > 0); st.n_anim = numel(ua);
% animal-level LME on session values
tbl = table(z, categorical(an), 'VariableNames', {'y', 'animal'});
st.n_trans = NaN;
try
    lme = fitlme(tbl, 'y ~ 1 + (1|animal)');
    ci = coefCI(lme);
    st.lme_est = lme.Coefficients.Estimate(1); st.lme_p = lme.Coefficients.pValue(1); st.lme_ci = ci(1, :);
catch
    st.lme_est = NaN; st.lme_p = NaN; st.lme_ci = [NaN NaN];
end
end

%% ------------------------------------------------------------------------
function row = stat_row(dec, region, analysis, measure, st)
am = sprintf('%.3f ', st.animal_mean);
row = {string(dec), string(region), string(analysis), string(measure), st.mean, st.sem, st.p_signrank, ...
    st.n_sess_pos, st.n_sess, st.n_anim_pos, st.n_anim, st.lme_est, st.lme_ci(1), st.lme_ci(2), st.lme_p, ...
    st.n_trans, string(strtrim(am))};
end

function plot_stat(x, st, acol, colr)
if nargin < 4, colr = 'k'; end
jit = 0.08 * randn(numel(st.session_mean), 1);
scatter(x + jit, st.session_mean, 12, acol(st.session_animal, :), 'filled', 'MarkerFaceAlpha', 0.6);
scatter(x + 0.25 + zeros(size(st.animal_mean)), st.animal_mean, 25, acol(st.animals, :), 'd', 'filled');
errorbar(x, st.mean, st.sem, 'o', 'Color', colr, 'MarkerFaceColor', colr, 'CapSize', 0, 'LineWidth', 1.5);
if st.lme_p < 0.05, s = '*'; else, s = ''; end
text(x, max(st.session_mean) + 0.02, sprintf('%s\n%.2g', s, st.lme_p), 'HorizontalAlignment', 'center', 'FontSize', 7);
end

function finish_axes(labels, ylab, ttl)
yline(0, 'k:');
set(gca, 'XTick', 1:numel(labels), 'XTickLabel', labels, 'TickDir', 'out', 'Box', 'off', 'TickLabelInterpreter', 'none');
xlim([0.5 numel(labels) + 0.5]);
ylabel(ylab); title(ttl, 'FontWeight', 'normal', 'Interpreter', 'none');
end

function c = redblue()
c = [linspace(0, 1, 32)' linspace(0, 1, 32)' ones(32, 1); ones(32, 1) linspace(1, 0, 32)' linspace(1, 0, 32)'];
end
