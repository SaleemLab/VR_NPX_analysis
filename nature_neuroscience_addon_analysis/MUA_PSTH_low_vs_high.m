%% Dominant / non-dominant V1 MUA PSTH around ripple onset: Low vs High (Peak vs Trough), POST
% Firing-rate check for the condition comparisons (R1.2, R4.2/4/6): is dominant-
% (or non-dominant-) hemisphere V1 MUA around ripples different between the two
% event groups? All events, no MUA matching.
%
% Dominant hemisphere from the V1 bias sign in the dominance window (as
% MUA_equalised_geomean): V1 > 0 -> V1 R, V1 < 0 -> V1 L; the non-dominant
% hemisphere is the other one. Done separately for dominance from the pre
% (-0.2-0 s), post (0-0.2 s) and full (-0.2-0.2 s) V1 window.
% Conditions (within session, as MUA_equalised_geomean):
%   ripple_power : ripple power percentile <= 25 (Low) vs > 75 (High)
%   spindle      : dominant-hemisphere spindle percentile <= 25 (Low) vs > 75 (High)
%   SO_phase     : dominant-hemisphere SO peak [-pi/2, pi/2] vs trough
% MUA: z-scored, 10 ms bins, -2 to 2 s around ripple onset (ripple_MUA_V1_L/R).
% Per session: mean PSTH per condition, and mean MUA in -0.2-0 s and 0-0.2 s.
% Across sessions: mean +- SE PSTH; paired signed-rank test High vs Low
% (Trough vs Peak) on the session means, for each MUA window.
% Sessions need >= min_events events in both conditions.
%   hemisphere : 'dominant' (default) -> MUA_equalised_event_control\MUA_PSTH
%                'non_dominant' -> MUA_equalised_event_control_with_nondominant\MUA_PSTH
% Output: nature_neuroscience_revision_analysis\<folder>\MUA_PSTH
function MUA_PSTH_low_vs_high(hemisphere)
    if nargin < 1, hemisphere = 'dominant'; end
    switch hemisphere
        case 'dominant'
            v = MUA_equalised_variant('dominant'); hemi_label = 'Dominant V1 MUA';
        case 'non_dominant'
            v = MUA_equalised_variant('with_nondominant'); hemi_label = 'Non-dominant V1 MUA';
        otherwise
            error('MUA_PSTH_low_vs_high:hemisphere', 'Unknown hemisphere %s', hemisphere);
    end
    rev_dir = addon_setup();
    out_dir = fullfile(rev_dir, v.folder, 'MUA_PSTH');
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end
    min_events = 10;          % per condition per session
    plot_range = [-0.5 0.5];  % s, PSTH x-axis

    % --- Load non-duplicated POST data (rows aligned to ripple_info) ---
    load(fullfile(rev_dir, 'ripple_info.mat'), 'ripple_info');
    load(fullfile(rev_dir, 'KDE_reactivation_ripples_PSTH_POST_unique.mat'), 'KDE_reactivation_ripples_PSTH');
    V1 = (KDE_reactivation_ripples_PSTH.V1_z_logodds_ripples + KDE_reactivation_ripples_PSTH.nan_mask)';
    clear KDE_reactivation_ripples_PSTH
    load(fullfile(rev_dir, 'ripple_MUA_V1_L.mat'), 'ripple_MUA_V1_L');
    load(fullfile(rev_dir, 'ripple_MUA_V1_R.mat'), 'ripple_MUA_V1_R');
    n_ev = numel(ripple_info.session_id);
    assert(size(ripple_MUA_V1_L, 1) == n_ev && size(V1, 2) == n_ev, 'MUA / PSTH rows do not match ripple_info');

    % --- Windows ---
    timebin = 0.01;
    bin_edges = -1:timebin:1;
    bin_centers = bin_edges(1:end-1) + timebin/2;
    mua_centers = -2 + timebin/2:timebin:2 - timebin/2;
    in_plot = mua_centers > plot_range(1) & mua_centers < plot_range(2);
    t = mua_centers(in_plot);
    dom_windows = struct( ...
        'name',  {'pre', 'post', 'full'}, ...
        'label', {'dominance -0.2-0 s', 'dominance 0-0.2 s', 'dominance -0.2-0.2 s'}, ...
        'range', {[-0.2 0], [0 0.2], [-0.2 0.2]});
    mua_windows = struct( ...
        'name',  {'pre', 'post'}, ...
        'label', {'MUA -0.2-0 s', 'MUA 0-0.2 s'}, ...
        'range', {[-0.2 0], [0 0.2]});

    sess = ripple_info.session_id(:);
    sessions = unique(sess)';
    n_sess = numel(sessions);
    ripple_pct = within_session_percentile(ripple_info.ripple_power(:), sess);
    spindle_pct = within_session_percentile(ripple_info.spindle_percentile, sess);

    analyses = {'ripple_power', 'spindle', 'SO_phase'};
    titles = {'Ripple power', 'Spindle (dom. hemi)', 'SO phase (dom. hemi)'};

    R = struct();   % per analysis x dominance window: session PSTHs and window means
    T = table();    % per session x dominance window x analysis x condition x MUA window
    for w = 1:numel(dom_windows)
        wn = dom_windows(w).name;
        in_wb = bin_centers > dom_windows(w).range(1) & bin_centers < dom_windows(w).range(2);
        V1_b = mean(V1(in_wb, :), 'omitnan')';
        hemi = NaN(n_ev, 1); hemi(V1_b < 0) = 1; hemi(V1_b > 0) = 2;
        MUA = NaN(n_ev, numel(mua_centers), 'single');
        if strcmp(hemisphere, 'dominant')
            MUA(hemi == 1, :) = ripple_MUA_V1_L(hemi == 1, :);
            MUA(hemi == 2, :) = ripple_MUA_V1_R(hemi == 2, :);
        else
            MUA(hemi == 1, :) = ripple_MUA_V1_R(hemi == 1, :);
            MUA(hemi == 2, :) = ripple_MUA_V1_L(hemi == 2, :);
        end
        valid = ~isnan(hemi) & all(isfinite(MUA), 2);
        win_MUA = NaN(n_ev, numel(mua_windows));
        for m = 1:numel(mua_windows)
            in_m = mua_centers > mua_windows(m).range(1) & mua_centers < mua_windows(m).range(2);
            win_MUA(:, m) = double(mean(MUA(:, in_m), 2));
        end

        for iA = 1:numel(analyses)
            [masks, labels] = condition_masks(analyses{iA}, ripple_info, ripple_pct, spindle_pct, hemi);
            psth = NaN(n_sess, numel(t), 2);   % session x time x condition
            psth_se = NaN(n_sess, numel(t), 2);
            wm = NaN(n_sess, numel(mua_windows), 2);
            n_c = zeros(n_sess, 2);
            for si = 1:n_sess
                s = sessions(si);
                for c = 1:2
                    i = masks{c} & valid & sess == s;
                    n_c(si, c) = sum(i);
                    if n_c(si, c) == 0, continue; end
                    x = double(MUA(i, in_plot));
                    psth(si, :, c) = mean(x, 1);
                    psth_se(si, :, c) = std(x, 0, 1) / sqrt(size(x, 1));
                    wm(si, :, c) = mean(win_MUA(i, :), 1);
                    for m = 1:numel(mua_windows)
                        T = [T; table(s, string(wn), string(analyses{iA}), string(labels{c}), ...
                            string(mua_windows(m).name), n_c(si, c), wm(si, m, c), ...
                            'VariableNames', {'session_id', 'dominance_window', 'analysis', 'condition', ...
                            'MUA_window', 'n_events', 'mean_MUA'})]; %#ok<AGROW>
                    end
                end
            end
            keep = all(n_c >= min_events, 2);
            R.(sprintf('%s_%s', analyses{iA}, wn)) = struct('labels', {labels}, 'sessions', sessions, ...
                'keep', keep, 'n_events', n_c, 'psth', psth, 'psth_se', psth_se, 'window_MUA', wm, 't', t);
        end
    end

    % --- Session-level paired tests (B vs A) ---
    S = table();
    for w = 1:numel(dom_windows)
        for iA = 1:numel(analyses)
            D = R.(sprintf('%s_%s', analyses{iA}, dom_windows(w).name));
            for m = 1:numel(mua_windows)
                a = D.window_MUA(D.keep, m, 1); b = D.window_MUA(D.keep, m, 2);
                [p2, p1l, p1r] = deal(NaN);
                if numel(a) >= 3
                    p2 = signrank(b, a);
                    p1r = signrank(b, a, 'tail', 'right');
                    p1l = signrank(b, a, 'tail', 'left');
                end
                S = [S; table(string(analyses{iA}), string(dom_windows(w).name), string(mua_windows(m).name), ...
                    string(sprintf('%s vs %s', D.labels{2}, D.labels{1})), numel(a), ...
                    mean(a), std(a) / sqrt(numel(a)), mean(b), std(b) / sqrt(numel(b)), mean(b - a), ...
                    p2, p1r, p1l, ...
                    'VariableNames', {'analysis', 'dominance_window', 'MUA_window', 'comparison', 'n_sessions', ...
                    'mean_A', 'SE_A', 'mean_B', 'SE_B', 'mean_diff', ...
                    'p_signrank_two', 'p_signrank_B_gt_A', 'p_signrank_B_lt_A'})]; %#ok<AGROW>
            end
        end
    end
    disp(S);

    save(fullfile(out_dir, 'MUA_PSTH_low_vs_high.mat'), 'R', 'T', 'S', 'min_events', 'dom_windows', 'mua_windows', 'hemisphere');
    writetable(T, fullfile(out_dir, 'MUA_PSTH_session_window_means.csv'));
    writetable(S, fullfile(out_dir, 'MUA_PSTH_session_stats.csv'));

    for w = 1:numel(dom_windows)
        for iA = 1:numel(analyses)
            key = sprintf('%s_%s', analyses{iA}, dom_windows(w).name);
            plot_session_PSTHs(R.(key), sprintf('%s, %s', titles{iA}, dom_windows(w).label), ...
                mua_windows, hemi_label, fullfile(out_dir, ['MUA_PSTH_sessions_' key]));
        end
    end
    plot_mean_PSTH(R, S, analyses, titles, dom_windows, mua_windows, hemi_label, fullfile(out_dir, 'MUA_PSTH_mean_SE'));
    plot_window_means(R, S, analyses, titles, dom_windows, mua_windows, hemi_label, fullfile(out_dir, 'MUA_PSTH_window_means'));
    fprintf('Saved to %s\n', out_dir);
end

function pct = within_session_percentile(X, session_ids)
    % Percentile (0-100) of each column of X within each session
    pct = NaN(size(X));
    for s = unique(session_ids)'
        idx = session_ids == s;
        for c = 1:size(X, 2)
            x = X(idx, c);
            ok = isfinite(x);
            p = NaN(size(x));
            p(ok) = (tiedrank(x(ok)) - 0.5) ./ sum(ok) * 100;
            pct(idx, c) = p;
        end
    end
end

function [masks, labels] = condition_masks(analysis, info, ripple_pct, spindle_pct, hemi)
    % Logical event masks {A, B}; spindle / SO from the dominant hemisphere
    n = numel(hemi);
    switch analysis
        case 'ripple_power'
            masks = {ripple_pct <= 25, ripple_pct > 75};
            labels = {'Low', 'High'};
        case 'spindle'
            f = dominant_feature(spindle_pct, hemi, n);
            masks = {f <= 25, f > 75};
            labels = {'Low', 'High'};
        case 'SO_phase'
            f = dominant_feature(info.SO_phase, hemi, n);
            masks = {f >= -pi/2 & f <= pi/2, (f >= -pi & f < -pi/2) | (f > pi/2 & f <= pi)};
            labels = {'Peak', 'Trough'};
    end
end

function f = dominant_feature(F, hemi, n)
    f = NaN(n, 1);
    f(hemi == 1) = F(hemi == 1, 1);
    f(hemi == 2) = F(hemi == 2, 2);
end

function shade_windows(ax, mua_windows)
    yl = ylim(ax);
    win_colors = {[0.3 0.6 1], [1 0.7 0.2]};
    for m = 1:numel(mua_windows)
        r = mua_windows(m).range;
        patch(ax, [r(1) r(2) r(2) r(1)], [yl(1) yl(1) yl(2) yl(2)], win_colors{m}, ...
            'FaceAlpha', 0.12, 'EdgeColor', 'none', 'HandleVisibility', 'off', 'Tag', 'window_shade');
    end
    uistack(findobj(ax, 'Tag', 'window_shade'), 'bottom');
    ylim(ax, yl);
end

function plot_band(ax, t, mu, se, color)
    fill(ax, [t fliplr(t)], [mu + se fliplr(mu - se)], color, 'FaceAlpha', 0.25, 'EdgeColor', 'none', ...
        'HandleVisibility', 'off');
    plot(ax, t, mu, '-', 'Color', color, 'LineWidth', 1.5);
end

function save_figure(fig, save_path)
    savefig(fig, [save_path '.fig']);
    exportgraphics(fig, [save_path '.png'], 'Resolution', 300);
    exportgraphics(fig, [save_path '.pdf'], 'ContentType', 'vector');
end

function plot_session_PSTHs(D, title_str, mua_windows, hemi_label, save_path)
    % One tile per session: mean +- SE (across events) PSTH per condition
    cond_colors = {[0.5 0.5 0.5], [231, 41, 138] / 256};
    n_sess = numel(D.sessions);
    n_col = ceil(sqrt(n_sess));
    n_row = ceil(n_sess / n_col);
    fig = figure('Name', title_str, 'Position', [50 50 260 * n_col 200 * n_row]);
    tiledlayout(n_row, n_col, 'TileSpacing', 'compact', 'Padding', 'compact');
    for si = 1:n_sess
        ax = nexttile; hold on;
        for c = 1:2
            if D.n_events(si, c) > 0
                plot_band(ax, D.t, D.psth(si, :, c), D.psth_se(si, :, c), cond_colors{c});
            end
        end
        xline(0, 'k--', 'HandleVisibility', 'off');
        shade_windows(ax, mua_windows);
        excl = '';
        if ~D.keep(si), excl = ' (excluded)'; end
        title(sprintf('Session %d: n = %d / %d%s', D.sessions(si), D.n_events(si, 1), D.n_events(si, 2), excl), ...
            'FontSize', 8, 'FontWeight', 'normal');
        set(ax, 'TickDir', 'out', 'Box', 'off', 'Color', 'none', 'FontSize', 8);
        xlim(D.t([1 end]));
        if si == 1
            legend(D.labels, 'Box', 'off', 'Location', 'best');
            ylabel([hemi_label ' (z)']);
            xlabel('Time from ripple onset (s)');
        end
    end
    sgtitle([title_str ': ' hemi_label ' PSTH per session (mean \pm SE across events)'], 'FontSize', 11);
    save_figure(fig, save_path);
end

function plot_mean_PSTH(R, S, analyses, titles, dom_windows, mua_windows, hemi_label, save_path)
    % Rows: dominance window; columns: analysis. Mean +- SE across sessions.
    cond_colors = {[0.5 0.5 0.5], [231, 41, 138] / 256};
    fig = figure('Name', [hemi_label ' PSTH'],'Position', [80 80 380 * numel(analyses) 340 * numel(dom_windows)]);
    tiledlayout(numel(dom_windows), numel(analyses), 'TileSpacing', 'compact');
    for w = 1:numel(dom_windows)
        for iA = 1:numel(analyses)
            D = R.(sprintf('%s_%s', analyses{iA}, dom_windows(w).name));
            ax = nexttile; hold on;
            n = sum(D.keep);
            for c = 1:2
                P = D.psth(D.keep, :, c);
                plot_band(ax, D.t, mean(P, 1), std(P, 0, 1) / sqrt(n), cond_colors{c});
            end
            xline(0, 'k--', 'HandleVisibility', 'off');
            shade_windows(ax, mua_windows);
            Sr = S(S.analysis == analyses{iA} & S.dominance_window == dom_windows(w).name, :);
            p_str = arrayfun(@(m) sprintf('%s: p = %.3g', mua_windows(m).label, ...
                Sr.p_signrank_two(Sr.MUA_window == mua_windows(m).name)), 1:numel(mua_windows), 'UniformOutput', false);
            title([{sprintf('%s, %s', titles{iA}, dom_windows(w).label), sprintf('n = %d sessions', n)}, p_str], ...
                'FontSize', 8, 'FontWeight', 'normal');
            set(ax, 'TickDir', 'out', 'Box', 'off', 'Color', 'none', 'FontSize', 10);
            xlim(D.t([1 end]));
            xlabel('Time from ripple onset (s)'); ylabel([hemi_label ' (z)']);
            legend(D.labels, 'Box', 'off', 'Location', 'northwest');
        end
    end
    sgtitle([hemi_label ' around ripple onset: mean \pm SE across sessions (signed-rank on session means)'], 'FontSize', 11);
    save_figure(fig, save_path);
end

function plot_window_means(R, S, analyses, titles, dom_windows, mua_windows, hemi_label, save_path)
    % Session-mean MUA per window, paired lines across conditions
    cond_colors = {[0.6 0.6 0.6], [231, 41, 138] / 256};
    n_col = numel(analyses) * numel(mua_windows);
    fig = figure('Name', [hemi_label ' window means'],'Position', [80 80 230 * n_col 330 * numel(dom_windows)]);
    tiledlayout(numel(dom_windows), n_col, 'TileSpacing', 'compact');
    for w = 1:numel(dom_windows)
        for iA = 1:numel(analyses)
            D = R.(sprintf('%s_%s', analyses{iA}, dom_windows(w).name));
            for m = 1:numel(mua_windows)
                ax = nexttile; hold on;
                Y = squeeze(D.window_MUA(D.keep, m, :));
                plot([1 2], Y', '-', 'Color', [0 0 0 0.25], 'LineWidth', 1);
                for c = 1:2
                    scatter(c * ones(size(Y, 1), 1), Y(:, c), 25, cond_colors{c}, 'filled', 'MarkerFaceAlpha', 0.7);
                    errorbar(c + 0.25 * (2 * c - 3), mean(Y(:, c)), std(Y(:, c)) / sqrt(size(Y, 1)), 'o', 'Color', 'k', ...
                        'MarkerFaceColor', cond_colors{c}, 'LineWidth', 1.5, 'CapSize', 6);
                end
                yline(0, 'k--');
                Sr = S(S.analysis == analyses{iA} & S.dominance_window == dom_windows(w).name & ...
                    S.MUA_window == mua_windows(m).name, :);
                title({titles{iA}, dom_windows(w).label, mua_windows(m).label, ...
                    sprintf('p = %.3g, n = %d', Sr.p_signrank_two, Sr.n_sessions)}, 'FontSize', 8, 'FontWeight', 'normal');
                set(ax, 'XTick', [1 2], 'XTickLabel', D.labels, 'TickDir', 'out', 'Box', 'off', 'Color', 'none', 'FontSize', 10);
                xlim([0.4 2.6]);
                ylabel([hemi_label ' (z, session mean)']);
            end
        end
    end
    sgtitle([hemi_label ' per session: Low vs High (Peak vs Trough), two-sided signed-rank'], 'FontSize', 11);
    save_figure(fig, save_path);
end
