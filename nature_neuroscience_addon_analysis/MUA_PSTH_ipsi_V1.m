%% Ipsilateral V1 MUA PSTH around ripple onset: Low vs High (Peak vs Trough), POST
% Simplified version of MUA_PSTH_low_vs_high without the track-bias based
% dominant hemisphere: each V1 hemisphere's MUA is split by its own
% (ipsilateral) spindle power / SO phase.
% Conditions (within session):
%   ripple_power : ripple power percentile <= 25 (Low) vs > 75 (High)
%   spindle      : ipsilateral spindle percentile <= 25 (Low) vs > 75 (High),
%                  re-ranked within session
%   SO_phase     : ipsilateral SO peak [-pi/2, pi/2] vs trough
% Per session: mean PSTH per condition for V1 L and V1 R, then the mean of the
% two hemispheres. Across sessions: mean +- SE PSTH and paired signed-rank test
% (High vs Low, Trough vs Peak) on the session mean MUA in -0.1-0 s and 0-0.1 s.
% MUA: z-scored, 10 ms bins, -2 to 2 s around ripple onset (ripple_MUA_V1_L/R).
% Output: nature_neuroscience_revision_analysis\MUA_PSTH_ipsi_V1
function MUA_PSTH_ipsi_V1(plot_range, out_dir)
    % plot_range: PSTH x-axis (s), default [-0.6 0.6]; out_dir: default
    % nature_neuroscience_revision_analysis\MUA_PSTH_ipsi_V1
    rev_dir = addon_setup();
    if nargin < 1 || isempty(plot_range)
        plot_range = [-0.6 0.6];
    end
    if nargin < 2 || isempty(out_dir)
        out_dir = fullfile(rev_dir, 'MUA_PSTH_ipsi_V1');
    end
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end

    load(fullfile(rev_dir, 'ripple_info.mat'), 'ripple_info');
    load(fullfile(rev_dir, 'ripple_MUA_V1_L.mat'), 'ripple_MUA_V1_L');
    load(fullfile(rev_dir, 'ripple_MUA_V1_R.mat'), 'ripple_MUA_V1_R');
    MUA = {ripple_MUA_V1_L, ripple_MUA_V1_R};   % hemisphere 1 = L, 2 = R
    clear ripple_MUA_V1_L ripple_MUA_V1_R
    n_ev = numel(ripple_info.session_id);
    assert(size(MUA{1}, 1) == n_ev && size(MUA{2}, 1) == n_ev, 'MUA rows do not match ripple_info');

    timebin = 0.01;
    mua_centers = -2 + timebin/2:timebin:2 - timebin/2;
    in_plot = mua_centers > plot_range(1) & mua_centers < plot_range(2);
    t = mua_centers(in_plot);
    mua_windows = struct('name', {'pre', 'post'}, 'label', {'-0.1-0 s', '0-0.1 s'}, 'range', {[-0.1 0], [0 0.1]});
    in_win = arrayfun(@(m) t > m.range(1) & t < m.range(2), mua_windows, 'UniformOutput', false);

    sess = ripple_info.session_id(:);
    sessions = unique(sess)';
    n_sess = numel(sessions);
    ripple_pct = within_session_percentile(ripple_info.ripple_power(:), sess);
    spindle_pct = within_session_percentile(ripple_info.spindle_percentile, sess);

    analyses = {'ripple_power', 'spindle', 'SO_phase'};
    titles = {'Ripple power', 'Spindle (ipsilateral)', 'SO phase (ipsilateral)'};

    R = struct();
    S = table();
    T = table();
    for iA = 1:numel(analyses)
        psth = NaN(n_sess, numel(t), 2, 2);   % session x time x condition x hemisphere
        n_c = zeros(n_sess, 2, 2);
        for h = 1:2
            [masks, labels] = condition_masks(analyses{iA}, ripple_info, ripple_pct, spindle_pct, h);
            valid = all(isfinite(MUA{h}(:, in_plot)), 2);
            for si = 1:n_sess
                for c = 1:2
                    i = masks{c} & valid & sess == sessions(si);
                    n_c(si, c, h) = sum(i);
                    if n_c(si, c, h) > 0
                        psth(si, :, c, h) = mean(double(MUA{h}(i, in_plot)), 1);
                    end
                end
            end
        end
        psth_sess = mean(psth, 4, 'omitnan');    % mean of the two hemispheres
        wm = NaN(n_sess, numel(mua_windows), 2);
        for m = 1:numel(mua_windows)
            wm(:, m, :) = mean(psth_sess(:, in_win{m}, :), 2);
        end
        R.(analyses{iA}) = struct('labels', {labels}, 'sessions', sessions, 't', t, 'psth', psth_sess, ...
            'psth_hemi', psth, 'n_events', n_c, 'window_MUA', wm);

        for m = 1:numel(mua_windows)
            a = wm(:, m, 1); b = wm(:, m, 2);
            ok = isfinite(a) & isfinite(b);
            S = [S; table(string(analyses{iA}), string(mua_windows(m).name), ...
                string(sprintf('%s vs %s', labels{2}, labels{1})), sum(ok), ...
                mean(a(ok)), std(a(ok)) / sqrt(sum(ok)), mean(b(ok)), std(b(ok)) / sqrt(sum(ok)), ...
                mean(b(ok) - a(ok)), signrank(b(ok), a(ok)), ...
                'VariableNames', {'analysis', 'MUA_window', 'comparison', 'n_sessions', ...
                'mean_A', 'SE_A', 'mean_B', 'SE_B', 'mean_diff', 'p_signrank'})]; %#ok<AGROW>
            for c = 1:2
                T = [T; table(sessions', repmat(string(analyses{iA}), n_sess, 1), repmat(string(labels{c}), n_sess, 1), ...
                    repmat(string(mua_windows(m).name), n_sess, 1), n_c(:, c, 1), n_c(:, c, 2), wm(:, m, c), ...
                    'VariableNames', {'session_id', 'analysis', 'condition', 'MUA_window', ...
                    'n_events_L', 'n_events_R', 'mean_MUA'})]; %#ok<AGROW>
            end
        end
    end
    disp(S);

    save(fullfile(out_dir, 'MUA_PSTH_ipsi_V1.mat'), 'R', 'S', 'T', 'mua_windows');
    writetable(S, fullfile(out_dir, 'MUA_PSTH_ipsi_V1_stats.csv'));
    writetable(T, fullfile(out_dir, 'MUA_PSTH_ipsi_V1_session_means.csv'));
    plot_summary(R, S, analyses, titles, mua_windows, fullfile(out_dir, 'MUA_PSTH_ipsi_V1'));
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

function [masks, labels] = condition_masks(analysis, info, ripple_pct, spindle_pct, h)
    % Logical event masks {A, B}; spindle / SO from hemisphere h
    switch analysis
        case 'ripple_power'
            masks = {ripple_pct <= 25, ripple_pct > 75};
            labels = {'Low', 'High'};
        case 'spindle'
            f = spindle_pct(:, h);
            masks = {f <= 25, f > 75};
            labels = {'Low', 'High'};
        case 'SO_phase'
            f = info.SO_phase(:, h);
            masks = {f >= -pi/2 & f <= pi/2, (f >= -pi & f < -pi/2) | (f > pi/2 & f <= pi)};
            labels = {'Peak', 'Trough'};
    end
end

function plot_summary(R, S, analyses, titles, mua_windows, save_path)
    % Columns: analysis. Row 1: mean +- SE PSTH across sessions; rows 2-3: session
    % mean MUA per window, paired across conditions
    style = AUC_low_high_style();
    n_row = 1 + numel(mua_windows);
    fig = figure('Name', 'Ipsilateral V1 MUA PSTH', 'Position', [80 80 330 * numel(analyses) 260 * n_row]);
    tiledlayout(n_row, numel(analyses), 'TileSpacing', 'compact');
    for iA = 1:numel(analyses)
        D = R.(analyses{iA});
        cond_colors = num2cell(style.(analyses{iA}).colors, 2);
        if ismember(analyses{iA}, {'ripple_power', 'spindle'})
            % Pale bar-plot Low colour is too faint for lines: light vs dark shade of the High colour
            base = cond_colors{2};
            cond_colors = {0.75 * base + 0.25 * [1 1 1], 0.85 * base};
        end
        ok = all(isfinite(D.window_MUA(:, :)), 2);
        n = sum(ok);

        ax = nexttile(iA); hold on;
        win_colors = {[0.3 0.6 1], [1 0.7 0.2]};
        for m = 1:numel(mua_windows)
            r = mua_windows(m).range;
            patch(ax, [r(1) r(2) r(2) r(1)], [-10 -10 10 10], win_colors{m}, 'FaceAlpha', 0.12, ...
                'EdgeColor', 'none', 'HandleVisibility', 'off');
        end
        h = gobjects(1, 2);
        y_lim = [Inf -Inf];
        for c = 1:2
            P = D.psth(ok, :, c);
            mu = mean(P, 1); se = std(P, 0, 1) / sqrt(n);
            fill(ax, [D.t fliplr(D.t)], [mu + se fliplr(mu - se)], cond_colors{c}, 'FaceAlpha', 0.3, ...
                'EdgeColor', 'none', 'HandleVisibility', 'off');
            h(c) = plot(ax, D.t, mu, '-', 'Color', cond_colors{c}, 'LineWidth', 1.5);
            y_lim = [min(y_lim(1), min(mu - se)) max(y_lim(2), max(mu + se))];
        end
        y_lim = y_lim + [-0.1 0.25] * diff(y_lim);   % extra top margin above the SE bands
        ylim(ax, y_lim);
        yt = ax.YTick;
        ax.YTick = yt(yt < y_lim(2) - 0.05 * diff(y_lim));   % no tick label against the title
        xline(ax, 0, 'k--', 'HandleVisibility', 'off');
        lg = legend(h, D.labels, 'Box', 'off', 'Location', 'northwest');
        lg.ItemTokenSize = [12 18];
        xlim(ax, D.t([1 end]));
        xlabel('Time relative to ripple onset (s)'); ylabel('V1 MUA (z)');
        title(sprintf('%s (n = %d sessions)', titles{iA}, n), 'FontSize', 10, 'FontWeight', 'bold');
        set(ax, 'TickDir', 'out', 'Box', 'off', 'Color', 'none', 'FontSize', 10);

        for m = 1:numel(mua_windows)
            ax = nexttile(m * numel(analyses) + iA); hold on;
            Y = squeeze(D.window_MUA(ok, m, :));
            plot([1 2], Y', '-', 'Color', [0 0 0 0.25], 'LineWidth', 1);
            for c = 1:2
                scatter(c * ones(n, 1), Y(:, c), 40, cond_colors{c}, 'filled', 'MarkerFaceAlpha', 0.6);
                errorbar(c + 0.25 * (2 * c - 3), mean(Y(:, c)), std(Y(:, c)) / sqrt(n), 'o', 'Color', 'k', ...
                    'MarkerFaceColor', cond_colors{c}, 'LineWidth', 1.5, 'CapSize', 6);
            end
            yline(0, 'k--');
            p = S.p_signrank(S.analysis == analyses{iA} & S.MUA_window == mua_windows(m).name);
            title(sprintf('MUA %s: p = %.3g', mua_windows(m).label, p), 'FontSize', 10, 'FontWeight', 'normal');
            set(ax, 'XTick', [1 2], 'XTickLabel', D.labels, 'TickDir', 'out', 'Box', 'off', 'Color', 'none', 'FontSize', 10);
            xlim([0.4 2.6]);
            ylabel('Mean V1 MUA (z)');
        end
    end
    sgtitle('Ipsilateral V1 MUA around ripple onset (mean of hemispheres; mean \pm SE across sessions; signed-rank)', 'FontSize', 10);
    savefig(fig, [save_path '.fig']);
    exportgraphics(fig, [save_path '.png'], 'Resolution', 300);
    exportgraphics(fig, [save_path '.pdf'], 'ContentType', 'vector');
end
