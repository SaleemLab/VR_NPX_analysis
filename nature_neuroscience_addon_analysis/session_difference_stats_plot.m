function S = session_difference_stats_plot(T, value_field, y_label, out_dir, prefix)
% Condition x sleep interaction across sessions (main PRE vs POST test).
% Per session and sleep: Delta = condition 2 - condition 1 (High - Low,
% Trough - Peak, Unilateral Trough - Bilateral Trough; see session_condition_labels).
% POST Delta vs PRE Delta, one-sided POST Delta > PRE Delta (two-sided also saved):
%   paired (sessions with both PRE and POST Delta): Wilcoxon signed-rank and exact
%     sign-flip permutation test on the mean difference
%   unpaired (all sessions with a Delta in each sleep): Wilcoxon rank-sum
% T: table with session_id, sleep ("PRE"/"POST"), V1_window, analysis, condition
%    and the per-session value in T.(value_field). Single-condition analyses are skipped.
% Plots show the one-sided paired permutation and signed-rank p.
% Saves <prefix>stats.csv, <prefix><analysis>.png/.fig (all sessions) and
% <prefix>matched_<analysis>.png/.fig (only sessions with both PRE and POST Delta).

windows = unique(T.V1_window, 'stable');
analyses = unique(T.analysis, 'stable');
window_labels = containers.Map({'pre', 'post'}, {'V1 -0.2-0 s (pre)', 'V1 0-0.2 s (post)'});
sleep_colors = {[0.5 0.5 0.5], [231, 41, 138] / 256};

S = table();
for iA = 1:numel(analyses)
    labels = session_condition_labels(T, analyses(iA));
    if numel(labels) ~= 2
        continue
    end
    delta_label = sprintf('%s - %s', labels(2), labels(1));
    nW = numel(windows);
    fig = figure('Name', char(analyses(iA)), 'Position', [100 100 380*nW 560]);
    tl = tiledlayout(fig, 1, nW, 'TileSpacing', 'compact');
    fig_m = figure('Name', [char(analyses(iA)) ' matched'], 'Position', [100 100 380*nW 560]);
    tl_m = tiledlayout(fig_m, 1, nW, 'TileSpacing', 'compact');
    for w = 1:nW
        Ypre = session_condition_matrix(T, value_field, "PRE", windows(w), analyses(iA), labels);
        Ypost = session_condition_matrix(T, value_field, "POST", windows(w), analyses(iA), labels);
        d_pre = Ypre(:, 2) - Ypre(:, 1);
        d_post = Ypost(:, 2) - Ypost(:, 1);
        row = run_difference_tests(analyses(iA), windows(w), delta_label, d_pre, d_post);
        S = [S; row]; %#ok<AGROW>

        % --- Matched plot: only sessions with both PRE and POST Delta (paired tests) ---
        ax = nexttile(tl_m);
        n = plot_matched_pre_post(ax, d_pre, d_post, [1 2]);
        yline(ax, 0, 'k--');
        set(ax, 'XTick', [1 2], 'XTickLabel', {'PRE', 'POST'}, 'TickDir', 'out', 'Box', 'off', 'Color', 'none', 'FontSize', 10);
        xlim(ax, [0.4 2.6]);
        ylabel(ax, sprintf('%s: %s', y_label, delta_label));
        title(ax, {window_labels(char(windows(w))), ...
            sprintf('POST > PRE (n = %d): perm p = %.3g (two-sided %.3g)', n, row.p_perm_one_sided, row.p_perm_two_sided), ...
            sprintf('signrank p = %.3g (two-sided %.3g)', row.p_signrank_one_sided, row.p_signrank_two_sided)}, ...
            'FontSize', 8, 'FontWeight', 'normal');

        % --- Plot: PRE Delta vs POST Delta, all sessions, paired lines for sessions with both ---
        ax = nexttile(tl); hold(ax, 'on');
        ok = isfinite(d_pre) & isfinite(d_post);
        plot(ax, [1 2], [d_pre(ok) d_post(ok)]', '-', 'Color', [0 0 0 0.25], 'LineWidth', 1);
        d_all = {d_pre(isfinite(d_pre)), d_post(isfinite(d_post))};
        for iS = 1:2
            y = d_all{iS};
            scatter(ax, iS * ones(size(y)), y, 25, sleep_colors{iS}, 'filled', 'MarkerFaceAlpha', 0.7);
            errorbar(ax, iS + 0.3, mean(y), std(y) / sqrt(numel(y)), 'o', 'Color', 'k', ...
                'MarkerFaceColor', sleep_colors{iS}, 'LineWidth', 1.5, 'CapSize', 6);
        end
        yline(ax, 0, 'k--');
        set(ax, 'XTick', [1 2], 'XTickLabel', {'PRE', 'POST'});
        xlim(ax, [0.4 2.6]);
        ylabel(ax, sprintf('%s: %s', y_label, delta_label));
        title(ax, {window_labels(char(windows(w))), ...
            sprintf('POST > PRE (paired, n = %d): perm p = %.3g, signrank p = %.3g', ...
                row.n_paired, row.p_perm_one_sided, row.p_signrank_one_sided), ...
            sprintf('rank-sum (n = %d vs %d): p = %.3g', row.n_post_all, row.n_pre_all, row.p_ranksum_one_sided)}, ...
            'FontSize', 8, 'FontWeight', 'normal');
        set(ax, 'TickDir', 'out', 'Box', 'off', 'Color', 'none', 'FontSize', 10);
    end
    linkaxes(findobj(fig, 'Type', 'axes'), 'y');
    sgtitle(fig, sprintf('%s: (%s) POST vs PRE (one-sided)', strrep(analyses(iA), '_', ' '), delta_label), 'FontSize', 12);
    savefig(fig, fullfile(out_dir, [prefix char(analyses(iA)) '.fig']));
    exportgraphics(fig, fullfile(out_dir, [prefix char(analyses(iA)) '.png']), 'Resolution', 300);
    exportgraphics(fig, fullfile(out_dir, [prefix char(analyses(iA)) '.pdf']), 'ContentType', 'vector');
    linkaxes(findobj(fig_m, 'Type', 'axes'), 'y');
    sgtitle(fig_m, sprintf('%s: (%s) POST vs PRE, matched sessions (one-sided)', strrep(analyses(iA), '_', ' '), delta_label), ...
        'FontSize', 12);
    savefig(fig_m, fullfile(out_dir, [prefix 'matched_' char(analyses(iA)) '.fig']));
    exportgraphics(fig_m, fullfile(out_dir, [prefix 'matched_' char(analyses(iA)) '.png']), 'Resolution', 300);
    exportgraphics(fig_m, fullfile(out_dir, [prefix 'matched_' char(analyses(iA)) '.pdf']), 'ContentType', 'vector');
end
disp(S);
writetable(S, fullfile(out_dir, [prefix 'stats.csv']));
end

function row = run_difference_tests(analysis, window, delta_label, d_pre, d_post)
% Paired (signed-rank, permutation) and unpaired (rank-sum) POST Delta vs PRE Delta
ok = isfinite(d_pre) & isfinite(d_post);
n = sum(ok);
ps2 = NaN; ps1 = NaN; pp2 = NaN; pp1 = NaN; pr2 = NaN; pr1 = NaN;
if n >= 3
    ps2 = signrank(d_post(ok), d_pre(ok));
    ps1 = signrank(d_post(ok), d_pre(ok), 'tail', 'right');
    [pp2, pp1] = paired_perm_test(d_post(ok) - d_pre(ok), 'right');
end
pre_all = d_pre(isfinite(d_pre)); post_all = d_post(isfinite(d_post));
if numel(pre_all) >= 3 && numel(post_all) >= 3
    pr2 = ranksum(post_all, pre_all);
    pr1 = ranksum(post_all, pre_all, 'tail', 'right');
end
row = table(string(analysis), string(window), string(delta_label), n, mean(d_pre(ok)), mean(d_post(ok)), ...
    sum(d_post(ok) > d_pre(ok)), ps2, ps1, pp2, pp1, ...
    numel(pre_all), numel(post_all), mean(pre_all), mean(post_all), pr2, pr1, ...
    'VariableNames', {'analysis', 'V1_window', 'delta', 'n_paired', 'mean_delta_PRE', 'mean_delta_POST', ...
    'n_POST_gt_PRE', 'p_signrank_two_sided', 'p_signrank_one_sided', 'p_perm_two_sided', 'p_perm_one_sided', ...
    'n_pre_all', 'n_post_all', 'mean_delta_PRE_all', 'mean_delta_POST_all', 'p_ranksum_two_sided', 'p_ranksum_one_sided'});
end
