function S = session_vs_shuffle_stats_plot(T, actual_field, ref_field, ref_label, out_dir, prefix)
% Paired Wilcoxon signed-rank across sessions: per-session actual AUC vs a
% per-session shuffle reference (e.g. shuffle mean, 95th prctile, CI upper bound).
% Saves two-sided and one-sided (actual > reference) p; plots show one-sided.
% T: table with session_id, sleep, V1_window, analysis, condition, T.(actual_field), T.(ref_field).
% Saves <prefix>stats.csv and <prefix><analysis>.png/.fig in out_dir.

S = table();
groups = unique(T(:, {'sleep', 'V1_window', 'analysis', 'condition'}), 'stable');
for g = 1:height(groups)
    Tg = T(T.sleep == groups.sleep(g) & T.V1_window == groups.V1_window(g) & ...
           T.analysis == groups.analysis(g) & T.condition == groups.condition(g), :);
    a = Tg.(actual_field); r = Tg.(ref_field);
    ok = isfinite(a) & isfinite(r);
    p2 = NaN; p1 = NaN;
    if sum(ok) >= 3
        p2 = signrank(a(ok), r(ok));
        p1 = signrank(a(ok), r(ok), 'tail', 'right');
    end
    S = [S; table(groups.sleep(g), groups.V1_window(g), groups.analysis(g), groups.condition(g), sum(ok), ...
        mean(a(ok)), mean(r(ok)), sum(a(ok) > r(ok)), p2, p1, ...
        'VariableNames', {'sleep', 'V1_window', 'analysis', 'condition', 'n_sessions', 'mean_actual', ...
        'mean_reference', 'n_sessions_above', 'p_two_sided', 'p_one_sided_actual_gt_ref'})]; %#ok<AGROW>
end
disp(S);
writetable(S, fullfile(out_dir, [prefix 'stats.csv']));

window_labels = containers.Map({'pre', 'post'}, {'V1 -0.2-0 s (pre)', 'V1 0-0.2 s (post)'});
analyses = unique(T.analysis, 'stable');
windows = unique(T.V1_window, 'stable');
sleeps = ["PRE", "POST"];
sleep_colors = {[0.5 0.5 0.5], [231, 41, 138] / 256};
for iA = 1:numel(analyses)
    Ta = T(T.analysis == analyses(iA), :);
    labels = unique(Ta.condition, 'stable');
    nCond = numel(labels); nW = numel(windows);
    fig = figure('Name', char(analyses(iA)), 'Position', [100 100 330*nW*nCond 420*2]);
    tiledlayout(2, nW, 'TileSpacing', 'compact');
    for iS = 1:2
        for w = 1:nW
            ax = nexttile; hold on;
            xt = []; xl = {}; p_lines = {};
            for iCond = 1:nCond
                Tc = Ta(Ta.sleep == sleeps(iS) & Ta.V1_window == windows(w) & Ta.condition == labels(iCond), :);
                x0 = (iCond - 1) * 2.5;
                Y = [Tc.(ref_field) Tc.(actual_field)];
                plot(x0 + [1 2], Y', '-', 'Color', [0 0 0 0.25], 'LineWidth', 1);
                scatter(x0 + ones(height(Tc), 1), Y(:, 1), 25, [0.6 0.6 0.6], 'filled', 'MarkerFaceAlpha', 0.7);
                scatter(x0 + 2 * ones(height(Tc), 1), Y(:, 2), 25, sleep_colors{iS}, 'filled', 'MarkerFaceAlpha', 0.7);
                for c = 1:2
                    y = Y(isfinite(Y(:, c)), c);
                    errorbar(x0 + c + (2*c - 3) * 0.3, mean(y), std(y) / sqrt(numel(y)), 'o', 'Color', 'k', ...
                        'MarkerFaceColor', 'k', 'LineWidth', 1.5, 'CapSize', 6);
                end
                xt = [xt, x0 + [1 2]]; %#ok<AGROW>
                xl = [xl, {ref_label, char(labels(iCond))}]; %#ok<AGROW>
                st = S(S.sleep == sleeps(iS) & S.V1_window == windows(w) & S.analysis == analyses(iA) & S.condition == labels(iCond), :);
                p_lines{end+1} = sprintf('%s > %s: p = %.3g (%d/%d above)', labels(iCond), ref_label, ...
                    st.p_one_sided_actual_gt_ref, st.n_sessions_above, st.n_sessions); %#ok<AGROW>
            end
            yline(0, 'k--');
            set(ax, 'XTick', xt, 'XTickLabel', xl, 'XTickLabelRotation', 30);
            xlim([0.4, max(xt) + 0.6]);
            ylabel('Reactivation coherence AUC');
            title([{sprintf('%s, %s', sleeps(iS), window_labels(char(windows(w))))}, p_lines], 'FontSize', 9, 'FontWeight', 'normal');
            set(ax, 'TickDir', 'out', 'Box', 'off', 'Color', 'none', 'FontSize', 10);
        end
    end
    linkaxes(findobj(fig, 'Type', 'axes'), 'y');
    sgtitle(sprintf('Session AUC vs %s: %s (paired signed-rank, one-sided)', ref_label, strrep(analyses(iA), '_', ' ')), 'FontSize', 12);
    savefig(fig, fullfile(out_dir, [prefix char(analyses(iA)) '.fig']));
    exportgraphics(fig, fullfile(out_dir, [prefix char(analyses(iA)) '.png']), 'Resolution', 300);
end
end
