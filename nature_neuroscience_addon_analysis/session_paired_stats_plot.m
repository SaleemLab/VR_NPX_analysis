function S = session_paired_stats_plot(T, value_field, y_label, out_dir, prefix, test_method)
% Tests across sessions on a per-session value, and plots.
% T: table with session_id, sleep ("PRE"/"POST"), V1_window, analysis, condition
%    and the per-session value in T.(value_field).
% test_method (default 'signrank'):
%   'signrank'    - paired Wilcoxon signed-rank (sessions with both values)
%   'permutation' - exact paired sign-flip permutation test on the mean difference
%                   (all 2^n sign assignments; sessions with both values)
%   'ranksum'     - unpaired Wilcoxon rank-sum, all sessions available in each group
%                   (ignores pairing; PRE and POST share sessions so not independent)
% Comparisons (both two-sided and one-sided p saved):
%   within sleep (all sessions of that sleep): condition 1 vs 2 of each analysis,
%     one-sided condition 2 > condition 1 (High > Low, Trough > Peak,
%     Unilateral Trough > Bilateral Trough; see session_condition_labels)
%   PRE vs POST: each condition, one-sided POST > PRE for
%     High / Trough / Unilateral Trough / All, POST < PRE for Low / Peak / Bilateral Trough
% Plots show the one-sided p. Saves <prefix>stats.csv and <prefix><analysis>.png/.fig,
% and for paired tests <prefix>PREvsPOST_matched_<analysis>.png/.fig (PRE vs POST per
% condition, only sessions with both values, i.e. the points the test uses).
if nargin < 6
    test_method = 'signrank';
end
test_titles = containers.Map({'signrank', 'permutation', 'ranksum'}, ...
    {'paired signed-rank', 'paired permutation', 'unpaired rank-sum'});

sleeps = ["PRE", "POST"];
windows = unique(T.V1_window, 'stable');
analyses = unique(T.analysis, 'stable');
window_labels = containers.Map({'pre', 'post'}, {'V1 -0.2-0 s (pre)', 'V1 0-0.2 s (post)'});

S = table();
for iA = 1:numel(analyses)
    labels = session_condition_labels(T, analyses(iA));
    for w = 1:numel(windows)
        % Within sleep: condition 2 vs condition 1
        if numel(labels) == 2
            for iS = 1:2
                Y = session_condition_matrix(T, value_field, sleeps(iS), windows(w), analyses(iA), labels);
                S = [S; run_test(analyses(iA), windows(w), sprintf('%s: %s vs %s', sleeps(iS), labels(2), labels(1)), ...
                    Y(:, 2), Y(:, 1), 'right', sprintf('%s > %s', labels(2), labels(1)), test_method)]; %#ok<AGROW>
            end
        end
        % PRE vs POST per condition
        Ypre = session_condition_matrix(T, value_field, "PRE", windows(w), analyses(iA), labels);
        Ypost = session_condition_matrix(T, value_field, "POST", windows(w), analyses(iA), labels);
        for iCond = 1:numel(labels)
            if is_reference_condition(labels(iCond))
                tail = 'left'; direction = 'POST < PRE';
            else
                tail = 'right'; direction = 'POST > PRE';
            end
            S = [S; run_test(analyses(iA), windows(w), sprintf('%s: POST vs PRE', labels(iCond)), ...
                Ypost(:, iCond), Ypre(:, iCond), tail, direction, test_method)]; %#ok<AGROW>
        end
    end
end
disp(S);
writetable(S, fullfile(out_dir, [prefix 'stats.csv']));

% --- Plots: one figure per analysis, one tile per V1 window ---
sleep_colors = {[0.5 0.5 0.5], [231, 41, 138] / 256};
for iA = 1:numel(analyses)
    labels = session_condition_labels(T, analyses(iA));
    nCond = numel(labels); nW = numel(windows);
    fig = figure('Name', char(analyses(iA)), 'Position', [100 100 420*nW 560]);
    tiledlayout(1, nW, 'TileSpacing', 'compact');
    for w = 1:nW
        ax = nexttile; hold on;
        if nCond == 2
            xpos = {[1 2], [3.5 4.5]};
            xlabels = {sprintf('PRE %s', labels(1)), sprintf('PRE %s', labels(2)), ...
                       sprintf('POST %s', labels(1)), sprintf('POST %s', labels(2))};
        else
            xpos = {1, 2};
            xlabels = {'PRE', 'POST'};
        end
        Y = cell(1, 2);
        for iS = 1:2
            Y{iS} = session_condition_matrix(T, value_field, sleeps(iS), windows(w), analyses(iA), labels);
            if nCond == 2
                plot(xpos{iS}, Y{iS}', '-', 'Color', [0 0 0 0.25], 'LineWidth', 1);
            end
        end
        if nCond == 1
            plot([1 2], [Y{1} Y{2}]', '-', 'Color', [0 0 0 0.25], 'LineWidth', 1);
        end
        for iS = 1:2
            for iCond = 1:nCond
                y = Y{iS}(:, iCond); y = y(isfinite(y));
                x = xpos{iS}(iCond);
                scatter(x * ones(size(y)), y, 25, sleep_colors{iS}, 'filled', 'MarkerFaceAlpha', 0.7);
                if nCond == 2, dx = (2 * iCond - 3) * 0.3; else, dx = 0.3; end
                errorbar(x + dx, mean(y), std(y) / sqrt(numel(y)), 'o', 'Color', 'k', ...
                    'MarkerFaceColor', sleep_colors{iS}, 'LineWidth', 1.5, 'CapSize', 6);
            end
        end
        yline(0, 'k--');
        set(ax, 'XTick', [xpos{:}], 'XTickLabel', xlabels, 'XTickLabelRotation', 30);
        xlim([min([xpos{:}]) - 0.6, max([xpos{:}]) + 0.6]);
        ylabel(y_label);
        Sw = S(S.analysis == analyses(iA) & S.V1_window == windows(w), :);
        if strcmp(test_method, 'ranksum')
            p_lines = arrayfun(@(i) sprintf('%s (%s): p = %.3g, n = %d vs %d', Sw.comparison(i), Sw.one_sided_direction(i), ...
                Sw.p_one_sided(i), Sw.n_test(i), Sw.n_ref(i)), 1:height(Sw), 'UniformOutput', false);
        else
            p_lines = arrayfun(@(i) sprintf('%s (%s): p = %.3g, n = %d', Sw.comparison(i), Sw.one_sided_direction(i), ...
                Sw.p_one_sided(i), Sw.n_sessions(i)), 1:height(Sw), 'UniformOutput', false);
        end
        title([{window_labels(char(windows(w)))}, p_lines], 'FontSize', 8, 'FontWeight', 'normal');
        set(ax, 'TickDir', 'out', 'Box', 'off', 'Color', 'none', 'FontSize', 10);
    end
    linkaxes(findobj(fig, 'Type', 'axes'), 'y');
    sgtitle(sprintf('%s (%s, one-sided)', strrep(analyses(iA), '_', ' '), test_titles(test_method)), 'FontSize', 12);
    savefig(fig, fullfile(out_dir, [prefix char(analyses(iA)) '.fig']));
    exportgraphics(fig, fullfile(out_dir, [prefix char(analyses(iA)) '.png']), 'Resolution', 300);
    exportgraphics(fig, fullfile(out_dir, [prefix char(analyses(iA)) '.pdf']), 'ContentType', 'vector');
end

% --- Matched PRE vs POST plots (paired tests only): sessions with both values ---
if ~strcmp(test_method, 'ranksum')
    for iA = 1:numel(analyses)
        labels = session_condition_labels(T, analyses(iA));
        nCond = numel(labels); nW = numel(windows);
        fig = figure('Name', [char(analyses(iA)) ' matched'], 'Position', [100 100 380*nW*max(nCond, 1.3) 560]);
        tiledlayout(1, nW, 'TileSpacing', 'compact');
        for w = 1:nW
            ax = nexttile; hold on;
            Ypre = session_condition_matrix(T, value_field, "PRE", windows(w), analyses(iA), labels);
            Ypost = session_condition_matrix(T, value_field, "POST", windows(w), analyses(iA), labels);
            xt = []; xl = {}; p_lines = {};
            for iCond = 1:nCond
                xpos = 3 * (iCond - 1) + [1 2];
                n = plot_matched_pre_post(ax, Ypre(:, iCond), Ypost(:, iCond), xpos);
                xt = [xt, xpos]; %#ok<AGROW>
                xl = [xl, {sprintf('PRE %s', labels(iCond)), sprintf('POST %s', labels(iCond))}]; %#ok<AGROW>
                st = S(S.analysis == analyses(iA) & S.V1_window == windows(w) & ...
                    S.comparison == sprintf('%s: POST vs PRE', labels(iCond)), :);
                p_lines{end+1} = sprintf('%s (%s): p = %.3g (two-sided %.3g), n = %d', st.comparison, ...
                    st.one_sided_direction, st.p_one_sided, st.p_two_sided, n); %#ok<AGROW>
            end
            yline(0, 'k--');
            set(ax, 'XTick', xt, 'XTickLabel', xl, 'XTickLabelRotation', 30);
            xlim([0.3, max(xt) + 0.7]);
            ylabel(y_label);
            title([{window_labels(char(windows(w)))}, p_lines], 'FontSize', 8, 'FontWeight', 'normal');
            set(ax, 'TickDir', 'out', 'Box', 'off', 'Color', 'none', 'FontSize', 10);
        end
        linkaxes(findobj(fig, 'Type', 'axes'), 'y');
        sgtitle(sprintf('%s: PRE vs POST, matched sessions (%s, one-sided)', strrep(analyses(iA), '_', ' '), ...
            test_titles(test_method)), 'FontSize', 12);
        savefig(fig, fullfile(out_dir, [prefix 'PREvsPOST_matched_' char(analyses(iA)) '.fig']));
        exportgraphics(fig, fullfile(out_dir, [prefix 'PREvsPOST_matched_' char(analyses(iA)) '.png']), 'Resolution', 300);
        exportgraphics(fig, fullfile(out_dir, [prefix 'PREvsPOST_matched_' char(analyses(iA)) '.pdf']), 'ContentType', 'vector');
    end
end
end

function row = run_test(analysis, window, comparison, y, x, tail, direction, test_method)
% Test y vs x; one-sided tail on y - x (y > x for 'right')
p2 = NaN; p1 = NaN;
if strcmp(test_method, 'ranksum')
    % Unpaired: every session with a value in each group
    xs = x(isfinite(x)); ys = y(isfinite(y));
    n_ref = numel(xs); n_test = numel(ys);
    n_paired = sum(isfinite(x) & isfinite(y));
    n_gt = NaN;
    if n_ref >= 3 && n_test >= 3
        p2 = ranksum(ys, xs);
        p1 = ranksum(ys, xs, 'tail', tail);
    end
else
    ok = isfinite(x) & isfinite(y);
    xs = x(ok); ys = y(ok);
    n_ref = sum(ok); n_test = sum(ok); n_paired = sum(ok);
    n_gt = sum(ys > xs);
    if n_paired >= 3
        switch test_method
            case 'signrank'
                p2 = signrank(ys, xs);
                p1 = signrank(ys, xs, 'tail', tail);
            case 'permutation'
                [p2, p1] = paired_perm_test(ys - xs, tail);
        end
    end
end
row = table(string(analysis), string(window), string(comparison), n_paired, mean(xs), mean(ys), ...
    n_gt, p2, p1, string(direction), n_ref, n_test, string(test_method), ...
    'VariableNames', {'analysis', 'V1_window', 'comparison', 'n_sessions', 'mean_ref', 'mean_test', ...
    'n_test_gt_ref', 'p_two_sided', 'p_one_sided', 'one_sided_direction', 'n_ref', 'n_test', 'test'});
end
