function plot_AUC_low_vs_high(folder_name)
% Simplified view of reactivation_AUC_coherence_PRE output: for each analysis
% only the two contrasting conditions (ripple power and spindle Q1 vs Q4, SO
% peak vs trough, and dominant SO trough with non-dominant peak vs trough),
% for the pre (-0.2-0 s) and post (0-0.2 s) V1 windows.
% Reads the saved AUC_*.mat in nature_neuroscience_revision_analysis\<folder_name>
% and saves AUC_low_vs_high.fig/.png there.
% Per window: Shuffle (mean, 95% CI), condition A (pale bar) and B (full bar)
% with 95% bootstrap CI. '*' : condition 95% CI does not overlap the shuffled
% 95% CI; 'n.s.' otherwise.
    if exist('D:\corticohippocampal_replay', 'dir') > 0
        analysis_folder = 'D:\corticohippocampal_replay';
    elseif exist('P:\corticohippocampal_replay', 'dir') > 0
        analysis_folder = 'P:\corticohippocampal_replay';
    end
    in_dir = fullfile(analysis_folder, 'nature_neuroscience_revision_analysis', folder_name);

    green = [76 187 135] / 255; purple = [140 105 185] / 255;
    pale = @(c) 0.35 * c + 0.65 * [1 1 1];   % pale version for the Low condition
    sets = struct( ...
        'file',   {'AUC_ripple_percentile', 'AUC_spindle_percentile', 'AUC_SO_phase', 'AUC_SO_trough_bilateral'}, ...
        'title',  {'Ripple power (HPC->V1)', 'Spindle (V1->HPC)', 'SO phase (V1->HPC)', {'Dominant SO trough;', 'non-dominant SO (V1->HPC)'}}, ...
        'conds',  {[1 4], [1 4], [1 2], [1 2]}, ...
        'labels', {{'Low', 'High'}, {'Low', 'High'}, {'Peak', 'Trough'}, {'Peak', 'Trough'}}, ...
        'colors', {[pale(green); green], [pale(purple); purple], ...
                   [80 175 220; 230 100 165] / 255, [245 190 215; 200 40 120] / 255});
    y_lim = [-0.2 0.44];   % fixed so PRE and POST figures are comparable; top leaves room for labels
    windows = {'pre', 'post'};
    window_labels = {'V1 -0.2-0 s', 'V1 0-0.2 s'};
    group_x = [1 2 3];     % Shuffle, A, B within a window
    group_gap = 4;

    fig = figure('Name', [folder_name ' low vs high'], 'Position', [200 200 1330 440]);
    tiledlayout(1, numel(sets), 'TileSpacing', 'compact');
    for iSet = 1:numel(sets)
        L = load(fullfile(in_dir, [sets(iSet).file '.mat']), 'res');
        bar_cols = sets(iSet).colors;
        ax = nexttile; hold on;
        xt = []; xtl = {};
        for w = 1:numel(windows)
            x = group_x + (w - 1) * group_gap;
            r = arrayfun(@(c) L.res(c).(windows{w}), sets(iSet).conds);
            % Shuffle: pooled over the two conditions' shuffled CIs
            ci_s = [min(arrayfun(@(q) q.AUC_CI_shuffled(1), r)) max(arrayfun(@(q) q.AUC_CI_shuffled(2), r))];
            m_s = mean([r.AUC_mean_shuffled]);
            errorbar(x(1), m_s, m_s - ci_s(1), ci_s(2) - m_s, 'k', 'LineWidth', 1.5, 'CapSize', 8);
            for k = 1:2
                bar(x(k + 1), r(k).AUC_mean, 0.7, 'FaceColor', bar_cols(k, :), 'EdgeColor', 'none');
                errorbar(x(k + 1), r(k).AUC_mean, r(k).AUC_mean - r(k).AUC_CI(1), r(k).AUC_CI(2) - r(k).AUC_mean, ...
                    'k', 'LineWidth', 1.5, 'CapSize', 8);
                if r(k).AUC_CI(1) > ci_s(2) || r(k).AUC_CI(2) < ci_s(1)
                    text(x(k + 1), r(k).AUC_CI(2) + 0.008, '*', 'HorizontalAlignment', 'center', ...
                        'VerticalAlignment', 'bottom', 'FontSize', 16);
                else
                    text(x(k + 1), r(k).AUC_CI(2) + 0.008, 'n.s.', 'HorizontalAlignment', 'center', ...
                        'VerticalAlignment', 'bottom', 'FontSize', 11);
                end
            end
            text(mean(x), y_lim(2), window_labels{w}, 'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'top', 'FontSize', 11);
            xt = [xt x]; %#ok<AGROW>
            xtl = [xtl {'Shuffle'} sets(iSet).labels]; %#ok<AGROW>
        end
        yline(0, 'k--');
        set(ax, 'XTick', xt, 'XTickLabel', xtl, 'XTickLabelRotation', 45);
        xlim([0.3, xt(end) + 0.7]);
        ylim(y_lim);
        yticks(-0.2:0.1:0.3);
        ylabel('Reactivation coherence AUC');
        title(sets(iSet).title, 'FontWeight', 'normal');
        set(ax, 'TickDir', 'out', 'Box', 'off', 'Color', 'none', 'FontSize', 12);
    end
    sgtitle(strrep(folder_name, '_', ' '), 'FontSize', 13);
    savefig(fig, fullfile(in_dir, 'AUC_low_vs_high.fig'));
    exportgraphics(fig, fullfile(in_dir, 'AUC_low_vs_high.png'), 'Resolution', 300);
end
