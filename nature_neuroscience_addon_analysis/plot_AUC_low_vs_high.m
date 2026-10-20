function plot_AUC_low_vs_high(folder_name)
% Simplified view of reactivation_AUC_coherence_PRE output: for each analysis
% only the two contrasting conditions (ripple power and spindle Q1 vs Q4, SO
% peak vs trough, and dominant SO trough with non-dominant peak vs trough),
% for the pre (-0.2-0 s) and post (0-0.2 s) V1 windows.
% Reads the saved AUC_*.mat in nature_neuroscience_revision_analysis\<folder_name>
% and saves AUC_low_vs_high.fig/.png/.pdf there.
% Per window: Shuffle (mean, 95% CI), condition A (pale bar) and B (full bar)
% with 95% bootstrap CI. '*' : condition 95% CI does not overlap the shuffled
% 95% CI; 'n.s.' otherwise. Panel drawing: AUC_low_high_panel; colours:
% AUC_low_high_style (shared with plot_AUC_low_vs_high_controls).
    if exist('D:\corticohippocampal_replay', 'dir') > 0
        analysis_folder = 'D:\corticohippocampal_replay';
    elseif exist('P:\corticohippocampal_replay', 'dir') > 0
        analysis_folder = 'P:\corticohippocampal_replay';
    end
    in_dir = fullfile(analysis_folder, 'nature_neuroscience_revision_analysis', folder_name);

    style = AUC_low_high_style();
    sets = struct( ...
        'file',  {'AUC_ripple_percentile', 'AUC_spindle_percentile', 'AUC_SO_phase', 'AUC_SO_trough_bilateral'}, ...
        'style', {'ripple_power', 'spindle', 'SO_phase', 'SO_trough_uni_bi'}, ...
        'conds', {[1 4], [1 4], [1 2], [1 2]});
    y_lim = [-0.2 0.44];   % fixed so PRE and POST figures are comparable; top leaves room for labels
    windows = {'pre', 'post'};
    window_labels = {'V1 -0.2-0 s', 'V1 0-0.2 s'};

    fig = figure('Name', [folder_name ' low vs high'], 'Position', [200 200 1330 440]);
    tiledlayout(1, numel(sets), 'TileSpacing', 'compact');
    for iSet = 1:numel(sets)
        L = load(fullfile(in_dir, [sets(iSet).file '.mat']), 'res');
        st = style.(sets(iSet).style);
        R = struct('m', {}, 'ci', {}, 'ms', {}, 'cis', {});
        for w = 1:numel(windows)
            r = arrayfun(@(c) L.res(c).(windows{w}), sets(iSet).conds);
            R(w).m = [r.AUC_mean];
            R(w).ci = [r(1).AUC_CI(:) r(2).AUC_CI(:)];
            R(w).ms = [r.AUC_mean_shuffled];
            R(w).cis = [r(1).AUC_CI_shuffled(:) r(2).AUC_CI_shuffled(:)];
        end
        ax = nexttile;
        AUC_low_high_panel(ax, R, st.labels, st.colors, window_labels, y_lim);
        title(ax, st.title, 'FontWeight', 'normal');
    end
    sgtitle(strrep(folder_name, '_', ' '), 'FontSize', 13);
    savefig(fig, fullfile(in_dir, 'AUC_low_vs_high.fig'));
    exportgraphics(fig, fullfile(in_dir, 'AUC_low_vs_high.png'), 'Resolution', 300);
    exportgraphics(fig, fullfile(in_dir, 'AUC_low_vs_high.pdf'), 'ContentType', 'vector');
end
