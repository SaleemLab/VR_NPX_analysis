% plotTimebinnedSignalCorrScatter
% This function generates a scatter plot showing the correlation between timebinned spiking responses from two clusters across 
% multiple stimulus orientations. Each orientation is plotted in a different color, and a best-fit regression line is overlaid. 
% Figure is auto-saved in the cd.
% ERB 2025
% INPUTS:
%   x_vals            - Vector of timebinned responses for cluster1 
%   y_vals            - Vector of timebinned responses for cluster2 
%   cluster1_id       - Cluster ID for x_vals
%   cluster2_id       - Cluster ID for y_vals
%   r_val             - Precomputed signal correlation between the clusters (from matrix)
%   ordered_oris      - Vector of orientations of presented stimlui (e.g. [0, 45, 90, ...])
%   subject_number    - String/ID of the subject or experiment
%   Stimulus_type     - Stimulus type (e.g. 'ABCD')
%   depth_for_analysis- Depth label or range of cortical layer (e.g. 'L4')
%   temporal_structure- Whether the analysis uses a simple spikecount response figure or timebinned spikes
%   stim_window       - Time window [start, end] after stimulus onset (in seconds)
%   gOSI_threshold    - gOSI cutoff used in filtering clusters
%   z_score_period    - Label for z-scoring baseline (entire_session is appropriate here)

function plotTimebinnedSignalCorrScatter(x_vals, y_vals, cluster1_id, cluster2_id, r_val, ori_labels, ordered_oris,...
    subject_number, Stimulus_type, depth_for_analysis, temporal_structure, stim_window, gOSI_threshold, z_score_period)
    % Compute Pearson correlation across all timebinned points
    [R, P] = corr(x_vals(:), y_vals(:), 'rows', 'complete');

    n_timebins = length(x_vals) / length(ordered_oris);

    fig = figure;
    cmap = lines(length(ordered_oris));  % distinct colors for each orientation
    hold on;

    jitter_amount = 0.01 * range(x_vals);
    % Plot timebin points for each orientation with a color
    for i = 1:length(ordered_oris)
        idx_range = (i-1)*n_timebins + (1:n_timebins);
        
        x_jittered = x_vals(idx_range) + jitter_amount * randn(size(idx_range));
        y_jittered = y_vals(idx_range) + jitter_amount * randn(size(idx_range));
        
        scatter(x_jittered, y_jittered, 100, ...
            'MarkerFaceColor', cmap(i,:), 'MarkerEdgeColor', 'k', ...
            'DisplayName', sprintf('Ori %s', ori_labels{i}));
    end

    % Fit overall regression line
    coeffs = polyfit(x_vals, y_vals, 1);
    x_fit = linspace(min(x_vals), max(x_vals), 100);
    y_fit = polyval(coeffs, x_fit);
    plot(x_fit, y_fit, 'r--', 'LineWidth', 2, 'HandleVisibility', 'off');

    xlabel(sprintf('Cluster %d mean timebinned response (%.2f–%.2fs after stim onset) Z-scored over %s', cluster1_id, stim_window(1), stim_window(2), z_score_period), 'Interpreter', 'none');
    ylabel(sprintf('Cluster %d mean timebinned response (%.2f–%.2fs after stim onset) Z-scored over %s', cluster2_id, stim_window(1), stim_window(2), z_score_period), 'Interpreter', 'none');
    title(sprintf('Signal corr (this plot) r = %.2f, p = %.3g\n (matrix r = %.2f)', R, P, r_val));

    legend('show');
    axis square;
    grid on;

    min_margin = 0.01;
    x_margin = max(0.05 * range(x_vals), min_margin);
    y_margin = max(0.05 * range(y_vals), min_margin);
    xlim([min(x_vals) - x_margin, max(x_vals) + x_margin]);
    ylim([min(y_vals) - y_margin, max(y_vals) + y_margin]);

    % Save figure
    fig_filename = sprintf('%s_%s_%s_TimebinnedSignalScatter_%dvs%d - %s (%.2f–%.2fs after stim onset), gOSI over %.2f.fig', ...
        subject_number, Stimulus_type, depth_for_analysis, cluster1_id, cluster2_id, ...
        temporal_structure, stim_window(1), stim_window(2), gOSI_threshold);
    savefig(fig, fullfile(pwd, fig_filename));

end