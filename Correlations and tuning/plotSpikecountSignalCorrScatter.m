% This function creates a scatter plot comparing the mean trial spike counts of two neuronal clusters in response to different
% stimuli. Each point represents a stimulus orientation condition, labeled accordingly.
% Also fits and plots a regression line and calculates the Pearson correlation coefficient.
% Figure is auto-saved in the cd.
% ERB 2025
%
% Inputs:
%   x_vals            - Vector of mean trial spike counts for cluster 1 (z-scored).
%   y_vals            - Vector of mean trial spike counts for cluster 2 (z-scored).
%   cluster1_id       - Numeric identifier for the first cluster.
%   cluster2_id       - Numeric identifier for the second cluster.
%   r_val             - Correlation coefficient value from a matrix or prior calculation (for title reference).
%   ori_labels        - Cell array of strings labeling each data point's stimulus orientation.
%   subject_number    - Identifier string for the subject (animal/session).
%   Stimulus_type     - String describing the stimulus type used in the experiment.
%   depth_for_analysis- String indicating the depth or brain region analyzed.
%   temporal_structure- String describing the temporal structure of the stimulus presentation.
%   stim_window       - Two-element vector [start_time, end_time] in seconds relative to stimulus onset.
%   gOSI_threshold    - Numeric value indicating the threshold for orientation selectivity index used.
%   z_score_period    - String describing the baseline period used for z-scoring spike counts (entire_session is appropriate here).

function plotSpikecountSignalCorrScatter(x_vals, y_vals, cluster1_id, cluster2_id, r_val, ori_labels,...
    subject_number, Stimulus_type, depth_for_analysis, temporal_structure, stim_window, gOSI_threshold, z_score_period)
    
    [R, P] = corr(x_vals(:), y_vals(:), 'rows', 'complete');

    fig=figure;
    scatter(x_vals, y_vals, 100, 'filled');
    hold on;

    % Fit line
    coeffs = polyfit(x_vals, y_vals, 1);
    x_fit = linspace(min(x_vals), max(x_vals), 100);
    y_fit = polyval(coeffs, x_fit);
    plot(x_fit, y_fit, 'r--', 'LineWidth', 2);


    xlabel(sprintf('Cluster %d mean Trial Spikecounts (%.2f–%.2fs after stim onset) Z-scored over %s', cluster1_id, stim_window(1), stim_window(2), z_score_period), 'Interpreter', 'none');
    ylabel(sprintf('Cluster %d mean Trial Spikecounts (%.2f–%.2fs after stim onset) Z-scored over %s', cluster2_id, stim_window(1), stim_window(2), z_score_period), 'Interpreter', 'none');
    title(sprintf('Signal corr (from this plot) r = %.2f, p = %.3g\n (r = %.2f matrix plot)', R, P, r_val));
    
    axis square;
    grid on;
    
    for i = 1:length(x_vals)
        text(x_vals(i) + 0.01*range(x_vals), y_vals(i) + 0.01*range(y_vals), ori_labels{i}, 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'black');
    end
    
    min_margin = 0.01;
    x_margin = max(0.05 * range(x_vals), min_margin);
    y_margin = max(0.05 * range(y_vals), min_margin);
    xlim([min(x_vals) - x_margin, max(x_vals) + x_margin]);
    ylim([min(y_vals) - y_margin, max(y_vals) + y_margin]);

    hold off;


    % Save figure
    fig_filename = sprintf('%s_%s_%s_SignalScatter_%dvs%d - %s (%.2f–%.2fs after stim onset), gOSI over %0.2f.fig', ...
        subject_number, Stimulus_type, depth_for_analysis, cluster1_id, cluster2_id, temporal_structure, stim_window(1), stim_window(2), ...
        gOSI_threshold);
    savefig(fig, fullfile(pwd, fig_filename));
    
end