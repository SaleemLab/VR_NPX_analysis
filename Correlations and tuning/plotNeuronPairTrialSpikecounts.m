function plotNeuronPairTrialSpikecounts(responses_this_ori, neuron_pair, sorted_cluster_ids, ori_letter, znoise_corr_value, ori_angle, subject_number, Stimulus_type, depth_for_analysis, temporal_structure, stim_window, gOSI_threshold)
% plotNeuronPairTrialSpikecounts  Create a scatter plot of spike counts for a neuron pair
% This function plots the trial-by-trial spike counts of two neurons for a given stimulus orientation, and fits a linear regression line
% to visualize their correlation.
% Figure is saved in the cd
% ERB 2025
% Inputs:
%   responses_this_ori   - [nTrials x nNeurons] matrix of z-scored spike counts for the current orientation
%   neuron_pair           - [1 x 2] vector containing the indices of the two neurons (columns in responses_this_ori) to be compared
%   sorted_cluster_ids    - [1 x nNeurons] vector of cluster IDs corresponding to each column of responses_this_ori (already sorted)
%   ori_letter            - Character representing the stimulus identity (e.g., 'A', 'B', 'C', 'D')
%   znoise_corr_value            - Scalar Pearson correlation coefficient between the two neurons from z-scored data (z-scored by trial-type, used in noise correlation matrices)
%   ori_angle             - Orientation angle in radians (used for display)
%

    x = responses_this_ori(:, neuron_pair(1));
    y = responses_this_ori(:, neuron_pair(2));

    % Add jitter to x and y to see overlapping datapoints
    jitter_amount = 0.1;  % adjust as needed 
    x_jittered = x + randn(size(x)) * jitter_amount;
    y_jittered = y + randn(size(y)) * jitter_amount;

    % Compute correlation and significance
    [R, P] = corr(x, y);
    
    fig=figure;
    scatter(x_jittered, y_jittered, 20, 'filled');
    hold on;
    
    % Fit line
    coeffs = polyfit(x, y, 1);
    x_fit = linspace(min(x), max(x), 100);
    y_fit = polyval(coeffs, x_fit);
    plot(x_fit, y_fit, 'r--', 'LineWidth', 2);
    
    % Annotate
    xlabel(sprintf('Neuron %d (cluster ID %d) Trial Spikecount', neuron_pair(1), sorted_cluster_ids(neuron_pair(1))));
    ylabel(sprintf('Neuron %d (cluster ID %d) Trial Spikecount', neuron_pair(2), sorted_cluster_ids(neuron_pair(2))));
    if contains(Stimulus_type, 'TRAIN')
        title(sprintf(['Orientation %s (%d°): r = %.2f (raw), p = %.3g\n' ...
                   '(r = %.2f znoise-scored)'], ori_letter, round(rad2deg(ori_angle)), R, P, znoise_corr_value));
    elseif contains(Stimulus_type, 'GAVNIK')
        title(sprintf(['Orientation %s (%d°): r = %.2f (raw), p = %.3g\n' ...
                   '(r = %.2f znoise-scored)'], ori_letter, ori_angle, R, P, znoise_corr_value));
    end

    grid on;
    
    % Save figure
    fig_filename = sprintf('%s_%s_%s_%s_NoiseScatter_%dvs%d - %s (%.2f–%.2fs after stim onset), gOSI over %0.2f.fig', ...
        subject_number, Stimulus_type, depth_for_analysis, ori_letter, ...
        sorted_cluster_ids(neuron_pair(1)), sorted_cluster_ids(neuron_pair(2)),temporal_structure, stim_window(1), stim_window(2), gOSI_threshold);
    
    savefig(fig, fullfile(pwd, fig_filename));

end