% For making CDF plots and performing the Kolmogorov Smirnov test on signal
% and noise correlation matrix data across days. Run this code after running
% main_PCA_of_PSTH_Ellie.m for all the days required (that code saves the
% necessary data in the struct used by this code).
% ERB 2025
clear all
% === User parameters ===
base_path = 'V:\Ellie\DATA\SUBJECTS\M00013\analysis';
cd(base_path);
day_folders = {'20250203', '20250204', '20250205', '20250206', '20250207'};  % your session dates
%'20250203', '20250204', '20250205', '20250206', '20250207', '20250210', '20250211'
%'20250217', '20250218', '20250219', '20250220', '20250221'
Stimulus_type = 'TRAIN';  % your stimulus type folder name here 'TRAIN'  'GAVNIK_ABCD'

% === Load correlation data across days ===
signal_data = cell(1, length(day_folders));
noise_data = cell(length(day_folders), 1);

for d = 1:length(day_folders)
    day_path = fullfile(base_path, day_folders{d}, Stimulus_type);
    corr_file = fullfile(day_path, 'correlation_summary.mat');
    
    if exist(corr_file, 'file')
        loaded = load(corr_file);
        if isfield(loaded, 'correlation_summary')
            signal_data{d} = loaded.correlation_summary.unique_signal_values;
            % On first day, load metadata like ordered_oris, labels, etc.
            if d == 1
                ordered_oris = loaded.correlation_summary.ordered_oris;
                letters = loaded.correlation_summary.letters;
                colors = loaded.correlation_summary.colors;
                ori_labels = loaded.correlation_summary.ori_labels;
                subject = loaded.correlation_summary.subject;
                stimulus_type = loaded.correlation_summary.stimulus_type;
                depth_for_analysis = loaded.correlation_summary.depth;
                stim_window = loaded.correlation_summary.stim_window;
                gOSI_threshold = loaded.correlation_summary.gOSI_threshold;
                n_noise_types = length(ordered_oris);  % Now we can define this
                noise_data = cell(length(day_folders), n_noise_types);  % Resize with known value
            end    

            % Load noise data for this session
            for o = 1:length(loaded.correlation_summary.unique_noise_values)
                noise_data{d, o} = loaded.correlation_summary.unique_noise_values{o};
            end
        else
            error('correlation_summary struct not found in %s', corr_file);
        end
    else
        error('File not found: %s', corr_file);
    end
end


% === Plot Signal correlations figure ===
figure; hold on;
signal_color = [0.2 0.2 0.2]; % dark grey base
for d = 1:length(day_folders)
    shade = 1 - (d-1)*0.15;  % lighter for later days
    c = signal_color * shade + (1 - shade); % blend towards white
    [f, x] = ecdf(signal_data{d});
    
    % Use dashed lines for first and last day, solid otherwise
    if d == 1 || d == length(day_folders)
        ls = '--';
        lw = 3;  % thicker
    else
        ls = '-';
        lw = 2;
    end
    
    % Plot the CDF
    plot(x, f, 'LineWidth', lw, 'Color', c, 'LineStyle', ls);
    
    % Plot the median
    med_val = median(signal_data{d});
    plot([med_val med_val], [0 0.5], 'LineStyle', ls, 'Color', c, 'LineWidth', lw - 0.5, 'HandleVisibility', 'off');
end

xlabel('Pairwise Correlation (r)');
ylabel('CDF');
xlim([-1 1]);
legend(day_folders, 'Location', 'Best');
box on;

% KS test for signal correlations
day1_idx = 1;
dayN_idx = length(day_folders);
[h_sig, p_sig, D_sig] = kstest2(signal_data{day1_idx}, signal_data{dayN_idx});
title(sprintf('%s %s %s Signal Correlations (%.2f–%.2fs after stim onset, gOSI > %0.2f), Across Days (%s vs %s KS D=%.3f, p=%.3g)', subject, stimulus_type, depth_for_analysis,...
    stim_window(1), stim_window(2), gOSI_threshold, day_folders{day1_idx}, day_folders{dayN_idx}, D_sig, p_sig), 'Interpreter', 'none');
fig_filename_signal = sprintf('%s %s %s Signal Correlations Across Days (%s to %s).fig', ...
                        subject, stimulus_type, depth_for_analysis, day_folders{day1_idx}, day_folders{dayN_idx});
%savefig(gcf, fig_filename_signal); 


% === Plot Noise correlations, one figure per orientation ===
fig_handles_noise = cell(n_noise_types, 1);         % To store figure handles
fig_filenames_noise = cell(n_noise_types, 1);       % To store corresponding filenames

for ori_idx = 1:n_noise_types
    fig_handles_noise{ori_idx} = figure; hold on;   % Store figure handle
    base_color = colors(ori_idx, :);
    for d = 1:length(day_folders)
        shade = 1 - (d-1)*0.15;  % lighter for later days
        c = base_color * shade + (1 - shade); % blend towards white
        [f, x] = ecdf(noise_data{d, ori_idx});
        % Use dashed lines for first and last day, solid otherwise
        if d == 1 || d == length(day_folders)
            ls = '--';
            lw = 3;  % thicker
        else
            ls = '-';
            lw = 2;  
        end
        
        % Plot the CDF
        plot(x, f, 'LineWidth', lw, 'Color', c, 'LineStyle', ls);

        % Plot matching median line
        med_val = median(noise_data{d, ori_idx});
        plot([med_val med_val], [0 0.5], 'LineStyle', ls, 'Color', c, 'LineWidth', lw - 0.5, 'HandleVisibility', 'off');
    end
    xlabel('Pairwise Correlation (r)');
    ylabel('CDF');
    xlim([-0.5 0.7]);
    legend(day_folders, 'Location', 'Best');
    box on;

    % KS test
    [h_noise, p_noise, D_noise] = kstest2(noise_data{day1_idx, ori_idx}, noise_data{dayN_idx, ori_idx});
    title(sprintf('%s %s %s Noise Correlations %s (%.2f–%.2fs after stim onset, gOSI > %0.2f), Across Days (%s vs %s KS D=%.3f, p=%.3g)', subject, stimulus_type, depth_for_analysis,...
        ori_labels{ori_idx}, stim_window(1), stim_window(2), gOSI_threshold, day_folders{day1_idx}, day_folders{dayN_idx}, D_noise, p_noise), 'Interpreter', 'none');
    fig_filenames_noise{ori_idx} = sprintf('%s %s %s Noise Correlations %s Across Days (%s to %s).fig', ...
                        subject, stimulus_type, depth_for_analysis, ori_labels{ori_idx}, day_folders{day1_idx}, day_folders{dayN_idx});

end
%for ori_idx = 1:n_noise_types
%    savefig(fig_handles_noise{ori_idx}, fig_filenames_noise{ori_idx});
%end



% === Plot Unit Firing Rates (trial averages) across Days for Each Orientation ===
fig_handles_fr = cell(n_noise_types, 1);         % To store figure handles
fig_filenames_fr = cell(n_noise_types, 1);       % To store filenames
stim_duration = diff(stim_window);               % Duration in seconds

for ori_idx = 1:n_noise_types
    fig_handles_fr{ori_idx} = figure; hold on;
    base_color = colors(ori_idx, :);
    
    all_rates = [];     % Collect all mean firing rates
    all_groups = [];    % Corresponding group indices
    all_colors = [];    % Colors for each day’s violins

    for d = 1:length(day_folders)
        % Reload responses
        day_path = fullfile(base_path, day_folders{d}, Stimulus_type);
        corr_file = fullfile(day_path, 'correlation_summary.mat');
        loaded = load(corr_file);

        responses = loaded.correlation_summary.responses_this_ori{ori_idx};  % trials × neurons
        mean_spike_counts = mean(responses, 1);  % average across trials
        mean_firing_rates = mean_spike_counts / stim_duration;

        % Accumulate data
        all_rates = [all_rates; mean_firing_rates(:)];
        all_groups = [all_groups; repmat(d, size(mean_firing_rates(:)))];

        % Color styling
        shade = 1 - (d-1)*0.15;  % lighter for later days
        this_color = base_color * shade + (1 - shade);
        all_colors = [all_colors; repmat(this_color, length(mean_firing_rates), 1)];
    end
    
    % Use violinplot.m (File Exchange version)
    h = violinplot(all_rates, all_groups, 'ShowData', false);  % thick central vertical range is the interquartile range (25th-75th percentile). Thinner vertical line is the data range

    % Apply color to each violin
    for d = 1:length(h)
        shade = 1 - (d-1)*0.15;
        c = base_color * shade + (1 - shade);
        h(d).ViolinColor = {c};
        h(d).EdgeColor = c * 0.6;
        h(d).BoxColor = c * 0.5;
        h(d).MedianColor = [0 0 0];  % central black circle is the median firing rate
    end

    % Add legend for median marker
    hold on;
    h_legend = scatter(NaN, NaN, 50, 'ko', 'filled'); % black filled circle marker for median
    legend(h_legend, 'Median unit firing rate', 'Location', 'best');
    hold off;


    % Label and style
    xticks(1:length(day_folders));
    xticklabels(day_folders);
    xlabel('Day');
    ylabel('V1 unit firing rates, averaged across trials (Hz)');
    title(sprintf('%s %s %s unit firing rates – %s (%.2f–%.2fs after stim onset, gOSI > %0.2f)', subject, stimulus_type, depth_for_analysis, ori_labels{ori_idx}, stim_window(1), stim_window(2), gOSI_threshold), 'Interpreter', 'none');
    box on;

    % Save filename
    fig_filenames_fr{ori_idx} = sprintf('%s %s %s unit FR violin plots %s (%s to %s) (%.2f to %.2fs after stim onset, gOSI over %0.2f).fig', ...
        subject, stimulus_type, depth_for_analysis, letters(ori_idx), day_folders{1}, day_folders{end}, stim_window(1), stim_window(2), gOSI_threshold);
end

% === (Optional) Save figures ===
for ori_idx = 1:n_noise_types
    savefig(fig_handles_fr{ori_idx}, fig_filenames_fr{ori_idx});
end






% === Scatter plots of signal vs noise correlation, one per orientation ===
figure_handles = cell(n_noise_types, 1);

for ori_idx = 1:n_noise_types
    figure_handles{ori_idx} = figure; hold on;
    yline(0, 'k-', 'LineWidth', 0.5, 'HandleVisibility', 'off'); % line at zero noise correlation
    base_color = colors(ori_idx, :);
    legend_entries = cell(length(day_folders), 1);

    for d = 1:length(day_folders)
        % Load signal and noise correlations
        signal_vals = signal_data{d};
        noise_vals = noise_data{d, ori_idx};

        if length(signal_vals) ~= length(noise_vals)
            warning('Mismatch in signal/noise data length for day %s, orientation %s', ...
                    day_folders{d}, ori_labels{ori_idx});
            continue;
        end

        % Shade based on day
        shade = 1 - (d-1)*0.15;
        c = base_color * shade + (1 - shade); % blend toward white

        % Scatter plot of signal vs noise correlations
        scatter(signal_vals, noise_vals, 10, 'MarkerEdgeColor', c, ...
            'MarkerFaceColor', c, 'MarkerFaceAlpha', 0.4, 'MarkerEdgeAlpha', 0.4, 'HandleVisibility', 'off');
        
        % Fit line of best fit (least squares linear regression)
        coeffs = polyfit(signal_vals, noise_vals, 1);
        x_fit = linspace(-1, 1, 100);
        y_fit = polyval(coeffs, x_fit);
        % Set dashed for first and last day
        if d == 1 || d == length(day_folders)
            line_style = ':';
            line_width = 3;
        else
            line_style = '-';
            line_width = 2;
        end


        plot(x_fit, y_fit, 'LineStyle', line_style, 'Color', c, 'LineWidth', line_width);

        % Compute Pearson correlation coefficient
        [r_val, p_val] = corr(signal_vals, noise_vals, 'Type', 'Pearson');

        % Store legend entry
        legend_entries{d} = sprintf('Day %s (r=%.2f, p=%.3f)', day_folders{d}, r_val, p_val);
    end

    % Labeling and formatting
    xlabel('Signal Correlation (r)');
    ylabel('Noise Correlation (r)');
    title(sprintf('%s %s - %s - %s Unit Pair Signal vs Noise Correlations (%.2f–%.2fs after stim onset, gOSI > %0.2f) by Day ', ...
        subject, stimulus_type, ori_labels{ori_idx}, depth_for_analysis, stim_window(1), stim_window(2), gOSI_threshold), ...
        'Interpreter', 'none');
    
    xlim([-1 1]); ylim([-0.5 0.6]);
    legend(legend_entries, 'Location', 'Best');
    box on;
    axis square;

    % Optionally save the figure
    fig_filename_scatter = sprintf('%s %s %s Unit Pair Signal_vs_Noise_Corrs Orientation_%s (%s to %s).fig', ...
        subject, stimulus_type, depth_for_analysis, ori_labels{ori_idx}, ...
        day_folders{1}, day_folders{end});
    savefig(figure_handles{ori_idx}, fig_filename_scatter);
end





% Store proportion of opposite-sign correlations across days and orientations
opposite_sign_props = zeros(length(day_folders), n_noise_types);  % days × orientations

% Store the cluster ID pairs for each day and orientation where signs are opposite
opposite_sign_pairs = cell(length(day_folders), n_noise_types);  % each cell: n_pairs × 2
for d = 1:length(day_folders)
    day_path = fullfile(base_path, day_folders{d}, Stimulus_type);
    corr_file = fullfile(day_path, 'correlation_summary.mat');
    
    if exist(corr_file, 'file')
        loaded = load(corr_file);
        if isfield(loaded, 'correlation_summary')
            sig_vals = loaded.correlation_summary.unique_signal_values;
            cluster_pairs = loaded.correlation_summary.cluster_id_pairs;

            for ori_idx = 1:n_noise_types
                noise_vals = loaded.correlation_summary.unique_noise_values{ori_idx};

                if length(sig_vals) ~= length(noise_vals)
                    error('Mismatch in pair counts for day %s, orientation %s', ...
                        day_folders{d}, ori_labels{ori_idx});
                end

                % Logical index of opposite-sign cases
                opposite_sign = (sig_vals > 0.01 & noise_vals < -0.01) | (sig_vals < -0.01 & noise_vals > 0.01);
                opposite_sign_props(d, ori_idx) = mean(opposite_sign);

                % Save the corresponding cluster ID pairs for inspection
                opposite_sign_pairs{d, ori_idx} = cluster_pairs(opposite_sign, :);
            end
        else
            error('correlation_summary struct not found in %s', corr_file);
        end
    else
        error('File not found: %s', corr_file);
    end
end

% === Plot bar chart of opposite-sign proportions ===
figure;
bar(opposite_sign_props, 'grouped');
ylim([0 1]);
ylabel('Proportion of Cluster Pairs with Opposite Sign');
xlabel('Day');
set(gca, 'XTickLabel', day_folders);
legend(ori_labels, 'Location', 'Best');  % A, B, C, D
title(sprintf('%s %s %s Signal vs Noise Correlation: Opposite Sign Proportions (%.2f–%.2fs after stim onset, gOSI > %0.2f), (%s to %s)', subject, stimulus_type, depth_for_analysis,...
        stim_window(1), stim_window(2), gOSI_threshold, day_folders{day1_idx}, day_folders{dayN_idx}), 'Interpreter', 'none');
box on;
fig_filename_bars = sprintf('%s %s %s Signal vs Noise Correlation: Opposite Sign Proportions (%s to %s).fig', ...
                        subject, stimulus_type, depth_for_analysis, day_folders{day1_idx}, day_folders{dayN_idx});
%savefig(gcf, fig_filename_bars); 

%%% For auditing/sense-checking. Example: show opposite-sign pairs for day 2 and orientation B
% day_idx = 1;
% ori_idx = 2;  % B
% disp(['Opposite-sign cluster pairs on ', day_folders{day_idx}, ...
%       ' for orientation ', letters(ori_idx), ':']);
% disp(opposite_sign_pairs{day_idx, ori_idx});


