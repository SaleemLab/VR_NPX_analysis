%% plot_V1_ripple_regression_highlight.m
% Re-plots the track difference vs ripple difference regression for ALL,
% LOW, and HIGH ripple powers across 3 time windows.
% Allows highlighting specific cells and plotting their PSTHs.

clear all;
analysis_folder = 'C:\Users\masah\OneDrive\Documents\corticohippocampal_replay';
if ~exist(analysis_folder, 'dir')
    analysis_folder = 'D:\corticohippocampal_replay';
end
if ~exist(analysis_folder, 'dir')
    analysis_folder = 'P:\corticohippocampal_replay';
end

%% CONFIGURATION
% Specify cells to highlight as [session_id, cell_id]
cells_to_highlight = []; % e.g., [16, 54; 20, 102]

%% Load Data
% load_path = fullfile(analysis_folder, 'V1-HPC sleep reactivation', 'reactivation coherence SUA', 'V1_ripple_spatial_PSTH.mat');
load_path = fullfile(analysis_folder, 'V1-HPC sleep reactivation', 'reactivation coherence SUA', 'V1_ripple_spatial_PSTH_z.mat');
fprintf('Loading data from %s...\n', load_path);
if exist(load_path, 'file')
    load(load_path, 'V1_cells_data', 'x_bins');
else
    error('Extracted PSTH data not found. Please run extract_V1_ripple_PSTH.m first.');
end

n_cells = length(V1_cells_data);
windows = {'ALL', 'PRE', 'POST'};
win_ranges = {[-0.2, 0.2], [-0.2, 0], [0, 0.2]};

% Arrays to hold regression data
x_track = zeros(n_cells, 1);
session_ids = zeros(n_cells, 1);
cell_ids = zeros(n_cells, 1);
animal_ids = zeros(n_cells, 1);

y_all = zeros(n_cells, 3);
y_high = zeros(n_cells, 3);
y_low = zeros(n_cells, 3);

for i = 1:n_cells
    data = V1_cells_data(i);
    x_track(i) = data.z_track_diff;
    session_ids(i) = data.session_id;
    cell_ids(i) = data.cell_id;
    animal_ids(i) = data.animal_id;
    
    for w = 1:3
        idx = x_bins >= win_ranges{w}(1) & x_bins <= win_ranges{w}(2);
        
        y_all(i, w) = mean(data.z_all_T1_m(idx)) - mean(data.z_all_T2_m(idx));
        y_high(i, w) = mean(data.z_high_T1_m(idx)) - mean(data.z_high_T2_m(idx));
        y_low(i, w) = mean(data.z_low_T1_m(idx)) - mean(data.z_low_T2_m(idx));
    end
end

%% Plotting Regressions
conditions = {'All ripples', 'Low ripples', 'High ripples'};
y_data_all = {y_all, y_low, y_high};
colors = [175, 141, 195; 127, 191, 123] / 256; % Purple and Green

for cond = 1:3
    Fig = figure('Name', sprintf('V1 context ripple modulation (%s)', conditions{cond}), 'Position', [60 60 1270 400]);
    y_cond = y_data_all{cond};
    
    for w = 1:3
        subplot(1, 3, w);
        x = x_track;
        y = y_cond(:, w);
        
        % Filter NaNs
        valid = ~isnan(x) & ~isnan(y);
        x_val = x(valid); y_val = y(valid);
        s_ids = session_ids(valid); c_ids = cell_ids(valid);
        
        % Fit line
        if length(x_val) > 1
            fitted = polyfit(x_val, y_val, 1);
            y_fit = polyval(fitted, [min(x_val) max(x_val)]);
            
            scatter(x_val, y_val, 20, 'k', 'filled', 'MarkerFaceAlpha', 0.1, 'MarkerEdgeAlpha', 0.1);
            hold on;
            plot([min(x_val) max(x_val)], y_fit, 'r-', 'LineWidth', 2);
        end
        xline(0, '--'); yline(0, '--');
        
        % Highlight specific cells
        for k = 1:size(cells_to_highlight, 1)
            hi_session = cells_to_highlight(k, 1);
            hi_cell = cells_to_highlight(k, 2);
            idx_hi = find(s_ids == hi_session & c_ids == hi_cell);
            if ~isempty(idx_hi)
                col = colors(mod(k-1, size(colors, 1)) + 1, :);
                scatter(x_val(idx_hi), y_val(idx_hi), 100, col, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
            end
        end
        
        ylim([-0.2 0.2]); xlim([-2.5 2.5]);
        xlabel('Track difference (z)'); ylabel('Ripple difference (z)');
        title(sprintf('%s - %s', conditions{cond}, windows{w}));
        set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
    end
end

%% Plot Highlighted Cells PSTH
for k = 1:size(cells_to_highlight, 1)
    hi_session = cells_to_highlight(k, 1);
    hi_cell = cells_to_highlight(k, 2);
    
    idx = find(session_ids == hi_session & cell_ids == hi_cell, 1);
    if isempty(idx)
        warning('Cell %d in session %d not found.', hi_cell, hi_session);
        continue;
    end
    
    data = V1_cells_data(idx);
    fig = figure('Name', sprintf('Highlighted Cell (Session %d, Cell %d)', hi_session, hi_cell), 'Position', [100 100 1500 400]);
    
    % Spatial PSTH
    subplot(1, 4, 1);
    pos_bins = 1:1:140;
    plot_psth(pos_bins, data.spatial_1_m, data.spatial_1_se, [0 0 1]); hold on;
    plot_psth(pos_bins, data.spatial_2_m, data.spatial_2_se, [1 0 0]);
    xlabel('Position (cm)'); ylabel('Firing Rate (Hz)');
    title(sprintf('Spatial PSTH\nZ-Diff: %.3f', data.z_track_diff));
    set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
    
    % Ripple PSTH ALL
    subplot(1, 4, 2);
    plot_psth(x_bins, data.raw_all_T1_m, data.raw_all_T1_se, [0 0 1]); hold on;
    plot_psth(x_bins, data.raw_all_T2_m, data.raw_all_T2_se, [1 0 0]);
    xlim([-0.5 0.5]); xline(0,'--');
    xlabel('Time from ripple (s)'); ylabel('Firing Rate (Hz)');
    title('Ripple PSTH (All)');
    set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
    
    % Ripple PSTH HIGH
    subplot(1, 4, 3);
    plot_psth(x_bins, data.raw_high_T1_m, data.raw_high_T1_se, [0 0 1]); hold on;
    plot_psth(x_bins, data.raw_high_T2_m, data.raw_high_T2_se, [1 0 0]);
    xlim([-0.5 0.5]); xline(0,'--');
    xlabel('Time from ripple (s)'); ylabel('Firing Rate (Hz)');
    title('Ripple PSTH (High)');
    set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
    
    % Ripple PSTH LOW
    subplot(1, 4, 4);
    plot_psth(x_bins, data.raw_low_T1_m, data.raw_low_T1_se, [0 0 1]); hold on;
    plot_psth(x_bins, data.raw_low_T2_m, data.raw_low_T2_se, [1 0 0]);
    xlim([-0.5 0.5]); xline(0,'--');
    xlabel('Time from ripple (s)'); ylabel('Firing Rate (Hz)');
    title('Ripple PSTH (Low)');
    set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
end

% Helper function for plotting PSTH with shaded error manually
function plot_psth(x_vals, m, se, color)
    if isempty(m) || all(isnan(m(:))), return; end
    m = m(:)'; se = se(:)'; x_vals = x_vals(:)';
    fill([x_vals fliplr(x_vals)], [m-se fliplr(m+se)], color, 'EdgeColor', 'none', 'FaceAlpha', 0.15);
    hold on;
    plot(x_vals, m, 'Color', color, 'LineWidth', 1.2);
end
