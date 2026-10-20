%% Code to 
analysis_folder = pwd;
load(fullfile(analysis_folder,'processed_data','spatial_map_all.mat'),'spatial_map_all');
% Initialize containers to gather data across all sessions
num_sessions =22;

%% Spatial tuning across all sessions
num_sessions = length(spatial_map_all.spatial_response);

% Initialize containers to gather data across all sessions
all_maps_trackL_odd  = [];
all_maps_trackL_even = [];
all_maps_trackR_odd  = [];
all_maps_trackR_even = [];
all_cell_labels      = {};

for nsession = 1:num_sessions
    % Fetch data structures for this session
    spatial_data = spatial_map_all.spatial_response{nsession};
    labels_data  = spatial_map_all.region{nsession};
    
    num_cells = size(spatial_data, 1);
    if num_cells == 0, continue; end
    
    % Store brain region labels for this session's cells
    all_cell_labels = [all_cell_labels; labels_data(:, 1)];
    
    % Temporary holders for this session's calculated averages
    sess_L_odd  = [];
    sess_L_even = [];
    sess_R_odd  = [];
    sess_R_even = [];
    
    for ncell = 1:num_cells
        % Track L (Column 1) & Track R (Column 2)
        map_L = spatial_data{ncell, 1}; 
        map_R = spatial_data{ncell, 2}; 
        
        % Track L Mean Maps (Odd vs. Even Laps)
        if ~isempty(map_L)
            sess_L_odd(ncell, :)  = mean(map_L(1:2:end, :), 1, 'omitnan');
            sess_L_even(ncell, :) = mean(map_L(2:2:end, :), 1, 'omitnan');
        else
            sess_L_odd(ncell, :)  = zeros(1, 140); 
            sess_L_even(ncell, :) = zeros(1, 140);
        end
        
        % Track R Mean Maps (Odd vs. Even Laps)
        if ~isempty(map_R)
            sess_R_odd(ncell, :)  = mean(map_R(1:2:end, :), 1, 'omitnan');
            sess_R_even(ncell, :) = mean(map_R(2:2:end, :), 1, 'omitnan');
        else
            sess_R_odd(ncell, :)  = zeros(1, 140);
            sess_R_even(ncell, :) = zeros(1, 140);
        end
    end
    
    % Append current session data to global matrix
    all_maps_trackL_odd  = [all_maps_trackL_odd;  sess_L_odd];
    all_maps_trackL_even = [all_maps_trackL_even; sess_L_even];
    all_maps_trackR_odd  = [all_maps_trackR_odd;  sess_R_odd];
    all_maps_trackR_even = [all_maps_trackR_even; sess_R_even];
end

% Replace remaining NaNs with 0
all_maps_trackL_odd(isnan(all_maps_trackL_odd))   = 0;
all_maps_trackL_even(isnan(all_maps_trackL_even)) = 0;
all_maps_trackR_odd(isnan(all_maps_trackR_odd))   = 0;
all_maps_trackR_even(isnan(all_maps_trackR_even)) = 0;

num_bins = size(all_maps_trackL_odd, 2); % Usually 140


% Loop Through the Target Regions and Generate Plots
% Define the list of regions you want to cycle through
target_regions = {'HPC', 'V1_L', 'V1_R'};

for r = 1:length(target_regions)
    target_region = target_regions{r};
    
    % Find cells matching the current region keyword
    target_indices = find(contains(all_cell_labels, target_region));
    
    if isempty(target_indices)
        warning('No cells found matching the region filter: %s. Skipping...', target_region);
        continue;
    end
    
    % Extract target sub-matrices for this region
    region_L_odd  = all_maps_trackL_odd(target_indices, :);
    region_L_even = all_maps_trackL_even(target_indices, :);
    region_R_odd  = all_maps_trackR_odd(target_indices, :);
    region_R_even = all_maps_trackR_even(target_indices, :);
    
    % --- Cross-Track Normalization & Sorting (ODD Laps) ---
    combined_odd = [region_L_odd, region_R_odd];
    norm_combined_odd = normalize(combined_odd, 2, 'range');
    
    norm_L_odd = norm_combined_odd(:, 1:num_bins);
    norm_R_odd = norm_combined_odd(:, num_bins+1:end);
    
    % Find independent peak locations for sorting reference
    [~, peak_idx_L] = max(norm_L_odd, [], 2);
    [~, peak_idx_R] = max(norm_R_odd, [], 2);
    
    [~, sort_order_L] = sort(peak_idx_L);
    [~, sort_order_R] = sort(peak_idx_R);
    
    % --- Cross-Track Normalization (EVEN Laps for Plotting) ---
    combined_even = [region_L_even, region_R_even];
    norm_combined_even = normalize(combined_even, 2, 'range');
    
    norm_L_even = norm_combined_even(:, 1:num_bins);
    norm_R_even = norm_combined_even(:, num_bins+1:end);
    
    % --- Generate Figures (2x2 Layout) ---
    fig = figure('Position', [50 + (r*40), 50 + (r*40), 700, 600]);
    sgtitle(sprintf('Population Maps for %s Cells (Cross-Track Normalized)', target_region), ...
        'FontSize', 15, 'FontWeight', 'bold');
    
    % 1. Track L (Even) sorted by Track L (Odd)
    subplot(2, 2, 1);
    imagesc(norm_L_even(sort_order_L, :));
    axis xy; colormap(flip(gray)); colorbar;
    title('Track L (Even Laps) | Sorted by Track L (Odd)');
    xlabel('Position Bins'); ylabel('Cells'); clim([0.1 1]);
    
    % 2. Track L (Even) sorted by Track R (Odd)
    subplot(2, 2, 2);
    imagesc(norm_L_even(sort_order_R, :));
    axis xy; colormap(flip(gray)); colorbar;
    title('Track L (Even Laps) | Sorted by Track R (Odd)');
    xlabel('Position Bins'); ylabel('Cells'); clim([0.1 1]);
    
    % 3. Track R (Even) sorted by Track L (Odd)
    subplot(2, 2, 3);
    imagesc(norm_R_even(sort_order_L, :));
    axis xy; colormap(flip(gray)); colorbar;
    title('Track R (Even Laps) | Sorted by Track L (Odd)');
    xlabel('Position Bins'); ylabel('Cells'); clim([0.1 1]);
    
    % 4. Track R (Even) sorted by Track R (Odd)
    subplot(2, 2, 4);
    imagesc(norm_R_even(sort_order_R, :));
    axis xy; colormap(flip(gray)); colorbar;
    title('Track R (Even Laps) | Sorted by Track R (Odd)');
    xlabel('Position Bins'); ylabel('Cells'); clim([0.1 1]);
end


%%%%
%%%%
%%%%
%% visualise Mean FR diff between Track L and Track R for V1 and HC
figure
subplot(2,2,1)
hold on
V1_L = [];V1_L_z=[];
V1_R = [];V1_R_z=[];
for nsession = 1:22
    V1_R_this = spatial_map_all.mean_FR{nsession}(contains(spatial_map_all.region{nsession},'V1_R'),:);
    V1_L_this = spatial_map_all.mean_FR{nsession}(contains(spatial_map_all.region{nsession},'V1_L'),:);
    scatter(V1_R_this(:,1),V1_R_this(:,2),'b','filled', ...
        'MarkerFaceAlpha', 0.1, ...  % Set face transparency (0 = clear, 1 = opaque)
        'MarkerEdgeAlpha', 0.1);     % Set edge transparency
    scatter(V1_L_this(:,1),V1_L_this(:,2), 'r','filled', ...
        'MarkerFaceAlpha', 0.1, ...  % Set face transparency (0 = clear, 1 = opaque)
        'MarkerEdgeAlpha', 0.1);     % Set edge transparency

    V1_L = [V1_L; V1_L_this];
    V1_R = [V1_R; V1_R_this];

    V1_R_this = spatial_map_all.mean_FR_z{nsession}(contains(spatial_map_all.region{nsession},'V1_R'),:);
    V1_L_this = spatial_map_all.mean_FR_z{nsession}(contains(spatial_map_all.region{nsession},'V1_L'),:);
    V1_L_z = [V1_L_z; V1_L_this];
    V1_R_z = [V1_R_z; V1_R_this];
end
xlim([0 50]);ylim([0 50])
plot([0 50],[0 50],'k--')
xlabel('Track L mean firing rate (Hz)')
ylabel('Track R mean firing rate (Hz)')
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);

subplot(2,2,2)
histogram(V1_L(:,1)-V1_L(:,2),-20:1:20,'Normalization','probability')
hold on;
histogram(V1_R(:,1)-V1_R(:,2),-20:1:20,'Normalization','probability')
xlabel('Mean Track L - R firing rate difference (Hz)')
ylabel('Proportion of cells')
legend('V1 L','V1 R','box','off')
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);

subplot(2,2,3)
histogram(V1_L_z(:,1)-V1_L_z(:,2),-2:0.1:2,'Normalization','probability')
hold on;
histogram(V1_R_z(:,1)-V1_R_z(:,2),-2:0.1:2,'Normalization','probability')
xlabel('Mean Track L - R firing rate difference (z)')
ylabel('Proportion of cells')
sgtitle('V1')
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);

figure
subplot(2,2,1)
hold on
HC_L = [];HC_L_z=[];
HC_R = [];HC_R_z=[];
for nsession = 1:22
    HC_R_this = spatial_map_all.mean_FR{nsession}(contains(spatial_map_all.region{nsession},'HPC_R'),:);
    HC_L_this = spatial_map_all.mean_FR{nsession}(contains(spatial_map_all.region{nsession},'HPC_L'),:);
    scatter(HC_R_this(:,1),HC_R_this(:,2),'b','filled', ...
        'MarkerFaceAlpha', 0.1, ...  % Set face transparency (0 = clear, 1 = opaque)
        'MarkerEdgeAlpha', 0.1);     % Set edge transparency
    scatter(HC_L_this(:,1),HC_L_this(:,2), 'r','filled', ...
        'MarkerFaceAlpha', 0.1, ...  % Set face transparency (0 = clear, 1 = opaque)
        'MarkerEdgeAlpha', 0.1);     % Set edge transparency

    HC_L = [HC_L; HC_L_this];
    HC_R = [HC_R; HC_R_this];

    HC_R_this = spatial_map_all.mean_FR_z{nsession}(contains(spatial_map_all.region{nsession},'HPC_R'),:);
    HC_L_this = spatial_map_all.mean_FR_z{nsession}(contains(spatial_map_all.region{nsession},'HPC_L'),:);
    HC_L_z = [HC_L_z; HC_L_this];
    HC_R_z = [HC_R_z; HC_R_this];
end
xlim([0 50]);ylim([0 50])
plot([0 50],[0 50],'k--')
xlabel('Track L mean firing rate (Hz)')
ylabel('Track R mean firing rate (Hz)')
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);

subplot(2,2,2)
histogram(HC_L(:,1)-HC_L(:,2),-20:1:20,'Normalization','probability')
hold on;
histogram(HC_R(:,1)-HC_R(:,2),-20:1:20,'Normalization','probability')
xlabel('Mean Track L - R firing rate difference (Hz)')
ylabel('Proportion of cells')
legend('HC L','HC R','box','off')

subplot(2,2,3)
histogram(HC_L_z(:,1)-HC_L_z(:,2),-2:0.1:2,'Normalization','probability')
hold on;
histogram(HC_R_z(:,1)-HC_R_z(:,2),-2:0.1:2,'Normalization','probability')
xlabel('Mean Track L - R firing rate difference (z)')
ylabel('Proportion of cells')
sgtitle('HC')
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);


%%%%
%%%%
%%%%
%% Track preference of one-track vs both-track tuned cells, split by hemisphere
% A cell is "significant" on a track if odd/even stability > 0.95 (and, if
% use_peak_percentile, peak percentile > 0.95). Cells significant on exactly
% one track are "one-track tuned". Expectation: V1_L prefers Track R and
% V1_R prefers Track L (contralateral), while HPC_L/HPC_R show no hemispheric
% bias for either one-track or both-track tuned cells.
sig_threshold       = 0.95;
use_peak_percentile = true;
min_cells_session   = 3;   % min one-track cells per hemisphere for session-level fraction
out_dir = fullfile(analysis_folder, 'track_tuning_hemisphere');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

% Pool cells across sessions
cell_FR = []; cell_FR_z = []; cell_sig = []; cell_region = strings(0,1); cell_session = [];
for nsession = 1:num_sessions
    region_this = string(spatial_map_all.region{nsession});
    if isempty(region_this), continue; end
    sig_this = spatial_map_all.odd_even_stability_RUN1{nsession} > sig_threshold;
    if use_peak_percentile
        sig_this = sig_this & spatial_map_all.peak_percentile_RUN1{nsession} > sig_threshold;
    end
    cell_FR      = [cell_FR;      spatial_map_all.mean_FR{nsession}];
    cell_FR_z    = [cell_FR_z;    spatial_map_all.mean_FR_z{nsession}];
    cell_sig     = [cell_sig;     sig_this];
    cell_region  = [cell_region;  region_this(:,1)];
    cell_session = [cell_session; nsession*ones(numel(region_this),1)];
end
is_L_only = cell_sig(:,1) & ~cell_sig(:,2);
is_R_only = ~cell_sig(:,1) & cell_sig(:,2);
is_both   = cell_sig(:,1) & cell_sig(:,2);
is_one    = is_L_only | is_R_only;
% Track selectivity index: +1 = fires only on Track L, -1 = only on Track R
sel_index = (cell_FR(:,1) - cell_FR(:,2)) ./ (cell_FR(:,1) + cell_FR(:,2));

two_sided_binom = @(k, n) min(1, 2*min(binocdf(k, n, 0.5), 1 - binocdf(k-1, n, 0.5)));

area_names = {'V1', 'HPC'};
hemi_names = {'L', 'R'};
hemi_cols  = {'r', 'b'};   % same colours as the scatter plots above (left = red, right = blue)
fprintf('\n=== One-track vs both-track tuning by hemisphere (threshold %.2f, peak percentile used: %d) ===\n', ...
    sig_threshold, use_peak_percentile);

for a = 1:numel(area_names)
    area = area_names{a};
    in_hemi = cell(1,2);
    for h = 1:2
        in_hemi{h} = cell_region == sprintf('%s_%s', area, hemi_names{h});
    end

    fig = figure('Position', [100 + a*40, 100 + a*40, 1300, 750], 'Color', 'w');

    % --- (1) Proportion of tuned cells: L-only / R-only / both ---
    subplot(2,3,1)
    frac = zeros(2,3); n_tuned = zeros(1,2);
    for h = 1:2
        n_tuned(h) = sum(in_hemi{h} & (is_one | is_both));
        frac(h,:)  = [sum(in_hemi{h} & is_L_only), sum(in_hemi{h} & is_R_only), ...
            sum(in_hemi{h} & is_both)] / max(n_tuned(h), 1);
    end
    b = bar(frac, 'stacked');
    b(1).FaceColor = [0.85 0.4 0.1]; b(2).FaceColor = [0.2 0.5 0.8]; b(3).FaceColor = [0.6 0.6 0.6];
    set(gca, 'XTickLabel', {sprintf('%s_L', area), sprintf('%s_R', area)}, 'TickLabelInterpreter', 'none');
    for h = 1:2, text(h, 1.03, sprintf('n=%d', n_tuned(h)), 'HorizontalAlignment', 'center'); end
    ylabel('Proportion of tuned cells'); ylim([0 1.1]); xlim([0.4 3.6]);
    legend({'Track L only', 'Track R only', 'Both'}, 'Box', 'off', 'Location', 'east');
    set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);

    % --- (2) Fraction preferring Track L among one-track cells (pooled + per session) ---
    subplot(2,3,2)
    hold on
    k_L = zeros(1,2); n_one = zeros(1,2); p_binom = zeros(1,2);
    for h = 1:2
        k_L(h)     = sum(in_hemi{h} & is_L_only);
        n_one(h)   = sum(in_hemi{h} & is_one);
        p_binom(h) = two_sided_binom(k_L(h), n_one(h));
        bar(h, k_L(h)/max(n_one(h),1), 0.6, 'FaceColor', hemi_cols{h}, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    end
    sessions = unique(cell_session);
    sess_frac = nan(numel(sessions), 2);
    for s = 1:numel(sessions)
        for h = 1:2
            idx = cell_session == sessions(s) & in_hemi{h} & is_one;
            if sum(idx) >= min_cells_session
                sess_frac(s,h) = sum(idx & is_L_only) / sum(idx);
            end
        end
    end
    paired = all(~isnan(sess_frac), 2);
    plot([1 2], sess_frac(paired,:)', '-', 'Color', [0.7 0.7 0.7]);
    for h = 1:2
        scatter(h*ones(sum(~isnan(sess_frac(:,h))),1), sess_frac(~isnan(sess_frac(:,h)),h), 25, hemi_cols{h}, 'filled', ...
            'MarkerFaceAlpha', 0.6);
    end
    plot([0.5 2.5], [0.5 0.5], 'k--');
    % Pooled 2x2 test (hemisphere x preferred track) and session-level paired test
    [~, p_fisher] = fishertest([k_L(1), n_one(1)-k_L(1); k_L(2), n_one(2)-k_L(2)]);
    if sum(paired) >= 3
        p_sess = signrank(sess_frac(paired,1), sess_frac(paired,2));
    else
        p_sess = NaN;
    end
    title({sprintf('Fisher p=%.2g', p_fisher), sprintf('session signrank p=%.2g (n=%d)', p_sess, sum(paired))}, 'FontWeight', 'normal');
    set(gca, 'XTick', [1 2], 'XTickLabel', {sprintf('%s_L', area), sprintf('%s_R', area)}, 'TickLabelInterpreter', 'none');
    for h = 1:2, text(h, 1.03, sprintf('n=%d', n_one(h)), 'HorizontalAlignment', 'center'); end
    xlim([0.5 2.5]); ylim([0 1.1]);
    ylabel('Frac. one-track cells on Track L');
    set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);

    % --- (3) Mean FR on Track L vs Track R for one-track cells ---
    subplot(2,3,3)
    hold on
    for h = 1:2
        idx = in_hemi{h} & is_one;
        scatter(cell_FR(idx,1), cell_FR(idx,2), 15, hemi_cols{h}, 'filled', 'MarkerFaceAlpha', 0.3, 'MarkerEdgeAlpha', 0.3);
    end
    max_FR = max([prctile(cell_FR(is_one & (in_hemi{1} | in_hemi{2}), :), 99, 'all'), 1]);
    plot([0 max_FR], [0 max_FR], 'k--');
    xlim([0 max_FR]); ylim([0 max_FR]); axis square
    xlabel('Track L mean firing rate (Hz)'); ylabel('Track R mean firing rate (Hz)');
    legend({sprintf('%s_L', area), sprintf('%s_R', area)}, 'Box', 'off', 'Interpreter', 'none');
    title('One-track tuned cells', 'FontWeight', 'normal');
    set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);

    % --- (4-6) Selectivity index / z-scored FR difference distributions ---
    dist_panels = {
        sel_index,                    is_one,  -1:0.1:1, 'Selectivity (L-R)/(L+R)',             'One-track tuned';
        sel_index,                    is_both, -1:0.1:1, 'Selectivity (L-R)/(L+R)',             'Both-track tuned';
        cell_FR_z(:,1)-cell_FR_z(:,2), is_one,  -2:0.1:2, 'Mean Track L - R firing rate diff (z)', 'One-track tuned'};
    for d = 1:size(dist_panels,1)
        subplot(2,3,3+d)
        hold on
        vals = cell(1,2);
        for h = 1:2
            vals{h} = dist_panels{d,1}(in_hemi{h} & dist_panels{d,2});
            vals{h} = vals{h}(~isnan(vals{h}));
            histogram(vals{h}, dist_panels{d,3}, 'Normalization', 'probability', ...
                'FaceColor', hemi_cols{h}, 'EdgeColor', 'none', 'FaceAlpha', 0.4);
        end
        yl = ylim;
        for h = 1:2
            plot(median(vals{h})*[1 1], yl, '-', 'Color', hemi_cols{h}, 'LineWidth', 1.5);
        end
        plot([0 0], yl, 'k--');
        if all(cellfun(@numel, vals) > 0)
            p_rs = ranksum(vals{1}, vals{2});
        else
            p_rs = NaN;
        end
        title(sprintf('%s | ranksum p=%.2g', dist_panels{d,5}, p_rs), 'FontWeight', 'normal');
        xlabel(dist_panels{d,4}); ylabel('Proportion of cells');
        legend({sprintf('%s_L (n=%d)', area, numel(vals{1})), sprintf('%s_R (n=%d)', area, numel(vals{2}))}, ...
            'Box', 'off', 'Interpreter', 'none', 'Location', 'northwest');
        set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
    end
    sgtitle(sprintf('%s: track preference by hemisphere', area), 'FontWeight', 'bold');

    % Summary to command window
    for h = 1:2
        fprintf('%s_%s: tuned=%d | L-only=%d, R-only=%d, both=%d | frac L among one-track=%.2f (binomial vs 0.5 p=%.2g)\n', ...
            area, hemi_names{h}, n_tuned(h), k_L(h), n_one(h)-k_L(h), sum(in_hemi{h} & is_both), ...
            k_L(h)/max(n_one(h),1), p_binom(h));
    end
    fprintf('%s hemisphere x preferred track: Fisher p=%.2g; session-level signrank p=%.2g (n=%d sessions)\n', ...
        area, p_fisher, p_sess, sum(paired));

    savefig(fig, fullfile(out_dir, sprintf('%s_track_preference_by_hemisphere.fig', area)));
    exportgraphics(fig, fullfile(out_dir, sprintf('%s_track_preference_by_hemisphere.pdf', area)), 'ContentType', 'vector');
end