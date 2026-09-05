function boot_output = plot_UP_survival_counting_process(T, feature_names, varargin)
% PLOT_UP_SURVIVAL_COUNTING_PROCESS Performs multivariable Cox proportional 
% hazards regression on counting-process interval data and plots survival outputs.
%
% Inputs:
%   T             - Table created by build_UP_counting_process_intervals
%   feature_names - Cell array of string column names in T to include as covariates
%                   e.g., {'inRipple', 'ripplePower', 'rippleHPC_MUA_sum', 'cumRippleHPC_MUA'}
%
% Optional Parameters:
%   'title_name'      - Title string for figures (default: 'UP Transition Cox Model')
%   'feature_labels'  - Cell array of display labels for covariates
%   'strata_var'      - Column name for session/subject stratification (default: 'session_id')
%   'nBoot'           - Number of bootstrap iterations (default: 1000)
%   'timebin'         - Bin size for survival curve plotting (default: 0.015)
%   'max_time'        - Upper limit of time axis in seconds (default: 0.3)
%   'stratify_feature'- Column name used to split survival curves into Low/High groups
%   'zscore_option'   - Boolean flag to z-score covariates before Cox fit (default: true)

p = inputParser;
addParameter(p, 'title_name', 'UP Transition Cox Model', @ischar);
addParameter(p, 'feature_labels', {}, @iscell);
addParameter(p, 'strata_var', 'session_id', @ischar);
addParameter(p, 'nBoot', 1000, @isnumeric);
addParameter(p, 'timebin', 0.015, @isnumeric);
addParameter(p, 'max_time', 0.3, @isnumeric);
addParameter(p, 'stratify_feature', '', @ischar);
addParameter(p, 'zscore_option', true, @(x) islogical(x) || isnumeric(x));

parse(p, varargin{:});

title_name       = p.Results.title_name;
feature_labels   = p.Results.feature_labels;
strata_var       = p.Results.strata_var;
nBoot            = p.Results.nBoot;
timebin          = p.Results.timebin;
max_time         = p.Results.max_time;
stratify_feature = p.Results.stratify_feature;
zscore_option    = logical(p.Results.zscore_option);

if isempty(feature_labels)
    feature_labels = feature_names;
end

% Extract matrices from table
X      = T{:, feature_names};
T_time = [T.start, T.stop];
C      = T.censoring; % 1 = censored, 0 = event (MATLAB coxphfit convention)
S      = T{:, strata_var};
upID   = T.upID;

% Filter valid rows without NaNs
valid_idx = find(all(~isnan(X), 2) & all(~isnan(T_time), 2) & ~isnan(C) & ~isnan(S));

X      = X(valid_idx, :);
T_time = T_time(valid_idx, :);
C      = C(valid_idx);
S      = S(valid_idx);
upID   = upID(valid_idx);

nVars = length(feature_names);
uniqueUPs = unique(upID);
nUPs = length(uniqueUPs);

% Z-score continuous/numeric covariates to standardize Cox coefficients (b per 1 SD)
mean_X = mean(X, 1, 'omitnan');
std_X  = std(X, 0, 1, 'omitnan');

if zscore_option
    for v = 1:nVars
        if std_X(v) > 0
            X(:, v) = (X(:, v) - mean_X(v)) / std_X(v);
        end
    end
end

% 1. Full Dataset Cox Model Fit
try
    [b_full, ~, ~, stats_full] = coxphfit(X, T_time, 'Censoring', C, 'Strata', S);
    p_full = stats_full.p;
catch ME
    warning('Full Cox model fit failed: %s', ME.message);
    b_full = nan(nVars, 1);
    p_full = nan(nVars, 1);
end

% 2. Bootstrap Cox Regression (Resampling unique UP events)
boot_b = nan(nBoot, nVars);
boot_p = nan(nBoot, nVars);

parfor iBoot = 1:nBoot
    s = RandStream('mrg32k3a', 'Seed', iBoot);
    
    % Resample UP events with replacement
    sampled_up_indices = datasample(s, 1:nUPs, nUPs, 'Replace', true);
    
    % Map sampled UP IDs back to interval rows
    boot_row_indices = [];
    for u = 1:nUPs
        targetUP = uniqueUPs(sampled_up_indices(u));
        rows = find(upID == targetUP);
        boot_row_indices = [boot_row_indices; rows];
    end
    
    X_boot      = X(boot_row_indices, :);
    T_time_boot = T_time(boot_row_indices, :);
    C_boot      = C(boot_row_indices);
    S_boot      = S(boot_row_indices);
    
    try
        [b_tmp, ~, ~, stats_tmp] = coxphfit(X_boot, T_time_boot, ...
            'Censoring', C_boot, 'Strata', S_boot);
        boot_b(iBoot, :) = b_tmp';
        boot_p(iBoot, :) = stats_tmp.p';
    catch
        boot_b(iBoot, :) = NaN;
        boot_p(iBoot, :) = NaN;
    end
end

% Clean bootstrap outputs
valid_boots = ~any(isnan(boot_b), 2);
boot_b_clean = boot_b(valid_boots, :);
boot_p_clean = boot_p(valid_boots, :);

barData = median(boot_b_clean, 1);
lowerCI = barData - prctile(boot_b_clean, 2.5, 1);
upperCI = prctile(boot_b_clean, 97.5, 1) - barData;
p50     = median(boot_p_clean, 1);

boot_output.b_full = b_full;
boot_output.p_full = p_full;
boot_output.b_boot = boot_b_clean;
boot_output.p_boot = boot_p_clean;
boot_output.mean_X = mean_X;
boot_output.std_X  = std_X;
boot_output.zscored = zscore_option;
boot_output.feature_names = feature_names;
boot_output.feature_labels = feature_labels;

% 3. Plot Figure 1: Cox Coefficients Bar Chart
fig1 = figure('Name', [title_name ' - Standardized Cox Coefficients']);
fig1.Position = [350, 59, min(900, 220 * nVars + 200), 550];

x_pos = 1:nVars;
hold on;
barColors = lines(nVars);

for v = 1:nVars
    bar(x_pos(v), barData(v), 0.5, 'FaceColor', barColors(v, :), ...
        'EdgeColor', 'none', 'FaceAlpha', 0.6);
    errorbar(x_pos(v), barData(v), lowerCI(v), upperCI(v), ...
        'k', 'linestyle', 'none', 'linewidth', 1.5, 'CapSize', 8);
    
    % Annotate p-value
    y_text = barData(v) + sign(barData(v)) * (upperCI(v) + 0.02 * max(abs(barData)));
    if barData(v) < 0
        y_text = barData(v) - lowerCI(v) - 0.05 * max(abs(barData));
    end
    text(x_pos(v), y_text, sprintf('p_{50%%} = %.3e', p50(v)), ...
        'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold');
end

yline(0, '--k', 'LineWidth', 1);
xlim([0.3, nVars + 0.7]);
xticks(x_pos);
xticklabels(feature_labels);
xtickangle(25);
if zscore_option
    ylabel('Standardized Cox Coefficient (b per SD)');
else
    ylabel('Raw Cox Coefficient (b)');
end
title([title_name ' (Multivariable Cox Model)']);
set(gca, 'TickDir', 'out', 'Box', 'off', 'Color', 'none', 'FontSize', 12);

% 4. Plot Figure 2: Stratified Survival Curves (if stratify_feature is specified)
if isempty(stratify_feature)
    stratify_feature = feature_names{1};
end

if ismember(stratify_feature, T.Properties.VariableNames)
    fig2 = figure('Name', [title_name ' - Survival Curves']);
    fig2.Position = [700, 120, 800, 600];
    
    % Extract UP event-level feature value and termination time
    up_feats = nan(nUPs, 1);
    up_times = nan(nUPs, 1);
    
    for u = 1:nUPs
        uID = uniqueUPs(u);
        rows = find(T.upID == uID);
        termRow = rows(end);
        up_times(u) = T.stop(termRow);
        
        vals = T{rows, stratify_feature};
        if any(vals > 0)
            up_feats(u) = max(vals(vals > 0));
        else
            up_feats(u) = max(vals);
        end
    end
    
    % Percentiles calculated among non-zero entries if applicable
    nonzero_feats = up_feats(up_feats > 0);
    if length(nonzero_feats) > 10
        low_thresh  = prctile(nonzero_feats, 25);
        high_thresh = prctile(nonzero_feats, 75);
    else
        low_thresh  = prctile(up_feats, 25);
        high_thresh = prctile(up_feats, 75);
    end
    
    binEdges   = 0:timebin:max_time;
    binCenters = binEdges(1:end-1) + diff(binEdges)/2;
    nGrid      = length(binCenters);
    
    low_mask  = (up_feats <= low_thresh);
    high_mask = (up_feats >= high_thresh);
    
    groups = {low_mask, high_mask};
    group_labels = {'Low (<=25th pct)', 'High (>=75th pct)'};
    group_colors = [0, 90, 50; 74, 20, 134] / 256;
    
    P = zeros(2, 1);
    for g = 1:2
        g_times = up_times(groups{g});
        
        if isempty(g_times)
            continue;
        end
        
        y_boot = zeros(nBoot, nGrid);
        for iBoot = 1:nBoot
            s = RandStream('mrg32k3a', 'Seed', iBoot);
            samp_times = g_times(datasample(s, 1:length(g_times), length(g_times), 'Replace', true));
            
            [temp_y, temp_x] = ecdf(samp_times, 'Function', 'survivor');
            temp_x(isnan(temp_y)) = []; temp_y(isnan(temp_y)) = [];
            [temp_x, uniq_idx] = unique(temp_x);
            temp_y = temp_y(uniq_idx);
            
            y_boot(iBoot, :) = interp1(temp_x, temp_y, binCenters, 'previous', 'extrap');
        end
        
        col = group_colors(g, :);
        x_plot = [0, binCenters];
        y_plot_boot = [ones(nBoot, 1), y_boot];
        
        y_mean = mean(y_plot_boot, 1)';
        LCI    = prctile(y_plot_boot, 2.5, 1)';
        UCI    = prctile(y_plot_boot, 97.5, 1)';
        
        % Subplot 1: Survival probability
        subplot(1, 2, 1);
        hold on;
        plot(x_plot, y_mean, 'Color', col, 'LineWidth', 2);
        P(g) = patch([x_plot, fliplr(x_plot)], [UCI', fliplr(LCI')], col, ...
            'FaceAlpha', 0.3, 'LineStyle', 'none');
        ylabel('Survival Probability of UP');
        xlabel('Time (s)');
        xlim([0, max_time]);
        ylim([0, 1.05]);
        set(gca, 'TickDir', 'out', 'Box', 'off', 'Color', 'none', 'FontSize', 12);
        
        % Subplot 2: Slope / Hazard Derivative
        subplot(1, 2, 2);
        hold on;
        y_diff_boot = diff(y_plot_boot, 1, 2);
        y_diff_mean = [0; mean(y_diff_boot, 1)'];
        LCI_diff    = [0; prctile(y_diff_boot, 2.5, 1)'];
        UCI_diff    = [0; prctile(y_diff_boot, 97.5, 1)'];
        
        plot(x_plot, y_diff_mean, 'Color', col, 'LineWidth', 2);
        patch([x_plot, fliplr(x_plot)], [UCI_diff', fliplr(LCI_diff')], col, ...
            'FaceAlpha', 0.3, 'LineStyle', 'none');
        ylabel('Slope of Survival Probability (\Delta S / \Delta t)');
        xlabel('Time (s)');
        xlim([0, max_time]);
        set(gca, 'TickDir', 'out', 'Box', 'off', 'Color', 'none', 'FontSize', 12);
    end
    
    subplot(1, 2, 1);
    legend(P(1:2), group_labels, 'Box', 'off', 'Location', 'southwest');
    subplot(1, 2, 2);
    legend(P(1:2), group_labels, 'Box', 'off', 'Location', 'northeast');
end

end
