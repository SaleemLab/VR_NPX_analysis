function boot_output = plot_last_ripple_UP_survival(T, feature_names, varargin)
% PLOT_LAST_RIPPLE_UP_SURVIVAL Evaluates and plots Cox proportional hazard
% models and Kaplan-Meier survival curves for UP state termination after the last ripple.
%
% Inputs:
%   T             - Table created by build_last_ripple_to_DOWN_summary_table.m
%   feature_names - Cell array of string feature names in T to evaluate
%                   (e.g., {'last_ripple_HPC_MUA_sum', 'past_ripples_HPC_MUA_sum', 'non_ripple_HPC_MUA_sum'})
%
% Options:
%   'feature_labels'  - Cell array of user-friendly feature display labels
%   'title_name'      - Character array for figure title prefix
%   'timebin'         - Bin size for survival curve interpolation (default 0.015 s)
%   'nBoot'           - Number of bootstrap iterations (default 1000)
%   'strata_var'      - Column name in T for Cox stratification (default 'subject_id')
%   'max_time'        - Maximum time on survival x-axis (default 0.3 s)
%   'zscore_features' - Logical, whether to z-score features for comparable beta coefficients (default true)
%   'model_type'      - 'multivariable' (default) or 'univariable' / 'univariate'
%   'plot_ecdf'       - Logical, whether to plot Kaplan-Meier / ECDF curves (default true)
%   'plot_bar'        - Logical, whether to plot Cox coefficient bar chart (default true)

p = inputParser;
addParameter(p, 'feature_labels', {}, @iscell);
addParameter(p, 'title_name', '', @ischar);
addParameter(p, 'timebin', 0.015, @isnumeric);
addParameter(p, 'nBoot', 1000, @isnumeric);
addParameter(p, 'strata_var', 'subject_id', @ischar);
addParameter(p, 'max_time', 0.3, @isnumeric);
addParameter(p, 'zscore_features', true, @islogical);
addParameter(p, 'model_type', 'multivariable', @(x) ischar(x) || isstring(x));
addParameter(p, 'plot_ecdf', false, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'plot_bar', true, @(x) islogical(x) || isnumeric(x));

parse(p, varargin{:});

feature_labels = p.Results.feature_labels;
title_name     = p.Results.title_name;
timebin        = p.Results.timebin;
nBoot          = p.Results.nBoot;
strata_var     = p.Results.strata_var;
max_time       = p.Results.max_time;
zscore_feat    = p.Results.zscore_features;
model_type_str = lower(char(p.Results.model_type));
plot_ecdf      = logical(p.Results.plot_ecdf);
plot_bar       = logical(p.Results.plot_bar);

is_multivariable = ~ismember(model_type_str, {'univariable', 'univariate'});

if isempty(feature_labels)
    feature_labels = feature_names;
end

% 1. Filter valid rows where outcome and features are non-NaN
if plot_bar
    valid_idx = ~isnan(T.last_ripple_to_UP_term) & T.last_ripple_to_UP_term > 0;
    for f = 1:length(feature_names)
        valid_idx = valid_idx & ~isnan(T.(feature_names{f}));
    end

    T_sub = T(valid_idx, :);
else
    valid_idx = ~isnan(T.last_ripple_to_UP_term) & T.last_ripple_to_UP_term > 0;
    T_sub = T(valid_idx, :);
end

cox_time = T_sub.last_ripple_to_UP_term;
X_feat   = zeros(height(T_sub), length(feature_names));
for f = 1:length(feature_names)
    vals = T_sub.(feature_names{f});
    if zscore_feat
        std_val = std(vals, 'omitnan');
        if std_val > 0
            vals = (vals - mean(vals, 'omitnan')) / std_val;
        end
    end
    X_feat(:, f) = vals;
end

if ismember(strata_var, T_sub.Properties.VariableNames)
    strata_used = T_sub.(strata_var);
else
    strata_used = ones(height(T_sub), 1);
end

nFeat = length(feature_names);
nObs  = length(cox_time);

%% 2. Bootstrapped Cox Regression (Multivariable or Univariable)
if is_multivariable
    fprintf('Running %d bootstrap iterations for multivariable Cox regression...\n', nBoot);
else
    fprintf('Running %d bootstrap iterations for univariable Cox regression (per feature)...\n', nBoot);
end

if plot_bar
boot_b = nan(nBoot, nFeat);
boot_p = nan(nBoot, nFeat);

parfor iBoot = 1:nBoot
    s = RandStream('mrg32k3a', 'Seed', iBoot);
    boot_idx = datasample(s, 1:nObs, nObs);
    
    time_samp   = cox_time(boot_idx);
    X_samp      = X_feat(boot_idx, :);
    strata_samp = strata_used(boot_idx);
    
    if is_multivariable
        try
            [b_tmp, ~, ~, stats_tmp] = coxphfit(X_samp, time_samp, ...
                'Censoring', zeros(size(time_samp)), ...
                'Strata', strata_samp);
            boot_b(iBoot, :) = b_tmp';
            boot_p(iBoot, :) = stats_tmp.p';
        catch
            % Fallback without stratification if singular matrix occurs
            try
                [b_tmp, ~, ~, stats_tmp] = coxphfit(X_samp, time_samp, ...
                    'Censoring', zeros(size(time_samp)));
                boot_b(iBoot, :) = b_tmp';
                boot_p(iBoot, :) = stats_tmp.p';
            catch
                boot_b(iBoot, :) = NaN;
                boot_p(iBoot, :) = NaN;
            end
        end
    else
        % Univariable fit per feature individually
        b_vec = nan(1, nFeat);
        p_vec = nan(1, nFeat);
        for f = 1:nFeat
            try
                [b_tmp, ~, ~, stats_tmp] = coxphfit(X_samp(:, f), time_samp, ...
                    'Censoring', zeros(size(time_samp)), ...
                    'Strata', strata_samp);
                b_vec(f) = b_tmp;
                p_vec(f) = stats_tmp.p;
            catch
                try
                    [b_tmp, ~, ~, stats_tmp] = coxphfit(X_samp(:, f), time_samp, ...
                        'Censoring', zeros(size(time_samp)));
                    b_vec(f) = b_tmp;
                    p_vec(f) = stats_tmp.p;
                catch
                    b_vec(f) = NaN;
                    p_vec(f) = NaN;
                end
            end
        end
        boot_b(iBoot, :) = b_vec;
        boot_p(iBoot, :) = p_vec;
    end
end

% Compute bootstrap summary statistics
barData = nan(1, nFeat);
lowerCI = nan(1, nFeat);
upperCI = nan(1, nFeat);
p50     = nan(1, nFeat);

for f = 1:nFeat
    vals_b = boot_b(~isnan(boot_b(:, f)), f);
    vals_p = boot_p(~isnan(boot_p(:, f)), f);
    if ~isempty(vals_b)
        barData(f) = median(vals_b);
        lowerCI(f) = median(vals_b) - prctile(vals_b, 2.5);
        upperCI(f) = prctile(vals_b, 97.5) - median(vals_b);
        p50(f)     = prctile(vals_p, 50);
    end
end

boot_output.b = boot_b;
boot_output.p = boot_p;
boot_output.feature_names  = feature_names;
boot_output.feature_labels = feature_labels;
boot_output.barData    = barData;
boot_output.lowerCI    = lowerCI;
boot_output.upperCI    = upperCI;
boot_output.p50        = p50;
boot_output.model_type = model_type_str;
end

%% 3. Plot Cox Hazard Coefficients Bar Chart (Narrower & Tighter)
if plot_bar
    px_per_bar     = 30;   
    left_margin    = 85;    
    right_margin   = 45;
    bottom_margin  = 90;  
    top_margin     = 55;    
    fig_height     = 400;   

    total_width = left_margin + right_margin + (nFeat * px_per_bar);

    % Create Figure
    fig_bar = figure('Units', 'pixels');
    fig_bar.Position = [300, 300, total_width, fig_height];
    if ~isempty(title_name)
        fig_bar.Name = [title_name, ' - Cox Coefficients'];
    end

    % Set Axes with Fixed Pixel Margins
    ax = axes('Parent', fig_bar, 'Units', 'pixels', ...
        'Position', [left_margin, bottom_margin, nFeat * px_per_bar, fig_height - bottom_margin - top_margin]);

    hold(ax, 'on');
    colors = lines(nFeat);
    x_pos  = 1:nFeat;

    % Zero reference line
    yline(ax, 0, 'Color', [0.4 0.4 0.4], 'LineStyle', ':', 'LineWidth', 1);

    % Global offset for p-value labels based on overall scale
    y_scale_ref = max(abs([barData(:) + upperCI(:); barData(:) - lowerCI(:)]));
    if isempty(y_scale_ref) || isnan(y_scale_ref) || y_scale_ref == 0
        y_scale_ref = 1;
    end
    y_offset = 0.04 * y_scale_ref;

    for f = 1:nFeat
        % Bar: width set to 0.4 (narrower profile)
        bar(ax, x_pos(f), barData(f), 0.4, 'FaceColor', colors(f, :), ...
            'EdgeColor', 'none', 'FaceAlpha', 0.7);
        
        % Errorbar
        errorbar(ax, x_pos(f), barData(f), lowerCI(f), upperCI(f), ...
            'Color', 'k', 'LineStyle', 'none', 'LineWidth', 1.3, 'CapSize', 4);
        
        % Adaptive label positioning above/below CI whiskers
        if barData(f) >= 0
            y_text = barData(f) + upperCI(f) + y_offset;
            va = 'bottom';
        else
            y_text = barData(f) - lowerCI(f) - y_offset;
            va = 'top';
        end
        
        text(ax, x_pos(f), y_text, sprintf('p = %.2e', p50(f)), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', va, ...
            'FontSize', 8.5);
    end

    % Formatting
    xlim(ax, [0.45, nFeat + 0.55]);
    xticks(ax, x_pos);
    xticklabels(ax, feature_labels);
    xtickangle(ax, 30);

    if zscore_feat
        ylabel(ax, 'Standardized Cox Hazard Coefficient (\beta per SD)');
    else
        ylabel(ax, 'Cox Hazard Coefficient (\beta)');
    end

    if is_multivariable
        title(ax, ['Multivariable Cox Regression: ', title_name], 'Interpreter', 'none');
    else
        title(ax, ['Univariable Cox Regression: ', title_name], 'Interpreter', 'none');
    end
    set(ax, 'TickDir', 'out', 'Box', 'off', 'Color', 'none', 'FontSize', 10);
end

%% 4. Plot Kaplan-Meier Survival Curves per Feature (Optional)
if plot_ecdf
    binEdges   = 0:timebin:max_time;
    binCenters = binEdges(1:end-1) + diff(binEdges)/2;

    for f = 1:nFeat
        feat_vals = X_feat(:, f);
        low_thresh  = prctile(feat_vals, 20);
        high_thresh = prctile(feat_vals, 80);
        % low_thresh  = prctile(feat_vals, 50);
        % high_thresh = prctile(feat_vals, 50);
        
        fig_km = figure;
        fig_km.Position = [400 200 700 320];
        if ~isempty(title_name)
            fig_km.Name = sprintf('%s - %s KM Survival', title_name, feature_labels{f});
        end
        
        colour_lines = [44, 162, 95; 215, 48, 39] / 255; % Low = green, High = red
        
        P = [];
        for n = 1:2 % 1 = Low (<=20th pct), 2 = High (>=80th pct)
            if n == 1
                grp_idx = find(feat_vals <= low_thresh);
            else
                grp_idx = find(feat_vals >= high_thresh);
            end
            
            y_boot = zeros(1000, length(binCenters));
            for iBoot = 1:1000
                s = RandStream('mrg32k3a', 'Seed', iBoot);
                sampled_times = cox_time(datasample(s, grp_idx, length(grp_idx)));
                
                [temp_y, temp_x] = ecdf(sampled_times, 'Function', 'survivor');
                temp_x(isnan(temp_y)) = []; temp_y(isnan(temp_y)) = [];
                [temp_x, uIdx] = unique(temp_x);
                temp_y = temp_y(uIdx);
                y_boot(iBoot, :) = interp1(temp_x, temp_y, binCenters, 'previous', 'extrap');
            end
            
            x_plot = [0, binCenters];
            y_boot_plot = [ones(1000, 1), y_boot];
            
            y_mean = mean(y_boot_plot, 1)';
            LCI    = prctile(y_boot_plot, 2.5, 1)';
            UCI    = prctile(y_boot_plot, 97.5, 1)';
            
            % Plot Survival Curve S(t)
            subplot(1, 2, 1);
            hold on;
            plot(x_plot, y_mean, 'Color', colour_lines(n, :), 'LineWidth', 1.8);
            patch([x_plot, fliplr(x_plot)], [UCI', fliplr(LCI')], colour_lines(n, :), ...
                'FaceAlpha', 0.25, 'LineStyle', 'none');
            ylabel('Survival Probability of UP');
            xlabel('Time from Last Ripple (s)');
            xlim([0 max_time]);
            ylim([0 1.05]);
            title(feature_labels{f});
            set(gca, 'TickDir', 'out', 'Box', 'off', 'Color', 'none', 'FontSize', 11);
            
            % Plot Derivative dS(t)/dt
            subplot(1, 2, 2);
            hold on;
            y_diff_boot = diff(y_boot_plot, 1, 2);
            y_diff_mean = [0; mean(y_diff_boot, 1)'];
            LCI_diff    = [0; prctile(y_diff_boot, 2.5, 1)'];
            UCI_diff    = [0; prctile(y_diff_boot, 97.5, 1)'];
            
            plot(x_plot, y_diff_mean, 'Color', colour_lines(n, :), 'LineWidth', 1.8);
            P(n) = patch([x_plot, fliplr(x_plot)], [UCI_diff', fliplr(LCI_diff')], colour_lines(n, :), ...
                'FaceAlpha', 0.25, 'LineStyle', 'none');
            ylabel('Slope of Survival Probability');
            xlabel('Time from Last Ripple (s)');
            xlim([0 max_time]);
            set(gca, 'TickDir', 'out', 'Box', 'off', 'Color', 'none', 'FontSize', 11);
        end
        
        legend(P(1:2), {'Low', 'High'}, 'Box', 'off', 'Location', 'best');
    end
end

end
