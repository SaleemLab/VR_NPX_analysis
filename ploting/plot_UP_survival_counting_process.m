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
addParameter(p, 'timebin', 0.02, @isnumeric);
addParameter(p, 'max_time', 1, @isnumeric);
addParameter(p, 'min_time', 0.1, @isnumeric);
addParameter(p, 'stratify_feature', '', @ischar);
addParameter(p, 'zscore_option', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'bootstrap', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'plot_survival', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'plot_bar', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'unadjusted_by_session', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'is_multivariate', true, @(x) islogical(x) || isnumeric(x));

parse(p, varargin{:});

title_name            = p.Results.title_name;
feature_labels        = p.Results.feature_labels;
strata_var            = p.Results.strata_var;
nBoot                 = p.Results.nBoot;
timebin               = p.Results.timebin;
max_time              = p.Results.max_time;
min_time              = p.Results.min_time;

stratify_feature      = p.Results.stratify_feature;
zscore_option         = logical(p.Results.zscore_option);
bootstrap             = logical(p.Results.bootstrap);
plot_survival         = logical(p.Results.plot_survival);
plot_bar              = logical(p.Results.plot_bar);
unadjusted_by_session = logical(p.Results.unadjusted_by_session);
is_multivariate       = logical(p.Results.is_multivariate);


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
b_full  = nan(nVars, 1);
p_full  = nan(nVars, 1);
se_full = nan(nVars, 1);
ci_full = nan(nVars, 2);

if is_multivariate
    try
        [b_full, ~, ~, stats_full] = coxphfit(X, T_time, 'Censoring', C, 'Strata', S);
        p_full  = stats_full.p;
        se_full = stats_full.se;
        if isfield(stats_full, 'ci') && ~isempty(stats_full.ci)
            ci_full = stats_full.ci;
        else
            ci_full = [b_full - 1.96 * se_full, b_full + 1.96 * se_full];
        end
    catch ME
        warning('Full multivariable Cox model fit failed: %s', ME.message);
    end
else
    for v = 1:nVars
        try
            [b_v, ~, ~, stats_v] = coxphfit(X(:, v), T_time, 'Censoring', C, 'Strata', S);
            b_full(v)  = b_v;
            p_full(v)  = stats_v.p;
            se_full(v) = stats_v.se;
            if isfield(stats_v, 'ci') && ~isempty(stats_v.ci)
                ci_full(v, :) = stats_v.ci;
            else
                ci_full(v, :) = [b_v - 1.96 * stats_v.se, b_v + 1.96 * stats_v.se];
            end
        catch ME
            warning('Univariate Cox model fit failed for feature %s: %s', feature_names{v}, ME.message);
        end
    end
end

% 2. Bootstrap Cox Regression (Resampling unique UP events)
if bootstrap
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

        if is_multivariate
            try
                [b_tmp, ~, ~, stats_tmp] = coxphfit(X_boot, T_time_boot, ...
                    'Censoring', C_boot, 'Strata', S_boot);
                boot_b(iBoot, :) = b_tmp';
                boot_p(iBoot, :) = stats_tmp.p';
            catch
                boot_b(iBoot, :) = NaN;
                boot_p(iBoot, :) = NaN;
            end
        else
            b_row = nan(1, nVars);
            p_row = nan(1, nVars);
            for v = 1:nVars
                try
                    [b_tmp, ~, ~, stats_tmp] = coxphfit(X_boot(:, v), T_time_boot, ...
                        'Censoring', C_boot, 'Strata', S_boot);
                    b_row(v) = b_tmp;
                    p_row(v) = stats_tmp.p;
                catch
                end
            end
            boot_b(iBoot, :) = b_row;
            boot_p(iBoot, :) = p_row;
        end
    end

    % Clean bootstrap outputs
    valid_boots = ~any(isnan(boot_b), 2);
    boot_b_clean = boot_b(valid_boots, :);
    boot_p_clean = boot_p(valid_boots, :);

end

if bootstrap
    barData = median(boot_b_clean, 1);
    lowerCI = barData - prctile(boot_b_clean, 2.5, 1);
    upperCI = prctile(boot_b_clean, 97.5, 1) - barData;
    p_plot  = median(boot_p_clean, 1);
    p_label_prefix = 'p_{50%}';
else
    barData = b_full(:)';
    lowerCI = (b_full - ci_full(:, 1))';
    upperCI = (ci_full(:, 2) - b_full)';
    p_plot  = p_full(:)';
    p_label_prefix = 'p';
end

boot_output.b_full = b_full;
boot_output.se_full = se_full;
boot_output.p_full = p_full;
boot_output.ci_full = ci_full;
if bootstrap
boot_output.b_boot = boot_b_clean;
boot_output.p_boot = boot_p_clean;
% boot_output.bootstrap = bootstrap;
end
boot_output.mean_X = mean_X;
boot_output.std_X  = std_X;
boot_output.zscored = zscore_option;
boot_output.feature_names = feature_names;
boot_output.feature_labels = feature_labels;

% 3. Plot Figure 1: Cox Coefficients Bar Chart (Narrower & Fixed Size)
if plot_bar
    px_per_bar     = 50;   
    left_margin    = 85;    
    right_margin   = 45;
    bottom_margin  = 90;  
    top_margin     = 55;    
    fig_height     = 400;   

    total_width = left_margin + right_margin + (nVars * px_per_bar);

    fig1 = figure('Units', 'pixels');
    fig1.Position = [300, 300, total_width, fig_height];
    if is_multivariate
        fig1.Name = [title_name ' (Multivariable Cox Model) - Cox Coefficients'];
    else
        fig1.Name = [title_name ' (Univariate Cox Models) - Cox Coefficients'];
    end

    ax1 = axes('Parent', fig1, 'Units', 'pixels', ...
        'Position', [left_margin, bottom_margin, nVars * px_per_bar, fig_height - bottom_margin - top_margin]);

    hold(ax1, 'on');
    barColors = lines(nVars);
    x_pos  = 1:nVars;

    % Zero reference line
    yline(ax1, 0, 'Color', [0.4 0.4 0.4], 'LineStyle', ':', 'LineWidth', 1);

    % Global offset for p-value labels based on overall scale
    y_scale_ref = max(abs([barData(:) + upperCI(:); barData(:) - lowerCI(:)]));
    if isempty(y_scale_ref) || isnan(y_scale_ref) || y_scale_ref == 0
        y_scale_ref = 1;
    end
    y_offset = 0.04 * y_scale_ref;

    for v = 1:nVars
        % Bar: width set to 0.4 (fixed narrow bar profile matching plot_last_ripple_UP_survival)
        bar(ax1, x_pos(v), barData(v), 0.4, 'FaceColor', barColors(v, :), ...
            'EdgeColor', 'none', 'FaceAlpha', 0.7);
        
        % Errorbar
        errorbar(ax1, x_pos(v), barData(v), lowerCI(v), upperCI(v), ...
            'Color', 'k', 'LineStyle', 'none', 'LineWidth', 1.3, 'CapSize', 4);
        
        % Adaptive label positioning above/below CI whiskers
        if barData(v) >= 0
            y_text = barData(v) + upperCI(v) + y_offset;
            va = 'bottom';
        else
            y_text = barData(v) - lowerCI(v) - y_offset;
            va = 'top';
        end
        
        text(ax1, x_pos(v), y_text, sprintf('%s = %.2e', p_label_prefix, p_plot(v)), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', va, ...
            'FontSize', 8.5);
    end

    xlim(ax1, [0.45, nVars + 0.55]);
    xticks(ax1, x_pos);
    xticklabels(ax1, feature_labels);
    xtickangle(ax1, 30);

    if zscore_option
        ylabel(ax1, 'Standardized Cox Hazard Coefficient (\beta per SD)');
    else
        ylabel(ax1, 'Cox Hazard Coefficient (\beta)');
    end

    if is_multivariate
        title(ax1, ['Multivariable Cox Model: ' title_name], 'Interpreter', 'none');
    else
        title(ax1, ['Univariate Cox Models: ' title_name], 'Interpreter', 'none');
    end
    set(ax1, 'TickDir', 'out', 'Box', 'off', 'Color', 'none', 'FontSize', 10);
end

% 4. Plot Figure 2: Stratified Survival Curves for tested variables
if isempty(stratify_feature) || strcmpi(stratify_feature, 'all')
    features_to_plot = feature_names;
else
    features_to_plot = {stratify_feature};
end

if plot_survival

for fIter = 1:length(features_to_plot)
    currFeat = features_to_plot{fIter};
    currFeatLabel = feature_labels{strcmp(feature_names, currFeat)};
    if isempty(currFeatLabel)
        currFeatLabel = currFeat;
    end
    
    if ismember(currFeat, T.Properties.VariableNames)
        featRaw = T{valid_idx, currFeat};

        % Representative Low/High values (percentiles among non-zero entries if
        % there are enough of them, else among all entries)
        nonzero_feat = featRaw(featRaw > 0);
        if sum(~isnan(nonzero_feat)) > 10
            low_val  = prctile(nonzero_feat, 25);
            high_val = prctile(nonzero_feat, 75);
        else
            low_val  = prctile(featRaw, 25);
            high_val = prctile(featRaw, 75);
        end

        if high_val > low_val
            rowGroup = nan(size(featRaw));
            rowGroup(featRaw <= low_val)  = 0;
            rowGroup(featRaw >= high_val) = 1;

            binEdges = min_time:timebin:max_time;
            if binEdges(end) < max_time
                binEdges = [binEdges, max_time];
            end
            nBins  = numel(binEdges) - 1;
            x_plot = [0, binEdges(2:end)];

            group_colors = [0, 90, 50; 74, 20, 134] / 256;
            legend_labels = {sprintf('Low (\\leq%.3g)', low_val), sprintf('High (\\geq%.3g)', high_val)};

            fig2 = figure('Name', [title_name ' - Survival Curves (' currFeatLabel ')']);
            fig2.Position = [700, 120, 1100, 500];

            % --- Left panel: model-based (adjusted) predicted survival
            subplot(1, 2, 1); hold on;
            featIdx = find(strcmp(feature_names, currFeat), 1);
            if isempty(featIdx)
                text(0.5, 0.5, sprintf('"%s" is not a fitted covariate\n(no Cox coefficient available)', currFeat), ...
                    'Units', 'normalized', 'HorizontalAlignment', 'center');
            elseif any(isnan(b_full))
                text(0.5, 0.5, 'Cox model fit failed; no baseline hazard available', ...
                    'Units', 'normalized', 'HorizontalAlignment', 'center');
            else
                if is_multivariate
                    lp = X * b_full(:);
                    model_subtitle = 'Multivariable Breslow baseline hazard \times Cox coefficient';
                else
                    lp = X(:, featIdx) * b_full(featIdx);
                    model_subtitle = 'Univariate Breslow baseline hazard \times Cox coefficient';
                end
                [Lambda0, stratWeight] = local_breslow_grouped(T_time, C, S, lp, binEdges);
                w = stratWeight / sum(stratWeight);

                if zscore_option && std_X(featIdx) > 0
                    low_z  = (low_val  - mean_X(featIdx)) / std_X(featIdx);
                    high_z = (high_val - mean_X(featIdx)) / std_X(featIdx);
                else
                    low_z  = low_val;
                    high_z = high_val;
                end

                hr_low  = exp(b_full(featIdx) * low_z);
                hr_high = exp(b_full(featIdx) * high_z);
                S_low_model  = [1, w' * exp(-Lambda0 * hr_low)];
                S_high_model = [1, w' * exp(-Lambda0 * hr_high)];

                if bootstrap
                    coefSamples = boot_b_clean(:, featIdx);
                else
                    coefSamples = [b_full(featIdx); ci_full(featIdx, 1); ci_full(featIdx, 2)];
                end
                nSamp = numel(coefSamples);
                S_low_boot  = zeros(nSamp, nBins);
                S_high_boot = zeros(nSamp, nBins);
                for ib = 1:nSamp
                    S_low_boot(ib, :)  = w' * exp(-Lambda0 * exp(coefSamples(ib) * low_z));
                    S_high_boot(ib, :) = w' * exp(-Lambda0 * exp(coefSamples(ib) * high_z));
                end
                if bootstrap
                    LCI_low  = [1, prctile(S_low_boot,  2.5, 1)];
                    UCI_low  = [1, prctile(S_low_boot,  97.5, 1)];
                    LCI_high = [1, prctile(S_high_boot, 2.5, 1)];
                    UCI_high = [1, prctile(S_high_boot, 97.5, 1)];
                else
                    LCI_low  = [1, min(S_low_boot,  [], 1)];
                    UCI_low  = [1, max(S_low_boot,  [], 1)];
                    LCI_high = [1, min(S_high_boot, [], 1)];
                    UCI_high = [1, max(S_high_boot, [], 1)];
                end

                patch([x_plot, fliplr(x_plot)], [UCI_low, fliplr(LCI_low)], group_colors(1, :), 'FaceAlpha', 0.25, 'LineStyle', 'none');
                patch([x_plot, fliplr(x_plot)], [UCI_high, fliplr(LCI_high)], group_colors(2, :), 'FaceAlpha', 0.25, 'LineStyle', 'none');
                p1 = plot(x_plot, S_low_model,  'Color', group_colors(1, :), 'LineWidth', 2);
                p2 = plot(x_plot, S_high_model, 'Color', group_colors(2, :), 'LineWidth', 2);
                legend([p1, p2], legend_labels, 'Box', 'off', 'Location', 'southwest');

                % field_name = matlab.lang.makeValidName(currFeat);
                % boot_output.survival.(field_name).model.time    = x_plot;
                % boot_output.survival.(field_name).model.S_low   = S_low_model;
                % boot_output.survival.(field_name).model.S_high  = S_high_model;
                % boot_output.survival.(field_name).model.CI_low  = [LCI_low; UCI_low];
                % boot_output.survival.(field_name).model.CI_high = [LCI_high; UCI_high];
            end
            title({['Model-based: ' currFeatLabel], model_subtitle});
            ylabel('Predicted survival probability of UP'); xlabel('Time (s)');
            xlim([0, max_time]); ylim([0, 1.05]);
            set(gca, 'TickDir', 'out', 'Box', 'off', 'Color', 'none', 'FontSize', 12);

            % --- Right panel: Simon-Makuch time-dependent-covariate Kaplan-Meier
            subplot(1, 2, 2); hold on;

            if unadjusted_by_session
                unique_sessions = unique(S);
                nSess = length(unique_sessions);
                sess_S_low  = nan(nSess, numel(x_plot));
                sess_S_high = nan(nSess, numel(x_plot));
                
                for iSess = 1:nSess
                    s_mask = (S == unique_sessions(iSess));
                    T_time_s = T_time(s_mask, :);
                    C_s      = C(s_mask);
                    feat_s   = featRaw(s_mask);
                    
                    % Compute within-session 25th / 75th percentiles
                    nonzero_s = feat_s(feat_s > 0);
                    if sum(~isnan(nonzero_s)) > 10
                        low_s  = prctile(nonzero_s, 20);
                        high_s = prctile(nonzero_s, 80);
                    else
                        low_s  = prctile(feat_s, 20);
                        high_s = prctile(feat_s, 80);
                    end
                    
                    if high_s > low_s
                        rowGroup_s = nan(size(feat_s));
                        rowGroup_s(feat_s <= low_s)  = 0;
                        rowGroup_s(feat_s >= high_s) = 1;
                        
                        [sl, sh] = local_simon_makuch_grouped(T_time_s, C_s, rowGroup_s, binEdges);
                        sess_S_low(iSess, :)  = [1, sl];
                        sess_S_high(iSess, :) = [1, sh];
                    % else
                    %     1
                    end


                end
                
                mean_low_sm  = mean(sess_S_low, 1, 'omitnan');
                std_low_sm   = std(sess_S_low, 0, 1, 'omitnan');
                n_low_sm     = sum(~isnan(sess_S_low), 1);
                se_low_sm    = std_low_sm ./ sqrt(max(1, n_low_sm));
                
                mean_high_sm = mean(sess_S_high, 1, 'omitnan');
                std_high_sm  = std(sess_S_high, 0, 1, 'omitnan');
                n_high_sm    = sum(~isnan(sess_S_high), 1);
                se_high_sm   = std_high_sm ./ sqrt(max(1, n_high_sm));
                
                se_low_sm(1)  = 0;
                se_high_sm(1) = 0;
                
                LSE_low_sm  = mean_low_sm - se_low_sm;
                USE_low_sm  = mean_low_sm + se_low_sm;
                LSE_high_sm = mean_high_sm - se_high_sm;
                USE_high_sm = mean_high_sm + se_high_sm;
                
                patch([x_plot, fliplr(x_plot)], [USE_low_sm, fliplr(LSE_low_sm)], group_colors(1, :), 'FaceAlpha', 0.25, 'LineStyle', 'none');
                patch([x_plot, fliplr(x_plot)], [USE_high_sm, fliplr(LSE_high_sm)], group_colors(2, :), 'FaceAlpha', 0.25, 'LineStyle', 'none');
                p3 = plot(x_plot, mean_low_sm,  'Color', group_colors(1, :), 'LineWidth', 2);
                p4 = plot(x_plot, mean_high_sm, 'Color', group_colors(2, :), 'LineWidth', 2);
                % 
                % boot_output.survival.(field_name).simon_makuch.time         = x_plot;
                % boot_output.survival.(field_name).simon_makuch.S_low        = mean_low_sm;
                % boot_output.survival.(field_name).simon_makuch.S_high       = mean_high_sm;
                % boot_output.survival.(field_name).simon_makuch.SE_low       = se_low_sm;
                % boot_output.survival.(field_name).simon_makuch.SE_high      = se_high_sm;
                % boot_output.survival.(field_name).simon_makuch.sess_S_low   = sess_S_low;
                % boot_output.survival.(field_name).simon_makuch.sess_S_high  = sess_S_high;
                % boot_output.survival.(field_name).simon_makuch.n_sessions   = nSess;
                
                title_suffix = 'Mean \pm SE across sessions (within-session percentiles)';
            else
                [S_low_sm, S_high_sm] = local_simon_makuch_grouped(T_time, C, rowGroup, binEdges);
                S_low_sm  = [1, S_low_sm];
                S_high_sm = [1, S_high_sm];

                if bootstrap
                    S_low_boot_sm  = zeros(nBoot, nBins);
                    S_high_boot_sm = zeros(nBoot, nBins);
                    parfor iBoot = 1:nBoot
                        s = RandStream('mrg32k3a', 'Seed', iBoot);
                        sampled_up_indices = datasample(s, 1:nUPs, nUPs, 'Replace', true);
                        boot_row_indices = [];
                        for u = 1:nUPs
                            targetUP = uniqueUPs(sampled_up_indices(u));
                            boot_row_indices = [boot_row_indices; find(upID == targetUP)]; %#ok<AGROW>
                        end
                        [sl, sh] = local_simon_makuch_grouped(T_time(boot_row_indices, :), C(boot_row_indices), rowGroup(boot_row_indices), binEdges);
                        S_low_boot_sm(iBoot, :)  = sl;
                        S_high_boot_sm(iBoot, :) = sh;
                    end
                    LCI_low_sm  = [1, prctile(S_low_boot_sm,  2.5, 1)];
                    UCI_low_sm  = [1, prctile(S_low_boot_sm,  97.5, 1)];
                    LCI_high_sm = [1, prctile(S_high_boot_sm, 2.5, 1)];
                    UCI_high_sm = [1, prctile(S_high_boot_sm, 97.5, 1)];

                    patch([x_plot, fliplr(x_plot)], [UCI_low_sm, fliplr(LCI_low_sm)], group_colors(1, :), 'FaceAlpha', 0.25, 'LineStyle', 'none');
                    patch([x_plot, fliplr(x_plot)], [UCI_high_sm, fliplr(LCI_high_sm)], group_colors(2, :), 'FaceAlpha', 0.25, 'LineStyle', 'none');
                    % 
                    % boot_output.survival.(field_name).simon_makuch.CI_low  = [LCI_low_sm; UCI_low_sm];
                    % boot_output.survival.(field_name).simon_makuch.CI_high = [LCI_high_sm; UCI_high_sm];
                end
                p3 = plot(x_plot, S_low_sm,  'Color', group_colors(1, :), 'LineWidth', 2);
                p4 = plot(x_plot, S_high_sm, 'Color', group_colors(2, :), 'LineWidth', 2);

                % boot_output.survival.(field_name).simon_makuch.time    = x_plot;
                % boot_output.survival.(field_name).simon_makuch.S_low   = S_low_sm;
                % boot_output.survival.(field_name).simon_makuch.S_high  = S_high_sm;
                title_suffix = 'time-dependent covariate group';
            end

            legend([p3, p4], legend_labels, 'Box', 'off', 'Location', 'southwest');
            title({['Simon-Makuch (unadjusted): ' currFeatLabel], title_suffix});
            ylabel('Survival probability of UP'); xlabel('Time (s)');
            xlim([0, max_time]); ylim([0, 1.05]);
            set(gca, 'TickDir', 'out', 'Box', 'off', 'Color', 'none', 'FontSize', 12);

            % boot_output.survival.(field_name).low_val  = low_val;
            % boot_output.survival.(field_name).high_val = high_val;
        else
            warning('Feature "%s" has no separation between Low/High percentile thresholds; skipping survival curve plots.', currFeat);
        end
    end
end
end

end

function [Lambda0, stratWeight] = local_breslow_grouped(T_time, C, S, lp, binEdges)
% Grouped-time (life-table) Breslow cumulative baseline hazard per stratum,
% as a step function evaluated at binEdges(2:end).
starts = T_time(:, 1);
stops  = T_time(:, 2);
stratList = unique(S);
nStrat = numel(stratList);
nBins  = numel(binEdges) - 1;
Lambda0 = zeros(nStrat, nBins);
stratWeight = zeros(nStrat, 1);

for si = 1:nStrat
    sMask = (S == stratList(si));
    stratWeight(si) = sum(sMask & C == 0);
    cumHaz = 0;
    for k = 1:nBins
        t0 = binEdges(k);
        t1 = binEdges(k + 1);
        riskMask = sMask & (starts < t0) & (stops >= t0);
        d = sum(sMask & C == 0 & stops >= t0 & stops < t1);
        denom = sum(exp(lp(riskMask)));
        if denom > 0
            cumHaz = cumHaz + d / denom;
        end
        Lambda0(si, k) = cumHaz;
    end
end
end

function [S_low, S_high] = local_simon_makuch_grouped(T_time_s, C_s, rowGroup_s, binEdges)
% Grouped-time (life-table) Simon-Makuch product-limit estimator: at each
% bin, risk-set and event-group membership come from whichever interval is
% active at the start of the bin, so a UP event can contribute to different
% groups' risk sets across its own duration.
starts = T_time_s(:, 1);
stops  = T_time_s(:, 2);
nBins  = numel(binEdges) - 1;
S_low  = nan(1, nBins);
S_high = nan(1, nBins);
sLow  = 1;
sHigh = 1;
hasLow = false;
hasHigh = false;

for k = 1:nBins
    t0 = binEdges(k);
    t1 = binEdges(k + 1);
    activeMask = (starts < t0) & (stops >= t0);
    riskLow  = sum(activeMask & rowGroup_s == 0);
    riskHigh = sum(activeMask & rowGroup_s == 1);
    dEventMask = (C_s == 0) & (stops >= t0) & (stops < t1);
    dLow  = sum(dEventMask & rowGroup_s == 0);
    dHigh = sum(dEventMask & rowGroup_s == 1);
    
    if riskLow > 0
        sLow = sLow * (1 - dLow / riskLow);
        hasLow = true;
    end
    if hasLow
        S_low(k) = sLow;
    end
    
    if riskHigh > 0
        sHigh = sHigh * (1 - dHigh / riskHigh);
        hasHigh = true;
    end
    if hasHigh
        S_high(k) = sHigh;
    end
end
end
