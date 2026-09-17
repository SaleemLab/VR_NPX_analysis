function [model, h_fig, h_scatter, c_vals] = plot_2D_scatter(x, y, varargin)
% plot_2D_scatter Creates a 2D scatter plot where point colors reflect 
% local point density (normalized 2D kernel/histogram density), or a specified 3rd variable.
% Supports Linear Mixed-Effects (LME) statistical modeling for p-values.
% Supports native log, log2, and log10 scaling with transformed regression fits.
%
% Inputs:
%   x - Vector of X-coordinates
%   y - Vector of Y-coordinates
%
% Optional Name-Value Pairs:
%   'xscale'                  - X-axis scale: 'linear', 'log', 'log2', 'log10' (default: 'linear')
%   'yscale'                  - Y-axis scale: 'linear', 'log', 'log2', 'log10' (default: 'linear')
%   'z' or 'color_var'        - Vector of 3rd variable values to color points (default: [])
%   'use_mixed_effects'      - Logical, whether to fit LME for p-value (default: false)
%   'subject_id'/'animal_id' - Vector of subject/animal IDs for random intercept
%   'session_id'              - Vector of session IDs for random intercept
%   'normalize_density'      - Logical, whether to normalize density to [0, 1] (default: true)
%   'xlabel_text'             - X-axis label string (default: 'X')
%   'ylabel_text'             - Y-axis label string (default: 'Y')
%   'clabel_text'             - Colorbar label string (default: 'Normalized Density' or '3rd Variable')
%   'title_text'              - Plot title string (default: '')
%   'marker_size'             - Size of scatter markers (default: 18)
%   'marker'                  - Marker type string (default: 'o')
%   'alpha'                   - Marker transparency alpha in [0, 1] (default: 0.7)
%   'colormap_name'           - Colormap string ('turbo', 'parula', 'viridis', etc., default: 'turbo')
%   'add_fit_line'            - Logical, whether to add linear trend line (default: false)
%   'show_stats'              - Logical, whether to display statistical metrics (default: false)
%   'nbins'                   - Number of grid bins for 2D density estimation (default: 200)
%   'parent_ax'               - Target axes handle. If empty, creates new figure (default: [])
%   'title_name'              - Title name for figure window
%
% Outputs:
%   model     - Fitted LME or LinearModel object
%   h_fig     - Figure handle (or parent axes handle if provided)
%   h_scatter - Scatter graphics object handle
%   c_vals    - Calculated density or 3rd variable color values for each point

%% 1. Parse Input Arguments
p = inputParser;
addRequired(p, 'x', @isnumeric);
addRequired(p, 'y', @isnumeric);
addParameter(p, 'xscale', 'linear', @(s) any(validatestring(s, {'linear', 'log', 'log2', 'log10'})));
addParameter(p, 'yscale', 'linear', @(s) any(validatestring(s, {'linear', 'log', 'log2', 'log10'})));
addParameter(p, 'z', [], @isnumeric);
addParameter(p, 'color_var', [], @isnumeric);
addParameter(p, 'use_mixed_effects', false, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'subject_id', []);
addParameter(p, 'animal_id', []);
addParameter(p, 'session_id', []);
addParameter(p, 'smooth_color', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'smooth_sigma', 2.0, @isnumeric);
addParameter(p, 'normalize_density', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'xlabel_text', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'ylabel_text', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'clabel_text', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'title_text', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'marker_size', 18, @isnumeric);
addParameter(p, 'marker', 'o', @(x) ischar(x) || isstring(x));
addParameter(p, 'alpha', 0.7, @(x) isnumeric(x) && x >= 0 && x <= 1);
addParameter(p, 'colormap_name', 'turbo', @(x) ischar(x) || isstring(x));
addParameter(p, 'add_fit_line', false, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'show_stats', false, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'nbins', 200, @isnumeric);
addParameter(p, 'parent_ax', [], @(x) isempty(x) || isgraphics(x, 'axes'));
addParameter(p, 'title_name', '', @(x) ischar(x) || isstring(x));
parse(p, x, y, varargin{:});

xscale_opt  = lower(char(p.Results.xscale));
yscale_opt  = lower(char(p.Results.yscale));
z_val       = p.Results.z;
if isempty(z_val)
    z_val   = p.Results.color_var;
end
subj_id     = p.Results.subject_id;
if isempty(subj_id)
    subj_id = p.Results.animal_id;
end
sess_id     = p.Results.session_id;
use_lme     = logical(p.Results.use_mixed_effects);
if ~isempty(subj_id) || ~isempty(sess_id)
    use_lme = true;
end
smooth_param = p.Results.smooth_color;
if islogical(smooth_param)
    do_smooth = smooth_param;
    smooth_sigma = p.Results.smooth_sigma;
elseif isnumeric(smooth_param)
    do_smooth = (smooth_param > 0);
    smooth_sigma = smooth_param;
else
    do_smooth = true;
    smooth_sigma = p.Results.smooth_sigma;
end
norm_dens   = logical(p.Results.normalize_density);
xlabel_str  = char(p.Results.xlabel_text);
ylabel_str  = char(p.Results.ylabel_text);
clabel_str  = char(p.Results.clabel_text);
title_str   = char(p.Results.title_text);
sz          = p.Results.marker_size;
mk          = char(p.Results.marker);
alpha_val   = p.Results.alpha;
cmap_str    = char(p.Results.colormap_name);
do_fit      = logical(p.Results.add_fit_line);
do_stats    = logical(p.Results.show_stats);
nbins       = p.Results.nbins;
ax_target   = p.Results.parent_ax;
title_name  = p.Results.title_name;

%% 2. Data Cleaning & Transformations
x = x(:);
y = y(:);
valid_mask = ~isnan(x) & ~isinf(x) & ~isnan(y) & ~isinf(y);

% Enforce positivity for log-transformed dimensions
if ~strcmp(xscale_opt, 'linear')
    valid_mask = valid_mask & (x > 0);
end
if ~strcmp(yscale_opt, 'linear')
    valid_mask = valid_mask & (y > 0);
end

if ~isempty(z_val)
    z_val = z_val(:);
    if length(z_val) ~= length(x)
        error('Length of 3rd variable (z) must match length of x and y.');
    end
    valid_mask = valid_mask & ~isnan(z_val) & ~isinf(z_val);
end
if ~isempty(subj_id)
    subj_id = subj_id(:);
    if length(subj_id) ~= length(x)
        error('Length of subject_id must match length of x and y.');
    end
    if isnumeric(subj_id)
        valid_mask = valid_mask & ~isnan(subj_id) & ~isinf(subj_id);
    elseif iscategorical(subj_id)
        valid_mask = valid_mask & ~isundefined(subj_id);
    elseif iscell(subj_id) || isstring(subj_id)
        valid_mask = valid_mask & ~cellfun(@isempty, cellstr(subj_id));
    end
end
if ~isempty(sess_id)
    sess_id = sess_id(:);
    if length(sess_id) ~= length(x)
        error('Length of session_id must match length of x and y.');
    end
    if isnumeric(sess_id)
        valid_mask = valid_mask & ~isnan(sess_id) & ~isinf(sess_id);
    elseif iscategorical(sess_id)
        valid_mask = valid_mask & ~isundefined(sess_id);
    elseif iscell(sess_id) || isstring(sess_id)
        valid_mask = valid_mask & ~cellfun(@isempty, cellstr(sess_id));
    end
end

x = x(valid_mask);
y = y(valid_mask);
if ~isempty(z_val), z_val = z_val(valid_mask); end
if ~isempty(subj_id), subj_id = subj_id(valid_mask); end
if ~isempty(sess_id), sess_id = sess_id(valid_mask); end

if isempty(x)
    warning('No valid finite data points to plot.');
    model = []; h_fig = []; h_scatter = []; c_vals = [];
    return;
end

% Apply coordinate transformations
x_plot = x;
y_plot = y;
switch xscale_opt
    case 'log',   x_plot = log(x);
    case 'log2',  x_plot = log2(x);
    case 'log10', x_plot = log10(x);
end
switch yscale_opt
    case 'log',   y_plot = log(y);
    case 'log2',  y_plot = log2(y);
    case 'log10', y_plot = log10(y);
end

%% 3. Calculate Color Values (Density vs 3rd Variable)
x_min = min(x_plot); x_max = max(x_plot);
y_min = min(y_plot); y_max = max(y_plot);
dx = (x_max - x_min) * 0.05; if dx == 0, dx = 1e-5; end
dy = (y_max - y_min) * 0.05; if dy == 0, dy = 1e-5; end

x_edges = linspace(x_min - dx, x_max + dx, nbins + 1);
y_edges = linspace(y_min - dy, y_max + dy, nbins + 1);
x_centers = (x_edges(1:end-1) + x_edges(2:end)) / 2;
y_centers = (y_edges(1:end-1) + y_edges(2:end)) / 2;

[~, ~, x_bin] = histcounts(x_plot, x_edges);
[~, ~, y_bin] = histcounts(y_plot, y_edges);
valid_grid = (x_bin > 0 & x_bin <= nbins & y_bin > 0 & y_bin <= nbins);

if smooth_sigma > 0
    k_half = max(3, ceil(3 * smooth_sigma));
    k1d = exp((-k_half:k_half).^2 / (-2 * smooth_sigma^2));
    k2d = (k1d' * k1d); 
    k2d = k2d / sum(k2d(:));
else
    k2d = 1;
end

if isempty(z_val)
    density_grid = accumarray([y_bin(valid_grid), x_bin(valid_grid)], 1, [nbins, nbins]);
    
    if smooth_sigma > 0
        density_smooth = conv2(density_grid, k2d, 'same');
    else
        density_smooth = density_grid;
    end
    
    if do_smooth
        c_vals = interp2(x_centers, y_centers, density_smooth, x_plot, y_plot, 'linear', 0);
        c_vals = max(0, c_vals);
    else
        c_vals = zeros(size(x_plot));
        for i = 1:length(x_plot)
            if valid_grid(i)
                c_vals(i) = density_smooth(y_bin(i), x_bin(i));
            end
        end
    end
    
    if norm_dens
        c_min = min(c_vals);
        c_max = max(c_vals);
        if c_max > c_min
            c_vals = (c_vals - c_min) / (c_max - c_min);
        end
        if isempty(clabel_str)
            clabel_str = 'Normalized Density';
        end
    else
        if isempty(clabel_str)
            clabel_str = 'Point Density';
        end
    end
else
    if do_smooth && smooth_sigma > 0 && ~any(strcmp('smooth_color', p.UsingDefaults))
        z_sum_grid = accumarray([y_bin(valid_grid), x_bin(valid_grid)], z_val(valid_grid), [nbins, nbins]);
        cnt_grid   = accumarray([y_bin(valid_grid), x_bin(valid_grid)], 1, [nbins, nbins]);
        
        z_sum_smooth = conv2(z_sum_grid, k2d, 'same');
        cnt_smooth   = conv2(cnt_grid, k2d, 'same');
        
        z_spatial_grid = z_sum_smooth ./ max(eps, cnt_smooth);
        
        c_vals = interp2(x_centers, y_centers, z_spatial_grid, x_plot, y_plot, 'linear', NaN);
        nan_mask = isnan(c_vals);
        c_vals(nan_mask) = z_val(nan_mask);
    else
        c_vals = z_val;
    end
    
    if isempty(clabel_str)
        clabel_str = '3rd Variable';
    end
end

%% 4. Sort Points by Color Value
[c_vals_sorted, s_idx] = sort(c_vals, 'ascend');
x_sorted = x_plot(s_idx);
y_sorted = y_plot(s_idx);

%% 5. Setup Plot & Axes
if isempty(ax_target)
    h_fig = figure('Color', 'w');
    h_fig.Position = [300, 300, 600, 500];
    ax = gca;
    if ~isempty(title_name)
        h_fig.Name = title_name;
    end
else
    ax = ax_target;
    h_fig = ancestor(ax, 'figure');
end
hold(ax, 'on');

h_scatter = scatter(ax, x_sorted, y_sorted, sz, c_vals_sorted, mk, 'filled', ...
    'MarkerEdgeAlpha', min(1, alpha_val + 0.1), ...
    'MarkerFaceAlpha', alpha_val);

try
    colormap(ax, cmap_str);
catch
    colormap(ax, 'parula');
end
cbar = colorbar(ax);
cbar.Label.String = clabel_str;
cbar.Label.FontSize = 10.5;
try
    cbar.TickDirection = 'out';
catch
end

% Linear trend line fitted on transformed coordinates
if do_fit && length(x_sorted) >= 2
    p_fit = polyfit(x_sorted, y_sorted, 1);
    x_line = linspace(min(x_sorted), max(x_sorted), 200);
    y_line = polyval(p_fit, x_line);
    plot(ax, x_line, y_line, 'k--', 'LineWidth', 1.8, 'DisplayName', 'Linear Fit');
end

% Statistical modeling on transformed coordinates
model = [];
if do_stats && length(x_sorted) >= 3
    if use_lme && (~isempty(subj_id) || ~isempty(sess_id))
        try
            tbl_lme = table(zscore(x_plot), zscore(y_plot), 'VariableNames', {'x', 'y'});
            formula_str = 'y ~ x';
            
            if ~isempty(subj_id) && ~isempty(sess_id)
                tbl_lme.subject_id = categorical(subj_id);
                tbl_lme.session_id = categorical(sess_id);
                formula_str = 'y ~ x + (1|subject_id) + (1|session_id)';
            elseif ~isempty(subj_id)
                tbl_lme.subject_id = categorical(subj_id);
                formula_str = 'y ~ x + (1|subject_id)';
            elseif ~isempty(sess_id)
                tbl_lme.session_id = categorical(sess_id);
                formula_str = 'y ~ x + (1|session_id)';
            end
            
            lme = fitlme(tbl_lme, formula_str);
            beta_lme = lme.Coefficients.Estimate(2);
            p_val    = lme.Coefficients.pValue(2);
            r_sq     = lme.Rsquared.Ordinary;
            model    = lme;
            stats_str = sprintf('LME \\beta = %.2f \nR^2 = %.3f \np = %.2e', ...
                beta_lme, r_sq, p_val);
        catch ME
            warning('LME fit failed (%s). Falling back to Pearson correlation.', ME.message);
            mdl   = fitlm(zscore(x_sorted), zscore(y_sorted));
            model = mdl;
            beta  = mdl.Coefficients.Estimate(2);
            r_sq  = mdl.Rsquared.Ordinary;
            p_val = mdl.ModelFitVsNullModel.Pvalue;
            stats_str = sprintf('\\beta = %.2f \nR^2 = %.3f \np = %.2e', ...
                beta, r_sq, p_val);
        end
    else
        mdl   = fitlm(zscore(x_sorted), zscore(y_sorted));
        model = mdl;
        beta  = mdl.Coefficients.Estimate(2);
        r_sq  = mdl.Rsquared.Ordinary;
        p_val = mdl.ModelFitVsNullModel.Pvalue;
        stats_str = sprintf('\\beta = %.2f \nR^2 = %.3f \np = %.2e', ...
            beta, r_sq, p_val);
    end
    
    xl = xlim(ax); yl = ylim(ax);
    text(ax, xl(1) + 0.04*(xl(2)-xl(1)), yl(2) - 0.08*(yl(2)-yl(1)), stats_str, ...
        'FontSize', 10, 'BackgroundColor', [1 1 1 0.75], 'EdgeColor', [0.7 0.7 0.7]);
end

%% 6. Axes Styling & Tick Formatting
if ~isempty(xlabel_str), xlabel(ax, xlabel_str, 'FontSize', 11); end
if ~isempty(ylabel_str), ylabel(ax, ylabel_str, 'FontSize', 11); end
if ~isempty(title_str),  title(ax, title_str, 'FontSize', 12, 'Interpreter', 'none'); end
set(ax, 'TickDir', 'out', 'Box', 'off', 'Color', 'none', 'FontSize', 12);
grid(ax, 'on');
ax.GridAlpha = 0.15;

% Map tick labels back to exponential representation for transformed axes
% Set clean, round raw-value tick marks mapped onto transformed space
if ~strcmp(xscale_opt, 'linear')
    raw_min = 0.8*min(x);
    raw_max = 1.2*max(x);
    
    % Pick clean, standard whole-number tick values spanning the data range
    all_raw_candidates = [0.01, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000];
    clean_raw_x = all_raw_candidates(all_raw_candidates >= raw_min & all_raw_candidates <= raw_max);
    
    % If data range is too narrow, fall back to linear spacing
    if length(clean_raw_x) < 3
        clean_raw_x = round(linspace(raw_min, raw_max, 5), 1);
    end
    
    % Map the clean raw values to where they lie on the transformed axis
    switch xscale_opt
        case 'log10', tick_pos_x = log10(clean_raw_x);
        case 'log2',  tick_pos_x = log2(clean_raw_x);
        case 'log',   tick_pos_x = log(clean_raw_x);
    end
    
    set(ax, 'XTick', tick_pos_x, 'XTickLabel', clean_raw_x);
end

if ~strcmp(yscale_opt, 'linear')
    raw_min = 0.8*min(y);
    raw_max = 1.2*max(y);
    
    % all_raw_candidates = [0.1, 0.5, 1, 2, 5, 6, 8, 10, 12, 15, 20, 30, 50, 100];
    all_raw_candidates = [0.01, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000];
    clean_raw_y = all_raw_candidates(all_raw_candidates >= raw_min & all_raw_candidates <= raw_max);
    
    if length(clean_raw_y) < 3
        clean_raw_y = round(linspace(raw_min, raw_max, 5));
    end
    
    switch yscale_opt
        case 'log10', tick_pos_y = log10(clean_raw_y);
        case 'log2',  tick_pos_y = log2(clean_raw_y);
        case 'log',   tick_pos_y = log(clean_raw_y);
    end
    
    set(ax, 'YTick', tick_pos_y, 'YTickLabel', clean_raw_y);
end

end