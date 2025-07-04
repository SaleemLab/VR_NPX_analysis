function smart_rasterise_scatter(fig_handle, output_filename, options)
% SMART_RASTERISE_SCATTER
% Rasterizes only scatter dots in selected subplots and exports the rest as vector.
%
% Usage:
%   smart_rasterise_scatter(gcf, 'output.pdf');
%   smart_rasterise_scatter(gcf, 'output.pdf', struct('subplot_idx', [1 2]));

if nargin < 1 || isempty(fig_handle), fig_handle = gcf; end
if nargin < 2 || isempty(output_filename), output_filename = 'rasterised_output.pdf'; end
if nargin < 3, options = struct(); end

% Optionally select specific subplots
if isfield(options, 'subplot_idx')
    raster_idx = options.subplot_idx;
else
    raster_idx = [];  % Auto-detect all with scatter
end

% Find and sort axes visually (top-down, left-right)
ax_all = findall(fig_handle, 'Type', 'Axes');
ax_pos = arrayfun(@(a) a.Position, ax_all, 'UniformOutput', false);
ax_pos = cat(1, ax_pos{:});
[~, sort_idx] = sortrows(ax_pos, [-2, 1]);
ax_all = ax_all(sort_idx);

% Process axes
for i = 1:numel(ax_all)
    ax = ax_all(i);
    scatterList = findall(ax, 'Type', 'Scatter');

    if isempty(scatterList)
        continue;
    end
    if ~isempty(raster_idx) && ~ismember(i, raster_idx)
        continue;
    end

    % Rasterize only the scatter
    rasterise_scatter_via_background_axes(ax, scatterList);
end

% Export final figure
exportgraphics(fig_handle, output_filename, 'ContentType', 'vector');
disp(['Exported to: ', output_filename]);
end

%% ------------------------------------------------------------------------
function rasterise_scatter_via_background_axes(ax, scatterList)
    % Get axis limits and parent figure
    xlim_ = ax.XLim;
    ylim_ = ax.YLim;
    fig = ancestor(ax, 'figure');

    % Ensure normalized position
    original_units = ax.Units;
    ax.Units = 'normalized';
    ax_pos = ax.Position;

    % Create offscreen figure for rasterization
    tempFig = figure('Visible','off', 'Position', [100 100 800 600]);
    tempAx = axes(tempFig, 'Position', [0 0 1 1]);
    hold(tempAx, 'on'); axis(tempAx, 'off');
    xlim(tempAx, xlim_);
    ylim(tempAx, ylim_);
    set(tempAx, 'YDir', 'normal');

    % Copy scatter only
    for s = 1:numel(scatterList)
        copyobj(scatterList(s), tempAx);
    end

    % Export raster
    tmpfile = [tempname, '.png'];
    exportgraphics(tempFig, tmpfile, 'Resolution', 300);
    close(tempFig);
    img = imread(tmpfile);
    delete(tmpfile);

    % Create background axis in original figure
    rasterAx = axes('Parent', fig, ...
        'Position', ax_pos, ...
        'Units', 'normalized', ...
        'Color', 'none', ...
        'XLim', xlim_, ...
        'YLim', ylim_, ...
        'YDir', 'normal', ...
        'Visible', 'off');

    % Insert image aligned to data space (pixel centers)
    nX = size(img, 2);
    nY = size(img, 1);
    dx = diff(xlim_) / nX;
    dy = diff(ylim_) / nY;

    xImg = linspace(xlim_(1) + dx/2, xlim_(2) - dx/2, nX);
    yImg = linspace(ylim_(1) + dy/2, ylim_(2) - dy/2, nY);

    image(rasterAx, 'XData', xImg, 'YData', yImg, 'CData', img);

    % Delete original scatter
    delete(scatterList);

    % Push background axis to bottom
    uistack(rasterAx, 'bottom');

    % Link axes for consistent zoom/pan
    linkaxes([ax rasterAx]);

    % Restore axis units
    ax.Units = original_units;
end
