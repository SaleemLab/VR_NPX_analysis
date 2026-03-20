function plotBarWithCI(eventIndexCell, data, colors, varargin)
% plotBarWithCI plots bar chart with error bars based on indexed events
%
% INPUTS:
%   eventIndexCell : {n x 1} cell array where each cell contains event indices
%   data           : matrix of data (size [1 x nevent] or [nbin x nevent])
%   colors         : n x 3 matrix of RGB values (one per bar)
%   Optional Name-Value pairs:
%       'XNames'   : cell array of x-tick labels
%       'TimeBin'  : integer, specifying which time bin to use (only if nbin > 1)
%
% Example usage:
%   data = randn(5, 300); % 5 timebins x 300 events
%   eventIndexCell = {find(rand(1,300)>0.5), find(rand(1,300)<0.5)};
%   colors = [0.2 0.6 0.8; 0.8 0.3 0.3];
%   plotBarWithError_v2(eventIndexCell, data, colors, 'XNames', {'Group A','Group B'}, 'TimeBin', 3);

% Parse inputs
p = inputParser;
addParameter(p, 'XNames', []);
addParameter(p, 'YNames', []);
addParameter(p, 'TimeBin', []);
parse(p, varargin{:});

xNames = p.Results.XNames;
timeBin = p.Results.TimeBin;

nBars = numel(eventIndexCell);

means = zeros(1, nBars);
lowerError = zeros(1, nBars);
upperError = zeros(1, nBars);

% Determine if data is timeseries
isTimeseries = (size(data,1) > 1);

for i = 1:nBars
    indices = eventIndexCell{i};
    if isTimeseries
        if isempty(timeBin)
            error('Data is time series, but TimeBin was not specified.');
        end
        selectedData = data(timeBin, indices);
    else
        selectedData = data(indices);
    end

    means(i) = mean(selectedData);
    pct = prctile(selectedData, [2.5 97.5]);
    lowerError(i) = means(i) - pct(1);
    upperError(i) = pct(2) - means(i);
end

figure;
hold on;

% Plot bars
for i = 1:nBars
    bar(i, means(i), 'FaceColor', colors(i,:), 'EdgeColor', 'none');
end

% Plot error bars
errorbar(1:nBars, means, lowerError, upperError, ...
    'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize',10);

% Aesthetics
if ~isempty(xNames)
    set(gca, 'XTick', 1:nBars, 'XTickLabel', xNames);
else
    set(gca, 'XTick', 1:nBars);
end

xlim([0 nBars+1]);
box on;
ylabel('Bootstrapped ');
hold off;

end