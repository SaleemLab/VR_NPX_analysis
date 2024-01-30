function  place_cell_reliability_analysis(clusters,Behaviour)
    
% Assuming spikeTimes and positions are given
binSize = 2; % cm
windowSize = 250; % ms
gaussianWidth = 5; % arbitrary, should be optimized
reliabilityThreshold = 0.01;
numFolds = 10;

% Discretize the position
discretePositions = discretize(positions, 0:binSize:max(positions));

% Initialize reliability
reliability = zeros(1, numFolds);

% Cross-validation
cv = cvpartition(length(spikeTimes), 'KFold', numFolds);
for i = 1:cv.NumTestSets
    % Training data
    spikeTimesTrain = spikeTimes(cv.training(i));
    discretePositionsTrain = discretePositions(cv.training(i));
    
    % Test data
    spikeTimesTest = spikeTimes(cv.test(i));
    discretePositionsTest = discretePositions(cv.test(i));
    
    % Calculate spike count map and occupancy map
    spikeCountMap = histcounts(spikeTimesTrain(discretePositionsTrain), 0:windowSize:max(spikeTimesTrain));
    occupancyMap = histcounts(discretePositionsTrain, 0:binSize:max(discretePositionsTrain));
    
    % Smooth maps
    gaussianWindow = fspecial('gaussian', [1, length(spikeCountMap)], gaussianWidth);
    smoothedSpikeCountMap = conv(spikeCountMap, gaussianWindow, 'same');
    smoothedOccupancyMap = conv(occupancyMap, gaussianWindow, 'same');
    
    % Calculate response profile
    responseProfile = smoothedSpikeCountMap ./ smoothedOccupancyMap;
    
    % Calculate reliability on test data
    y = histcounts(spikeTimesTest(discretePositionsTest), 0:windowSize:max(spikeTimesTest)) / windowSize; % firing rate
    yHat = interp1(1:length(responseProfile), responseProfile, discretePositionsTest); % predicted firing rate
    mu = mean(y); % mean firing rate
    
    reliability(i) = 1 - sum((y - yHat).^2) / sum((y - mu).^2);
end

% Average reliability
avgReliability = mean(reliability);

% Filter neurons based on average reliability
if avgReliability > reliabilityThreshold
    % This neuron is considered for further analysis
end