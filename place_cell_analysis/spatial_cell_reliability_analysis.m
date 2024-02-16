function  [bestReliability bestReliability_shuffled] = spatial_cell_reliability_analysis(spikeTimes,Behaviour,position_shuffled,Task_info,x_window,x_bin_width,shuffle_option)
    
windowSize = 0.20; % ms
% gaussianWidth = 5; % arbitrary, should be optimized
reliabilityThreshold = 0.01;
numFolds = 10;
windowWidth = [2 5 10 15 20]; % gaussian windows for smoothing rate map


% Spike position when speed above 5
Behaviour.position(abs(Behaviour.speed)<5)=nan;
Behaviour.position(abs(Behaviour.speed)>100)=nan;
% discretePositions = discretize(spikePositions, x_window(1):x_bin_width:x_window(2));

% Behaviour variables
timevec = Behaviour.tvec';
timevec_edge = (timevec(1)-(timevec(2)-timevec(1))/2....
    :mean(diff(timevec)):...
    timevec(end)+(timevec(end)-timevec(end-1))/2)';

bestReliability = calculate_ratemap_explained_variance(spikeTimes,Behaviour,x_window,x_bin_width,numFolds,windowWidth);


% Initialize reliability for shuffles
bestReliability_shuffled = zeros(1000,2,numFolds);
if shuffle_option == 1
    tic
    parfor nshuffle = 1:1000
        Behaviour_temp = Behaviour;
        Behaviour_temp.position = position_shuffled{nshuffle};

        bestReliability_shuffled(nshuffle,:,:) = calculate_ratemap_explained_variance(spikeTimes,Behaviour_temp,x_window,x_bin_width,numFolds,windowWidth)
    end
    toc
end


function bestReliability = calculate_ratemap_explained_variance(spikeTimes,Behaviour,x_window,x_bin_width,numFolds,windowWidth)

% Initialize reliability
reliability = zeros(1, numFolds);
bestReliability = nan(2, numFolds);

for track_id = 1:max(Behaviour.track_ID)
    % Cross-validation
    trackPosition = x_bin_width*discretize(Behaviour.position(Behaviour.track_ID == track_id),x_window(1):x_bin_width:x_window(2));
    trackTime = Behaviour.tvec(Behaviour.track_ID == track_id);

    %     spikePositions = interp1(Behaviour.tvec(Behaviour.track_ID == track_id),Behaviour.position(Behaviour.track_ID == track_id),spikeTimes,'nearest');
    %     spikeCount = histcounts(spikeTimes, timevec_edge)';
    %     spikeCount = spikeCount((Behaviour.track_ID == track_id))';

    % 10 fold cross validation in chunks
    foldSize = floor(length(trackTime) / numFolds);

    % Create a cell array to hold the indices for each fold
    foldIndices = cell(numFolds, 1);

    % Fill the cell array with the indices for each fold
    for i = 1:numFolds
        if i ~= numFolds
            foldIndices{i} = (1+(i-1)*foldSize):(i*foldSize);
        else
            % If it's the last fold, include the remaining data points
            foldIndices{i} = (1+(i-1)*foldSize):length(trackTime);
        end
    end
    
    for nWin = 1:length( windowWidth)
        for i = 1:numFolds
            testIndices = foldIndices{i};
            trainIndices = setdiff(1:length(trackTime), testIndices);

            % Training data
            spikePositionsTrain = interp1(trackTime(trainIndices),trackPosition(trainIndices),spikeTimes,'nearest');
            spikeTimesTrain = interp1(trackTime(trainIndices),trackTime(trainIndices),spikeTimes,'nearest');
            trackPositionTrain = trackPosition(trainIndices);

            % Test data
%             spikePositionsTest = interp1(trackTime(testIndices),trackPosition(testIndices),spikeTimes,'nearest');
            spikeTimesTest = interp1(trackTime(testIndices),trackTime(testIndices),spikeTimes,'nearest');
            trackPositionTest = trackPosition(testIndices);

            % Calculate spike count map and occupancy map
            spikeCountMap = histcounts(spikePositionsTrain, x_window(1):x_bin_width:x_window(2));
            occupancyMap = mean(diff(trackTime))*histcounts(trackPositionTrain,x_window(1):x_bin_width:x_window(2));

            % Define Gaussian window for smoothing
            gaussianWindow = gausswin(windowWidth(nWin));

            % Normalize to have an area of 1 (i.e., to be a probability distribution)
            gaussianWindow = gaussianWindow / sum(gaussianWindow);

            smoothedSpikeCountMap = conv(spikeCountMap, gaussianWindow, 'same');
            smoothedOccupancyMap = conv(occupancyMap, gaussianWindow, 'same');

            % Calculate response profile
            responseProfile = smoothedSpikeCountMap ./ smoothedOccupancyMap;

            % Calculate reliability on test data
%             trackTimeTest = trackTime;
%             trackTimeTest(cv.test(i)) = nan;
            trackTimeTest = trackTime(testIndices);
            y = histcounts(spikeTimesTest, trackTimeTest)./(diff(trackTimeTest));
            yHat = interp1(1:length(responseProfile), responseProfile, trackPositionTest(1:end-1)/x_bin_width); % predicted firing rate
            yHat(isnan(yHat))=0;
            mu = mean(histcounts(spikeTimesTrain,  trackTime(trainIndices))./(diff(trackTime(trainIndices))),'omitnan'); % mean firing rate

            reliability(i) = 1 - sum((y - yHat).^2) / sum((y - mu).^2); % percentage variance explained
            
            if sum(isnan(bestReliability(track_id,i))) == 1
                bestReliability(track_id,i) = reliability(i);
            elseif mean(reliability) > mean(bestReliability(track_id,:))
                bestReliability(track_id,:) = reliability;
%                 nWin
            end

        end
    end
end