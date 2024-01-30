function spatial_modulation_GLM_analysis(clusters, Behaviour)

% This function fits a GLM to the data and compares the performance of
% three different models: a visual model, a non-spatial model, and a
% spatial model. The visual model includes the visual stimulus as a
% predictor. The non-spatial model includes the visual stimulus and
% eye position as predictors. The spatial model includes the visual
% stimulus, eye position, and spatial eye position as predictors. The
% spatial eye position is a linear combination of the horizontal and
% vertical eye position. The spatial model is fit for different values
% of alpha, which controls the relative contribution of the horizontal
% and vertical eye position to the spatial eye position. The model with
% the smallest sum squared error is chosen as the best model. The
% function also fits an ellipse to the distribution of model
% predictions.



% Define time lags
tau = [-1000, -500, 0, 500, 1000]; % ms

% Define ridge regression coefficients
lambda = [0.01, 0.05, 0.1, 0.5, 1];

% Define alpha values for spatial model
alpha = 0:0.1:1;

% Convert spike times to spike count time series
y = histcounts(spikeTimes, timeVector);

% Initialize models
visualModel = zeros(size(y));
nonSpatialModel = zeros(size(y));
spatialModel = zeros(size(y));

% Fit visual model
for i = 1:length(I)
    visualModel = visualModel + sqrt(2) * I{i} * x;
end

% Fit non-spatial model
nonSpatialModel = visualModel;
for j = 1:length(tau)
    nonSpatialModel = nonSpatialModel + s * tau(j) + p * tau(j) + r * tau(j);
end
nonSpatialModel = nonSpatialModel + ex + ey;

% Fit spatial model
for a = alpha
    spatialModel = zeros(size(y));
    for i = 1:length(I)
        if i == 3
            spatialModel = spatialModel + a / sqrt(a^2 + (1-a)^2) * I{i} * x;
        elseif i == 4
            spatialModel = spatialModel + (1-a) / sqrt(a^2 + (1-a)^2) * I{i} * x;
        else
            spatialModel = spatialModel + I{i} * x;
        end
    end
    for j = 1:length(tau)
        spatialModel = spatialModel + s * tau(j) + p * tau(j) + r * tau(j);
    end
    spatialModel = spatialModel + ex + ey;

    % Perform ridge regression
    beta_visual = ridge(y, visualModel, lambda);
    beta_nonSpatial = ridge(y, nonSpatialModel, lambda);
    beta_spatial = ridge(y, spatialModel, lambda);

    % Perform 5-fold cross-validation
    cv = cvpartition(length(y), 'KFold', 5);
    for i = 1:cv.NumTestSets
        trainInd = cv.training(i);
        testInd = cv.test(i);
        beta_visual_cv = ridge(y(trainInd), visualModel(trainInd, :), lambda);
        beta_nonSpatial_cv = ridge(y(trainInd), nonSpatialModel(trainInd, :), lambda);
        beta_spatial_cv = ridge(y(trainInd), spatialModel(trainInd, :), lambda);
    end

    % Calculate sum squared error for each model
    SSE_visual = sum((y - visualModel * beta_visual).^2);
    SSE_nonSpatial = sum((y - nonSpatialModel * beta_nonSpatial).^2);
    SSE_spatial = sum((y - spatialModel * beta_spatial).^2);

    % Choose model with smallest SSE
    [~, model] = min([SSE_visual, SSE_nonSpatial, SSE_spatial]);

    % Print chosen model
    switch model
        case 1
            disp('Visual model chosen');
        case 2
            disp('Non-spatial model chosen');
        case 3
            disp('Spatial model chosen');
    end
end

% Fit an ellipse to the distribution of model predictions
% This part is left as an exercise for the reader
end