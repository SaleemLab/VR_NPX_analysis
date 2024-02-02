function spatial_modulation_GLM_analysis(clusters,Behaviour,Task_info)
% Main function for spatial modulatiom GLM analysis

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

% Define Gaussian window for smoothing
windowWidth = 0.2; % seconds
gaussianWindow = gausswin(windowWidth/mean(diff(timevec)));

% Normalize to have an area of 1 (i.e., to be a probability distribution)
gaussianWindow = gaussianWindow / sum(gaussianWindow);


% Behaviour variables
timevec = clusters_ks3(nprobe).timevec';
timevec_edge = (timevec(1)-(timevec(2)-timevec(1))/2....
    :mean(diff(timevec)):...
    timevec(end)+(timevec(end)-timevec(end-1))/2)';

face_motion_energy = Behaviour.face_motion_enegy;
pupil_size = Behaviour.pupil_size;
lick_count_L = Behaviour.lick_count(1,:);
lick_count_R = Behaviour.lick_count(2,:);
speed = Behaviour.speed;
speed(speed>200) = nan;
speed(speed<-200) = nan; % remove extreme speed

reward = histcounts(Task_info.reward_delivery_time,timevec_edge);
reward_index = find(reward > 0);
filler_bins = 0.25*1/mean(diff(timevec)); % 250ms window

for event = 1:length(reward_index)
    reward(reward_index(event)- filler_bins) = 1; % add 1 to 250ms before reward
     reward(reward_index(event) + filler_bins) = 1;% add 1 to 250ms after reward
end


% Scale pupilDiameter and eyePosition to [-1, 1]
pupil_size = 2 * (pupil_size - min(pupil_size)) / (max(pupil_size) - min(pupil_size)) - 1;
% eyePosition = 2 * (eyePosition - min(eyePosition)) / (max(eyePosition) - min(eyePosition)) - 1;

% Scale speed, face motion energy and lick counts to [0, 1]
speed = (speed - min(speed)) / (max(speed) - min(speed));
face_motion_energy = (face_motion_energy- min(face_motion_energy))/(max(face_motion_energy)-min(face_motion_energy));
lick_count_L = (lick_count_L - min(lick_count_L)) / (max(lick_count_L) - min(lick_count_L));
lick_count_R = (lick_count_R - min(lick_count_R)) / (max(lick_count_R) - min(lick_count_R)); % Not sure if it is best way

% Define time lags
tau = [-1,-0.5,0,0.5,1]; % s

% Apply 5 time lags to behavioural variables
pupil_size_lagged = nan(length(timevec),5);
speed_lagged = nan(length(timevec),5);
reward_lagged = nan(length(timevec),5);
lick_L_lagged = nan(length(timevec),5);
lick_R_lagged = nan(length(timevec),5);
face_energy_lagged = nan(length(timevec),5);

for t = 1:length(tau)
    lag_frames = fix(tau(t)/mean(diff(timevec)));
    pupil_size_lagged(:,t) = create_lag(pupil_size, lag_frames);
    speed_lagged(:, t) = create_lag(speed, lag_frames); 
    reward_lagged(:,t) = create_lag(reward, lag_frames); 
    lick_L_lagged(:,t) = create_lag(lick_count_L, lag_frames);
    lick_R_lagged(:,t) = create_lag(lick_count_R, lag_frames);
    face_energy_lagged(:,t) = create_lag(face_motion_energy, lag_frames);
end

% Define ridge regression coefficients
lambda = [0.01, 0.05, 0.1, 0.5, 1];

% Define alpha values for spatial model
alpha = 0:0.1:1;

% Define spatial basis functions
% I = @(x, i, a, b) i == 1 * (x >= 0 & x < 5) + ...
%                   i == 2 * (x >= 5 & x < 10) + ...
%                   i >= 3 & i <= 10 * ((a * (x >= (i-2)*10 & x < (i-1)*10)) + (b * (x >= (i+3)*10 & x < (i+4)*10))) + ...
%                   i == 11 * (x >= 90 & x < 95) + ...
%                   i == 12 * (x >= 95 & x <= 100);

% Define spatial basis functions
% I = @(x, i, a, b) (i <= 7) .* (x >= (i-1)*5 & x <= i*5) + ... % Before visually repeating segments
%                   (i >= 8 & i <= 13) .* ((a * (x >= (i-1)*5 & x < (i)*5)) + (b * (x >= (i+7)*5 & x < (i+8)*10))) + ... % Visually repeating segments
%                   (i >= 14 & i <= 15) .* (x >= (i-1)*5 & x <= i*5) + ...   % segment between two repeating segments
%                   (i >= 16 & i <= 22) .* (x >= (i+5)*5 & x <= (i+6)*5); % after visually repeating segments

I = @(x, i, a, b) (i <= 6) .* (x > (i-1)*5 & x <= i*5) + ... % Before visually repeating segments
    (i >= 7 & i <= 14) .* ((a * (x > (i-1)*5 & x <= (i)*5)) + (b * (x > (i+7)*5 & x <= (i+8)*5))) + ... % Visually repeating segments
    (i >= 15 & i <= 20) .* (x > (i+7)*5 & x <= (i+8)*5); % after visually repeating segments

fig = figure;
fig.Position = [680,259,1111,719];
count= 1;
% for ncell = 1:length(good_unit)
for ncell = 1:4:length(good_unit)
    % Convert spike times to spike count time series
    SpikeTimes = clusters_ks3(nprobe).spike_times(clusters.spike_id == clusters.cluster_id(ncell));

    % Convolve spike count time series with Gaussian window
    y = histcounts(SpikeTimes, timevec_edge)';
    y = conv(y, gaussianWindow, 'same');

    for track_id = 1:max(Behaviour.track_ID)
        % Position data (speed filtered and specific to one track)
        x = Behaviour.position';
        x(Behaviour.speed < 5) = nan;
        x(Behaviour.track_ID ~= track_id) = nan;
        %     TimeVector(isnan(x)) = [];
        %     y(isnan(x)) = [];
        %     x(isnan(x)) = [];


        % Initialize models and model performance
        BestVisualModelError = [];
        BestVisualModelVar = [];

        BestSpatialModelError = [];
        BestSpatialModelVar = [];

        BestNonSpatialModelError = [];
        BestNonSpatialModelVar = [];

        BestSpatialModel = [];
        BestNonSpatialModel = [];
        BestVisualModel = [];

        parfor nshuffle = 1:1000
            s = RandStream('mrg32k3a','Seed',nshuffle+10000*track_id); % Set random seed for resampling
            % Fit models
            for l = 1:length(lambda)
                % Initialize models
                visualModel = zeros(size(y));
                nonSpatialModel = zeros(size(y));
                spatialModel = zeros(size(y));

                %%%%% Fit non-spatial model
                for t = 1:length(tau)
                    X = [pupil_size_lagged(:,t), speed_lagged(:,t), reward_lagged(:,t),...
                        lick_L_lagged(:,t), lick_R_lagged(:,t),face_energy_lagged(:,t)];
                    beta = (X'*X + lambda(l)*eye(size(X,2))) \ X'*y;
                    nonSpatialModel = nonSpatialModel + beta * X;
                    beta_nonSpatial = ridge(y, nonSpatialModel, lambda(l));

                    % Perform 5-fold cross-validation
                    cv = cvpartition(length(y), 'KFold', 5);
                    for i = 1:cv.NumTestSets
                        trainInd = cv.training(i);
                        testInd = cv.test(i);
                        beta_nonSpatial_cv = ridge(y(trainInd), nonSpatialModel(trainInd, :), lambda(l));
                    end
                end

                for a = 1:length(alpha)
                    for i = 1:20
                        %%%%% Fit visual model
                        X = I(x, i, sqrt(2), sqrt(2));
                        beta = (X'*X + lambda(l)*eye(size(X,2))) \ X'*y;
                        visualModel = visualModel + beta * X;

                        %%%%% Fit spatial model
                        X = I(x, i, alpha(a)/sqrt(alpha(a)^2+(1-alpha(a))^2), (1-alpha(a))/sqrt(alpha(a)^2+(1-alpha(a))^2));
                        beta = (X'*X + lambda(l)*eye(size(X,2))) \ X'*y;

                        %  beta_hat = (X'*X + lambda*eye(size(X,2))) \ X'*y;

                        spatialModel = spatialModel + beta * X;
                    end

                    % Perform ridge regression
                    beta_visual = ridge(y, visualModel, lambda(l));
                    beta_spatial = ridge(y, spatialModel, lambda(l));

                    % Perform 5-fold cross-validation
                    cv = cvpartition(length(y), 'KFold', 5);
                    for i = 1:cv.NumTestSets
                        trainInd = cv.training(i);
                        testInd = cv.test(i);
                        beta_visual_cv = ridge(y(trainInd), visualModel(trainInd, :), lambda(l));
                        beta_nonSpatial_cv = ridge(y(trainInd), nonSpatialModel(trainInd, :), lambda(l));
                        beta_spatial_cv = ridge(y(trainInd), spatialModel(trainInd, :), lambda(l));
                    end

                    % Calculate sum squared error for each model
                    SSE_visual = sum((y - visualModel * beta_visual).^2);
                    %         SSE_nonSpatial = sum((y - nonSpatialModel * beta_nonSpatial).^2);
                    SSE_spatial = sum((y - spatialModel * beta_spatial).^2);
                    SST = sum((y - mean(y)).^2);

                    if isempty(BestVisualModelError)
                        BestVisualModelError = SSE_visual;
                        BestVisualModel = visualModel;
                    else
                        if BestVisualModelError > SSE_visual
                            BestVisualModelError = SSE_visual;
                            BestVisualModel = visualModel;
                        end
                    end

                    %         if spatialModelError > sum((y - spatialModel).^2)
                    %             BestSpatialModelError = sum((y - spatialModel).^2);
                    %         end
                    if isempty(BestSpatialModelError)
                        BestSpatialModelError = SSE_spatial;
                        BestSpatialModel = spatialModel;
                        BestAlpha = alpha(a);
                    else
                        if BestSpatialModelError > SSE_spatial
                            BestSpatialModelError = SSE_spatial;
                            BestSpatialModel = spatialModel;
                            BestAlpha = alpha(a);
                        end
                    end
                end
            end
        end

        % Fit models
        for l = 1:length(lambda)
            % Initialize models
            visualModel = zeros(size(y));
            nonSpatialModel = zeros(size(y));
            spatialModel = zeros(size(y));

            %%%%% Fit non-spatial model
            for t = 1:length(tau)
                X = [pupil_size_lagged(:,t), speed_lagged(:,t), reward_lagged(:,t),...
                    lick_L_lagged(:,t), lick_R_lagged(:,t),face_energy_lagged(:,t)];
                beta = (X'*X + lambda(l)*eye(size(X,2))) \ X'*y;
                nonSpatialModel = nonSpatialModel + beta * X;
                beta_nonSpatial = ridge(y, nonSpatialModel, lambda(l));

                % Perform 5-fold cross-validation
                cv = cvpartition(length(y), 'KFold', 5);
                for i = 1:cv.NumTestSets
                    trainInd = cv.training(i);
                    testInd = cv.test(i);
                    beta_nonSpatial_cv = ridge(y(trainInd), nonSpatialModel(trainInd, :), lambda(l));
                end
            end

            for a = 1:length(alpha)
                for i = 1:20
                    %%%%% Fit visual model
                    X = I(x, i, sqrt(2), sqrt(2));
                    beta = (X'*X + lambda(l)*eye(size(X,2))) \ X'*y;
                    visualModel = visualModel + beta * X;

                    %%%%% Fit spatial model
                    X = I(x, i, alpha(a)/sqrt(alpha(a)^2+(1-alpha(a))^2), (1-alpha(a))/sqrt(alpha(a)^2+(1-alpha(a))^2));
                    beta = (X'*X + lambda(l)*eye(size(X,2))) \ X'*y;

                    %  beta_hat = (X'*X + lambda*eye(size(X,2))) \ X'*y;

                    spatialModel = spatialModel + beta * X;
                end

                % Perform ridge regression
                beta_visual = ridge(y, visualModel, lambda(l));
                beta_spatial = ridge(y, spatialModel, lambda(l));

                % Perform 5-fold cross-validation
                cv = cvpartition(length(y), 'KFold', 5);
                for i = 1:cv.NumTestSets
                    trainInd = cv.training(i);
                    testInd = cv.test(i);
                    beta_visual_cv = ridge(y(trainInd), visualModel(trainInd, :), lambda(l));
                    beta_nonSpatial_cv = ridge(y(trainInd), nonSpatialModel(trainInd, :), lambda(l));
                    beta_spatial_cv = ridge(y(trainInd), spatialModel(trainInd, :), lambda(l));
                end

                % Calculate sum squared error for each model
                SSE_visual = sum((y - visualModel * beta_visual).^2);
                %         SSE_nonSpatial = sum((y - nonSpatialModel * beta_nonSpatial).^2);
                SSE_spatial = sum((y - spatialModel * beta_spatial).^2);
                SST = sum((y - mean(y)).^2);

                if isempty(BestVisualModelError)
                    BestVisualModelError = SSE_visual;
                    BestVisualModel = visualModel;
                else
                    if BestVisualModelError > SSE_visual
                        BestVisualModelError = SSE_visual;
                        BestVisualModel = visualModel;
                    end
                end

                %         if spatialModelError > sum((y - spatialModel).^2)
                %             BestSpatialModelError = sum((y - spatialModel).^2);
                %         end
                if isempty(BestSpatialModelError)
                    BestSpatialModelError = SSE_spatial;
                    BestSpatialModel = spatialModel;
                    BestAlpha = alpha(a);
                else
                    if BestSpatialModelError > SSE_spatial
                        BestSpatialModelError = SSE_spatial;
                        BestSpatialModel = spatialModel;
                        BestAlpha = alpha(a);
                    end
                end
            end
        end

        %     % Fit ellipse to the distribution of model predictions
        %     %     covMatrix = cov([visualModel, nonSpatialModel, spatialModel]);
        %     covMatrix = cov([visualModel, spatialModel]);
        %     [eigVec, eigVal] = eig(covMatrix);
        %     % Find the index of the largest eigenvalue
        %     [~, largestEigInd] = max(diag(eigVal));
        %
        %     % Get the largest eigenvector (corresponding to the major axis of the ellipse)
        %     majorAxis = eigVec(:, largestEigInd);
        %     angle = rad2deg(atan2(majorAxis(2), majorAxis(2)));

        spatial_modulation(nprobe).(ncell)

        if count >= 5
            fig = figure
            fig.Position = [680,259,1111,719];
            count= 1;
        end
        subplot(2,2,count)
        plot(timevec,BestSpatialModel./max(BestSpatialModel),'r');hold on
        plot(timevec,x./max(x),'m');hold on
        plot(timevec,y./max(y),'k');
        legend('spatial model','X','Cell spiking activity')
        count = 1 + count;
        title(sprintf('Cell %i %ium Best alpha %.2f SEM %.2f ',good_unit(ncell),clusters_ks3(nprobe).peak_depth(good_unit(ncell)),BestAlpha,BestSpatialModelError))

        subplot(2,2,count)
        %     plot(TimeVector,BestSpatialModel./max(BestSpatialModel),'r');hold on
        plot(timevec,BestVisualModel./max(BestVisualModel),'b');hold on
        plot(timevec,x./max(x),'m');hold on
        plot(timevec,y./max(y),'k');
        legend('visual model','X','Cell spiking activity')
        count = 1 + count;
        title(sprintf('Cell %i %ium Best alpha %.2f SEM %.2f ',good_unit(ncell),clusters_ks3(nprobe).peak_depth(good_unit(ncell)),BestAlpha,BestVisualModelError))
    end
end

