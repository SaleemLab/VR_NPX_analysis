function dimension_reduction_two_tracks(tvec_template,position,speed,lap_times,track_ID,spikes_template,event_times,spike_target)


% Define time bin edges
time_bin_size=0.1;
timevec_edge_RUN= interp1(tvec_template,tvec_template,tvec_template(1):time_bin_size:tvec_template(end));
speed_RUN= interp1(tvec_template,speed,tvec_template(1):time_bin_size:tvec_template(end));
timevec_edge_RUN(speed_RUN<1)=nan;

bins = [];
samples = timevec_edge_RUN';

bins(:,1)=samples;
bins(:,2)=samples+time_bin_size;


position_interp1 = interp1(tvec_template,position,timevec_edge_RUN);
position_interp1(isnan(position_interp1))=0;
position_interp1 = discretize(position_interp1,28)*5;


[T1_bins,~,index] = InIntervals(timevec_edge_RUN',lap_times(track_ID==1,:));
[T2_bins,~,index] = InIntervals(timevec_edge_RUN',lap_times(track_ID==2,:));
[~,~,~,~,~,n] = ActivityTemplates(spikes_template,'bins',bins);


time_bin_size=0.05;
% timevec_edge= interp1(tvec_target,tvec_target,tvec_target(1):time_bin_size:tvec_target(end));
% speed_RUN= interp1(tvec_template,speed,tvec_template(1):time_bin_size:tvec_template(end));
% timevec_edge_RUN(speed_RUN<1)=nan;
% tvec_template

ripples_time_edges=[];
ripples_id=[];
for i = 1: length(event_times)
    event_duration = event_times(i,2) - event_times(i,1);
%     event_duration = event_times(i,2) - event_times(i,1);
    if event_duration <0.1
        event_duration = 0.1;
    end
    num_bins = floor(event_duration / time_bin_size);
    timebins_edges = linspace(event_times(i,1), event_times(i,1) + num_bins *  time_bin_size, num_bins+1);
    ripples_time_edges = [ripples_time_edges timebins_edges];
    ripples_id = [ripples_id i*ones(1,length(timebins_edges))];
end

ripple_bins = [];
ripple_bins(:,1)=ripples_time_edges;
ripple_bins(:,2)=ripples_time_edges+time_bin_size;
[~,~,~,~,~,n_ripples] = ActivityTemplates(spike_target,'bins',ripple_bins);

% Calculate mean activity across all neurons for each ripple
mean_ripple_activity = mean(n_ripples, 2); % Average firing rate per ripple

% Regress out the mean ripple activity
coeff = (mean_ripple_activity' * n_ripples) / (mean_ripple_activity' * mean_ripple_activity);
ripple_global_pattern = mean_ripple_activity * coeff; % Global dynamics component
n_ripples_residuals = zscore(n_ripples - ripple_global_pattern,0,1); % Residuals


% 
% [status,interval,index] = InIntervals(values,event_times)

%%%%% PCA or ICA
[templates,correlations,projection1,weights,variance] = ActivityTemplates(spikes_template,'bins',bins(T1_bins,:),'mode','pca');
%           [templates,correlations,projection1,weights,variance] = ActivityTemplates(spikes_template,'bins',bins,'mode','pca');
%           scatter3(n(T1_bins,:) * weights(:, 1),n(T1_bins,:) * weights(:, 2),n(T1_bins,:) * weights(:, 3),20,position_interp1(T1_bins)'/max(position_interp1),'filled');hold on
%           scatter3(n(T2_bins,:) * weights(:, 1),n(T2_bins,:) * weights(:, 2),n(T2_bins,:) * weights(:, 3),20,position_interp1(T2_bins)'/max(position_interp1),'filled');hold on
%           colorbar
%           colormap(gray)

scatter3(n(T1_bins,:) * weights(:, 1),n(T1_bins,:) * weights(:, 2),n(T1_bins,:) * weights(:, 3),15,'r','filled','MarkerFaceAlpha',0.1);hold on
scatter3(n(T2_bins,:) * weights(:, 1),n(T2_bins,:) * weights(:, 2),n(T2_bins,:) * weights(:, 3),15,'b','filled','MarkerFaceAlpha',0.1);hold on
scatter3(n_ripples* weights(:, 1),n_ripples* weights(:, 2),n_ripples* weights(:, 3),15,'magenta','filled','MarkerFaceAlpha',0.1);hold on



for nshuffle = 1:1000
    s = RandStream('mrg32k3a','Seed',nshuffle); % Set random seed for resampling
    random_cell_index = randperm(s,size(n_ripples,2));
    n_ripples_shuffled = n_ripples(:,random_cell_index);
%     scatter3(n_ripples_shuffled* weights(:, 1),n_ripples_shuffled* weights(:, 2),n_ripples_shuffled* weights(:, 3),15,'k','filled','MarkerFaceAlpha',0.1);hold on

end

% Combine data and labels
X = zscore([n(T1_bins,:); n(T2_bins,:)]);
Y = [ones(sum(T1_bins),1); 2 * ones(sum(T2_bins),1)];
position_Z = [position_interp1(T1_bins) position_interp1(T2_bins)];

% Create 10-fold stratified partition
rng(200)
numFolds = 10;
cv = cvpartition(Y, 'KFold', numFolds);

for fold = 1:numFolds
    tic
    % Get training and test indices for this fold
    trainIdx = training(cv, fold); % Logical vector for training samples
    testIdx = test(cv, fold);     % Logical vector for testing samples


    %%%%% Perform PLS regression
    [XL, YL, XS, YS, beta,PCTVAR] = plsregress(X(trainIdx,:),Y(trainIdx),3); % 3 components

    % Project data onto the PLS components

    [W, ~] = qr(XL, 0); % Orthonormalize XL using QR decomposition
    plsScores_Track1 = n(testIdx & Y == 1,:) * W;
    plsScores_Track2 = n(testIdx & Y == 2,:) * W;

%     plsScores_train = X*XL; % The PLS latent scores (reduced dimensions)
%     svmModel = fitcsvm(plsScores_train, Y, 'KernelFunction', 'linear', 'KernelScale', 'auto', 'BoxConstraint', 1);
%     CVMdl = crossval(svmModel);

    %            plsScores_Sleep = spikeCounts_Sleep * XL; % Project sleep events

    %     % Visualize first PLS component
%     figure;
%     scatter3(plsScores_Track2(:, 1), plsScores_Track2(:, 2), plsScores_Track2(:, 3),10,'b','filled','MarkerFaceAlpha',0.1);  hold on;
%     scatter3(plsScores_Track1(:, 1), plsScores_Track1(:, 2), plsScores_Track1(:, 3),10,'r','filled','MarkerFaceAlpha',0.1);hold on;
%     scatter3(plsScores_ripples(:, 1), plsScores_ripples(:, 2), plsScores_ripples(:, 3),10,'m','filled','MarkerFaceAlpha',0.1);  hold on;

    % Train SVM on PLS scores
    plsScores_train = X(trainIdx,:)*W; % The PLS latent scores (reduced dimensions)
    %     mdl = fitcsvm(plsScores_train, Y(trainIdx), 'KernelFunction', 'linear', 'KernelScale', 'auto', 'BoxConstraint', 1);
    mdl = fitclinear(plsScores_train, Y(trainIdx), 'Learner', 'logistic', 'Regularization', 'ridge');
    plsScores_test = X(testIdx,:)* W;
    predictedTrack = predict(mdl, plsScores_test); % Predict track identity
%     plot(cumsum(predictedTrack-1));hold on;plot(cumsum(Y(testIdx)-1));
    
%     sum(predictedTrack == Y(testIdx))/sum(testIdx)

    temp = position_Z(testIdx);
    correct_position1 = temp(Y(testIdx)==predictedTrack  & Y(testIdx) == 1);
    correct_position2 = temp(Y(testIdx)==predictedTrack  & Y(testIdx) == 2);
    incorrect_position1 = temp(Y(testIdx)~=predictedTrack  & Y(testIdx) == 1);
    incorrect_position2 = temp(Y(testIdx)~=predictedTrack  & Y(testIdx) == 2);

    % Plot 3D PLS scores with labels
%     figure;
%     scatter3(plsScores_test(Y(testIdx)==1, 1), plsScores_test(Y(testIdx)==1, 2), plsScores_test(Y(testIdx)==1, 3), 10,'b','filled','MarkerFaceAlpha',0.1); hold on;% Track 1
%     scatter3(plsScores_test(Y(testIdx)==2, 1), plsScores_test(Y(testIdx)==2, 2), plsScores_test(Y(testIdx)==2, 3), 10,'r','filled','MarkerFaceAlpha',0.1); % Track 2


    PLS.XL{fold} = XL; % Weight for projection
    PLS.trainIdx{fold} = trainIdx;
    PLS.testIdx{fold} = testIdx;
    PLS.correct_position1{fold} = correct_position1;
    PLS.incorrect_position1{fold} = incorrect_position1;
    PLS.correct_position2{fold} = correct_position2;
    PLS.incorrect_position2{fold} = incorrect_position2;
    PLS.hit_rate(fold) = sum(predictedTrack == Y(testIdx))/sum(testIdx);
    PLS.XL{fold} = XL; % The weights for the linear combinations of the predictors ( 𝑋 X) to form the latent components.
    PLS.W{fold} = W; % XL was produced in a way that maximizes covariacne between X and Y, orthonormalized weights for projection.
    PLS.PCTVAR{fold} = PCTVAR; % variance explaiend

    PLS.X = X;
    PLS.Y = Y;
    toc
end


for nshuffle = 1:1000
    PLS_shuffled.correct_position1{nshuffle} = [];
    PLS_shuffled.incorrect_position1{nshuffle} = [];
    PLS_shuffled.correct_position2{nshuffle} = [];
    PLS_shuffled.incorrect_position2{nshuffle} = [];
%     tic
    for fold = 1:numFolds

        % Get training and test indices for this fold
        trainIdx = training(cv, fold); % Logical vector for training samples
        testIdx = test(cv, fold);     % Logical vector for testing samples

        s = RandStream('mrg32k3a','Seed',nshuffle+1000000); % Set random seed for resampling
        random_cell_index = randperm(s,size(X,2));
        X_shuffled = X(:,random_cell_index);

        %%%%% Perform PLS regression
        [XL, YL, XS, YS, beta,PCTVAR] = plsregress(X_shuffled(trainIdx,:),Y(trainIdx),3); % 3 components

        % Project data onto the PLS components

        [W, ~] = qr(XL, 0); % Orthonormalize XL using QR decomposition
        plsScores_Track1 = n(testIdx & Y == 1,:) * W;
        plsScores_Track2 = n(testIdx & Y == 2,:) * W;

        %     plsScores_train = X*XL; % The PLS latent scores (reduced dimensions)
        %     svmModel = fitcsvm(plsScores_train, Y, 'KernelFunction', 'linear', 'KernelScale', 'auto', 'BoxConstraint', 1);
        %     CVMdl = crossval(svmModel);

        %            plsScores_Sleep = spikeCounts_Sleep * XL; % Project sleep events

        %     % Visualize first PLS component
        %     figure;
        %     scatter3(plsScores_Track2(:, 1), plsScores_Track2(:, 2), plsScores_Track2(:, 3),10,'b','filled','MarkerFaceAlpha',0.1);  hold on;
        %     scatter3(plsScores_Track1(:, 1), plsScores_Track1(:, 2), plsScores_Track1(:, 3),10,'r','filled','MarkerFaceAlpha',0.1);hold on;
        %     scatter3(plsScores_ripples(:, 1), plsScores_ripples(:, 2), plsScores_ripples(:, 3),10,'m','filled','MarkerFaceAlpha',0.1);  hold on;

        % Train SVM on PLS scores
        plsScores_train = X_shuffled(trainIdx,:)*W; % The PLS latent scores (reduced dimensions)
        %     mdl = fitcsvm(plsScores_train, Y(trainIdx), 'KernelFunction', 'linear', 'KernelScale', 'auto', 'BoxConstraint', 1);
        mdl = fitclinear(plsScores_train, Y(trainIdx), 'Learner', 'logistic', 'Regularization', 'ridge');

        plsScores_test = X(testIdx,:)* W;
        predictedTrack = predict(mdl, plsScores_test); % Predict track identity
        %     plot(cumsum(predictedTrack-1));hold on;plot(cumsum(Y(testIdx)-1));

        %     sum(predictedTrack == Y(testIdx))/sum(testIdx)

        temp = position_Z(testIdx);
        correct_position1 = temp(Y(testIdx)==predictedTrack  & Y(testIdx) == 1);
        correct_position2 = temp(Y(testIdx)==predictedTrack  & Y(testIdx) == 2);
        incorrect_position1 = temp(Y(testIdx)~=predictedTrack  & Y(testIdx) == 1);
        incorrect_position2 = temp(Y(testIdx)~=predictedTrack  & Y(testIdx) == 2);

        % Plot 3D PLS scores with labels
        %     figure;
        %     scatter3(plsScores_test(Y(testIdx)==1, 1), plsScores_test(Y(testIdx)==1, 2), plsScores_test(Y(testIdx)==1, 3), 10,'b','filled','MarkerFaceAlpha',0.1); hold on;% Track 1
        %     scatter3(plsScores_test(Y(testIdx)==2, 1), plsScores_test(Y(testIdx)==2, 2), plsScores_test(Y(testIdx)==2, 3), 10,'r','filled','MarkerFaceAlpha',0.1); % Track 2


        PLS_shuffled.correct_position1{nshuffle} = [PLS_shuffled.correct_position1{nshuffle} correct_position1];
        PLS_shuffled.incorrect_position1{nshuffle} = [PLS_shuffled.incorrect_position1{nshuffle} incorrect_position1];
        PLS_shuffled.correct_position2{nshuffle} = [PLS_shuffled.correct_position2{nshuffle} correct_position2];
        PLS_shuffled.incorrect_position2{nshuffle} = [PLS_shuffled.incorrect_position2{nshuffle} incorrect_position2];
    end

    %     PLS.hit_rate{nshuffle} = sum(predictedTrack == Y(testIdx))/sum(testIdx);
%     toc
end

colour_lines = [215,25,28;44,123,182]/256;

x = unique(position_interp1);
for fold = 1:numFolds
    correct_position_distribution1(fold,:) = histcounts(PLS.correct_position1{fold},length(x))./((histcounts(horzcat(PLS.incorrect_position1{fold}),length(x)))+(histcounts(horzcat(PLS.correct_position1{fold}),length(x))));
end
y = mean(correct_position_distribution1);
SE = std(y)./sqrt(numFolds);
PLOT = plot(x,y,'Color',colour_lines(1,:));hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[y+SE fliplr(y-SE)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');


for fold = 1:numFolds
    correct_position_distribution2(fold,:) = histcounts(PLS.correct_position2{fold},length(x))./((histcounts(horzcat(PLS.incorrect_position2{fold}),length(x)))+(histcounts(horzcat(PLS.correct_position2{fold}),length(x))));
end
y = mean(correct_position_distribution2);
SE = std(y)./sqrt(numFolds);
PLOT = plot(x,y,'Color',colour_lines(2,:));hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[y+SE fliplr(y-SE)],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');

correct_position_distribution_shuffled=[];
for nshuffle = 1:500
    correct_position = [PLS_shuffled.correct_position1{nshuffle} PLS_shuffled.correct_position2{nshuffle}];
    incorrect_position = [PLS_shuffled.incorrect_position1{nshuffle} PLS_shuffled.incorrect_position2{nshuffle}];
    correct_position_distribution_shuffled(nshuffle,:) = histcounts(correct_position,length(x))./(histcounts(correct_position,length(x))+histcounts(incorrect_position,length(x)));
end
y = mean(correct_position_distribution_shuffled);
CI_U = prctile(correct_position_distribution_shuffled,97.5,1);
CI_L = prctile(correct_position_distribution_shuffled,2.5,1);
PLOT = plot(x,y,'Color','k');hold on;
ERROR_SHADE(3) = patch([x fliplr(x)],[CI_U fliplr(CI_L)],'k','FaceAlpha','0.3','LineStyle','none');
legend(ERROR_SHADE(1:3),{'Track 1','Track 2','Cell ID shuffle'},'Box','off')
xlabel('Position')
ylabel('Proportion of correct track classification')
title('Logistic Ridge regression of PLS latent components')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


%%%%%  Remove data point when track id classification based on PLS
%%%%%  regression

X = [n(T1_bins,:); n(T2_bins,:)];
Y = [ones(sum(T1_bins),1); 2 * ones(sum(T2_bins),1)];
position_Z = [position_interp1(T1_bins) position_interp1(T2_bins)];

position_bins = unique(position_interp1);
[index,Locb] = ismember(position_Z,position_bins(mean(correct_position_distribution1) > CI_U & mean(correct_position_distribution2) > CI_U));

% Find Bins with 'Good' positions
position_Z = position_Z(index);
X = zscore(X(index,:));
Y = Y(index);

[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(X,Y,3,'CV',10,'Options',statset('UseParallel',true)); % 3 components

[W, ~] = qr(XL, 0); % Orthonormalize XL using QR decomposition
plsScores_Track1 = X(Y == 1,:) * W;
plsScores_Track2 = X(Y == 2,:) * W;
plsScores_shuffled1 = X(Y == 1,random_cell_index) * W;
plsScores_shuffled2 = X(Y == 2,random_cell_index) * W;
plsScores_ripples = n_ripples_residuals * W;
figure;
scatter3(plsScores_Track2(:, 1), plsScores_Track2(:, 2), plsScores_Track2(:, 3),10,'b','filled','MarkerFaceAlpha',0.1);  hold on;
scatter3(plsScores_Track1(:, 1), plsScores_Track1(:, 2), plsScores_Track1(:, 3),10,'r','filled','MarkerFaceAlpha',0.1);hold on;
% scatter3(plsScores_shuffled1(:, 1), plsScores_shuffled1(:, 2), plsScores_shuffled1(:, 3),10,'k','filled','MarkerFaceAlpha',0.1);hold on;

scatter3(plsScores_ripples(:, 1), plsScores_ripples(:, 2), plsScores_ripples(:, 3),10,'m','filled','MarkerFaceAlpha',0.1);  hold on;
s = RandStream('mrg32k3a','Seed',nshuffle); % Set random seed for resampling
random_cell_index = randperm(s,size(n_ripples_residuals,2));
n_ripples_shuffled = n_ripples_residuals(:,random_cell_index);

ripple_bins = [];
ripple_bins(:,1)=ripples_time_edges+0.5;
ripple_bins(:,2)=ripples_time_edges+0.5+time_bin_size;
[~,~,~,~,~,n_ripples_time_shifted] = ActivityTemplates(spike_target,'bins',ripple_bins);
% scatter3(n_ripples_time_shifted* W(:, 1),n_ripples_time_shifted* W(:, 2),n_ripples_time_shifted* W(:, 3),15,'k','filled','MarkerFaceAlpha',0.1);hold on
scatter3(n_ripples_shuffled* W(:, 1),n_ripples_shuffled* W(:, 2),n_ripples_shuffled* W(:, 3),15,'k','filled','MarkerFaceAlpha',0.1);hold on
xlabel('PLS Component 1');
ylabel('PLS Component 2');
legend('Track 1', 'Track 2','Ripples','Ripple (shuffled)');
title('PLS-DA Separation');

for nEvent = 1:max(ripples_id);
    this_event_bins = find(ripples_id==nEvent);
    for nshuffle = 1:500
        s = RandStream('mrg32k3a','Seed',nshuffle); % Set random seed for resampling
        random_cell_index = randperm(s,size(n_ripples_residuals,2));
        n_ripples_shuffled = n_ripples_residuals(this_event_bins,random_cell_index);
        scatter3(n_ripples_shuffled* W(:, 1),n_ripples_shuffled* W(:, 2),n_ripples_shuffled* W(:, 3),15,'k','filled','MarkerFaceAlpha',0.3);hold on
    end
end
scatter3(plsScores_ripples(this_event_bins, 1), plsScores_ripples(this_event_bins, 2), plsScores_ripples(this_event_bins, 3),10,'m','filled','MarkerFaceAlpha',1);  hold on;
xlabel('PLS Component 1');
ylabel('PLS Component 2');
legend('Track 1', 'Track 2');
title('PLS-DA Separation');


%%%%% Representational similarity
% Determine similarity between ripple representation to 
% Track 1 and Track 2 representations

%%%% Mahalanobis distances
plsScores_ripples = n_ripples * W;
dist_to_Track1 = mahal(plsScores_ripples, plsScores_Track1); % MATLAB's mahal function
dist_to_Track2 = mahal(plsScores_ripples, plsScores_Track2);

% Normalize to calculate similarity
similarity_timeseries = dist_to_Track1 ./ (dist_to_Track1 + dist_to_Track2);
% 
% for nEvent = 1:max(ripples_id);
%     this_event_bins = find(ripples_id==nEvent);
%     if length(this_event_bins>3)
%     similarity(this_event_bins);
% 
%     end
% end

similarity_timeseries_shuffled=[];
for nshuffle = 1:1000
    s = RandStream('mrg32k3a','Seed',nshuffle); % Set random seed for resampling
    random_cell_index = randperm(s,size(n_ripples,2));
    n_ripples_shuffled = n_ripples(:,random_cell_index);

    plsScores_ripples = n_ripples_shuffled * W;
    dist_to_Track1 = mahal(plsScores_ripples, plsScores_Track1); % MATLAB's mahal function
    dist_to_Track2 = mahal(plsScores_ripples, plsScores_Track2);
    similarity_timeseries_shuffled(nshuffle,:) = dist_to_Track1 ./ (dist_to_Track1 + dist_to_Track2);
end

similarity_shuffled=[];
similarity = [];
for nEvent = 1:max(ripples_id);
    this_event_bins = find(ripples_id==nEvent);
    similarity_shuffled(:,nEvent) = mean(similarity_timeseries_shuffled(:,this_event_bins),2);
    similarity(nEvent) = mean(similarity_timeseries(this_event_bins));
end

% Visualize similarity
figure;
% similarity_timeseries_shuffled = reshape(similarity_timeseries_shuffled,1,[]);
histogram([similarity_shuffled],100,'Normalization','probability');hold on;
histogram(similarity,100,'Normalization','probability');
xlabel('Ripple Time Bin');
ylabel('Track 2 Similarity');
title('Mahalanobis-Based Similarity');




%%% Nearest-Neighbor-Based Method
dist_to_all_Track1 = pdist2(plsScores_ripples, plsScores_Track1); % Pairwise distances to Track 1
dist_to_all_Track2 = pdist2(plsScores_ripples, plsScores_Track2); % Pairwise distances to Track 2

% Find minimum distances
min_dist_Track1 = min(dist_to_all_Track1, [], 2); % Closest distance to Track 1
min_dist_Track2 = min(dist_to_all_Track2, [], 2); % Closest distance to Track 2

% Normalize to calculate similarity
similarity_timeseries = min_dist_Track2 ./ (min_dist_Track1 + min_dist_Track2);

similarity_timeseries_shuffled=[];
for nshuffle = 1:1000
    s = RandStream('mrg32k3a','Seed',nshuffle); % Set random seed for resampling
    random_cell_index = randperm(s,size(n_ripples,2));
    n_ripples_shuffled = n_ripples(:,random_cell_index);

    plsScores_ripples = n_ripples_shuffled * W;

    %%% Nearest-Neighbor-Based Method
    dist_to_all_Track1 = pdist2(plsScores_ripples, plsScores_Track1); % Pairwise distances to Track 1
    dist_to_all_Track2 = pdist2(plsScores_ripples, plsScores_Track2); % Pairwise distances to Track 2

    % Find minimum distances
    min_dist_Track1 = min(dist_to_all_Track1, [], 2); % Closest distance to Track 1
    min_dist_Track2 = min(dist_to_all_Track2, [], 2); % Closest distance to Track 2

    similarity_timeseries_shuffled(nshuffle,:) = min_dist_Track1 ./ (min_dist_Track1 + min_dist_Track2);
    %     for nEvent = 1:max(ripples_id);
    %         this_event_bins = find(ripples_id==nEvent);
    %         n_ripples_shuffled(this_event_bins,:)
    %     end
end

similarity_shuffled=[];
similarity = [];
for nEvent = 1:max(ripples_id);
    this_event_bins = find(ripples_id==nEvent);
    similarity_shuffled(:,nEvent) = mean(similarity_timeseries_shuffled(:,this_event_bins),2);
    similarity(nEvent) = mean(similarity_timeseries(this_event_bins));

    similarity_T1_sig = find(similarity<prctile(similarity_shuffled,2.5));
    similarity_T2_sig = find(similarity>prctile(similarity_shuffled,97.5));
end
similarity_T1_sig = similarity<similarity_shuffled

similarity;

% Visualize similarity
figure;
% similarity_timeseries_shuffled = reshape(similarity_timeseries_shuffled,1,[]);
% histogram(similarity_shuffled(:,[similarity_T1_sig similarity_T2_sig]),200,'Normalization','probability');hold on;
% histogram(similarity([similarity_T1_sig similarity_T2_sig]),200,'Normalization','probability');
histogram(similarity_shuffled,200,'Normalization','probability');hold on;
histogram(similarity,200,'Normalization','probability');
xlabel('Ripple Time Bin');
ylabel('Track 2 Similarity');
title('Nearest Neighbour Similarity');

