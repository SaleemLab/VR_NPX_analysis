function [PLS KDE_RUN KDE_reactivation] = PLS_KDE_Reactivation_decoding(tvec_template,position,speed,lap_times,track_ID,spikes_template,event_times,spike_target,options)
% PLS_KDE_Reactivation_decoding - 1. Perform PLS regression to obtain PLS latent components that are maximised for covariance between track identity and spiking data during RUN, which can be projected to RUN data or Ripple data
%  2. Use KDE to find probability of track reactivation during ripple based on PLS projection
%
% Inputs:
%   tvec_template   - Time vector template for interpolation.
%   position        - Position data corresponding to the time vector.
%   speed           - Speed data corresponding to the time vector.
%   lap_times       - Lap times for different tracks.
%   track_ID        - Track identifiers for each lap.
%   spikes_template - Spike data template for activity templates.
%   event_times     - Event times for ripple events.
%   spike_target    - Spike data for target neurons.
%   options         - Struct containing additional options (e.g., SUBJECT, SESSION).
%
% Outputs:
%   PLS             - Struct containing PLS regression results, including weights, indices, and hit rates.
%   KDE_RUN         - Struct containing KDE analysis results during running periods, including bandwidth and AUC.
%   KDE_reactivation- Struct containing KDE reactivation analysis results, including ripple bias and probabilities.
%
% Description:
%   This function performs Partial Least Squares (PLS) regression and Kernel Density Estimation (KDE) analysis on neural data.
%   It first interpolates the input data to create time bins, then calculates activity templates and performs PLS regression.
%   The function also performs cross-validation and logistic regression to classify track identity based on PLS components.
%   Additionally, it calculates KDE for ripple events and evaluates reactivation strength based on PLS projections.
%
% Example:
%   [PLS, KDE_RUN, KDE_reactivation] = PLS_KDE_Reactivation_decoding(tvec_template, position, speed, lap_times, track_ID, spikes_template, event_times, spike_target, options);
%

    % Initialize output variables
PLS = struct();
KDE_RUN = struct();
KDE_reactivation = struct();

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

time_bin_size=0.02;
time_bin_size_moving = 0.01;
% timevec_edge= interp1(tvec_target,tvec_target,tvec_target(1):time_bin_size:tvec_target(end));
% speed_RUN= interp1(tvec_template,speed,tvec_template(1):time_bin_size:tvec_template(end));
% timevec_edge_RUN(speed_RUN<1)=nan;
% tvec_template

ripples_time_edges=[];
ripples_id=[];
for i = 1: length(event_times)
    event_duration = event_times(i,2) - event_times(i,1)-0.02;
%     event_duration = event_times(i,2) - event_times(i,1);
    if event_duration <0.1
        event_duration = 0.1;
    end
    num_bins = floor(event_duration / time_bin_size_moving);
    timebins_edges = linspace(event_times(i,1), event_times(i,1) + num_bins *  time_bin_size_moving, num_bins+1);
    ripples_time_edges = [ripples_time_edges timebins_edges];
    ripples_id = [ripples_id i*ones(1,length(timebins_edges))];
end

ripple_bins = [];
ripple_bins(:,1)=ripples_time_edges;
ripple_bins(:,2)=ripples_time_edges+time_bin_size;
[~,~,~,~,~,n_ripples] = ActivityTemplates(spike_target,'bins',ripple_bins);


% Calculate ripple residual
% Calculate mean activity across all neurons for each ripple
mean_ripple_activity = mean(n_ripples, 2); % Average firing rate per ripple
% Regress out the mean ripple activity
coeff = (mean_ripple_activity' * n_ripples) / (mean_ripple_activity' * mean_ripple_activity);
ripple_global_pattern = mean_ripple_activity * coeff; % Global dynamics component
n_ripples_residuals = zscore(n_ripples - ripple_global_pattern,0,1); % Residuals

% n_ripples_residuals = n_ripples;
% n_ripples_residuals = zscore(n_ripples - mean_ripple_activity); % Residuals
% [status,interval,index] = InIntervals(values,event_times)

%%%%% PCA or ICA
% [templates,correlations,projection1,weights,variance] = ActivityTemplates(spikes_template,'bins',bins(T1_bins,:),'controlbins',bins(T2_bins,:),'mode','pca');
%           [templates,correlations,projection1,weights,variance] = ActivityTemplates(spikes_template,'bins',bins,'mode','pca');
%           scatter3(n(T1_bins,:) * weights(:, 1),n(T1_bins,:) * weights(:, 2),n(T1_bins,:) * weights(:, 3),20,position_interp1(T1_bins)'/max(position_interp1),'filled');hold on
%           scatter3(n(T2_bins,:) * weights(:, 1),n(T2_bins,:) * weights(:, 2),n(T2_bins,:) * weights(:, 3),20,position_interp1(T2_bins)'/max(position_interp1),'filled');hold on
%           colorbar
%           colormap(gray)

% scatter3(n(T1_bins,:) * weights(:, 1),n(T1_bins,:) * weights(:, 2),n(T1_bins,:) * weights(:, 3),15,'r','filled','MarkerFaceAlpha',0.1);hold on
% scatter3(n(T2_bins,:) * weights(:, 1),n(T2_bins,:) * weights(:, 2),n(T2_bins,:) * weights(:, 3),15,'b','filled','MarkerFaceAlpha',0.1);hold on
% scatter3(n_ripples* weights(:, 1),n_ripples* weights(:, 2),n_ripples* weights(:, 3),15,'magenta','filled','MarkerFaceAlpha',0.1);hold on

% for nshuffle = 1:1000
%     s = RandStream('mrg32k3a','Seed',nshuffle); % Set random seed for resampling
%     random_cell_index = randperm(s,size(n_ripples,2));
%     n_ripples_shuffled = n_ripples(:,random_cell_index);
% %     scatter3(n_ripples_shuffled* weights(:, 1),n_ripples_shuffled* weights(:, 2),n_ripples_shuffled* weights(:, 3),15,'k','filled','MarkerFaceAlpha',0.1);hold on
% 
% end
disp('Caclulate PLS during RUN')
% Combine data and labels
X = [n(T1_bins,:); n(T2_bins,:)];
Y = [ones(sum(T1_bins),1); 2 * ones(sum(T2_bins),1)];
position_Z = [position_interp1(T1_bins) position_interp1(T2_bins)];

[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(X,Y,3); % 3 components
[W, ~] = qr(XL, 0); % Orthonormalize XL using QR decomposition

% figure
% plsScores_Track1 = X(Y == 1,:) * W;
% plsScores_Track2 = X(Y == 2,:) * W;
% scatter3(plsScores_Track2(:, 1), plsScores_Track2(:, 2), plsScores_Track2(:, 3),10,'b','filled','MarkerFaceAlpha',0.1);  hold on;
% scatter3(plsScores_Track1(:, 1), plsScores_Track1(:, 2), plsScores_Track1(:, 3),10,'r','filled','MarkerFaceAlpha',0.1);hold on;

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

region_name = {'Left','Right'};

colour_lines = [215,25,28;44,123,182]/256;
title_text = sprintf('%s %s Logistic Ridge regression of PLS components %s',options.SUBJECT,options.SESSION);
nfigure = 1;
fig(nfigure)=figure;
fig(nfigure).Name=title_text;
fig(nfigure).Position = [680 482 630 500]
% fig(nfigure).Name=sprintf('Logistic Ridge regression of PLS latent components %s',region_name{options.probe_hemisphere});

x = unique(position_interp1);
for fold = 1:numFolds
    correct_position_distribution1(fold,:) = histcounts(PLS.correct_position1{fold},length(x))./((histcounts(horzcat(PLS.incorrect_position1{fold}),length(x)))+(histcounts(horzcat(PLS.correct_position1{fold}),length(x))));
end
y = mean(correct_position_distribution1);
SE = std(correct_position_distribution1)./sqrt(numFolds);
PLOT = plot(x,y,'Color',colour_lines(1,:));hold on;
ERROR_SHADE(1) = patch([x fliplr(x)],[y+SE fliplr(y-SE)],colour_lines(1,:),'FaceAlpha','0.3','LineStyle','none');


for fold = 1:numFolds
    correct_position_distribution2(fold,:) = histcounts(PLS.correct_position2{fold},length(x))./((histcounts(horzcat(PLS.incorrect_position2{fold}),length(x)))+(histcounts(horzcat(PLS.correct_position2{fold}),length(x))));
end
y = mean(correct_position_distribution2);
SE = std(correct_position_distribution2)./sqrt(numFolds);
PLOT = plot(x,y,'Color',colour_lines(2,:));hold on;
ERROR_SHADE(2) = patch([x fliplr(x)],[y+SE fliplr(y-SE)],colour_lines(2,:),'FaceAlpha','0.3','LineStyle','none');

correct_position_distribution_shuffled=[];
for nshuffle = 1:1000
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
xlabel('Position (cm)')
ylabel('Proportion of correct track classification')
title(title_text)
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

PLS.good_idx = index;
PLS.good_position = unique(position_Z);


%%%%% Project PLS weights to ripple data
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(X,Y,3,'CV',10,'Options',statset('UseParallel',true)); % 3 components
[W, ~] = qr(XL, 0); % Orthonormalize XL using QR decomposition

%%%%%%%%%%%%%% Reactivation strength
% % calculate projection matrix 
% P_template = zeros(size(W,1),size(W,1),size(W,2));
% for i = 1:size(W,2)
%     P_template(:,:,i) = W(:,i)*W(:,i)';
%     P_template(:,:,i) = P_template(:,:,i) - diag(diag(P_template(:,:,i))); % remove the diagonal
% end
% % figure
% % for i = 1:size(W,2)
% %     subplot(2,2,i)
% %     imagesc(squeeze(P_template(:,:,i)))
% %     colorbar
% % end
% % 
% % Compute reactivation strengths
% nTemplates = size(P_template,3);
% strength = zeros(size(n_ripples,1),nTemplates);
% for i = 1:nTemplates
%     template = P_template(:,:,i);
%     % Set the diagonal to zero to not count coactivation of i and j when i=j
%     template = template - diag(diag(template));
%     strength(:,i) = nansum(n_ripples*(template).*n_ripples,2);
% end
% 
% % Compute reactivation strengths
% nTemplates = size(P_template,3);
% strength_shuffled = zeros(size(n_ripples,1),nTemplates,1000);
% for nshuffle = 1:1000
%     for i = 1:nTemplates
%         template = P_template(:,:,i);
%         % Set the diagonal to zero to not count coactivation of i and j when i=j
%         template = template - diag(diag(template));
%         s = RandStream('mrg32k3a','Seed',nshuffle); % Set random seed for resampling
%         random_cell_index = randperm(s,size(n_ripples_residuals,2));
% 
%         n_ripples_shuffled = n_ripples(:,random_cell_index);
%         strength_shuffled(:,i,nshuffle) = nansum(n_ripples_shuffled*(template).*n_ripples_shuffled,2);
%     end
% end
% 
% 
% for i = 1:nTemplates
% %     prctile(strength(:,i),0.5)
% %     prctile(strength(:,i),99.5)
%     subplot(2,2,i)
%     histogram(strength(:,i), prctile(strength(:,i),2.5):0.01:prctile(strength(:,i),97.5),'Normalization','probability','EdgeAlpha',0);hold on
%     histogram(strength_shuffled(:,i,:),prctile(strength(:,i),2.5):0.01:prctile(strength(:,i),97.5),'Normalization','probability','EdgeAlpha',0)
% end



%%%%% Visualise projection

region_name = {'Left','Right'};

colour_lines = [215,25,28;44,123,182]/256;
title_text = sprintf('%s %s PLS latent components visualisation %s',options.SUBJECT,options.SESSION);
nfigure = 2;
fig(nfigure)=figure;
fig(nfigure).Name=title_text;
fig(nfigure).Position = [300 380 1300 600]

% s = RandStream('mrg32k3a','Seed',2); % Set random seed for resampling
% random_cell_index = randperm(s,size(n_ripples,2));
plsScores_Track1 = X(Y == 1,:) * W;
plsScores_Track2 = X(Y == 2,:) * W;
% plsScores_shuffled1 = X(Y == 1,random_cell_index) * W;
% plsScores_shuffled2 = X(Y == 2,random_cell_index) * W;
% plsScores_ripples = n_ripples_residuals * W;
% figure;
subplot(2,2,1)
scatter3(plsScores_Track2(:, 1), plsScores_Track2(:, 2), plsScores_Track2(:, 3),10,'b','filled','MarkerFaceAlpha',0.1);  hold on;
scatter3(plsScores_Track1(:, 1), plsScores_Track1(:, 2), plsScores_Track1(:, 3),10,'r','filled','MarkerFaceAlpha',0.1);hold on;
xlabel('PLS Component 1');
ylabel('PLS Component 2');
zlabel('PLS Component 3');
legend('Track 1', 'Track 2','Box','off');
title('T1-T2 PLS components');
% legend(ERROR_SHADE(1:3),{'Track 1','Track 2','Cell ID shuffle'},'Box','off')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)


subplot(2,2,2)
scatter3(plsScores_Track2(:, 1), plsScores_Track2(:, 2), plsScores_Track2(:, 3),10,'b','filled','MarkerFaceAlpha',0.1);  hold on;
scatter3(plsScores_Track1(:, 1), plsScores_Track1(:, 2), plsScores_Track1(:, 3),10,'r','filled','MarkerFaceAlpha',0.1);hold on;
plsScores_ripples = n_ripples * W;
scatter3(plsScores_ripples(:, 1), plsScores_ripples(:, 2), plsScores_ripples(:, 3),10,'m','filled','MarkerFaceAlpha',0.1);  hold on;
legend('Track 1', 'Track 2','Ripples','Box','off')
title('ripple PLS projection');
% legend(ERROR_SHADE(1:3),{'Track 1','Track 2','Cell ID shuffle'},'Box','off')
title(title_text)
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
xlabel('PLS Component 1');
ylabel('PLS Component 2');
zlabel('PLS Component 3');


subplot(2,2,4)
scatter3(plsScores_Track2(:, 1), plsScores_Track2(:, 2), plsScores_Track2(:, 3),10,'b','filled','MarkerFaceAlpha',0.1);  hold on;
scatter3(plsScores_Track1(:, 1), plsScores_Track1(:, 2), plsScores_Track1(:, 3),10,'r','filled','MarkerFaceAlpha',0.1);hold on;
plsScores_ripples = n_ripples_residuals * W;
scatter3(plsScores_ripples(:, 1), plsScores_ripples(:, 2), plsScores_ripples(:, 3),10,'m','filled','MarkerFaceAlpha',0.1);  hold on;
legend('Track 1', 'Track 2','Ripples residuals','Box','off')
title('ripple PLS projection');
% legend(ERROR_SHADE(1:3),{'Track 1','Track 2','Cell ID shuffle'},'Box','off')
title(title_text)
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
xlabel('PLS Component 1');
ylabel('PLS Component 2');
zlabel('PLS Component 3');

sgtitle(title_text)

% scatter3(plsScores_shuffled2(:, 1), plsScores_Track2(:, 2), plsScores_Track2(:, 3),10,'b','filled','MarkerFaceAlpha',0.01);  hold on;
% scatter3(plsScores_shuffled1(:, 1), plsScores_Track1(:, 2), plsScores_Track1(:, 3),10,'r','filled','MarkerFaceAlpha',0.01);hold on;
% 
% % scatter3(plsScores_shuffled1(:, 1), plsScores_shuffled1(:, 2), plsScores_shuffled1(:, 3),10,'k','filled','MarkerFaceAlpha',0.1);hold on;

% scatter3(plsScores_ripples(:, 1), plsScores_ripples(:, 2), plsScores_ripples(:, 3),10,'m','filled','MarkerFaceAlpha',0.1);  hold on;
% 
% 
% s = RandStream('mrg32k3a','Seed',nshuffle); % Set random seed for resampling
% random_cell_index = randperm(s,size(n_ripples_residuals,2));
% [XL_shuffled] = plsregress(X(:,random_cell_index),Y,3,'CV',10,'Options',statset('UseParallel',true)); % 3 components
% [W_shuffled, ~] = qr(XL_shuffled, 0); % Orthonormalize XL using QR decomposition
% n_ripples_residuals_shuffled = n_ripples(:,random_cell_index);
% % ripple_bins = [];
% % ripple_bins(:,1)=ripples_time_edges+0.5;
% % ripple_bins(:,2)=ripples_time_edges+0.5+time_bin_size;
% % [~,~,~,~,~,n_ripples_time_shifted] = ActivityTemplates(spike_target,'bins',ripple_bins);
% % scatter3(n_ripples_time_shifted* W(:, 1),n_ripples_time_shifted* W(:, 2),n_ripples_time_shifted* W(:, 3),15,'k','filled','MarkerFaceAlpha',0.1);hold on
% scatter3(n_ripples_residuals_shuffled* W(:, 1),n_ripples_residuals_shuffled* W(:, 2),n_ripples_residuals_shuffled* W(:, 3),15,'k','filled','MarkerFaceAlpha',0.05);hold on
% % scatter3(n_ripples_residuals* W_shuffled(:, 1),n_ripples_residuals* W_shuffled(:, 2),n_ripples_residuals* W_shuffled(:, 3),15,'k','filled','MarkerFaceAlpha',0.1);hold on
% xlabel('PLS Component 1');
% ylabel('PLS Component 2');
% legend('Track 1', 'Track 2','Ripples','Ripple (shuffled)');
% title('PLS-DA Separation');

disp('Find best bandwidth for KDE of PLS projection')
%%%%% Create 10-fold partition for finding best bandwidth for KDE PLS
rng(200)
numFolds = 10;
cv = cvpartition(Y, 'KFold', numFolds);

% KDE for Track 1 and Track 2 with optimized bandwidth
Bandwith_list = 0.1:0.1:2;
track_separation=[];
KDE_RUN = [];
% title_text = sprintf('%s %s PLS components KDE RUN ',options.SUBJECT,options.SESSION);

for i = 1:20
    pdf_T1_CV=[];
    pdf_T2_CV=[];
    Y_CV = [];
    Z_CV = [];
    for fold = 1:numFolds
        s = RandStream('mrg32k3a','Seed',fold); % Set random seed for resampling
        %         random_cell_index = randperm(s,size(n_ripples_residuals,2));
        %         random_cell_index = randperm(s,size(n_ripples_residuals,2));
        %         n_ripples_shuffled = n_ripples(:,random_cell_index);

        % Get training and test indices for this fold
        trainIdx = training(cv, fold); % Logical vector for training samples
        testIdx = test(cv, fold);     % Logical vector for testing samples

        [pdf_T1] = mvksdensity( X(Y == 1 & trainIdx==1,:) * W,X(testIdx,:) * W, 'Bandwidth',Bandwith_list(i)); % Cross-validated bandwidth
        [pdf_T2] = mvksdensity( X(Y == 2 & trainIdx==1,:) * W,X(testIdx,:) * W, 'Bandwidth', Bandwith_list(i));

        pdf_T1_CV = [pdf_T1_CV;pdf_T1];
        pdf_T2_CV = [pdf_T2_CV;pdf_T2];
        Y_CV = [Y_CV; Y(Y == 1 & testIdx==1) ;Y(Y == 2 & testIdx==1)];% Track ID
        Z_CV = [Z_CV; position_Z(Y == 1 & testIdx==1)' ;position_Z(Y == 2 & testIdx==1)']; % Position
        %     [pdf_T1] = mvksdensity( X(Y == 1,:) * W,n_ripples * W, 'Bandwidth',Bindiwith_list(i)); % Cross-validated bandwidth
        %     [pdf_T2] = mvksdensity( X(Y == 2,:) * W,n_ripples * W, 'Bandwidth', Bindiwith_list(i));
    end
    %     pdf_T2_CV = movmedian(pdf_T2_CV,3);
    %     pdf_T1_CV = movmedian(pdf_T1_CV,3);

    % Compute separability metric (e.g., Bhattacharyya distance)
%     separability(i) = -log(sum(sqrt(pdf_T1_CV .* pdf_T2_CV))); % Bhattacharyya distance
%     [emd_value] = emd(pdf_T1_CV, pdf_T2_CV, distance_metric);

    RUN_bias = pdf_T1_CV./(pdf_T1_CV+pdf_T2_CV);
    for nboot = 1:1000
        s = RandStream('mrg32k3a','Seed',nboot); % Set random seed for resampling

        index = datasample(s,1:length(RUN_bias),length(RUN_bias));
        data_resampled(nboot,:) = RUN_bias(index);
        track_label_resampled(nboot,:) = Y_CV(index)-1;
        [x,y,T,A] = perfcurve(Y_CV(index),RUN_bias(index),1,'XVals',0:0.05:1,'NBoot',1,'BootType','per');
        %             [X,Y,T,A] = perfcurve(track_label(index),z_log_odds(index),1,'XVals',0:0.05:1);

        FPR = x;
        TPR(nboot,:) = y(:,1);
        AUC(nboot) = A(1);
    end

    KDE_RUN.Bandwidth(i) = Bandwith_list(i);
    KDE_RUN.AUC(i,:) = AUC;
    KDE_RUN.FPR = FPR;
    KDE_RUN.TPR(i,:,:) = TPR;


%     fig(10) = figure(10);
%     fig(10).Position = [300 70 1330 895];
%     fig(10).Name = [title_text,'T1-T2 bias distribution']
% 
% 
%     subplot(4,5,i)
%     histogram(pdf_T2_CV./(pdf_T1_CV+pdf_T2_CV),0:0.005:1,'EdgeAlpha',0)
%     title(sprintf('AUC %.4f Bandwith %.1f',median(AUC),Bandwith_list(i)))
%     set(gca,"TickDir","out",'box', 'off','Color','none')
%     fontsize(gcf,12,"points")
% 
%     fig(11) = figure(11);
%     fig(11).Position = [300 70 1330 895];
%     fig(11).Name = [title_text,'T1-T2 bias']
% 
%     subplot(4,5,i)
%     plot(pdf_T2_CV./(pdf_T1_CV+pdf_T2_CV));hold on;plot(Y_CV-1)
%     title(sprintf('AUC %.4f Bandwith %.1f',median(AUC),Bandwith_list(i)))
%     set(gca,"TickDir","out",'box', 'off','Color','none')
%     fontsize(gcf,12,"points")
% 
%     fig(12) = figure(12);
%     fig(12).Position = [300 70 1330 895];
%     fig(12).Name = [title_text,'T1-T2 bias AUC']
% 
%     subplot(4,5,i)
%     hold on
%     x = FPR';
%     CI_shuffle = prctile(TPR,[2.5 97.5]);
% %     plot(x, CI_shuffle(2,:), 'r--', 'LineWidth', 1);
% %     plot(x, CI_shuffle(1,:), 'r--', 'LineWidth', 1);
%     x2 = [x, fliplr(x)];
%     inBetween = [CI_shuffle(1,:), fliplr(CI_shuffle(2,:))];
%     h(2) = fill(x2, inBetween, 'r','FaceAlpha',0.2,'EdgeColor','r');
% 
%     h(1) = plot([0 1],[0 1],'k--')
%     title(sprintf('AUC %.4f Bandwith %.1f',median(AUC),Bandwith_list(i)))
%     set(gca,"TickDir","out",'box', 'off','Color','none')
%     %     legend([h(2) h(1)],{'Real','chance'})
%     xlabel('False positive rate')
%     ylabel('True positive rate')
%     %                 sgtitle('lap z log odds ROC two track discrimination in V1 for left probe')
%     fontsize(gcf,12,"points")
end

%%% Best bandwith for KDE based on RUN Track 1 and Track 2 separation
mean_AUC = median(KDE_RUN.AUC,2);
AUC_UCI = prctile(KDE_RUN.AUC,97.5,2);
AUC_LCI = prctile(KDE_RUN.AUC,2.5,2);

sig_difference=[];
CI_narrowness = [];
for i = 1:size(KDE_RUN.AUC,1)
    for j = 1:size(KDE_RUN.AUC,1)

        sig_difference(i,j) = AUC_LCI(i) > AUC_UCI(j) ;
    end
    CI_narrowness(i) = AUC_UCI(i) - AUC_LCI(i) ;
end
[~,index] = max(mean_AUC);% Bandwith that leads to maximum separation between tracks
Bandwidth = Bandwith_list(index);


% n_ripples_residuals = n_ripples;
%%% Reactivation analysis based on KDE of PLS projection to ripple data
[pdf_T1] = mvksdensity( X(Y == 1,:) * W,n_ripples_residuals * W, 'Bandwidth',Bandwidth); % Cross-validated bandwidth
[pdf_T2] = mvksdensity( X(Y == 2,:) * W,n_ripples_residuals * W, 'Bandwidth', Bandwidth);

pdf_T1_shuffled = zeros(1000,size(n_ripples_residuals,1));
pdf_T2_shuffled = zeros(1000,size(n_ripples_residuals,1));
disp('Caclulate Reactivation based on KDE of PLS projection to ripple data')
% tic
parfor nshuffle = 1:1000
    tic
    s = RandStream('mrg32k3a','Seed',nshuffle); % Set random seed for resampling
    random_cell_index = randperm(s,size(n_ripples_residuals,2));
    n_ripples_shuffled = n_ripples_residuals(:,random_cell_index);
    pdf_T1_shuffled(nshuffle,:) = mvksdensity( X(Y == 1,:) * W,n_ripples_shuffled * W, 'Bandwidth',Bandwidth); % Cross-validated bandwidth
    pdf_T2_shuffled(nshuffle,:)= mvksdensity( X(Y == 2,:) * W,n_ripples_shuffled * W, 'Bandwidth', Bandwidth);
    toc
end
% toc

KDE_bias_shuffled_event=[];
KDE_bias_event = [];
KDE_bias = pdf_T1./(pdf_T1+pdf_T2);
KDE_bias_shuffled = pdf_T1_shuffled./(pdf_T1_shuffled+pdf_T2_shuffled);

count=0;
for nEvent = 1:max(ripples_id)
    count = count + 1;
    this_event_bins = find(ripples_id==nEvent);
%     KDE_bias_shuffled_event(:,nEvent) = mean(KDE_bias_shuffled(:,this_event_bins),2);
    [KDE_bias_T1_event(nEvent) tidx] = max(KDE_bias(this_event_bins));
    KDE_bias_shuffled_T1_event(:,nEvent) = KDE_bias_shuffled(:,this_event_bins(tidx));
    KDE_T1_event_pvalue(nEvent) = sum(KDE_bias_T1_event(nEvent)>KDE_bias_shuffled_T1_event(:,nEvent))/size(KDE_bias_shuffled_T1_event,1);
    KDE_T1_peak_tidx(nEvent) = this_event_bins(tidx);

    this_event_bins = find(ripples_id==nEvent);
    %     KDE_bias_shuffled_event(:,nEvent) = mean(KDE_bias_shuffled(:,this_event_bins),2);
    [KDE_bias_T2_event(nEvent) tidx] = min(KDE_bias(this_event_bins));
    KDE_bias_shuffled_T2_event(:,nEvent) = KDE_bias_shuffled(:,this_event_bins(tidx));
    KDE_T2_event_pvalue(nEvent) = sum(KDE_bias_T2_event(nEvent)<KDE_bias_shuffled_T2_event(:,nEvent))/size(KDE_bias_shuffled_T2_event,1);
    KDE_T2_peak_tidx(nEvent) = this_event_bins(tidx);

%     KDE_event_bias(nEvent) = sum(KDE_bias(this_event_bins));
%     KDE_event_bias_shuffled(:,nEvent) = sum(KDE_bias_shuffled(:,this_event_bins),2);
%      sum(KDE_event_bias(nEvent)<KDE_bias_shuffled_T2_event(:,nEvent))/size(KDE_bias_shuffled_T2_event,1);
end
KDE_reactivation = [];
KDE_reactivation.ripple_id = ripples_id;
KDE_reactivation.ripple_spike_counts = n_ripples;
KDE_reactivation.ripple_bins = ripple_bins;
KDE_reactivation.ripple_bias = KDE_bias;
KDE_reactivation.ripple_T1_probability = pdf_T1;
KDE_reactivation.ripple_T2_probability = pdf_T2;

KDE_reactivation.ripple_proj = W;
KDE_reactivation.T1_peak_tidx = KDE_T1_peak_tidx;
KDE_reactivation.T2_peak_tidx = KDE_T2_peak_tidx;
KDE_reactivation.T1_peak_bias = KDE_bias(KDE_T1_peak_tidx);
KDE_reactivation.T2_peak_bias = KDE_bias(KDE_T2_peak_tidx);
KDE_reactivation.T1_pvalue = KDE_T1_event_pvalue;
KDE_reactivation.T2_pvalue = KDE_T2_event_pvalue;


KDE_reactivation.ripple_T1_probability_shuffled = pdf_T1_shuffled;
KDE_reactivation.ripple_T2_probability_shuffled = pdf_T1_shuffled;



%%%%% Visualisation of Peak bias distribution
nfigure = nfigure+ 1;
fig(nfigure) = figure;
fig(nfigure).Name=sprintf('KDE peak bias distribution %s ripples',region_name{options.probe_hemisphere});

subplot(2,2,1)
histogram(KDE_bias_T1_event,100,'FaceAlpha',0.5,'Normalization','probability');hold on
histogram(KDE_bias_shuffled_T1_event,100,'FaceAlpha',0.5,'Normalization','probability')
set(gca,"TickDir","out",'box', 'off','Color','none')

subplot(2,2,2)
histogram(KDE_bias_T2_event,100,'FaceAlpha',0.5,'Normalization','probability');hold on
histogram(KDE_bias_shuffled_T2_event,100,'FaceAlpha',0.5,'Normalization','probability')
set(gca,"TickDir","out",'box', 'off','Color','none')
fontsize(gcf,12,"points")

% sum(KDE_T1_event_pvalue >= 0.95 & KDE_T2_event_pvalue >= 0.95)
% sum(KDE_T1_event_pvalue >= 0.975)
% sum(KDE_T2_event_pvalue >= 0.975)


%%%%% Visualisation of significant T1 and T2 ripple events
nfigure = nfigure+ 1;
fig(nfigure) = figure;
fig(nfigure).Position = [279,55,1800,933];
fig(nfigure).Name=sprintf('T1 events KDE bias %s ripples',region_name{options.probe_hemisphere});
count = 0
for nEvent = find(KDE_T1_event_pvalue >=0.95)
    count = count + 1;
    if count < 201
        this_event_bins = find(ripples_id==nEvent);

        x = 0:time_bin_size_moving:time_bin_size_moving*(length(this_event_bins)-1);
        subplot(10,20,count)
        plot(x,(KDE_bias(this_event_bins)),'r')
        hold on
        y = mean(KDE_bias_shuffled(:,this_event_bins));
        CI_U = prctile(KDE_bias_shuffled(:,this_event_bins),97.5,1);
        CI_L = prctile(KDE_bias_shuffled(:,this_event_bins),2.5,1);
        PLOT = plot(x,y,'Color','k');hold on;
        ERROR_SHADE = patch([x fliplr(x)],[CI_U fliplr(CI_L)],'k','FaceAlpha','0.3','LineStyle','none');
        ylim([0 1])
        set(gca,"TickDir","out",'box', 'off','Color','none')
    end
end
fontsize(gcf,10,"points")
% Add a common x-label and y-label
% Use 'Position' to adjust placement
han = axes(fig(nfigure), 'Visible', 'off'); % Create invisible axes
han.XLabel.Visible = 'on'; % Turn on visibility for x-label
han.YLabel.Visible = 'on'; % Turn on visibility for y-label
xlabel(han, 'Time (s)');
ylabel(han, 'T1/T2 bias');


fig(nfigure) = figure
fig(nfigure).Position = [279,55,1800,933];
fig(nfigure).Name=sprintf('T2 events KDE bias %s ripples',region_name{options.probe_hemisphere});
count = 0
for nEvent = find(KDE_T2_event_pvalue >=0.95)
    count = count + 1;
    if count < 201
        this_event_bins = find(ripples_id==nEvent);

        x = 0:time_bin_size_moving:time_bin_size_moving*(length(this_event_bins)-1);
        subplot(10,20,count)
        plot(x,(KDE_bias(this_event_bins)),'r')
        hold on
        y = mean(KDE_bias_shuffled(:,this_event_bins));
        CI_U = prctile(KDE_bias_shuffled(:,this_event_bins),97.5,1);
        CI_L = prctile(KDE_bias_shuffled(:,this_event_bins),2.5,1);
        PLOT = plot(x,y,'Color','k');hold on;
        ERROR_SHADE = patch([x fliplr(x)],[CI_U fliplr(CI_L)],'k','FaceAlpha','0.3','LineStyle','none');
        ylim([0 1])
        set(gca,"TickDir","out",'box', 'off','Color','none')
    end
end
fontsize(gcf,10,"points")
% Add a common x-label and y-label
% Use 'Position' to adjust placement
han = axes(fig(nfigure), 'Visible', 'off'); % Create invisible axes
han.XLabel.Visible = 'on'; % Turn on visibility for x-label
han.YLabel.Visible = 'on'; % Turn on visibility for y-label
xlabel(han, 'Time (s)');
ylabel(han, 'T1/T2 bias');


disp('PLS KDE decoding finished')
% %%%%%% Visualising single event
% figure;
% scatter3(plsScores_Track2(:, 1), plsScores_Track2(:, 2), plsScores_Track2(:, 3),10,'b','filled','MarkerFaceAlpha',0.1);  hold on;
% scatter3(plsScores_Track1(:, 1), plsScores_Track1(:, 2), plsScores_Track1(:, 3),10,'r','filled','MarkerFaceAlpha',0.1);hold on;
% 
% for nEvent = 1:max(ripples_id);
%     this_event_bins = find(ripples_id==nEvent);
%     for nshuffle = 1:1000
%         s = RandStream('mrg32k3a','Seed',nshuffle); % Set random seed for resampling
%         random_cell_index = randperm(s,size(n_ripples_residuals,2));
%         n_ripples_shuffled = n_ripples_residuals(this_event_bins,random_cell_index);
%         scatter3(n_ripples_shuffled* W(:, 1),n_ripples_shuffled* W(:, 2),n_ripples_shuffled* W(:, 3),15,'k','filled','MarkerFaceAlpha',0.05);hold on
%     end
% end
% for tidx = 1:length(this_event_bins)
% scatter3(plsScores_ripples(this_event_bins(tidx), 1), plsScores_ripples(this_event_bins(tidx), 2), plsScores_ripples(this_event_bins(tidx), 3),10,'','filled','MarkerFaceAlpha',1);  hold on;
% end
% xlabel('PLS Component 1');
% ylabel('PLS Component 2');
% legend('Track 1', 'Track 2');
% title('PLS-DA Separation');
% 

% 
% %%%%% Representational similarity
% % Determine similarity between ripple representation to 
% % Track 1 and Track 2 representations
% 
% %%%% Mahalanobis distances
% [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(X,Y,3,'CV',10,'Options',statset('UseParallel',true)); % 3 components
% [W, ~] = qr(XL, 0); % Orthonormalize XL using QR decomposition
% plsScores_ripples = n_ripples_residuals * W;
% dist_to_Track1 = mahal(plsScores_ripples, plsScores_Track1); % MATLAB's mahal function
% dist_to_Track2 = mahal(plsScores_ripples, plsScores_Track2);
% 
% % Normalize to calculate similarity
% similarity_timeseries = dist_to_Track1 ./ (dist_to_Track1 + dist_to_Track2);
% % 
% % for nEvent = 1:max(ripples_id);
% %     this_event_bins = find(ripples_id==nEvent);
% %     if length(this_event_bins>3)
% %     similarity(this_event_bins);
% % 
% %     end
% % end
% dist_to_Track1_shuffled=[];
% dist_to_Track2_shuffled=[];
% similarity_timeseries_shuffled=[];
% for nshuffle = 1:1000
%     s = RandStream('mrg32k3a','Seed',nshuffle); % Set random seed for resampling
%     random_cell_index = randperm(s,size(n_ripples,2));
% % 
% %     % Calculate mean activity across all neurons for each ripple
% %     mean_ripple_activity = mean(n_ripples, 2); % Average firing rate per ripple
% %     % Regress out the mean ripple activity
% %     coeff = (mean_ripple_activity' * n_ripples(:,random_cell_index)) / (mean_ripple_activity' * mean_ripple_activity);
% %     ripple_global_pattern = mean_ripple_activity * coeff; % Global dynamics component
% %     n_ripples_residuals_shuffled = zscore(n_ripples(:,random_cell_index) - ripple_global_pattern,0,1); % Residuals
% % plsScores_ripples = n_ripples_residuals_shuffled * W;
%     [XL_shuffled] = plsregress(X(:,random_cell_index),Y,3,'CV',10,'Options',statset('UseParallel',true)); % 3 components
%     [W_shuffled, ~] = qr(XL_shuffled, 0); % Orthonormalize XL using QR decomposition
% 
% %     n_ripples_shuffled = n_ripples_residuals(:,random_cell_index);
% 
%     plsScores_ripples = n_ripples_residuals * W_shuffled;
%     dist_to_Track1_shuffled(nshuffle,:) = mahal(plsScores_ripples, plsScores_Track1); % MATLAB's mahal function
%     dist_to_Track2_shuffled(nshuffle,:) = mahal(plsScores_ripples, plsScores_Track2);
%     similarity_timeseries_shuffled(nshuffle,:) = dist_to_Track1_shuffled(nshuffle,:) ./ (dist_to_Track1_shuffled(nshuffle,:) + dist_to_Track2_shuffled(nshuffle,:));
% end
% 
% similarity_shuffled=[];
% similarity = [];
% max(ripples_id);
% figure
% for nEvent = 1:100
%     this_event_bins = find(ripples_id==nEvent);
%     %     similarity_shuffled(:,nEvent) = mean(similarity_timeseries_shuffled(:,this_event_bins),2);
%     %     similarity(nEvent) = mean(similarity_timeseries(this_event_bins));
% 
%     [similarity(nEvent) tidx] = min(similarity_timeseries(this_event_bins));
%     similarity_shuffled(:,nEvent) = similarity_timeseries_shuffled(:,this_event_bins(tidx));
%     x = 0:time_bin_size_moving:time_bin_size_moving*(length(this_event_bins)-1);
% 
%     subplot(10,10,nEvent)
%     plot(x,(similarity_timeseries(this_event_bins)),'r')
%     hold on
%     y = mean(similarity_timeseries_shuffled(:,this_event_bins));
%     CI_U = prctile(similarity_timeseries_shuffled(:,this_event_bins),97.5,1);
%     CI_L = prctile(similarity_timeseries_shuffled(:,this_event_bins),2.5,1);
%     PLOT = plot(x,y,'Color','k');hold on;
%     ERROR_SHADE = patch([x fliplr(x)],[CI_U fliplr(CI_L)],'k','FaceAlpha','0.3','LineStyle','none');
% end
% 
% KDE_bias_shuffled_event=[];
% KDE_bias_event = [];
% KDE_bias = pdf_T1./(pdf_T1+pdf_T2);
% KDE_bias_shuffled = pdf_T1_shuffled./(pdf_T1_shuffled+pdf_T2_shuffled);
% 
% 
% % Visualize similarity
% figure;
% % similarity_timeseries_shuffled = reshape(similarity_timeseries_shuffled,1,[]);
% histogram([similarity_shuffled],100,'Normalization','probability');hold on;
% histogram(similarity,100,'Normalization','probability');
% xlabel('Ripple Time Bin');
% ylabel('Track 2 Similarity');
% title('Mahalanobis-Based Similarity');
% 
% 
% similarity_shuffled=[];
% similarity = [];
% for nEvent = 1:max(ripples_id);
%     this_event_bins = find(ripples_id==nEvent);
%     Track1_dist_shuffled(:,nEvent) = mean(dist_to_Track1_shuffled(:,this_event_bins),2);
%     Track2_dist_shuffled(:,nEvent) = mean(dist_to_Track2_shuffled(:,this_event_bins),2);
% 
%     Track1_dist(:,nEvent) = mean(dist_to_Track1(this_event_bins));
%     Track2_dist(:,nEvent) = mean(dist_to_Track2(this_event_bins));
% end
% 
% 
% histogram(Track1_dist,0:0.01:5,'Normalization','probability','EdgeAlpha',0);hold on;
% histogram(Track1_dist_shuffled,0:0.01:5,'Normalization','probability','EdgeAlpha',0);hold on;
% 
% 
% histogram(Track2_dist,0:0.01:5,'Normalization','probability','EdgeAlpha',0);hold on;
% histogram(Track2_dist_shuffled,0:0.01:5,'Normalization','probability','EdgeAlpha',0);
% 
% 
% similarity_shuffled=[];
% similarity = [];
% for nEvent = 1:max(ripples_id);
%     this_event_bins = find(ripples_id==nEvent);
%     similarity_shuffled(:,nEvent) = mean(similarity_timeseries_shuffled(:,this_event_bins),2);
%     similarity(nEvent) = mean(similarity_timeseries(this_event_bins));
% 
%     similarity_T1_sig = find(similarity<prctile(similarity_shuffled,2.5));
%     similarity_T2_sig = find(similarity>prctile(similarity_shuffled,97.5));
% end
% similarity_T1_sig = similarity<similarity_shuffled
% 
% similarity_T1_sig = find(similarity<prctile(similarity_shuffled,2.5));
% similarity_T2_sig = find(similarity>prctile(similarity_shuffled,97.5));
% 
% % Visualize similarity
% figure;
% % similarity_timeseries_shuffled = reshape(similarity_timeseries_shuffled,1,[]);
% histogram(similarity_shuffled(:,[similarity_T1_sig similarity_T2_sig]),200,'Normalization','probability');hold on;
% histogram(similarity([similarity_T1_sig similarity_T2_sig]),200,'Normalization','probability');
% % histogram(similarity_shuffled,200,'Normalization','probability');hold on;
% % histogram(similarity,200,'Normalization','probability');
% xlabel('Ripple Time Bin');
% ylabel('Track 2 Similarity');
% title('Nearest Neighbour Similarity');

