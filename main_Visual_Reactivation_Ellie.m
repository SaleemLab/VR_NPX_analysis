

clear all
addpath(genpath('C:\Users\eleanor.benoit\Documents\GitHub\VR_NPX_analysis'))

SUBJECTS = {'M00069'};  %% set this - 1/4
option = 'V1-HPC';
experiment_info = subject_session_stimuli_mapping_Ellie(SUBJECTS, option);
params = create_cluster_selection_params('sorting_option','ellie');
psthBinSize = 0.01; % this script divides this by 10 (to 0.001s) for raster plots

%% 2/4
Stimulus_type = 'Sleep_Box_2'; % 'Sleep_Box' 'Sleep_Box_1' 'Sleep_Box_2' 'Sleep_Box_3'
cd('V:\Ellie\DATA\SUBJECTS\M00069\analysis\20250929\GAVNIK_ABCD_1')
load('pls_V1_20250929_GAVNIK_ABCD_1.mat')
pls_model.X

cd('V:\Ellie\DATA\SUBJECTS\M00069\analysis\20250929\Sleep_Box_2')
load('SleepData.mat')
SleepData.V1.X_sleep

[qr_W,~] = qr(pls_model.XL,0);
%%% Full PLS model during visual stimulation
X = pls_model.X;
popMean = mean(X,2);
Xreg = [ones(size(popMean)) popMean];
X_residual = zeros(size(X));

for ncell = 1:size(X,2)
    beta = Xreg \ X(:,ncell);
    X_residual(:,ncell) = X(:,ncell) - Xreg*beta;
end
XS_train_plot = X_residual * qr_W(:,1:PLS_for_centroids);

stimulus_template = [];
figure;
hold on;
for nstimulus = 1:size(pls_model.Y,2)
    stimulus_template{nstimulus} = XS_train_plot(pls_model.Y(:,nstimulus) == 1,:);
    scatter3(stimulus_template{nstimulus}(:,1),stimulus_template{nstimulus}(:,2),stimulus_template{nstimulus}(:,3))
end
% axis equal

%%% ripple data projection
[nEvents,nTimebins,nNeurons] = size(SleepData.V1.X_sleep);
X_sleep_2D = reshape(SleepData.V1.X_sleep,...
                     nEvents*nTimebins,...
                     nNeurons);
popMean = mean(X_sleep_2D,2);
Xreg = [ones(size(popMean)) popMean];
X_residual = zeros(size(X_sleep_2D));
for ncell = 1:nNeurons
    beta = Xreg \ X_sleep_2D(:,ncell);
    X_residual(:,ncell) = X_sleep_2D(:,ncell) - Xreg*beta;
end

XS_sleep_plot = X_residual * qr_W(:,1:PLS_for_centroids);

% 
% XS_sleep_plot_3D = reshape(XS_sleep_plot,...
%                            nEvents,...
%                            nTimebins,...
%                            PLS_for_centroids);



%%%% KDE cross validation to get optimal bandwidth

nFold = 5;
bandwidth_list = linspace(0.1,1,20);
nStim = 4;
nPLS = pls_model.n_components;

acc = zeros(nFold,length(bandwidth_list));
rng(200);
cv = cvpartition(size(X,1),'KFold',nFold);

for f = 1:nFold

    trainIdx = training(cv,f);
    testIdx = test(cv,f);

    X_train = X_residual(trainIdx,:);
    X_test  = X_residual(testIdx,:);
    Y_train = Y(trainIdx,:);   % one-hot
    Y_test  = Y(testIdx,:);    % one-hot
    [~,Y_train_lab] = max(Y_train,[],2);
    [~,Y_test_lab]  = max(Y_test,[],2);

    [XL,~,~,~,~,~,~,~] = plsregress(X_train,Y_train,nPLS);

    [qr_W,~] = qr(XL,0);

    Z_train = X_train * qr_W;
    Z_test  = X_test  * qr_W;

    stim_data = cell(nStim,1);

    for s = 1:nStim
        stim_data{s} = Z_train(Y_train_lab == s,:);
    end

    for b = 1:length(bandwidth_list)

        h = bandwidth_list(b);

        nTest = size(Z_test,1);
        logLik = zeros(nTest,nStim);

        for s = 1:nStim

            Xs = stim_data{s};

            for i = 1:nTest
                x = Z_test(i,:);
                logLik(i,s) = log(mvksdensity(Xs,x,'Bandwidth',h));
            end
        end

        logLik = logLik - max(logLik,[],2);

        P = exp(logLik);
        P = P ./ sum(P,2);

        [~,pred] = max(P,[],2);

        acc(f,b) = mean(pred == Y_test_lab);
    end
end

mean_acc = mean(acc,1);
[~,best_idx] = max(mean_acc);
best_bandwidth = bandwidth_list(best_idx);
figure
plot(bandwidth_list,mean_acc);


%%% Now project ripple data
nStim = numel(stimulus_template);
nSleep = size(XS_sleep_plot,1);

logLik_sleep = zeros(nSleep,nStim);

h = best_bandwidth;

for s = 1:nStim

    Xs = stimulus_template{s};

    for i = 1:nSleep

        x = XS_sleep_plot(i,:);

        logLik_sleep(i,s) = log(mvksdensity(Xs,x,'Bandwidth',h));

    end
end

% convert to probabilities
logLik_sleep = logLik_sleep - max(logLik_sleep,[],2);

P_sleep = exp(logLik_sleep);
P_sleep = P_sleep ./ sum(P_sleep,2);

% winner stimulus per timebin
[~,sleep_stim_pred] = max(P_sleep,[],2);


figure;


nStimulus = size(pls_model.Y,2);

cmap = lines(nStimulus);

h = gobjects(nStimulus,2);

for nstimulus = 1:nStimulus
    nexttile
    hold on;
    col = cmap(nstimulus,:);

    h(nstimulus,1) = scatter3( ...
        stimulus_template{nstimulus}(:,1), ...
        stimulus_template{nstimulus}(:,2), ...
        stimulus_template{nstimulus}(:,3), ...
        30, col, 'o', 'filled');

    idx = (sleep_stim_pred == nstimulus);

    h(nstimulus,2) = scatter3( ...
        XS_sleep_plot(idx,1), ...
        XS_sleep_plot(idx,2), ...
        XS_sleep_plot(idx,3), ...
        10, col, 'd');
end

axis equal;
grid on;
xlabel('PC1');
ylabel('PC2');
zlabel('PC3');

legend(h(:,1), strcat("Stim ", string(1:nStimulus)), 'Location','best');


figure;
hold on;
for nstimulus = 1:size(pls_model.Y,2)
    % stimulus_template{nstimulus} = XS_train_plot(pls_model.Y(:,nstimulus) == 1,:);
    scatter3(stimulus_template{nstimulus}(:,1),stimulus_template{nstimulus}(:,2),stimulus_template{nstimulus}(:,3))
end

cmap = lines(6);
col = [0 0 0];
scatter3(XS_sleep_plot(:,1),XS_sleep_plot(:,2),XS_sleep_plot(:,3),2, col, 'd','MarkerFaceAlpha', 0.1)

sgtitle(sprintf(['%s Stimulus-Evoked Activity and Peri-Ripple %s Activity in PLS Latent Space Day %d'], ...
                                 SleepData.subject, Stimulus_type, experiment_info(nsession).date), ...
                        'Interpreter', 'none');