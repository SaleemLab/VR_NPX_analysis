function [RRR, RUN, RIP] = RRR_BlockDiag_Reactivation_decoding( ...
    tvec_template, position, speed, lap_times, track_ID, ...
    spikes_template, event_times, spike_target, options, varargin)
% RRR_BlockDiag_Reactivation_decoding
%
%   Ridge-regularised reduced-rank regression (RRR) with block-diagonal
%   position bases to decode track identity and position during RUN and to
%   project ripple activity into the same space to obtain track log-odds.
%
% Inputs:
%   tvec_template   - Time vector template for interpolation.
%   position        - Position data corresponding to the time vector.
%   speed           - Speed data corresponding to the time vector.
%   lap_times       - Lap times for different tracks (Nx2).
%   track_ID        - Track identifiers for each lap (N x 1, values 1 or 2).
%   spikes_template - Spike data template for activity templates (RUN cells).
%   event_times     - Event times for ripple events (Mx2).
%   spike_target    - Spike data for target neurons (ripple decoding).
%   options         - Struct containing additional options (e.g., SUBJECT, SESSION).
%
% Optional name/value pairs:
%   'event_type'          - default 'ripples'
%   'time_bin_size_RUN'   - default 0.1
%   'time_bin_size'       - default 0.02
%   'time_bin_size_moving'- default 0.01
%   'shuffle_option'      - 1 = do cell-ID shuffle CV (default), 0 = skip
%   'plot_option'         - 1 = basic plots, 0 = no plots
%
% Outputs:
%   RRR - struct with model and hyperparameters:
%       .lambda_opt, .rank_opt
%       .pos_bins
%       .B               (C x 2*nPos) final weight matrix
%       .muX, .sdX       (1 x C) z-scoring parameters
%       .cv_track_AUC    (nLambda x nRank)
%       .cv_pos_err_cm   (nLambda x nRank)
%       .shuffle_track_AUC  (nShuffle x 1)
%       .shuffle_pos_err_cm (nShuffle x 1)
%
%   RUN - struct with RUN decoding performance:
%       .track_true, .pos_true          (N_run x 1)
%       .CV.logodds, .CV.track_hat      (N_run x 1)
%       .CV.pos_hat_T1, .CV.pos_hat_T2
%       .CV.pos_hat_true
%       .CV.track_confmat, .CV.pos_confmat
%       .CV.track_AUC, .CV.mean_abs_pos_err_cm
%       .FULL.* same fields for model fit on all bins
%
%   RIP - struct with ripple projections:
%       .event_id        (N_rippleBins x 1)
%       .event_bins      (N_rippleBins x 2)
%       .logodds         (N_rippleBins x 1)
%       .track_hat       (N_rippleBins x 1)
%       .pos_hat_T1, .pos_hat_T2, .pos_hat
%       .event_peak_logodds (nEvents x 1)
%       .event_peak_idx      (nEvents x 1)
%
% Example:
%   [RRR, RUN, RIP] = RRR_BlockDiag_Reactivation_decoding( ...
%       tvec_template, position, speed, lap_times, track_ID, ...
%       spikes_template, event_times, spike_target, options);

%% ------------------------------------------------------------------------
%  Parse inputs
% -------------------------------------------------------------------------
p = inputParser;
addParameter(p,'event_type','ripples',@ischar);
addParameter(p,'time_bin_size_RUN',0.1,@isnumeric);
addParameter(p,'time_bin_size',0.02,@isnumeric);
addParameter(p,'time_bin_size_moving',0.01,@isnumeric);
addParameter(p,'shuffle_option',1,@isnumeric);
addParameter(p,'plot_option',1,@isnumeric);
% addParameter(p,'analysis_option',1,@isnumeric);% 1 is reactivation only, 2 is sequence

parse(p,varargin{:});
event_type        = p.Results.event_type;
time_bin_size_RUN = p.Results.time_bin_size_RUN;
time_bin_size     = p.Results.time_bin_size;
time_bin_size_moving = p.Results.time_bin_size_moving;
shuffle_option    = p.Results.shuffle_option;
plot_option       = p.Results.plot_option;

RRR = struct();
RUN = struct();
RIP = struct();

%% ------------------------------------------------------------------------
%  RUN binning and spike counts (same as your existing pipeline)
% -------------------------------------------------------------------------
% Time bins for RUN
timevec_edge_RUN = interp1(tvec_template, tvec_template, ...
    tvec_template(1):time_bin_size_RUN:tvec_template(end));

speed_RUN = interp1(tvec_template, speed, ...
    tvec_template(1):time_bin_size_RUN:tvec_template(end));
timevec_edge_RUN(speed_RUN < 1) = nan; % exclude immobility

samples = timevec_edge_RUN';
bins    = [samples, samples + time_bin_size_RUN];

% Position (discretised to 5cm bins over 28 bins)
bin_size = 10;
no_bins = 140/bin_size;

position_interp1 = interp1(tvec_template, position, timevec_edge_RUN);
position_interp1(isnan(position_interp1)) = 0;
position_interp1 = discretize(position_interp1, no_bins) * bin_size;

% Track-specific RUN bins
[T1_bins,~,~] = InIntervals(timevec_edge_RUN', lap_times(track_ID==1,:));
[T2_bins,~,~] = InIntervals(timevec_edge_RUN', lap_times(track_ID==2,:));

% RUN spike counts: n is (N_runBins x nCells_template)
[~,~,~,~,~,n] = ActivityTemplates(spikes_template,'bins',bins);

%% ------------------------------------------------------------------------
%  Ripple binning and spike counts 
% -------------------------------------------------------------------------
ripples_time_edges = [];
ripples_id = [];

for iEvt = 1:size(event_times,1)
    event_duration = event_times(iEvt,2) - event_times(iEvt,1);

    if contains(event_type,'ripples')
        if iEvt < size(event_times,1)
            if event_duration < 0.1 && (event_times(iEvt,1) + 0.1 < event_times(iEvt+1,1))
                event_duration = 0.1;
            end
        elseif event_duration < 0.1
            event_duration = 0.1;
        end
    end

    num_bins = floor(event_duration / time_bin_size_moving);
    timebins_edges = linspace(event_times(iEvt,1), ...
                              event_times(iEvt,1) + num_bins*time_bin_size_moving, ...
                              num_bins+1);

    ripples_time_edges = [ripples_time_edges, timebins_edges];
    ripples_id = [ripples_id, iEvt*ones(1,length(timebins_edges))];
end

ripple_bins = [];
ripple_bins(:,1) = ripples_time_edges(:);
ripple_bins(:,2) = ripples_time_edges(:) + time_bin_size;

% Spike counts for ripple bins
if sum(event_times(2:end,1) - event_times(1:end-1,2) < 0) < 1
    [~,~,~,~,~,n_ripples] = ActivityTemplates(spike_target,'bins',ripple_bins);
else
    % In case larger window / overlapping bins: manual counts
    spike_times = spike_target(:,1);
    unit_ids    = spike_target(:,2);
    nNeurons    = max(unit_ids);

    spike_counts_all = cell(1,size(event_times,1));

    parfor iEvt = 1:size(event_times,1)
        this_bins = ripple_bins(ripples_id == iEvt,:);
        this_counts = zeros(size(this_bins,1), nNeurons);

        for unit = 1:nNeurons
            unit_spikes = spike_times(unit_ids == unit);
            this_counts(:,unit) = CountInIntervals(unit_spikes, this_bins);
        end

        spike_counts_all{iEvt} = this_counts;
    end

    spike_counts_all = cat(1,spike_counts_all{:});
    n_ripples = zscore(spike_counts_all, 0, 1);
end

% Ripple residuals: regress out global pattern
mean_ripple_activity = mean(n_ripples, 2); % N_rippleBins x 1
coeff = (mean_ripple_activity' * n_ripples) / ...
        (mean_ripple_activity' * mean_ripple_activity);
ripple_global_pattern = mean_ripple_activity * coeff;
n_ripples_residuals = zscore(n_ripples - ripple_global_pattern, 0, 1);

%% ------------------------------------------------------------------------
%  Build RUN design matrices for RRR 
% -------------------------------------------------------------------------
X_run_full = [n(T1_bins,:); n(T2_bins,:)]; % N_run x C
Y_track    = [ones(sum(T1_bins),1); 2*ones(sum(T2_bins),1)]; % N_run x 1

pos_run_full = [position_interp1(T1_bins) position_interp1(T2_bins)]';
pos_run_full = pos_run_full(:); % N_run x 1

N_run   = size(X_run_full,1);
nCells  = size(X_run_full,2);

% Z-score spikes across RUN bins (per neuron)
muX = mean(X_run_full,1);
sdX = std(X_run_full,[],1);
sdX(sdX==0) = 1;
X_run_z = (X_run_full - muX) ./ sdX;

% Position bins (used to build one-hot basis)
pos_bins = unique(pos_run_full);
nPos     = numel(pos_bins);

% Position indices per RUN bin
[~,pos_idx] = ismember(pos_run_full,pos_bins);

% Build block-diagonal bases Φ_L, Φ_R
Phi_L = zeros(N_run, nPos);
Phi_R = zeros(N_run, nPos);

idx_T1 = (Y_track == 1);
idx_T2 = (Y_track == 2);

Phi_L(sub2ind(size(Phi_L), find(idx_T1), pos_idx(idx_T1))) = 1;
Phi_R(sub2ind(size(Phi_R), find(idx_T2), pos_idx(idx_T2))) = 1;

Y_basis = [Phi_L Phi_R]; % N_run x (2*nPos)

%% ------------------------------------------------------------------------
%  Hyperparameter search: ridge-regularised RRR over (lambda, rank)
% -------------------------------------------------------------------------
lambda_grid = logspace(-2,4,30);      % you can tweak this
max_rank    = min([nCells, 2*nPos 10]);
rank_grid   = 1:max_rank;

nLambda = numel(lambda_grid);
nRank   = numel(rank_grid);

cv_track_AUC  = nan(nLambda, nRank);
cv_pos_err_cm = nan(nLambda, nRank);

rng(200);
Kfold = 10;
cv = cvpartition(Y_track, 'KFold', Kfold);
tic
for il = 1:nLambda
    lambda = lambda_grid(il);
    for ir = 1:nRank
        r = rank_grid(ir);

        AUC_fold = nan(Kfold,1);
        err_fold = nan(Kfold,1);
        err_fold_false = nan(Kfold,1);

        for f = 1:Kfold
            trIdx = training(cv,f);
            teIdx = test(cv,f);

            Xtr = X_run_z(trIdx,:);
            Xte = X_run_z(teIdx,:);
            Ytr = Y_basis(trIdx,:);
            pos_true   = pos_run_full(teIdx);
            track_true = Y_track(teIdx);

            % Fit ridge+RRR
            B = fit_ridge_rrr(Xtr, Ytr, lambda, r);

            % Predict
            Yhat = Xte * B;

            [logodds, track_hat, pos_hat_T1, pos_hat_T2, ~] = ...
                decode_track_pos_from_Yhat(Yhat,pos_bins);

            % ROC for log odds 
            fpr = 0:0.01:1;
            log_odds = logodds;
  
            %                 [~, ~, ~, AUC_real(i)] = perfcurve(track_true(idx), log_odds(idx),1,'XVals', fpr);
            [x,y,~,AUC_real] = perfcurve(track_true, log_odds,1, 'XVals', fpr);

            
            AUC_fold(f) = mean(AUC_real);

            % Position error using true track
            pos_pred_true = pos_hat_T1;
            pos_pred_true(track_true==2) = pos_hat_T2(track_true==2);
            err_fold(f) = mean(abs(pos_pred_true - pos_true));

            pos_pred_false = pos_hat_T2;
            pos_pred_false(track_true==2) = pos_hat_T1(track_true==2);
            err_fold_false(f) = mean(abs(pos_pred_false - pos_true));
        end

        cv_track_AUC(il,ir)  = mean(AUC_fold);
        cv_pos_err_cm(il,ir) = mean(err_fold);
        cv_pos_err_cm_false(il,ir) = mean(err_fold_false);
    end
end

% Choose best (lambda, rank): maximise track AUC, then minimize pos error
[~,idx_best] = max(cv_track_AUC(:));
[best_il, best_ir] = ind2sub(size(cv_track_AUC), idx_best);
lambda_opt = lambda_grid(best_il);
rank_opt   = rank_grid(best_ir);

fprintf('RRR: chosen lambda = %.3g, rank = %d, CV track AUC = %.3f, CV pos err = %.1f cm\n', ...
    lambda_opt, rank_opt, cv_track_AUC(best_il,best_ir), cv_pos_err_cm(best_il,best_ir));
toc



%% ------------------------------------------------------------------------
%  Cross-validated RUN predictions using (lambda_opt, rank_opt)
% -------------------------------------------------------------------------
logodds_cv      = nan(N_run,1);
% track_hat_cv    = nan(N_run,1);
pos_hat_T1_cv   = nan(N_run,1);
pos_hat_T2_cv   = nan(N_run,1);
pos_hat_true_cv = nan(N_run,1);

for f = 1:Kfold
    trIdx = training(cv,f);
    teIdx = test(cv,f);

    Xtr = X_run_z(trIdx,:);
    Xte = X_run_z(teIdx,:);
    Ytr = Y_basis(trIdx,:);

    B = fit_ridge_rrr(Xtr, Ytr, lambda_opt, rank_opt);
    Yhat = Xte * B;

    [logodds, track_hat, pos_hat_T1, pos_hat_T2, ~] = ...
        decode_track_pos_from_Yhat(Yhat,pos_bins);

    logodds_cv(teIdx)    = logodds;
%     track_hat_cv(teIdx)  = track_hat;
    pos_hat_T1_cv(teIdx) = pos_hat_T1;
    pos_hat_T2_cv(teIdx) = pos_hat_T2;

    pos_hat_true = pos_hat_T1;
    pos_hat_true(Y_track(teIdx)==2) = pos_hat_T2(Y_track(teIdx)==2);
    pos_hat_true_cv(teIdx) = pos_hat_true;

    pos_hat_false = pos_hat_T1;
    pos_hat_false(Y_track(teIdx)==1) = pos_hat_T2(Y_track(teIdx)==1);
    pos_hat_false_cv(teIdx) = pos_hat_false;
end


%% ------------------------------------------------------------------------
%  Cell-ID shuffle CV on RUN (using lambda_opt, rank_opt)
% -------------------------------------------------------------------------
if shuffle_option
    nShuffle = 1000;
else
    nShuffle = 0;
end


shuffle_pos_err_cm = nan(nShuffle,1);
logodds_shuffled      = nan(nShuffle,N_run);

for nsh = 1:nShuffle
    s = RandStream('mrg32k3a','Seed',nsh+1e6);
    perm = randperm(s,nCells);
    X_sh = X_run_z(:,perm);

    AUC_fold = nan(Kfold,1);
    err_fold = nan(Kfold,1);

    for f = 1:Kfold
        trIdx = training(cv,f);
        teIdx = test(cv,f);

        Xtr = X_run_z(trIdx,:);
        Xte = X_sh(teIdx,:);
        Ytr = Y_basis(trIdx,:);

        B = fit_ridge_rrr(Xtr, Ytr, lambda_opt, rank_opt);
        Yhat = Xte * B;

        [logodds, track_hat, pos_hat_T1, pos_hat_T2, ~] = ...
            decode_track_pos_from_Yhat(Yhat,pos_bins);

        track_true = Y_track(teIdx);
        pos_true   = pos_run_full(teIdx);

        fpr = 0:0.01:1;
        logodds_shuffled(nsh,teIdx) = logodds;

        pos_pred_true = pos_hat_T1;
        pos_pred_true(track_true==2) = pos_hat_T2(track_true==2);
        err_fold(f) = mean(abs(pos_pred_true - pos_true));
    end

%     shuffle_track_AUC(nsh)  = mean(AUC_fold);
    shuffle_pos_err_cm(nsh) = mean(err_fold);
end


% ROC for log odds
fpr = 0:0.01:1;
TPR_real = zeros(1000, length(fpr));
track_AUC_cv = zeros(1,1000);
log_odds = logodds_cv;
parfor i = 1:1000
    s = RandStream('mrg32k3a','Seed',i); % Set random seed for resampling
    idx = randsample(s,length(Y_track), length(Y_track), true);
    %                 [~, ~, ~, AUC_real(i)] = perfcurve(track_true(idx), log_odds(idx),1,'XVals', fpr);
    [x,y,~,track_AUC_cv(i)] = perfcurve(Y_track(idx), log_odds(idx),1, 'XVals', fpr);
    [C,ia,ic] = unique(x);
    y = y(ia);
    x = x(ia);

    TPR_interp = interp1(x, y, fpr, 'linear', 'extrap');  % force to match fpr
    TPR_real(i,:) = TPR_interp;
end


% ROC for log odds cell id shuffled
fpr = 0:0.01:1;
TPR_cellID_shuf = zeros(1000, length(fpr));
cellID_shuffle_track_AUC = zeros(1,1000);
log_odds = logodds_shuffled;
% log_odds = logodds_cv;
parfor i = 1:1000
    s = RandStream('mrg32k3a','Seed',i); % Set random seed for resampling
    idx = randsample(s,length(Y_track), length(Y_track), true);
%     [x,y,~,cellID_shuffle_track_AUC(i)] = perfcurve(Y_track(idx), log_odds,1, 'XVals', fpr);
    [x,y,~,cellID_shuffle_track_AUC(i)] = perfcurve(Y_track(idx), log_odds(idx),1, 'XVals', fpr);
    [C,ia,ic] = unique(x);
    y = y(ia);
    x = x(ia);

    TPR_interp = interp1(x, y, fpr, 'linear', 'extrap');  % force to match fpr
    TPR_cellID_shuf(i,:) = TPR_interp;
end

% ROC for log odds track label shuffled
fpr = 0:0.01:1;
TPR_shuf = zeros(1000, length(fpr));
shuffle_track_AUC = zeros(1,1000);
% log_odds = logodds_shuffled;
log_odds = logodds_cv;
parfor i = 1:1000
    s = RandStream('mrg32k3a','Seed',i); % Set random seed for resampling
    idx = randsample(s,length(Y_track), length(Y_track), true);
    [x,y,~,shuffle_track_AUC(i)] = perfcurve(Y_track(idx), log_odds,1, 'XVals', fpr);
%     [x,y,~,shuffle_track_AUC(i)] = perfcurve(Y_track(idx), log_odds(idx),1, 'XVals', fpr);
    [C,ia,ic] = unique(x);
    y = y(ia);
    x = x(ia);

    TPR_interp = interp1(x, y, fpr, 'linear', 'extrap');  % force to match fpr
    TPR_shuf(i,:) = TPR_interp;
end



mean_abs_pos_err_cv = mean(abs(pos_hat_true_cv - pos_run_full));

% Confusion matrices
[~, true_pos_idx] = ismember(pos_run_full,pos_bins);
[~, pred_pos_idx_true] = ismember(pos_hat_true_cv, pos_bins);

% track_confmat_cv = confusionmat(Y_track, track_hat_cv);

pos_confmat_T1   = confusionmat(true_pos_idx(track_true ==1), pos_hat_T1_cv(track_true ==1)/bin_size);
pos_confmat_T2   = confusionmat(true_pos_idx(track_true ==2), pos_hat_T2_cv(track_true ==2)/bin_size);


if plot_option
    title_text = sprintf('%s %s %s RRR ROC validation %.0f ms bins', ...
        options.SUBJECT, options.SESSION, event_type, time_bin_size_RUN*1000);

    fig = figure('Name', title_text, 'Position', [200 100 866 741]); hold on;

    nexttile
    CI_real = prctile(TPR_real, [0.5 99.5]);
    CI_shuf = prctile(TPR_shuf, [0.5 99.5]);
    hold on
    plot([0 1],[0 1],'--k')
    PLOT(1) = fill([fpr fliplr(fpr)], [CI_real(2,:) fliplr(CI_real(1,:))], ...
        [231,41,138]/256, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    PLOT(2) = fill([fpr fliplr(fpr)], [CI_shuf(2,:) fliplr(CI_shuf(1,:))], ...
        'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    plot(fpr, mean(TPR_real,1), 'Color', [231,41,138]/256);
    plot(fpr, mean(TPR_shuf,1), 'k', 'LineWidth', 1.5);
    xlabel('True positive rate')
    ylabel('False Positive rate')
    %     xline(0.5, '--k'); xlabel('False Positive Rate'); ylabel('True Positive Rate');
    title(title_text); legend(PLOT(1:2),{ 'Real', 'Shuffle'},'box', 'off');
    set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);

    nexttile

    AUC = [track_AUC_cv;shuffle_track_AUC];
    bar_colors = [231,41,138; 0, 0, 0]/255;
    x_pos=[1 2];
    for i = 1:2

        mean_AUC = mean(AUC(i,:),'omitnan');
        auc_errors = [mean_AUC-prctile(AUC(i,:), [0.5]),...  % lower error
            prctile(AUC(i,:), [99.5])-mean_AUC];   % upper error
        hold on
        bar(x_pos(i), mean_AUC, 0.4, 'FaceAlpha',0.5, 'FaceColor', bar_colors(i,:), 'EdgeColor', 'none');
        errorbar(x_pos(i), mean_AUC, auc_errors(1), auc_errors(2), 'k', 'linestyle', 'none', 'linewidth', 1);
    end


    xticks([1 2]);
    xticklabels({'Real', 'Shuffled'});
    ylabel('AUC');
    ylim([0 1]);
    title('Mean AUC – Real vs Shuffled');
    set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);


    nexttile
    CI_real = prctile(TPR_real, [0.5 99.5]);
    CI_shuf = prctile(TPR_cellID_shuf, [0.5 99.5]);
    hold on
    plot([0 1],[0 1],'--k')
    PLOT(1) = fill([fpr fliplr(fpr)], [CI_real(2,:) fliplr(CI_real(1,:))], ...
        [231,41,138]/256, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    PLOT(2) = fill([fpr fliplr(fpr)], [CI_shuf(2,:) fliplr(CI_shuf(1,:))], ...
        'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    plot(fpr, mean(TPR_real,1), 'Color', [231,41,138]/256);
    plot(fpr, mean(TPR_cellID_shuf,1), 'k', 'LineWidth', 1.5);
    xlabel('True positive rate')
    ylabel('False Positive rate')
    %     xline(0.5, '--k'); xlabel('False Positive Rate'); ylabel('True Positive Rate');
    title(title_text); legend(PLOT(1:2),{ 'Real', 'Shuffle'},'box', 'off');
    set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);

    nexttile

    AUC = [track_AUC_cv;cellID_shuffle_track_AUC];
    bar_colors = [231,41,138; 0, 0, 0]/255;
    x_pos=[1 2];
    for i = 1:2

        mean_AUC = mean(AUC(i,:),'omitnan');
        auc_errors = [mean_AUC-prctile(AUC(i,:), [0.5]),...  % lower error
            prctile(AUC(i,:), [99.5])-mean_AUC];   % upper error
        hold on
        bar(x_pos(i), mean_AUC, 0.4, 'FaceAlpha',0.5, 'FaceColor', bar_colors(i,:), 'EdgeColor', 'none');
        errorbar(x_pos(i), mean_AUC, auc_errors(1), auc_errors(2), 'k', 'linestyle', 'none', 'linewidth', 1);
    end


    xticks([1 2]);
    xticklabels({'Real', 'Cell ID Shuffled'});
    ylabel('AUC');
    ylim([0 1]);
    title('Mean AUC – Real vs Shuffled');
    set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);



end


%% ------------------------------------------------------------------------
%  Final RUN model fit on all bins (FULL)
% -------------------------------------------------------------------------
B_final = fit_ridge_rrr(X_run_z, Y_basis, lambda_opt, rank_opt);
Yhat_full = X_run_z * B_final;

%% ------------------------------------------------------------------------
%  Project new data with RRR model to get log-odds and positions
% -------------------------------------------------------------------------
X_rip_z = (n_ripples_residuals - muX) ./ sdX;
Yhat_rip = X_rip_z * B_final;

[logodds_rip, track_hat_rip, pos_hat_T1_rip, pos_hat_T2_rip, pos_hat_rip] = ...
    decode_track_pos_from_Yhat(Yhat_rip,pos_bins);

event_id_vec = ripples_id(:);
event_bins   = ripple_bins;

nEvents = max(event_id_vec);
event_peak_logodds = nan(nEvents,1);
event_peak_idx     = nan(nEvents,1);

for e = 1:nEvents
    idx = find(event_id_vec == e);
    if isempty(idx), continue; end

    % Peak by absolute log-odds
    [~,rel] = max(abs(logodds_rip(idx)));
    event_peak_logodds(e) = logodds_rip(idx(rel));
    event_peak_idx(e)     = idx(rel);
end

%% ------------------------------------------------------------------------
%  Pack outputs
% -------------------------------------------------------------------------
RRR.method            = 'ridge_RRR_blockdiag_position';
RRR.lambda_opt        = lambda_opt;
RRR.rank_opt          = rank_opt;
RRR.pos_bins          = pos_bins;
RRR.B                 = B_final;
RRR.muX               = muX;
RRR.sdX               = sdX;
RRR.cv_track_AUC      = cv_track_AUC;
RRR.cv_pos_err_cm     = cv_pos_err_cm;
RRR.cv_pos_err_cm_false     = cv_pos_err_cm_false;
RRR.lambda_grid       = lambda_grid;
RRR.rank_grid         = rank_grid;

RUN.track_true        = Y_track;
RUN.pos_true          = pos_run_full;
RUN.pos_bins          = pos_bins;

RUN.logodds        = logodds_cv;
RUN.track_hat      = track_hat_cv;
RUN.pos_hat_T1     = pos_hat_T1_cv;
RUN.pos_hat_T2     = pos_hat_T2_cv;
RUN.pos_hat_true   = pos_hat_true_cv;
RUN.track_confmat  = track_confmat_cv;
RUN.pos_confmat    = pos_confmat_cv;
RUN.track_AUC      = track_AUC_cv;
RUN.mean_abs_pos_err_cm = mean_abs_pos_err_cv;
RUN.track_AUC_shuffled = shuffle_track_AUC;
RUN.mean_pos_err_cm_shuffled = shuffle_pos_err_cm;


RIP.event_id          = event_id_vec;
RIP.event_bins        = event_bins;
RIP.logodds           = logodds_rip;
RIP.track_hat         = track_hat_rip;
RIP.pos_hat_T1        = pos_hat_T1_rip;
RIP.pos_hat_T2        = pos_hat_T2_rip;
RIP.pos_hat           = pos_hat_rip;
RIP.event_peak_logodds= event_peak_logodds;
RIP.event_peak_idx    = event_peak_idx;
RIP.pos_bins          = pos_bins;

%% ------------------------------------------------------------------------
%  (Optional) minimal plotting
% -------------------------------------------------------------------------
if plot_option
    try
        % Track decoding performance vs shuffle
        figure('Name','RRR RUN decoding performance','Color','w');
        subplot(1,2,1);
        histogram(shuffle_track_acc,30,'FaceAlpha',0.5,'EdgeColor','none');
        hold on;
        yl = ylim;
        plot([track_AUC track_acc_full], yl, 'r','LineWidth',2);
        xlabel('Track decoding accuracy'); ylabel('Count');
        title('Track decoding (shuffle vs real)');
        box off; set(gca,'TickDir','out');

        subplot(1,2,2);
        histogram(shuffle_pos_err_cm,30,'FaceAlpha',0.5,'EdgeColor','none');
        hold on;
        yl = ylim;
        plot([mean_abs_pos_err_full mean_abs_pos_err_full], yl, 'r','LineWidth',2);
        xlabel('Mean |pos error| (cm)'); ylabel('Count');
        title('Position error (shuffle vs real)');
        box off; set(gca,'TickDir','out');

        % Ripple log-odds examples
        figure('Name','RRR ripple log-odds','Color','w');
        ev_ids = unique(event_id_vec);
        nShow = min(20, numel(ev_ids));
        for k = 1:nShow
            e = ev_ids(k);
            idx = find(event_id_vec==e);
            t = (0:numel(idx)-1) * time_bin_size_moving;
            subplot(4,5,k);
            plot(t, logodds_rip(idx));
            yline(0,'k--');
            xlabel('Time (s)');
            ylabel('log P(T1) - log P(T2)');
            title(sprintf('Event %d',e));
            box off; set(gca,'TickDir','out');
        end
    catch
        warning('Plotting failed, but decoding completed successfully.');
    end
end

fprintf('RRR BlockDiag decoding finished.\n');

end

%% ------------------------------------------------------------------------
%  Subfunctions
% -------------------------------------------------------------------------
function B = fit_ridge_rrr(X, Y, lambda, r)
% X: N x C, Y: N x D
% Ridge solution + reduced-rank approximation via SVD.

[~,C] = size(X);
Sxx = (X' * X) + lambda * eye(C);
B_ridge = Sxx \ (X' * Y);   % C x D

[U,S,V] = svd(B_ridge,'econ');
max_r = min([r, size(S,1), size(S,2)]);
B = U(:,1:max_r) * S(1:max_r,1:max_r) * V(:,1:max_r)';

end

function [track_logodds, track_hat, pos_hat_T1, pos_hat_T2, pos_hat_all] = ...
         decode_track_pos_from_Yhat(Yhat,pos_bins)
% Decode track log-odds and track-conditional positions from Yhat = X*B.

nPos = numel(pos_bins);
Y_L  = Yhat(:,1:nPos);
Y_R  = Yhat(:,nPos+1:end);

% Rectify to get non-negative pseudo-probabilities
Y_L(Y_L<0) = 0;
Y_R(Y_R<0) = 0;

sumL = sum(Y_L,2) + eps;
sumR = sum(Y_R,2) + eps;

pL = sumL ./ (sumL + sumR);
pR = sumR ./ (sumL + sumR);

track_logodds = log(pL+eps) - log(pR+eps);
track_hat = ones(size(track_logodds));
track_hat(track_logodds < 0) = 2;

% Track-conditional position distributions
p_pos_L = Y_L ./ sumL; % row-normalised
p_pos_R = Y_R ./ sumR;

pos_bins_col = pos_bins(:);
% pos_hat_T1 = p_pos_L * pos_bins_col;
% pos_hat_T2 = p_pos_R * pos_bins_col;

% maximum-probability position estimate (MAP)

[~, idx] = max(p_pos_L');        % index of highest probability bin Track 1
pos_hat_T1 = pos_bins_col(idx);  % corresponding position estimate

[~, idx] = max(p_pos_R');        % index of highest probability bin Track 1
pos_hat_T2 = pos_bins_col(idx);  % corresponding position estimate


pos_hat_all = pos_hat_T1;
pos_hat_all(track_hat==2) = pos_hat_T2(track_hat==2);

end
