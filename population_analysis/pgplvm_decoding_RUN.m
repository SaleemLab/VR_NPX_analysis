function [model,folds] = pgplvm_decoding_RUN( ...
    spike_times, spike_ids, ...
    tvec, position, speed, track_id, opts)
% P-GPLVM (LMT-style) training + constrained test inference + KNN decode
%  - AR1/OU prior on time for X(t) (latentTYPE=1 in LMT)
%  - SE prior on tuning f_n(X) with Poisson likelihood (Laplace)
%  - 5-fold temporally contiguous CV
%  - Final model trained on all RUN bins (for decoding new data)
%
% Inputs:
%   spike_times [nSp x 1] (s)
%   spike_ids   [nSp x 1] unit IDs
%   tvec        [T x 1]   continuous timebase (s)
%   position    [T x 1]
%   speed       [T x 1]
%   track_id    [T x 1]
%   opts        struct: binSize_train(0.1), binSize_test(0.02),
%                      runSpeedThresh(1), latentDim(1),
%                      nIter_train(60), nIter_infer(50), KNN_K(25),
%                      seed(1)
%
% Output:
%   out.cv.folds(f): per-fold decoding (t, X, pos_true/hat, trk_aff, rmse, mae)
%   out.cv.summary:  averaged CV metrics (rmse, mae)
%   out.model:       trained-on-all RUN model + banks for decoding
%   out.meta:        metadata (unit_ids, edges, isRun, opts)

%% ---- options
if nargin < 7, opts = struct; end
opts.binSize_train    = fget(opts,'binSize_train',0.1);
opts.binSize_test     = fget(opts,'binSize_test',0.020);
opts.runSpeedThresh   = fget(opts,'runSpeedThresh',1);
opts.latentDim        = fget(opts,'latentDim',1);
% opts.KNN_K            = fget(opts,'KNN_K',25);
opts.nIter_train      = fget(opts,'nIter_train',30);
opts.nIter_infer      = fget(opts,'nIter_infer',5);
opts.seed             = fget(opts,'seed',1);
rng(opts.seed);

%% ---- binning & RUN selection (mean speed per bin; require >50% valid samples)
[tvec, position, speed, track_id] = makecol(tvec, position, speed, track_id);
edges   = (tvec(1):opts.binSize_train:tvec(end))';
centers = edges(1:end-1) + opts.binSize_train/2;
contBinIdx = discretize(tvec, edges);

[Y, unit_ids] = bin_counts(spike_times(:), spike_ids(:), edges);   % [T x N]
[T, N] = size(Y); %#ok<NASGU>

isValid = ~isnan(position) & ~isnan(speed) & ismember(track_id,[1 2]);
bin_frac_valid = accum_by_bin(contBinIdx, isValid, @(x) mean(x)>0.5, true);
bin_speed_mean = accum_by_bin(contBinIdx, speed, @mean, false);
isRun = bin_frac_valid & (bin_speed_mean >= opts.runSpeedThresh);

pos_bin = accum_by_bin(contBinIdx, position, @mean, false);
trk_bin = mode_by_bin(contBinIdx, track_id);

Y_run   = Y(isRun,:);      % [T_run x N]
t_run   = centers(isRun);  % [T_run x 1]
pos_run = pos_bin(isRun);  % [T_run x 1]
trk_run = trk_bin(isRun);  % [T_run x 1]
T_run   = size(Y_run,1);


%% ---- 10 contiguous folds (train on training bins, infer latent on test bins)

nFolds     = 5;
idx_all    = (1:T_run)';
fold_edges = round(linspace(1, T_run+1, nFolds+1));
folds      = struct([]);

K_list = [10 15 20 25 30 50 75 100];

nK     = numel(K_list);

for f = 1:nFolds
    te = idx_all(fold_edges(f):fold_edges(f+1)-1);
    tr = setdiff(idx_all, te);

    Ytr   = Y_run(tr,:);  Yte   = Y_run(te,:);
    ttr   = t_run(tr);    tte   = t_run(te);
    postr = pos_run(tr);  poste = pos_run(te);
    trktr = trk_run(tr);  trkte = trk_run(te);

    % ---- Train LMT-style P-GPLVM on TRAIN -------------------------------
    model_tr = train_lmt_pgplvm(Ytr, opts.latentDim, opts.nIter_train);
    % %
    % figure;
    % hold on;
    % ax1 = axes;
    % scatter3(ax1,model_tr.X_train(trktr==1,1),model_tr.X_train(trktr==1,2),model_tr.X_train(trktr==1,3),20,postr(trktr==1),'filled');
    % colormap(ax1, abyss);          % colormap for first scatter
    % set(ax1, 'Color', 'none');
    % xlim([-0.1 0.1])
    % ylim([-0.1 0.1])
    % zlim([-0.1 0.1])
    % 
    % % figure
    % 
    % ax2 = axes;
    % scatter3(ax2,model_tr.X_train(trktr==2,1),model_tr.X_train(trktr==2,2),model_tr.X_train(trktr==2,3),20,postr(trktr==2),'filled');
    % hold on;
    % colormap(ax2, "autumn");          % colormap for first scatter
    % set(ax2, 'Color', 'none');
    % xlim([-0.1 0.1])
    % ylim([-0.1 0.1])
    % zlim([-0.1 0.1])

    % 
    % figure
    % scatter3(model_tr.X_train(trktr==1,1),model_tr.X_train(trktr==1,2),model_tr.X_train(trktr==1,3),10,'filled');
    % hold on
    % scatter3(model_tr.X_train(trktr==2,1),model_tr.X_train(trktr==2,2),model_tr.X_train(trktr==2,3),10,'filled');

    % % ---- Constrained inference of TEST latent trajectory ----------------
    Xte1 = infer_new_traj_lmt(Yte, model_tr, opts.nIter_infer);
    Xte = infer_new_traj_epgplvm(Yte, model_tr, opts.nIter_infer);
    % 
    z_Xte =Xte;
    z_Xtr = model_tr.X_train;


    figure
    subplot(2,2,1)
    scatter3(z_Xte(trkte==1,1),z_Xte(trkte==1,2),z_Xte(trkte==1,3),10,'blue','MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.2);
    hold on;
    scatter3(z_Xtr(trktr==1,1),z_Xtr(trktr==1,2),z_Xtr(trktr==1,3),10,'k','MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.2);


    subplot(2,2,2)
    scatter3(z_Xte(trkte==1,1),z_Xte(trkte==1,2),z_Xte(trkte==1,3),10,'blue','MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.2);
    hold on;
    scatter3(z_Xtr(trktr==2,1),z_Xtr(trktr==2,2),z_Xtr(trktr==2,3),10,'k','MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.2);


    % ---- Decode TEST via KNN on TRAIN latent bank -----------------------
   % ---- Preallocate error containers ----------------------------------
    % We'll fill these as [nK x nTrk] and [nK x 1]
    nTrk = numel(1:max(trktr));
    rmse_trk = nan(nK, nTrk);
    mae_trk  = nan(nK, nTrk);

    decoded_position  = nan(nK, size(Xte,1), nTrk);
    track_probability  = nan(nK, size(Xte,1), nTrk);


    parfor kIdx = 1:nK
        K = K_list(kIdx);

        % Multitrack KNN decode (track-conditional position)
        [pos_by_trk_local, trk_aff_local] = knn_decode( ...
            Xte, model_tr.X_train, postr, trktr, K);

        decoded_position(kIdx,:,:) = pos_by_trk_local;
        track_probability(kIdx,:,:) = trk_aff_local;

        % ----- Track-specific errors ------------------------------------
        for tIdx = 1:nTrk
            mask = (trkte == tIdx);      % time bins belonging to this track

            if sum(mask)>0
                pos_hat_trk = pos_by_trk_local(mask, tIdx);
                pos_true_trk = poste(mask);

                rmse_trk(kIdx, tIdx) = sqrt(mean((pos_hat_trk - pos_true_trk).^2, 'omitnan'));
                mae_trk(kIdx,  tIdx) = mean(abs(pos_hat_trk - pos_true_trk), 'omitnan');
            end
        end
    end

    % ----- Store everything for this fold -------------------------------
    folds(f).fold        = f;
    folds(f).model = model_tr;
    folds(f).t           = tte;
    folds(f).X           = Xte;
    folds(f).pos_true    = poste;
    folds(f).trk_true    = trkte;

    folds(f).K_list      = K_list;
    folds(f).rmse    = rmse_trk;   % [nK x nTrk]
    folds(f).mae     = mae_trk;    % [nK x nTrk]
    folds(f).decoded_position = decoded_position;
    folds(f).track_probability = track_probability;


    % Best K for this fold (by overall RMSE)
    [~, bestK_idx] = min(mean(rmse_trk','omitnan'));
    folds(f).bestK_idx = bestK_idx;
    folds(f).bestK     = K_list(bestK_idx);
end


% Collect per-fold best K indices
bestK_idx_all = [folds.bestK_idx];               % [1 x nFolds]
bestK_all     = [folds.bestK];                   % the actual K values

[bestK,F,C] = mode(bestK_all);

if F <= 5 % if less than 5 folds out of 10 prefer different K
    [~,bestK_idx] = min(abs(median(bestK_all) - K_list));
    bestK     = K_list(bestK_idx);
else % mode
    [~,bestK_idx] = min(abs(bestK - K_list));
end

folds.bestK_all = bestK;


%%%%%%%% Repeat with cell id shuffles
for f = 1:nFolds
    te = idx_all(fold_edges(f):fold_edges(f+1)-1);
    tr = setdiff(idx_all, te);

    Ytr   = Y_run(tr,:);  Yte   = Y_run(te,:);
    ttr   = t_run(tr);    tte   = t_run(te);
    postr = pos_run(tr);  poste = pos_run(te);
    trktr = trk_run(tr);  trkte = trk_run(te);

    %
    % % ---- Constrained inference of TEST latent trajectory ----------------
    nTrk = numel(1:max(trktr));
    rmse_trk_shuf = nan(1000, nTrk);
    mae_trk_shuf  = nan(1000, nTrk);
    trk_aff_shuf  = nan(1000, size(Yte,1), nTrk);

    parfor nshuffle = 1:1000
        s = RandStream('mrg32k3a','Seed',nshuffle);  % deterministic per shuffle

        cell_index  = randperm(s, size(Yte,2));
        Y_run_shuf = Yte(:, cell_index);

        Xte_shuf = infer_new_traj_lmt(Y_run_shuf, folds(f).model, opts.nIter_infer);

        % ---- Decode TEST via KNN on TRAIN latent bank -----------------------
        % ---- Preallocate error containers ----------------------------------

        K = bestK;

        % Multitrack KNN decode (track-conditional position)
        [pos_by_trk, trk_aff] = knn_decode( ...
            Xte_shuf, model_tr.X_train, postr, trktr, K);

        trk_aff_shuf(nshuffle,:,:) = trk_aff;

        % ----- Track-specific errors ------------------------------------
        for tIdx = 1:nTrk
            mask = (trkte == tIdx);      % time bins belonging to this track
            if sum(mask)>0
                pos_hat_trk = pos_by_trk(mask, tIdx);
                pos_true_trk = poste(mask);

                rmse_trk_shuf(nshuffle, tIdx) = sqrt(mean((pos_hat_trk - pos_true_trk).^2, 'omitnan'));
                mae_trk_shuf(nshuffle,  tIdx) = mean(abs(pos_hat_trk - pos_true_trk), 'omitnan');
            end
        end
    end


    % ----- Store shuffled results for this fold -------------------------------
    folds(f).rmse_shuffled    = rmse_trk_shuf;   % [nshuffle x nTrk]
    folds(f).mae_shuffled     = mae_trk_shuf;    % [nshuffle x nTrk]
    folds(f).track_probability_shuffled     = trk_aff_shuf;    % [nshuffle x nTrk]

end

% 
% %%%% ROC plot for track discriminiabiity 
% labels=decoded_Bayesian_RUN.track_id(good_idx);
% log_odds=decoded_Bayesian_RUN.position_log_odds(good_idx);
% 
% fpr = 0:0.01:1;
% TPR_real = zeros(1000, length(fpr));
% AUC_real = zeros(1,1000);
% parfor i = 1:1000
%     s = RandStream('mrg32k3a','Seed',i); % Set random seed for resampling
%     idx = randsample(s,length(labels), length(labels), true);
%     %         [~, ~, ~, AUC_real(i)] = perfcurve(labels(idx), log_odds(idx),1,'XVals', fpr);
%     [x,y,~,AUC_real(i)] = perfcurve(labels(idx), log_odds(idx),1, 'XVals', fpr);
%     [C,ia,ic] = unique(x);
%     y = y(ia);
%     x = x(ia);
% 
%     TPR_interp = interp1(x, y, fpr, 'linear', 'extrap');  % force to match fpr
%     TPR_real(i,:) = TPR_interp;
% end
% Bayesian_log_odds_RUN_validation.AUC_real_good_bins = AUC_real;
% 
% TPR_shuf = zeros(1000, length(fpr));
% AUC_shuf = zeros(1,1000);
% for i = 1:1000
%     log_odds_shuffled = decoded_Bayesian_RUN_shuffled{i}.position_log_odds(good_idx);
%     s = RandStream('mrg32k3a','Seed',i); % Set random seed for resampling
%     idx = randsample(s,length(labels), length(labels), true);
%     [x, y, ~, AUC_shuf(i)] = perfcurve(labels, log_odds_shuffled, 1, 'XVals', fpr);
% 
%     %     [x,y,~,AUC_shuf(i)] = perfcurve(labels(idx), log_odds(:),1, 'XVals', fpr);
%     [C,ia,ic] = unique(x);
%     y = y(ia);
%     x = x(ia);
% 
%     TPR_interp = interp1(x, y, fpr, 'linear', 'extrap');  % force to match fpr
%     TPR_shuf(i,:) = TPR_interp;
% end
% Bayesian_log_odds_RUN_validation.AUC_cell_ID_shuffle_good_bins = AUC_shuf;


%% ---- Final model trained on ALL RUN (for decoding new data) -------------
model_all = train_lmt_pgplvm(Y_run, opts.latentDim, opts.nIter_train);

model = struct('lmt',          model_all, ...
               'X_bank',       model_all.X_train, ...
               'pos_bank',     pos_run(:), ...
               'trk_bank',     trk_run(:), ...
               'unit_ids',     unit_ids, ...
               'binSize_train',opts.binSize_train, ...
               't_run',        t_run);
model.bestK =bestK;




end % === main function end ===

% 
% %% ========================================================================
% %  LMT-style trainer (simplified reference implementation)
% %  - Uses prior_kernel_sp, covariance_fun, StateSpaceModelsofSpikeTrains_ref,
% %    logmargli_gplvm_se_block_ref_nogrid from the LMT repo.
% %  - Single continuous segment (no multi-trial handling).
% % ========================================================================
% function model = train_lmt_pgplvm_run(Y, latentDim, nIter)
% 
% % Y: [T x N] spike counts
% % latentDim: number of latent dimensions
% % nIter: number of outer iterations (alternating tuning + latent)
% 
% [nt0, nneur] = size(Y);
% nf = latentDim;
% 
% % Single segment: ntr = 1, nt = T
% ntr = 1;
% nt  = nt0;
% 
% yy = reshape(Y, [ntr, nt, nneur]);   % [1 x nt x N]
% 
% % ---- initial latent from PLDS (LMT PLDS folder) ----
% xinit = run_plds(Y, nf)';            % [T x nf]  (PLDS returns nf x T)
% if all(abs(xinit(:)) < eps)
%     xinit = randn(size(xinit))*1e-4;
% end
% xinit = reshape(xinit, [ntr, nt, nf]);  % [1 x nt x nf]
% 
% % ---- hyperparameters / options (LMT-like defaults) ----
% latentTYPE  = 1;      % 1: AR1/OU time kernel
% ffTYPE      = 2;      % 2: SE tuning kernel
% rhoxx       = 100;    % latent variance
% lenxx       = 100;    % latent length-scale (in bins)
% rhoff       = 1;      % tuning variance
% lenff_ratio = 0.1;    % tuning lengthscale = median(range(x))/lenff_ratio
% sigma2      = 3;      % small noise on tuning kernel (LMT anneals this; we keep fixed)
% 
% % ---- build time prior operators in "u" space (LMT utilities) ------------
% tgrid = (1:nt)';      % discrete time index
% 
% [Bbase, BTbase, nu] = prior_kernel_sp(rhoxx, lenxx, nt, latentTYPE, tgrid);
% % Wrap B and BT to act on [ntr x nt x nf] arrays (same as demo2_*_ref)
% Bfun  = @(x,flag) permute( ...
%             reshape(Bbase(reshape(permute(x,[2,1,3]), size(x,2), []), flag), ...
%                     [], size(x,1), size(x,3)), ...
%             [2,1,3]);
% BTfun = @(x,flag) permute( ...
%             reshape(BTbase(reshape(permute(x,[2,1,3]), size(x,2), []), flag), ...
%                     [], size(x,1), size(x,3)), ...
%             [2,1,3]);
% 
% % ---- initialize u and x from xinit --------------------------------------
% uu     = Bfun(xinit, 1);    % [1 x nu x nf] "whitened" latent
% xxsamp = Bfun(uu, 0);       % [1 x nt x nf] latent trajectory
% xxsamp_mt = reshape(xxsamp, [], nf);   % [T x nf]
% 
% % ---- grid over latent space (inducing / tuning grid) --------------------
% ng     = choose_ng(nf);
% ranges = [min(xxsamp_mt); max(xxsamp_mt)]';   % [nf x 2]
% xgrid  = gen_grid(ranges, ng, nf);            % [ng^nf x nf] (LMT utility)
% 
% % ---- tuning GP hyperparams and initial fgrid ----------------------------
% lenff  = median(range(xxsamp_mt,1)) / lenff_ratio;
% covfun = covariance_fun(rhoff, lenff, ffTYPE);    % LMT tuning kernel
% cxx    = covfun(xgrid, xgrid);
% invcxx = pdinv(cxx + cxx(1,1)*sigma2*eye(size(cxx)));
% 
% % rough initial fgrid using GP regression from averaged latent & rates
% xavg = squeeze(mean(xxsamp,1));  % [nt x nf]
% if nf == 1
%     xavg = xavg(:);
% end
% K_xx  = covfun(xavg, xavg);
% K_xg  = covfun(xavg, xgrid);
% mean_rate = squeeze(mean(yy,1));  % [nt x N]
% fgrid = (K_xg' / (K_xx + 1e-6*eye(size(K_xx)))) * mean_rate;  % [ng x N]
% fgrid = fgrid / max(abs(fgrid(:)) + eps);  % scale
% 
% % ---- optimization settings (minFunc) ------------------------------------
% opt_f             = [];
% opt_f.Method      = 'scg';
% opt_f.TolFun      = 1e-4;
% opt_f.MaxIter     = 1e1;
% opt_f.maxFunEvals = 1e1;
% opt_f.Display     = 'off';
% 
% opt_u = opt_f;
% 
% % matrix view of spikes: [T x N]
% yymat = reshape(permute(yy, [2,1,3]), [], nneur);  % [nt x N]
% 
% % ---- alternating updates: tuning, then latent ---------------------------
% for it = 1:nIter
%     % --- (1) optimize tuning functions fgrid given current latent --------
%     xxsamp_mt = reshape(xxsamp, [], nf);         % [T x nf]
%     ctx       = covfun(xxsamp_mt, xgrid);        % [T x Ng]
%     kk        = ctx * invcxx;                    % [T x Ng]
% 
%     obj_f = @(ff) StateSpaceModelsofSpikeTrains_ref(ff, yymat, invcxx, kk);
%     ff0   = fgrid(:);
%     ffnew = minFunc(obj_f, ff0, opt_f);
%     fgrid = reshape(ffnew, size(xgrid,1), nneur);
% 
%     % --- (2) optimize latent uu given tuning (new trajectory) ------------
%     invkf = invcxx * fgrid;                      % [Ng x N]
% 
%     u0 = uu(:);
%     obj_u = @(uvec) logmargli_gplvm_se_block_ref_nogrid( ...
%                         uvec, xgrid, invkf, yymat, ...
%                         Bfun, covfun, nf, BTfun, ntr);
% 
%     u_new = minFunc(obj_u, u0, opt_u);
%     uu    = reshape(u_new, [ntr, nu, nf]);
% 
%     % update latent trajectory
%     xxsamp = Bfun(uu, 0);
% 
%     % (optional) convergence check
%     if it > 1
%         delta = max(abs(u_new - u0));
%         if delta < 1e-4
%             break;
%         end
%     end
% end
% 
% % flatten latent into [T x nf]
% X_train = reshape(xxsamp, [], nf);
% 
% % pack model
% model = struct();
% model.X_train    = X_train;
% model.rhoxx      = rhoxx;
% model.lenxx      = lenxx;
% model.rhoff      = rhoff;
% model.lenff      = lenff;
% model.sigma2     = sigma2;
% model.latentTYPE = latentTYPE;
% model.ffTYPE     = ffTYPE;
% model.nneur      = nneur;
% model.xgrid      = xgrid;        % [Ng x nf]
% model.fgrid      = fgrid;        % [Ng x N] (typically log-rate)
% model.invcxx     = invcxx;       % [Ng x Ng] inverse grid kernel
% 
% end


%% ========================================================================
%  LMT-style trainer using pgplvm_la (from LMT repo)
%  - Single continuous segment (no multi-trial handling).
%  - Uses PLDS init (run_plds)
% ========================================================================
function model = train_lmt_pgplvm(Y, latentDim, nIter)
% Train P-GPLVM using pgplvm_la and build a latent grid + tuning functions
% for later use (tuning curves + consistent test-time inference).
%
% Inputs:
%   Y         [T x N] spike counts
%   latentDim scalar, number of latent dimensions (nf)
%   nIter     number of outer iterations for pgplvm_la (setopt.niter)
%
% Output:
%   model struct with fields:
%     X_train    [T x nf]   latent trajectory (MAP)
%     ff_time    [T x N]   log-rate at each time bin (from pgplvm_la)
%     xgrid      [Ng x nf] latent grid
%     fgrid      [Ng x N] tuning values on grid (log-rate)
%     invcxx     [Ng x Ng] inverse kernel on grid
%     rhoxx,lenxx,rhoff,lenff,latentTYPE,ffTYPE,sigma2,nneur

[T, nneur] = size(Y);
nf = latentDim;

%% --- 1. Initial latent from PLDS (same idea as before) --------------
xplds = run_plds(Y, nf)';   % [T x nf]
if all(abs(xplds(:)) < eps)
    xplds = randn(size(xplds)) * 1e-4;
end

%% --- 2. Set options for pgplvm_la -----------------------------------
setopt = struct();
setopt.sepx_flag    = 0;
setopt.sigma2_init  = 3; % small noise term on the tuning-kernel (K_ff) to keep inverses stable
setopt.sigma2_end   = min([0.1, setopt.sigma2_init]);
setopt.lr           = 0.95; % Increase lr (e.g., 0.97–0.99) so sigma2 drops faster,
setopt.latentTYPE   = 1;  % 1: AR1/OU
setopt.ffTYPE       = 2;  % 2: SE
setopt.initTYPE     = 1;  % use PLDS init
setopt.la_flag      = 3;  % 2 standard Laplace. 3 Decoupled
setopt.rhoxx        = 100;
setopt.lenxx        = 100;
setopt.rhoff        = 1;
setopt.lenff        = median(range(xplds,1));  % rough init
setopt.lenff_ratio  = 1;
setopt.b            = 0;
setopt.r            = 1;
setopt.nsevar       = 1;
% hyperparameters to optimize (like demo2)
setopt.hypid        = [2,3,4];  % lenxx, rhoff, lenff
setopt.xinit        = xplds;    % initialization
setopt.xpldsmat = xplds;   % 
setopt.xplds    = xplds;   % 
setopt.niter        = nIter;
setopt.opthyp_flag  = 0;        % (set to 1 if you want hyper-opt)

% pgplvm_la expects xplds (for plotting/diagnostics); we just pass it
result = pgplvm_la(Y, xplds, [], setopt);

% result typically has:
%   result.xxsamp : [T x nf] MAP latent
%   result.ffmat  : [T x N] log-rates per time bin
X_train = result.xxsamp;
ff_time = result.ffmat;

%% --- 3. Build latent grid and tuning GP on grid ---------------------
latentTYPE = setopt.latentTYPE;
ffTYPE     = setopt.ffTYPE;
rhoxx      = result.rhoxx;
lenxx      = result.lenxx;
rhoff      = result.rhoff;
lenff      = result.lenff;

% Small noise term for tuning GP on grid (not the same as Poisson noise,
% just a numerical regularizer in kernel inversion)
sigma2 = 0.1;

% Build latent grid (like demo2, but slightly cleaned)
nf = size(X_train, 2);
ng = choose_ng(nf);   % you already have this helper

ranges = [min(X_train); max(X_train)]';   % [nf x 2]
xgrid  = gen_grid(ranges, ng, nf);        % [Ng x nf]

% GP over tuning as a function of latent
covfun = covariance_fun(rhoff, lenff, ffTYPE);
cxx    = covfun(xgrid, xgrid);                            % [Ng x Ng]
invcxx = pdinv(cxx + sigma2 * cxx(1,1) * eye(size(cxx))); % [Ng x Ng]

% Smooth time-varying log-rates ff_time over latent using GP regression.
% We treat (X_train, ff_time) as training data; xgrid as query.
K_xx   = covfun(X_train, X_train);                        % [T x T]
K_gx   = covfun(xgrid,  X_train);                         % [Ng x T]
% add tiny jitter for numerical stability
A      = K_xx + 1e-6 * eye(T);
% fgrid: [Ng x N], log-rate at each latent grid point
fgrid  = K_gx / A * ff_time;

%% --- 4. Pack model ---------------------------------------------------
model = struct();
model.X_train    = X_train;      % [T x nf]
model.ff_time    = ff_time;      % [T x N]
model.xgrid      = xgrid;        % [Ng x nf]
model.fgrid      = fgrid;        % [Ng x N]
model.invcxx     = invcxx;       % [Ng x Ng]

model.rhoxx      = rhoxx;
model.lenxx      = lenxx;
model.rhoff      = rhoff;
model.lenff      = lenff;
model.sigma2     = sigma2;       % GP noise for grid inversion
model.latentTYPE = latentTYPE;
model.ffTYPE     = ffTYPE;
model.nneur      = nneur;
end

%% ========================================================================
%  Constrained inference of a NEW latent trajectory given a trained model
%  - Fixes hyperparameters from training.
% - Decoupled 
% - Based on Luo et al 2024
% ========================================================================

function Xnew = infer_new_traj_epgplvm(Ynew, model, nIter)
% Infer a new latent trajectory using the Extended P-GPLVM test-time
% procedure (Luo et al. 2024) as closely as possible within the LMT code.
%
% Two stages:
%   (1) Per-bin Bayesian initialization on the latent grid xgrid using
%       the learned tuning curves fgrid.
%   (2) Constrained refinement using training tuning curves (TCs)
%       (X_train, ff_time) as the GP support set, with fixed
%       hyperparameters and temporal smoothness prior.
%
% Inputs:
%   Ynew   [T_new x N] spike counts for new segment
%   model  struct from train_lmt_pgplvm (must have:
%           X_train, ff_time, xgrid, fgrid, invcxx,
%           rhoxx, lenxx, rhoff, lenff, ffTYPE, latentTYPE,
%           nneur, sigma2)
%   nIter  number of outer iterations for refinement (e.g. 10–15)
%
% Output:
%   Xnew   [T_new x nf] inferred latent trajectory

    % -------------------------------------------------------------
    % Basic size checks
    % -------------------------------------------------------------
    [T_new, N] = size(Ynew);
    nf         = size(model.X_train, 2);

    if N ~= model.nneur
        error('infer_new_traj_epgplvm: neuron count mismatch (Ynew has %d, model has %d).', ...
              N, model.nneur);
    end

    % -------------------------------------------------------------
    % Stage 1: Bayesian per-bin init on latent grid (Eqs. 9–12)
    % -------------------------------------------------------------
    % Use the training latent grid + tuning curves as the discrete
    % support {z_k}, and compute:
    %   log p(z_k | y_t) ∝ log p(y_t | z_k)
    %                   = sum_n [ y_{t,n} log λ_{k,n} - λ_{k,n} ]
    % (factorial terms are constant in z_k and are dropped)

    xgrid  = model.xgrid;   % [Ng x nf]
    fgrid  = model.fgrid;   % [Ng x N] (log firing rates per bin)
    Ng     = size(xgrid, 1);

    % Convert to mean counts per bin
    lambda_grid      = exp(fgrid);                % [Ng x N]
    log_lambda_grid  = log(lambda_grid + 1e-12);  % [Ng x N]
    sum_lambda_grid  = sum(lambda_grid, 2);       % [Ng x 1]

    xinit = zeros(T_new, nf);

    for t = 1:T_new
        y = double(Ynew(t, :));  % [1 x N]

        % log-likelihood over grid points (up to const in y)
        % ll_k = sum_n y_n * log λ_{k,n} - λ_{k,n}
        ll_mat  = bsxfun(@times, log_lambda_grid, y);  % [Ng x N]
        loglik  = sum(ll_mat, 2) - sum_lambda_grid;    % [Ng x 1]

        % MAP grid point
        [~, k_star]   = max(loglik);
        xinit(t, :)   = xgrid(k_star, :);
    end

    % -------------------------------------------------------------
    % Stage 2: temporal GP prior + TC-based mapping (Eqs. 13–15)
    % -------------------------------------------------------------

    % Time prior hyperparameters (fixed from training)
    latentTYPE = model.latentTYPE;
    rhoxx      = model.rhoxx;
    lenxx      = model.lenxx;

    % Tuning GP hyperparameters (fixed from training)
    rhoff      = model.rhoff;
    lenff      = model.lenff;
    ffTYPE     = model.ffTYPE;
    sigma2     = model.sigma2;   % small noise on tuning kernel

    % Wrap Ynew into [ntr x nt x N]
    ntr = 1;
    nt  = T_new;
    yy  = reshape(Ynew, [ntr, nt, N]);  % [1 x T_new x N]

    % Build time prior operators for this segment
    tgrid = (1:nt)';

    [Bbase, BTbase, nu] = prior_kernel_sp(rhoxx, lenxx, nt, latentTYPE, tgrid);

    % Wrap B and BT to act on [ntr x nt x nf] arrays (LMT convention)
    Bfun  = @(x,flag) permute( ...
                reshape(Bbase(reshape(permute(x,[2,1,3]), size(x,2), []), flag), ...
                        [], size(x,1), size(x,3)), ...
                [2,1,3]);

    BTfun = @(x,flag) permute( ...
                reshape(BTbase(reshape(permute(x,[2,1,3]), size(x,2), []), flag), ...
                        [], size(x,1), size(x,3)), ...
                [2,1,3]);

    % Initialize uu (whitened latent) from xinit
    xinit_3d = reshape(xinit, [ntr, nt, nf]);  % [1 x T_new x nf]
    uu       = Bfun(xinit_3d, 1);             % [1 x nu x nf]
    xxsamp   = Bfun(uu, 0);                   % [1 x nt x nf]

    % -------------------------------------------------------------
    % Tuning curves from training as TC support set
    % -------------------------------------------------------------
    % In the paper, the training trajectories and their log tuning
    % curves (TCs) form the support set {X_TC, f_TC} for the GP
    % mapping; the mapping for new data is constrained by this set.
    %
    % Here we use:
    %   X_TC  = model.X_train   [T_train x nf]
    %   f_TC  = model.ff_time   [T_train x N]
    % and precompute invK_TC * f_TC.

    X_TC  = model.X_train;   % [T_train x nf]
    f_TC  = model.ff_time;   % [T_train x N]

    covfun = covariance_fun(rhoff, lenff, ffTYPE);

    K_TC   = covfun(X_TC, X_TC);   % [T_train x T_train]
    invK_TC = pdinv(K_TC + sigma2 * K_TC(1,1) * eye(size(K_TC)));  % [T_train x T_train]

    invkf_TC = invK_TC * f_TC;     % [T_train x N]

    % Matrix view of spikes for objective
    yymat = reshape(permute(yy, [2,1,3]), [], N);  % [T_new x N]

    % -------------------------------------------------------------
    % Optimization over uu (whitened latent), TCs fixed
    % -------------------------------------------------------------
    % Objective is the LMT log-marginal-likelihood with:
    %   - Bfun/BTfun from time prior
    %   - covfun from tuning GP
    %   - support set {X_TC, invkf_TC}
    %
    % This corresponds to Eqs. 13–15: smoothness parameters frozen,
    % mapping constrained by TCs, and we only refine X_new.

    opt_u             = [];
    opt_u.Method      = 'scg';
    opt_u.TolFun      = 1e-4;
    opt_u.MaxIter     = 1e1;
    opt_u.maxFunEvals = 1e1;
    opt_u.Display     = 'off';

    for it = 1:nIter
        u0 = uu(:);

        obj_u = @(uvec) logmargli_gplvm_se_block_ref_nogrid( ...
                            uvec, X_TC, invkf_TC, yymat, ...
                            Bfun, covfun, nf, BTfun, ntr);

        u_new = minFunc(obj_u, u0, opt_u);
        uu    = reshape(u_new, [ntr, nu, nf]);
        xxsamp = Bfun(uu, 0);

        % Simple convergence check
        if it > 1
            delta = max(abs(u_new - u0));
            if delta < 1e-4
                break;
            end
        end
    end

    % Flatten latent trajectory into [T_new x nf]
    Xnew = reshape(xxsamp, [], nf);
end





%% ========================================================================
%  Constrained inference of a NEW latent trajectory given a trained model
%  - Fixes hyperparameters from training and reuses LMT objectives.
% - Decoupled 
% ========================================================================
function Xnew = infer_new_traj_lmt(Ynew, model, nIter)
% Infer a new latent trajectory given a trained P-GPLVM model
% (from train_lmt_pgplvm).
%
% Inputs:
%   Ynew   [T_new x N]  spike counts for new segment
%   model  struct from training (must contain xgrid, fgrid, invcxx, hyperparams)
%   nIter  number of outer iterations (usually small, e.g. 5–10)
%
% Output:
%   Xnew   [T_new x nf] inferred latent trajectory for Ynew

    [T_new, N] = size(Ynew);
    nf         = size(model.X_train, 2);

    if N ~= model.nneur
        error('infer_new_traj_lmt_pgplvm: neuron count mismatch (Ynew has %d, model has %d).', ...
              N, model.nneur);
    end

    %% --- 1. PLDS init on new data --------------------------------------
    xinit = run_plds(Ynew, nf)';   % [T_new x nf]
    if all(abs(xinit(:)) < eps)
        xinit = randn(size(xinit)) * 1e-4;
    end

    xinit_mt = reshape(xinit, [], nf);

    % training stats
    mu_tr  = mean(model.X_train, 1);
    std_tr = std(model.X_train, [], 1) + 1e-6;

    % init stats
    mu_in  = mean(xinit_mt, 1);
    std_in = std(xinit_mt, [], 1) + 1e-6;

    % rescale PLDS init to training stats
    xinit_mt = bsxfun(@minus,  xinit_mt, mu_in);
    xinit_mt = bsxfun(@times,  xinit_mt, std_tr ./ std_in);
    xinit_mt = bsxfun(@plus,   xinit_mt, mu_tr);

    ntr = 1;
    nt  = T_new;
    yy  = reshape(Ynew, [ntr, nt, N]);     % [1 x T_new x N]

    xinit = reshape(xinit_mt, [ntr, nt, nf]);
    % xinit = reshape(xinit, [ntr, nt, nf]); % [1 x T_new x nf]

    %% --- 2. Time prior operators for new time grid ---------------------
    latentTYPE = model.latentTYPE;
    rhoxx      = model.rhoxx;
    lenxx      = model.lenxx;
    rhoff      = model.rhoff;
    lenff      = model.lenff;

    tgrid = (1:nt)';

    [Bbase, BTbase, nu] = prior_kernel_sp(rhoxx, lenxx, nt, latentTYPE, tgrid);

    % Wrap B and BT to act on [ntr x nt x nf] arrays (same pattern as training)
    Bfun  = @(x,flag) permute( ...
                reshape(Bbase(reshape(permute(x,[2,1,3]), size(x,2), []), flag), ...
                        [], size(x,1), size(x,3)), ...
                [2,1,3]);

    BTfun = @(x,flag) permute( ...
                reshape(BTbase(reshape(permute(x,[2,1,3]), size(x,2), []), flag), ...
                        [], size(x,1), size(x,3)), ...
                [2,1,3]);

    % init u and x for new segment
    uu     = Bfun(xinit, 1);   % [1 x nu x nf] whitened latent
    xxsamp = Bfun(uu, 0);      % [1 x nt x nf] latent trajectory

    %% --- 3. Fixed tuning GP from training ------------------------------
    covfun = covariance_fun(rhoff, lenff, model.ffTYPE);

    xgrid  = model.xgrid;       % [Ng x nf]
    fgrid  = model.fgrid;       % [Ng x N]
    invcxx = model.invcxx;      % [Ng x Ng]

    % Precompute invKF = K_grid^{-1} * fgrid (Ng x N)
    invkf = invcxx * fgrid;

    % Matrix view of spikes: [T_new x N]
    yymat = reshape(permute(yy, [2,1,3]), [], N);

    %% --- 4. Optimize latent trajectory with fixed tuning ---------------
    % Options for minFunc (same as training style)
    opt_u             = [];
    opt_u.Method      = 'scg';
    opt_u.TolFun      = 1e-4;
    opt_u.MaxIter     = 1e1;
    opt_u.maxFunEvals = 1e1;
    opt_u.Display     = 'off';

    for it = 1:nIter
        u0 = uu(:);

        % Objective: log marginal likelihood given fixed tuning (invkf, xgrid)
        obj_u = @(uvec) logmargli_gplvm_se_block_ref_nogrid( ...
                            uvec, xgrid, invkf, yymat, ...
                            Bfun, covfun, nf, BTfun, ntr);

        u_new = minFunc(obj_u, u0, opt_u);
        uu    = reshape(u_new, [ntr, nu, nf]);
        xxsamp = Bfun(uu, 0);

        % Simple convergence check
        if it > 1
            delta = max(abs(u_new - u0));
            if delta < 1e-4
                break;
            end
        end
    end

    % flatten latent into [T_new x nf]
    Xnew = reshape(xxsamp, [], nf);
end



%% ========================================================================
%  KNN decoder in latent space
% ========================================================================
function [pos_by_trk, trk_aff]  = knn_decode(Xq, Xbank, pos_bank, trk_bank, K)

% Inputs:
%   Xq        [Tq x m]   query latent trajectory
%   Xbank     [Tb x m]   training latent trajectory
%   pos_bank  [Tb x 1]   positions for training bins
%   trk_bank  [Tb x 1]   track IDs for training bins (can be 1,2,3,... or any integers)
%   K         scalar     number of nearest neighbours
%
% Outputs:
%   pos_by_trk [Tq x nTrk]  decoded position for EACH track hypothesis
%                           (column i corresponds to track_ids(i))
%   trk_aff    [Tq x nTrk]  affinity / soft evidence for each track
%                           (fraction of neighbours from that track)
%   track_ids  [1 x nTrk]   unique track labels (sorted)
%
% Example use:
%   [pos_by_trk, trk_aff, track_ids] = knn_decode_multitrack(Xte, Xtrain, pos_train, trk_train, 25);
%   % Hard decision per bin:
%   [~, best_idx] = max(trk_aff, [], 2);
%   pos_hat = pos_by_trk(sub2ind(size(pos_by_trk), (1:size(pos_by_trk,1))', best_idx));

    % ---- KNN in latent space -------------------------------------------
    [idx, D] = knnsearch(Xbank, Xq, 'K', K);   % idx: [Tq x K], D: [Tq x K]

    % distance weights
    w = 1 ./ (D + 1e-9);       % [Tq x K]
    w = w ./ sum(w, 2);        % rows sum to 1

    % neighbours' positions and tracks
    pos_nn = pos_bank(idx);    % [Tq x K]
    nn     = trk_bank(idx);    % [Tq x K]

    % ---- unique tracks -------------------------------------------------
    track_ids = unique(trk_bank(:))';   % row vector
    nTrk = numel(track_ids);
    Tq   = size(Xq, 1);

    pos_by_trk = nan(Tq, nTrk);
    trk_aff    = nan(Tq, nTrk);

    % ---- loop over tracks ----------------------------------------------
    for i = 1:nTrk
        trk = track_ids(i);

        % mask of neighbours from this track
        mask = (nn == trk);                    % [Tq x K] logical

        % affinity: fraction of neighbours from this track
        trk_aff(:, i) = mean(mask, 2);         % [Tq x 1]

        % weights restricted to this track
        w_trk = w;
        w_trk(~mask) = 0;

        denom = sum(w_trk, 2);                 % total weight from this track

        % distance-weighted position using ONLY this track's neighbours
        pos_i = sum(w_trk .* pos_nn, 2) ./ max(denom, 1e-9);

        % if denom == 0 (no neighbours of this track), this will give NaN
        pos_i(denom == 0) = NaN;

        pos_by_trk(:, i) = pos_i;
    end

end


%% ========================================================================
%  Binning & helper utilities (unchanged from your original code)
% ========================================================================
function [Y, unit_ids] = bin_counts(st, sid, edges)
unit_ids = unique(sid(:)'); 
nU       = numel(unit_ids); 
nB       = numel(edges)-1;
Y        = zeros(nB,nU,'double');
b        = discretize(st, edges);
for u = 1:nU
    bu = b(sid==unit_ids(u));
    bu = bu(isfinite(bu) & bu>=1 & bu<=nB);
    if ~isempty(bu)
        Y(:,u) = accumarray(bu,1,[nB 1]);
    end
end
end

function y = accum_by_bin(binIdx, v, fun, logicalOut)
if nargin<4, logicalOut=false; end
K = max(binIdx(~isnan(binIdx)));
y = nan(K,1);
for k=1:K
    ii = (binIdx==k) & ~isnan(v);
    if any(ii)
        y(k) = fun(v(ii));
    end
end
if logicalOut, y = logical(y); end
end

function m = mode_by_bin(binIdx, v)
K = max(binIdx(~isnan(binIdx))); 
m = nan(K,1);
for k=1:K
    ii = (binIdx==k) & ~isnan(v);
    if any(ii), m(k)=mode(v(ii)); end
end
end

function [a,b,c,d]=makecol(a,b,c,d), a=a(:); b=b(:); c=c(:); d=d(:); end

function v = fget(s,f,d), v=d; if isfield(s,f)&&~isempty(s.(f)), v=s.(f); end, end


%% ========================================================================
%  Small helpers for the LMT bits
% ========================================================================
function ng = choose_ng(nf)
% Roughly match LMT demo grid sizes
switch nf
    case 1, ng = 25;
    case 2, ng = 10;
    case 3, ng = 10;
    case 4, ng = 5;
    case 5, ng = 5;
    case 6, ng = 4;
    case 7, ng = 3;
    otherwise, ng = 2;
end
end
