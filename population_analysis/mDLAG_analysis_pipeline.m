function [estParams,trackedParams,flags] = mDLAG_analysis_pipeline(seqTrue,binWidth,yDims,options)

% save_all_figures('C:\Users\masahiro.takigawa\Documents\GitHub\mDLAG\mDLAG\demo\results',[])
%% =====================================
% 0) spike counts
% ======================================
if ~isfield(options,'load_fitted_model')
    options.load_fitted_model=[];
end

if isempty(options.load_fitted_model)
%% =====================================
% 1a) Initialize mDLAG model parameters
% ======================================

xDim_fit = 10;           % Number of latent variables to be fitted
randomSeed = 0;          % Seed the random number generator for reproducibility
prior_val = 1e-12;       % Set to very small value to set uninformative prior
prior.d.beta = prior_val;
prior.phi.a = prior_val;
prior.phi.b = prior_val;
prior.alpha.a = prior_val;
prior.alpha.b = prior_val;
saveCcov = false;        % Set to false to save memory when saving initParams
initParams = init_mdlag(seqTrue, ...
                        yDims, ...
                        xDim_fit, ...
                        binWidth, ...
                        'prior', prior, ...
                        'randomSeed', randomSeed, ...
                        'saveCcov', saveCcov);

%% =======================
% 1b) Set fitting parameters
% ===========================

% Let's explicitly define relevant optional arguments, for 
% the sake of demonstration:
tol = 1e-8;           % Tolerance to determine fitting convergence
maxIters = 5000;      % Maximum fitting iterations
freqLB = 10;          % Check for convergence of lower bound every freqLB iterations
freqParam = 100;      % Store intermediate delay and timescale estimates every freqParam iterations
learnDelays = true;   % Toggle whether to learn delay parameters
verbose = true;       % Print fitting progress
minVarFrac = 0.001;   % Private noise variances will not be allowed to go below this value
maxDelayFrac  = 0.5;  % Maximum delay magnitude (unit: fraction of trial length)
maxTauFrac = 1.0;     % Maximum timescale magnitude (unit: fraction of trial length)
pruneX = true;        % For speed-up, remove latents that become inactive in all groups
saveXcov = false;     % Set to false to save memory when saving final results
saveCcov = false;     % Set to false to save memory when saving final results
segLength = Inf;      % Optional speedup to cut trials into smaller segments during fitting

%% ========================================================================
% 1c) Fit mDLAG model
%     NOTE: You can skip to Section 2a if a trained model already exists.
% =========================================================================

% Optional speedup to cut trials into smaller segments during fitting
seqTrueCut = cutTrials(seqTrue, 'segLength', segLength);

[estParams,~,trackedParams,flags] ...
    = em_mdlag(initParams, ...
               seqTrueCut, ...
               xDim_fit, ...
               'prior', prior, ...
               'tol', tol, ...
               'maxIters', maxIters, ...
               'freqLB', freqLB, ...
               'freqParam', freqParam, ...
               'learnDelays', learnDelays, ...
               'verbose', verbose, ...
               'minVarFrac', minVarFrac, ...
               'maxDelayFrac', maxDelayFrac, ...
               'maxTauFrac', maxTauFrac, ...
               'pruneX', pruneX, ...
               'saveXcov', saveXcov, ...
               'saveCcov', saveCcov);
           

%% ========================
% 1d) Save fitting results
% =========================
           
% save('demo/results/demo_mdlag_results.mat', ...
%      'estParams', 'trackedParams', 'flags');
if exist(fullfile(options.ANALYSIS_DATAPATH,'mDLAG_result'))==0
    mkdir(fullfile(options.ANALYSIS_DATAPATH,'mDLAG_result'))
end
DIR = dir(fullfile(options.ANALYSIS_DATAPATH,'mDLAG_result','*mDLAG*.mat'));
all_files= [];
run_num = [];
if ~isempty(DIR)
    for  nfile = 1:length(DIR)

        all_files{nfile} = DIR(nfile).name;

        run_num(nfile) =  str2num(cell2mat(extractBetween( DIR(nfile).name,'run','.mat')));
    end
else
    all_files= '';
end

% fullfile(options.ANALYSIS_DATAPATH,'mDLAG_result'
if sum(contains(all_files,'mDLAG_run1.mat'))==0
    save(fullfile(options.ANALYSIS_DATAPATH,'mDLAG_result','mDLAG_run1.mat'),'estParams','trackedParams','flags');
else
    save(fullfile(options.ANALYSIS_DATAPATH,'mDLAG_result',sprintf('mDLAG_run%i.mat',max(run_num)+1)),'estParams','trackedParams','flags');
end

else
    % Let's explicitly define relevant optional arguments, for
    % the sake of demonstration:
    tol = 1e-8;           % Tolerance to determine fitting convergence
    maxIters = 5000;      % Maximum fitting iterations
    freqLB = 10;          % Check for convergence of lower bound every freqLB iterations
    freqParam = 100;      % Store intermediate delay and timescale estimates every freqParam iterations
    learnDelays = true;   % Toggle whether to learn delay parameters
    verbose = true;       % Print fitting progress
    minVarFrac = 0.001;   % Private noise variances will not be allowed to go below this value
    maxDelayFrac  = 0.5;  % Maximum delay magnitude (unit: fraction of trial length)
    maxTauFrac = 1.0;     % Maximum timescale magnitude (unit: fraction of trial length)
    pruneX = true;        % For speed-up, remove latents that become inactive in all groups
    saveXcov = false;     % Set to false to save memory when saving final results
    saveCcov = false;     % Set to false to save memory when saving final results
    segLength = Inf;      % Optional speedup to cut trials into smaller segments during fitting

    run_id = options.load_fitted_model; % 1 means run1
    load(fullfile(options.ANALYSIS_DATAPATH,'mDLAG_result',sprintf('mDLAG_run%i.mat',run_id)));
end
%% ========================
% 2b) Check fitting results
% =========================

% Display flags indicating fitting procedure status
% flags.Convergence -- Indicates that fitting converged according to 'tol'.
%                      Else 'maxIters' was reached before convergence.
% flags.DecreasingLowerBound -- Indicates that the lower bound (objective
%                               function) decreased at some point during
%                               fitting, which suggests there's something
%                               wrong with the fitting or data.
% flags.PrivateVarianceFloor -- Indicates that the variance floor was used
%                               on one or more observed dimensions. 
%                               See em_mdlag.m header for more info.
% flags.xDimsRemoved         -- Number of latent dimensions removed (if 
%                               pruneX is true) due to low variance in all 
%                               groups.
fprintf('\n');
disp(flags)
units = 'sec';%units of time of binWidth (for labels)
plotFittingProgress(trackedParams, binWidth, ...
                    'freqLB', freqLB, ...
                    'freqParam', freqParam, ...
                    'units', units);
%% ==============================================================
% 3) Visualize recovery of the loading matrix and ARD parameters
% ===============================================================

% Loadings matrices

% Ground truth
Ctrue = vertcat(paramsTrue.Cs{:});
hinton(Ctrue);

% Estimate
% NOTE: In general, the columns of Cest are unordered, and will not match
%       the order in Ctrue. Here, we can reorder and flip the sign of each
%       dimension to facilitate comparison with the ground truth.
reorder = [ 6  4  1  2  3  7  5];
rescale = [-1 -1 -1  1  1  1  1];
Cest = vertcat(estParams.C.means{:});
hinton(Cest(:,reorder).*rescale);

% Alpha parameters
% The following plots visualize the shared variance explained by each
% latent variable in each area.

% Ground truth
alpha_inv_true = 1./paramsTrue.alphas;
% Normalize by the shared variance in each area
alpha_inv_rel_true = alpha_inv_true ./ sum(alpha_inv_true,2);

figure;
hold on;
b = bar(alpha_inv_rel_true');
for groupIdx = 1:numGroups
    b(groupIdx).FaceColor = groupColors{groupIdx};
end
xlabel('Latent variable');
ylabel('Frac. shared var. exp.');
ylim([0 0.5]);
title('Ground truth')

% Estimate
% NOTE: In general, the columns of alpha_est are unordered, and will not
%       match the order in alpha_true.
alpha_inv_est = 1./estParams.alpha.mean;
% Normalize by the shared variance in each area
alpha_inv_rel_est = alpha_inv_est ./ sum(alpha_inv_est,2);

figure;
hold on;
b = bar(alpha_inv_rel_est(:,reorder)');
for groupIdx = 1:numGroups
    b(groupIdx).FaceColor = groupColors{groupIdx};
end
xlabel('Latent variable');
ylabel('Frac. shared var. exp.');
ylim([0 0.5]);
title('Estimate')

%% ========================================================================
% 4) Explore leave-group-out prediction performance and SNRs on train data
% =========================================================================

% Leave-group-out predictive performance
[R2, MSE] = pred_mdlag(seqTrue, estParams);
fprintf('Leave-group-out R^2:\n    %1.4f\n', R2);

% Signal-to-noise ratio of each group, according to estimated mDLAG model
% parameters
snr = computeSNR(estParams.C, estParams.phi);
fprintf('Signal-to-noise ratio of each group:\n    %1.4f  %1.4f  %1.4f\n', snr);

%% =================================================================
% 5) Determine and then visualize the dimensionalities of all types
% ==================================================================

cutoff_sharedvar = 0.02; % Minimum shared variance within a group that a latent must explain
cutoff_snr = 0.001;      % Minimum SNR that a group must have for latents to be significant
[dims,sigDims,varExp,dimTypes] = computeDimensionalities(estParams, ...
                                                         cutoff_sharedvar, ...
                                                         cutoff_snr);
                                              
% Visualize the number of each type of dimension
plotDimensionalities(dims, dimTypes, ...
                     'groupNames', groupNames, ...
                     'plotZeroDim', false);

% Visualize the shared variance explained by each dimension type in each
% group
plotVarExp(varExp, dimTypes, ...
           'groupNames', groupNames, ...
           'plotZeroDim', false);

%% ======================================
% 6) Visualize recovery of GP parameters
% =======================================

% Ground truth
plotGPparams_mdlag(paramsTrue,binWidth,'sigDims',sigDimsTrue,'units',units);

% Estimate
plotGPparams_mdlag(estParams,binWidth,'sigDims',sigDims,'units',units);

%% ============================================
% 7) Visualize recovery of latent time courses
% =============================================

xspec = 'xve'; % 'xve' gives latent time courses scales by shared variance
               % 'xsm' gives latent time courses with normalized variances
                 
% Ground truth
[seqTrue, sortParams] = scaleByVarExp(seqTrue, paramsTrue, alpha_inv_rel_true, ...
                                     'sortDims', false, ...
                                     'sortGroup', 1, ...
                                     'numDim', 10, ...
                                     'indatafield', 'xsm', ...
                                     'outdatafield', 'xve');
plotDimsVsTime_mdlag(seqTrue, xspec, sortParams, binWidth, ...
                     'nPlotMax', 1, ...
                     'plotSingle', false, ...
                     'plotMean', true, ...
                     'units', units, ...
                     'trialGroups', {}, ...
                     'trialColors', {});
                 
% Estimate
[seqEst,~,~] = inferX(seqTrue, estParams);
% NOTE: In general, latent time courses are unordered, and will not match
%       the order in seqTrue.
[seqEst, sortParams] = scaleByVarExp(seqEst, estParams, alpha_inv_rel_est, ...
                                     'sortDims', false, ...
                                     'sortGroup', 1, ...
                                     'numDim', 10, ...
                                     'indatafield', 'xsm', ...
                                     'outdatafield', 'xve');
plotDimsVsTime_mdlag(seqEst, xspec, sortParams, binWidth, ...
                     'nPlotMax', 1, ...
                     'plotSingle', false, ...
                     'plotMean', true, ...
                     'units', units, ...
                     'trialGroups', {}, ...
                     'trialColors', {});
                 
% Overlay estimate on ground truth, for an example trial                 
trialIdx = 1;
seqEst(trialIdx).xve = seqEst(trialIdx).xve([reorder reorder+estParams.xDim reorder+estParams.xDim*2],:) .* repmat(rescale, 1, 3)';                
plotDimsVsTime_mdlag([seqTrue(trialIdx); seqEst(trialIdx)], xspec, paramsTrue, binWidth, ...
                     'nPlotMax', 1, ...
                     'plotSingle', false, ...
                     'plotMean', true, ...
                     'units', units, ...
                     'trialGroups', {[1], [2]}, ...
                     'trialColors', {'k', '#D35FBC'});

%% =============================================================
% 8) Perform a pairwise analysis of interactions between groups
% ==============================================================

% Compute pairwise dimensionalities and shared variances explained
[pairDims,pairVarExp,pairs] = computeDims_pairs(dims,dimTypes,varExp);

% Visualize results
plotDims_pairs(pairDims, pairs, numGroups, 'groupNames', groupNames);
plotVarExp_pairs(pairVarExp, pairs, numGroups, 'groupNames', groupNames);