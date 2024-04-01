function mDLAG_analysis_pipeline(seqTrue,binWidth,yDims)


%% =====================================
% 0) spike counts
% ======================================


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
