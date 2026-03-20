function [estParams,trackedParams,flags] = mDLAG_fitting(seqTrue,binWidth,yDims,options)

% save_all_figures('C:\Users\masahiro.takigawa\Documents\GitHub\mDLAG\mDLAG\demo\results',[])
%% =====================================
% 0) spike counts
% ======================================
if ~isfield(options,'load_fitted_model')
    options.load_fitted_model=[];
end

if isempty(options.load_fitted_model)
    disp('run id not found. check if options.load_fitted_model is empty')
else
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

save(fullfile(options.ANALYSIS_DATAPATH,'mDLAG_result',sprintf('mDLAG_%s_run%i.mat',options.eventname,options.load_fitted_model)),'estParams','trackedParams','flags');

% DIR = dir(fullfile(options.ANALYSIS_DATAPATH,'mDLAG_result',sprintf('*mDLAG*_%s.mat',options.eventname)));
% 
% all_files= [];
% run_num = [];
% if ~isempty(DIR)
%     for  nfile = 1:length(DIR)
% 
%         all_files{nfile} = DIR(nfile).name;
% 
%         run_num(nfile) =  str2num(cell2mat(extractBetween( DIR(nfile).name,'run','.mat')));
%     end
% else
%     all_files= '';
% end
% 
% % fullfile(options.ANALYSIS_DATAPATH,'mDLAG_result'
% if sum(contains(all_files,'run1.mat'))==0
%     save(fullfile(options.ANALYSIS_DATAPATH,'mDLAG_result',sprintf('mDLAG_%s_run1.mat',options.eventname)),'estParams','trackedParams','flags');
% else
%     save(fullfile(options.ANALYSIS_DATAPATH,'mDLAG_result',sprintf('mDLAG_%s_run%i.mat',options.eventname,max(run_num)+1)),'estParams','trackedParams','flags');
% end

end

