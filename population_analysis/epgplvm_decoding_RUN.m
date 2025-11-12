function out = epgplvm_decoding_RUN( ...
    spike_times, spike_ids, ...
    tvec, position, speed, track_id, opts)
% Extended P-GPLVM (LMT-style) training + constrained test inference + KNN decode
%  - OU prior on time for X(t)
%  - SE prior on tuning f_n(X) with Poisson likelihood (Laplace)
%  - 5-fold temporally contiguous CV
%  - Final model trained on all RUN bins
%
% Inputs:
%   spike_times [nSp x 1] (s)
%   spike_ids   [nSp x 1] unit IDs
%   tvec        [T x 1]   continuous timebase (s)
%   position    [T x 1]
%   speed       [T x 1]
%   track_id    [T x 1]
%   opts        struct: binSize_run(0.02), runSpeedThresh(5), latentDim(1),
%                      nIter_train(60), nIter_infer(50), KNN_K(25),
%                      seed(1), theta_init (rho,delta,r,ell)
%
% Output:
%   out.cv.folds(f): per-fold decoding (t, X, pos_true/hat, trk_aff)
%   out.cv.summary:  averaged CV metrics (rmse, mae)
%   out.model:       trained-on-all RUN model
%   out.meta:        metadata

%% ---- options
if nargin < 7, opts = struct; end
opts.binSize_train    = fget(opts,'binSize_train',0.1);
opts.binSize_test    = fget(opts,'binSize_test',0.020);
opts.runSpeedThresh = fget(opts,'runSpeedThresh',1);
opts.latentDim      = fget(opts,'latentDim',1);
opts.KNN_K          = fget(opts,'KNN_K',25);
opts.nIter_train    = fget(opts,'nIter_train',60);
opts.nIter_infer    = fget(opts,'nIter_infer',50);
opts.theta_init     = fget(opts,'theta_init',struct('rho',5,'delta',10,'r',2,'ell',20));
opts.seed           = fget(opts,'seed',1);
rng(opts.seed);

%% ---- binning & RUN selection (mean speed per bin; require >50% valid samples)
[tvec, position, speed, track_id] = makecol(tvec, position, speed, track_id);
edges   = (tvec(1):opts.binSize_train:tvec(end))';
centers = edges(1:end-1) + opts.binSize_train/2;
contBinIdx = discretize(tvec, edges);

[Y, unit_ids] = bin_counts(spike_times(:), spike_ids(:), edges);   % [T x N]
[T, N] = size(Y);

isValid = ~isnan(position) & ~isnan(speed) & ismember(track_id,[1 2]);
bin_frac_valid = accum_by_bin(contBinIdx, isValid, @(x) mean(x)>0.5, true);
bin_speed_mean = accum_by_bin(contBinIdx, speed, @mean, false);
isRun = bin_frac_valid & (bin_speed_mean >= opts.runSpeedThresh);

pos_bin = accum_by_bin(contBinIdx, position, @mean, false);
trk_bin = mode_by_bin(contBinIdx, track_id);

Y_run   = Y(isRun,:);   t_run   = centers(isRun);
pos_run = pos_bin(isRun); trk_run = trk_bin(isRun);
T_run = size(Y_run,1);


%% ---- 5 contiguous folds
nFolds = 5;
idx_all = (1:T_run)';
fold_edges = round(linspace(1, T_run+1, nFolds+1));
folds = struct([]);

for f = 1:nFolds
    te = idx_all(fold_edges(f):fold_edges(f+1)-1);
    tr = setdiff(idx_all, te);

    Ytr = Y_run(tr,:); Yte = Y_run(te,:);
    ttr = t_run(tr);   tte = t_run(te);
    postr = pos_run(tr);  poste = pos_run(te);
    trktr = trk_run(tr);  trkte = trk_run(te);

    % ---- Train on TRAIN: decoupled Laplace (tuning GP + time GP) --------------
    m = opts.latentDim;
    X0 = init_latents(Ytr, m);
    theta = opts.theta_init;  % {rho,delta} (tuning) & {r,ell} (time OU)

    [fit, ~] = train_epgplvm_decoupled(Ytr, ttr, X0, theta, opts.nIter_train);

    % ---- Constrained inference on TEST (fix theta + tuning) -------------------
    Xte = infer_constrained_new_trajectory(Yte, tte, fit, opts.nIter_infer);

    % ---- Decode via KNN on TRAIN latent bank ---------------------------------
    [pos_hat, trk_aff] = knn_decode(Xte, fit.X_train, postr, trktr, opts.KNN_K);

    folds(f).fold = f;
    folds(f).t = tte;
    folds(f).X = Xte;
    folds(f).pos_true = poste;
    folds(f).pos_hat  = pos_hat;
    folds(f).trk_true = trkte;
    folds(f).trk_aff  = trk_aff;
    folds(f).rmse = sqrt(mean((pos_hat - poste).^2,'omitnan'));
    folds(f).mae  = mean(abs(pos_hat - poste),'omitnan');
end

cv_rmse = mean([folds.rmse]);
cv_mae  = mean([folds.mae]);

%% ---- Final model trained on ALL RUN ------------------------------------------
X0_all = init_latents(Y_run, opts.latentDim);
[fit_all, train_diag] = train_epgplvm_decoupled(Y_run, t_run, X0_all, opts.theta_init, opts.nIter_train); %Train using all data

model = struct('fit',fit_all, ...
               'X_bank',fit_all.X_train, ...
               'pos_bank',pos_run(:), ...
               'trk_bank',trk_run(:), ...
               'unit_ids',unit_ids, ...
               'binSize_train',opts.binSize_train, ...
               't_run',t_run);

out = struct();
out.cv.folds   = folds;
out.cv.summary = struct('rmse',cv_rmse,'mae',cv_mae);
out.model      = model;
out.meta       = struct('unit_ids',unit_ids,'edges',edges,'isRun',isRun,'opts',opts);
end

%% ============================ Core training (LMT/2017 + 2024) ============================
function [fit, train_diag] = train_epgplvm_decoupled(Y, tvec, Xinit, theta0, nIter)
% Decoupled Laplace:
%   (1) given X, for each neuron n, find Laplace mode f_n(X) under Poisson
%       with SE kernel Kx(X,X | rho,delta)
%   (2) given {f_n}, update X with OU time prior via Gauss–Newton step
%   (3) ML over theta = {rho,delta,r,ell} using marginal-Laplace objective (optional light)

X = Xinit; theta = theta0;
[T, N] = size(Y); m = size(X,2);

for it = 1:nIter
    % --- (1) per-neuron Laplace (Poisson) to get f_hat and alpha = y - mu at the mode
    [f_all, alpha_all, Wdiag_all] = laplace_tuning_all(Y, X, theta.rho, theta.delta);

    % --- (2) update X with OU prior using correct gradient & GN diagonal
    % g_like = sum_n J_n(X)' * (y_n - mu_n),   H_like ≈ sum_n J_n.^2 * mu_n  (per element)
    [g_like, H_like] = like_grad_hdiag_X(X, f_all, alpha_all, Wdiag_all, theta.rho, theta.delta);

    % prior over X: OU kernel in time (per latent dim, independent)
    Kt  = ou_from_t(tvec, theta.r, theta.ell);         % [T x T]
    iKt = chol_inv(Kt);                                 % stable inverse via chol
    g_prior   = iKt * X;                               % [T x m]
    H_prior_d = repmat(diag(iKt),1,m);                 % diag approx

    step = (g_like + g_prior) ./ max(H_like + H_prior_d, 1e-8);
    X = X - step;

    % --- (3) quick theta refinement (small #f-evals)
    theta = refine_theta(Y, X, tvec, theta);

    % early stop if tiny movement
    if mean(abs(step(:))) < 1e-5, break; end
end

fit = struct('theta',theta, ...
             'X_train',X, ...
             'tvec_train',tvec, ...
             'rho',theta.rho,'delta',theta.delta, ...
             'r',theta.r,'ell',theta.ell);
train_diag = struct('iters',it,'theta',theta);
end

function Xnew = infer_constrained_new_trajectory(Ynew, t_new, fit, nIter)
% Constrained test-time MAP:
%   fix theta and the fitted tuning GPs (via Laplace on current X),
%   then optimize X_new with OU prior + correct Poisson gradient.

[T,~] = size(Ynew); m = size(fit.X_train,2);
X = init_latents(Ynew, m);   % warm-start

for it=1:nIter
    % Laplace per neuron, given current X and fixed {rho,delta}
    [f_all, alpha_all, Wdiag_all] = laplace_tuning_all(Ynew, X, fit.rho, fit.delta);

    % gradient/diag Hessian w.r.t X using Jacobians of tuning GPs
    [g_like, H_like] = like_grad_hdiag_X(X, f_all, alpha_all, Wdiag_all, fit.rho, fit.delta);

    % OU prior with fixed r,ell on NEW time grid
    Kt  = ou_from_t(t_new, fit.r, fit.ell);
    iKt = chol_inv(Kt);
    g_prior   = iKt * X;
    H_prior_d = repmat(diag(iKt),1,m);

    step = (g_like + g_prior) ./ max(H_like + H_prior_d, 1e-8);
    X = X - step;

    if mean(abs(step(:))) < 1e-5, break; end
end
Xnew = X;
end

%% ============================ Likelihood pieces ==============================
function [f_all, alpha_all, Wdiag_all] = laplace_tuning_all(Y, X, rho, delta)
% Per-neuron Laplace to get:
%  f_n (T x 1)  : mode under Poisson with GP prior over X (SE kernel)
%  alpha_n      : y - mu at the mode (useful: f = K * alpha)
%  Wdiag_n      : diag(exp(f))  (for Gauss–Newton)
[T,N] = size(Y);
f_all      = zeros(T,N);
alpha_all  = zeros(T,N);
Wdiag_all  = zeros(T,N);

K = se_XK(X, X, rho, delta);            % shared across neurons at given X
iK = chol_inv(K);

for n = 1:N
    y = Y(:,n);
    f = log1p(y);                       % initialize at log(count+1)
    for it=1:30
        mu = exp(f);
        g  = y - mu - iK*f;             % gradient wrt f
        W  = mu;                        % diag(W) = mu under Poisson link
        H  = iK + diag(W);              % negative Hessian: -(iK + W)
        % Solve H * step = g   (Newton step)
        step = solve_pd(H, g);
        f = f + step;
        if max(abs(step)) < 1e-6, break; end
    end
    mu = exp(f);
    alpha = y - mu;                     % from K^{-1} f = y - mu at optimum
    f_all(:,n) = f;
    alpha_all(:,n) = alpha;
    Wdiag_all(:,n) = mu;
end
end

function [g_like, H_like] = like_grad_hdiag_X(X, f_all, alpha_all, Wdiag_all, rho, delta)
% Computes gradient and Gauss–Newton diagonal Hessian of the Poisson log-likelihood wrt X.
% IMPORTANT: jacobian_f_wrt_X returns J = ∂f/∂X (already alpha-weighted),
% so we should NOT multiply by alpha again here.

[T,m] = size(X);
N = size(f_all,2);
g_like = zeros(T,m);
H_like = zeros(T,m);

for n = 1:N
    alpha = alpha_all(:,n);        %#ok<NASGU> % not used directly in this corrected form
    mu    = Wdiag_all(:,n);        % mu = exp(f) at the Laplace mode
    % J(t,d) = ∂f(t)/∂X(t,d) = sum_j alpha_j * ∂k(x_t,x_j)/∂x_t(d)
    J = jacobian_f_wrt_X(X, alpha_all(:,n), rho, delta);

    % Gradient: ∂L/∂X = (∂f/∂X) * (y - μ) ; our J already includes (y - μ) via alpha, so just add J.
    g_like = g_like + J;

    % GN diag Hessian: sum_n μ .* (∂f/∂X)^2  (since ∂μ/∂X = μ * ∂f/∂X for exp link)
    H_like = H_like + (J.^2) .* mu;
end
end

function J = jacobian_f_wrt_X(X, alpha, rho, delta)
% X: [T x m], alpha: [T x 1]
% J(t,:) = sum_j alpha_j * ∂k(x_t, x_j)/∂x_t
T = size(X,1); m = size(X,2);
J = zeros(T,m);
for t=1:T
    dif = X(t,:) - X;                 % [T x m]
    sq = sum(dif.^2,2);
    k  = rho * exp(-0.5*sq/max(delta,1e-9)^2);            % [T x 1]
    % ∂k/∂x_t = k .* (-(dif)/delta^2)  (note dif = x_t - x_j)
    dk = (k ./ max(delta,1e-9)^2) .* (-dif);              % [T x m]
    J(t,:) = alpha' * dk;                                  % 1 x m
end
end

%% ============================ Priors / kernels ===============================
function Kt = ou_from_t(tvec, r, ell)
tvec=tvec(:);
Dt = abs(tvec - tvec');               % |t_i - t_j|
Kt = r * exp(-Dt/max(ell,1e-9));
Kt = Kt + 1e-6*eye(size(Kt));         % jitter
end

function K = se_XK(X1, X2, rho, delta)
D2 = pdist2(X1,X2,'euclidean').^2;
K = rho * exp(-0.5 * D2 / max(delta,1e-9)^2);
K = K + 1e-6*(size(K,1)==size(K,2))*eye(size(K,1)); % jitter if square
end

%% ============================ Theta refinement ===============================
function theta = refine_theta(Y, X, tvec, theta0)
% A light ML refinement: a few fmincon steps on theta = [rho,delta,r,ell]
% using marginal Laplace objective; keep small to stay robust.
lb=[1e-4,1e-4,1e-4,1e-4]; ub=[1e3,1e3,1e3,1e3];
th0 = log([theta0.rho,theta0.delta,theta0.r,theta0.ell]);
obj = @(th) - laplace_marglik_theta(exp(th),Y,X,tvec);   % negative log-like
opts = optimoptions('fmincon','Display','off','Algorithm','interior-point', ...
    'SpecifyObjectiveGradient','off','MaxFunctionEvaluations',40, ...
    'OptimalityTolerance',1e-3,'StepTolerance',1e-3);
th = fmincon(obj, th0, [],[],[],[], log(lb), log(ub), [], opts);
th = exp(th);
theta = struct('rho',th(1),'delta',th(2),'r',th(3),'ell',th(4));
end

function L = laplace_marglik_theta(th, Y, X, tvec)
% Same structure as your original, but using our tuning Laplace.
rho=th(1); delta=th(2); r=th(3); ell=th(4);
[T,N]=size(Y); m=size(X,2);

Kx = se_XK(X,X,rho,delta); iKx = chol_inv(Kx);
Kt = ou_from_t(tvec,r,ell); iKt = chol_inv(Kt);

loglike_sum=0; logdet_sum=0; quad_sum=0;
for n=1:N
    y=Y(:,n); f=log1p(y);
    for it=1:30
        mu=exp(f); g=y - mu - iKx*f; H=iKx + diag(mu);
        step = solve_pd(H,g); f=f+step;
        if max(abs(step))<1e-6, break; end
    end
    mu=exp(f); alpha=y-mu;
    % Laplace terms
    ll = sum(y.*f - mu);
    SW = sqrt(mu); D = diag(SW); A = eye(T) + D*(Kx*D);
    logdetA = logdet_pd(A);
    quad = -0.5*(f'*(iKx*f));
    loglike_sum = loglike_sum + ll;
    logdet_sum  = logdet_sum  - 0.5*logdetA;
    quad_sum    = quad_sum + quad;
end

prior_quad=0;
for d=1:m
    xd=X(:,d); prior_quad = prior_quad - 0.5*(xd'*(iKt*xd));
end
prior_logdet = -0.5*m*logdet_pd(Kt);

L = loglike_sum + logdet_sum + quad_sum + prior_quad + prior_logdet;
end

%% ============================ Decode / metrics ===============================
function [pos_hat, trk_aff] = knn_decode(Xq, Xbank, pos_bank, trk_bank, K)
[idx,D] = knnsearch(Xbank, Xq, 'K',K);
w = 1./(D + 1e-9); w = w ./ sum(w,2);
pos_hat = sum(w .* pos_bank(idx),2);
nn = trk_bank(idx);
trk_aff = [mean(nn==1,2), mean(nn==2,2)];
end

%% ============================ Utils =========================================
function [Y, unit_ids] = bin_counts(st, sid, edges)
unit_ids = unique(sid(:)'); nU=numel(unit_ids); nB=numel(edges)-1;
Y = zeros(nB,nU,'double'); b = discretize(st, edges);
for u=1:nU
    bu = b(sid==unit_ids(u));
    bu = bu(isfinite(bu) & bu>=1 & bu<=nB);
    if ~isempty(bu), Y(:,u) = accumarray(bu,1,[nB 1]); end
end
end

function y = accum_by_bin(binIdx, v, fun, logicalOut)
if nargin<4, logicalOut=false; end
K = max(binIdx(~isnan(binIdx)));
y = nan(K,1);
for k=1:K
    ii = (binIdx==k) & ~isnan(v);
    if any(ii)
        if logicalOut, y(k) = fun(v(ii)); else, y(k)=fun(v(ii)); end
    end
end
if logicalOut, y = logical(y); end
end

function m = mode_by_bin(binIdx, v)
K = max(binIdx(~isnan(binIdx))); m = nan(K,1);
for k=1:K
    ii=(binIdx==k) & ~isnan(v);
    if any(ii), m(k)=mode(v(ii)); end
end
end

function X0 = init_latents(Y, d)
% INIT_LATENTS  Returns [T x d] latent init from RUN data Y [T x N].
% Prefer PLDS if available; else PCA(log1p(Y)).
% Columns of X0 are z-scored.

T = size(Y,1);   % time bins
N = size(Y,2);   % No of Neurons

if exist('run_plds','file')==2 % if LMT repo is present 
    % run_plds expects Time x Neuron
    Xraw = run_plds(Y, d);

    % Unwrap common output formats
    if isstruct(Xraw)
        % Try common field names
        for fn = ["X","xs","latent","Z","z","x","factors","states"]
            if isfield(Xraw, fn)
                Xraw = Xraw.(fn);
                break;
            end
        end
    end

    % Now coerce to [T x d]
    if ismatrix(Xraw)
        % If it's [d x T], transpose
        if size(Xraw,2) == T && size(Xraw,1) ~= T
            X0 = Xraw';
        elseif size(Xraw,1) == T
            X0 = Xraw;
        else
            error('init_latents: PLDS output has incompatible size [%d x %d], expected T=%d in one dimension.', ...
                   size(Xraw,1), size(Xraw,2), T);
        end
    else
        error('init_latents: Unrecognized PLDS return type.');
    end

    % Trim/pad dims
    if size(X0,2) > d
        X0 = X0(:,1:d);
    elseif size(X0,2) < d
        X0(:,end+1:d) = 0;
    end

else
    % ---------- GPFA-lite fallback (PCA + OU smoothing) ----------
    % 1) Variance stabilization (Anscombe) and neuron-wise z-score
    Z = 2*sqrt(max(Y,0) + 3/8);
    mu = mean(Z,1,'omitnan');
    sg = std(Z,0,1,'omitnan'); sg = max(sg,1e-9);
    Z = (Z - mu) ./ sg;

    % 2) Get top-d temporal scores via PCA (memory-light if svds exists)
    got = false;
    if exist('svds','file')==2
        try
            [U,S,~] = svds(Z, min(d, min(T,size(Z,2))));
            Scores = U * S;                 % [T x d]
            got = true;
        catch
        end
    end
    if ~got
        try
            [U,S,~] = svd(Z,'econ');
            r = min([d, size(U,2), size(S,1)]);
            Scores = U(:,1:r) * S(1:r,1:r); % [T x r]
            if r < d, Scores(:,end+1:d) = 0; end
            got = true;
        catch
            % Final fallback to stats toolbox PCA if available
            [~,Smat,~] = pca(Z,'Algorithm','svd','Centered',false, ...
                               'NumComponents', min(d, min(size(Z))));
            Scores = Smat(:,1:min(d,size(Smat,2)));
            if size(Scores,2) < d, Scores(:,end+1:d) = 0; end
            got = true;
        end
    end

    % 3) OU GP smoothing of each score (PLDS-ish temporal prior)
    %    (Skip smoothing if T very large to avoid T×T memory.)
    doSmooth = (T <= 20000);         % adjust threshold if needed
    ellBins  = max(5, round(T/50));  % OU timescale in bins
    lamRidge = 1e-2;                 % small observation noise
    if doSmooth
        t = (1:T)'; D = abs(t - t');
        K = exp(-D / max(ellBins,1e-9));   % OU kernel
        I = speye(T);
        % Posterior mean under GP prior with tiny noise:
        % X0(:,j) = (K + lam*I) \ (K * yj)
        X0 = zeros(T, d);
        for j = 1:d
            yj = Scores(:,j);
            X0(:,j) = (K + lamRidge*I) \ (K * yj);
        end
    else
        X0 = Scores;
    end
end

% Standardize columns
X0 = (X0 - mean(X0,1)) ./ max(std(X0,0,1), 1e-9);

% Final safety check
assert(size(X0,1) == T && size(X0,2) == d, ...
    'init_latents: X0 must be [T x d]; got [%d x %d], expected [%d x %d].', ...
    size(X0,1), size(X0,2), T, d);
end


function Ainv = chol_inv(A)
L = chol((A+A')/2 + 1e-9*eye(size(A)),'lower');
Ainv = L'\(L\eye(size(A)));
end

function x = solve_pd(H, b)
% robust solve for SPD-ish H
[L,p] = chol((H+H')/2 + 1e-12*eye(size(H)),'lower');
if p==0
    x = L'\(L\b);
else
    x = pinv(H)*b;
end
end

function v = iff(c,a,b), if c, v=a; else, v=b; end
end

function val = logdet_pd(A)
% log|A| for SPD-ish A
[L,p] = chol((A+A')/2 + 1e-12*eye(size(A)),'lower');
if p==0
    val = 2*sum(log(diag(L)));
else
    e = eig((A+A')/2); e = max(e,1e-12); val = sum(log(e));
end
end

function [a,b,c,d]=makecol(a,b,c,d), a=a(:); b=b(:); c=c(:); d=d(:); end

function v = fget(s,f,d), v=d; if isfield(s,f)&&~isempty(s.(f)), v=s.(f); end, end
