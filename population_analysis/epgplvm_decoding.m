function sleep = epgplvm_decoding( ...
    spike_times, spike_ids, ...
    replay_onset, replay_offset, ...
    model, opts)
% Decode replay/PBE events using trained EP-GPLVM model
%
% Inputs:
%   spike_times, spike_ids : spikes (s, IDs)
%   replay_onset/offset    : [K x 1] event times (s)
%   model                  : output of epgplvm_decoding_RUN
%   opts                   : (optional) decoding options
%
% Output:
%   sleep(k): decoded latent trajectory and metrics per event

if nargin<6, opts = struct; end
opts.pbe_binSize = fget(opts,'pbe_binSize',0.005);
opts.nIter_infer = fget(opts,'nIter_infer',50);
opts.KNN_K       = fget(opts,'KNN_K',25);
opts.nShuffles   = fget(opts,'nShuffles',200);

tmin = min(replay_onset); tmax = max(replay_offset);
edges = (tmin:opts.pbe_binSize:tmax)';
Y_all = bin_spike_timewise(spike_times(:), spike_ids(:), edges, model.unit_ids);
tvec_pbe = edges(1:end-1) + opts.pbe_binSize/2;

% Adjust tuning for bin size difference
fit_pbe = rescale_fit_for_pbe(model.fit, model.binSize_run, opts.pbe_binSize);

K = numel(replay_onset);
sleep = repmat(struct('event_t',[],'X',[],'pos_hat',[],'trk_aff',[], ...
                      'LLH',[],'zLLH_circ',[],'zLLH_time',[], ...
                      'spatial_consistency',[],'step_distance',[]),K,1);

for k=1:K
    on = replay_onset(k); off = replay_offset(k);
    in_ev = tvec_pbe>=on & tvec_pbe<=off;
    if ~any(in_ev), continue; end
    Yev = double(Y_all(in_ev,:)); tvec_ev = tvec_pbe(in_ev);

    Xev = infer_new_latents_constrained(Yev, fit_pbe, tvec_ev, opts.nIter_infer);
    [pos_ev, trk_aff_ev] = knn_manifold_decode(Xev, model.X_bank, model.pos_bank, model.trk_bank, opts.KNN_K);

    LLH = path_loglik_poisson(Yev, Xev, fit_pbe);
    [zC,zT] = llh_shuffle_zscores(Yev, fit_pbe, tvec_ev, opts.nShuffles);
    Scons = spatial_consistency(Xev, model.X_bank);
    stepD = mean(vecnorm(diff(Xev,1,1),2,2),'omitnan');

    sleep(k).event_t = [on off];
    sleep(k).X = Xev;
    sleep(k).pos_hat = pos_ev;
    sleep(k).trk_aff = trk_aff_ev;
    sleep(k).LLH = LLH;
    sleep(k).zLLH_circ = zC;
    sleep(k).zLLH_time = zT;
    sleep(k).spatial_consistency = Scons;
    sleep(k).step_distance = stepD;
end
end


%% ============================ Helpers ============================
function v = fget(s,f,d), v=d; if isfield(s,f)&&~isempty(s.(f)), v=s.(f); end, end

function Y = bin_spike_timewise(st, sid, edges, unit_ids)
nU = numel(unit_ids); nB = numel(edges)-1;
Y = zeros(nB,nU,'single'); b = discretize(st, edges);
for i=1:nU
    mask = (sid==unit_ids(i)); bu = b(mask);
    bu = bu(isfinite(bu) & bu>=1 & bu<=nB);
    if ~isempty(bu), Y(:,i) = Y(:,i) + accumarray(bu,1,[nB 1]); end
end
Y = double(Y);
end

function fit2 = rescale_fit_for_pbe(fit, bin_run, bin_pbe)
fit2 = fit;
log_scale = log(max(bin_run/bin_pbe,1e-12));
old_pred = fit.tuning.pred; old_jac = fit.tuning.jac;
fit2.tuning.pred = @(Xq) old_pred(Xq) + log_scale;
fit2.tuning.jac  = @(Xq) old_jac(Xq);
end

function LLH = path_loglik_poisson(Y, X, fit)
Mu = exp( fit.tuning.pred(X) );
LLH = sum(Y(:).*log(Mu(:)+1e-12) - Mu(:));
end

function [zCirc, zTime] = llh_shuffle_zscores(Y, fit, tvec, S)
LL = zeros(S,2);
for s=1:S
    sh = randi(size(Y,1)-1);
    Yc = circshift(Y,sh,1);
    Xc = infer_new_latents_constrained(Yc, fit, tvec, 10);
    LL(s,1) = path_loglik_poisson(Yc, Xc, fit);

    ord = randperm(size(Y,1));
    Yr = Y(ord,:); Xt = tvec(ord);
    Xr = infer_new_latents_constrained(Yr, fit, Xt, 10);
    LL(s,2) = path_loglik_poisson(Yr, Xr, fit);
end
X0 = infer_new_latents_constrained(Y, fit, tvec, 10);
LL0 = path_loglik_poisson(Y, X0, fit);
zCirc = (LL0 - mean(LL(:,1))) / (std(LL(:,1))+1e-9);
zTime = (LL0 - mean(LL(:,2))) / (std(LL(:,2))+1e-9);
end

function Xnew = infer_new_latents_constrained(Ynew, fit, tvec_new, nIter)
[T,~]=size(Ynew); d=fit.dim; X=pca_init(Ynew,d);
for it=1:nIter
    iKt = pinv_stable(ou_kernel_from_t(tvec_new, fit.theta.r, fit.theta.ell));
    [~, gX_like, Hdiag_like] = poisson_ll_grad_hdiag(Ynew, X, fit.tuning);
    g_prior = iKt * X; Hdiag_p = repmat(diag(iKt),1,d);
    step = (gX_like + g_prior) ./ max(Hdiag_like + Hdiag_p,1e-6);
    X = X - step;
end
Xnew=X;
end

function X0 = pca_init(Y, d)
Y=double(Y); Y(~isfinite(Y))=0; Y(Y<0)=0; Ylog=log1p(Y);
sd=std(Ylog,0,1,'omitnan'); keep=sd>0 & isfinite(sd); Z=Ylog(:,keep)-mean(Ylog(:,keep),1,'omitnan');
try
    [~,score]=pca(Z,'Algorithm','svd','Centered',false,'NumComponents',min(d,min(size(Z))));
    X0=score(:,1:min(d,size(score,2)));
catch
    [U,S,~]=svd(Z,'econ'); r=min([d,size(U,2),size(S,1)]); X0=U(:,1:r)*S(1:r,1:r);
end
X0=(X0-mean(X0,1))./max(std(X0,0,1),1e-9);
end

function [ll, gX, Hdiag] = poisson_ll_grad_hdiag(Y, X, F)
Mu=exp(F.pred(X)); ll=sum(Y(:).*log(Mu(:)+1e-12)-Mu(:)); %#ok<NASGU>
J=F.jac(X); gX=J; Hdiag=repmat(mean(Mu,2),1,size(X,2));
end

function K = ou_kernel_from_t(tvec, r, ell)
tvec=tvec(:); Dt=abs(tvec - tvec'); K=r*exp(-Dt/max(ell,1e-8)); K=K+1e-6*eye(size(K));
end
function iK = pinv_stable(K), iK=(K+1e-6*eye(size(K)))\eye(size(K)); end

function [pos_hat,trk_aff]=knn_manifold_decode(Xq,Xbank,pos_bank,trk_bank,K)
idx=knnsearch(Xbank,Xq,'K',K);
D=pdist2(Xq,Xbank); [~,ord]=sort(D,2,'ascend'); ord=ord(:,1:K);
w=1./(D(sub2ind(size(D),(1:size(D,1))',ord))+1e-6); w=w./sum(w,2);
pos_hat=sum(w.*pos_bank(ord),2);
trk_neighbors=trk_bank(idx); trk_aff=[mean(trk_neighbors==1,2), mean(trk_neighbors==2,2)];
end
