% Generic fitting function for speed tuning 
function fitSpeedTuning(forfit_sfs,forfit_tfs,forfit_f0array)

% If we included a SF of 0, the fitting function doesn't handle that, so
% need to make it a small number
forfit_sfs(forfit_sfs <0.001) = 0.001;

% Get preferred
[A_f0 A_f0_ix] = max(forfit_f0array(:));
[A_f1 A_f1_ix] = max(forfit_f1array(:));
forfit_prefSF_f0 = forfit_sfs(A_f0_ix);
forfit_prefTF_f0 = forfit_tfs(A_f0_ix);

% To visualise fit predictions, interpolate axis use meshgrid (note, should put this in log2 space not
% matlabs generic log10)
x_axis_tfs = 2.^linspace(log2(min(forfit_tfs(:))),log2(max(forfit_tfs(:))),128);
tp = unique(forfit_sfs(:));
if tp(1) < tp(2)/8
    % If we included an SF of near 0, we need to make a log space that
    % doesn't include that...
    tpl = tp(2)/2;
    x_axis_sfs = 2.^linspace(log2(tpl),log2(max(tp)),128);
else
    x_axis_sfs = 2.^linspace(log2(min(tp)),log2(max(tp)),128);
end
[interp_x_axis_sfs,interp_x_axis_tfs] = meshgrid(x_axis_sfs,x_axis_tfs);

% Get starting guess for speed model F0
initparamsF0 = [spont_f0array A_f0 forfit_prefSF_f0 forfit_prefTF_f0 1/forfit_prefSF_f0 2 0];
sftf = [forfit_sfs(:),forfit_tfs(:)];
obsresp = forfit_f0array(:);

% Get lsq fit
[fitparams,pve] = getSpeedTuningFitLsq(sftf,obsresp,initparamsF0);

% Get predictions for actual measurement points
pred_f0 = speedmodelInt(fitparams,sftf);
pred_f0 = reshape(pred_f0,6,7)'; % to keep plotting consistent with F0 and F1 matrices in main plots

% Get predictions for interpolated measurement pointss
interp_f0 = speedmodelInt(fitparams,[interp_x_axis_sfs(:),interp_x_axis_tfs(:)]);
interp_f0 = reshape(interp_f0,sqrt(length(interp_f0)),sqrt(length(interp_f0)))';  % to keep plotting consistent with F0 and F1 matrices in main plots
end

function [fitparams,pve] = getSpeedTuningFitLsq(sftf,obsresp,initparams)
% Set up fit including bounds
lb = [0 0 0.0001 0.1 0.1 0.1 -2];
ub = [initparams(1)*2 1000 0.5 30 5 5 2];
options = optimset('TolFun',1e-5,'TolX',1e-4,'Maxiter',100000,'MaxFunEvals',100000, ...
    'Display','off');

% Fit
x = initparams;
xdata = sftf;
ydata = obsresp;
x = lsqcurvefit(@(x,xdata) speedmodelInt2(x,xdata),x,xdata,ydata,lb,ub,options);
fitparams = x;

% Get predicted response
B = fitparams(1);
A = fitparams(2);
prefSF = fitparams(3);
prefTF = fitparams(4);
SFsd = fitparams(5);
TFsd = fitparams(6);
slope = fitparams(7);
TFs = xdata(:,2);
SFs = xdata(:,1);

% generate preferred TF for each specific SF
bestTFs = slope.*(log2(SFs)-log2(prefSF))+log2(prefTF);
term1 = exp(-1.*(   (log2(SFs)-log2(prefSF)).^2  ./ (2.*(SFsd.^2)) ));
term2 = exp(-1.*(   (log2(TFs)-bestTFs).^2 ./ (2.*(TFsd.^2)) ));
predresp = B+A.*term1.*term2;

tpredresp = predresp;
tobsresp = obsresp;
pve = 100*(1-mean((tpredresp-tobsresp).^2)/mean((mean(tobsresp)-tobsresp).^2));


    function [modelout,SFs,TFs] = speedmodelInt2(param,xvals)
        
        if nargin < 2
            [SFs, TFs] = meshgrid([0 0.0094 0.0187 0.0375 0.0750 0.1500 0.3000],[0.4688 0.9375 1.8751 3.7502 7.5003 15.0006]);
            TFs = TFs(:);
            SFs = SFs(:);
        else
            TFs = xvals(:,2);
            SFs = xvals(:,1);
        end
        if nargin < 1
            param = [2 10 0.05 3 1 2 1];
        end
        
        B  = param(1);
        A = param(2);
        prefSF = param(3);
        prefTF = param(4);
        SFsd = param(5);
        TFsd = param(6);
        slope = param(7);
        
        % generate preferred TF for each specific SF
        bestTFs = slope.*(log2(SFs)-log2(prefSF))+log2(prefTF);
        
        term1 = exp(-1.*(   (log2(SFs)-log2(prefSF)).^2  ./ (2.*(SFsd.^2)) ));
        term2 = exp(-1.*(   (log2(TFs)-bestTFs).^2 ./ (2.*(TFsd.^2)) ));
        
        modelout = B+A.*term1.*term2;
    end
end

% function modelout = speedmodel(param,sfs,tfs)
% Model of speed tuning, as in Gale & Murphy 2014
% Inputs:
%		Parameters:	max Response (A)
%                   preferred SF over all TFs (prefSF)
%					preferred TF over all SFs (prefTF)
%					SFsd
%                   TFsd
%                   linear relationship of prefTF and stimSF in log SF space (slope)
%
% Calculate a two dimensional gaussian function in which the preferred TF
% can depend on the SF of the stimulus
% History: GdF 22 March 2014 Wrote it

function [modelout,SFs,TFs] = speedmodelInt(param,xvals)

if nargin < 2
    [SFs, TFs] = meshgrid([0 0.0094 0.0187 0.0375 0.0750 0.1500 0.3000],[0.4688 0.9375 1.8751 3.7502 7.5003 15.0006]);
    TFs = TFs(:);
    SFs = SFs(:);
else
    TFs = xvals(:,2);
    SFs = xvals(:,1);
end
if nargin < 1
    param = [2 10 0.05 3 1 2 1];
end

B  = param(1);
A = param(2);
prefSF = param(3);
prefTF = param(4);
SFsd = param(5);
TFsd = param(6);
slope = param(7);

% generate preferred TF for each specific SF
bestTFs = slope.*(log2(SFs)-log2(prefSF))+log2(prefTF);

term1 = exp(-1.*(   (log2(SFs)-log2(prefSF)).^2  ./ (2.*(SFsd.^2)) ));
term2 = exp(-1.*(   (log2(TFs)-bestTFs).^2 ./ (2.*(TFsd.^2)) ));

modelout = B+A.*term1.*term2;
end