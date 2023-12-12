% function [fitparams, modelout, options] = ibn_fit2DGaussLsq(forfit_f0,forfit_xval,forfit_yval,options)
% Function to fit 2D gaussian to generic input data
%
% Inputs -  all inputs are 2D matrices (may have unequal numbers of rows
%           and columns)
%   forfit_f0 - this is the amplitude (of LFP, spike rate, ca2+ signal etc)
%   forfit_xval - these are x positions of the centre of the stimulus
%       (usually in degrees, but whatever the parameter, that will  be the unit
%        returned)
%   forfit_yval - these are the y positions of the above 
%   options - 
%           rangeBg % eg = 2; ie from -2:2 x the minimum
%           sdInit % eg 3; ie the initial sd in multiples of the spacing of
%                       stimuli
%           interpModelOutMultiplier- multiplier to make interpolated model
%           for output (comes out in options)
% Outputs
%   fitparams [centre x, centre y, gain, 1 sd, maintained signal]
%   modelout [same dimensions as input]
%
% Dependencies
%   lsqcurvefit
%
% History 
%   TW and SGS 11th Feb 2020 wrote it (adapted from : )
%   TW added lines in the centre of mass function
%   SZ added the R2 as an output and at line 70 to use it for evaluating
%   how good the gaussianf fit is

function [fitparams, modelout, options,R2,initparamsBl] = ibn_fit2DGaussLsq(forfit_f0,forfit_xval,forfit_yval,options)
if ~exist('options','var')
    sdInit = 2;
    rangeBg = 2; % eg = 2; ie from half to twice the minimum
    interpModelOutMultiplier = 10;
else
    rangeBg = options.rangeBg;
    sdInit = options.sdInit;
    interpModelOutMultiplier = options.interpModelOutMultiplier;
end

% Calculate spacing of stimuli
spacing = mean(diff(forfit_xval(1,:)));  % ie the everage spacing in x

% First use centre of mass to estimate likely position
% [col,row] = centreOfMass2DInt(forfit_f0);
StartPosition=[min(forfit_xval(:))+2*spacing, min(forfit_yval(:))+2*spacing; ...
               min(forfit_xval(:))+2*spacing, max(forfit_yval(:))-2*spacing; ...
    max(forfit_xval(:))-2*spacing, min(forfit_yval(:))+2*spacing; ...
    max(forfit_xval(:))-2*spacing, max(forfit_yval(:))-2*spacing];

% For Gaussian fit
initparamsBl(1) = StartPosition(1,1);         % centre pos, x
initparamsBl(2) = StartPosition(1,2);         % centre pos, y
initparamsBl(3) = max(forfit_f0(:));           % peak gain (gain of centre)
initparamsBl(4) = sdInit*spacing;               % centre size 1*sd (in units of stimulus spacing)
initparamsBl(5) = median(forfit_f0(:));            % background activity 


% Set up fit including bounds
% Generate a vector for the lower and upper bounds (lb and ub respectively)
% of different parameters. [ centreX,centreY,gain of centre,centre size, background activity]
options.lb = [min(forfit_xval(:))-initparamsBl(4), min(forfit_yval(:))-initparamsBl(4),...
    -1*max(forfit_f0(:))*10, 0.5*spacing, initparamsBl(5)-initparamsBl(5)*rangeBg+0.01]; %
options.ub = [max(forfit_xval(:))+initparamsBl(4), max(forfit_yval(:))+initparamsBl(4),...
    max(forfit_f0(:))*10, 4*spacing, initparamsBl(5)+initparamsBl(5)*rangeBg+0.01]; % Just centre, no surround
fitFuncoptions = optimset('TolFun',1e-5,'TolX',1e-4,'Maxiter',100000,'MaxFunEvals',100000, ...
    'Display','off');

% Fit
xinputs = [forfit_xval(:) forfit_yval(:)];
ydata = forfit_f0(:);

for thisstart=1:size(StartPosition,1)
% For Gaussian fit
initparamsBl(1) = StartPosition(thisstart,1);         % centre pos, x
initparamsBl(2) = StartPosition(thisstart,2);         % centre pos, y

% The going through different parameters based on the pre-defined starting
% conditions and bounds is all done within lsqcurvefit. The outputs are the
% parameters of the best fit.
[Allfitparams(thisstart,:), AllR2(thisstart,:)] = lsqcurvefit(@(x,xdata) xyposGaussianInt(x,xinputs),...
    initparamsBl, xinputs, ydata, options.lb, options.ub, fitFuncoptions);
end

[bestfit,bestindx]=min(mean(AllR2,2));

fitparams=Allfitparams(bestindx,:);
R2=AllR2(bestindx,:)/range(forfit_f0(:));


% Get predicted response
% Create distance metric from center pixel
[~,gaussFiltD] = cart2pol(forfit_xval-fitparams(1),forfit_yval-fitparams(2));
% baseline + peak gain * gaussian function. 
modelout = fitparams(5)+ fitparams(3).*exp(-(gaussFiltD.^2./(2*fitparams(4).^2)));

% Get interpolated prediction
[options.interpx,options.interpy] = meshgrid( ...
    linspace(min(xinputs(:,1)),max(xinputs(:,1)),size(forfit_xval,2)*interpModelOutMultiplier), ...
    linspace(min(xinputs(:,2)),max(xinputs(:,2)),size(forfit_yval,1)*interpModelOutMultiplier));
% Create distance metric from center pixel
[~,gaussFiltD] = cart2pol(options.interpx-fitparams(1),options.interpy-fitparams(2));
% baseline + peak gain * gaussian function. 
options.interpmodelout = fitparams(5)+ fitparams(3).*exp(-(gaussFiltD.^2./(2*fitparams(4).^2)));


    function modelout2 = xyposGaussianInt(param,xinputs)
        % Create interpolated distance metric from center pixel
        [~,gaussFiltDInt] = cart2pol(xinputs(:,1) - param(1),xinputs(:,2) - param(2));
        % Create center mechanism (in 2D matrix)
        column_model = param(5)+ param(3).*exp(-(gaussFiltDInt.^2./(2*param(4).^2)));
        % Columnise for lsq to fit properly
        modelout2 = column_model(:);
    end

end


function [col,row] = centreOfMass2DInt(data)
xProj = sum(data,1);
yProj = sum(data,2)';
col = fix(sum(xProj.*(1:length(xProj)))/sum(xProj));
col = min([col,length(xProj)]);
col = max([col,1]);
row = fix(sum(yProj.*(1:length(yProj)))/sum(yProj));
row = min([row,length(yProj)]);
row = max([row,1]);
end