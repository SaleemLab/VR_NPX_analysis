function [fitparams, modelout] = fit2DGauss(forfit_f0,forfit_xval,forfit_yval)

% Inputs:
% forfit_f0: your ratemap; i.e. a 2d matrix 

[col,row] = centreOfMass2DInt(forfit_f0);

% For Gaussian fit
initparamsBl(1) = forfit_xval(col,row);         % centre pos, x
initparamsBl(2) = forfit_yval(col,row);           % centre pos, y
initparamsBl(3) = forfit_f0(row,col);              % peak gain (gain of centre)
initparamsBl(4) = 5;                            % centre size half of best diameter
initparamsBl(5) = min(forfit_f0(:));        % background activity

% Set up fit including bounds
% Generate a vector for the lower and upper bounds (lb and ub respectively)
% of different parameters. [ centreX,centreY,gain of centre,centre size, background activity]
lb = [min(forfit_xval(:))-initparamsBl(4), min(forfit_yval(:))-initparamsBl(4),...
    0, 1, initparamsBl(5)/2+0.01]; %
ub = [max(forfit_xval(:))+initparamsBl(4), max(forfit_yval(:))+initparamsBl(4),...
    1000, 100, initparamsBl(5)*2+0.01]; % Just centre, no surround
options = optimset('TolFun',1e-5,'TolX',1e-4,'Maxiter',100000,'MaxFunEvals',100000, ...
    'Display','off');

% Fit
% x = initparamsBl;
xinputs = [forfit_xval(:) forfit_yval(:)];
ydata = forfit_f0(:);
% The going through different parameters based on the pre-defined starting
% conditions and bounds is all done within lsqcurvefit. The outputs are the
% parameters of the best fit.
fitparams = lsqcurvefit(@(x,xdata) xyposGaussianInt(x,xinputs),...
    initparamsBl, xinputs, ydata, lb, ub, options);

% Get predicted response
% Create interpolated distance metric from center pixel
[tt,gaussFiltD] = cart2pol([xinputs(:,2) - fitparams(1)],[xinputs(:,1) - fitparams(2)]);

% Create center mechanism
% baseline + peak gain * gaussian function. It doesn't seem like there is
% any opportunity for the centre to be offset in this line
modelout = fitparams(5)+ fitparams(3).*exp(-(gaussFiltD./(fitparams(4))).^2);
% the below line must be correct as the matrix comes out the right size
modelout = reshape(modelout,size(forfit_xval,1),size(forfit_xval,2));

    function [modelout,modelout_full] = xyposGaussianInt(param,xinputs)
        % Create interpolated distance metric from center pixel
        [tt,gaussFiltD] = cart2pol([xinputs(:,1) - param(1)],[xinputs(:,2) - param(2)]);
        % Create center mechanism
        column_model = param(5)+ param(3).*exp(-(gaussFiltD./(param(4))).^2);
        % Retrieve points that should align with data comparisons
        modelout = column_model(:);
    end

end


function [col,row] = centreOfMass2DInt(data)
xProj = sum(data,1);
yProj = sum(data,2)';
col = fix(sum(xProj.*(1:length(xProj)))/sum(xProj));
row = fix(sum(yProj.*(1:length(yProj)))/sum(yProj));
end