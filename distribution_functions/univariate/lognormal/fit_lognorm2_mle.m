% FIT_LOGNORM2_MLE Two-parameter lognormal, maximum likelihood parameter estimation
%
%SYNOPSYS
% [parmhat, output] = FIT_LOGNORM2_MLE(data, def_type)
%
% if def_type == 'par' (default)
% parmhat = [scale_hat (mu), shape_hat (sigma)]
%
% if def_type == 'mom'
% parmhat = [mean_hat, std_hat]
%
% See also
% lognfit, lognormpdf, fit_lognorm3_mle

function [parmhat, max_nLL, output] = fit_lognorm2_mle(data, def_type)

if nargin < 2
    def_type = 'par';
end

% starting point is the method of moment value
meanX       = mean(data);
stdX        = std(data);

switch lower(def_type)
    case {'p', 'par', 'parm', 'parameter', 'parameters'}
        scale       = log(meanX.^2./sqrt((stdX).^2 + meanX.^2));
        shape       = sqrt(log((stdX).^2./meanX.^2 + 1));
        mom_par     = [scale, shape];
        
        % refine options if needed
        options = optimset('MaxFunEvals', 2e3, ...
            'MaxIter', 2e3, ...
            'TolFun', abs(nloglike_par(mom_par))*1e-8, ...
            'TolX', min(abs(mom_par))*1e-8);
       
        [parmhat, max_nLL, ~, output] = fminsearch(@nloglike_par, mom_par, options);
        
    case {'m', 'mom', 'moment', 'moments'}
        mom_par     = [meanX, stdX];
        
        % refine options if needed
        options = optimset('MaxFunEvals', 4e3, ...
            'MaxIter', 2e3, ...
            'TolFun', abs(nloglike_mom(mom_par))*1e-8, ...
            'TolX', min(abs(mom_par))*1e-8);
        
        [parmhat, max_nLL, ~, output] = fminsearch(@nloglike_mom, mom_par, options);
        
    otherwise
        error(['Unknown def_type: ', def_type])
end

% negative log-likelihood function
    function nLL = nloglike_par(param)
        scale_ = param(1);
        shape_ = param(2);
        
        nLL = -sum(log(lognpdf(data, scale_, shape_)));
    end

    function nLL = nloglike_mom(param)
        mean_ = param(1);
        std_  = param(2);
        nLL = -sum(log(lognormpdf(data, mean_, std_/mean_)));
    end
end