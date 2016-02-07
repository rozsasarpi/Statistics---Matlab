% FIT_WEIBULL3MIN_MLE Three-parameter weibull maximum likelihood parameter estimation
% using GEV parametrization with constrain on shape parameter (< 0)
%
%SYNOPSYS
% [parmhat, max_nLL, output] = FIT_WEIBULLM3MIN_MLE(data)
%
%INPUT
% data      sample vector
%
%OUTPUT
% parmhat   MLE estimation of parameters /shape, scale, threshold/
% output    output argument of the optimization algorithm
%
% See also
% gevfit, gevstat, gev_stat 

function [parmhat, max_nLL, output] = fit_weibull3min_mle(data)

% starting point is the method of moment value
mom_par = gev_stat(mean(data), std(data)/mean(data), skewness(data, 0), 'mom');
mom_par(1) = min(-1e-4, mom_par(1));

loc_ub = (max(data) - mom_par(2)/-mom_par(1))*0.9;
mom_par(3) = max(loc_ub, mom_par(3));
% if isfinite(nloglike(mom_par))
%     error('..')
% end
% refine options if needed
% options = optimset('MaxFunEvals', 8e3, ... 
%                    'MaxIter', 4e3, ...
%                    'TolFun', max(abs(nloglike(mom_par))*1e-4, 1e-10), ...
%                    'TolX', max(mean(abs(mom_par))*1e-4), 1e-10);
options = optimoptions('fmincon',...
                   'MaxFunEvals', 8e3, ... 
                   'MaxIter', 4e3, ...
                   'TolFun', max(abs(nloglike(mom_par))*1e-4, 1e-10), ...
                   'TolX', max(mean(abs(mom_par))*1e-4, 1e-10), 'Display', 'notify');

lb = [-Inf, 1e-6, -Inf];
ub = [-1e-6, Inf, Inf];
[parmhat, max_nLL, ~, output] = fmincon(@nloglike, mom_par, [],[],[],[], lb, ub, [], options);

% negative log-likelihood function
    function nLL = nloglike(param)
        shape = param(1);
        scale = param(2);
        location = param(3);
        
        nLL = -sum(log(gevpdf(data, shape, scale, location)));
    end

end