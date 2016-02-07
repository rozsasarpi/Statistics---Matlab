% FIT_LOGNORM3_MLE Three-parameter lognormal maximum likelihood parameter estimation
%
%SYNOPSYS
% [parmhat, max_nLL, output] = FIT_LOGNORM3_MLE(data)
%
%INPUT
% data      sample vector
%
%OUTPUT
% parmhat   MLE estimation of parameters /shape, scale, threshold/
% output    output argument of the optimization algorithm
%
% See also
% logn3fit, lognorm3pdf, lognorm3inv, lognorm3stat 

function [parmhat, max_nLL, output] = fit_lognorm3_mle(data)

% starting point is the method of moment value
mom_par = lognorm3stat(mean(data), std(data)/mean(data), skewness(data, 0), 'mom');

if ~isfinite(nloglike(mom_par))
    
    scale_mom = log(mean(data).^2./sqrt(std(data).^2 + mean(data).^2));
    shape_mom = sqrt(log(std(data).^2./mean(data).^2 + 1));
    
    mom_par = [shape_mom, scale_mom, 0];
end

% refine options if needed
options = optimset('MaxFunEvals', 8e3, ... 
                   'MaxIter', 4e3, ...
                   'TolFun', abs(nloglike(mom_par))*1e-4, ...
                   'TolX', mean(abs(mom_par))*1e-4);
               
[parmhat, max_nLL, ~, output] = fminsearch(@nloglike, mom_par, options);

% negative log-likelihood function
    function nLL = nloglike(param)
        shape = param(1);
        scale = param(2);
        thres = param(3);
        
        nLL = -sum(log(lognorm3pdf(data, shape, scale, thres, 'par')));
    end

end