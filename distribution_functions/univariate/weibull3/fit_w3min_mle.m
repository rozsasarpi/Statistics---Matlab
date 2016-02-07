% FIT_W3MIN_MLE Three-parameter minima Weibull, maximum likelihood parameter estimation
% using STRUREL's parametrization
%
%SYNOPSYS
% [parmhat, max_nLL, output] = FIT_W3MIN_MLE(data)
%
%INPUT
% data      sample vector
%
%OUTPUT
% parmhat   MLE estimation of parameters /tau, w, k/
% output    output argument of the optimization algorithm
%
%NOTES:
%  parametrization from STRUREL manual, Appendix A
% (1) F(x) = 1-exp(-((x-tau)/(w-tau))^k)
%
% (2) Reference(s):
%  [1] STRUREL 2003. COMREL & SYSREL: Users Manual. Componental and Systel Reliability Anlysis. Using Built-in Symbolic Processor. München, Germany: RCP GmbH.
%
% See also
% gevfit, gevstat, gev_stat, fit_w3max_mle 

function [parmhat, max_nLL, output] = fit_w3min_mle(data)

% starting point is the method of moment value
mom_par = w3minstat(mean(data), std(data)/mean(data), skewness(data, 0), 'mom');

%WARNING!
mom_par(1) = min(min(data)*0.9, mom_par(1));

% refine options if needed
options = optimset('MaxFunEvals', 8e3, ... 
                   'MaxIter', 4e3, ...
                   'TolFun', abs(nloglike(mom_par))*1e-4, ...
                   'TolX', mean(abs(mom_par))*1e-4);
% options = optimset('MaxFunEvals', 8e3, ... 
%                    'MaxIter', 4e3, ...
%                    'TolFun', max(abs(nloglike(mom_par))*1e-4, 1e-10), ...
%                    'TolX', max(mean(abs(mom_par))*1e-4, 1e-10));

[parmhat, max_nLL, ~, output] = fminsearch(@nloglike, mom_par, options);
               
% negative log-likelihood function
    function nLL = nloglike(param)
        tau = param(1);
        w   = param(2);
        k   = param(3);
        
        nLL = -sum(log(w3minpdf(data, tau, w, k)));
    end

end