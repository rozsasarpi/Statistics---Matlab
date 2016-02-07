% FIT_W3MAX_MLE Three-parameter Weibull maximum likelihood parameter estimation
% using STRUREL's parametrization
%
%SYNOPSYS
% [parmhat, max_nLL, output] = FIT_W3MAX_MLE(data)
%
%INPUT
% data      sample vector
%
%OUTPUT
% parmhat   MLE estimation of parameters /omega, w, k/
% output    output argument of the optimization algorithm
%
%NOTES:
%  parametrization from STRUREL manual, Appendix A
% (1) F(x) = exp(-((omega-x)/(omega-w))^k)
%
% (2) Reference(s):
%  [1] STRUREL 2003. COMREL & SYSREL: Users Manual. Componental and Systel Reliability Anlysis. Using Built-in Symbolic Processor. München, Germany: RCP GmbH.
%
% See also
% gevfit, gevstat, gev_stat, fit_w3min_mle 

function [parmhat, max_nLL, output] = fit_w3max_mle(data)

% starting point is the method of moment value
mom_par = w3maxstat(mean(data), std(data)/mean(data), skewness(data, 0), 'mom');

mom_par(1) = max(max(data)*1.1, mom_par(1));

% refine options if needed
% options = optimset('MaxFunEvals', 8e3, ... 
%                    'MaxIter', 4e3, ...
%                    'TolFun', abs(nloglike(mom_par))*1e-4, ...
%                    'TolX', mean(abs(mom_par))*1e-4);
% % options = optimset('MaxFunEvals', 8e3, ... 
% %                    'MaxIter', 4e3, ...
% %                    'TolFun', max(abs(nloglike(mom_par))*1e-4, 1e-10), ...
% %                    'TolX', max(mean(abs(mom_par))*1e-4, 1e-10));    

% [parmhat, max_nLL, ~, output] = fminsearch(@nloglike, mom_par, options);

% WARNING - the constrained optimization often finds the initial point as optimum solution
% sometimes it is considerably different than that obtained by the unconstrained optimization
options = optimoptions('fmincon',...
                        'Algorithm', 'sqp',...
                        'MaxFunEvals', 8e3, ... 
                        'MaxIter', 4e3, ...
                        'TolFun', abs(nloglike(mom_par))*1e-6, ...
                        'TolX', mean(abs(mom_par))*1e-6, 'Display', 'notify'); 

A = [-1, 1, 0;
      0, 1, 0;
      0, 0, -1];
b = [0, Inf, 0];
lb = [max(data)*(1+sqrt(eps)), -Inf, sqrt(eps)];
ub = [Inf, Inf, Inf];
[parmhat, max_nLL, ~, output] = fmincon(@nloglike, mom_par,A,b,[],[],lb,ub,[],options);
               
% negative log-likelihood function
    function nLL = nloglike(param)
        omega   = param(1);
        w       = param(2);
        k       = param(3);
        
        nLL = -sum(log(w3maxpdf(data, omega, w, k)));
    end

end