%W3MAXSTAT Conversion between parameters and ~moments 
% of three-parameter maxima Weibull probability distribution (STRUREL)
%
%SYNOPSYS
% a = W3MAXSTAT(meanx_omega, covx_w, skewx_k, def_type)
%
%INPUT
% meanx_omega       - mean value of sample x; if def_type='moment' /scalar/
%                   - omega parameter; if def_type='par' /scalar/
%
% covx_w            - coefficient of variation of sample x; if def_type='moment' /scalar/
%                   - w parameter; if def_type='par' /scalar/
%
% skewx_k           - skewness of sample x (see Notes (1)!); if def_type='moment' /scalar/
%                   - k (shape) parameter; if def_type='par' /scalar/
%
% def_type          - how the distribution function is defined: 'moment' (default) or 'par' /string/
%
%OUTPUT
% a                 - parameters of the distribution [omega, w,    k]; if def_type='moment' /vector/
%                   - ~moments of the distribution   [meanx, covx, skewx]; if def_type='par' /vector/
%
%NOTES:
%  parametrization from STRUREL manual, Appendix A
% (1) F(x) = exp(-((omega-x)/(omega-w))^k)
%
% (2) If the distribution is given by ~moments, then it is recommended to use the following formula:
%   skewx = g_1 * sqrt(n(n-1)) / (n-2)
%
% using:
%   m_r   = sum_i (x_i - mu)^r / n for the sample moments of order r
%   g_1   = m_3 / m_2^(3/2)
%
% This corresponds to matlab built-in skewness(.) function with flag = 0 (bias corrected)!
%
% (3) Formulas of mean, variance and skewness is from:
%  [1] STRUREL 2003. COMREL & SYSREL: Users Manual. Componental and Systel Reliability Anlysis. Using Built-in Symbolic Processor. München, Germany: RCP GmbH.
%  [2] Muraleedharan G. 2013. Characteristic and Moment Generating Functions of Three Parameter Weibull Distribution-an Independent Approach. Research Journal of Mathematical and Statistical Sciences. Vol 1(8). p.25-27.
%
%See also:
% gevstat, gev_stat, w3minstat

function a = w3maxstat(meanx_omega, covx_w, skewx_k, def_type)

if nargin < 4
    def_type = 'moment';
end
if nargin < 3
    error('Too few input arguments!');
end

switch def_type
    case {'mom', 'moment'}

        % following the notations of [1] for convinience
        mu_x    = meanx_omega;
        sigma_x = meanx_omega*covx_w;
        G_x     = skewx_k;
        
        xini    = 2;
        
        %[xsol, fval] = fzero(@(x) skew_fun(x) - G_x, xini);
        [xsol, fval] = fminsearch(@(x) (skew_fun(x) - G_x).^2, xini);
%         disp(fval)
        
        k = xsol;
        
        % check the domain of parameters
        if k < 0
            error('With the given ~moments the ''k'' (shape) parameter is smaller than equal to 0!')
        end
        
        % scale = omega - w
        % shift = omega
        scale = sigma_x/sqrt(gamma(1+2/k) - gamma(1+1/k)^2);
        shift = mu_x + scale*gamma(1+1/k);
        
        omega = shift;
        w     = omega - scale; 
        
        a = [omega, w, k];
        
    case {'par', 'param', 'parameter'}

        omega   = meanx_omega;
        w       = covx_w;
        k       = skewx_k;
        
        % check the domain of parameters
        if k < 0
            error('The ''k'' parameter should be non-negative!')
        end
        if omega < k
            error('Parameter ''omega'' should be smaller than ''k'' parameter!')
        end

        
        mu_x    = omega - (omega-w)*gamma(1+1/k);
        sigma_x = (omega-w)*sqrt(gamma(1+2/k) - gamma(1+1/k)^2);
        
        K       = 1:3;
        g       = gamma(1+K/k);
        G_x     = -((g(3) - 3*g(1)*g(2) + 2*g(1)^3)/(g(2) - g(1)^2)^(3/2));
        
        meanx   = mu_x;
        covx    = sigma_x/meanx;
        skewx   = G_x;
        
        a = [meanx, covx, skewx];
            
    otherwise
        error(['Unknown def_type: ', def_type, ' !'])
end

    function S = skew_fun(x)
        KK = 1:3;
        gg = gamma(1+KK/x);
        S  = -((gg(3) - 3*gg(1)*gg(2) + 2*gg(1)^3)/(gg(2) - gg(1)^2)^(3/2)); 
    end

end