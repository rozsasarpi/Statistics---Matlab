%LOGNORM3INV Three-parameter lognormal probability density function
%
%SYNOPSYS
% x = LOGNORM3INV(P, meanx_shape, covx_scale, skewx_thres, def_type)
%
%INPUT
% P                 - P where we are interested in the value of random variable /vector/
%
% meanx_shape       - mean value of sample x; if def_type='moment' /scalar/
%                   - shape parameter; if def_type='par' /scalar/
%                     scale parameter (sigma) of the corresponding normal distribution
%
% covx_scale        - coefficient of variation of sample x; if def_type='moment' /scalar/
%                   - scale parameter; if def_type='par' /scalar/
%                     location parameter (mu) of the corresponding normal distribution
%
% skewx_thres       - skewness of sample x (see Notes (1)!); if def_type='moment' /scalar/
%                   - threshold (shift) parameter; if def_type='par' /scalar/
%
% def_type          - how the distribution function is defined: 'moment' (default) or 'par' /string/
%
%OUTPUT
% x                 - x value corresponding to P /vector/
%
%NOTES:
%
% (1) If the distribution is given by ~moments, then the skewness is defined as followed:
%   skewx = g_1 * sqrt(n(n-1)) / (n-2)
%
% using:
%   m_r   = sum_i (x_i - mu)^r / n for the sample moments of order r
%   g_1   = m_3 / m_2^(3/2)
%
% This corresponds to matlab built-in skewness(.) function with flag = 0 (bias corrected)!
%
% (2) Formulas of mean, variance and skewness is from:
%  [1] Singh. V.P. (1998). Entropy-Based Parameter Estimation in Hydrology
%
%See also:
% lognorm3pdf, logn3fit, skewness

function x = lognorm3inv(P, meanx_shape, covx_scale, skewx_thres, def_type)

if nargin < 5
    def_type = 'moment';
end
if nargin < 4
    error('Too few input arguments!');
end

switch def_type
    case {'mom', 'moment'}
        % following the notations of [1] for convinience
        mu_x    = meanx_shape;
        sigma_x = meanx_shape*covx_scale;
        G_x     = skewx_thres;
        
        B       = 1/2*(-G_x + (G_x^2 + 4)^(1/2));
        theta   = (1 - B^(2/3))/B^(1/3);
        w       = theta^2 + 1;
        sigma_y = sqrt(log(w));
        mu_y    = 1/2*log(sigma_x^2/(w*(w-1)));
        
        a       = mu_x - sigma_x/theta;
        b       = mu_y;
        c       = sigma_y^2;      
        shape   = sqrt(c);
        scale   = b;
        thres   = a;

    case {'par', 'param', 'parameter'}
        shape   = meanx_shape;
        scale   = covx_scale;
        thres   = skewx_thres;
    otherwise
        error(['Unknow def_type: ', def_type, ' !'])
end
% keyboard
% check the domain of parameters
if shape < 0
    error('The shape parameter should be non-negative!')
end

if any(P < 0 | P > 1) 
    error('Probabilities should be in the interval [0,1]!')
end

x       = zeros(size(P));
idx     = P == 1;
x(idx)  = Inf;

P(idx) = [];

% inverse cdf function
x(~idx) = logninv(P, scale, shape) + thres;

end