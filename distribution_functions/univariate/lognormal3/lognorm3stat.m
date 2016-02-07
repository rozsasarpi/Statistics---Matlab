%LOGNORM3STAT Conversion between parameters and ~moments 
% of three-parameter lognormal probability distribution
%
%SYNOPSYS
% a = LOGNORM3STAT(meanx_shape, covx_scale, skewx_thres, def_type)
%
%INPUT
% meanx_shape       - mean value of sample x; if def_type='moment' /scalar/
%                   - shape parameter; if def_type='par' /scalar/
%
% covx_scale        - coefficient of variation of sample x; if def_type='moment' /scalar/
%                   - scale parameter; if def_type='par' /scalar/
%
% skewx_thres       - skewness of sample x (see Notes (1)!); if def_type='moment' /scalar/
%                   - threshold (shift) parameter; if def_type='par' /scalar/
%
% def_type          - how the distribution function is defined: 'moment' (default) or 'par' /string/
%
%OUTPUT
% a                 - parameters of the distribution [shape, scale, thres]; if def_type='moment' /vector/
%                   - ~moments of the distribution   [meanx, covx, skewx]; if def_type='par' /vector/
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
% lognstat, lognorm3pdf, logn3fit

function a = lognorm3stat(meanx_shape, covx_scale, skewx_thres, def_type)

if nargin < 4
    def_type = 'moment';
end
if nargin < 3
    error('Too few input arguments!');
end

switch def_type
    case {'mom', 'moment'}
        
%         if skewx_thres < 0
%             error('No admissible moment estimate exists for negative skewness! [http://www.jstor.org/stable/2287466]')
%         end
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
        
        % check the domain of parameters
        if shape < 0
            error('With the given ~moments the shape parameter is smaller than equal to 0!')
        end
        
        a = [shape, scale, thres];
        
    case {'par', 'param', 'parameter'}

        shape   = meanx_shape;
        scale   = covx_scale;
        thres   = skewx_thres;
        
        % check the domain of parameters
        if shape < 0
            error('The shape parameter should be non-negative!')
        end
        
        c       = shape^2;
        b       = scale;
        a       = thres;
        
        sigma_y = sqrt(c);
        mu_y    = b;
        w       = exp(sigma_y^2);
        theta   = sqrt(w - 1);
        sigma_x = sqrt(w*(w-1)*exp(2*mu_y));
        mu_x    = a + sigma_x/theta;
        G_x     = theta^3 + 3*theta;
        
        meanx   = mu_x;
        covx    = sigma_x/meanx;
        skewx   = G_x;
        
%         if skewx < 0
%             error('No admissible moment estimate exists for negative skewness! [http://www.jstor.org/stable/2287466]')
%         end
        
        a = [meanx, covx, skewx];
            
    otherwise
        error(['Unknown def_type: ', def_type, ' !'])
end

end