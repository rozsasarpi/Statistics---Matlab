%GEV_STAT Conversion between parameters and ~moments
% of three-parameter lognormal probability distribution
%
%SYNOPSYS
% a = GEV_STAT(meanx_shape, covx_scale, skewx_loc, def_type)
%
%INPUT
% meanx_shape       - mean value of sample x; if def_type='moment' /scalar/
%                   - shape parameter; if def_type='par' /scalar/
%
% covx_scale        - coefficient of variation of sample x; if def_type='moment' /scalar/
%                   - scale parameter; if def_type='par' /scalar/
%
% skewx_location    - skewness of sample x (see Notes (1)~); if def_type='moment' /scalar/
%                   - threshold (shift) parameter; if def_type='par' /scalar/
%
% def_type          - how the distribution function is defined: 'moment' (default) or 'par' /string/
%
%OUTPUT
% A                 - parameters of the distribution [shape, scale, thres]; if def_type='moment' /vector/
%                   - ~moments of the distribution   [meanx, covx, skewx]; if def_type='par' /vector/
%
%NOTES:
% (1) If the distribution is given by ~moments, then the skewness is defined as followed:
%   skewx = g_1 * sqrt(n(n-1)) / (n-2)
%
% using:
%   m_r   = sum_i (x_i - location)^r / n for the sample moments of order r
%   g_1   = m_3 / m_2^(3/2)
%
% This corresponds to matlab built-in skewness(.) function with flag = 0 (bias corrected)~
%
% (2) def = 'par'; Notation follows that of Wikipedia and Coles(2001)
%     formulas are also from wikipedia (October 2014)
% (3) def = 'par'; Accept vectorized input as well
% (4) def = 'mom'; formulas of mean, variance and skewness is from [1998_Entropy-Based Parameter Estimation in Hydrology]
%     following the notation of the referenced book
%     connection to wikipedia notation:
%     sigma = a (scale)
%     ksi = -b (shape)
%     mu = c (location)
%
%See also:
% gevstat

function A = gev_stat(meanx_shape, covx_scale, skewx_location, def_type)

if nargin < 4
    def_type = 'moment';
end
if nargin < 3
    error('Too few input arguments~');
end

% Check input consistency
lengths     = [length(meanx_shape), length(covx_scale), length(skewx_location)];
length_diff = diff(lengths);
flag        = any(length_diff ~= 0);
if flag == 1
    error('The input parameters must have the same size!')
end

switch def_type
    case {'mom', 'moment'}
        meanx   = meanx_shape;
        covx    = covx_scale;
        skewx   = skewx_location;
        
        varx    = (covx*meanx)^2;
        
        if skewx > 1.139547 %Fréchet
            
            % numerical solution for parameter 'b'
            [b, fval] = fminsearch(@fun_f, -0.1); % a<b !!
            
            a = sqrt(varx*b^2/(gamma(1 + 2*b) - gamma(1 + b)^2));
            c = meanx - a/b*(1 - gamma(1 + b));
            
            A = [-b, a, c];
            
        else % Weibull distribution
            
            % following the wikipedia article notation (same as Coles.2001.Introduction the extreme value theory)
            % the notations should be unified..
            % numerical solution for parameter 'xi'
            [xi, fval] = fminsearch(@fun_w, -0.1);
%             disp(fval)
            sigma = sqrt(varx*xi^2/(gamma(1 - 2*xi) - gamma(1 - xi)^2));
            mu    = meanx - sigma/xi*(gamma(1 - xi) - 1);
            
            A = [xi, sigma, mu];
        end
        
        
    case {'par', 'param', 'parameter'}
        shape   = meanx_shape;
        scale   = covx_scale;
        location= skewx_location;
        
        n = length(shape);
        M = zeros(n, 1);
        V = M;
        S = V;
        
        index_0       = shape == 0;
        index_lt_0    = shape < 0;
        index_gt_0    = shape > 0;
        index_lt_05   = logical((shape < 0.5).*(shape ~= 0));
        index_gte_05  = shape >= 0.5;
        index_lt_1    = logical((shape < 1).*(shape ~= 0));
        index_gte_1   = shape >= 1;
        
        g1 = gamma(1 - 1*shape);
        g2 = gamma(1 - 2*shape);
        g3 = gamma(1 - 3*shape);
        
        % mean
        i    = index_0;
        M(i) = location(i) + scale(i)*0.5772157;
        i    = index_lt_1;
        M(i) = location(i) + scale(i).*(gamma(1-shape(i))-1)./shape(i);
        i    = index_gte_1;
        M(i) = Inf;
        
        % variance
        i    = index_0;
        V(i) = scale(i).^2*pi^2/6;
        i    = index_lt_05;
        V(i) = scale(i).^2.*(g2(i) - g1(i).^2)./shape(i).^2;
        i    = index_gte_05;
        V(i) = Inf;
        
        % skewness
        i    = index_0;
        S(i) = 1.139547099404649;
        i    = index_lt_0;
        S(i) = -(g3(i) - 3*g1(i).*g2(i) + 2*g1(i).^3)./(g2(i) - g1(i).^2).^(3/2);
        i    = index_gt_0;
        S(i) = (g3(i) - 3*g1(i).*g2(i) + 2*g1(i).^3)./(g2(i) - g1(i).^2).^(3/2);
        
        cov_ = sqrt(V)./M;
        
        A    = [M, cov_, S];
    otherwise
        error(['Unknown def_type: ', def_type, ' ~'])
end

    function f = fun_f(b)
        M_2_ = gamma(1 + 2*b) - gamma(1 + b).^2;
        M_3_ = -gamma(1 + 3*b) + 3*gamma(1 + b).*gamma(1 + 2*b) - 2*gamma(1 + b).^3;
        
        skew_ = -M_3_./M_2_.^(3/2); % WARNING, - sign added, probably a typo in the book
        
        f = (skewx - skew_).^2;
    end

    function f = fun_w(xi)
        % following the wikipedia article notation (same as Coles.2001.Introduction the extreme value theory)
        % the notations should be unified..
        % numerical solution for parameter 'xi'
        M_2_ = gamma(1 - 2*xi) - gamma(1 - xi).^2;
        M_3_ = -gamma(1 - 3*xi) + 3*gamma(1 - xi).*gamma(1 - 2*xi) - 2*gamma(1 - xi).^3;
        
        skew_ = M_3_./M_2_.^(3/2);
        
        f = (skewx - skew_).^2;
        
    end
end