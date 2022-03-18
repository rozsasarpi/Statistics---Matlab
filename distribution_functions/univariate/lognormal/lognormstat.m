%LOGNORMSTAT Conversion between parameters and ~moments
% of two-parameter lognormal probability distribution
%
%SYNOPSYS
% a = LOGNORMSTAT(meanx_shape, covx_scale, def_type)
%
%INPUT
% meanx_shape       - mean value of sample x; if def_type='moment' /scalar/
%                   - shape parameter; if def_type='par' /scalar/
%
% covx_scale        - coefficient of variation of sample x; if def_type='moment' /scalar/
%                   - scale parameter; if def_type='par' /scalar/
%
%
% def_type          - how the distribution function is defined: 'moment' (default) or 'par' /string/
%
%OUTPUT
% a                 - parameters of the distribution [shape, scale]; if def_type='moment' /vector/
%                   - ~moments of the distribution   [meanx, covx]; if def_type='par' /vector/
%
%NOTES:
% X = exp(Y)        - X is lognormally distributed, Y is normally distrib.
% shape             - mean of Y
% scale             - standard deviation of Y
%
%See also:
% lognstat, lognorm3stat

function [meanx_shape, covx_scale] = lognormstat(meanx_shape, covx_scale, def_type)

if nargin < 3
    def_type = 'moment';
end
if nargin < 2
    error('Too few input arguments!');
end

switch def_type
    case {'mom', 'moment'}
        mu_x    = meanx_shape;
        sigma_x = meanx_shape.*covx_scale;
        
        shape   = log(mu_x.^2./sqrt(sigma_x.^2 + mu_x.^2));
        scale   = sqrt(log(sigma_x.^2./mu_x.^2 + 1));
        
        meanx_shape = shape;
        covx_scale = scale;
    case {'par', 'param', 'parameter'}
        shape   = meanx_shape;
        scale   = covx_scale;
        
        meanx   = exp(shape + scale.^2/2);
        covx    = sqrt(exp(scale.^2) - 1);
        
        meanx_shape = meanx;
        covx_scale = covx;
    otherwise
        error(['Unknown def_type: ', def_type, ' !'])
end

end