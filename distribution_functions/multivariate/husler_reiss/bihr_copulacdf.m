% Bivariate Hüsler-Reiss copula cumulative distribution function (pdf)
% (uniform marginals)
%
%SYNOPSYS
%   c = BIHR_COPULACDF(u, delta)
%
%INPUT
% u         [0,1]x[0,1] /vector; nx2/
% delta     [0,Inf[ /scalar/
%
%OUTPUT
% c         probability density at u /vector; nx1/
%
%SEE ALSO
% binorm_copulapdf, bit_copulapdf

function C = bihr_copulacdf(u, delta)

%==========================================================================
% INPUT CHECK & INITIALIZATION
%==========================================================================
if delta < 0
    error('delta should be the element of [0,Inf]!')
end

if any(any(u > 1)) || any(any(u < 0))
    error('u should be the element of [0,1]!')
end

if size(u,2) ~=2
    error('Wrong input format! u should be an nx2 matrix!')
end

%==========================================================================
% CALCULATION
%==========================================================================
if isinf(delta)
    C = nan(size(u));
else
    u1      = u(:,1);
    u2      = u(:,2);
    u1p     = -log(u1);
    u2p     = -log(u2);
    
    % copula cdf
    % taken from R copula package huslerReissCopula.R phuslerReissCopula
%     C = exp(-u1p.* normcdf(1/delta + 0.5 * delta * log(u1p./u2p))...
%         -u2p.* normcdf(1/delta + 0.5 * delta * log(u2p./u1p)));
    
end
end