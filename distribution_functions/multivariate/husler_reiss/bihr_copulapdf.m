% Bivariate Hüsler-Reiss copula probability density function (pdf)
% (uniform marginals)
%
%SYNOPSYS
%   c = BIHR_COPULAPDF(u, delta)
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

function c = bihr_copulapdf(u, delta)

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
    c = nan(size(u));
else
    u1      = u(:,1);
    u2      = u(:,2);
    u1p     = -log(u(:,1));
    u2p     = -log(u(:,2));
    z       = u1p./u2p;
    
    % copula cdf
    C = exp(-u1p.* normcdf(1/delta + 0.5 * delta * log(u1p./u2p))...
        -u2p.* normcdf(1/delta + 0.5 * delta * log(u2p./u1p)));
%     keyboard
    % copula pdf
    % taken from R copula package huslerReissCopula.R dhuslerReissCopula
    % function; verified from [Joe (1997). Multivariate Models and Dependence Concepts, p.142]
    c = 1./(u1.* u2).* C.*...
        (normcdf(1/delta - 0.5 * delta * log(z)).*...
         normcdf(1/delta + 0.5 * delta * log(z)) +...
        0.5 * delta./ u2p.* normpdf(1/delta + 0.5 * delta * log(z)));
end
end