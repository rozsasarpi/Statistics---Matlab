% Bivariate extreme t copula cumulative probability density function (pdf)
% (uniform marginals)
%
%SYNOPSYS
%   c = BITEV_COPULACDF(u, rho, nu)
%
%INPUT
% u         [0,1]x[0,1] /vector; nx2/
% rho       correlation coefficient ]-1,1[ /scalar/
% nu        degree of freedom /scalar/ (default = 2)
%
%OUTPUT
% c         probability density at u /vector; nx1/
%
%SEE ALSO
% bitev_copulapdf, binorm_copulapdf, bit_copulapdf

% taken from R copula package tevCopula.R ptevCopula
function C = bitev_copulacdf(u, rho, nu)

%==========================================================================
% INPUT CHECK & INITIALIZATION
%==========================================================================
if abs(rho) > 1
    error('rho should be the element of [-1,1]!')
end

if any(any(u > 1)) || any(any(u < 0))
    error('u should be the element of [0,1]!')
end

% default degrees of freedom
if nargin < 3
    nu = 2;
end

if nu < 1
    error('nu should be greater than or equal 1!')
end

%==========================================================================
% CALCULATION
%==========================================================================
u1      = u(:,1);
u2      = u(:,2);
uu      = u1.* u2;
r       = uu;
p       = r > 0;
nna     = ~isnan(p);
p       = p & nna; % p: positive uu
logu    = log(uu(p));
r(p)    = exp(logu.* A_tev(log(u2(p))./ logu, rho, nu)) ;
r(~p & nna) = 0; % when one u is zero

C = r;

end