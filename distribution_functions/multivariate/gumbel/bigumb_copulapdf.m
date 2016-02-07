% Bivariate Gumbel copula probability density function (pdf)
% (uniform marginals)
%
%SYNOPSYS
%   c = BIGUMB_COPULAPDF(u, theta)
%
%INPUT
% u         [0,1]x[0,1] /vector; nx2/
% theta       correlation coefficient [1,Inf] /scalar/
%
%OUTPUT
% c         probability density at u /vector; nx1/
%
%NOTES
% faster than the general, multidimensional matlab built-in copulapdf
% however for highly correlated rvs numerical problems arise, 0/0 -> NaN
%
%SEE ALSO
% bihr_copulapdf, bit_copulapdf, bnorm_copulapdf

function c = bigumb_copulapdf(u, theta)

%==========================================================================
% INPUT CHECK & INITIALIZATION
%==========================================================================
if theta < 1 
    error('theta should be the element of [1,Inf]!')
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
u1  = u(:,1);
u2  = u(:,2);
u1p = -log(u1);
u2p = -log(u2);

% Joe. 1997: Multivariate Models and Dependence Concepts. p.142
% cdf
C = exp(-(u1p.^theta + u2p.^theta).^(1/theta));

%pdf
c = C.*(u1.*u2).^(-1).*(u1p.*u2p).^(theta-1)./...
    ((u1p.^theta + u2p.^theta).^(2-1/theta)).*...
    ((u1p.^theta + u2p.^theta).^(1/theta) + theta - 1);
end