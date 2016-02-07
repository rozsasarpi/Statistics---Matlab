% Bivariate Guass copula probability density function (pdf)
% (uniform marginals)
%
%SYNOPSYS
%   c = BINORM_COPULAPDF(u, rho)
%
%INPUT
% u         [0,1]x[0,1] /vector; nx2/
% rho       correlation coefficient ]-1,1[ /scalar/
%
%OUTPUT
% c         probability density at u /vector; nx1/
%
%NOTES
% faster than the general, multidimensional matlab built-in copulapdf
%
%SEE ALSO
% bihr_copulapdf, bit_copulapdf

function c = binorm_copulapdf(u, rho)

%==========================================================================
% INPUT CHECK & INITIALIZATION
%==========================================================================
if abs(rho) >= 1
    error('rho should be the element of ]-1,1[!')
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
u1 = u(:,1);
u2 = u(:,2);

c = 1/sqrt(1 - rho^2) * exp((2*rho*norminv(u1).*norminv(u2) -...
    rho^2*(norminv(u1).^2 + norminv(u2).^2))/(2*(1-rho^2)));

end