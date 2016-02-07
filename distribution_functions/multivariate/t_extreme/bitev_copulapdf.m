% Bivariate extreme t copula probability density function (pdf)
% (uniform marginals)
%
%SYNOPSYS
%   c = BITEV_COPULAPDF(u, rho, nu)
%
%INPUT
% u         [0,1]x[0,1] /vector; nx2/
% rho       [-1,1] /scalar/
% nu        degrees of freedom
%
%OUTPUT
% c         probability density at u /vector; nx1/
%
%SEE ALSO
% binorm_copulapdf, bit_copulapdf

function c = bitev_copulapdf(u, rho, nu)

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

% taken from R copula package tevCopula.R dtevCopula
u1      = u(:,1);
u2      = u(:,2);
C       = bitev_copulacdf(u, rho, nu);
logu    = log(u1.* u2);
w       = log(u2)./ logu;

dwdu1   = -log(u2)./(u1.*log(u1.*u2).^2);
dwdu2   = (log(u1.*u2)-log(u2))./(u2.*log(u1.*u2).^2);
d2wdu1du2 = (2*log(u2)-log(u1.*u2))./(u1.*u2.*log(u1.*u2).^3);

A       = A_tev(w, rho, nu);
Ader    = dAdu_tev(w, rho, nu);

Ader1   = Ader(:,1);
Ader2   = Ader(:,2);
dCdu1   = C.* (1./ u1.* A + logu.* Ader1.* dwdu1); % unused?! same in copula package
dCdu2   = C.* (1./ u2.* A + logu.* Ader1.* dwdu2);

% copula pdf
c = dCdu2.* (1./ u1.* A + logu.* Ader1.* dwdu1) +...
    C.* (1./ u1.* Ader1.* dwdu2 + 1./ u2.* Ader1.* dwdu1 +...
    logu.* Ader2.* dwdu2.* dwdu1 + logu.* Ader1.* d2wdu1du2);

end