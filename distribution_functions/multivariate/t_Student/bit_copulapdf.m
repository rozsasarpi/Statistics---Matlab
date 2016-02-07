% Bivariate t (Student) copula probability density function (pdf)
% (uniform marginals)
%
%SYNOPSYS
%   c = BIT_COPULAPDF(u, rho, nu)
%
%INPUT
% u         [0,1]x[0,1] /vector; nx2/
% rho       correlation coefficient ]-1,1[ /scalar/
% nu        degree of freedom /scalar/ (default = 2)
%
%OUTPUT
% c         probability density at u /vector; nx1/
%
%NOTES
% faster than the general, multidimensional matlab built-in copulapdf
%
%SEE ALSO
% bihr_copulapdf, binorm_copulapdf

function c = bit_copulapdf(u, rho, nu)

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
tinv_u  = tinv(u,nu);
tinv_u1 = tinv_u(:,1);
tinv_u2 = tinv_u(:,2);

c = bitpdf([tinv_u1, tinv_u2], rho, nu)./(tpdf(tinv_u1,nu).*tpdf(tinv_u2,nu));

    function Y = bitpdf(x, rho, nu)
        x1 = x(:,1);
        x2 = x(:,2);
        Y = gamma((nu+2)/2)/(gamma(nu/2)*pi*nu*sqrt(1-rho^2)).*...
            (1 + (x1.^2 + x2.^2 - 2*rho*x1.*x2)/(nu*(1-rho^2))).^(-(nu+2)/2);
        
    end
end