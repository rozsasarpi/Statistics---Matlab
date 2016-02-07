% Bivariate Clayton copula cumulative distribution function (cdf)
% (uniform marginals)
%
%SYNOPSYS
%   C = BICLAY_COPULACDF(u, theta)
%
%INPUT
% u         [0,1]x[0,1] /vector; nx2/
% theta     copula parameter [-1,Inf]/{0} /scalar/
%
%OUTPUT
% C         probability at u /vector; nx1/
%
%NOTES
% faster than the general, multidimensional matlab built-in copulapdf
% theta and Kendall tau connection:
%   theta = 2*k_tau/(1-k_tau);
%
%SEE ALSO
% biclay_copulapdf

function C = biclay_copulacdf(u, theta)

%==========================================================================
% INPUT CHECK & INITIALIZATION
%==========================================================================
if theta < -1 || theta == 0
    error('theta should be the element of [-1,Inf]/{0}!')
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

C = (u1.^-theta + u2.^-theta -1).^-(1/theta);
C = arrayfun(@(y) max([y, 0]), C);

end