% taken from R copula package tevCopula.R ATev
%
% A = A_TEV(w, rho, nu)
%

function A = A_tev(w, rho, nu)

A   = nan(size(w));
idx = w == 0 | w == 1;

wnu = (w./(1-w)).^(1/nu);
x   = (wnu - rho)./ sqrt(1 - rho.^2).* sqrt(nu + 1);
y   = (1./ wnu - rho)./ sqrt(1 - rho.^2).* sqrt(nu + 1);
A(~idx) = w.* tcdf(x, nu + 1) + (1 - w).* tcdf(y, nu + 1);

end