% taken from R copula package tevCopula.R dAduTev
%
% Ader = DADU_TEV(w, rho, nu)

function Ader = dAdu_tev(w, rho, nu)

wnu = (w./(1-w)).^(1/nu);
x   = (wnu - rho)./ sqrt(1 - rho.^2).* sqrt(nu + 1);
y   = (1./ wnu - rho)./ sqrt(1 - rho.^2).* sqrt(nu + 1);

% prepare dx, dy
dxdw = (sqrt(nu+1).*(w./(1-w).^2+1./(1-w)).*...
    (w./(1-w)).^(1./nu-1))./(nu.*sqrt(1-rho.^2));
d2xdw2 = ((1./nu-1).*sqrt(nu+1).*(w./(1-w).^2+1./(1-w)).^2.*...
    (w./(1-w)).^(1./nu-2))./(nu.*sqrt(1-rho.^2))+...
    (sqrt(nu+1).*((2*w)./(1-w).^3+2./(1-w).^2).*...
    (w./(1-w)).^(1./nu-1))./(nu*sqrt(1-rho.^2));

dydw = -(sqrt(nu+1).*(w./(1-w).^2+1./(1-w)).*...
    (w./(1-w)).^(-1./nu-1))./(nu.*sqrt(1-rho.^2));
d2ydw2 = -((-1./nu-1).*sqrt(nu+1).*(w./(1-w).^2+1./(1-w)).^2.*...
    (w./(1-w)).^(-1./nu-2))./(nu.*sqrt(1-rho.^2))-...
    (sqrt(nu+1).*((2*w)./(1-w).^3+2./(1-w).^2).*...
    (w./(1-w)).^(-1./nu-1))./(nu.*sqrt(1-rho.^2));

% prepare ddtx, ddty, derivative of the t-(nu) density
ddtx    = ddens(x, nu + 1);
ddty    = ddens(y, nu + 1);

% now collect the results
der1 = tcdf(x, nu + 1) + w.* tpdf(x, nu + 1).* dxdw -...
    tcdf(y, nu + 1) + (1 - w).* tpdf(y, nu + 1).* dydw;
der2 = tpdf(x, nu + 1).* dxdw +...
    tpdf(x, nu + 1).* dxdw + w.* ddtx.* dxdw.^2 + w.* tpdf(x, nu + 1).* d2xdw2 +...
    (- tpdf(y, nu + 1).* dydw ) +...
    (- tpdf(y, nu + 1).* dydw) + (1 - w).* ddty.* dydw.^2 + (1 - w).* tpdf(y, nu + 1).* d2ydw2;

Ader = [der1, der2];

% Nested function
    function d = ddens(u, nu)
        d = gamma(0.5 * (nu + 1)) / gamma(0.5 * nu) / sqrt(nu * pi)*...
            -((nu+1).*u.*((nu+u.^2)/nu).^(-0.5*nu-1.5))/nu;
    end
end