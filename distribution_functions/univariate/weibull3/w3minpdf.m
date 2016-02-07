function p = w3minpdf(x, tau, w, k)
% p = W3MINPDF(x, tau, w, k)
% Probability distribution of three-parameter minima Weibull distribution (STRUREL)

% parameter range check
if any(k < 0)
    error('The ''k'' parameter should be non-negative!')
end
if any(tau > w)
    error('Parameter ''tau'' should be smaller than ''w'' parameter!')
end

idx = ~(x < tau);

p = nan(length(x), 1);

p(~idx) = 0;
p(idx)  = k./(w-tau).*((x(idx)-tau)./(w-tau)).^(k-1).*exp(-((x(idx)-tau)./(w-tau)).^k);

end
