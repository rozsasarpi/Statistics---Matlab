function P = w3mincdf(x, tau, w, k)
% P = W3MAXCDF(x, tau, w, k)
% Cumulative probability distribution of three-parameter minima Weibull distribution (STRUREL)

% parameter range check
if any(k < 0)
    error('The ''k'' parameter should be non-negative!')
end
if any(tau > w)
    error('Parameter ''tau'' should be smaller than ''w'' parameter!')
end

idx = ~(x < tau);

P = nan(length(x), 1);

P(~idx) = 0;

P(idx) = 1-exp(-((x-tau)./(w-tau)).^k);

end
