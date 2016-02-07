function P = w3maxcdf(x, omega, w, k)
% P = W3MAXCDF(x, omega, w, k)

% parameter range check
if any(k < 0)
    error('The ''k'' parameter should be non-negative!')
end
if any(omega < k)
    error('Parameter ''omega'' should be smaller than ''k'' parameter!')
end

idx = x < omega;

P = nan(length(x), 1);

P(~idx) = 1;

P(idx) = exp(-((omega-x(idx))./(omega-w)).^k);

end
