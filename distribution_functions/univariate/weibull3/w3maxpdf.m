function p = w3maxpdf(x, omega, w, k)
% p = W3MAXPDF(x, omega, w, k)
% Probability distribution of three-parameter maxima Weibull distribution (STRUREL)

% parameter range check
if any(k < 0)
    error('The ''k'' parameter should be non-negative!')
end
if any(omega < k)
    error('Parameter ''omega'' should be smaller than ''k'' parameter!')
end

idx = x < omega;

p = nan(length(x), 1);

p(~idx) = 0;
p(idx)  = k./(omega-w).*((omega-x(idx))./(omega-w)).^(k-1).*exp(-((omega-x(idx))./(omega-w)).^k);

end
