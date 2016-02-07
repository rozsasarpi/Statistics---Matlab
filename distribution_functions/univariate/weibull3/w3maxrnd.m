function R = w3maxrnd(omega, w, k, m, n)
% R = W3MAXRND(omega, w, k, m, n)
% Random three-parameter maxima Weibull distributed numbers (STRUREL)

if nargin < 5
    n = 1;
end
if nargin < 4
    m = 1;
else

% parameter range check
if any(k < 0)
    error('The ''k'' parameter should be non-negative!')
end
if any(omega < k)
    error('Parameter ''omega'' should be smaller than ''k'' parameter!')
end

RU = rand(m, n);
R = w3maxinv(RU, omega, w, k);

end
