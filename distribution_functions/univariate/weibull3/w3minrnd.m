function R = w3minrnd(tau, w, k, m, n)
% R = W3MINRND(tau, w, k, m, n)
% Random three-parameter minima Weibull distributed numbers (STRUREL)

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
if any(tau > w)
    error('Parameter ''tau'' should be smaller than ''w'' parameter!')
end

RU = rand(m, n);
R = w3mininv(RU, tau, w, k);

end
