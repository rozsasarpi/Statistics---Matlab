function x = w3maxinv(P, omega, w, k)
% x = W3MAXINV(P, omega, w, k)

% parameter range check
if any(k < 0)
    error('The ''k'' parameter should be non-negative!')
end
if any(omega < k)
    error('Parameter ''omega'' should be smaller than ''k'' parameter!')
end
if any(P < 0 | P > 1) 
    error('Probabilities should be in the interval [0,1]!')
end

x = omega - (omega - w).*(-log(P)).^(1/k);

end
