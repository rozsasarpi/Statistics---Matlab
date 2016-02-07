function x = w3mininv(P, tau, w, k)
% x = W3MININV(P, tau, w, k)

% parameter range check
if any(k < 0)
    error('The ''k'' parameter should be non-negative!')
end
if any(tau > w)
    error('Parameter ''tau'' should be smaller than ''w'' parameter!')
end
if any(P < 0 | P > 1) 
    error('Probabilities should be in the interval [0,1]!')
end

x = tau + (w - tau).*(-log(1-P)).^(1/k);

end
