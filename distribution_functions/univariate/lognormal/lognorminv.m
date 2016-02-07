% Lognormal inverse cumulative distribution functioon
%
% X = LOGNORMINV(P, meanX, covX)
%
%INPUT
% P     - probability where we are interested in the value of X
% meanX - mean value of random variable X
% covX  - coefficient of variation of random variable X
% 
%OUTPUT
% X     - corresponding to P probability

function X = lognorminv(P, meanX, covX)

    mu_lognorm = log(meanX.^2./sqrt((meanX.*covX).^2 + meanX.^2));
    sigma_lognorm = sqrt(log((meanX.*covX).^2./meanX.^2 + 1));

    X = logninv(P, mu_lognorm, sigma_lognorm);
    
end
