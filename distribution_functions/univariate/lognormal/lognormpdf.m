% Lognormal probability density functioon
%
% P = LOGNORMPDF(X, meanX, covX)
%
%INPUT
% X     - value where we are interested in the probability density
% meanX - mean value of random variable X
% covX  - coefficient of variation of random variable X
% 
%OUTPUT
% P     - probability corresponding to value X

function P = lognormpdf(X, meanX, covX)

    mu_lognorm = log(meanX.^2./sqrt((meanX.*covX).^2 + meanX.^2));
    sigma_lognorm = sqrt(log((meanX.*covX).^2./meanX.^2 + 1));

    P = lognpdf(X, mu_lognorm, sigma_lognorm);
    
end