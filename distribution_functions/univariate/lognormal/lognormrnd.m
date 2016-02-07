% Lognormaly distributed random numbers
%
% R = LOGNORMRND(meanX, covX, m, n)
%
%INPUT
% meanX - mean value of random variable X
% covX  - coefficient of variation of random variable X
% mxn     - number of random values (rows, columns)
% 
%OUTPUT
% R     - lognormaly distributed random variables /vector/

function R = lognormrnd(meanX, covX, m, n)

    mu_lognorm = log(meanX.^2./sqrt((meanX.*covX).^2 + meanX.^2));
    sigma_lognorm = sqrt(log((meanX.*covX).^2./meanX.^2 + 1));

    R = lognrnd(mu_lognorm, sigma_lognorm, m, n);
    
end