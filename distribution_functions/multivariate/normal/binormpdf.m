% Bivariate normal probability density function (pdf)
%
%SYNOPSYS
%   y = BINORMPDF(X, MU, COV)
%
%INPUT
%
%
%OUTPUT
%
%

function y = binormpdf(X, MU, COV)

if COV(1,2) ~= COV(2,1)
    error('The covariance matrix /COV/ must be symmetric!')
end

if size(X,2) ~= 2
    if numel(X) == 2
        X = X';
    else
        error('Wrong dimension(s) for X!')
    end
end

x1      = X(:,1);
x2      = X(:,2);
mu1     = MU(1);
mu2     = MU(2);

sigma1  = sqrt(COV(1,1));
sigma2  = sqrt(COV(2,2));
rho     = COV(1,2)/(sigma1*sigma2);

y = 1/(2*pi*sigma1*sigma2*sqrt(1-rho^2))*exp(-1/(2*(1-rho^2))*...
    ((x1-mu1).^2/sigma1^2 + (x2-mu2).^2/sigma2^2 - 2*rho*(x1-mu1).*(x2-mu2)/(sigma1*sigma2)));
end