% Log Pearson 3 probability density functioon
%
% f = LP3PDF(X, mu, sigma, gam)
%
%INPUT
% X     - value where we are interested in the probability density
% mu    - location parameter
% sigma - scale parameter
% gam   - shape parameter
% 
%OUTPUT
% f     - probability corresponding to value X
%
%REF: Notations follow: Hosking and Wallis. (2013). Regional Frequency Analysis: An Approach Based on L-Moments
%
%NOTE: gam = 0     normal distribution
%      gam = 2     exponential distribution
%      gam = -2    reverse exponential distribution

function f = lp3pdf(X, mu, sigma, gam)

if gam == 0
    %[-inf, inf]
    f       = normpdf(X, mu, sigma);
else
    alpha   = 4/gam^2;
    beta    = 1/2*sigma*abs(gam);
    ksi     = mu - 2*sigma/gam;
    
    if gam > 0
        %[ksi, inf]
        f       = (X-ksi).^(alpha-1).*exp(-(X-ksi)/beta)/(beta^alpha*gamma(alpha));
    elseif gam < 0
        %[-inf, ksi]
        f       = (ksi-X).^(alpha-1).*exp(-(ksi-X)/beta)/(beta^alpha*gamma(alpha));
    end
end

end