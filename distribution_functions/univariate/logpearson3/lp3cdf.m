% Log Pearson 3 cumulative distribution function
%
% F = LP3CDF(X, mu, sigma, gam)
%
%INPUT
% X     - value where we are interested in the CDF
% mu    - location parameter
% sigma - scale parameter
% gam   - shape parameter
% 
%OUTPUT
% F     - probability corresponding to that x is smaller than X
%
%REF: Notations follow: Hosking and Wallis. (2013). Regional Frequency Analysis: An Approach Based on L-Moments
%
%NOTE: gam = 0     normal distribution
%      gam = 2     exponential distribution
%      gam = -2    reverse exponential distribution

function F = lp3cdf(X, mu, sigma, gam)

if gam == 0
    %[-inf, inf]
    F       = normcdf(X, mu, sigma);
else
    alpha   = 4/gam^2;
    beta    = 1/2*sigma*abs(gam);
    ksi     = mu - 2*sigma/gam;
    
    if gam > 0
        %[ksi, inf]
        F       = gammainc((X - ksi)/beta, alpha, 'lower');
    elseif gam < 0
        %[-inf, ksi]
        F       = gammainc((ksi - X)/beta, alpha, 'upper');
    end
end

end