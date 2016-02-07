% Log Pearson 3 inverse cumulative distribution function
%
% X = LP3INV(P, mu, sigma, gam)
%
%INPUT
% P     - probability where we are interested in the value of X
% mu    - location parameter
% sigma - scale parameter
% gam   - shape parameter
% 
%OUTPUT
% X     - corresponding to P probability
%
%REF: Notations follow: Hosking and Wallis. (2013). Regional Frequency Analysis: An Approach Based on L-Moments
%
%NOTE: gam = 0     normal distribution
%      gam = 2     exponential distribution
%      gam = -2    reverse exponential distribution

function X = lp3inv(P, mu, sigma, gam)

if gam == 0
    %[-inf, inf]
    X       = norminv(P, mu, sigma);
else
    alpha   = 4/gam^2;
    beta    = 1/2*sigma*abs(gam);
    ksi     = mu - 2*sigma/gam;
    
    if gam > 0
        %[ksi, inf]
        X       = gammaincinv(P, alpha, 'lower')*beta + ksi;
    elseif gam < 0
        %[-inf, ksi]
        X       = -gammaincinv(P, alpha, 'upper')*beta + ksi;
    end
end

end