% Log Pearson 3 distributed random numbers
%
% R = LP3RND(mu, sigma, gam, n, m)
%
%INPUT
% mu    - location parameter
% sigma - scale parameter
% gam   - shape parameter
% n     - number of random number generated (row)
% m     - number of random number generated (column)
% 
%OUTPUT
% n X m - matrix of random number per maximum Log Pearson 3 distribution
%
%REF: Notations follow: Hosking and Wallis. (2013). Regional Frequency Analysis: An Approach Based on L-Moments
%
%NOTE: gam = 0     normal distribution
%      gam = 2     exponential distribution
%      gam = -2    reverse exponential distribution

function R = lp3rnd(mu, sigma, gam, n, m)

if gam == 0
    %[-inf, inf]
    R       = normrnd(mu, sigma, n, m);
else
    alpha   = 4/gam^2;
    beta    = 1/2*sigma*abs(gam);
    ksi     = mu - 2*sigma/gam;
    
    RU      = rand(n, m);
    if gam > 0
        %[ksi, inf]
        RU      = rand(n, m);
        R       = gammaincinv(RU, alpha, 'lower')*beta + ksi;
    elseif gam < 0
        %[-inf, ksi]
        R       = -gammaincinv(RU, alpha, 'upper')*beta + ksi;
    end
end

end