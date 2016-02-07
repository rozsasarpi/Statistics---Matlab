% Gumbel (max) cumulative distribution function
%
% P = GUMBELCDF(X, meanX, covX)
%
%INPUT
% X     - value where we are interested in the probability
% meanX - mean value of random variable X
% covX  - coefficient of variation of random variable X
% 
%OUTPUT
% P     - probability corresponding to value X
%
% maximum Gumbel distribution
% cumulative distribution function
% matlab built-in is a minimumum distribution!

function P = gumbelcdf(x, meanX, covX)


gamma = 0.5772156649015328606065120900824024310421;
beta = sqrt(6)./pi*covX.*meanX;
mupar = meanX - beta.*gamma;

P = exp(-exp(-(x-mupar)./beta));

end