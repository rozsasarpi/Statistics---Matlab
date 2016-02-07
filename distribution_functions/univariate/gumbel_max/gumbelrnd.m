%maximum Gumbel distribution
%random numbers with Gumbel(max) distribution
%
% R = GUMBELRND(meanX, stdX, n, m)
%
%INPUT
% meanX - mean value
% stdX  - standard deviation
% n     - number of random number generated (row)
% m     - number of random number generated (column)
%
%OUTPUT
% n X m - matrix of random number per maximum Gumbel distribution

function R = gumbelrnd(meanX, stdX, n, m)

gamma   = 0.5772156649015328606065120900824024310421;
beta    = sqrt(6)/pi*stdX;
mupar   = meanX - beta*gamma;

P = rand(n,m);

R = -beta*log(-log(P)) + mupar;

end