%Maximum Gumbel probability density function
%
% fG = GUMBELPDF(x, meanX, stdX)
%
% WARNING! should be changed to this:
%SYNOPSYS
% fG = GUMBELPDF(P, p1, p2, method)
%
%INPUT
% x         
%
%OPTIONAL
% method    paramter definition type
%
% case method = 'moments' (default)
% p1        mean
% p2        standard deviation
%
% case method = 'parameter'
% p1        mu paramter
% p2        beta parameter
%
%OUTPUT
% fG        probability density at x      
%
%NOTE(S)
% Matlab built in is a minimumum distribution, defined with parameters

function fG = gumbelpdf(x, meanX, stdX)

%maximum Gumbel distribution
%matlab built in is a minimumum distribution
gamma = 0.5772156649015328606065120900824024310421;
beta = sqrt(6)/pi*stdX;
mupar = meanX - beta*gamma;

z = (x-mupar)./beta;

fG = 1./beta.*exp(-z-exp(-z));

end