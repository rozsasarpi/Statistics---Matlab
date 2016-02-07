%TRUNCPPDF Truncated Pareto probability density function
%
%SYNOPSYS
% f = TRUNCPPDF(x, exponent, lowerlim, upperlim)
%
%INPUT
% x                 - value where we are interested in the probability density /vector/
%
% exponent          - a
%
% lowerlim          - L
%
% upperlim          - U
%
%OUTPUT
% f                 - probability corresponding to value x /vector/
%
%NOTES:
%
% (1) model from 2010_Estimation of Pareto Distribution Functions from Samples Contaminated by Measurement Errors. MSc thesis
%     + http://en.wikipedia.org/wiki/Pareto_distribution#Bounded_Pareto_distribution

function f = truncppdf(x, exponent, lowerlim, upperlim)

if nargin < 4
    error('Too few input arguments!');
end

a = exponent;
L = lowerlim;
U = upperlim;

if any(x < L) || any(x > U)
    error('Values of ''x'' should be between lowerlim (L) and upperlim (U)!')
end

f = a.*x.^(-a-1)./(L.^-a - U^-a);

end