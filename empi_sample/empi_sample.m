% Random sampling from the empirical distribution function of a sample
%
%SYNOPSYS
% R = EMPI_SAMPLE(obs, n)
%
%INPUT
% obs   sample, vector of observations
% n     number of random numbers required
%
%OUTPUT
% R     random sample, vector of random elements from empirical dist.
%
%NOTE(S)
% (1)The empirical distribution function is not unique, herein
% empiF_i = (i - 2/5)/(k + 1/5)
% formula is used, where 'k' is the length of 'obs'
% (2) Between the points determined with the above formula, linear interpolation is used.
% (3) Linear extrapolation also added!
%
%TO DO:
% adding other methods, e.g. kernel distr. appr to sample from
%
%See also
% ecdf, rand, interp1

function R = empi_sample(obs, n)

xx      = sort(obs);
k       = length(obs);

empiF   = ((1:k) - 2/5)/(k + 1/5);
r_uni   = rand(n,1);
R       = interp1(empiF, xx ,r_uni, 'linear', 'extrap');

end
