 function mste = st_error(n_sample, statistics)
 % Estimate standard error of a statistic based on the number of 
 % observations and normality assumption
 %
 %SYNOPSYS
 % smte = ST_ERROR(n_sample, statistics)
 %
 %INPUT
 % n_sample     sample size /numeric vector/
 % statistics   statistics for which standard error is to be estimated
 %              'mean', 'std', 'var'
 %
 %OUTPUT
 % mste         std_error = mste*statistics_hat
 %
 % REF:
 % [http://web.eecs.umich.edu/~fessler/papers/files/tr/stderr.pdf]


if nargin < 2
    statistics = 'std';
end

n = n_sample(:);

switch lower(statistics)
    case 'mean'
        mste = 1/sqrt(n_sample);
        
    case 'std'
        K_n = sqrt((n-1)/2).*exp(log(gamma((n-1)/2)) - log(gamma(n/2)));
        V_n = 2*((n-1)/2 - gamma(n/2).^2/(gamma((n-1)/2)).^2);
        mste = K_n.*sqrt(V_n./(n-1));
        
    case {'var', 'variance'}
        mste = sqrt(2/(n-1));
        
    otherwise
        error(['Unknown statistics: ', statistics])
end

end