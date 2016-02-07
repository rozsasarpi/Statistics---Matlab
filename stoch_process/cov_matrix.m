% Generates covariance matrix points given with their coordinates (1D)  
% (random field or stochastic process, same rationale)
%
%SYNOPSIS
% CM = COV_MATRIX(X, Corr)
%
%INPUT
% X             - coordinates of the nodes /vector; nx1/
% Corr.         - structure defining the autocorrelation function
%   .length     - correlation length of the random variable/field /scalar/
%   .type       - type of the autocorrelation function /char; 'exp',
%                 'cauchy'/ (default = 'exp')
%
%OPTIONAL
% Corr.
%   .pow        - power/exponent, not applicabel for exp and gaussian types
%                 (larger means faster decrease in
%                 correlation with distance), /scalar/ (default = 2)
%
%OUTPUT
% CM            - covariance matrix of the random field in the nodes'
%                 position
%                 /matrix, nxn/
%COMMENTS
% (1) NOTE: The script uses Statistical Toolbox's functions (pdist, squareform).

function CM = cov_matrix(X, Corr)

%==========================================================================
% INPUT CHECK & INITIALIZATION
%==========================================================================
if nargin < 2
    error('Correlation structure (Corr) should be provided!')
end

if ~isfield(Corr, 'length')
    error('Correlation length (Corr.length) should be specified!')
end

% set default value if not specified
if ~isfield(Corr, 'pow')
    Corr.pow  = 2;
end

% set default value if not specified
if ~isfield(Corr, 'type')
    Corr.type  = 'exponential';
end

corr_length = Corr.length;
corr_pow    = Corr.pow;
corr_type   = Corr.type;

X = X(:);

%==========================================================================
% CALCULATION
%==========================================================================
% Eucledian distance of the nodes
dx = squareform(pdist(X,'euclidean'));

% covariance matrix (https://cran.r-project.org/web/packages/geoR/geoR.pdf)
switch lower(corr_type)
    case {'e', 'exp', 'exponential'}
        CM = exp(-dx./corr_length);
    case {'c', 'cauchy'}
        CM = (1 + (dx./corr_length).^2).^(-corr_pow);
    case {'g', 'gauss', 'gaussian'}
        CM = exp(-(dx./corr_length).^2);
    otherwise
        error(['Unknown correlation type /corr_type/: ', corr_type])
end

end