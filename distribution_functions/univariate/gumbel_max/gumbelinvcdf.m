%Maximum Gumbel inverse cumulative distribution function
%
%SYNOPSYS
% IFG = GUMBELINVCDF(P, scale_mean, loc_std,, def_type)
%
%INPUT
% P         probability
%
%OPTIONAL
% def_type    paramter definition type
%
% case def_type = 'moments' (default)
% scale_mean    mean
% loc_std       standard deviation
%
% case def_type = 'parameter'
% scale_mean    scale (beta) parameter    
% loc_std       location (mu) parameter
%
%OUTPUT
% IFG
%
%NOTE(S)
% Matlab built in is a minimumum distribution, defined only with parameters

function IFG = gumbelinvcdf(P, scale_mean, loc_std, def_type)
if nargin < 4
    def_type = 'moment';
end

if any(P>1) || any(P<0)
    error('P should be in [0,1]!')
end

gamma = 0.5772156649015328606065120900824024310421;

switch lower(def_type)
    case {'moment', 'moments', 'mom'}
        std_    = loc_std;
        mean_   = scale_mean;
        
        scale   = sqrt(6)/pi*std_;
        loc     = mean_ - scale*gamma;
        
    case {'param', 'par', 'parameter'}
        loc     = loc_std;
        scale   = scale_mean;
end

IFG = -scale.*log(-log(P)) + loc;

end