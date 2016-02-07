% GUMBEL_STAT Conversion between parameters and ~moments
% of max Gumbel probability distribution
%
%SYNOPSYS
% parm = GUMBEL_STAT(scale_mean, loc_std, def_type)
%
% if def_type == 'par' (default)
% parm = [mean, std]
%
% if def_type == 'mom'
% parm = [scale, loc]
%
% See also
% gev_stat, gevstat

function parm = gumbel_stat(scale_mean, loc_std, def_type)

gamma = 0.5772156649015328606065120900824024310421;

switch lower(def_type)
    case {'p', 'par', 'parm', 'parameter', 'parameters'}
        loc     = loc_std;
        scale   = scale_mean;
        
        std_    = pi/sqrt(6)*scale;
        mean_   = loc + scale*gamma;
        parm    = [mean_, std_];
        
    case {'m', 'mom', 'moment', 'moments'}
        std_    = loc_std;
        mean_   = scale_mean;
        
        scale   = sqrt(6)/pi*std_;
        loc     = mean_ - scale*gamma;
        parm    = [scale, loc];
        
    otherwise
        error(['Unknown def_type: ', def_type])
end