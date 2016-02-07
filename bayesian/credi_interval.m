function ci = credi_interval(sample_vec, credi_mass, type)
% Credible interval from a representative sample
%
%USAGE:
% ci = CREDI_INTERVAL(sample_vec, credi_mass, type)
%
%INPUT:
%   sample_vec      Representative values from a probability distribution
%                   (typically MCMC output) /vector/.
%
%   credi_mass      Indicating the mass within the credible interval
%                   that is to be estimated, /real/ from range ]0,1[.
%OPTIONAL:
%   type            type of credible interval (the mass does not uniquely
%                   determine it, [default='hd']
%       'hd'        highest density interval; narrowest credible interval
%                   with the specified mass.
%       'equal'     Probability of being above or below the interval is equal.
%
%OUTPUT:
%   ci              Endpoints of the highest density credible interval.
%                   /vector [1x2]/
%
% inpired by Kurschke's HDI script written in R
%
% See also
% credi_contour, quant_tile

%CUSTOM FUNCTIONS:
% quan_tile.m

if nargin < 3
    type = 'hd';
end

if numel(credi_mass) ~=1 || credi_mass <= 0 || credi_mass >= 1
    error('''credi_max'' should be a real number between 0 and 1!')
end

sample_vec = sample_vec(:);

switch lower(type)
    case {'hd', 'hdi'}
        sorted_sample   = sort(sample_vec);
        ci_idx_inc      = floor(credi_mass.*length(sorted_sample));
        nci             = length(sorted_sample) - ci_idx_inc;
        ci_width        = sorted_sample((1:nci) + ci_idx_inc) - sorted_sample(1:nci);
        
        [~, min_idx]    = min(ci_width);
        ci              = sorted_sample([min_idx, min_idx+ci_idx_inc])';
        
    case {'equal', 'eqi', 'eq'}
        ci = quan_tile(sample_vec, [(1-credi_mass)/2, (1+credi_mass)/2]);
        ci = ci';
    otherwise
        error(['Unknown credible interval type: ', type])
end
