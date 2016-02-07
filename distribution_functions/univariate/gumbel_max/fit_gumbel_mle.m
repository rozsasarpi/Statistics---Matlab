% FIT_GUMBEL_MLE Max Gumbel, maximum likelihood parameter estimation
%
%SYNOPSYS
% [parmhat, max_nLL, output] = FIT_GUMBEL_MLE(data, def_type)
%
% if def_type == 'par' (default)
% parmhat = [scale_hat, loc_hat]
%
% if def_type == 'mom'
% parmhat = [mean_hat, std_hat]
%
% See also
% gevfit, evdfit, gumbelpdf

function [parm_hat, max_nLL, output] = fit_gumbel_mle(data, def_type)

if nargin < 2
    def_type = 'par';
end

gamma = 0.5772156649015328606065120900824024310421;

% starting point is the method of moment value
meanX   = mean(data);
stdX    = std(data);
m       = length(data);

switch lower(def_type)
    case {'p', 'par', 'parm', 'parameter', 'parameters'}
        
        scale   = sqrt(6)/pi*stdX;
        loc     = meanX - scale*gamma;
        mom_par = [scale, loc];
        % refine options if needed
        options = optimset('MaxFunEvals', 2e3, ...
            'MaxIter', 2e3, ...
            'TolFun', abs(nloglike_par(mom_par))*1e-6, ...
            'TolX', min(abs(mom_par))*1e-6);
        
        [parm_hat, max_nLL, ~, output] = fminsearch(@nloglike_par, mom_par, options);
        
    case {'m', 'mom', 'moment', 'moments'}
        
        mom_par = [meanX, stdX];
        % refine options if needed
        options = optimset('MaxFunEvals', 2e3, ...
            'MaxIter', 2e3, ...
            'TolFun', abs(nloglike_mom(mom_par))*1e-6, ...
            'TolX', min(abs(mom_par))*1e-6);
        
        [parm_hat, max_nLL, ~, output] = fminsearch(@nloglike_mom, mom_par, options);
        
    otherwise
        error(['Unknown def_type: ', def_type])
end

% negative log-likelihood function
    function nLL = nloglike_par(param)
        scale_  = param(1);
        loc_    = param(2);
        
        z       = (data - loc_)/scale_;
        
        LL      = m*log(1/scale_) - sum(z + exp(-z));
        nLL     = -LL;
    end

    function nLL = nloglike_mom(param)
        mean_   = param(1);
        std_    = param(2); 
        scale_  = sqrt(6)/pi*std_;
        loc_    = mean_ - scale_*gamma;
        
        z       = (data - loc_)/scale_;
        
        LL      = m*log(1/scale_) - sum(z + exp(-z));
        nLL     = -LL;
    end

end