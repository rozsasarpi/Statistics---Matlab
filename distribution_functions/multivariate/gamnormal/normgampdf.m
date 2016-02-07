% Gamma-normal probability density function
%
% P = NORMGAMPDF(X, mu, lambda, alpha, beta)
%
%INPUT - OUT OF DATE!
% X             - value pair where we are interested in the probability ([Xg, Xn])
% alpha, beta   - gamma distr. parameters (shape and rate parameters, repsectively)
% mu, lambda    - normal distribution parameters (mean value and precision parameters, respectively)
%
%
%OUTPUT
% P     - probability corresponding to value pair [Xg, Xn]
%
%NOTE(S):
% following the notation of [http://en.wikipedia.org/wiki/Normal-gamma_distribution]
% http://www.cs.ubc.ca/~murphyk/Papers/bayesGauss.pdf

function P = normgampdf(X, param, parametrization)
    
    if nargin < 3
        parametrization = 'normgamma';
    end

    Xn = X(:,1);
    Xg = X(:,2);
    
    switch lower(parametrization)
        case {'normgamma', 'normalgamma', 'normgam', 'normalgam'} 
            % notation after [Rackwitz (1983) Predictive distribution of strength under control]
            
            mu  = Xn;
            %h   = Xg;            % if given with precision !!!
            h   = 1./Xg.^2;      % if given with standard deviation !!!
            x_  = param(1);
            s   = param(2);
            n   = param(3);
            nu  = param(4);
            
            P = sqrt(h*n)/sqrt(2*pi).*exp(-1/2*((mu-x_)./(1./sqrt(h*n))).^2).*...
                (1/2*nu*s^2*h).^(1/2*nu-1).*exp(-1/2*nu*s^2*h)/(2*gamma(nu/2))*nu*s^2;
            
%             % parametrization
%             % [http://www.cs.ubc.ca/~murphyk/Teaching/CS340-Fall07/reading/NG.pdf] Section 3.1
%             alpha = 1;
%             beta = 1;
%             mu0 = 0;
%             kappa = 2;
% 
%             x_  = mu0;
%             n   = kappa;
%             s   = sqrt(beta/alpha);
%             nu  = 2*alpha;
            
        case {'lnormgamma', 'lognormalgamma', 'logngamma', 'lognormgamma', 'lnormgam', 'lognormalgam', 'logngam', 'lognormgam'} 
            % notation after [Rackwitz (1983) Predictive distribution of strength under control]
            % SZAAAAR!
            
            mu_lognorm = log(Xn.^2./sqrt((Xg).^2 + Xn.^2));
            sigma_lognorm = sqrt(log((Xg).^2./Xn.^2 + 1));
            
            mu  = mu_lognorm;
            h   = 1./sigma_lognorm.^2;
            x_  = param(1);
            s   = param(2);
            n   = param(3);
            nu  = param(4);
            
            P = sqrt(h*n)/sqrt(2*pi).*exp(-1/2*((mu-x_)./(1./sqrt(h*n))).^2).*...
                (1/2*nu*s^2*h).^(1/2*nu-1).*exp(-1/2*nu*s^2*h)/(2*gamma(nu/2))*nu*s^2;    

        case '??'

            mu = param(1);
            lambda = param(2);
            alpha = param(3);
            beta = param(4);

            %ZNG = gamma(alpha)/(beta^alpha)*(2*pi/kappa
            P = beta.^alpha.*sqrt(lambda)./(gamma(alpha).*sqrt(2*pi)).*Xg.^(alpha-1/2).*exp(-beta.*Xg).*exp(-lambda.*Xg.*(Xn-mu).^2/2);
        otherwise
            error('Unknown normgam parametrization!')
    
    end
    
end