% Probability distribution function (pdf) of the sum of two independent random variables (X1, X2)
%
% Y = X1 + X2
% Convolution of two random variables -> performed using discrete Fourier
% transformation
%
%SYNOPSYS:
% pdy = SUM_2RV(pd1, pd2, n, sign)
%
%INPUT:
% pd1           structure corresponding to rv X1
%  .fx_fun      handler of X1's pdf
%  .mean        mean value of X1
%  .min         'min' of X1 (to cover almost entirely the pdf's area)
%  .max         'max' of X1 (to cover almost entirely the pdf's area) 
%
% pd2           structure corresponding to rv X2, defined the same way as pd1
%
%OPTIONAL:
% n             number of point used in discrete Fourier transformation (default=1e3)
% sign          'addition' or 'subtraction' of X1, X2 (default='addition')
%
%OUTPUT:
% pdy
%
%EXAMPLE:
%
% pd1.fx_fun    = @(x) normpdf(x,1,1);
% pd1.mean      = 1;
% pd1.min       = 1-6*1;
% pd1.max       = 1+6*1;
%
% pd2.fx_fun    = @(x) normpdf(x,2,1);
% pd2.mean      = 2;
% pd2.min       = 2-6*1;
% pd2.max       = 2+6*1;
%
% pdy           = SUM_2RV(pd1, pd2);
% pdy           = SUM_2RV(pd1, pd2, 1e2);
% pdy           = SUM_2RV(pd1, pd2, [], 'subtraction');
%
%NOTE(S):
% Y = X1 + X2
% the pdf of Y is shifted to the correct position based on the equivalence of mean values:
% E[Y] = E[X1 + X2]
% this requires that the rvs are given in a sufficently wide domain (covers almost all the densities)
% it would be great to be able to extract the shift without req of 'full' domain
% the sufficiency of the [min, max] domain is checked internally!
%
% FAR FROM BEING GENERAL OR COMPLETE!
%
%See also
% prod_2rv, conv

function pdy = sum_2rv(pd1, pd2, n, sign)
    
    % pre-processing, input check
    if nargin < 4
        sign = 'add';
    end
    if nargin < 3
        n = 1e3;
    end
    
    if pd1.min > pd1.max || pd2.min > pd2.max
        error('The lower bound (pd.min) of the support should be smaller than the upper bound (pd.max)!')
    end
    
    switch lower(sign)
        case {'a','add','addition'}
            % do nothing
        case {'s','sub','subtract','subtraction'} %in case of sign = 'subtract' the order of rvs is important! X1-X2
            pd2.fx_fun  = @(x) pd2.fx_fun(-x);
            pd2.mean    = -pd2.mean;
            temp        = pd2.max;
            pd2.max     = -pd2.min;
            pd2.min     = -temp;
        otherwise
            error(['Unknown ''sign'': ', sign])
    end
    
    % discretize the wider min-max range with n points - might not be sufficient for the smaller one!
    if abs(pd1.max - pd1.min) > abs(pd2.max - pd2.min)
        x1      = linspace(pd1.min, pd1.max, n);
        dx      = abs(diff(x1(1:2)));
        x2      = pd2.min:dx:pd2.max;
    else
        x2      = linspace(pd2.min, pd2.max, n);
        dx      = abs(diff(x2(1:2)));
        x1      = pd1.min:dx:pd1.max;        
    end

    % discrete points of the pdfs
    fx1     = pd1.fx_fun(x1);
    fx2     = pd2.fx_fun(x2);
    
    % convolution
    fy      = dx*conv(fx1, fx2, 'full');
    xx      = (0:numel(fy)-1)*dx;    

    % check the area under Y's pdf
    A       = trapz(xx,fy);
    if A < 1-1e-3
        warning(['The area under the pdf of the sum is considerably smaller than 1! A = ' num2str(A),...
            '. The [min, max] range or number of points (n) should be modified!'])
    end
    
    % shift to the correct position of Y's pdf
    mean_y  = trapz(xx,fy.*xx)/A;
    shift   = mean_y - (pd1.mean + pd2.mean);
    
    %WARNING!
    %####################################################
    if any(isnan(fy)) || any(isnan(xx)) || isnan(shift)
        warning('NaNs produced during convolution! The output is a ~0, flat pdf.')
        pdy.fx_fun  = @(y) ones(size(y))*1e-18;
    else
    %####################################################
    
    % function of Y's pdf
    fy_fun  = @(y) interp1(xx-shift, fy, y);
    
    % output /the min max might be refined, with small quantiles/
    pdy.mean    = mean_y-shift;
    pdy.min     = min(xx-shift);
    pdy.max     = max(xx-shift);
    pdy.fx_fun  = fy_fun;
    
    end
    
end