% Numerically constructs the observed information matrix
% equivalent of the curvature matrix (second derivate approx.) of a function in a given point
%
%SYNOPSYS
% OI = OIMATRIX(fun, x)
%
%INPUT
% fun       function (with vector input) handler /function handler with n-element vector input argument/, loglikelihood
% x         vector containing the maximum likelihood estimates /n-element vector/
%           (coordinates of the point where we are interested in the curvature)
%
%OUTPUT
% OI        observed information matrix (curvature) /nxn matrix/     
%
%CUSTOM FUNCTIONS
%
%
%NOTES
% general formulation, can handle function with arbitrary number of inputs
% num2str and eval with 40 digits precision, quite lame implementation
% the derivatives are calculated using central finite differences
% truncation error ~1e-11 (O(h^2))
%
% ref e.g. Finite difference lecture notes. http://www.mathematik.uni-dortmund.de/~kuzmin/cfdintro/lecture4.pdf
%
%See also
% oimatrix_3par

function OI = oimatrix(fun, x)

% x = [2, 2, 2];
% fun  = @(x) x(1)*x(2)^2 + x(1)*x(3)^2;
% keyboard
% WARNING!
prec    = 30;

% INPUT CHECK, PRE-PROCESS
d       = numel(x);
OI      = zeros(d,d);

% Construct the obsered information matrix
% loop through the rows
for i = 1:d
    % loop through the columns
    for j = i:d
        % Second, 'pure' derivative (diagonal)
        if i == j 
            g_pos = x(i);
            
            %..............................................................
            % lame construction of the input argument for reduced dimension function
            g_arg = '[';
            for k = 1:d
                if k == i
                    g_arg = [g_arg, 'u'];
                else
                    g_arg = [g_arg, num2str(x(k), prec)];
                end
                
                if k ~= d
                    g_arg = [g_arg, ','];
                end
            end
            g_arg = [g_arg, ']'];
            %..............................................................
            
            OI(i,j) = cfd_ii(@(u) fun(eval(g_arg)), g_pos);
            
        % Second, mixed derivate (off-diagonal)    
        else
            g_pos = [x(i), x(j)];
            
            %..............................................................
            % lame construction of the input argument for reduced dimension function
            g_arg = '[';
            for k = 1:d
                if k == i
                    g_arg = [g_arg, 'u'];
                elseif k == j
                    g_arg = [g_arg, 'v'];
                else
                    g_arg = [g_arg, num2str(x(k), prec)];
                end
                
                if k ~= d
                    g_arg = [g_arg, ','];
                end
            end
            g_arg = [g_arg, ']'];
            %..............................................................
                       
            OI(i,j) = cfd_ij(@(u,v) fun(eval(g_arg)), g_pos);
            OI(j,i) = OI(i,j);
        end
    end
end

OI = -OI;

% if any(eig(OI)) < 0
%     warning('The correspondign covariance matrix is negative definite!')
% end

%% NESTED FUNCTIONS
% second, 'pure' derivative
    function nd = cfd_ii(g_fun, pos)
        
        u   = pos;
        
        %h   = (eps)^(1/3)*arrayfun(@(y) max([y, 1]), abs(x));
        h   = (eps)^(1/4)*max(abs(u), 1);
        uph = u + h;
        umh = u - h;
        du  = uph - umh; % to account for rounding error
        
        nd  = (g_fun(uph) - 2*g_fun(u) + g_fun(umh))./(1/2*du)^2;
        
    end

% second, mixed derivate
    function nd = cfd_ij(g_fun, pos)
        
        u   = pos(1);
        v   = pos(2);
        
        hu  = (eps)^(1/4)*max(abs(u), 1);
        uph = u + hu;
        umh = u - hu;
        du  = uph - umh; % to account for rounding error
        
        hv  = (eps)^(1/4)*max(abs(v), 1);
        vph = v + hv;
        vmh = v - hv;
        dv  = vph - vmh; % to account for rounding error
        
        nd  = (g_fun(uph, vph) - g_fun(uph, vmh) - g_fun(umh, vph) + g_fun(umh, vmh))./(4*1/2*du*1/2*dv);
        
    end

end