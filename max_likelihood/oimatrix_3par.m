% Numerically constructs the observed information matrix
%
% not general but could be easily extended
% works only with a single point and with trivariate function (currently)
% the derivatives are calculated using central finite differences
% truncation error ~1e-11 (O(h^2))
%
% ref e.g. Finite difference lecture notes. http://www.mathematik.uni-dortmund.de/~kuzmin/cfdintro/lecture4.pdf
%
%See also
% oimatrix

function OI = oimatrix_3par(fun, cx, cy, cz)

% cx = 2;
% cy = 2;
% cz = 2;
% fun3  = @(x,y,z) x*y^2 + x*z^2;

% derivatives
f_xx = cfd_11(@(x) fun(x,cy, cz), cx);
f_yy = cfd_11(@(y) fun(cx,y, cz), cy);
f_zz = cfd_11(@(z) fun(cx,cy, z), cz);

f_xy = cfd_12(@(x,y) fun(x, y, cz), cx, cy);
f_xz = cfd_12(@(x,z) fun(x, cy, z), cx, cz);
f_yz = cfd_12(@(y,z) fun(cx, y, z), cy, cz);

% observed information matrix
OI = -[f_xx, f_xy, f_xz; 
      f_xy, f_yy, f_yz;
      f_xz, f_yz, f_zz];


%% NESTED FUNCTIONS
% second 'pure' derivative
    function nd = cfd_11(fun, x)
        
        %h   = (eps)^(1/3)*arrayfun(@(y) max([y, 1]), abs(x));
        h   = (eps)^(1/4)*max(abs(x), 1);
        xph = x + h;
        xmh = x - h;
        dx  = xph - xmh; % to account for rounding error
        
        nd  = (fun(xph) - 2*fun(x) + fun(xmh))./(1/2*dx)^2;
        
    end

% second mixed derivate
    function nd = cfd_12(fun, x, y)
        
        hx  = (eps)^(1/4)*max(abs(x), 1);
        xph = x + hx;
        xmh = x - hx;
        dx  = xph - xmh; % to account for rounding error
        
        hy  = (eps)^(1/4)*max(abs(y), 1);
        yph = y + hy;
        ymh = y - hy;
        dy  = yph - ymh; % to account for rounding error
        
        nd  = (fun(xph, yph) - fun(xph, ymh) - fun(xmh, yph) + fun(xmh, ymh))./(4*1/2*dx*1/2*dy);
        
    end

end