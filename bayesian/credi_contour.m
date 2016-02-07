function CI = credi_contour(sample_vec, credi_mass)
% Highest density credible countour (2D credible interval) 
% from a representative 2D sample
%
%USAGE:
% CI = CREDI_CONTOUR(sample_vec, credi_mass)
%
%INPUT:
% sample_vec    Nx2 matrix, N is the number of observation pairs
% credi_mass    scalar [0,1]
%
%OUTPUT:
% CI            contour matrix, /Mx2 matrix/
%
%NOTE:
% Uses 2D kernel estimate for finding the isocontour with _credi_mass_
% probability mass.
%
% Probability mass of the kernel density function above level _p_ is 
% estimated using a simple, grid based, numerical integration is sufficient
% for the mainly illustrative purpose (considering other uncertainties 
% it is reasonable)
%
%EXAMPLE:
% example.1:
%   credi_mass = 0.9;
%   sample_vec = mvnrnd([10,1],[1,0.5; 0.5,1],1e3);
%   CI = credi_contour(sample_vec, credi_mass);
%   plot(sample_vec(:,1),sample_vec(:,2),'r.','MarkerSize',5)
%   hold on
%   plot(CI(:,1), CI(:,2), 'black')
%
% example.2:
%   credi_mass = 0.9;
%	sample_vec = mvnrnd([1,1],[1,0.5; 0.5,1],1e3) + mvnrnd([10,10],[1,0.5; 0.5,1],1e3);
%   CI = credi_contour(sample_vec, credi_mass);
%   plot(sample_vec(:,1),sample_vec(:,2),'r.','MarkerSize',5)
%   hold on
%   plot(CI(:,1), CI(:,2), 'black')
%
% TODO: more roboust _p0_ estimate
% TODO: improve input check
%
%CUSTOM FUNCTIONS:
% kde2d.m

if nargin < 2
    error('Two input arguments are required!')
end

if size(sample_vec,1) < 3 || size(sample_vec,2) ~= 2
    error('''sample_vec'' should be in the form of Nx2 !')
end

% Very fast 2D kernel density estimation by Botev
[~,density,X,Y] = kde2d(sample_vec,2^8);

% grid cell dimensions
dx = diff(X(1,1:2));
dy = diff(Y(1:2,1));

% initial value
p0  = 0.5*max(max(density));

opt = optimset('fzero');
opt = optimset(opt, 'TolX',p0*1e-3);
opt = optimset(opt, 'Display','notify');

x = fzero(@(p) sum(sum(density(density > p)))*dx*dy - credi_mass, p0, opt);

C = contourc(X(1,:), Y(:,1), density,[x, x]);

CI = [C(1,2:end)', C(2,2:end)'];

end



