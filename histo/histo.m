function [nelements,centers] = histo(data,nbins)
%Creates histogram from sample
%
% Accomplishes the same task as _hist(data,nbins)_
% but accepts only vector input!
%
%USAGE:
% [nelements,centers] = HISTO(data,nbins)
%
%INPUT:
% data          Data to distribute among bins. /vector/
%
%OPTIONAL:
% nbins         Number of bins. /integer/ [default=10]
%
%OUTPUT:
% nelements     Number of elments in each bin. /vector/
% centers       Bin centers. /vector/
%
%
%NOTES:
% Simple implementation to remove the dependence on toolboxes.
% Contrary to _hist_ this function should be closed by _;_ even if not output
% arguments are specified.
%
% See also
% bar

if all(size(data) ~= 1) || numel(size(data)) > 2
    error('''data'' should be a vector!')
end

if nargin < 2
    nbins = 10;%min(10,numel(data));
end

% reshape to column vector
data = data(:);

% Find max and min
idx  = isfinite(data);
xmin = min(data(idx));
xmax = max(data(idx));

% bin width
wbin = (xmax-xmin)/(nbins-1);

% xmax and xmin is set to be the center of the edge bins
centers = xmin + (0:nbins-1)*wbin;
edges   = [centers-wbin/2; centers+wbin/2];

% ]a, b] type counting
% (a, b]
nelements = zeros(1,nbins);
% loop through the bins; thanks JIT
for ii = 1:nbins
    nelements(ii) = sum(data >= edges(1,ii) & data < edges(2,ii));
end

% plot histogram if there is no output argument
if nargout == 0
    bar(centers, nelements,1);
end

end