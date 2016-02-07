% transform numbers form [0,1] to *^10 scale symmetrically
%
%SYNOPSYS
% tx = LOG10_TRANS(x)
%
% stretches(zoom) the vicintiy of 1 and 0, useful for visualizing the tails
% of cumulative distribution function

function tx = log10_trans(x)

tx          = nan(size(x));
idx         = x < 0.5;

tx(idx)     = log10(x(idx)*2);
tx(~idx)    = -log10(-2*(x(~idx) - 1));

end