function P = lognorm3cdf(x, shape, scale, thres)
% P = LOGNORM3CDF(x, shape, scale, thres)
%
 
P = logncdf(x-thres, scale, shape);

end
