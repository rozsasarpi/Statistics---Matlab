function R = lognorm3rnd(shape, scale, thres, m, n)
% R = LOGNORM3RND(shape, scale, thres, m, n)

RN = normrnd(0,1,m,n);

R = exp(RN*shape+scale) + thres;

end