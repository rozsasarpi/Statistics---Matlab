% Probability distribution function (pdf) of the product of two independent random variables
%
% Y = X1*X2
%
%SYNOPSYS:
% pdy = PROD_2RV(pd1, pd2)
%
%INPUT:
% pd1           structure referring to rv X1
%  .fx_fun      handler of X1's pdf
%
% pd2           structure referring to rv X2, defined the same way as pd1
%
%OUTPUT:
% pdy           
%  .fx_fun      handler of Y's pdf
%
% Only calculates a single point, hence for many points it may be ineffective

function pdy = prod_2rv(pd1, pd2)

fy = @(y) integral(@(x) pd1.fx_fun(x).*pd2.fx_fun(y./x).*1./abs(x), pd1.min, pd1.max, 'AbsTol',1e-30, 'RelTol',1e-12); % WARNING! limits!

pdy.fx_fun = fy;

% Mellin transformation could be the solution (analogous to Fourier transform for sum), My = Mx1*Mx2!
% need a fast algorithm that can calculate discrete Mellin transformation

end

%EXAMPLE 1
% clc
% close all
% sigma = 0.1;
% pd1.fx_fun = @(x) lognpdf(x,1,0.1);
% pd1.min = 0;
% pd1.max = 20;
% pd2.fx_fun = @(x) lognpdf(x,1,sigma);
% 
% pdy = prod_2rv(pd1, pd2);
% 
% x = 0:0.1:20;
% 
% tic
% y = zeros(numel(x),1);
% for i = 1:numel(x)
%     y(i) = pdy.fx_fun(x(i));
% end
% toc
% 
% pd3.fx_fun = @(x) lognpdf(x,2,sqrt(0.1^2+sigma^2));
% plot(x, pd3.fx_fun(x))
% hold on
% plot(x,y,'--r')
% plot(x, pd1.fx_fun(x),'--g')

%EXAMPLE 2
% clc
% close all
% sigma = 0.001;
% pd1.fx_fun = @(x) lognpdf(x,1,0.1);
% pd1.min = 0;
% pd1.max = 20;
% pd2.fx_fun = @(x) normpdf(x,1,sigma);
% 
% pdy = prod_2rv(pd1, pd2);
% 
% x = 0:0.05:5;
% 
% tic
% y = zeros(numel(x),1);
% for i = 1:numel(x)
%     y(i) = pdy.fx_fun(x(i));
% end
% toc
% 
% plot(x,y,'r')
% hold on
% plot(x, pd1.fx_fun(x),'--b')