% Maximum likelihood diagnostic plots, to check if real optimum has found
%
%SYNOPSYS:
% MLE_DIAGNOSTIC_PLOT(param, nloglike)
%
%INPUT:
% param     - maximum likleihood paramter estimates
% nloglike  - handler of the negative loglikelihood function

function mle_diagnostic_plot(param, nloglike)

figure('Position',[100, 200, 1400, 400]);

nLL = nloglike(param);
varnum = numel(param);
% loop through the parameters
for jj = 1:varnum
    subplot(1,varnum,jj)
    xx = 0.98*param(jj):param(jj)/1000:1.02*param(jj);
    par_arg = param;
    
    yy = zeros(size(xx));
    %loop through the ppoint in around the ml estimate
    for ii = 1:numel(xx)
        par_arg(jj) = xx(ii);
        yy(ii) = nloglike(par_arg);
    end
    
    plot(xx,yy)
    hold on
    plot(param(jj), nLL, 'ro')
    xlabel(['param',num2str(jj)])
    
end

end