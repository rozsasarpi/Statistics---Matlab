% Convert Kendall tau to rho parameter for tev copula
%
% delta = TEV_KTAU2RHO(ktauq)

function rho = tev_ktau2rho(ktauq)

load('tev_tau_rho.mat', 'ktau', 'rho')

rho = interp1(ktau, rho, ktauq);

end