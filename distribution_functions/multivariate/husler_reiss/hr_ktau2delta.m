% Convert Kendall tau to delta parameter for Hüsler-Reiss copula
%
% delta = HR_KTAU2DELTA(ktauq)

function delta = hr_ktau2delta(ktauq)

load('hr_tau_delta.mat', 'ktau', 'delta')

delta = interp1(ktau, delta, ktauq);

end