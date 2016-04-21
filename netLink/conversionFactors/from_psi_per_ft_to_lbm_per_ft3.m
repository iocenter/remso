function output = from_psi_per_ft_to_lbm_per_ft3(input)
% output = from_psi_per_ft_to_lbm_per_ft3(input)
%
% Performs conversion from gradient in psi/ft to
% density in lbm/ft^3 (exact).
%
% JDJ, 03-12-01, corrected 06-05-03 (thanks to JRR).

output =  input * 144;