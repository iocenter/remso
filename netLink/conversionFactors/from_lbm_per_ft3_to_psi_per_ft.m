function output = from_lbm_per_ft3_to_psi_per_ft(input)
% output = from_lbm_per_ft3_to_psi_per_ft(input)
%
% Performs conversion from density in lbm/ft^3 to
% gradient in psi/ft (exact).
%
% JDJ, 03-12-01, corrected 06-05-03 (thanks to JRR).

output =  input / 144;