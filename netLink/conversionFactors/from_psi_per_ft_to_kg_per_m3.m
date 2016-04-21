function output = from_psi_per_ft_to_kg_per_m3(input)
% output = from_psi_per_ft_to_kg_per_m3(input)
%
% Performs unit conversion from psi/ft to kg/m^3.
%
% JDJ, 29-11-01.

output =  input * 2.262059e4 / 9.80665;