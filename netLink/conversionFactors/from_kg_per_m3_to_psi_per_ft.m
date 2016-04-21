function output = from_kg_per_m3_to_psi_per_ft(input)
% output = from_kg_per_m3_to_psi_per_ft(input)
%
% Performs unit conversion from kg/m^3 to psi/ft.
%
% JDJ, 29-11-01.

output =  input * 9.80665 / 2.262059e4;