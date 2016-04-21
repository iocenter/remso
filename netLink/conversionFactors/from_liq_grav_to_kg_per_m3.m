function output = from_liq_grav_to_kg_per_m3(input)
% output = from_liq_grav_to_kg_per_m3(input)
%
% Performs conversion from specific liquid gravity (density
% relative to density of water at standard conditions = 999 kg/m^3 =
% 62.4 lbm/ft^3) to kg/m^3.
% 
% JDJ, 29-11-01.

output =  input * 999.;