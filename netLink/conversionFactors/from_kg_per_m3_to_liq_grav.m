function output = from_kg_per_m3_to_liq_grav(input)
% output = from_kg_per_m3_to_liq_grav(input)
%
% Performs conversion from kg/m^3 to specific liquid gravity (density
% relative to density of water at standard conditions = 999 kg/m^3 =
% 62.4 lbm/ft^3).
%  
% JDJ, 29-11-01.

output =  input / 999.;