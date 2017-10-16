function output = from_kg_per_m3_to_gas_grav(input)
% output = from_kg_per_m3_to_gas_grav(input)
%
% Performs conversion from kg/m^3 to specific gas gravity (density
% relative to density of air at standard conditions = 1.23 kg/m^3 =
% 76.3 * 10^-3 lbm/ft^3).
%  
% JDJ, 29-11-01.

output =  input / 1.23;